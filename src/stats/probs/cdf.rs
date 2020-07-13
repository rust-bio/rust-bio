// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Support for discrete probability distributions in terms of cumulative distribution
//! functions (CDF).
//!
//! # Examples
//!
//! Example usage of all CDF functions:
//! ```
//! use approx::assert_relative_eq;
//! use bio::stats::probs::cdf::{Entry, CDF};
//! use bio::stats::probs::{LogProb, Prob};
//! use ordered_float::NotNan;
//! use std::ops::Range;
//! // pmf1 is an example PMF with `LogProb(0.0)` at `0`, with `LogProb(0.1)`
//! // at `{1, 2, ..., 8}` and LogProb(0.2) at `10`
//! let mut pmf1 = vec![Entry::new(0, LogProb((0.0 as f64).ln()))];
//! for i in 1..=8 {
//!     pmf1.push(Entry::new(i, LogProb((0.1 as f64).ln())));
//! }
//! pmf1.push(Entry::new(10, LogProb((0.2 as f64).ln())));
//!
//! // create the cumulative distribution function from the probability mass function
//! let cdf = CDF::from_pmf(pmf1.clone());
//! assert_relative_eq!(*cdf.get(&0).unwrap(), (0.0 as f64).ln(), epsilon = 0.0);
//! assert_relative_eq!(
//!     *cdf.get(&3).unwrap(),
//!     (0.3 as f64).ln(),
//!     epsilon = 0.0000000001
//! );
//!
//! // get back the original probability mass value at 7
//! assert_relative_eq!(
//!     *cdf.get_pmf(&7).unwrap(),
//!     (0.1 as f64).ln(),
//!     epsilon = 0.00001
//! );
//!
//! // Check that cdf sums up to 1.0
//! assert_relative_eq!(
//!     f64::from(cdf.total_prob()),
//!     (1.0 as f64).ln(),
//!     epsilon = 0.0
//! );
//!
//! // copy a CDF via its iter() function
//! let mut cdf_copy = CDF::from_cdf(cdf.iter().cloned());
//! assert_eq!(cdf.len(), cdf_copy.len());
//!
//! // get the maximum a posteriori probability estimate
//! assert_eq!(cdf_copy.map().unwrap(), &10);
//!
//! // get the 50% credible interval
//! assert_eq!(cdf_copy.credible_interval(0.5).unwrap(), &2..&8);
//!
//! // cdf_vec is an example Entry vector with `LogProb(0.0)` at `ordered_float::NotNan`
//! // values `{0.0, 1.0, 2.0}` and increasing by `LogProb(0.2)` at each to `{3.0, 4.0, ..., 7.0}`
//! let mut cdf_vec = Vec::new();
//! for i in 0..=2 {
//!     cdf_vec.push(Entry::new(
//!         NotNan::new(i as f64).unwrap(),
//!         LogProb::ln_zero(),
//!     ))
//! }
//! for i in 3..=7 {
//!     cdf_vec.push(Entry::new(
//!         NotNan::new(i as f64).unwrap(),
//!         LogProb(((i - 2) as f64 * 0.2f64).ln()),
//!     ));
//! }
//!
//! // create cdf from vector of `Entry`s
//! let mut cdf_from_vec = CDF::from_cdf(cdf_vec.into_iter());
//!
//! assert_relative_eq!(
//!     *cdf_from_vec.get(&NotNan::new(2.0).unwrap()).unwrap(),
//!     LogProb::ln_zero(),
//!     epsilon = 0.0
//! );
//! assert_relative_eq!(
//!     *cdf_from_vec.get(&NotNan::new(4.0).unwrap()).unwrap(),
//!     LogProb((0.4 as f64).ln()),
//!     epsilon = 0.0
//! );
//!
//! // get the number of `Entry`s in cdf_from_vec
//! assert_eq!(cdf_from_vec.len(), 8);
//!
//! // remove three zero values at `{0.0, 1.0, 2.0}` with `CDF::reduce()`
//! cdf_from_vec = CDF::reduce(cdf_from_vec);
//! assert_eq!(cdf_from_vec.len(), 5);
//! ```

use std::f64;
use std::iter;
use std::ops::Range;
use std::slice;

use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::stats::LogProb;

/// An `Entry` associates a `LogProb` with a value on an ordered axis. It can for example be
/// used to set up probability mass functions or cumulative distribution functions ([CDF](struct.CDF)).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Entry<T: Ord> {
    /// A `value` on the ordered axis, which has to have the Trait [`std::cmp::Ord`](https://doc.rust-lang.org/std/cmp/trait.Ord.html) implemented.
    pub value: T,
    /// A probability at that `value` / point x on the x-axis.
    pub prob: LogProb,
}

impl<T: Ord> Entry<T> {
    /// Create a new `Entry` for `prob` at `value`.
    ///
    /// `value` needs to have the Trait [`std::cmp::Ord`](https://doc.rust-lang.org/std/cmp/trait.Ord.html)
    /// implemented. As `f64` only has PartialOrd, use something like the [`ordered_float` crate](https://docs.rs/ordered-float/1.0.2/ordered_float/)
    /// if you want to use floating point numbers.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::stats::probs::cdf::Entry;
    /// use bio::stats::LogProb;
    /// let entry = Entry::new(5, LogProb(0.6));
    /// assert_eq!(entry.value, 5);
    /// assert_eq!(entry.prob, LogProb(0.6));
    /// ```
    pub fn new(value: T, prob: LogProb) -> Self {
        Entry { value, prob }
    }
}

/// Implementation of a cumulative distribution function as a vector of `Entry`s.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CDF<T: Ord> {
    inner: Vec<Entry<T>>,
}

impl<T: Ord> CDF<T> {
    /// Create CDF from a vector representing a probability mass function (PMF).
    /// The PMF may contain duplicate values the probabilities of which are summed
    /// during generation of the CDF.
    ///
    /// Runtime complexity: O(n log n), where n is the number of `entries`.
    ///
    /// # Arguments
    ///
    /// * `entries` - The PMF as a vector of `Entry` objects (values with an associated `LogProb`).
    pub fn from_pmf(mut entries: Vec<Entry<T>>) -> Self {
        entries.sort_by(|a, b| a.value.cmp(&b.value));
        let mut inner: Vec<Entry<T>> = Vec::new();
        for mut e in entries {
            let p = inner
                .last()
                .map_or(LogProb::ln_zero(), |e| e.prob)
                .ln_add_exp(e.prob);
            if !inner.is_empty() && inner.last().unwrap().value == e.value {
                inner.last_mut().unwrap().prob = p;
            } else {
                e.prob = p;
                inner.push(e);
            }
        }
        let mut cdf = CDF { inner };

        // cap at prob=1.0 if there are slightly exceeding values due to numerical issues.
        for e in &mut cdf.inner {
            e.prob = e.prob.cap_numerical_overshoot(0.00001);
        }

        cdf
    }

    /// Create CDF from iterator. This can be used to replace the values of a CDF.
    ///
    /// Runtime complexity: O(n), where n is the number of `entries`.
    ///
    /// # Arguments
    ///
    /// * `entries` - An iterator over `Entry<T>` values, where T requires
    pub fn from_cdf<I: Iterator<Item = Entry<T>>>(entries: I) -> Self {
        CDF {
            inner: entries.collect_vec(),
        }
    }

    /// Reduce CDF by omitting values with zero probability.
    ///
    /// Runtime complexity: O(n), where n is the number of `value`s with `prob` of zero.
    pub fn reduce(self) -> Self {
        let mut inner = Vec::new();
        let mut last = LogProb::ln_zero();
        for e in self.inner {
            if last != e.prob {
                last = e.prob;
                inner.push(e);
            }
        }
        CDF { inner }
    }

    /// Downsample CDF to n entries. Panics if n <= 1 and returns identity if n is greater
    /// than the number of entries.
    ///
    /// Runtime complexity: O(m), where m is the original number of `Entry`s in `CDF`.
    ///
    /// # Arguments
    ///
    /// * `n` - Number of entries after downsampling.
    pub fn sample(mut self, n: usize) -> Self {
        assert!(n > 1);
        if self.inner.len() <= n {
            self
        } else {
            let s = self.inner.len() / (n - 1);
            let last = self.inner.pop().unwrap();
            let mut inner = self.inner.into_iter().step_by(s).collect_vec();
            inner.push(last);
            CDF { inner }
        }
    }

    /// Provide an iterator for the CDF.
    pub fn iter(&self) -> slice::Iter<'_, Entry<T>> {
        self.inner.iter()
    }

    /// Provide a mutable iterator over entries.
    ///
    /// This does not check for consistency. In other words, you
    /// should not change the order of the entries, nor the probabilities!
    pub fn iter_mut(&mut self) -> slice::IterMut<'_, Entry<T>> {
        self.inner.iter_mut()
    }

    /// Provide an iterator over the PMF corresponding to this CDF.
    pub fn iter_pmf(&self) -> CDFPMFIter<'_, T> {
        fn cdf_to_pmf<'a, G: Ord>(
            last_prob: &mut LogProb,
            e: &'a Entry<G>,
        ) -> Option<Entry<&'a G>> {
            let prob = e.prob.ln_sub_exp(*last_prob);
            *last_prob = e.prob;
            Some(Entry::new(&e.value, prob))
        }
        self.inner.iter().scan(LogProb::ln_zero(), cdf_to_pmf)
    }

    /// Get cumulative probability for a given value.
    ///
    /// If the value is not present, return the probability of the previous value.
    /// Time complexity: O(log n), where n is the number of `Entry`s in `CDF`.
    ///
    /// # Arguments
    ///
    /// * `value` - A value at which you're interested in the cumulative probability.
    pub fn get(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        } else {
            Some(match self.inner.binary_search_by(|e| e.value.cmp(value)) {
                Ok(i) => self.inner[i].prob,
                Err(i) => {
                    if i > 0 {
                        self.inner[i - 1].prob
                    } else {
                        LogProb::ln_zero()
                    }
                }
            })
        }
    }

    /// Get probability (i.e. probability mass) for a given `value`.
    ///
    /// Time complexity: O(log n), where n is the number of `Entry`s in `CDF`.
    pub fn get_pmf(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        } else {
            Some(match self.inner.binary_search_by(|e| e.value.cmp(value)) {
                Ok(i) => {
                    if i > 0 {
                        self.inner[i].prob.ln_sub_exp(self.inner[i - 1].prob)
                    } else {
                        self.inner[0].prob
                    }
                }
                Err(i) => {
                    if i > 0 {
                        self.inner[i - 1].prob
                    } else {
                        LogProb::ln_zero()
                    }
                }
            })
        }
    }

    /// Return total probability of the `CDF`.
    ///
    /// Time complexity: O(1).
    pub fn total_prob(&self) -> LogProb {
        self.inner.last().map_or(LogProb::ln_zero(), |e| e.prob)
    }

    /// Return maximum a posteriori probability estimate (MAP).
    ///
    /// Time complexity: O(n), where n is the number of `Entry`s in `CDF`.
    pub fn map(&self) -> Option<&T> {
        if let Some(mut max) = self.iter_pmf().next() {
            for e in self.iter_pmf() {
                if e.prob >= max.prob {
                    max = e;
                }
            }
            Some(max.value)
        } else {
            None
        }
    }

    /// Return w%-credible interval. The width w is a float between 0 and 1. Panics otherwise.
    /// E.g. provide `width=0.95` for the 95% credible interval.
    ///
    /// Runtime complexity: O(log n), where n is the number of `Entry`s in `CDF`.
    ///
    /// # Arguments
    ///
    /// * `width` - wanted width of the credible interval as a fraction of 1.
    pub fn credible_interval(&self, width: f64) -> Option<Range<&T>> {
        assert!(width >= 0.0 && width <= 1.0);

        if self.inner.is_empty() {
            return None;
        }

        let margin = 1.0 - width;
        let p_lower = OrderedFloat((margin / 2.0).ln());
        let p_upper = OrderedFloat((1.0 - margin / 2.0).ln());
        let lower = self
            .inner
            .binary_search_by(|e| OrderedFloat(*e.prob).cmp(&p_lower))
            .unwrap_or_else(|i| if i > 0 { i - 1 } else { 0 });
        let mut upper = self
            .inner
            .binary_search_by(|e| OrderedFloat(*e.prob).cmp(&p_upper))
            .unwrap_or_else(|i| i);
        if upper == self.inner.len() {
            upper -= 1;
        }

        Some(&self.inner[lower].value..&self.inner[upper].value)
    }

    /// Number of `Entry`s in the `CDF`.
    ///
    /// Time complexity: O(1)
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Returns `true` if a CDF is empty, false otherwise.
    ///
    /// Time complexity: O(1)
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

impl<T: Clone + Ord> CDF<T>
where
    f64: From<T>,
{
    /// Calculate expected value.
    ///
    /// Runtime complexity: O(n), where n is the number of `Entry`s in `CDF`.
    pub fn expected_value(&self) -> f64 {
        self.iter_pmf()
            .map(|e| f64::from(e.value.clone()) * e.prob.exp())
            .fold(0.0f64, |s, e| s + e)
    }

    /// Calculate variance.
    ///
    /// Runtime complexity: O(n), where n is the number of `Entry`s in `CDF`.
    pub fn variance(&self) -> f64 {
        let ev = self.expected_value();
        self.iter_pmf()
            .map(|e| (f64::from(e.value.clone()) - ev).powi(2) * e.prob.exp())
            .fold(0.0, |s, e| s + e)
    }

    /// Calculate standard deviation.
    ///
    /// Runtime complexity: O(n), where n is the number of `Entry`s in `CDF`.
    pub fn standard_deviation(&self) -> f64 {
        self.variance().sqrt()
    }
}

pub type CDFPMFIter<'a, T> = iter::Scan<
    slice::Iter<'a, Entry<T>>,
    LogProb,
    fn(&mut LogProb, &'a Entry<T>) -> Option<Entry<&'a T>>,
>;

#[cfg(test)]
mod test {
    use super::*;
    use crate::stats::LogProb;
    use ordered_float::NotNan;

    #[test]
    fn test_cdf() {
        let mut pmf = vec![Entry::new(NotNan::new(0.0).unwrap(), LogProb(0.1f64.ln()))];
        for i in 0..9 {
            pmf.push(Entry::new(
                NotNan::new(i as f64).unwrap(),
                LogProb(0.1f64.ln()),
            ));
        }
        println!("{:?}", pmf);

        let cdf = CDF::from_pmf(pmf.clone());
        println!("{:?}", cdf);
        for e in pmf.iter().skip(2) {
            assert_relative_eq!(*e.prob, *cdf.get_pmf(&e.value).unwrap(), epsilon = 0.000003);
        }
        assert_relative_eq!(*cdf.total_prob(), 1.0f64.ln());
        assert_relative_eq!(
            *cdf.get(&NotNan::new(1.0).unwrap()).unwrap(),
            0.3f64.ln(),
            epsilon = 0.00000001
        );
        {
            let ci = cdf.credible_interval(0.95).unwrap();
            assert_relative_eq!(**ci.start, 0.0);
            assert_relative_eq!(**ci.end, 8.0);
        }

        {
            for e in cdf.iter_pmf() {
                assert_relative_eq!(
                    e.prob.exp(),
                    if **e.value == 0.0 { 0.2 } else { 0.1 },
                    epsilon = 0.0001
                );
            }
        }

        assert_relative_eq!(cdf.sample(5).total_prob().exp(), 1.0);
    }
}
