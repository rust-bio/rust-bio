// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Support for discrete probability distributions in terms of cumulative distribution
//! functions (CDF).

use std::f64;
use std::iter;
use std::slice;
use std::ops::Range;

use num::traits::{cast, NumCast};
use itertools::Itertools;

use stats::LogProb;


#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde_macros", derive(Serialize, Deserialize))]
pub struct Entry<T: Ord> {
    pub value: T,
    pub prob: LogProb
}


impl<T: Ord> Entry<T> {
    pub fn new(value: T, prob: LogProb) -> Self {
        Entry { value: value, prob: prob }
    }
}



/// Implementation of a cumulative distribution function.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde_macros", derive(Serialize, Deserialize))]
pub struct CDF<T: Ord> {
    inner: Vec<Entry<T>>
}


impl<T: Ord> CDF<T> {
    /// Create CDF from given probability mass function (PMF). The PMF may contain duplicate values
    /// the probabilities of which are summed during generation of the CDF.
    ///
    /// # Arguments
    ///
    /// * `pmf` - the PMF as a vector of `Entry` objects
    pub fn from_pmf(mut entries: Vec<Entry<T>>) -> Self {
        entries.sort_by(|a, b| a.value.cmp(&b.value));
        let mut inner: Vec<Entry<T>> = Vec::new();
        for mut e in entries.into_iter() {
            let p = inner.last().map_or(LogProb::ln_zero(), |e| e.prob).ln_add_exp(e.prob);
            if !inner.is_empty() && inner.last().unwrap().value == e.value {
                inner.last_mut().unwrap().prob = p;
            }
            else {
                e.prob = p;
                inner.push(e);
            }
        }
        let mut cdf = CDF {
            inner: inner
        };

        if relative_eq!(*cdf.total_prob(), *LogProb::ln_one()) && *cdf.total_prob() > *LogProb::ln_one() {
            cdf.inner.last_mut().unwrap().prob = LogProb::ln_one();
        }

        cdf
    }

    /// Create CDF from iterator. This can be used to replace the values of a CDF.
    pub fn from_cdf<I: Iterator<Item = Entry<T>>>(entries: I) -> Self {
        CDF { inner: entries.collect_vec() }
    }

    /// Reduce CDF by omitting values with zero probability.
    pub fn reduce(self) -> Self {
        let mut inner = Vec::new();
        let mut last = LogProb::ln_zero();
        for e in self.inner.into_iter() {
            if last != e.prob {
                last = e.prob;
                inner.push(e);
            }
        }
        CDF { inner: inner }
    }

    /// Downsample CDF to n entries.
    pub fn sample(mut self, n: usize) -> Self {
        if self.inner.len() <= n {
            self
        }
        else {
            let s = self.inner.len() / (n - 1);
            let last = self.inner.pop().unwrap();
            let mut inner = self.inner.into_iter().step(s).collect_vec();
            inner.push(last);
            CDF {
                inner: inner
            }
        }
    }

    /// Provide iterator.
    pub fn iter(&self) -> slice::Iter<Entry<T>>{
        self.inner.iter()
    }

    /// Iterator over corresponding PMF.
    pub fn iter_pmf<'a>(&'a self) -> CDFPMFIter<'a, T> {
        fn cdf_to_pmf<'a, G: Ord>(last_prob: &mut LogProb, e: &'a Entry<G>) -> Option<Entry<&'a G>> {
            let prob = e.prob.ln_sub_exp(*last_prob);
            *last_prob = e.prob;
            Some(Entry::new(&e.value, prob))
        }
        self.inner.iter().scan(LogProb::ln_zero(), cdf_to_pmf)
    }

    /// Get cumulative probability for a given value. If the value is not present,
    /// return the probability of the previous value. Complexity O(log n).
    pub fn get(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        }
        else {
            Some(match self.inner.binary_search_by(|e| e.value.partial_cmp(value).unwrap()) {
                Ok(i) => self.inner[i].prob,
                Err(i) => if i > 0 { self.inner[i - 1].prob } else { LogProb::ln_zero() }
            })
        }
    }

    /// Get probability (i.e. probability mass) for a given value. Complexity O(log n).
    pub fn get_pmf(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        }
        else {
            Some(match self.inner.binary_search_by(|e| e.value.partial_cmp(value).unwrap()) {
                Ok(i) => if i > 0 { self.inner[i].prob.ln_sub_exp(self.inner[i - 1].prob) } else { self.inner[0].prob },
                Err(i) => if i > 0 { self.inner[i - 1].prob } else { LogProb::ln_zero() }
            })
        }
    }

    /// Return total probability.
    pub fn total_prob(&self) -> LogProb {
        self.inner.last().map_or(LogProb::ln_zero(), |e| e.prob)
    }

    /// Return maximum a posteriori probability estimate (MAP).
    pub fn map(&self) -> Option<&T> {
        if let Some(mut max) = self.iter_pmf().next() {
            for e in self.iter_pmf() {
                if e.prob >= max.prob {
                    max = e;
                }
            }
            Some(&max.value)
        } else {
            None
        }
    }

    /// Return 95% credible interval.
    pub fn credible_interval(&self) -> Range<&T> {
        let p_lower = LogProb(0.025f64.ln());
        let p_upper = LogProb(0.095f64.ln());
        let lower = self.inner.binary_search_by(|e| e.prob.partial_cmp(&p_lower).unwrap()).unwrap_or_else(|i| i);
        let upper = self.inner.binary_search_by(|e| e.prob.partial_cmp(&p_upper).unwrap()).unwrap_or_else(|i| i - 1);

        &self.inner[lower].value..&self.inner[upper].value
    }

    /// Number of entries in the CDF.
    pub fn len(&self) -> usize {
        self.inner.len()
    }
}


impl<T: NumCast + Clone + Ord> CDF<T> {
    /// Calculate expected value.
    pub fn expected_value(&self) -> f64 {
        self.iter_pmf().map(|e| {
            cast::<T, f64>(e.value.clone()).unwrap() * e.prob.exp()
        }).fold(0.0f64, |s, e| s + e)
    }

    /// Calculate variance.
    pub fn variance(&self) -> f64 {
        let ev = self.expected_value();
        self.iter_pmf().map(|e| {
                (cast::<T, f64>(e.value.clone()).unwrap() - ev).powi(2) * e.prob.exp()
        }).fold(0.0, |s, e| s + e)
    }

    /// Calculate standard deviation.
    pub fn standard_deviation(&self) -> f64 {
        self.variance().sqrt()
    }
}

pub type CDFPMFIter<'a, T> = iter::Scan<slice::Iter<'a, Entry<T>>, LogProb, fn(&mut LogProb, &'a Entry<T>) -> Option<Entry<&'a T>>>;


#[cfg(test)]
mod test {
    use super::*;
    use stats::LogProb;
    use ordered_float::NotNaN;


    #[test]
    fn test_cdf() {
        let mut pmf = vec![Entry::new(NotNaN::new(0.0).unwrap(), LogProb(0.1f64.ln()))];
        for i in 0..9 {
            pmf.push(Entry::new(NotNaN::new(i as f64).unwrap(), LogProb(0.1f64.ln())));
        }
        println!("{:?}", pmf);

        let cdf = CDF::from_pmf(pmf.clone());
        println!("{:?}", cdf);
        for e in pmf.iter().skip(2) {
            assert_ulps_eq!(*e.prob, *cdf.get_pmf(&e.value).unwrap(), epsilon = 0.0000000000001);
        }
        assert_relative_eq!(*cdf.total_prob(), 1.0f64.ln());
        assert_relative_eq!(*cdf.get(&NotNaN::new(1.0).unwrap()).unwrap(), 0.3f64.ln(), epsilon = 0.00000001);
    }
}
