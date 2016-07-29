// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f64;
use std::iter;
use std::slice;

use num::traits::{cast, NumCast};
use itertools::Itertools;

use stats::logprobs::{self, LogProb};


/// Implementation of a cumulative distribution function.
#[derive(Debug, Clone)]
pub struct CDF<T: PartialOrd> {
    inner: Vec<(T, LogProb)>
}


impl<T: PartialOrd> CDF<T> {
    /// Create CDF from given probability mass function (PMF). The PMF may contain duplicate values
    /// the probabilities of which are summed during generation of the CDF.
    pub fn from_pmf(mut entries: Vec<(T, LogProb)>) -> Self {
        entries.sort_by(|&(ref a, _), &(ref b, _)| a.partial_cmp(b).unwrap());
        let mut inner: Vec<(T, LogProb)> = Vec::new();
        for mut e in entries.into_iter() {
            let p = logprobs::add(inner.last().map_or(f64::NEG_INFINITY, |e| e.1), e.1);
            if !inner.is_empty() && inner.last().unwrap().0 == e.0 {
                inner.last_mut().unwrap().1 = p;
            }
            else {
                e.1 = p;
                inner.push(e);
            }
        }
        let mut cdf = CDF {
            inner: inner
        };

        if relative_eq!(cdf.total_prob(), 0.0) && cdf.total_prob() > 0.0 {
            cdf.inner.last_mut().unwrap().1 = 0.0;
        }

        cdf
    }

    /// Create CDF from iterator. This can be used to replace the values of a CDF.
    pub fn from_cdf<I: Iterator<Item = (T, LogProb)>>(entries: I) -> Self {
        CDF { inner: entries.collect_vec() }
    }

    pub fn reduce(self) -> Self {
        let mut inner = Vec::new();
        let mut last = f64::NEG_INFINITY;
        for e in self.inner.into_iter() {
            if last != e.1 {
                last = e.1;
                inner.push(e);
            }
        }
        CDF { inner: inner }
    }

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

    pub fn iter(&self) -> slice::Iter<(T, f64)>{
        self.inner.iter()
    }

    pub fn iter_pmf<'a>(&'a self) -> CDFPMFIter<'a, T> {
        fn cdf_to_pmf<'a, G>(last_prob: &mut LogProb, e: &'a (G, LogProb)) -> Option<(&'a G, LogProb)> {
            let &(ref value, cdf_prob) = e;
            let prob = logprobs::sub(cdf_prob, *last_prob);
            *last_prob = cdf_prob;
            Some((value, prob))
        }
        self.inner.iter().scan(f64::NEG_INFINITY, cdf_to_pmf)
    }

    /// Get cumulative probability for a given value. If the value is not present,
    /// return the probability of the previous value.
    pub fn get(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        }
        else {
            Some(match self.inner.binary_search_by(|e| e.0.partial_cmp(value).unwrap()) {
                Ok(i) => self.inner[i].1,
                Err(i) => if i > 0 { self.inner[i - 1].1 } else { f64::NEG_INFINITY }
            })
        }
    }

    pub fn get_pmf(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        }
        else {
            Some(match self.inner.binary_search_by(|e| e.0.partial_cmp(value).unwrap()) {
                Ok(i) => if i > 0 { logprobs::sub(self.inner[i].1, self.inner[i - 1].1) } else { self.inner[0].1 },
                Err(i) => if i > 0 { self.inner[i - 1].1 } else { f64::NEG_INFINITY }
            })
        }
    }

    pub fn total_prob(&self) -> LogProb {
        let &(_, prob) = self.inner.last().unwrap();
        prob
    }

    /// Return maximum a posteriori probability estimate (MAP).
    pub fn map(&self) -> &T {
        let mut max = self.iter_pmf().next().unwrap();
        for e in self.iter_pmf() {
            if e.1 >= max.1 {
                max = e;
            }
        }
        &max.0
    }

    pub fn credible_interval(&self) -> (&T, &T) {
        let lower = self.inner.binary_search_by(|&(_, p)| p.partial_cmp(&0.025f64.ln()).unwrap()).unwrap_or_else(|i| i);
        let upper = self.inner.binary_search_by(|&(_, p)| p.partial_cmp(&0.975f64.ln()).unwrap()).unwrap_or_else(|i| i - 1);

        (&self.inner[lower].0, &self.inner[upper].0)
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }
}


impl<T: NumCast + Clone + PartialOrd> CDF<T> {
    pub fn expected_value(&self) -> f64 {
        self.iter_pmf().map(|(value, prob)| {
            cast::<T, f64>(value.clone()).unwrap() * prob.exp()
        }).fold(0.0f64, |s, e| s + e)
    }

    pub fn variance(&self) -> f64 {
        let ev = self.expected_value();
        self.iter_pmf().map(|(value, prob)| {
                (cast::<T, f64>(value.clone()).unwrap() - ev).powi(2) * prob.exp()
        }).fold(0.0, |s, e| s + e)
    }

    pub fn standard_deviation(&self) -> f64 {
        self.variance().sqrt()
    }
}

pub type CDFPMFIter<'a, T> = iter::Scan<slice::Iter<'a, (T, LogProb)>, LogProb, fn(&mut LogProb, &'a (T, LogProb)) -> Option<(&'a T, LogProb)>>;


#[cfg(test)]
mod test {
    use super::*;

    use itertools::Itertools;


    #[test]
    fn test_cdf() {
        let mut pmf = vec![(0.0, 0.1f64.ln())];
        for i in 0..9 {
            pmf.push((i as f64, 0.1f64.ln()));
        }
        println!("{:?}", pmf);

        let cdf = CDF::from_pmf(pmf.clone());
        println!("{:?}", cdf);
        for &(value, prob) in pmf.iter().skip(2) {
            assert_ulps_eq!(prob, cdf.get_pmf(&value).unwrap(), epsilon = 0.0000000000001);
        }
        assert_relative_eq!(cdf.total_prob(), 1.0f64.ln());
        assert_relative_eq!(cdf.get(&1.0).unwrap(), 0.3f64.ln(), epsilon = 0.00000001);
    }
}
