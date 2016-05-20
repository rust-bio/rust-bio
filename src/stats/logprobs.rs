// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Handling log-probabilities.

use std::mem;
use std::f64;
use std::iter;

pub use stats::{Prob, LogProb};


/// A factor to convert log-probabilities to PHRED-scale (phred = p * LOG_TO_PHRED_FACTOR).
const LOG_TO_PHRED_FACTOR: f64 = -4.3429448190325175; // -10 * 1 / ln(10)


/// A factor to convert PHRED-scale to log-probabilities (p = phred * PHRED_TO_LOG_FACTOR).
const PHRED_TO_LOG_FACTOR: f64 = -0.23025850929940456; // 1 / (-10 * log10(e))


/// Calculate log(1 - p) with p given in log space without loss of precision as described in
/// http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf.
pub fn ln_1m_exp(p: LogProb) -> LogProb {
    assert!(p <= 0.0);
    if p < -0.693 {
        (-p.exp()).ln_1p()
    } else {
        (-p.exp_m1()).ln()
    }
}


/// Convert log scale probability to PHRED scale.
pub fn log_to_phred(p: LogProb) -> f64 {
    p * LOG_TO_PHRED_FACTOR
}


/// Convert PHRED scale probability to log scale.
pub fn phred_to_log(p: f64) -> LogProb {
    p * PHRED_TO_LOG_FACTOR
}


/// Calculate the sum of the given probabilities in a numerically stable way (Durbin 1998).
pub fn sum(probs: &[LogProb]) -> LogProb {
    if probs.is_empty() {
        f64::NEG_INFINITY
    } else {
        let mut pmax = probs[0];
        let mut imax = 0;
        for (i, &p) in probs.iter().enumerate().skip(1) {
            if p > pmax {
                pmax = p;
                imax = i;
            }
        }
        if pmax == f64::NEG_INFINITY {
            f64::NEG_INFINITY
        } else if pmax == f64::INFINITY {
            f64::INFINITY
        } else {
            // TODO use sum() once it has been stabilized: .sum::<usize>()
            pmax +
            (probs.iter()
                  .enumerate()
                  .filter_map(|(i, p)| {
                      if i == imax {
                          None
                      } else {
                          Some((p - pmax).exp())
                      }
                  })
                  .fold(0.0, |s, e| s + e))
                .ln_1p()
        }
    }
}


/// Add the given probabilities in a numerically stable way (Durbin 1998).
pub fn add(mut p0: LogProb, mut p1: LogProb) -> LogProb {
    if p1 > p0 {
        mem::swap(&mut p0, &mut p1);
    }
    if p0 == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    } else if p0 == f64::INFINITY {
        f64::INFINITY
    } else {
        p0 + (p1 - p0).exp().ln_1p()
    }
}


/// Subtract the given probabilities in a numerically stable way (Durbin 1998).
pub fn sub(p0: LogProb, p1: LogProb) -> LogProb {
    assert!(p0 >= p1, "Subtraction would lead to negative probability, which is undefined in log space.");
    if p0 == p1 || p0 == f64::NEG_INFINITY {
        // the first case leads to zero,
        // in the second case p0 and p1 are -inf, which is fine
        f64::NEG_INFINITY
    } else if p0 == f64::INFINITY {
        f64::INFINITY
    } else {
        p0 + ln_1m_exp(p1 - p0)
    }
}



fn scan_add(s: &mut LogProb, p: LogProb) -> Option<LogProb> {
    *s = add(*s, p);
    Some(*s)
}


/// Iterator returned by scans over logprobs.
pub type ScanIter<I> = iter::Scan<I, LogProb, fn(&mut LogProb, LogProb) -> Option<LogProb>>;


/// Calculate the cumulative sum of the given probabilities in a numerically stable way (Durbin 1998).
pub fn cumsum<I: Iterator<Item = LogProb>>(probs: I) -> ScanIter<I> {
    probs.scan(f64::NEG_INFINITY, scan_add)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64;
    use itertools::Itertools;

    #[test]
    fn test_sum() {
        let probs = [f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY];
        assert_eq!(sum(&probs), 0.0);
    }

    #[test]
    fn test_empty_sum() {
        assert_eq!(sum(&[]), f64::NEG_INFINITY);
    }

    #[test]
    fn test_cumsum() {
        let probs = vec![0.0f64.ln(), 0.01f64.ln(), 0.001f64.ln()];
        assert_eq!(cumsum(probs.into_iter()).collect_vec(),
                   [0.0f64.ln(), 0.01f64.ln(), 0.011f64.ln()]);
    }

    #[test]
    fn test_sub() {
        assert_eq!(sub(0.0f64.ln(), 0.0f64.ln()), f64::NEG_INFINITY);
        assert_relative_eq!(sub(1.0f64.ln(), 0.5f64.ln()), 0.5f64.ln());
    }

    #[test]
    fn test_ln_1m_exp() {
        assert_eq!(ln_1m_exp(f64::NEG_INFINITY), 0.0);
        assert_eq!(ln_1m_exp(-0.0), f64::NEG_INFINITY);
    }
}
