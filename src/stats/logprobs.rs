// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Handling log-probabilities.

use std::mem;
use std::f64;


/// An alias for `f64` to indicate a probability.
pub type Prob = f64;


/// An alias for `f64` to indicate a log-probability.
pub type LogProb = f64;


/// A factor to convert log-probabilities to PHRED-scale (phred = p * LOG_TO_PHRED_FACTOR).
const LOG_TO_PHRED_FACTOR: f64 = -4.3429448190325175; // -10 * 1 / ln(10)


/// A factor to convert PHRED-scale to log-probabilities (p = phred * PHRED_TO_LOG_FACTOR).
const PHRED_TO_LOG_FACTOR: f64 = -0.23025850929940456; // 1 / (-10 * log10(e))


/// Calculate log(1 - p) with p given in log space without loss of precision as described in
/// http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf.
pub fn ln_1m_exp(p: LogProb) -> LogProb {
    if p < -0.693 {
        (-p.exp()).ln_1p()
    }
    else {
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
pub fn log_prob_sum(probs: &[LogProb]) -> LogProb {
    if probs.is_empty() {
        f64::NEG_INFINITY
    }
    else {
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
        }
        else {
            pmax + (probs.iter().enumerate().filter_map(|(i, p)| if i != imax { Some((p - pmax).exp()) } else { None }).sum::<Prob>()).ln_1p()
        }
    }
}


/// Calcualte the sum of the given probabilities in a numerically stable way (Durbin 1998).
pub fn log_prob_add(mut p0: LogProb, mut p1: LogProb) -> LogProb {
    if p1 > p0 {
        mem::swap(&mut p0, &mut p1);
    }
    if p0 == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    }
    else {
        p0 + (p1 - p0).exp().ln_1p()
    }
}


/// Calcualte the cumulative sum of the given probabilities in a numerically stable way (Durbin 1998).
pub fn log_prob_cumsum(probs: &[LogProb]) -> Vec<LogProb> {
    probs.iter().scan(f64::NEG_INFINITY, |s, p| {
        *s = log_prob_add(*s, *p);
        Some(*s)
    }).collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64;

    #[test]
    fn test_log_prob_sum() {
        let probs = [f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY];
        assert_eq!(log_prob_sum(&probs), 0.0);
    }

    #[test]
    fn test_empty_log_prob_sum() {
        assert_eq!(log_prob_sum(&[]), f64::NEG_INFINITY);
    }
}
