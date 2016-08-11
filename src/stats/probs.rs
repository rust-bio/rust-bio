// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Handling log-probabilities.

use std::mem;
use std::f64;
use std::iter;
use std::ops::{Add, Sub, Mul, Div};
use num::traits::cast;
use num::NumCast;

use itertools::linspace;
use itertools::Itertools;
use itertools::misc::ToFloat;


/// A factor to convert log-probabilities to PHRED-scale (phred = p * `LOG_TO_PHRED_FACTOR`).
const LOG_TO_PHRED_FACTOR: f64 = -4.3429448190325175; // -10 * 1 / ln(10)


/// A factor to convert PHRED-scale to log-probabilities (p = phred * `PHRED_TO_LOG_FACTOR`).
const PHRED_TO_LOG_FACTOR: f64 = -0.23025850929940456; // 1 / (-10 * log10(e))


/// Calculate log(1 - p) with p given in log space without loss of precision as described in
/// http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf.
fn ln_1m_exp(p: f64) -> f64 {
    assert!(p <= 0.0);
    if p < -0.693 {
        (-p.exp()).ln_1p()
    } else {
        (-p.exp_m1()).ln()
    }
}


/// A newtype for probabilities.
custom_derive! {
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        NewtypeAdd(*),
        NewtypeSub(*),
        PartialEq,
        PartialOrd,
        Copy,
        Clone,
        Debug
    )]
    pub struct Prob(f64);
}


/// A newtype for log-scale probabilities.
custom_derive! {
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        NewtypeAdd(*),
        NewtypeSub(*),
        PartialEq,
        PartialOrd,
        Copy,
        Clone,
        Debug
    )]
    pub struct LogProb(f64);
}


/// A newtype for PHRED-scale probabilities.
custom_derive! {
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        NewtypeAdd(*),
        NewtypeSub(*),
        PartialEq,
        PartialOrd,
        Copy,
        Clone,
        Debug
    )]
    pub struct PHREDProb(f64);
}


/// Iterator returned by scans over logprobs.
pub type ScanIter<I> = iter::Scan<<I as IntoIterator>::IntoIter, LogProb, fn(&mut LogProb, LogProb) -> Option<LogProb>>;


impl LogProb {
    pub fn ln_zero() -> LogProb {
        LogProb(f64::NEG_INFINITY)
    }

    pub fn ln_one() -> LogProb {
        LogProb(0.0)
    }

    pub fn ln_one_minus_exp(&self) -> LogProb {
        LogProb(ln_1m_exp(**self))
    }

    pub fn ln_sum_exp(probs: &[LogProb]) -> LogProb {
        if probs.is_empty() {
            Self::ln_zero()
        } else {
            let mut pmax = probs[0];
            let mut imax = 0;
            for (i, &p) in probs.iter().enumerate().skip(1) {
                if p > pmax {
                    pmax = p;
                    imax = i;
                }
            }
            if pmax == Self::ln_zero() {
                Self::ln_zero()
            } else if *pmax == f64::INFINITY {
                LogProb(f64::INFINITY)
            } else {
                // TODO use sum() once it has been stabilized: .sum::<usize>()
                pmax + LogProb(
                    (probs.iter()
                          .enumerate()
                          .filter_map(|(i, p)| {
                              if i == imax {
                                  None
                              } else {
                                  Some((p - pmax).exp())
                              }
                          })
                          .fold(0.0, |s, e| s + e)
                     ).ln_1p()
                 )
            }
        }
    }

    pub fn ln_add_exp(self, other: LogProb) -> LogProb {
        let (mut p0, mut p1) = (self, other);
        if p1 > p0 {
            mem::swap(&mut p0, &mut p1);
        }
        if p0 == Self::ln_zero() {
            Self::ln_zero()
        } else if *p0 == f64::INFINITY {
            LogProb(f64::INFINITY)
        } else {
            p0 + LogProb((p1 - p0).exp().ln_1p())
        }
    }

    pub fn ln_sub_exp(self, other: LogProb) -> LogProb {
        let (p0, p1) = (self, other);
        assert!(p0 >= p1,
                "Subtraction would lead to negative probability, which is undefined in log space.");
        if relative_eq!(*p0, *p1) || p0 == Self::ln_zero() {
            // the first case leads to zero,
            // in the second case p0 and p1 are -inf, which is fine
            Self::ln_zero()
        } else if *p0 == f64::INFINITY {
            LogProb(f64::INFINITY)
        } else {
            p0 + (p1 - p0).ln_one_minus_exp()
        }
    }

    /// Calculate the cumulative sum of the given probabilities in a numerically stable way (Durbin 1998).
    pub fn ln_cumsum_exp<I: IntoIterator<Item = LogProb>>(probs: I) -> ScanIter<I> {
        probs.into_iter().scan(Self::ln_zero(), Self::scan_ln_add_exp)
    }

    /// Integrate numerically stable over given log-space density in the interval [a, b]. Uses the trapezoidal rule with n grid points.
    pub fn ln_integrate_exp<T, D>(density: &D, a: T, b: T, n: usize) -> LogProb
        where T: NumCast + Copy + Add<Output=T> + Sub<Output=T> + Div<Output=T> + Mul<Output=T>, D: Fn(T) -> LogProb, usize: ToFloat<T>
    {
        let mut probs = linspace(a, b, n).dropping(1).dropping_back(1).map(|v| LogProb(*density(v) + 2.0f64.ln())).collect_vec();
        probs.push(density(a));
        probs.push(density(b));
        let width: f64 = cast(b - a).unwrap();

        LogProb(*Self::ln_sum_exp(&probs) + width.ln() - (2.0 * (n - 1) as f64).ln())
    }

    fn scan_ln_add_exp(s: &mut LogProb, p: LogProb) -> Option<LogProb> {
        *s = s.ln_add_exp(p);
        Some(*s)
    }
}


impl From<LogProb> for Prob {
    fn from(p: LogProb) -> Prob {
        Prob(p.exp())
    }
}


impl From<PHREDProb> for Prob {
    fn from(p: PHREDProb) -> Prob {
        Prob(10.0f64.powf(-*p / 10.0))
    }
}


impl From<Prob> for LogProb {
    fn from(p: Prob) -> LogProb {
        LogProb(p.ln())
    }
}


impl From<PHREDProb> for LogProb {
    fn from(p: PHREDProb) -> LogProb {
        LogProb(*p * PHRED_TO_LOG_FACTOR)
    }
}


impl From<LogProb> for PHREDProb {
    fn from(p: LogProb) -> PHREDProb {
        PHREDProb(*p * LOG_TO_PHRED_FACTOR)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64;
    use itertools::Itertools;

    #[test]
    fn test_sum() {
        let probs = [LogProb::ln_zero(), LogProb::ln_one(), LogProb::ln_zero()];
        assert_eq!(LogProb::ln_sum_exp(&probs), LogProb::ln_one());
    }

    #[test]
    fn test_empty_sum() {
        assert_eq!(LogProb::ln_sum_exp(&[]), LogProb::ln_zero());
    }

    #[test]
    fn test_cumsum() {
        let probs = vec![LogProb::ln_zero(), LogProb(0.01f64.ln()), LogProb(0.001f64.ln())];
        assert_eq!(LogProb::ln_cumsum_exp(probs).collect_vec(),
                   [LogProb::ln_zero(), LogProb(0.01f64.ln()), LogProb(0.011f64.ln())]);
    }

    #[test]
    fn test_sub() {
        assert_eq!(LogProb::ln_one().ln_sub_exp(LogProb::ln_one()), LogProb::ln_zero());
        assert_relative_eq!(*LogProb::ln_one().ln_sub_exp(LogProb(0.5f64.ln())), *LogProb(0.5f64.ln()));
    }

    #[test]
    fn test_one_minus() {
        assert_eq!(LogProb::ln_zero().ln_one_minus_exp(), LogProb::ln_one());
        assert_eq!(LogProb::ln_one().ln_one_minus_exp(), LogProb::ln_zero());
    }

    #[test]
    fn test_integrate() {
        let density = |_| LogProb(0.1f64.ln());
        let prob = LogProb::ln_integrate_exp(&density, 0.0, 10.0, 5);
        assert_relative_eq!(*prob, *LogProb::ln_one(), epsilon=0.0000001);
    }
}
