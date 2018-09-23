// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Handling log-probabilities.

pub mod cdf;

use std::f64;
use std::iter;
use std::mem;
use std::ops::{Add, AddAssign, Div, Mul, Sub, SubAssign};

use itertools::Itertools;
use itertools_num::linspace;
use num_traits::{Float, Zero};
use ordered_float::NotNaN;

/// A factor to convert log-probabilities to PHRED-scale (phred = p * `LOG_TO_PHRED_FACTOR`).
const LOG_TO_PHRED_FACTOR: f64 = -4.342_944_819_032_517_5; // -10 * 1 / ln(10)

/// A factor to convert PHRED-scale to log-probabilities (p = phred * `PHRED_TO_LOG_FACTOR`).
const PHRED_TO_LOG_FACTOR: f64 = -0.230_258_509_299_404_56; // 1 / (-10 * log10(e))

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

custom_derive! {
    /// A newtype for probabilities.
    ///
    /// # Example
    ///
    /// ```
    /// #[macro_use]
    /// extern crate approx;
    /// # extern crate bio;
    /// # fn main() {
    /// use bio::stats::Prob;
    ///
    /// let p = Prob(0.5);
    /// let q = Prob(0.2);
    ///
    /// assert_relative_eq!(*(p + q), *Prob(0.7));
    /// # }
    /// ```
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        NewtypeAdd(*),
        NewtypeSub(*),
        NewtypeMul(*),
        NewtypeDiv(*),
        PartialEq,
        PartialOrd,
        Copy,
        Clone,
        Debug,
        Default
    )]
    #[derive(Serialize, Deserialize)]
    pub struct Prob(pub f64);
}

impl Prob {
    pub fn checked(p: f64) -> Result<Self, ProbError> {
        if p >= 0.0 && p <= 1.0 {
            Ok(Prob(p))
        } else {
            Err(ProbError::InvalidProb(p))
        }
    }
}

custom_derive! {
    /// A newtype for log-scale probabilities.
    ///
    /// # Example
    ///
    /// ```
    /// #[macro_use]
    /// extern crate approx;
    /// # extern crate bio;
    /// # fn main() {
    /// use bio::stats::{LogProb, Prob};
    ///
    /// // convert from probability
    /// let p = LogProb::from(Prob(0.5));
    /// // convert manually
    /// let q = LogProb(0.2f64.ln());
    /// // obtain zero probability in log-space
    /// let o = LogProb::ln_one();
    ///
    /// assert_relative_eq!(*Prob::from(p.ln_add_exp(q) + o), *Prob(0.7));
    /// # }
    /// ```
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
    #[derive(Serialize, Deserialize)]
    pub struct LogProb(pub f64);
}

custom_derive! {
    /// A newtype for PHRED-scale probabilities.
    ///
    /// # Example
    ///
    /// ```
    /// #[macro_use]
    /// extern crate approx;
    /// # extern crate bio;
    /// # fn main() {
    /// use bio::stats::{PHREDProb, Prob};
    ///
    /// let p = PHREDProb::from(Prob(0.5));
    ///
    /// assert_relative_eq!(*Prob::from(p), *Prob(0.5));
    /// # }
    /// ```
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
    #[derive(Serialize, Deserialize)]
    pub struct PHREDProb(pub f64);
}

/// Iterator returned by scans over logprobs.
pub type ScanIter<I> = iter::Scan<
    <I as IntoIterator>::IntoIter,
    LogProb,
    fn(&mut LogProb, LogProb) -> Option<LogProb>,
>;

static LOGPROB_LN_ZERO: LogProb = LogProb(f64::NEG_INFINITY);
static LOGPROB_LN_ONE: LogProb = LogProb(0.0);

impl LogProb {
    pub fn is_valid(&self) -> bool {
        !self.is_nan() && **self <= 0.0
    }

    /// Log-space representation of Pr=0
    pub fn ln_zero() -> LogProb {
        LOGPROB_LN_ZERO
    }

    /// Log-space representation of Pr=1
    pub fn ln_one() -> LogProb {
        LOGPROB_LN_ONE
    }

    /// sums of LogProbs, e.g. with `ln_sum_exp()` can end up
    /// slightly above the maximum of `LogProb <= 0` due to
    /// numerical imprecisions -- this function can rescue such
    /// values before panics due to asserts in other functions
    /// handling LogProbs, e.g. `ln_1m_exp`
    pub fn cap_numerical_overshoot(&self, epsilon: f64) -> LogProb {
        if *self <= LogProb::ln_one() {
            *self
        } else {
            let capped = **self - epsilon;
            if capped <= 0.0 {
                LogProb::ln_one()
            } else {
                panic!(
                    "Cannot correct LogProb {:?} -- not within given epsilon of 0.0 ({})",
                    **self, epsilon
                );
            }
        }
    }

    /// Numerically stable calculation of 1 - p in log-space.
    pub fn ln_one_minus_exp(&self) -> LogProb {
        LogProb(ln_1m_exp(**self))
    }

    /// Numerically stable sum of probabilities in log-space.
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
                    (probs
                        .iter()
                        .enumerate()
                        .filter_map(|(i, p)| {
                            if i == imax {
                                None
                            } else {
                                Some((p - pmax).exp())
                            }
                        }).fold(0.0, |s, e| s + e)).ln_1p(),
                )
            }
        }
    }

    /// Numerically stable addition probabilities in log-space.
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

    /// Numerically stable subtraction of probabilities in log-space.
    pub fn ln_sub_exp(self, other: LogProb) -> LogProb {
        let (p0, p1) = (self, other);
        assert!(
            p0 >= p1,
            "Subtraction would lead to negative probability, which is undefined in log space."
        );
        if *p1 == f64::NEG_INFINITY {
            p0
        } else if relative_eq!(*p0, *p1) || p0 == Self::ln_zero() {
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
        probs
            .into_iter()
            .scan(Self::ln_zero(), Self::scan_ln_add_exp)
    }

    /// Integrate numerically stable over given log-space density in the interval [a, b]. Uses the trapezoidal rule with n grid points.
    pub fn ln_trapezoidal_integrate_exp<T, D>(mut density: D, a: T, b: T, n: usize) -> LogProb
    where
        T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Mul<Output = T> + Float,
        D: FnMut(T) -> LogProb,
        f64: From<T>,
    {
        let mut probs = linspace(a, b, n)
            .dropping(1)
            .dropping_back(1)
            .map(|v| LogProb(*density(v) + 2.0f64.ln()))
            .collect_vec();
        probs.push(density(a));
        probs.push(density(b));
        let width = f64::from(b - a);

        LogProb(*Self::ln_sum_exp(&probs) + width.ln() - (2.0 * (n - 1) as f64).ln())
    }

    /// Integrate numerically stable over given log-space density in the interval [a, b]. Uses Simpson's rule with n (odd) grid points.
    pub fn ln_simpsons_integrate_exp<T, D>(mut density: D, a: T, b: T, n: usize) -> LogProb
    where
        T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Mul<Output = T> + Float,
        D: FnMut(T) -> LogProb,
        f64: From<T>,
    {
        assert_eq!(n % 2, 1, "n must be odd");
        let mut probs = linspace(a, b, n)
            .enumerate()
            .dropping(1)
            .dropping_back(1)
            .map(|(i, v)| {
                let weight = (2 + (i % 2) * 2) as f64;
                LogProb(*density(v) + weight.ln()) // factors alter between 2 and 4
            }).collect_vec();
        probs.push(density(a));
        probs.push(density(b));
        let width = f64::from(b - a);

        LogProb(*Self::ln_sum_exp(&probs) + width.ln() - ((n - 1) as f64).ln() - 3.0f64.ln())
    }

    fn scan_ln_add_exp(s: &mut LogProb, p: LogProb) -> Option<LogProb> {
        *s = s.ln_add_exp(p);
        Some(*s)
    }
}

impl<'a> iter::Sum<&'a LogProb> for LogProb {
    fn sum<I: Iterator<Item = &'a LogProb>>(iter: I) -> Self {
        iter.fold(LogProb(0.0), |a, b| a + *b)
    }
}

impl<'a> iter::Sum<LogProb> for LogProb {
    fn sum<I: Iterator<Item = LogProb>>(iter: I) -> Self {
        iter.fold(LogProb(0.0), |a, b| a + b)
    }
}

impl AddAssign for LogProb {
    fn add_assign(&mut self, other: LogProb) {
        *self = *self + other;
    }
}

impl SubAssign for LogProb {
    fn sub_assign(&mut self, other: LogProb) {
        *self = *self - other;
    }
}

impl From<NotNaN<f64>> for LogProb {
    fn from(p: NotNaN<f64>) -> LogProb {
        LogProb(*p)
    }
}

impl From<LogProb> for NotNaN<f64> {
    fn from(p: LogProb) -> NotNaN<f64> {
        NotNaN::from(*p)
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

impl From<Prob> for PHREDProb {
    fn from(p: Prob) -> PHREDProb {
        PHREDProb(-10.0 * p.log10())
    }
}

impl From<LogProb> for PHREDProb {
    fn from(p: LogProb) -> PHREDProb {
        PHREDProb(*p * LOG_TO_PHRED_FACTOR)
    }
}

impl Default for LogProb {
    fn default() -> LogProb {
        LogProb::ln_zero()
    }
}

impl Default for PHREDProb {
    fn default() -> PHREDProb {
        PHREDProb::from(Prob(0.0))
    }
}

impl Zero for Prob {
    fn zero() -> Self {
        Prob(0.0)
    }

    fn is_zero(&self) -> bool {
        **self == 0.0
    }
}

impl Zero for LogProb {
    fn zero() -> Self {
        LogProb::ln_zero()
    }

    fn is_zero(&self) -> bool {
        self.is_infinite() && self.signum() < 0.0
    }
}

impl Zero for PHREDProb {
    fn zero() -> Self {
        PHREDProb::from(Prob(0.0))
    }

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

quick_error! {
    #[derive(Debug)]
    pub enum ProbError {
        InvalidProb(value: f64) {
            description("invalid probability")
            display("probabilty {} not in interval [0,1]", value)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let probs = vec![
            LogProb::ln_zero(),
            LogProb(0.01f64.ln()),
            LogProb(0.001f64.ln()),
        ];
        assert_eq!(
            LogProb::ln_cumsum_exp(probs).collect_vec(),
            [
                LogProb::ln_zero(),
                LogProb(0.01f64.ln()),
                LogProb(0.011f64.ln()),
            ]
        );
    }

    #[test]
    fn test_sub() {
        assert_eq!(
            LogProb::ln_one().ln_sub_exp(LogProb::ln_one()),
            LogProb::ln_zero()
        );
        assert_relative_eq!(
            *LogProb::ln_one().ln_sub_exp(LogProb(0.5f64.ln())),
            *LogProb(0.5f64.ln())
        );
        assert_relative_eq!(
            *LogProb(-1.6094379124341).ln_sub_exp(LogProb::ln_zero()),
            *LogProb(-1.6094379124341)
        );
    }

    #[test]
    fn test_one_minus() {
        assert_eq!(LogProb::ln_zero().ln_one_minus_exp(), LogProb::ln_one());
        assert_eq!(LogProb::ln_one().ln_one_minus_exp(), LogProb::ln_zero());
    }

    #[test]
    fn test_trapezoidal_integrate() {
        let density = |_| LogProb(0.1f64.ln());
        let prob = LogProb::ln_trapezoidal_integrate_exp(density, 0.0, 10.0, 5);
        assert_relative_eq!(*prob, *LogProb::ln_one(), epsilon = 0.0000001);
    }

    #[test]
    fn test_simpsons_integrate() {
        let density = |_| LogProb(0.1f64.ln());
        let prob = LogProb::ln_simpsons_integrate_exp(density, 0.0, 10.0, 5);
        assert_relative_eq!(*prob, *LogProb::ln_one(), epsilon = 0.0000001);
    }

    #[test]
    fn test_cap_numerical_overshoot() {
        let under_top = LogProb(-0.00000005);
        assert_eq!(under_top, under_top.cap_numerical_overshoot(0.0000001));
        let over_top = LogProb(0.00000005);
        assert_eq!(
            LogProb::ln_one(),
            over_top.cap_numerical_overshoot(0.0000001)
        );
    }

    #[test]
    #[should_panic]
    fn test_cap_numerical_overshoot_panic() {
        let over_top = LogProb(0.00000005);
        over_top.cap_numerical_overshoot(0.00000001);
    }

    #[test]
    fn test_sum_one_zero() {
        assert_eq!(
            LogProb::ln_one().ln_add_exp(LogProb::ln_zero()),
            LogProb::ln_one()
        );
    }

    #[test]
    fn test_zero() {
        assert!(LogProb::zero().is_zero());
        assert!(Prob::zero().is_zero());
        assert!(PHREDProb::zero().is_zero());
    }
}
