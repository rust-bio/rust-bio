use std::mem;
use std::ops::{Add, AddAssign, Deref, Div, Mul, Sub};

use num_traits::identities::{One, Zero};
use num_traits::Pow;

use crate::stats::Prob;
use crate::utils::FastExp;
use std::iter::Sum;

const LOGPROB_ZERO: LogProb = LogProb(f64::NEG_INFINITY);
const LOGPROB_ONE: LogProb = LogProb(0.0);

#[derive(PartialOrd, PartialEq, Clone, Copy, Debug)]
pub struct LogProb(pub f64);

impl LogProb {
    pub fn exp(self) -> f64 {
        self.0.exp()
    }
    pub fn fastexp(self) -> f64 {
        self.0.fastexp()
    }
}

impl From<Prob> for LogProb {
    fn from(p: Prob) -> LogProb {
        assert!(p.0 >= 0.0);
        LogProb(p.ln())
    }
}

impl Into<Prob> for LogProb {
    fn into(self) -> Prob {
        Prob(self.0.fastexp())
    }
}

impl Zero for LogProb {
    fn zero() -> Self {
        LOGPROB_ZERO
    }

    fn is_zero(&self) -> bool {
        self.0.is_infinite() && self.0.signum() < 0.0
    }
}

impl One for LogProb {
    fn one() -> Self {
        LOGPROB_ONE
    }
}

impl Add for LogProb {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        if other == Self::zero() {
            self
        } else {
            let (mut x, mut y) = (self.0, other.0);
            if x < y {
                mem::swap(&mut x, &mut y)
            };
            if x == f64::NEG_INFINITY {
                Self::zero()
            } else if x == f64::INFINITY {
                LogProb(f64::INFINITY)
            } else {
                LogProb(x + (y - x).fastexp().ln_1p())
            }
        }
    }
}

impl AddAssign for LogProb {
    fn add_assign(&mut self, other: LogProb) {
        *self = *self + other;
    }
}

impl Sum<Self> for LogProb {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let probs: Vec<LogProb> = iter.collect();
        if probs.is_empty() {
            LogProb::zero()
        } else {
            let mut pmax = probs[0];
            let mut imax = 0;
            for (i, &p) in probs.iter().enumerate().skip(1) {
                if p > pmax {
                    pmax = p;
                    imax = i;
                }
            }
            if pmax == Self::zero() {
                Self::zero()
            } else if *pmax == f64::INFINITY {
                LogProb(f64::INFINITY)
            } else {
                LogProb(
                    pmax.0
                        + (probs
                            .iter()
                            .enumerate()
                            .filter_map(|(i, &p)| {
                                if i == imax || p == Self::zero() {
                                    None
                                } else {
                                    Some((p.0 - pmax.0).fastexp())
                                }
                            })
                            .sum::<f64>())
                        .ln_1p(),
                )
            }
        }
    }
}

impl Sub for LogProb {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        if other.is_zero() {
            self
        } else if self.is_zero() {
            LogProb(-other.0.exp())
        } else if self.0 - other.0 <= f64::EPSILON {
            LogProb::zero()
        } else {
            LogProb((self.0.exp() - other.0.exp()).ln())
        }
    }
}

impl Mul for LogProb {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: Self) -> Self::Output {
        LogProb(self.0 + rhs.0)
    }
}

// impl Mul<f64> for LogProb {
//     type Output = LogProb;
//
//     fn mul(self, other: f64) -> Self::Output {
//         LogProb((self.0.exp() * other).ln())
//     }
// }

impl Div for LogProb {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: LogProb) -> Self::Output {
        LogProb(self.0 - other.0)
    }
}

// impl Div<f64> for LogProb {
//     type Output = LogProb;
//
//     fn div(self, other: f64) -> Self::Output {
//         LogProb((self.0.exp() / other).ln())
//     }
// }

impl Pow<f64> for LogProb {
    type Output = LogProb;

    fn pow(self, rhs: f64) -> Self::Output {
        LogProb(self.0 * rhs)
    }
}

impl Pow<usize> for LogProb {
    type Output = LogProb;

    fn pow(self, rhs: usize) -> Self::Output {
        LogProb(self.0 * rhs as f64)
    }
}

impl Deref for LogProb {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
