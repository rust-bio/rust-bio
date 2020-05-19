use std::mem;
use std::ops::{Add, AddAssign, Deref, Div, Mul, Sub};

use num_traits::identities::{One, Zero};
use num_traits::Pow;

use crate::stats::Prob;
use crate::utils::FastExp;

const LOGPROB_ZERO: LogProb = LogProb(f64::NEG_INFINITY);
const LOGPROB_ONE: LogProb = LogProb(0.0);

#[derive(PartialOrd, PartialEq, Clone, Copy, Debug)]
pub struct LogProb(pub f64);

impl LogProb {
    pub fn exp(self) -> f64 {
        self.0.exp()
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
        LogProb((self.0.exp() + other.0.exp()).ln())
    }
}

impl AddAssign for LogProb {
    fn add_assign(&mut self, other: LogProb) {
        *self = *self + other;
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

impl Mul<f64> for LogProb {
    type Output = LogProb;

    fn mul(self, other: f64) -> Self::Output {
        LogProb((self.0.exp() * other).ln())
    }
}

impl Div for LogProb {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: LogProb) -> Self::Output {
        LogProb(self.0 - other.0)
    }
}

impl Div<f64> for LogProb {
    type Output = LogProb;

    fn div(self, other: f64) -> Self::Output {
        LogProb((self.0.exp() / other).ln())
    }
}

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

/// Numerically stable sum.
///
/// See http://www.cs.cmu.edu/~quake-papers/robust-arithmetic.ps for a proof of
/// correctness.
pub struct ExactSum {
    partials: Vec<f64>,
}

impl ExactSum {
    fn new() -> ExactSum {
        ExactSum { partials: vec![] }
    }

    fn add(&mut self, mut x: f64) {
        if x.is_infinite() || x.is_nan() {
            return;
        }
        let mut j = 0;
        // This inner loop applies `hi`/`lo` summation to each
        // partial so that the list of partial sums remains exact.
        for i in 0..self.partials.len() {
            let mut y: f64 = self.partials[i];
            if x.abs() < y.abs() {
                mem::swap(&mut x, &mut y);
            }
            // Rounded `x+y` is stored in `hi` with round-off stored in
            // `lo`. Together `hi+lo` are exactly equal to `x+y`.
            let hi = x + y;
            let lo = y - (hi - x);
            if lo != 0.0 {
                self.partials[j] = lo;
                j += 1;
            }
            x = hi;
        }
        if j >= self.partials.len() {
            self.partials.push(x);
        } else {
            self.partials[j] = x;
            self.partials.truncate(j + 1);
        }
    }

    fn sum(&self) -> f64 {
        let mut hi;
        let l = self.partials.len();
        if l > 0 {
            hi = self.partials[l - 1];
            let mut n = l - 1;
            let mut lo = 0.0;
            let (mut x, mut y, mut yr);
            for n_ in 1..l {
                n = l - n_ - 1;
                x = hi;
                y = self.partials[n];
                assert!(y.abs() < x.abs());
                hi = x + y;
                yr = hi - x;
                lo = y - yr;
                if lo != 0.0 {
                    break;
                }
            }
            if n > 0
                && ((lo < 0.0 && self.partials[n - 1] < 0.0)
                    || (lo > 0.0 && self.partials[n - 1] > 0.0))
            {
                y = lo * 2.0;
                x = hi + y;
                yr = x - hi;
                if y == yr {
                    hi = x;
                }
            }
            return hi;
        }
        0.0
    }
}

impl std::iter::FromIterator<f64> for ExactSum {
    fn from_iter<T>(iter: T) -> ExactSum
    where
        T: IntoIterator<Item = f64>,
    {
        let mut e = ExactSum::new();
        for i in iter {
            e.add(i);
        }
        e
    }
}

#[cfg(test)]
mod tests {
    use super::ExactSum;

    #[test]
    fn test_exact_sum() {
        let values = vec![1e-16, 1., 1e16];
        let s: ExactSum = values.iter().copied().collect();
        assert_eq!(s.sum(), 10000000000000002.0);
    }
}
