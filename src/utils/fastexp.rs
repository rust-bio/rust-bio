// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides a trait adding a fast approximation of the exponential function to f64.
//! This can be very useful if the exact value is not too important, for example when working with
//! `bio::stats::LogProb`.

use num_traits::Float;
use std::f64;
use std::ops;

const COEFF_0: f64 = 1.0;
const COEFF_1: f64 = 4.831_794_110;
const COEFF_2: f64 = 0.143_440_676;
const COEFF_3: f64 = 0.019_890_581;
const COEFF_4: f64 = 0.006_935_931;
#[allow(clippy::approx_constant)]
const ONEBYLOG2: f64 = 1.442_695_041;
const OFFSET_F64: i64 = 1023;
const FRACTION_F64: u32 = 52;
const MIN_VAL: f64 = -500.0;

/// This trait adds a fast approximation of exp to float types.
pub trait FastExp<V: Float + ops::MulAssign> {
    fn fastexp(&self) -> V;
}

impl FastExp<f64> for f64 {
    /// Fast approximation of exp() as shown by Kopczynski 2017:
    /// https://eldorado.tu-dortmund.de/bitstream/2003/36203/1/Dissertation_Kopczynski.pdf
    fn fastexp(&self) -> f64 {
        if *self > MIN_VAL {
            let mut x = ONEBYLOG2 * self;

            let mut bits = x as i64;

            x -= bits as f64;
            let mut f2 = x;
            let mut x_tmp = x;

            bits += OFFSET_F64;
            bits <<= FRACTION_F64;

            f2 *= COEFF_4;
            x_tmp += COEFF_1;
            f2 += COEFF_3;
            x_tmp *= x;
            f2 *= x;
            f2 += COEFF_2;
            f2 *= x_tmp;
            f2 += COEFF_0;

            f64::from_bits(bits as u64) * f2
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastexp() {
        let x = 1e-15_f64.ln();
        assert_relative_eq!(x.fastexp(), 1e-15);
        let x = 1e-8_f64.ln();
        assert_relative_eq!(x.fastexp(), 1e-8, epsilon = 0.00000000000002);
        let x = 0.5_f64.ln();
        assert_relative_eq!(x.fastexp(), 0.5, epsilon = 0.01);
        let x = -159.000_000_023_278_6_f64;
        assert_relative_eq!(x.fastexp(), (-159.000_000_023_278_6).exp());
    }
}
