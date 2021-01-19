// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Combinations with and without replacement.

use std::cmp;

/// Calculate the number of combinations when choosing
/// k elements from n elements without replacement, multiplied by a scaling factor.
/// Time complexity: O(min(k, n - k))
///
/// # Examples
/// ```
/// use approx::assert_relative_eq;
/// use bio::stats::combinatorics::scaled_combinations;
/// assert_relative_eq!(scaled_combinations(5, 3, 0.5), 5., epsilon = f64::EPSILON);
/// ```
pub fn scaled_combinations(n: u64, k: u64, scale: f64) -> f64 {
    if k > n {
        0.0
    } else {
        let mut comb = scale;
        for j in 0..cmp::min(k, n - k) {
            comb /= (j + 1) as f64;
            comb *= (n - j) as f64;
        }
        comb
    }
}

/// Calculate the number of combinations when choosing
/// k elements from n elements without replacement.
/// This is also known as n over k, or the binomial coefficient.
/// Time complexity: O(min(k, n - k))
///
/// # Examples
/// ```
/// use approx::assert_relative_eq;
/// use bio::stats::combinatorics::combinations;
/// assert_relative_eq!(combinations(5, 3), 10., epsilon = f64::EPSILON);
/// ```
pub fn combinations(n: u64, k: u64) -> f64 {
    scaled_combinations(n, k, 1.0)
}

/// Calculate the number of combinations when choosing
/// k elements from n elements with replacement.
/// Time complexity: O(min(k, n - k))
///
/// # Examples
/// ```
/// use approx::assert_relative_eq;
/// use bio::stats::combinatorics::combinations_with_repl;
/// assert_relative_eq!(combinations_with_repl(5, 3), 35., epsilon = f64::EPSILON);
/// ```
pub fn combinations_with_repl(n: u64, k: u64) -> f64 {
    combinations(n + k - 1, k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comb() {
        assert_relative_eq!(combinations(10, 3), 120.0, epsilon = f64::EPSILON);
        assert_relative_eq!(combinations_with_repl(10, 3), 220.0, epsilon = f64::EPSILON);
        assert_relative_eq!(
            combinations(200, 10),
            22451004309013280.0,
            epsilon = f64::EPSILON
        );
    }

    #[test]
    fn test_comb_scaled() {
        assert_relative_eq!(
            scaled_combinations(150, 80, 1e-5),
            6.664_393_816_347_938_4e38,
            epsilon = f64::EPSILON
        );
    }
}
