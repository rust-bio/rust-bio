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
/// use bio::stats::combinatorics::scaled_combinations;
/// assert_eq!(scaled_combinations(5, 3, 0.5), 5.);
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
/// use bio::stats::combinatorics::combinations;
/// assert_eq!(combinations(5, 3), 10.);
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
/// use bio::stats::combinatorics::combinations_with_repl;
/// assert_eq!(combinations_with_repl(5, 3), 35.);
/// ```
pub fn combinations_with_repl(n: u64, k: u64) -> f64 {
    combinations(n + k - 1, k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comb() {
        assert_eq!(combinations(10, 3), 120.0);
        assert_eq!(combinations_with_repl(10, 3), 220.0);
        assert_eq!(combinations(200, 10), 22451004309013280.0);
    }

    #[test]
    fn test_comb_scaled() {
        assert!((scaled_combinations(150, 80, 1e-5) - 6.6643938163479384e+38).abs() < 0.0000001);
    }
}
