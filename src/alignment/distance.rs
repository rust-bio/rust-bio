// Copyright 2015 Vadim Nazarov.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various subroutines for computing a distance between sequences.
//! Complexity: O(n) for strings of length n for the Hamming distance;
//! O(n * m) for strings of length n and m for the Levenshtein (or edit) distance.
//!
//! # Example
//!
//! ```
//! use bio::alignment::distance::*;
//!
//! let x = b"GTCTGCATGCG";
//! let y = b"TTTAGCTAGCG";
//! // GTCTGCATGCG
//! //  |  ||  |||
//! // TTTAGCTAGCG
//! match hamming(x, y) {
//!     Ok(s)  => println!("Score is: {}", s),  // Score is 5
//!     Err(e) => println!("Error: {}", e),     // No error in this example
//! }
//!
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! // ----ACCGTGGAT
//! //     ||||| |||
//! // AAAAACCGTTGAT
//! let l_score = levenshtein(x, y);  // Score is 5
//! assert_eq!(l_score, 5);
//! ```


use itertools::Zip;
use std::cmp::min;
use std::result::Result;


/// Compute the Hamming distance between two strings with `hamming`. If returns the `Result<u32, &str>` type
/// with the first element corresponding to the distance between two strings (a number of mismatches) and the second one to the error message
/// when two strings are not of equal sizes.
///
/// # Example
///
/// ```
/// use bio::alignment::distance::*;
///
/// let x = b"GTCTGCATGCG";
/// let y = b"TTTAGCTAGCG";
/// // GTCTGCATGCG
/// //  |  ||  |||
/// // TTTAGCTAGCG
/// match hamming(x, y) {
///     Ok(s)  => println!("Score is: {}", s),  // Score is 5
///     Err(e) => println!("Error: {}", e),     // No error in this example
/// }
/// assert!(hamming(x, y).is_ok());
/// assert_eq!(hamming(x, y).unwrap(), 5);
/// ```
///
/// In case of the error:
///
/// ```
/// use bio::alignment::distance::*;
///
/// let x = b"GACTATATCGA";
/// let y = b"TTTAGCTC";
/// match hamming(x, y) {
///     Ok(s)  => println!("Score is: {}", s),  // No score because strings are of unequal sizes.
///     Err(e) => println!("Error: {}", e),     // Will print "Error: hamming distance: strings are of unequal sizes!"
/// }
/// assert!(hamming(x, y).is_err());
/// ```
pub fn hamming(alpha: &[u8], beta: &[u8]) -> Result<u32, &'static str> {
    if alpha.len() == beta.len() {
        let mut score : u32 = 0;
        for (a, b) in Zip::new((alpha, beta)) {
            if a != b { score += 1; }
        }
        Ok(score)
    } else {
        Err("hamming distance: strings are of unequal sizes!")
    }
}


/// Compute the Levenshtein (or Edit) distance between two strings with `levenshtein`. It returns a distance between two strings,
/// i.e. minimal number of mismatches, insertions and deletions between two strings.
///
/// # Example
///
/// ```
/// use bio::alignment::distance::*;
///
/// let x = b"ACCGTGGAT";
/// let y = b"AAAAACCGTTGAT";
/// // ----ACCGTGGAT
/// //     ||||| |||
/// // AAAAACCGTTGAT
/// let l_score = levenshtein(x, y);  // Score is 5
/// assert_eq!(l_score, 5);
/// ```
#[allow(unused_assignments)]
pub fn levenshtein(alpha: &[u8], beta: &[u8]) -> u32 {
    let mut columns = [vec!(0u32; alpha.len() + 1), vec!(0u32; alpha.len() + 1)];
    let mut i_prev = 0;
    let mut i_cur = 1;

    for i in 0..columns[0].len() {
        columns[0][i] = i as u32;
    }

    for (j, item) in beta.iter().enumerate() {
        i_cur = i_cur % 2;
        i_prev = 1 - i_cur;

        columns[i_cur][0] = 1 + j as u32;
        for i in 1..columns[0].len() {
            columns[i_cur][i] = min(columns[i_prev][i - 1] +
                                    if alpha[i - 1] == *item {
                                        0
                                    } else {
                                        1
                                    },
                                    min(columns[i_cur][i - 1] + 1, columns[i_prev][i] + 1));
        }

        i_cur += 1;
    }

    columns[i_cur - 1][columns[0].len() - 1]
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_dist_good() {
        let x = b"GTCTGCATGCG";
        let y = b"TTTAGCTAGCG";
        // GTCTGCATGCG
        //  |  ||  |||
        // TTTAGCTAGCG
        assert!(hamming(x, y).is_ok());
        assert_eq!(hamming(x, y).unwrap(), 5);
    }

    #[test]
    fn test_hamming_dist_bad() {
        let x = b"GACTATATCGA";
        let y = b"TTTAGCTC";
        assert!(hamming(x, y).is_err());
    }

    #[test]
    fn test_levenshtein_dist() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        // ----ACCGTGGAT
        //     ||||| |||
        // AAAAACCGTTGAT
        assert_eq!(levenshtein(x, y), 5);
        assert_eq!(levenshtein(x, y), levenshtein(y, x));
        assert_eq!(levenshtein(b"AAA", b"TTTT"), 4);
        assert_eq!(levenshtein(b"TTTT", b"AAA"), 4);
    }
}
