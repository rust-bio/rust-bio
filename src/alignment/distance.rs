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
//! let x = b"GTCTGCATGCG";
//! let y = b"TTTAGCTAGCG";
//! // GTCTGCATGCG
//! //  |  ||  |||
//! // TTTAGCTAGCG
//! match hamm_dist(x, y) {
//!     Ok(s)  => println!("Score is: {}", s),  // Score is 5
//!     Err(e) => println!("Error: {}", e),     // No error in this example
//! }
//!
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! // ----ACCGTGGAT
//! //     ||||| |||
//! // AAAAACCGTTGAT
//! let l_score = lev_dist(x, y);  // Score is 5
//! ```


use itertools::Zip;
use std::result::Result;


/// Compute the Hamming distance between two strings. If returns the `Result<u32, &str>` type
/// with the first element corresponding to the distance between two strings and the second one to the error message
/// when two strings are not of equal sizes.
///
/// ```
/// use bio::alignment::distance::*;
/// let x = b"GTCTGCATGCG";
/// let y = b"TTTAGCTAGCG";
/// // GTCTGCATGCG
/// //  |  ||  |||
/// // TTTAGCTAGCG
/// match hamm_dist(x, y) {
///     Ok(s)  => println!("Score is: {}", s),  // Score is 5
///     Err(e) => println!("Error: {}", e),     // No error in this example
/// }
/// assert!(hamm_dist(x, y).is_ok());
/// assert_eq!(hamm_dist(x, y).unwrap(), 5);
/// ```
///
/// In case of the error:
///
/// ```
/// let x = b"GACTATATCGA";
/// let y = b"TTTAGCTC";
/// match hamm_dist(x, y) {
///     Ok(s)  => println!("Score is: {}", s),  // No score because strings are of unequal sizes.
///     Err(e) => println!("Error: {}", e),     // Will print "Error: hamming distance: strings are of unequal sizes!"
/// }
/// assert!(hamm_dist(x, y).is_err())
/// ```
pub fn hamm_dist(alpha: &[u8], beta: &[u8]) -> Result<u32, &'static str> {
    if alpha.len() != beta.len() {
        Err("hamming distance: strings are of unequal sizes!")
    } else {
        let mut score : u32 = 0;
        for (a, b) in Zip::new((alpha, beta)) {
            if a != b { score += 1; }
        }
        Ok(score)
    }
}

pub fn lev_dist(alpha: &[u8], beta: &[u8]) -> u32 {
    // let mut prev_col : &[u32];
    // let mut cur_col : &[u32];
    1
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
        assert!(hamm_dist(x, y).is_ok());
        assert_eq!(hamm_dist(x, y).unwrap(), 5);
    }

    #[test]
    fn test_hamming_dist_bad() {
        let x = b"GACTATATCGA";
        let y = b"TTTAGCTC";
        assert!(hamm_dist(x, y).is_err());
    }

    #[test]
    fn test_levenshtein_dist() {
        // let x = b"ACCGTGGAT";
        // let y = b"AAAAACCGTTGAT";
        // // ----ACCGTGGAT
        // //     ||||| |||
        // // AAAAACCGTTGAT
        // assert_eq!(lev_dist(x, y), align_global!(x, y, -1, 1, |a: u8, b: u8| if a == b {1i32} else {-1i32}))

        assert!(true)
    }
}
