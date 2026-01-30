// Copyright 2015-2017 Vadim Nazarov, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various subroutines for computing a distance between sequences. Features
//! both scalar and efficient vectorized distance functions with SIMD.

use crate::utils::TextSlice;

/// Compute the Hamming distance between two strings. Complexity: O(n).
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
/// assert_eq!(hamming(x, y), 5);
/// ```
pub fn hamming(alpha: TextSlice<'_>, beta: TextSlice<'_>) -> u64 {
    assert_eq!(
        alpha.len(),
        beta.len(),
        "hamming distance cannot be calculated for texts of different length ({}!={})",
        alpha.len(),
        beta.len()
    );
    let mut dist = 0;
    for (a, b) in alpha.iter().zip(beta) {
        if a != b {
            dist += 1;
        }
    }
    dist
}

/// Compute the Levenshtein (or Edit) distance between two strings. Complexity: O(n * m) with
/// n and m being the length of the given texts.
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
/// let ldist = levenshtein(x, y); // Distance is 5
/// assert_eq!(ldist, 5);
/// ```
#[allow(unused_assignments)]
pub fn levenshtein(alpha: TextSlice<'_>, beta: TextSlice<'_>) -> u32 {
    editdistancek::edit_distance(alpha, beta) as u32
}

pub mod simd {
    //! String distance routines accelerated with Single Instruction Multiple Data (SIMD)
    //! intrinsics.
    //!
    //! These routines will automatically fallback to scalar versions if AVX2 or SSE4.1 is
    //! not supported by the CPU.
    //!
    //! With AVX2, SIMD-accelerated Hamming distance can reach up to 40 times faster than
    //! the scalar version on strings that are long enough.
    //!
    //! The performance of SIMD-accelerated Levenshtein distance depends on the number of
    //! edits between two strings, so it can perform anywhere from 2 times to nearly 1000
    //! times faster than the scalar version. When the two strings are completely different,
    //! there could be no speedup at all. It is important to note that the algorithms work
    //! best when the number of edits is known to be small compared to the length of the
    //! strings (for example, 10% difference). This should be applicable in many situations.
    //!
    //! If AVX2 support is not available, there is a speed penalty for using SSE4.1 with
    //! smaller vectors.

    use crate::utils::TextSlice;
    use std::cmp::{max, min};

    /// SIMD-accelerated Hamming distance between two strings. Complexity: O(n / w), for
    /// SIMD vectors of length w (usually w = 16 or w = 32).
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alignment::distance::simd::*;
    ///
    /// let x = b"GTCTGCATGCG";
    /// let y = b"TTTAGCTAGCG";
    /// // GTCTGCATGCG
    /// //  |  ||  |||
    /// // TTTAGCTAGCG
    /// assert_eq!(hamming(x, y), 5);
    /// ```
    pub fn hamming(alpha: TextSlice<'_>, beta: TextSlice<'_>) -> u64 {
        assert_eq!(
            alpha.len(),
            beta.len(),
            "simd hamming distance cannot be calculated for texts of different length ({}!={})",
            alpha.len(),
            beta.len()
        );
        // triple_accel Hamming routine returns an u32
        triple_accel::hamming(alpha, beta) as u64
    }

    /// SIMD-accelerated Levenshtein (or Edit) distance between two strings. Complexity:
    /// O(k / w * (n + m)), with n and m being the length of the given texts, k being the
    /// number of edits, and w being the length of the SIMD vectors (usually w = 16 or
    /// w = 32).
    ///
    /// Uses exponential search, which is approximately two times slower than the usual
    /// O(n * m) implementation if the number of edits between the two strings is very large,
    /// but much faster for cases where the edit distance is low (when less than half of the
    /// characters in the strings differ).
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alignment::distance::simd::*;
    ///
    /// let x = b"ACCGTGGAT";
    /// let y = b"AAAAACCGTTGAT";
    /// // ----ACCGTGGAT
    /// //     ||||| |||
    /// // AAAAACCGTTGAT
    /// let ldist = levenshtein(x, y); // Distance is 5
    /// assert_eq!(ldist, 5);
    /// ```
    pub fn levenshtein(alpha: TextSlice<'_>, beta: TextSlice<'_>) -> u32 {
        triple_accel::levenshtein_exp(alpha, beta)
    }

    /// SIMD-accelerated bounded Levenshtein (or Edit) distance between two strings.
    /// Complexity: O(k / w * (n + m)), with n and m being the length of the given texts,
    /// k being the threshold on the number of edits, and w being the length of the SIMD vectors
    /// (usually w = 16 or w = 32).
    ///
    /// If the Levenshtein distance between two strings is greater than the threshold k, then
    /// `None` is returned. This is useful for efficiently calculating whether two strings
    /// are similar.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alignment::distance::simd::*;
    ///
    /// let x = b"ACCGTGGAT";
    /// let y = b"AAAAACCGTTGAT";
    /// // ----ACCGTGGAT
    /// //     ||||| |||
    /// // AAAAACCGTTGAT
    /// let ldist = bounded_levenshtein(x, y, 5); // Distance is 5
    /// assert_eq!(ldist, Some(5));
    ///
    /// let ldist = bounded_levenshtein(x, y, 4); // Threshold too low!
    /// assert_eq!(ldist, None);
    /// ```
    pub fn bounded_levenshtein(alpha: TextSlice<'_>, beta: TextSlice<'_>, k: u32) -> Option<u32> {
        editdistancek::edit_distance_bounded(
            alpha,
            beta,
            min(k as usize, max(alpha.len(), beta.len())),
        )
        .map(|x| x as u32)
    }
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
        assert_eq!(hamming(x, y), 5);
    }

    #[test]
    fn test_simd_hamming_dist_good() {
        let x = b"GTCTGCATGCG";
        let y = b"TTTAGCTAGCG";
        // GTCTGCATGCG
        //  |  ||  |||
        // TTTAGCTAGCG
        assert_eq!(simd::hamming(x, y), 5);
    }

    #[test]
    #[should_panic(
        expected = "hamming distance cannot be calculated for texts of different length (11!=8)"
    )]
    fn test_hamming_dist_bad() {
        let x = b"GACTATATCGA";
        let y = b"TTTAGCTC";
        hamming(x, y);
    }

    #[test]
    #[should_panic(
        expected = "simd hamming distance cannot be calculated for texts of different length (11!=8)"
    )]
    fn test_simd_hamming_dist_bad() {
        let x = b"GACTATATCGA";
        let y = b"TTTAGCTC";
        simd::hamming(x, y);
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

    #[test]
    fn test_simd_levenshtein_dist() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        // ----ACCGTGGAT
        //     ||||| |||
        // AAAAACCGTTGAT
        assert_eq!(simd::levenshtein(x, y), 5);
        assert_eq!(simd::levenshtein(x, y), simd::levenshtein(y, x));
        assert_eq!(simd::levenshtein(b"AAA", b"TTTT"), 4);
        assert_eq!(simd::levenshtein(b"TTTT", b"AAA"), 4);
    }

    #[test]
    fn test_simd_bounded_levenshtein_dist() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        // ----ACCGTGGAT
        //     ||||| |||
        // AAAAACCGTTGAT
        assert_eq!(simd::bounded_levenshtein(x, y, u32::MAX), Some(5));
        assert_eq!(
            simd::bounded_levenshtein(x, y, u32::MAX),
            simd::bounded_levenshtein(y, x, u32::MAX)
        );
        assert_eq!(
            simd::bounded_levenshtein(b"AAA", b"TTTT", u32::MAX),
            Some(4)
        );
        assert_eq!(
            simd::bounded_levenshtein(b"TTTT", b"AAA", u32::MAX),
            Some(4)
        );
    }
}
