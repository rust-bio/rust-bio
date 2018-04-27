// Copyright 2015-2017 Vadim Nazarov, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various subroutines for computing a distance between sequences.

use std::cmp::min;

use utils::TextSlice;

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
pub fn hamming(alpha: TextSlice, beta: TextSlice) -> u64 {
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
/// let ldist = levenshtein(x, y);  // Distance is 5
/// assert_eq!(ldist, 5);
/// ```
#[allow(unused_assignments)]
pub fn levenshtein(alpha: TextSlice, beta: TextSlice) -> u32 {
    let mut columns = [vec![0u32; alpha.len() + 1], vec![0u32; alpha.len() + 1]];
    let mut i_prev = 0;
    let mut i_cur = 1;

    for i in 0..columns[0].len() {
        columns[0][i] = i as u32;
    }

    for (j, item) in beta.iter().enumerate() {
        i_cur %= 2;
        i_prev = 1 - i_cur;

        columns[i_cur][0] = 1 + j as u32;
        for i in 1..columns[0].len() {
            columns[i_cur][i] = min(
                columns[i_prev][i - 1] + if alpha[i - 1] == *item { 0 } else { 1 },
                min(columns[i_cur][i - 1] + 1, columns[i_prev][i] + 1),
            );
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
        assert_eq!(hamming(x, y), 5);
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
