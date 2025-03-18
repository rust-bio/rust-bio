// Copyright 2014-2017 M. Rizky Luthfianto.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

lazy_static! {
    // taken from https://github.com/seqan/seqan/blob/master/include%2Fseqan%2Fscore%2Fscore_matrix_data.h#L327
    static ref MAT: ndarray::Array2<i32> = ndarray::Array::from_shape_vec((27, 27), vec![
         4, -2,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -1, -2,  0, -1, -1, -1,  1,  0,  0,
         0, -3, -2, -1,  0, -4,
        -2,  4, -3,  4,  1, -3, -1,  0, -3, -4,  0, -4, -3,  3, -1, -2,  0, -1,  0, -1, -1,
        -3, -4, -3,  1, -1, -4,
         0, -3,  9, -3, -4, -2, -3, -3, -1, -1, -3, -1, -1, -3, -2, -3, -3, -3, -1, -1, -2,
        -1, -2, -2, -3, -2, -4,
        -2,  4, -3,  6,  2, -3, -1, -1, -3, -4, -1, -4, -3,  1, -1, -1,  0, -2,  0, -1, -1,
        -3, -4, -3,  1, -1, -4,
        -1,  1, -4,  2,  5, -3, -2,  0, -3, -3,  1, -3, -2,  0, -1, -1,  2,  0,  0, -1, -1,
        -2, -3, -2,  4, -1, -4,
        -2, -3, -2, -3, -3,  6, -3, -1,  0,  0, -3,  0,  0, -3, -1, -4, -3, -3, -2, -2, -1,
        -1,  1,  3, -3, -1, -4,
         0, -1, -3, -1, -2, -3,  6, -2, -4, -4, -2, -4, -3,  0, -1, -2, -2, -2,  0, -2, -1,
        -3, -2, -3, -2, -1, -4,
        -2,  0, -3, -1,  0, -1, -2,  8, -3, -3, -1, -3, -2,  1, -1, -2,  0,  0, -1, -2, -1,
        -3, -2,  2,  0, -1, -4,
        -1, -3, -1, -3, -3,  0, -4, -3,  4,  3, -3,  2,  1, -3, -1, -3, -3, -3, -2, -1, -1,
         3, -3, -1, -3, -1, -4,
        -1, -4, -1, -4, -3,  0, -4, -3,  3,  3, -3,  3,  2, -3, -1, -3, -3, -3, -2, -1, -1,
         2, -3, -1, -3, -1, -4,
        -1,  0, -3, -1,  1, -3, -2, -1, -3, -3,  5, -2, -1,  0, -1, -1,  1,  2,  0, -1, -1,
        -2, -3, -2,  1, -1, -4,
        -1, -4, -1, -4, -3,  0, -4, -3,  2,  3, -2,  4,  2, -3, -1, -3, -2, -2, -2, -1, -1,
         1, -2, -1, -3, -1, -4,
        -1, -3, -1, -3, -2,  0, -3, -2,  1,  2, -1,  2,  5, -2, -1, -2,  0, -1, -1, -1, -1,
         1, -1, -1, -1, -1, -4,
        -2,  3, -3,  1,  0, -3,  0,  1, -3, -3,  0, -3, -2,  6, -1, -2,  0,  0,  1,  0, -1,
        -3, -4, -2,  0, -1, -4,
         0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1,
        -1, -2, -1, -1, -1, -4,
        -1, -2, -3, -1, -1, -4, -2, -2, -3, -3, -1, -3, -2, -2, -2,  7, -1, -2, -1, -1, -2,
        -2, -4, -3, -1, -2, -4,
        -1,  0, -3,  0,  2, -3, -2,  0, -3, -3,  1, -2,  0,  0, -1, -1,  5,  1,  0, -1, -1,
        -2, -2, -1,  3, -1, -4,
        -1, -1, -3, -2,  0, -3, -2,  0, -3, -3,  2, -2, -1,  0, -1, -2,  1,  5, -1, -1, -1,
        -3, -3, -2,  0, -1, -4,
         1,  0, -1,  0,  0, -2,  0, -1, -2, -2,  0, -2, -1,  1,  0, -1,  0, -1,  4,  1,  0,
        -2, -3, -2,  0,  0, -4,
         0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  1,  5,  0,
         0, -2, -2, -1,  0, -4,
         0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1,
        -1, -2, -1, -1, -1, -4,
         0, -3, -1, -3, -2, -1, -3, -3,  3,  2, -2,  1,  1, -3, -1, -2, -2, -3, -2,  0, -1,
         4, -3, -1, -2, -1, -4,
        -3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -3, -2, -1, -4, -2, -4, -2, -3, -3, -2, -2,
        -3, 11,  2, -3, -2, -4,
        -2, -3, -2, -3, -2,  3, -3,  2, -1, -1, -2, -1, -1, -2, -1, -3, -1, -2, -2, -2, -1,
        -1,  2,  7, -2, -1, -4,
        -1,  1, -3,  1,  4, -3, -2,  0, -3, -3,  1, -3, -1,  0, -1, -1,  3,  0,  0, -1, -1,
        -2, -3, -2,  4, -1, -4,
         0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1,
        -1, -2, -1, -1, -1, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4,  1
    ]).unwrap();
}

#[inline]
fn lookup(a: u8) -> usize {
    if a == b'Y' {
        23
    } else if a == b'Z' {
        24
    } else if a == b'X' {
        25
    } else if a == b'*' {
        26
    } else {
        (a - 65) as usize
    }
}

/// Return the BLOSUM62 substitution matrix score of [a, b]
///
/// # Example
///
/// ```
/// use bio::scores::blosum62;
/// assert_eq!(blosum62(b'H', b'A'), -2);
/// ```
pub fn blosum62(a: u8, b: u8) -> i32 {
    let a = lookup(a);
    let b = lookup(b);

    MAT[(a, b)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blosum62() {
        let score1 = blosum30(b'H', b'H');
        assert_eq!(score1, 8);
        let score2 = blosum62(b'O', b'*');
        assert_eq!(score2, -4);
        let score3 = blosum62(b'A', b'*');
        assert_eq!(score3, -4);
        let score4 = blosum62(b'*', b'*');
        assert_eq!(score4, 1);
        let score5 = blosum62(b'X', b'X');
        assert_eq!(score5, -1);
        let score6 = blosum62(b'X', b'Z');
        assert_eq!(score6, -1);
    }
}
