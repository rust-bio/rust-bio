// Copyright 2015 M. Rizky Luthfianto.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use ndarray;


lazy_static! {
	// Taken from https://github.com/seqan/seqan/blob/master/include%2Fseqan%2Fscore%2Fscore_matrix_data.h#L710
	// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
	static ref MAT: ndarray::Array2<i32> = ndarray::Array::from_shape_vec((27, 27), vec![
		 3,  0, -3,  0,  0, -4,  1, -2, -1, -2, -2, -2, -2,  0,  0,  1, -1, -2,  1,  1,  0,  0, -7, -4,  0,  0, -9,
		 0,  3, -5,  4,  3, -6,  0,  1, -3, -4,  0, -4, -3,  3, -1, -1,  1, -1,  1,  0, -1, -3, -6, -4,  2, -1, -9,
		-3, -5, 12, -6, -7, -6, -4, -4, -3, -5, -7, -7, -6, -5, -4, -4, -7, -4,  0, -3, -4, -2, -9,  0, -7, -4, -9,
		 0,  4, -6,  5,  4, -7,  0,  0, -3, -4,  0, -5, -4,  3, -1, -2,  2, -2,  0,  0, -1, -3, -8, -5,  3, -1, -9,
		 0,  3, -7,  4,  5, -7,  0,  0, -3, -4,  0, -4, -3,  2, -1, -1,  3, -2,  0, -1, -1, -2, -9, -5,  4, -1, -9,
		-4, -6, -6, -7, -7, 10, -6, -2,  1,  2, -7,  2,  0, -4, -3, -6, -6, -5, -4, -4, -3, -2,  0,  7, -6, -3, -9,
		 1,  0, -4,  0,  0, -6,  6, -3, -3, -4, -2, -5, -4,  0, -1, -1, -2, -4,  1,  0, -1, -2, -8, -6, -1, -1, -9,
		-2,  1, -4,  0,  0, -2, -3,  8, -3, -3, -1, -3, -3,  2, -1, -1,  3,  2, -1, -2, -1, -3, -3,  0,  2, -1, -9,
		-1, -3, -3, -3, -3,  1, -3, -3,  6,  4, -2,  2,  2, -2, -1, -3, -3, -2, -2,  0, -1,  4, -6, -2, -3, -1, -9,
		-2, -4, -5, -4, -4,  2, -4, -3,  4,  5, -3,  5,  3, -3, -2, -3, -3, -3, -3, -1, -2,  3, -4, -2, -3, -2, -9,
		-2,  0, -7,  0,  0, -7, -2, -1, -2, -3,  6, -4,  1,  1, -1, -2,  1,  4,  0,  0, -1, -3, -4, -5,  0, -1, -9,
		-2, -4, -7, -5, -4,  2, -5, -3,  2,  5, -4,  7,  4, -4, -2, -3, -2, -4, -4, -2, -2,  2, -2, -2, -3, -2, -9,
		-2, -3, -6, -4, -3,  0, -4, -3,  2,  3,  1,  4,  8, -2, -1, -3, -1, -1, -2, -1, -1,  2, -5, -3, -2, -1, -9,
		 0,  3, -5,  3,  2, -4,  0,  2, -2, -3,  1, -4, -2,  3,  0, -1,  1,  0,  1,  0,  0, -2, -5, -2,  1,  0, -9,
		 0, -1, -4, -1, -1, -3, -1, -1, -1, -2, -1, -2, -1,  0, -1, -1, -1, -1,  0,  0, -1, -1, -5, -3, -1, -1, -9,
		 1, -1, -4, -2, -1, -6, -1, -1, -3, -3, -2, -3, -3, -1, -1,  7,  0,  0,  1,  0, -1, -2, -7, -6, -1, -1, -9,
		-1,  1, -7,  2,  3, -6, -2,  3, -3, -3,  1, -2, -1,  1, -1,  0,  5,  1, -1, -1, -1, -3, -6, -5,  4, -1, -9,
		-2, -1, -4, -2, -2, -5, -4,  2, -2, -3,  4, -4, -1,  0, -1,  0,  1,  7, -1, -1, -1, -3,  2, -5,  0, -1, -9,
		 1,  1,  0,  0,  0, -4,  1, -1, -2, -3,  0, -4, -2,  1,  0,  1, -1, -1,  2,  2,  0, -1, -3, -3, -1,  0, -9,
		 1,  0, -3,  0, -1, -4,  0, -2,  0, -1,  0, -2, -1,  0,  0,  0, -1, -1,  2,  4,  0,  0, -6, -3, -1,  0, -9,
		 0, -1, -4, -1, -1, -3, -1, -1, -1, -2, -1, -2, -1,  0, -1, -1, -1, -1,  0,  0, -1, -1, -5, -3, -1, -1, -9,
		 0, -3, -2, -3, -2, -2, -2, -3,  4,  3, -3,  2,  2, -2, -1, -2, -3, -3, -1,  0, -1,  5, -8, -3, -2, -1, -9,
		-7, -6, -9, -8, -9,  0, -8, -3, -6, -4, -4, -2, -5, -5, -5, -7, -6,  2, -3, -6, -5, -8, 18, -1, -7, -5, -9,
		-4, -4,  0, -5, -5,  7, -6,  0, -2, -2, -5, -2, -3, -2, -3, -6, -5, -5, -3, -3, -3, -3, -1, 11, -5, -3, -9,
		 0,  2, -7,  3,  4, -6, -1,  2, -3, -3,  0, -3, -2,  1, -1, -1,  4,  0, -1, -1, -1, -2, -7, -5,  4, -1, -9,
		 0, -1, -4, -1, -1, -3, -1, -1, -1, -2, -1, -2, -1,  0, -1, -1, -1, -1,  0,  0, -1, -1, -5, -3, -1, -1, -9,
		-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9,  1
	]).unwrap();
}

#[inline]
fn lookup(a: u8) -> usize {
    if a == b'Y' {
        23 as usize
    } else if a == b'Z' {
        24 as usize
    } else if a == b'X' {
        25 as usize
    } else if a == b'*' {
        26 as usize
    } else {
        (a - 65) as usize
    }
}

pub fn pam200(a: u8, b: u8) -> i32 {
    let a = lookup(a);
    let b = lookup(b);

    MAT[(a, b)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pam200() {
        let score1 = pam200(b'A', b'A');
        assert_eq!(score1, 3);
        let score2 = pam200(b'*', b'*');
        assert_eq!(score2, 1);
        let score3 = pam200(b'A', b'*');
        assert_eq!(score3, -9);
        let score4 = pam200(b'Y', b'Z');
        assert_eq!(score4, -5);
        let score5 = pam200(b'X', b'X');
        assert_eq!(score5, -1);
        let score6 = pam200(b'X', b'Z');
        assert_eq!(score6, -1);
    }
}
