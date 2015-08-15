// Copyright 2015 Vadim Nazarov.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various subroutines for computing a distance between sequences.


use std::iter::Zip;


pub fn hamm_dist(alpha: &[u8], beta: &[u8]) -> Result<u32, bool> {
    if alpha.len() != beta.len() {
        Result(0, true)
    } else {
        let mut score : u32 = 0;
        for (a, b) in Zip::new(alpha, beta) {
            score += alpha == beta;
        }
        Result(score, false)
    }
}

pub fn lev_dist(alpha: &[u8], beta: &[u8]) -> u32 {

}

#[cfg(test)]
mod tests {

    #[test]
    fn test_hamming_dist() {
        // let x = b"GACTATATCG";
        // let y = b"TTTAGCTAGC";
        // // GTCTGCATGC
        // //  |  ||  ||
        // // TTTAGCTAGC
        // assert_eq!(hamm_dist(x, y), 5)

        assert!(true)
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
