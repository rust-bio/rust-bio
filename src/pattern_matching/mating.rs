// Copyright 2018 Jeff K
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Merge overlapping paired sequences
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::mating::{mate,merge};
//! let read1 = b"tacgattcgat";
//! let read2 = b"ttcgattacgt";
//! let offset = mate(read1, read2, 0).unwrap();
//! let contig = merge(read1, read2, offset);
//! assert_eq!(contig, b"tacgattcgattacgt");
//! ```
//!
//! Mate pair merging procedes two steps:
//!   mate: identifying the optimal index of overlap or rejecting the case
//!   merge: combines overlapping sequences 
//!
//! Mating is governed by an objective function and merging is guided by a
//! function that resolves conflicts between sequences.
//!
//! We also need to define the conditions under which mating fails. This can
//! be guessed with a k-mer distribution or targetted with
//! information about fragment size statistics

use std::cmp;

pub fn mate(r1: &[u8], r2: &[u8], min_overhang: usize) -> Option<usize> {
    let min_score = 0; // TODO: move this to score closure
    let min_offset = min_overhang;
    let max_offset = cmp::min(r1.len(), r2.len()) - min_overhang;
    let mut m: i16 = 0;
    let mut pos: usize = 0;
    for i in min_offset..max_offset {
        let h = score(&r2[0..i], &r1[max_offset-i..max_offset]);
        if h > m {
            m = h;
            pos = i;
        }
    }
    if pos > min_offset && m > min_score {
        return Some(pos)
    }
    None
}

fn score(r1: &[u8], r2: &[u8]) -> i16 {
    let mut s: i16 = 0;
    let len = r1.len();
    for i in 0..len {
        if r1[i] == r2[i] {
            s = s + 1;
        } else {
            s = s - 1;
        }
    }
    s
}

pub fn merge(r1: &[u8], r2: &[u8], overlap: usize) -> Vec<u8> {
    let r1_end = r1.len() - overlap; 
    let r2_end = r2.len() - overlap;
    let len = r1_end + r2_end + overlap;

    let mut seq = vec![0; len];
    seq[0..r1_end].copy_from_slice(&r1[0..r1_end]);
    seq[r1_end..len].copy_from_slice(&r2[0..r2.len()]);
    seq
}

#[cfg(test)]
mod tests {
    use super::{mate, merge};

    #[test]
    fn test_mate_pair() {
        let r1 = b"tacgattcgat";
        let r2 = b"ttcgattacgt";
        let offset = mate(r1, r2, 0).unwrap();
        assert_eq!(offset, 6);
    }

    #[test]
    fn test_merge_pair() {
        let r1 = b"tacgattcgat";
        let r2 = b"ttcgattacgt";
        assert_eq!(merge(r1, r2, 6), b"tacgattcgattacgt");
    }
}
