// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Backward nondeterministic DAWG matching (BNDM).
//! Best-case complexity: O(n / m) with pattern of length m <= 64 and text of length n.
//! Worst case complexity: O(n * m).
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::bndm;
//! let pattern = b"GAAAA";
//! let text = b"ACGGCTAGAAAAGGCTAGAAAA";
//! let bndm = bndm::BNDM::new(pattern);
//! let occ: Vec<usize> = bndm.find_all(text).collect();
//! assert_eq!(occ, [7, 17]);
//! ```


use pattern_matching::shift_and::masks;
use std::slice::SliceExt;


#[derive(Copy)]
pub struct BNDM {
    m: usize,
    masks: [u64; 256],
    accept: u64
}


impl BNDM {
    pub fn new(pattern: &[u8]) -> Self {
        let m = pattern.len();
        assert!(m <= 64, "Expecting a pattern of at most 64 symbols.");
        // take the reverse pattern and build nondeterministic
        // suffix automaton
        let mut rev = pattern.to_vec();
        rev.reverse();

        let (masks, accept) = masks(rev.as_slice());

        BNDM { m: m, masks: masks, accept: accept }
    }

    pub fn find_all<'a>(&'a self, text: &'a [u8]) -> BNDMMatches {
        BNDMMatches { bndm: self, window: self.m, text: text }
    }
}


pub struct BNDMMatches<'a> {
    bndm: &'a BNDM,
    window: usize,
    text: &'a [u8]
}


impl<'a> Iterator for BNDMMatches<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        while self.window <= self.text.len() {
            let mut occ = None;
            // bit mask of ones, all states active
            let mut active = (1u64 << self.bndm.m) - 1;
            let (mut j, mut lastsuffix) = (1, 0);
            // while not in fail state
            while active != 0 {
                // process j-th symbol from right
                active &= self.bndm.masks[self.text[self.window - j] as usize];
                if active & self.bndm.accept != 0 {
                    // reached accepting state
                    if j == self.bndm.m {
                        occ = Some(self.window - self.bndm.m);
                        break;
                    }
                    else {
                        // we reached the accepting state
                        // but not the end of the pattern
                        // hence, a suffix of the reverse pattern
                        // i.e. a prefix of the pattern of
                        // length j matches
                        // in case of a mismatch, we can shift 
                        // to this prefix
                        lastsuffix = j;
                    }
                }
                j += 1;
                active <<= 1;
            }
            // shift the window
            self.window += self.bndm.m - lastsuffix;
            if occ.is_some() {
                return occ;
            }
        }

        None
    }
}
