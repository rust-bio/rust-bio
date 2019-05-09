// Copyright 2014-2016 Johannes Köster.
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

use crate::pattern_matching::shift_and::masks;
use crate::utils::TextSlice;
use std::borrow::Borrow;

/// BNDM algorithm.
pub struct BNDM {
    m: usize,
    masks: [u64; 256],
    accept: u64,
}

impl BNDM {
    /// Create a new instance for a given pattern.
    pub fn new<C, P>(pattern: P) -> Self
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: DoubleEndedIterator + ExactSizeIterator,
    {
        let pattern = pattern.into_iter();
        let m = pattern.len();
        assert!(m <= 64, "Expecting a pattern of at most 64 symbols.");
        // take the reverse pattern and build nondeterministic
        // suffix automaton
        let (masks, accept) = masks(pattern.rev());

        BNDM { m, masks, accept }
    }

    /// Find all matches of pattern with a given text. Matches are returned as iterator over start positions.
    pub fn find_all<'a>(&'a self, text: TextSlice<'a>) -> Matches<'_> {
        Matches {
            bndm: self,
            window: self.m,
            text,
        }
    }
}

/// Iterator over start positions of matches.
pub struct Matches<'a> {
    bndm: &'a BNDM,
    window: usize,
    text: TextSlice<'a>,
}

impl<'a> Iterator for Matches<'a> {
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
                    } else {
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

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    #[test]
    fn test_find_all() {
        let text = b"dhjalkjwqnnnannanaflkjdklfj";
        let pattern = b"qnnnannan";
        let bndm = BNDM::new(pattern);
        assert_eq!(bndm.find_all(text).collect_vec(), [8]);
    }
}
