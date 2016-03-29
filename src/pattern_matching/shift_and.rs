// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! ShiftAnd algorithm for pattern matching.
//! Patterns may contain at most 64 symbols.
//! Complexity: O(n) with text length n.
//!
//! # Example
//!
//! ```rust
//! use bio::pattern_matching::shift_and;
//! let pattern = b"AAAA";
//! let text = b"ACGGCTAGAAAAGGCTAG";
//! let shiftand = shift_and::ShiftAnd::new(pattern);
//! let occ = shiftand.find_all(text).next().unwrap();
//! assert_eq!(occ, 8);
//! ```

use std::iter::Enumerate;

use utils::{TextSlice, IntoTextIterator, TextIterator};

/// ShiftAnd algorithm.
pub struct ShiftAnd {
    m: usize,
    masks: [u64; 256],
    accept: u64,
}


impl ShiftAnd {
    /// Create new ShiftAnd instance from a given pattern.
    pub fn new(pattern: TextSlice) -> ShiftAnd {
        assert!(pattern.len() <= 64,
                "Expecting a pattern of at most 64 symbols.");
        let (masks, accept) = masks(pattern);

        ShiftAnd {
            m: pattern.len(),
            masks: masks,
            accept: accept,
        }

    }

    /// Find all matches of pattern in the given text. Matches are returned as an iterator
    /// over start positions.
    pub fn find_all<'a, I: IntoTextIterator<'a>>(&'a self, text: I) -> Matches<I::IntoIter> {
        Matches {
            shiftand: self,
            active: 0,
            text: text.into_iter().enumerate(),
        }
    }
}


/// Calculate ShiftAnd masks. This function is called automatically when instantiating
/// a new ShiftAnd for a given pattern.
pub fn masks(pattern: &[u8]) -> ([u64; 256], u64) {
    let mut masks = [0; 256];

    let mut bit = 1;
    for &c in pattern.iter() {
        masks[c as usize] |= bit;
        bit *= 2;
    }

    (masks, bit / 2)
}


/// Iterator over start positions of matches.
pub struct Matches<'a, I: TextIterator<'a>> {
    shiftand: &'a ShiftAnd,
    active: u64,
    text: Enumerate<I>,
}


impl<'a, I: Iterator<Item = &'a u8>> Iterator for Matches<'a, I> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        for (i, &c) in self.text.by_ref() {
            self.active = ((self.active << 1) | 1) & self.shiftand.masks[c as usize];
            if self.active & self.shiftand.accept > 0 {
                return Some(i - self.shiftand.m + 1);
            }
        }

        None
    }
}
