// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! `ShiftAnd` algorithm for pattern matching.
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

use std::borrow::Borrow;
use std::iter::Enumerate;

/// `ShiftAnd` algorithm.
pub struct ShiftAnd {
    m: usize,
    masks: [u64; 256],
    accept: u64,
}

impl ShiftAnd {
    /// Create new ShiftAnd instance from a given pattern.
    pub fn new<C, P>(pattern: P) -> Self
    where
        P::IntoIter: ExactSizeIterator,
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
    {
        let pattern = pattern.into_iter();
        let m = pattern.len();
        assert!(m <= 64, "Expecting a pattern of at most 64 symbols.");
        let (masks, accept) = masks(pattern);

        ShiftAnd { m, masks, accept }
    }

    /// Find all matches of pattern in the given text. Matches are returned as an iterator
    /// over start positions.
    pub fn find_all<C, T>(&self, text: T) -> Matches<'_, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        Matches {
            shiftand: self,
            active: 0,
            text: text.into_iter().enumerate(),
        }
    }
}

/// Calculate ShiftAnd masks. This function is called automatically when instantiating
/// a new ShiftAnd for a given pattern.
pub fn masks<C, P>(pattern: P) -> ([u64; 256], u64)
where
    C: Borrow<u8>,
    P: IntoIterator<Item = C>,
{
    let mut masks = [0; 256];

    let mut bit = 1;
    for c in pattern {
        masks[*c.borrow() as usize] |= bit;
        bit *= 2;
    }

    (masks, bit / 2)
}

/// Iterator over start positions of matches.
pub struct Matches<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    shiftand: &'a ShiftAnd,
    active: u64,
    text: Enumerate<T>,
}

impl<'a, C, T> Iterator for Matches<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        for (i, c) in self.text.by_ref() {
            self.active = ((self.active << 1) | 1) & self.shiftand.masks[*c.borrow() as usize];
            if self.active & self.shiftand.accept > 0 {
                return Some(i - self.shiftand.m + 1);
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
        let shiftand = ShiftAnd::new(pattern);
        assert_eq!(shiftand.find_all(text).collect_vec(), [8]);
    }
}
