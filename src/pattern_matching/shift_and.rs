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
use std::slice;


#[derive(Copy)]
pub struct ShiftAnd {
    m: usize,
    masks: [u64; 256],
    accept: u64
}


impl ShiftAnd {
    /// Create new ShiftAnd instance.
    pub fn new(pattern: &[u8]) -> ShiftAnd {
        assert!(pattern.len() <= 64, "Expecting a pattern of at most 64 symbols.");
        let (masks, accept) = get_masks(pattern);

        ShiftAnd { m: pattern.len(), masks: masks, accept: accept }

    }

    /// Find all occurences of pattern in the given text.
    pub fn find_all<'a>(&'a self, text: &'a [u8]) -> ShiftAndMatches {
        ShiftAndMatches { shiftand: self, active: 0, text: text.iter().enumerate() }
    }
}


pub fn get_masks(pattern: &[u8]) -> ([u64; 256], u64) {
    let mut masks = [0; 256];

    let mut bit = 1;
    for &c in pattern.iter() {
        masks[c as usize] |= bit;
        bit *= 2;
    }

    (masks, bit / 2)
}


pub struct ShiftAndMatches<'a> {
    shiftand: &'a ShiftAnd,
    active: u64,
    text: Enumerate<slice::Iter<'a, u8>>
}


impl<'a> Iterator for ShiftAndMatches<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        for (i, &c) in self.text {
            self.active = ((self.active << 1) | 1) & self.shiftand.masks[c as usize];
            if self.active & self.shiftand.accept > 0 {
                return Some(i - self.shiftand.m + 1);
            }
        }

        None
    }
}
