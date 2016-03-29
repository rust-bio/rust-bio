// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Myers bit-parallel approximate pattern matching algorithm.
//! Finds all matches up to a given edit distance. The pattern has to fit into a bitvector,
//! and is here limited to 64 symbols.
//! Complexity: O(n)
//!
//! # Example
//!
//! ```
//! # extern crate itertools;
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers;
//! use itertools::Itertools;
//!
//! # fn main() {
//! let text = b"ACCGTGGATGAGCGCCATAG";
//! let pattern = b"TGAGCGT";
//!
//! let myers = Myers::new(pattern);
//! let occ = myers.find_all_end(text, 1).collect_vec();
//!
//! assert_eq!(occ, [(13, 1), (14, 1)]);
//! # }
//! ```


use std::iter;
use std::u64;

use utils::{TextSlice, IntoTextIterator, TextIterator};


/// Myers algorithm.
pub struct Myers {
    peq: [u64; 256],
    bound: u64,
    m: u8,
}


impl Myers {
    /// Create a new instance of Myers algorithm for a given pattern.
    pub fn new(pattern: TextSlice) -> Self {
        assert!(pattern.len() <= 64 && pattern.len() > 0);

        let mut peq = [0; 256];
        for (i, &a) in pattern.iter().enumerate() {
            peq[a as usize] |= 1 << i;
        }

        Myers {
            peq: peq,
            bound: 1 << (pattern.len() - 1),
            m: pattern.len() as u8,
        }
    }

    /// Create a new instance of Myers algorithm for a given pattern and a wildcard character
    /// that shall match any character.
    pub fn with_wildcard(pattern: TextSlice, wildcard: u8) -> Self {
        let mut myers = Self::new(pattern);
        // wildcard matches all symbols of the pattern.
        myers.peq[wildcard as usize] = u64::MAX;

        myers
    }

    fn step(&self, state: &mut State, a: u8) {
        let eq = self.peq[a as usize];
        let xv = eq | state.mv;
        let xh = (((eq & state.pv) + state.pv) ^ state.pv) | eq;

        let mut ph = state.mv | !(xh | state.pv);
        let mut mh = state.pv & xh;

        if ph & self.bound > 0 {
            state.dist += 1;
        } else if mh & self.bound > 0 {
            state.dist -= 1;
        }

        ph <<= 1;
        mh <<= 1;
        state.pv = mh | !(xv | ph);
        state.mv = ph & xv;
    }

    /// Calculate the global distance of the pattern to the given text.
    pub fn distance<'a, I: IntoTextIterator<'a>>(&self, text: I) -> u8 {
        let mut state = State::new(self.m);
        for &a in text {
            self.step(&mut state, a);
        }
        state.dist
    }

    /// Find all matches of pattern in the given text up to a given maximum distance.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn find_all_end<'a, I: IntoTextIterator<'a>>(&'a self,
                                                        text: I,
                                                        max_dist: u8)
                                                        -> Matches<I::IntoIter> {
        Matches {
            myers: self,
            state: State::new(self.m),
            text: text.into_iter().enumerate(),
            max_dist: max_dist,
        }
    }
}


/// The current algorithm state.
struct State {
    pv: u64,
    mv: u64,
    dist: u8,
}


impl State {
    /// Create new state.
    pub fn new(m: u8) -> Self {
        State {
            pv: (1 << m) - 1,
            mv: 0,
            dist: m,
        }
    }
}


/// Iterator over pairs of end positions and distance of matches.
pub struct Matches<'a, I: TextIterator<'a>> {
    myers: &'a Myers,
    state: State,
    text: iter::Enumerate<I>,
    max_dist: u8,
}


impl<'a, I: Iterator<Item = &'a u8>> Iterator for Matches<'a, I> {
    type Item = (usize, u8);

    fn next(&mut self) -> Option<(usize, u8)> {
        for (i, &a) in self.text.by_ref() {
            self.myers.step(&mut self.state, a);
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance() {
        let text = b"TGAGCNT";
        let pattern = b"TGAGCGT";

        let myers = Myers::new(pattern);
        assert_eq!(myers.distance(text), 1);

        let myers_wildcard = Myers::with_wildcard(pattern, b'N');
        assert_eq!(myers_wildcard.distance(text), 0);
    }
}
