// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Algorithm of Horspool.
//! Window-based, similar to but faster than Boyer-Moore.
//!
//! # Idea
//! Look at a search window m, match pattern backwards.
//! In case of a mismatch, you can jump behind that.
//! Best case time complexity: O(n / m)
//! Worst case time complexity: O(n * m)
//! With a large alphabet, you are likely
//! around the best case, and faster than the rather
//! complicated Boyer-Moore.
//!
//! The algorithm has two phases (let a be the last symbol in the window):
//!
//! 1. test phase: compare the last symbol of the window.
//!    If it matches, compare the whole pattern.
//!    If it does not match, continue with the shift phase.
//! 2. shift phase: let l[a] be the rightmost position of a in
//!    the pattern without the last symbol. If it does not occur
//!    let l[a] be -1. Shift the window by m - 1 - l[a]. I.e.
//!    we shift the window such that the rightmost a matches
//!    the a at the end of the last window.
//!    If a does not occur in the pattern, we shift by the whole length.
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::horspool::Horspool;
//! let text = b"ACGGCTAGGAAAAAGACTGAGGACTGAAAA";
//! let pattern = b"GAAAA";
//! let horspool = Horspool::new(pattern);
//! let occ: Vec<usize> = horspool.find_all(text).collect();
//! assert_eq!(occ, [8, 25]);
//! ```

use crate::utils::TextSlice;

/// Algorithm of Horspool.
pub struct Horspool<'a> {
    shift: Vec<usize>,
    m: usize,
    pattern: TextSlice<'a>,
}

impl<'a> Horspool<'a> {
    /// Create a new instance for a given pattern.
    pub fn new(pattern: TextSlice<'a>) -> Self {
        let m = pattern.len();
        let mut shift = vec![m; 256];
        // shift is m for all not occurring characters
        // and m - 1 - j for all others
        for (j, &a) in pattern[..m - 1].iter().enumerate() {
            shift[a as usize] = m - 1 - j;
        }

        Horspool { m, shift, pattern }
    }

    /// Find all matches with a given text. Matches are returned as an iterator over start
    /// positions.
    pub fn find_all<'b>(&'b self, text: TextSlice<'b>) -> Matches<'_> {
        Matches {
            horspool: self,
            text,
            n: text.len(),
            last: self.m - 1,
            pattern_last: self.pattern[self.m - 1],
        }
    }
}

/// Iterator over start positions of matches.
pub struct Matches<'a> {
    horspool: &'a Horspool<'a>,
    text: TextSlice<'a>,
    n: usize,
    last: usize,
    pattern_last: u8,
}

impl<'a> Iterator for Matches<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        loop {
            // shift until the last symbol matches
            while self.last < self.n && self.text[self.last] != self.pattern_last {
                self.last += self.horspool.shift[self.text[self.last] as usize];
            }
            // stop if end of text is reached
            if self.last >= self.n {
                return None;
            }

            // putative start position
            let i = self.last + 1 - self.horspool.m;
            let j = self.last;

            // shift again (after both match and mismatch, this makes sense)
            self.last += self.horspool.shift[self.pattern_last as usize];

            if self.text[i..j] == self.horspool.pattern[..self.horspool.m - 1] {
                return Some(i);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Horspool;
    use itertools::Itertools;

    #[test]
    fn test_shift() {
        let pattern = b"AACB";
        let horspool = Horspool::new(pattern);
        assert_eq!(horspool.shift[b'A' as usize], 2);
        assert_eq!(horspool.shift[b'C' as usize], 1);
        assert_eq!(horspool.shift[b'B' as usize], 4);
        assert_eq!(horspool.shift[b'X' as usize], 4);
    }

    #[test]
    fn test_find_all() {
        let text = b"dhjalkjwqnnnannanaflkjdklfj";
        let pattern = b"qnnnannan";
        let horspool = Horspool::new(pattern);
        assert_eq!(horspool.find_all(text).collect_vec(), [8]);
    }
}
