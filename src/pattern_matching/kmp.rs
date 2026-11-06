// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Algorithm of Knuth Morris and Pratt.
//! Constructs an automaton recognizing the pattern, and scans linearly over
//! a text of length n. Complexity: O(n).
//! The transition function delta is simulated via the lps-function, that assigns to each position
//! q in the pattern the longest prefix of the pattern that is suffix of pattern[..q+1].
//! Then, in the NFA for the pattern, active states after reading position q are
//! {q, lps(q), lps(lps(q)), ... 0}.
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::kmp::KMP;
//! let text = b"aaaaabbabbbbbbbabbab";
//! let pattern = b"abbab";
//! let kmp = KMP::new(pattern);
//! let occ: Vec<usize> = kmp.find_all(text).collect();
//! assert_eq!(occ, [4, 15]);
//! ```

use std::borrow::Borrow;
use std::iter::{repeat_n, Enumerate};

use crate::utils::TextSlice;

type Lps = Vec<usize>;

/// KMP algorithm.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct KMP<'a> {
    m: usize,
    lps: Lps,
    pattern: TextSlice<'a>,
}

impl<'a> KMP<'a> {
    /// Create a new instance for a given pattern.
    pub fn new(pattern: TextSlice<'a>) -> Self {
        let m = pattern.len();
        let lps = lps(pattern);

        KMP { lps, m, pattern }
    }

    fn delta(&self, mut q: usize, a: u8) -> usize {
        while q == self.m || (self.pattern[q] != a && q > 0) {
            q = self.lps[q - 1];
        }
        if self.pattern[q] == a {
            q += 1;
        }

        q
    }

    /// Find all matches of pattern in a given text. Matches are returned as iterator over start
    /// positions.
    pub fn find_all<C, T>(&self, text: T) -> Matches<'_, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        Matches {
            kmp: self,
            q: 0,
            text: text.into_iter().enumerate(),
        }
    }
}

fn lps(pattern: &[u8]) -> Lps {
    let (m, mut q) = (pattern.len(), 0);
    let mut lps: Lps = repeat_n(0, m).collect();
    for i in 1..m {
        while q > 0 && pattern[q] != pattern[i] {
            q = lps[q - 1];
        }
        if pattern[q] == pattern[i] {
            q += 1;
        }
        lps[i] = q;
    }

    lps
}

/// Iterator over start positions of matches.
#[derive(Clone, Debug)]
pub struct Matches<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    kmp: &'a KMP<'a>,
    q: usize,
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
            self.q = self.kmp.delta(self.q, *c.borrow());
            if self.q == self.kmp.m {
                return Some(1 + i - self.kmp.m);
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::{lps, KMP};
    use itertools::Itertools;

    #[test]
    fn test_find_all() {
        let text = b"dhjalkjwqnnnannanaflkjdklfj";
        let pattern = b"qnnnannan";
        let kmp = KMP::new(pattern);
        assert_eq!(kmp.find_all(text).collect_vec(), [8]);
    }

    #[test]
    fn test_find_all_at_start() {
        let text = b"dhjalkjwqnnnannanaflkjdklfj";
        let pattern = b"dhjalk";
        let kmp = KMP::new(pattern);
        assert_eq!(kmp.find_all(text).collect_vec(), [0]);
    }

    #[test]
    fn test_lps() {
        let pattern = b"ababaca";
        let lps = lps(pattern);
        assert_eq!(lps, [0, 0, 1, 2, 3, 0, 1]);
    }

    #[test]
    fn test_delta() {
        let pattern = b"abbab";
        let kmp = KMP::new(pattern);
        assert_eq!(kmp.delta(0, b'a'), 1);
        assert_eq!(kmp.delta(0, b'b'), 0);
        assert_eq!(kmp.delta(1, b'a'), 1);
        assert_eq!(kmp.delta(1, b'b'), 2);
        assert_eq!(kmp.delta(2, b'a'), 1);
        assert_eq!(kmp.delta(2, b'b'), 3);
        assert_eq!(kmp.delta(3, b'a'), 4);
        assert_eq!(kmp.delta(3, b'b'), 0);
        assert_eq!(kmp.delta(4, b'a'), 1);
        assert_eq!(kmp.delta(4, b'b'), 5);
        assert_eq!(kmp.delta(5, b'a'), 1);
        assert_eq!(kmp.delta(5, b'b'), 3);
    }
}
