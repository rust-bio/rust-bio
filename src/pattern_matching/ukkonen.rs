// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Bounded version of Ukkonens DP algorithm for approximate pattern matching.
//! Complexity: O(n * k) on random texts of length n.
//!
//! The algorithm finds all matches of a pattern in a text with up to k errors.
//! The idea is to use dynamic programming to column-wise explore the edit matrix, but to omit
//! parts of the matrix for which the error exceeds k. To achieve this, a value `lastk` is
//! maintained that provides the lower feasible boundary of the matrix.
//! Initially, lastk = min(k, m). In each iteration (over a column), lastk can increase by at most 1.
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::ukkonen::{unit_cost, Ukkonen};
//!
//! let mut ukkonen = Ukkonen::with_capacity(10, unit_cost);
//! let text = b"ACCGTGGATGAGCGCCATAG";
//! let pattern = b"TGAGCGA";
//! let occ: Vec<(usize, usize)> = ukkonen.find_all_end(pattern, text, 1).collect();
//! assert_eq!(occ, [(13, 1), (14, 1)]);
//! ```

use std::borrow::Borrow;
use std::cmp::min;
use std::iter::{self, repeat_n};

use crate::utils::TextSlice;

/// Default cost function (unit costs).
pub fn unit_cost(a: u8, b: u8) -> u32 {
    (a != b) as u32
}

/// Ukkonens algorithm.
#[allow(non_snake_case)]
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Ukkonen<F>
where
    F: Fn(u8, u8) -> u32,
{
    D: [Vec<usize>; 2],
    cost: F,
}

impl<F> Ukkonen<F>
where
    F: Fn(u8, u8) -> u32,
{
    /// Initialize algorithm with given capacity and cost function.
    pub fn with_capacity(m: usize, cost: F) -> Self {
        let get_vec = || Vec::with_capacity(m + 1);
        Ukkonen {
            D: [get_vec(), get_vec()],
            cost,
        }
    }

    /// Find all matches between pattern and text with up to k errors.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn find_all_end<'a, C, T>(
        &'a mut self,
        pattern: TextSlice<'a>,
        text: T,
        k: usize,
    ) -> Matches<'a, F, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let m = pattern.len();
        self.D[0].clear();
        self.D[0].extend(repeat_n(k + 1, m + 1));
        self.D[1].clear();
        self.D[1].extend(0..=m);
        Matches {
            ukkonen: self,
            pattern,
            text: text.into_iter().enumerate(),
            lastk: min(k, m),
            m,
            k,
        }
    }
}

/// Iterator over pairs of end positions and distance of matches.
#[derive(Debug)]
pub struct Matches<'a, F, C, T>
where
    F: Fn(u8, u8) -> u32,
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    ukkonen: &'a mut Ukkonen<F>,
    pattern: TextSlice<'a>,
    text: iter::Enumerate<T>,
    lastk: usize,
    m: usize,
    k: usize,
}

impl<'a, F, C, T> Iterator for Matches<'a, F, C, T>
where
    F: 'a + Fn(u8, u8) -> u32,
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    type Item = (usize, usize);

    fn next(&mut self) -> Option<(usize, usize)> {
        let cost = &self.ukkonen.cost;
        for (i, c) in &mut self.text {
            let col = i % 2;
            let prev = 1 - col;

            // start with zero edit distance (semi-global alignment)
            self.ukkonen.D[col][0] = 0;
            self.lastk = min(self.lastk + 1, self.m);
            // in each column, go at most one cell further than before
            // do not look at cells with too big k
            for j in 1..=self.lastk {
                self.ukkonen.D[col][j] = min(
                    min(self.ukkonen.D[prev][j] + 1, self.ukkonen.D[col][j - 1] + 1),
                    self.ukkonen.D[prev][j - 1] + (cost)(self.pattern[j - 1], *c.borrow()) as usize,
                );
            }

            // reduce lastk as long as k is exceeded: while lastk can increase by at most 1, it can
            // decrease more in one iteration.
            while self.ukkonen.D[col][self.lastk] > self.k {
                self.lastk -= 1;
            }

            if self.lastk == self.m {
                return Some((i, self.ukkonen.D[col][self.m]));
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_all_end() {
        let mut ukkonen = Ukkonen::with_capacity(10, unit_cost);
        let text = b"ACCGTGGATGAGCGCCATAG";
        let pattern = b"TGAGCGT";
        let occ: Vec<(usize, usize)> = ukkonen.find_all_end(pattern, text, 1).collect();
        assert_eq!(occ, [(13, 1), (14, 1)]);
    }

    #[test]
    fn test_find_start() {
        let mut u = Ukkonen::with_capacity(10, unit_cost);

        let pattern = b"ACCGT";
        // hit begins at 1st position
        let text1 = b"ACCGTGGATGAGCGCCATAG";
        // hit begins at 2nd position
        let text2 = b"AACCGTGGATGAGCGCCATAG";

        let occ: Vec<(usize, usize)> = u.find_all_end(pattern, text1, 1).collect();
        assert_eq!(occ, [(3, 1), (4, 0), (5, 1)]);
        let occ: Vec<(usize, usize)> = u.find_all_end(pattern, text2, 1).collect();
        assert_eq!(occ, [(4, 1), (5, 0), (6, 1)]);
    }
}
