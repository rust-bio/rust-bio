// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Bounded version of Ukkonens DP algorithm for approximate pattern matching.
//! Complexity: O(n * k) on random texts.
//!
//! The algorithm finds all matches of a pattern in a text with up to k errors.
//! Idea is to use dynamic programming to column-wise explore the edit matrix, but to omit
//! parts of the matrix for which the error exceeds k. To achieve this, a value `lastk` is
//! maintained that provides the lower feasible boundary of the matrix.
//! Initially, lastk = min(k, m). In each iteration (over a column), lastk can increase by at most 1.
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::ukkonen::{Ukkonen, unit_cost};
//!
//! let mut ukkonen = Ukkonen::with_capacity(10, unit_cost);
//! let text = b"ACCGTGGATGAGCGCCATAG";
//! let pattern = b"TGAGCGA";
//! let occ: Vec<(usize, usize)> = ukkonen.find_all_end(pattern, text, 1).collect();
//! assert_eq!(occ, [(13, 1), (14, 1)]);
//! ```

use std::cmp::min;
use std::iter::repeat;

use utils::TextSlice;


/// Default cost function (unit costs).
pub fn unit_cost(a: u8, b: u8) -> u32 {
    (a != b) as u32
}


/// Ukkonens algorithm.
#[allow(non_snake_case)]
pub struct Ukkonen<F>
    where F: Fn(u8, u8) -> u32
{
    D: [Vec<usize>; 2],
    cost: F,
}


impl<F> Ukkonen<F>
    where F: Fn(u8, u8) -> u32
{
    /// Initialize algorithm with given capacity and cost function.
    pub fn with_capacity(m: usize, cost: F) -> Self {
        let get_vec = || Vec::with_capacity(m + 1);
        Ukkonen {
            D: [get_vec(), get_vec()],
            cost: cost,
        }
    }

    /// Find all matches between pattern and text with up to k errors.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn find_all_end<'a>(&'a mut self,
                            pattern: TextSlice<'a>,
                            text: TextSlice<'a>,
                            k: usize)
                            -> Matches<F> {
        let m = pattern.len();
        self.D[0].clear();
        self.D[0].extend(repeat(k + 1).take(m + 1));
        self.D[1].clear();
        self.D[1].extend(0..m + 1);
        Matches {
            ukkonen: self,
            pattern: pattern,
            text: text,
            i: 1,
            lastk: min(k, m),
            m: m,
            n: text.len(),
            k: k,
        }
    }
}


/// Iterator over pairs of end positions and distance of matches.
pub struct Matches<'a, F>
    where F: 'a + Fn(u8, u8) -> u32
{
    ukkonen: &'a mut Ukkonen<F>,
    pattern: TextSlice<'a>,
    text: TextSlice<'a>,
    i: usize,
    lastk: usize,
    m: usize,
    n: usize,
    k: usize,
}


impl<'a, F> Iterator for Matches<'a, F>
    where F: 'a + Fn(u8, u8) -> u32
{
    type Item = (usize, usize);

    fn next(&mut self) -> Option<(usize, usize)> {
        let cost = &self.ukkonen.cost;
        while self.i <= self.n {
            let col = self.i % 2;
            let prev = 1 - col;

            // start with zero edit distance (semi-global alignment)
            self.ukkonen.D[col][0] = 0;
            self.lastk = min(self.lastk + 1, self.m);
            // in each column, go at most one cell further than before
            // do not look at cells with too big k
            for j in 1..self.lastk + 1 {
                self.ukkonen.D[col][j] =
                    min(min(self.ukkonen.D[prev][j] + 1, self.ukkonen.D[col][j - 1] + 1),
                        self.ukkonen.D[prev][j - 1] +
                        (cost)(self.pattern[j - 1], self.text[self.i - 1]) as usize);
            }

            // reduce lastk as long as k is exceeded: while lastk can increase by at most 1, it can
            // decrease more in one iteration.
            while self.ukkonen.D[col][self.lastk] > self.k {
                self.lastk -= 1;
            }

            println!("{} {}", self.i, self.lastk);

            self.i += 1;

            if self.lastk == self.m {
                return Some((self.i - 2, self.ukkonen.D[col][self.m]));
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
}
