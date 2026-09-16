// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Bounded version of Ukkonens DP algorithm for approximate pattern matching, extended
//! by optional support for swaps.
//! Complexity: O(n * k) on random texts of length n.
//!
//! The algorithm finds all matches of a pattern in a text with up to k errors.
//! The idea is to use dynamic programming to column-wise explore the edit matrix, but to omit
//! parts of the matrix for which the error exceeds k. To achieve this, a value `lastk` is
//! maintained that provides the lower feasible boundary of the matrix.
//! Initially, lastk = min(k, m). In each iteration (over a column), lastk can increase
//! by at most 1.
//!
//! Indels (cost 1 per inserted/deleted character, enabled by default) and swaps (CG -> GC with cost 1)
//! can be optionally enabled/disabled.
//! The cost function can be customized.
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
//! 
//! ukkonen.allow_swaps(true);
//! ukkonen.allow_indels(false);
//! let pattern = b"ABCD";
//! let text = b"ABDC";
//! let occ: Vec<(usize, usize)> = ukkonen.find_all_end(pattern, text, 1).collect();
//! assert_eq!(occ, vec![(3, 1)]);
//! ```

use std::borrow::Borrow;
use std::cmp::min;
use std::iter;
use std::iter::repeat_n;

use crate::utils::TextSlice;

/// Default cost function (unit costs).
pub fn unit_cost(a: u8, b: u8) -> u32 {
    (a != b) as u32
}

type ScoreFn<F> =
    fn(&Ukkonen<F>, u8, u8, Option<u8>, Option<u8>, usize, usize, usize, Option<usize>) -> usize;

/// Ukkonens algorithm.
#[allow(non_snake_case)]
#[derive(Default, Clone, Debug)]
pub struct Ukkonen<F>
where
    F: Fn(u8, u8) -> u32,
{
    D: [Vec<usize>; 3],
    cost: F,
    allow_swaps: bool,
    allow_indels: bool,
    score: Option<ScoreFn<F>>,
}

impl<F> Ukkonen<F>
where
    F: Fn(u8, u8) -> u32,
{
    /// Initialize algorithm with given capacity and cost function.
    pub fn with_capacity(m: usize, cost: F) -> Self {
        let get_vec = || Vec::with_capacity(m + 1);

        let mut instance = Ukkonen {
            D: [get_vec(), get_vec(), get_vec()],
            cost,
            allow_swaps: false,
            allow_indels: true,
            score: None,
        };
        instance.compile_score();

        instance
    }

    pub fn allow_swaps(&mut self, allow: bool) {
        self.allow_swaps = allow;
        self.compile_score();
    }

    pub fn allow_indels(&mut self, allow: bool) {
        self.allow_indels = allow;
        self.compile_score();
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
        let lastk = min(k, m);
        self.D[0].clear();
        self.D[0].extend(repeat_n(k + 1, m + 1));
        self.D[2].clear();
        if self.allow_indels {
            self.D[2].extend(0..=lastk);
            self.D[2].extend(repeat_n(k + 1, m - lastk));
        } else {
            self.D[2].push(0);
            self.D[2].extend(repeat_n(k + 1, m));
        }
        self.D[1].clear();
        self.D[1].extend(repeat_n(k + 1, m + 1));
        Matches {
            ukkonen: self,
            pattern,
            text: text.into_iter().enumerate(),
            lastk,
            m,
            k,
        }
    }

    fn score_substitution(&self, pattern_char: u8, text_char: u8, diag_score: usize) -> usize {
        let substitution_cost = (self.cost)(pattern_char, text_char) as usize;
        diag_score + substitution_cost
    }

    fn score_deletion(&self, left_score: usize) -> usize {
        left_score + 1
    }

    fn score_insertion(&self, up_score: usize) -> usize {
        up_score + 1
    }

    fn score_swap(
        &self,
        pattern_char: u8,
        text_char: u8,
        prev_pattern_char: Option<u8>,
        prev_text_char: Option<u8>,
        diag_score_2_hop: Option<usize>,
    ) -> usize {
        if let Some(diag_score_2_hop) = diag_score_2_hop {
            if let (Some(prev_pattern_char), Some(prev_text_char)) =
                (prev_pattern_char, prev_text_char)
            {
                if pattern_char == prev_text_char && text_char == prev_pattern_char {
                    return diag_score_2_hop + 1;
                }
            }
        }
        usize::MAX
    }

    fn score_only_substitution(
        &self,
        pattern_char: u8,
        text_char: u8,
        _prev_pattern_char: Option<u8>,
        _prev_text_char: Option<u8>,
        _up_score: usize,
        _left_score: usize,
        diag_score: usize,
        _diag_score_2_hop: Option<usize>,
    ) -> usize {
        self.score_substitution(pattern_char, text_char, diag_score)
    }

    fn score_with_swaps(
        &self,
        pattern_char: u8,
        text_char: u8,
        prev_pattern_char: Option<u8>,
        prev_text_char: Option<u8>,
        _up_score: usize,
        _left_score: usize,
        diag_score: usize,
        diag_score_2_hop: Option<usize>,
    ) -> usize {
        min(
            self.score_substitution(pattern_char, text_char, diag_score),
            self.score_swap(
                pattern_char,
                text_char,
                prev_pattern_char,
                prev_text_char,
                diag_score_2_hop,
            ),
        )
    }

    fn score_with_indels(
        &self,
        pattern_char: u8,
        text_char: u8,
        _prev_pattern_char: Option<u8>,
        _prev_text_char: Option<u8>,
        up_score: usize,
        left_score: usize,
        diag_score: usize,
        _diag_score_2_hop: Option<usize>,
    ) -> usize {
        min(
            min(
                self.score_substitution(pattern_char, text_char, diag_score),
                self.score_deletion(left_score),
            ),
            self.score_insertion(up_score),
        )
    }

    fn score_with_swaps_and_indels(
        &self,
        pattern_char: u8,
        text_char: u8,
        prev_pattern_char: Option<u8>,
        prev_text_char: Option<u8>,
        up_score: usize,
        left_score: usize,
        diag_score: usize,
        diag_score_2_hop: Option<usize>,
    ) -> usize {
        min(
            min(
                self.score_substitution(pattern_char, text_char, diag_score),
                self.score_deletion(left_score),
            ),
            min(
                self.score_insertion(up_score),
                self.score_swap(
                    pattern_char,
                    text_char,
                    prev_pattern_char,
                    prev_text_char,
                    diag_score_2_hop,
                ),
            ),
        )
    }

    fn compile_score(&mut self) {
        match (self.allow_swaps, self.allow_indels) {
            (true, true) => self.score = Some(Self::score_with_swaps_and_indels),
            (true, false) => self.score = Some(Self::score_with_swaps),
            (false, true) => self.score = Some(Self::score_with_indels),
            (false, false) => self.score = Some(Self::score_only_substitution),
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
        let mut prev_text_char = None;
        for (i, c) in &mut self.text {
            let col = i % 3;
            let prev = (i + 2) % 3;
            let prev_prev = (i + 1) % 3;
            let text_char = *c.borrow();

            // start with zero edit distance (semi-global alignment)
            self.ukkonen.D[col][0] = 0;
            self.lastk = min(self.lastk + 1, self.m);
            // in each column, go at most one cell further than before
            // do not look at cells with too big k
            let mut prev_pattern_char = None;
            for j in 1..=self.lastk {
                let pattern_char = self.pattern[j - 1];
                let diag_score = self.ukkonen.D[prev][j - 1];
                let diag_score_2_hop = if j > 1 && i >= 1 {
                    Some(self.ukkonen.D[prev_prev][j - 2])
                } else {
                    None
                };

                self.ukkonen.D[col][j] = self.ukkonen.score.unwrap()(
                    self.ukkonen,
                    pattern_char,
                    text_char,
                    prev_pattern_char,
                    prev_text_char,
                    self.ukkonen.D[col][j - 1],
                    self.ukkonen.D[prev][j],
                    diag_score,
                    diag_score_2_hop,
                );
                prev_pattern_char = Some(pattern_char);
            }

            if self.lastk < self.m {
                self.ukkonen.D[col][self.lastk + 1] = self.k + 1;
            }
            if self.ukkonen.allow_swaps && self.lastk + 1 < self.m {
                self.ukkonen.D[col][self.lastk + 2] = self.k + 1;
            }
            // reduce lastk as long as k is exceeded: while lastk can increase by at most 1, it can
            // decrease more in one iteration.
            while self.ukkonen.D[col][self.lastk] > self.k {
                self.lastk -= 1;
            }

            if self.lastk == self.m {
                return Some((i, self.ukkonen.D[col][self.m]));
            }
            prev_text_char = Some(text_char);
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

    #[test]
    fn test_allow_indels_toggle_changes_matches() {
        let pattern = b"GATTACA";
        let text = b"GATTTACA";

        let mut with_indels = Ukkonen::with_capacity(pattern.len(), unit_cost);
        with_indels.allow_indels(true);
        let occ_with_indels: Vec<(usize, usize)> =
            with_indels.find_all_end(pattern, text, 1).collect();
        assert!(occ_with_indels.contains(&(7, 1)));

        let mut without_indels = Ukkonen::with_capacity(pattern.len(), unit_cost);
        without_indels.allow_indels(false);
        let occ_without_indels: Vec<(usize, usize)> =
            without_indels.find_all_end(pattern, text, 1).collect();
        assert!(occ_without_indels.is_empty());
    }

    #[test]
    fn test_swaps_toggle_changes_matches() {
        let pattern = b"ABCD";
        let text = b"ABDC";

        let mut without_swaps = Ukkonen::with_capacity(pattern.len(), unit_cost);
        without_swaps.allow_indels(false);
        without_swaps.allow_swaps(false);
        let occ_without_swaps: Vec<(usize, usize)> =
            without_swaps.find_all_end(pattern, text, 1).collect();
        assert!(occ_without_swaps.is_empty());

        let mut with_swaps = Ukkonen::with_capacity(pattern.len(), unit_cost);
        with_swaps.allow_indels(false);
        with_swaps.allow_swaps(true);
        let occ_with_swaps: Vec<(usize, usize)> =
            with_swaps.find_all_end(pattern, text, 1).collect();
        assert_eq!(occ_with_swaps, vec![(3, 1)]);
    }

    #[test]
    fn test_swap_not_free_on_repeated_chars() {
        let pattern = b"AABA";
        let text = b"ABAB";

        let mut ukkonen = Ukkonen::with_capacity(pattern.len(), unit_cost);
        ukkonen.allow_indels(false);
        ukkonen.allow_swaps(true);

        // Optimal edit distance is 2 (one adjacent swap plus one substitution), not 1.
        let occ_k1: Vec<(usize, usize)> = ukkonen.find_all_end(pattern, text, 1).collect();
        assert!(occ_k1.is_empty());

        let occ_k2: Vec<(usize, usize)> = ukkonen.find_all_end(pattern, text, 2).collect();
        assert_eq!(occ_k2, vec![(3, 2)]);
    }
}
