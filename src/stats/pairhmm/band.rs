// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Bookkeeping for visiting only the cells inside the band of a banded forward algorithm.
//!
//! The forward algorithms keep two alternating columns. A cell can only be inside the band if
//! its top left or left neighbour (previous column) or its top neighbour (this column) is.
//! Hence, for every column, the cells to examine are the ranges of cells written in the
//! previous column, extended by one to the right, the first cell if an alignment may start in
//! this column, and the chain of top neighbours of every cell found inside the band. Cells
//! outside of those ranges are never written, so the caller only has to clear what a column
//! buffer still holds from two columns ago (see [`BandScan::stale`]).

use serde::{Deserialize, Serialize};

#[derive(Default, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub(crate) struct BandScan {
    // inclusive ranges of the cells (position in y + 1) written in each of the two columns
    live: [Vec<(usize, usize)>; 2],
    // candidate ranges of the current column
    candidates: Vec<(usize, usize)>,
    len_y: usize,
    // cursor over the candidates
    range: usize,
    range_end: usize,
    cell: usize,
    next: usize,
    in_range: bool,
}

impl BandScan {
    /// Forget everything (at the start of a computation).
    pub(crate) fn reset(&mut self) {
        self.live[0].clear();
        self.live[1].clear();
    }

    /// The ranges of cells that column buffer `k` still holds from two columns ago; the
    /// caller has to clear them before writing the column.
    #[inline]
    pub(crate) fn stale(&self, k: usize) -> &[(usize, usize)] {
        &self.live[k]
    }

    /// Start column buffer `curr`, whose predecessor is buffer `prev`: forget what `curr` held
    /// and derive the candidate cells. `may_start` tells whether an alignment may start in
    /// this column, i.e. whether the first cell is a candidate.
    #[inline]
    pub(crate) fn begin_column(&mut self, curr: usize, prev: usize, len_y: usize, may_start: bool) {
        self.live[curr].clear();
        self.candidates.clear();
        // An empty y has no cell to start in: without this guard the first cell would be
        // read past the end of the row.
        if may_start && len_y > 0 {
            self.candidates.push((1, 1));
        }
        self.candidates.extend(
            self.live[prev]
                .iter()
                .map(|&(a, b)| (a, (b + 1).min(len_y))),
        );
        self.len_y = len_y;
        self.range = 0;
        self.next = 1;
        self.in_range = false;
    }

    /// The next cell (position in y + 1) to examine in this column, or `None` when the column
    /// is done. `last_inside` tells whether the previously returned cell turned out to be
    /// inside the band, in which case its bottom neighbour has to be examined too.
    #[inline]
    pub(crate) fn next_cell(&mut self, last_inside: bool) -> Option<usize> {
        if self.in_range {
            self.cell += 1;
            if self.cell <= self.len_y && (self.cell <= self.range_end || last_inside) {
                return Some(self.cell);
            }
            self.in_range = false;
            self.next = self.cell;
        }
        while self.range < self.candidates.len() {
            let (a, b) = self.candidates[self.range];
            self.range += 1;
            let cell = a.max(self.next);
            if cell > b {
                continue;
            }
            self.range_end = b;
            self.cell = cell;
            self.in_range = true;
            return Some(cell);
        }
        None
    }

    /// Record that cell `j_` of column buffer `curr` has been written.
    #[inline]
    pub(crate) fn mark(&mut self, curr: usize, j_: usize) {
        match self.live[curr].last_mut() {
            Some((_, last)) if *last + 1 == j_ => *last = j_,
            _ => self.live[curr].push((j_, j_)),
        }
    }
}
