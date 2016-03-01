// Copyright 2014-2015 Johannes Köster, Vadim Nazarov.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Calculate alignments with a generalized variant of the Smith Waterman algorithm.
//! Complexity: O(n * m) for strings of length m and n.
//!
//! For quick computation of alignments and alignment scores there are 6 simple functions.
//!
//! # Example
//!
//! ```
//! use bio::alignment::pairwise::*;
//! use bio::alignment::AlignmentOperation::{Match, Subst};
//!
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
//! let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
//! let alignment = aligner.semiglobal(x, y);
//! assert_eq!(alignment.ystart, 4);
//! assert_eq!(alignment.xstart, 0);
//! assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
//!
//! // If you don't known sizes of future sequences, you could
//! // use Aligner::new().
//! // Global alignment:
//! let mut aligner = Aligner::new(-5, -1, &score);
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! let alignment = aligner.global(x, y);
//! assert_eq!(alignment.ystart, 0);
//! assert_eq!(alignment.xstart, 0);
//! assert_eq!(aligner.local(x, y).score, 7);
//! ```


use std::i32;
use std::iter::repeat;
use std::cmp::max;

use alignment::{Alignment, AlignmentOperation};
use data_structures::bitenc::BitEnc;


enum AlignmentType {
    Local,
    Semiglobal,
    Global,
}


/// Current internal state of alignment.
struct AlignmentState {
    m: usize,
    n: usize,
    best: i32,
    best_i: usize,
    best_j: usize,
    i: usize,
    j: usize,
    score: i32,
    col: usize,
}


macro_rules! align {
    (
        $aligner:ident, $x:ident, $y:ident, $state:ident,
        $init:block, $inner:block, $outer:block, $ret:block
    ) => (
        {
            let mut $state = AlignmentState {
                m: $x.len(), n: $y.len(),
                best: 0, best_i: 0, best_j: 0,
                i: 1, j: 1,
                score: 0,
                col: 0
            };

            while $state.i <= $state.n {
                $state.col = $state.i % 2;
                let prev = 1 - $state.col;

                // init code
                $init

                // read next y symbol
                let b = $y[$state.i - 1];
                $state.j = 1;

                while $state.j <= $state.m {
                    // read next x symbol
                    let a = $x[$state.j - 1];

                    // score for deletion
                    let d_score = max(
                        $aligner.S[prev][$state.j] + $aligner.gap_open,
                        $aligner.D[prev][$state.j] + $aligner.gap_extend
                    );
                    // score for insertion
                    let i_score = max(
                        $aligner.S[$state.col][$state.j-1] + $aligner.gap_open,
                        $aligner.I[$state.col][$state.j-1] + $aligner.gap_extend
                    );
                    // score for substitution
                    $state.score = $aligner.S[prev][$state.j-1] + ($aligner.score)(a, b);

                    if d_score > $state.score {
                        $state.score = d_score;
                        $aligner.traceback.del($state.i, $state.j);
                    }
                    else if i_score > $state.score {
                        $state.score = i_score;
                        $aligner.traceback.ins($state.i, $state.j);
                    }
                    else {
                        $aligner.traceback.subst($state.i, $state.j);
                    }

                    // inner code
                    $inner

                    $aligner.S[$state.col][$state.j] = $state.score;
                    $aligner.D[$state.col][$state.j] = d_score;
                    $aligner.I[$state.col][$state.j] = i_score;

                    $state.j += 1;
                }

                // outer code
                $outer

                $state.i += 1;
            }

            // return code
            $ret
        }
    );
}


/// A generalized Smith-Waterman aligner.
#[allow(non_snake_case)]
pub struct Aligner<'a, F>
    where F: 'a + Fn(u8, u8) -> i32
{
    S: [Vec<i32>; 2],
    I: [Vec<i32>; 2],
    D: [Vec<i32>; 2],
    traceback: Traceback,
    gap_open: i32,
    gap_extend: i32,
    score: &'a F,
}


const DEFAULT_ALIGNER_CAPACITY: usize = 200;


impl<'a, F> Aligner<'a, F>
    where F: Fn(u8, u8) -> i32
{
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `score` - function that returns the score for substitutions (also see bio::scores)
    ///
    pub fn new(gap_open: i32, gap_extend: i32, score: &'a F) -> Self {
        Aligner::with_capacity(DEFAULT_ALIGNER_CAPACITY,
                               DEFAULT_ALIGNER_CAPACITY,
                               gap_open,
                               gap_extend,
                               score)
    }

    /// Create new aligner instance. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `score` - function that returns the score for substitutions (also see bio::scores)
    ///
    pub fn with_capacity(m: usize, n: usize, gap_open: i32, gap_extend: i32, score: &'a F) -> Self {
        let get_vec = || Vec::with_capacity(m + 1);
        Aligner {
            S: [get_vec(), get_vec()],
            I: [get_vec(), get_vec()],
            D: [get_vec(), get_vec()],
            traceback: Traceback::with_capacity(m, n),
            gap_open: gap_open,
            gap_extend: gap_extend,
            score: score,
        }
    }

    /// Create new aligner instance with unit (equal to '-1') penalties for gap open and gap extend
    /// and unit score function ('1' if two letters are equal, '-1' if not). This is
    /// effectively equal to Levenshtein metric.
    // pub fn with_unit_cost() -> Self {
    //     let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    //     Aligner::new(-1, -1, *&score)
    // }

    fn init(&mut self, m: usize, alignment_type: AlignmentType) {
        // set minimum score to -inf, and allow to add gap_extend
        // without overflow
        let min_score = i32::MIN - self.gap_extend;
        for k in 0..2 {
            self.S[k].clear();
            self.I[k].clear();
            self.D[k].clear();
            self.I[k].extend(repeat(min_score).take(m + 1));
            self.D[k].extend(repeat(min_score).take(m + 1));
            match alignment_type {
                AlignmentType::Global | AlignmentType::Semiglobal => {
                    let mut s = &mut self.S[k];
                    // neutral start
                    s.push(0);
                    // other cells are gaps
                    let mut score = self.gap_open;
                    for _ in 1..m + 1 {
                        s.push(score);
                        score += self.gap_extend;
                    }
                }
                AlignmentType::Local => self.S[k].extend(repeat(0).take(m + 1)),
            }
        }
    }

    /// Calculate global alignment of x against y.
    pub fn global(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, AlignmentType::Global);
        self.traceback.init(m, n, AlignmentType::Global);

        align!(self,
               x,
               y,
               state,
               {
                   self.S[state.col][0] = self.gap_open + (state.i as i32 - 1) * self.gap_extend;
                   self.traceback.del(state.i, 0);
               },
               {},
               {},
               {
                   self.traceback.alignment(state.n, state.m, x, y, state.score)
               })
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, AlignmentType::Semiglobal);
        self.traceback.init(m, n, AlignmentType::Semiglobal);

        align!(self,
               x,
               y,
               state,
               {
                   self.S[state.col][0] = 0;
               },
               {},
               {
                   // the second condition ensures that score is overwritten if best
                   // does not reflect a full x-column (can happen in first iteration)
                   if state.score > state.best || state.best_j != state.m {
                       state.best = state.score;
                       state.best_i = state.i;
                       state.best_j = state.m;
                   }
               },
               {
                   self.traceback.alignment(state.best_i, state.best_j, x, y, state.best)
               })
    }

    /// Calculate local alignment of x against y.
    pub fn local(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, AlignmentType::Local);
        self.traceback.init(m, n, AlignmentType::Local);

        align!(self,
               x,
               y,
               state,
               {
                   self.S[state.col][0] = 0;
               },
               {
                   if state.score < 0 {
                       self.traceback.start(state.i, state.j);
                       state.score = 0;
                   } else if state.score > state.best {
                       state.best = state.score;
                       state.best_i = state.i;
                       state.best_j = state.j;
                   }
               },
               {},
               {
                   self.traceback.alignment(state.best_i, state.best_j, x, y, state.best)
               })
    }
}


/// Internal traceback.
struct Traceback {
    matrix: Vec<BitEnc>,
}


const TBSTART: u8 = 0b00;
const TBSUBST: u8 = 0b01;
const TBINS: u8 = 0b10;
const TBDEL: u8 = 0b11;


impl Traceback {
    fn with_capacity(m: usize, n: usize) -> Self {
        let mut matrix = Vec::with_capacity(n + 1);
        for _ in 0..n + 1 {
            matrix.push(BitEnc::with_capacity(2, m + 1));
        }
        Traceback { matrix: matrix }
    }

    fn init(&mut self, m: usize, n: usize, alignment_type: AlignmentType) {
        match alignment_type {
            AlignmentType::Global => {
                // set the first cell to start, the rest to insertions
                for i in 0..n + 1 {
                    self.matrix[i].clear();
                    self.matrix[i].push_values(m + 1, TBINS);
                }
                self.matrix[0].set(0, TBSTART);
            }
            AlignmentType::Semiglobal => {
                // set the first cell of each column to start, the rest to insertions
                for i in 0..n + 1 {
                    self.matrix[i].clear();
                    self.matrix[i].push_values(m + 1, TBINS);
                    self.matrix[i].set(0, TBSTART);
                }
            }
            AlignmentType::Local => {
                // set every cell to start
                for i in 0..n + 1 {
                    self.matrix[i].clear();
                    self.matrix[i].push_values(m + 1, TBSTART);
                }
            }
        }
    }

    fn start(&mut self, i: usize, j: usize) {
        self.matrix[i].set(j, TBSTART);
    }

    fn subst(&mut self, i: usize, j: usize) {
        self.matrix[i].set(j, TBSUBST);
    }

    fn del(&mut self, i: usize, j: usize) {
        self.matrix[i].set(j, TBDEL);
    }

    fn ins(&mut self, i: usize, j: usize) {
        self.matrix[i].set(j, TBINS);
    }

    fn get(&self, i: usize, j: usize) -> u8 {
        self.matrix[i].get(j).unwrap()
    }

    fn alignment(&self, mut i: usize, mut j: usize, x: &[u8], y: &[u8], score: i32) -> Alignment {
        let mut ops = Vec::with_capacity(x.len());

        loop {
            let (ii, jj, op) = match self.get(i, j) {
                TBSUBST => {
                    let op = if y[i - 1] == x[j - 1] {
                        AlignmentOperation::Match
                    } else {
                        AlignmentOperation::Subst
                    };
                    (i - 1, j - 1, op)
                }
                TBDEL => (i - 1, j, AlignmentOperation::Del),
                TBINS => (i, j - 1, AlignmentOperation::Ins),
                _ => {
                    break;
                }
            };

            ops.push(op);
            i = ii;
            j = jj;
        }

        ops.reverse();
        Alignment {
            ystart: i,
            xstart: j,
            xlen: x.len(),
            operations: ops,
            score: score,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use alignment::AlignmentOperation::{Match, Subst, Ins, Del};
    use scores::blosum62;

    #[test]
    fn test_semiglobal() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_local() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_global() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match,
                    Match]);
    }

    #[test]
    fn test_blosum62() {
        let x = b"AAAA";
        let y = b"AAAA";
        let score = &blosum62;
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.score, 16);
        assert_eq!(alignment.operations, [Match, Match, Match, Match]);
    }

    #[test]
    fn test_issue11() {
        let y = b"TACC";//GTGGAC";
        let x = b"AAAAACC";//GTTGACGCAA";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Ins, Ins, Ins, Subst, Match, Match, Match]);
    }


    #[test]
    fn test_issue12_1() {
        let x = b"CCGGCA";
        let y = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 1);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Subst, Subst, Subst]);
    }

    #[test]
    fn test_issue12_2() {
        let y = b"CCGGCA";
        let x = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.operations,
                   [Ins, Ins, Match, Subst, Ins, Ins, Ins, Ins, Subst, Match, Match]);
    }

    #[test]
    fn test_aligner_new() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::new(-5, -1, &score);

        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);

        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);

        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match,
                    Match]);
    }

    // #[test]
    // fn test_aligner_with_unit_cost() {
    //     let x = b"ACCGTGGAT";
    //     let y = b"AAAAACCGTTGAT";
    //     let mut aligner = Aligner::with_unit_cost();
    //     // ----ACCGTGGAT
    //     //     ||||| |||
    //     // AAAAACCGTTGAT
    //     assert_eq!(aligner.global(x, y).score, 5);
    //     assert_eq!(aligner.global(y, x).score, 5);
    //
    //     let x = b"AAA";
    //     let y = b"TTTT";
    //     assert_eq!(aligner.global(x, y).score, 4);
    //     assert_eq!(aligner.global(y, x).score, 4);
    // }
}
