// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
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

use alignment::{Alignment, AlignmentOperation};
use utils::TextSlice;

#[derive(Copy, Clone)]
enum AlignmentType {
    Local,
    Semiglobal,
    Global,
}


/// Value to use as a 'negative infinity' score. Should be close to i32::MIN,
/// but avoid underflow when used with reasonable scoring parameters. Use 0.9 * i32::MIN
const MIN_SCORE: i32 = -1932735283;

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

                    let mut tb = TracebackCell::new();

                    // score for deletion
                    let d_open = $aligner.S[prev][$state.j] + $aligner.gap_open;
                    let d_extend = $aligner.D[prev][$state.j] + $aligner.gap_extend;

                    if d_open > d_extend {
                        tb.set_d(TBSUBST);
                        $aligner.D[$state.col][$state.j] = d_open;
                    } else {
                        tb.set_d(TBDEL);
                        $aligner.D[$state.col][$state.j] = d_extend;
                    }

                    // score for insertion
                    let i_open = $aligner.S[$state.col][$state.j-1] + $aligner.gap_open;
                    let i_extend = $aligner.I[$state.col][$state.j-1] + $aligner.gap_extend;

                    if i_open > i_extend {
                        tb.set_i(TBSUBST);
                        $aligner.I[$state.col][$state.j] = i_open;
                    } else {
                        tb.set_i(TBINS);
                        $aligner.I[$state.col][$state.j] = i_extend;
                    };


                    // score for substitution
                    let match_score = ($aligner.score)(a, b);
                    $state.score = $aligner.S[prev][$state.j-1] + match_score;
                    tb.set_s(TBSUBST);

                    let from_d = $aligner.D[prev][$state.j-1] + match_score;
                    let from_i = $aligner.I[prev][$state.j-1] + match_score;

                    if from_d > $state.score {
                        $state.score = from_d;
                        tb.set_s(TBDEL);
                    }

                    if from_i > $state.score {
                        $state.score = from_i;
                        tb.set_s(TBINS);
                    }

                    $aligner.traceback.set($state.i, $state.j, tb);

                    // inner code
                    $inner

                    $aligner.S[$state.col][$state.j] = $state.score;
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

        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

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

    fn init(&mut self, m: usize, n: usize, alignment_type: AlignmentType) {

        self.traceback.init(m, n, alignment_type);

        // set minimum score to -inf, and allow to add gap_extend
        // without overflow
        let min_score = MIN_SCORE;
        for k in 0..2 {
            self.S[k].clear();
            self.I[k].clear();
            self.D[k].clear();

            self.D[k].extend(repeat(min_score).take(m + 1));

            match alignment_type {
                AlignmentType::Semiglobal |
                AlignmentType::Global => {
                    let mut i = &mut self.I[k];

                    // need one insertion to establish a gap
                    i.push(min_score);

                    // other cells are gaps
                    let mut score = self.gap_open;
                    for _ in 1..m + 1 {
                        i.push(score);
                        score += self.gap_extend;
                    }

                    self.S[k].push(0);
                    // Impossible to reach S state after first position in first column
                    self.S[k].extend(repeat(min_score).take(m));
                },

                AlignmentType::Local => {
                    self.S[k].extend(repeat(0).take(m + 1));
                    self.I[k].extend(repeat(min_score).take(m + 1));
                },
            }
        }
    }

    /// Calculate global alignment of x against y.
    pub fn global(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, n, AlignmentType::Global);

        align!(self,
               x,
               y,
               state,
               {
                   self.S[state.col][0] = MIN_SCORE;
                   self.I[state.col][0] = MIN_SCORE;
                   self.D[state.col][0] = self.gap_open + (state.i as i32 - 1) * self.gap_extend;

                   self.traceback.get_mut(state.i, 0).set_all(TBDEL);
               },
               {},
               {},
               {
                   self.alignment(state.n, state.m, x, y, state.score)
               })
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, n, AlignmentType::Semiglobal);

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
                   self.alignment(state.best_i, state.best_j, x, y, state.best)
               })
    }

    /// Calculate local alignment of x against y.
    pub fn local(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, n, AlignmentType::Local);

        align!(self,
               x,
               y,
               state,
               {
                   self.S[state.col][0] = 0;
               },
               {
                   if state.score < 0 {
                       self.traceback.get_mut(state.i, state.j).set_s(TBSTART);
                       state.score = 0;
                   } else if state.score > state.best {
                       state.best = state.score;
                       state.best_i = state.i;
                       state.best_j = state.j;
                   }
               },
               {},
               {
                   self.alignment(state.best_i, state.best_j, x, y, state.best)
               })
    }

    fn alignment(&self, yend: usize, xend: usize, x: TextSlice, y: TextSlice, score: i32) -> Alignment {

        let mut i = yend;
        let mut j = xend;

        let mut ops = Vec::with_capacity(x.len());

        let get = move |i,j,ty| {
                match ty {
                    TBDEL => self.traceback.get(i,j).get_d(),
                    TBINS => self.traceback.get(i,j).get_i(),
                    _ => self.traceback.get(i,j).get_s(),
                }
        };

        let mut which_mat = TBSUBST;

        loop {
            let tb = get(i, j, which_mat);

            if tb == TBSTART {
                break;
            }

            let (ii, jj, op) = match which_mat {
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
            which_mat = tb;
        }

        ops.reverse();
        Alignment {
            ystart: i,
            xstart: j,
            yend: yend,
            xend: xend,
            xlen: x.len(),
            operations: ops,
            score: score,
        }
    }

    
    // Debugging helper function for visualizing traceback matrices
    #[allow(dead_code)]
    fn print_traceback_matrices(&self, i: usize, j: usize)
    {
        for tb in &[TBSUBST, TBINS, TBDEL] {
            for jj in 0..(j+1) {
                let mut s = String::new();
                for ii in 0..(i+1) {
                    match self.traceback.get(ii,jj).get(*tb) {
                        TBSUBST => s.push_str(" M"),
                        TBDEL => s.push_str(" D"),
                        TBINS => s.push_str(" I"),
                        TBSTART => s.push_str(" S"),
                        _ => (),
                    }
                }
                println!("{}", s);
            }
        }
    }
}

/// Packed representation of one cell of a Smith-Waterman traceback matrix.
/// Stores the D, S, and I traceback matrix values in a single byte.
#[derive(Copy, Clone)]
struct TracebackCell {
    v: u8,
}

const DPOS: u8 = 0x3;
const SPOS: u8 = 0x3 << 2;
const IPOS: u8 = 0x3 << 4;  

impl TracebackCell {
    /// Initialize a blank traceback cell
    pub fn new() -> TracebackCell {
        TracebackCell { v: 0 }
    }

    pub fn set_d(&mut self, value: u8) {
        self.v = self.v & !DPOS | (value)
    }

    pub fn get_d(&self) -> u8 {
        self.v & 0x3
    }

    pub fn set_s(&mut self, value: u8) {
        self.v = self.v & !SPOS | (value << 2)
    }

    pub fn get_s(&self) -> u8 {
        (self.v >> 2) & 0x3
    }

    pub fn set_i(&mut self, value: u8) {
        self.v = self.v & !IPOS | (value << 4)
    }

    pub fn get_i(&self) -> u8 {
        (self.v >> 4) & 0x3
    }

    pub fn get(&self, which: u8) -> u8 {
        match which {
            TBINS => self.get_i(),
            TBDEL => self.get_d(),
            _ => self.get_s(),
        }
    }

    /// Set all matrices to the same value.
    pub fn set_all(&mut self, value: u8) {
        self.set_d(value);
        self.set_s(value);
        self.set_i(value);
    }
}

/// Internal traceback.
struct Traceback {
    matrix: Vec<Vec<TracebackCell>>,
}


const TBSTART: u8 = 0b00;
const TBSUBST: u8 = 0b01;
const TBINS: u8 = 0b10;
const TBDEL: u8 = 0b11;

impl Traceback {
    fn with_capacity(m: usize, n: usize) -> Self {
        let mut matrix = Vec::with_capacity(n + 1);
        for _ in 0..n + 1 {
            matrix.push(Vec::with_capacity(m+1));
        }
        Traceback { matrix: matrix }
    }

    fn init(&mut self, m: usize, n: usize, alignment_type: AlignmentType) {
        let mut ins = TracebackCell::new();
        ins.set_all(TBINS);

        match alignment_type {
            AlignmentType::Global => {
                // set the first cell to start, the rest to insertions
                for i in 0..n + 1 {
                    self.matrix[i].clear();
                    for _ in 0 .. m+1 {
                        self.matrix[i].push(ins)
                    }
                }

                self.matrix[0][0].set_all(TBSTART);
            }

            AlignmentType::Semiglobal => {
                // set the first cell of each column to start, the rest to insertions
                for i in 0..n + 1 {
                    self.matrix[i].clear();
                    for _ in 0 .. m+1 {
                        self.matrix[i].push(ins);
                    }

                    self.matrix[i][0].set_all(TBSTART);
                }
            }

            AlignmentType::Local => {
                // set every cell to start
                let mut start = TracebackCell::new();
                start.set_all(TBSTART);

                for i in 0..n + 1 {
                    self.matrix[i].clear();
                    for _ in 0 .. m + 1 {
                        self.matrix[i].push(start);
                    }
                }
            }
        }
    }

    fn set(&mut self, i: usize, j: usize, v: TracebackCell) {
        self.matrix[i][j] = v;
    }

    fn get(&self, i: usize, j: usize) -> &TracebackCell  {
        self.matrix[i].get(j).unwrap()
    }

    fn get_mut(&mut self, i: usize, j: usize) -> &mut TracebackCell  {
        self.matrix[i].get_mut(j).unwrap()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use alignment::AlignmentOperation::{Match, Subst, Ins, Del};
    use scores::blosum62;
    use std::iter::repeat;

    #[test]
    fn traceback_cell() {
        let mut tb = super::TracebackCell::new();

        tb.set_d(1);
        assert_eq!(tb.get_d(), 1);

        tb.set_s(2);
        assert_eq!(tb.get_s(), 2);
        assert_eq!(tb.get_d(), 1);

        tb.set_i(3);
        assert_eq!(tb.get_i(), 3);
        assert_eq!(tb.get_s(), 2);
        assert_eq!(tb.get_d(), 1);
        

    }

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


    // Test case for underflow of the SW score.
    #[test]
    fn test_semiglobal_gap_open_lt_mismatch() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -5i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -2, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Del, Match, Ins, Match, Match, Match]);
    }


    #[test]
    fn test_global_affine_ins() {
        let x = b"ACGAGAACA";
        let y = b"ACGACA";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -3i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);

        println!("aln:\n{}", alignment.pretty(x, y));
        assert_eq!(alignment.operations, [Match, Match, Match, Ins, Ins, Ins, Match, Match, Match]);
    }

    #[test]
    fn test_global_affine_ins2() {
        let x = b"AGATAGATAGATAGGGAGTTGTGTAGATGATCCACAGT";
        let y = b"AGATAGATAGATGTAGATGATCCACAGT";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);

        println!("aln:\n{}", alignment.pretty(x, y));

        let mut correct = Vec::new();
        correct.extend(repeat(Match).take(11));
        correct.extend(repeat(Ins).take(10));
        correct.extend(repeat(Match).take(17));

        assert_eq!(alignment.operations, correct);
    }


    #[test]
    fn test_local_affine_ins2() {
        let x = b"ACGTATCATAGATAGATAGGGTTGTGTAGATGATCCACAG";
        let y =  b"CGTATCATAGATAGATGTAGATGATCCACAGT";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.local(x, y);

        assert_eq!(alignment.xstart, 1);
        assert_eq!(alignment.ystart, 0);
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

        println!("\naln:\n{}", alignment.pretty(x, y));
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
                    [Subst, Match, Ins, Ins, Ins, Ins, Ins, Ins, Subst, Match, Match]);
    }


    #[test]
    fn test_issue12_3() {
        let y = b"CCGTCCGGCAA";
        let x = b"AAAAACCGTTGACGCAA";
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
        assert_eq!(alignment.operations,
                         [Ins, Ins, Ins, Ins, Ins, Ins, Match, Subst, Subst, Match, Subst, Subst, Subst, Match, Match, Match, Match]);


        let mut aligner = Aligner::with_capacity(y.len(), x.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(y, x);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                [Match, Subst, Subst, Match, Subst, Subst, Subst, Match, Match, Match, Match]);
    }


    #[test]
    fn test_left_aligned_del() {
        let x = b"GTGCATCATGTG";
        let y = b"GTGCATCATCATGTG";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        println!("\naln:\n{}", alignment.pretty(x, y));

        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Del, Del, Del, Match, Match, Match, Match, Match, Match, Match, Match, Match]);
    }


    #[test]
    fn test_left_aligned_ins() {

        let x = b"GTGCATCATCATGTG";
        let y = b"GTGCATCATGTG";
        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        println!("\naln:\n{}", alignment.pretty(x, y));

        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Ins, Ins, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match]);
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
