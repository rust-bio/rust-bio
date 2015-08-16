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
//! // Easy-to-use functions for alignment.
//! // Global alignment:
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! let alignment = align_global(x, y, -5, -1, &score);
//! assert_eq!(alignment.ystart, 0);
//! assert_eq!(alignment.xstart, 0);
//!
//! // Score of the local alignment:
//! let score = score_local(x, y, -5, -1, &score);
//! assert_eq!(score, 7);
//! ```


use std::i32;
use std::iter::repeat;
use std::cmp::max;

use alignment::{Alignment, AlignmentOperation};
use data_structures::bitenc::BitEnc;


enum AlignmentType {
    Local,
    Semiglobal,
    Global
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
    col: usize
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
pub struct Aligner<'a, F> where F: 'a + Fn(u8, u8) -> i32 {
    S: [Vec<i32>; 2],
    I: [Vec<i32>; 2],
    D: [Vec<i32>; 2],
    traceback: Traceback,
    gap_open: i32,
    gap_extend: i32,
    score: &'a F,
}


impl<'a, F> Aligner<'a, F> where F: Fn(u8, u8) -> i32 {
    /// Create new aligner instance. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `score` - function that returns the score for substitutions
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
            score: score
        }
    }

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
                    let ref mut s = self.S[k];
                    let mut score = self.gap_open;
                    for _ in 0..m+1 {
                        s.push(score);
                        score += self.gap_extend;
                    }
                },
                AlignmentType::Local => {
                    self.S[k].extend(repeat(0).take(m + 1))
                }
            }
        }
    }

    /// Calculate global alignment of x against y.
    pub fn global(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, AlignmentType::Global);
        self.traceback.init(m, n, AlignmentType::Global);

        align!(
            self, x, y, state,
            {
                self.S[state.col][0] = self.gap_open + (state.i as i32 - 1) * self.gap_extend;
                self.traceback.del(state.i, 0);
            }, {}, {},
            {
                self.traceback.alignment(state.n, state.m, x, y, state.score)
            }
        )
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, AlignmentType::Semiglobal);
        self.traceback.init(m, n, AlignmentType::Semiglobal);

        align!(
            self, x, y, state,
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
            { self.traceback.alignment(state.best_i, state.best_j, x, y, state.best) }
        )
    }

    /// Calculate local alignment of x against y.
    pub fn local(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, AlignmentType::Local);
        self.traceback.init(m, n, AlignmentType::Local);

        align!(
            self, x, y, state,
            { self.S[state.col][0] = 0; },
            {
                if state.score < 0 {
                    self.traceback.start(state.i, state.j);
                    state.score = 0;
                }
                else if state.score > state.best {
                    state.best = state.score;
                    state.best_i = state.i;
                    state.best_j = state.j;
                }
            },
            {},
            { self.traceback.alignment(state.best_i, state.best_j, x, y, state.best) }
        )
    }
}


/// Internal traceback.
struct Traceback {
    matrix: Vec<BitEnc>
}


const TBSTART: u8 = 0b00;
const TBSUBST: u8 = 0b01;
const TBINS: u8   = 0b10;
const TBDEL: u8   = 0b11;


impl Traceback {

    fn with_capacity(m: usize, n: usize) -> Self {
        let mut matrix = Vec::with_capacity(n+1);
        for _ in 0..n+1 {
            matrix.push(BitEnc::with_capacity(2, m + 1));
        }
        Traceback {
            matrix: matrix
        }
    }

    fn init(&mut self, m: usize, n: usize, alignment_type: AlignmentType) {
        match alignment_type {
            AlignmentType::Global => {
                // set the first cell to start, the rest to insertions
                for i in 0..n+1 {
                    self.matrix[i].clear();
                    self.matrix[i].push_values(m + 1, TBINS);
                }
                self.matrix[0].set(0, TBSTART);
            },
            AlignmentType::Semiglobal => {
                // set the first cell of each column to start, the rest to insertions
                for i in 0..n+1 {
                    self.matrix[i].clear();
                    self.matrix[i].push_values(m + 1, TBINS);
                    self.matrix[i].set(0, TBSTART);
                }
            },
            AlignmentType::Local => {
                // set every cell to start
                for i in 0..n+1 {
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
                    let op = if y[i-1] == x[j-1] {
                        AlignmentOperation::Match
                    }
                    else {
                        AlignmentOperation::Subst
                    };
                    (i - 1, j - 1, op)
                }
                TBDEL => {
                    (i - 1, j, AlignmentOperation::Del)
                }
                TBINS => {
                    (i, j - 1, AlignmentOperation::Ins)
                }
                _ => {
                    break;
                }
            };

            ops.push(op);
            i = ii;
            j = jj;
        }

        ops.reverse();
        Alignment { ystart: i, xstart: j, xlen: x.len(), operations: ops, score: score}
    }
}


/// Compute global alignment of the given sequences and return an `Alignment`.
pub fn align_global<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> Alignment
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).global(x, y)
}


/// Compute semiglobal alignment of the given sequences and return an `Alignment`.
pub fn align_semiglobal<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> Alignment
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).semiglobal(x, y)
}


/// Compute local alignment of the given sequences and return an `Alignment`.
pub fn align_local<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> Alignment
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).local(x, y)
}


/// Compute global alignment of the given sequences and return the alignment score.
pub fn score_global<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> i32
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).global(x, y).score
}


/// Compute semiglobal alignment of the given sequences and return the alignment score.
pub fn score_semiglobal<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> i32
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).semiglobal(x, y).score
}


/// Compute local alignment of the given sequences and return the alignment score.
pub fn score_local<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> i32
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).local(x, y).score
}


/// Compute global alignment of the given sequences and return the pretty-formatted alignment as a `String`.
pub fn pretty_global<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> String
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).global(x, y).pretty(x, y)
}


/// Compute semiglobal alignment of the given sequences and return the pretty-formatted alignment as a `String`.
pub fn pretty_semiglobal<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> String
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).semiglobal(x, y).pretty(x, y)
}


/// Compute local alignment of the given sequences and return the pretty-formatted alignment as a `String`.
pub fn pretty_local<'a, F>(x: &[u8], y: &[u8], gap_open: i32, gap_extend: i32, score: &'a F) -> String
    where F : Fn(u8, u8) -> i32 {
    Aligner::with_capacity(x.len(), y.len(), gap_open, gap_extend, score).local(x, y).pretty(x, y)
}


#[cfg(test)]
mod tests {
    use super::*;
    use alignment::AlignmentOperation::{Match, Subst, Ins, Del};
    use test::Bencher;


    // 5,000 random nucleotides, GC content = .55
    static STR_1: &'static [u8] = b"ATCTAACTATTCCCTGTGCCTTATGGGGGCCTGCGCTATCTGCCTGTCGAACCATAGGACTCGCGCCAGCGCGCAGGCTTGGATCGAGGTGAAATCTCCGGGGCCTAAGACCACGAGCGTCTGGCGTCTTGGCTAACCCCCCTACATGCTGTTATAGACAATCAGTGGAAACCCGGTGCCAGGGGGTGGAGTGACCTTAAGTCAGGGACGATATTAATCGGAAGGAGTATTCAACGCAATGAAGCCGCAGGGTTGGCGTGGGAATGGTGCTTCTGTCCAAGCAGGTAAGGGCATGAGGCCGCAACCGTCCCCCAAGCGTACAGGGTGCACTTTGCAACGATTTCGGAGTCCGGAGACTCGCTGTTTTCGAAATTTGCGCTCAAGGGCGGGTATTGAACCAGGCTTACGCCCAAGAACGTAGCAAGGTGACTCAAACAAGGTACATCTTGCCCGCGTTTCACACGAATCAAGTTGGAGGTTATGGAGCATAGTAACACGTGGGCGGCCAGTGGTCGGTTGCTACACCCCTGCCGCAACGTTGAAGGTCCCGGATTAGACTGGCTGGACCCATGCCGTGACACCCGTCACACTCCATTACCGTCTGCGGGTCACGGCTTGTTGTGGACTGGATTGCCATTCTCTCAGTGTATTACGCAGGCCGGCGCGCGGGTCCCATGTAAACCTGTCATAGCTTACCTGACTCTACTTGGAAGTGTGGCTAGGCCTTTGCCCACGCACCTGGTCGGTCCTCGTTTGCTTTTTAGGACCGGATGAACTACAGAGCGCTGCAAGAATCTCTACCTGCTTTACAAAGCGCTGGGTCCTACTCCAGCGGGATGTTTTATCTAAACACGATGAGAGGAGTATTCGTCAGGCCACATGGCTTTCTTGTCCTGGTCGGATCCATCGTTGGCGCCCGACCCCCCCACTCCGTAGTGAGTTCTTCGTCCGAGCCATTGCATGCCAGATCGGCAGACAGATAGCGGATCCAGTATATCCCTGGAAGCTATAGACGCACAGGTTGGAATCCTAAGCGAAGTCGCGCGTCCGAACCCAGCTCTACTTTAGTGGCCACGGGTTCTGGTCCCCCCGGGCCGCGGAACCGATTAGGGCCATGTACAACAATACTTATTAGTCACCTTTCAGACACGATCTCCCTGCTCAGTGGTATATGGTTCCTGCTATAATTAGCCACCCTCATAAGTTGCACTACTTCTGCGACCCAAGTGCACCCTTACCACGAAGACAGGATTGTCCGATCCCATACTGCGGCCTTGGCAGGGGGTTCGCAAGTCCCACCCCAAACGATGCTGAAGGCTCAGGTTACACAGGCACAAGTGCTATATACGCGAGTTCCCGCTCTTAACCTGGACCGAATGCGGGATCATGCATCGTACCACTGTGTTCGTGTCATCTAGGACGGGCGCAAAGGATACATAGTTCAATCAAGAATACCTTGTATTATTGTACACCTACCGGTCACCAGCCAACAATGTGCGGACGGCGTTGCGACTTGCTGGGCCTGATCTCACCGCCCTAGATACCGCACACTGGGCAATACGAGGTAAAGCCAGTCACCCAGTGTCGATCAACAGCTGACGTAACGGTAAGAGGCTCACAAAATCGCACCGCCGGCGTCCCCTGGGTATTTTACGTCAGCATCGGGTGGACTGGCATGAATCTTTACTCCCAGGCGGAAACGGGTGCGTGGACAAGCGAGCAGCAAACGAAAATTCCTGGCCTGCTTGGTGTCTCGTATCCCTCTTGGAGATCGAGGAAATGTTTCACGACCAAGGGAAAGGTCGCCCTACGAAATAGATTTGCGCTACTGTCCGCATAAGGAGTCCGGTGTAGCGAAGGATGAAGGCGACCCTAGGTAGCAACCGCCGGCTTCGGCGGTAAGGTATCACTCAGGAAGCAGGCACGGAAAGACACGGTCTAGCAGACCGTCTATCGGCTAGGTCAAATAGGGTGCTTTGATATCAGCATGTCCAGCCTTAGAATTCAGTTCAGCGCGCTGGTCTGGGTCGAGATAAAATCACCAGTACCCAAGACCAGGCGGGCTCGCCGCGTTGGCTAATCCTGGTACATCTTGTAATCAATGTTCAGAAGAAAATCTGTGTTAGAGGGACGAGTCACCACGTACCAATAGCGACAACGATCGGTCGGACTATTCATCGTGGTGGTGACGCTCGGATTACGCGGGAAAGGTGCTTGTGTCCCGACAGGCTAGGATATAATGCTGAGGCGCTGCCCCAACCGTTCAGCGTGGGGTTTGCTACAACTTCCGAGTGCTACGTGTGCGAGACCATGTTATGTATGCACAAGGCCGACAATAGGACGTAGCCTTCGAGTTAGTACGTAGCGTGGTCGCACAAGCACAGTAGATCCTCCCCGCGCATCCTATTTATTAAGTTAATTCTATAGCAATACGATCACATGCGGATGGGCAGTGGCCGGTAGTCACACGCCTACCGCGGTGCTCAATGACCGGGACTAGAGAGGCGAAGATTATGGCGTGTGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGATTGCTTTCCCAATCTCCGAGCGATTTAGCGTGACAGCCCCAGGGAACCCACAAAATGCGATCGCAGTCCACCCGATCGTACACAGAAAGGAGGGTCCCCATACGCCGACGCACCTGTTCGCACGTCGTATGCATAAACGAGCCGCACGAACCAGAGAGCATAAAGAGGACCTCTAGCTCCTTTACAAAGTACAGGTTCGCCTGCCGCCGGGATGCCTTACCTAGACGCAATGACGGACGTATTCCTCTGGCCTCAACGGTTCCTGCTTCCGCTGGGATCCAAGATTGGCGGCCGAAGCCGCCTTTCCAAAGTGAGTCCTTCGTCTGTGACTAACTGTGCCAGATCGTCTTGCAAACTCCCGATCCAGTTTAACTCACCAAACTATAGCCGTACAGACCCAAATCTTAAGTCATATCACGCGACTAGCCTCTGCTCAATTTCTGTGCTCAAGGGTTTTGGTCCGCCCGAGCGGTGCAGCCGATTAGGACCATCTAATGCACTTGTTACAAGACTTCTTTTAAACACTTTCTTCCTGCCCAGTGGCGGATGATAATGGTTGTTGCCAGCCGGCGTGGAAGGTAACAGCACCGGTGCGAGCCTAATGTGCCGTCTCCACCAACACAGGGCTGTCCGGTCGTATAATAGGACTCCGCAATGGGGTTAGCAAGTGGCAGCCTAAACGATGTCGGGGACTCGCGATGTACATGCTCTGGTTCAATACATACGTGACCCGGCAGTTATCCTGCATCGGAACGTCAATCGTGCATCGGGCCAGCGTAATCGTGTCATCTGGGAGGCGGCCGTAGGATAAATAATTCAATAAAGATGTCGTTTTGCTAGTATACGCCTAGGCGTCACCCGCCATCTCTGTGCAGGTGGGCCGACGAGACACTGCCCCTGATTTCTCCGCTACTAATAGCACACACGGGGCAATACCAGCACAAGCCAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTGAGGCACGCCAATTCGCACTGTCGGGGTCGCTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTCCTTCCTTCCAGGGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCCCGTGCTTGCTCTCACGAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACCGGACAGCGCATCAACGCAACGGATCTACGTTACAGCGTGCATAGCGAGAACGGAGTTGCCGACGACGAAAGCGACACTGGGATCTGTCCGTCGTCATTCGCGGAAAGCATCCGCTCACGAGGCGGACACTGATTGACACGGTTTTGCAGAAGGTTAGGGGAATAGGTCAAATTGAGTGGCTTAAAAACGCTATGTCTGGGATTAAAGTGTAGTAAACTGCGGTCAACGGAGACGGTTTTAAGACAGGAGTTCGCAAAACCAGGCGGGGTCGCCACGACGGCTATTCCTGGTGGTTTAGGCGTACAATGTCCTGAAGAATATTTAAGAAAGAAGCACCCCTCGTCGCCTAGAATTACCTACCGCGGTCGACCATACCTTCGATTGTCGCGCCCACCCTCCCATTAGTCGGCAGAGGTGGTTGTGTTGCGATAGCCCAGCATGATATCCTAAGGCGTTACGCCGATGGATATCCCACGGAATTGCCATAGGCGCTGAACGCTACACGGACGATACGAACTTATGTATGGAGCGGGTCATCGAAAGGTCATACCCTTGTAGTTAACATGTAGCCCGGCCCTATTAGTACAGCAGTGCCTTGAGCGGCATTCTCATTATTAAGTTTTCTCTACAGCCAAACGACCAAGTGCACTTCCGCGGAGCGCGGTGGAGACTCGTCCACCCGGCAGCTCTGTAATAGGGACTAAAAGAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGAACAGAGCGGTCTTACGGCCAGTCGTATCCCTTCTCGAGTTCCGTCCGGTTAAGCGTGACACTCCCAGTGTACCCGCAAACCGTGATGGCTGTGCTTGGGGTCAATCGCATGTAGGATGGTCTCCAGACACCGGGGCACCAGTTTTCACGCCCAAAGCATAAACGACGAGCAGTCATGAGAGTCTTAGAACTGGACGTGCCGTTTCTCTGCGAACAACACCTCGAGCTGTACCGTTGTTGCGCTGCCTAGATGCAGTGCCGCTCCTATCACATTTGCCTCGACGACTGCCGCCTTCGCTGTTTCCCTAGACACTCAACAGTAAGCGCCTTTTGTAGGCAGGGGCACCCCCTGTCAGTGGCTGCGCCAAAACGTCTTCGGATCCCCTTGTCCAATCAAACTGACCGAATTCTTTCATTTAAGACCCTAATATGACATCATTAGTGACTAAATGCCACTCCCAAAATTCTGCCCAGAAGCGTTTAAGTTCGCCCCACTAAAGTTGTCTAAAACGA";
    // 5,000 random nucleotides, GC content = .45
    static STR_2: &'static [u8] = b"CTAAAGTGGCGAAATTTATGGTGTGTGACCCGTTATGCTCCATTTCGGTCAGTGGGTCATTGCTAGTAGTCGATTGCATTGTCATTCTCCGAGTGATTTAGCGTGACAGCCGCAGGGAACCCATAAAATGTAATCGTAGTCCATCTGATCGTACTTAGAAATGAAGGTCCCCTTTTACCCACGCACCTGTTTACTCGTCGTTTGCTTTTAAGAACCGCACGAACCACAGAGCATAAAGAGAACCTCTAGTTCCTTTACAAAGTACTGGTTCCCTTTTCAGCAAGATGCCTTATCTAAATGCAATGACAGACGTATTCCTCAGGCCACATCGCTTCCTACTTTCGCTGGGATCCATCATTGGCAGCTGAAACCGCCATTCCATAGTGAGTCCTTCGTCTGTGTCTTTCTGTGCCAAATCGTCTAGCAAATTGCTGATCCAGTTTATCTCACGAAATTATAGTCATACAGACCGAAATTTTAAATCAAATCACGCGACTAGGCTCAGCTTTATTTTAGTGGTCATGGGTTTTGGTCCGCCCGAGCGGTGCAACCGATTAGGACCATGTAAAACATTTGTTACAAGTCTTCTTTTAAATACAATCTTCCTGCTCAGTAGCGCATGATTATCGTTGTTGCTAGCCAGTGTGGTAAGTAACAGCACCACTGCGAGCCTAATGTGCCCTTTCCACGAACACAAGGCTATCCGATCCTATATTAGGATTCCGCAATGGGGTTAGCAAATCGCACCCTAAACGATATTGAAGACTTGCGATGTACATGCTTTGGTACAATACATACGTGTTCCAGTTGTTATCCTGTATCGGAACTTCAATTATGCATCGCACCAGCATATTCATGTCATCTAGGAAGAGCGCGTAGGATAAATAATTCAATTAAGATGTCGTTATGCTAGTATACGTCTACCCGTCACCGGCCATCTGTGTGCAGATGGGGCGACGAGTTATTGACCCTGATTTCTCCACTTCTAATACCACACACTGGGCAATACGAGCTCAAGCTAGTCTCGCAGTAACGCTCATCAGCTAACGAAAGAGTTAAAGGCTCGCTAATTCGCACTGTCAGGGTCTCTTGGGTGTTTTGCACTAGCGTCAGGTAGGCTAGTATGTGTTTTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAATGTGCAGCATACGTATTTGCTCGACGTGTTTAGTCTCTCATACTTCTCCTGGAGATCAAGGAAATGTTTCTTGTCCAAGTGGACAACGGTTCTACGGAATGGATCTACGTTACTGCCTGCATAAAGAAAACGGAGTTGCTAAGGACGAAAGCGACTTTAGGTTCTAACTGTTGACTTTGGCGGAAAAGTTTCATTCAGGAAGCAGACACTGATTGACACGGTTTAGCAGAACGTTTGAGGATTAGGTTAAATTGAGTGGTTTAATATTGGTATGTCTGGGATTAAAATATAGTATAGTGTGTTAATCGGAGACGAATTAAAGACACGAGTTCCCAAAATCAAGCGGGCTCATTACAACGGTTAATCCTGGTAGTTTACGTGAACAATGTTCTGAAGAAAATTTATGAAAAAAGGACCCGTCATCGCCTACAATTACCTACAACGGTCGACCATACCTTCGATTATCGTGGCCACTCTCGGATTACACGGCAGAGGTGGTTGTGTTCCGATAGGCCAGTATATTATTCTAAGGCGTTACCCTAATCATTTTTCATCGGATTTGCTATAGCCCTTGAACGCTACATGCACGAAACCAAATTATGTATACACTGGGTCATCAATAGGATATAGTCTTGTAGTTAACATGTAGCCCGGCCGTATTAGTACAGTAGAGCCTTCATTGACATTCTGTTTATTAAATTATTTCTACAGCAAAACGATCATATGCAAATCCACAGTGCGCGATAGAGATACATTCACTCGGCTGCTCTGTAATAGGGACTAAAAAAGTGATGATTATCATGAGTGCCCCGTTATGGTCGTGTTCGATCAGAGCGCTCTTACGAGCAGTCGTATACTTTCTCGAATTCCGTGCAGTTAAGCGTGACAGTCCCAGTGAACCCACAAAACGTGATGGCAGTCCATGCAATCATACGCAAGAAGGATGGTCTCCAGACACCGGCGCACCAGTTTTCACGCCGAAAGCATAAACGAGGAGCACAAATGAAAGTGTTTGAACTGGACCTGTAGTTTCTCTACGAAAAATACCTTGAGCTGTTGCGTTGTTGCGCTGCCTAGATGCAGTGTTGCACATATCACTTTTGCTTCAACGACTGCTGCTTTCGCTGTAACCCTAGACAGACAACAATAAGCGCTTTTTGTAGGCAAGAGCTCCGCCTATGACTAACTGCGCCAAAACATCTTCCAATCCCCTTATCCAATTTAATTCATCGAATTCTTACAATTTAGACCCTAATATCACATCATTAGACATTAATTGCCTCTGCCAAAATTCTGTCTACAAATGTTTTAGTTCGCTCCAGTAAAGTTGTTAATAACGACTACTAAATCCGCATGTTACGGGATTTCTTATTAATTCTTTTTTCGTAAGGAACAGCGGATCTTAATGGATGGCGCCAGGTGGTATGGAAGCTAATAGCGCGGGTGAGAGGGTAATTAGCCGTCTTCACCAACACAACGCTATCGGGTCATACTATAAGATTCCACAATGCGACTACTTATAAGATGTCTTAACGGTATCCGCAACTTGTGATGTGCCTACTATGCTTAAATGCATATCTCGCTCAGTAACTTTCCAATATGAGAGCATCAATTGTAGATCGGGCCGAGATAATCATGTCGTCACGGAACTTATTGTAAGAGTAATAATTTAAAAGAGATGTCAGTTTGCTGGTTCACGTAAAGGTTCCTCACACTACCTCTAAATAAGTGAGCGGTCGTGACATTATCCCTGATTTTCTCACTACTATTAGTACTCACGACACAATTCTACCACAGCCTTGTTTCGCCAGAATGCCAGTCAGCATAAAGAAGAGCTCAAGGCAGGTCAACTCGCATTGTGAGAGTTACATGAACGTTCGGCACTACCGACACGAACCTCAGTTAGCGTACATCCTACCAGAGGTCTGTGGCCCCGTGGTCAAAAGTGCGGATTTCGTATTTGCTGCTCGTCAGTACTTTCAGAATCATGACCTGCACGGTAAAAAGACGCTTATTATGGAGTTCGACATGGCAATAACGCGACGAATCTACGTCATGACGAGAATAGTATAAACAAAACTGCTGACGGCAGAAGCGTCAAAGAAGTCTGTAAATTGTTATTCGCGAAAAACATCCGTCTCCGTGGGGGATAATCACCGACGCCATTTTATAGAAGCCTAGGGGAACAGATTGGTTTAATTAGCTTAAGAAAGTAAATTCTGGGATTATACTGTAGTAATCACTAATTTACGGTGAGGGTTTTATGGCGGATTTTTACAAATTCAAACCAGGTGATTTCAACAAATTTTGTTGACGATTTAGGCGCACTATCCCCTAAACTACAAATTAAAAAATAGCGTTCCTTGACGGCTAGAATTACTTACCGGCCTTCACCATACCTTCGATATTCGCGCCCACTCTCCCATTAATCCGTACAAGTGGATGTAATGCGATTGTCCGCTAAGATATTCTAACGTGTAACGTAGATAAGTATTTTACAGAGTTGCCGTACGCGTTGAACACTTCACAGATGATAGGAATTTGCGTATAGAGCGTGTTATTGAGGAGTTATACACCCGTAGACTACAATGGGCCCAACTCAATCAGAACTCGAGTGCCTTGAATAACATACTCATCACTAAACATTCTCAACAATCAATCGAGCAAGTCCATTATCAACGAGTGTGTTGCAGTTTTATTCTCTTGCCAGCATTGTAATAGGCACTAAAAGAATGATGATAGTCATGAGTACTGAGCTAAGACGGCGTCGATGCATAGCGGACTTTCGGTCAATCACAATTCCTCACGAGACTCGTCCTGTTGAGCGTATCACTCTCAATGTACAAGCAACCCAAGAAGGCTGTGCCTGGACTCAACTGGATGCAGGATGAACTCCAGACACGGGGTCACTACTCTTCATACATAAAGCAAGAACGTCGAACAGTCATGAAAGTCTTAGTACCGCACGTACCATCTTACTGTGAATATTGCTTGAAGCTGTACCGTTATTGGGGGGCAAAGATGAAGTTCTCTTCTTTTCATAATTGTACTGACGACAGTCGTGTTCTCGGTTTCTTCAAAGGTTAAAGAATAAAGGCTTATTGTAGGCAGAGGAACGCCCTTTTAGTGGCTGGCGTTAAGTATCTTCGGACCCCCTTGTCTATCCAGATTAATCGAATTCTCTCATTTAGGACCTTAGTAAGTCATCATTGGTATTTGAATGCGACCTTGAAGAAACCGCTTAAAAATGTCAATGGTTGATCCACTAAACTTCATTTAATTAACTCCTAAATCAGCGCGATAGGCTATTAGAGGTTTAATTTTGTATAGCAAGGTACTTCCGATCTTAATGAATGGCCGGAAAAGGTACGGACGCGATATGCGAGGGTGAAAGGGCAAATAGACAGGTTCGTCTTTGTCACGCTAGGAGGCAATTCTATAAGAATGCATATTGCATCGATACATAAAATGTCTCGATCGCATGCGCAATTTGTGAAGTGTCTATTATCCCTAAGCCCATTTCCCGCATAATAACCCCTGATTGTATCCGCATTTGATGCTACCCAGGTTGAGTTAGCGTCGAGCTCGCGGAACTTATTGCATGAGTAGAGTTGAGTAAGAGCTGTTAGATGGCTCGCTGAACTAATAGTTGTCCACAGAACGTCAAGATTAGAAAACGGTTGTAGCATTATCGGAGGTTCTCTAACTACTATCAATACCCGTGTCTTGACTCTGCTGCGGCTACCTATCGCCTGAAAACCAGTTGGTGTTAAGGGATGCTCTGTCCAGGACGCCACATGTAGTGAAACTTACATGTTCGTTGGGTTCACCCGACT";


    #[test]
    fn test_semiglobal() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_local() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_global() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations, [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_issue11() {
        let y = b"TACC";//GTGGAC";
        let x = b"AAAAACC";//GTTGACGCAA";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations, [Ins, Ins, Ins, Subst, Match, Match, Match]);
    }


    #[test]
    fn test_issue12_1() {
        let x = b"CCGGCA";
        let y = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 1);
        assert_eq!(alignment.operations, [Match, Match, Match, Subst, Subst, Subst]);
    }

    #[test]
    fn test_issue12_2() {
        let y = b"CCGGCA";
        let x = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.operations, [Ins, Ins, Match, Ins, Ins, Subst, Match, Ins, Ins, Match, Match]);
    }

    #[test]
    fn test_fn_semiglobal() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let alignment = align_semiglobal(x, y, -5, -1, &score);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_fn_local() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let alignment = align_local(x, y, -5, -1, &score);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_fn_global() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let alignment = align_global(x, y, -5, -1, &score);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations, [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_fn_semiglobal_score() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        assert_eq!(score_semiglobal(x, y, -5, -1, &score), align_semiglobal(x, y, -5, -1, &score).score);
    }

    #[test]
    fn test_fn_local_score() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        assert_eq!(score_local(x, y, -5, -1, &score), align_local(x, y, -5, -1, &score).score);
    }

    #[test]
    fn test_fn_global_score() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        assert_eq!(score_global(x, y, -5, -1, &score), align_global(x, y, -5, -1, &score).score);
    }

    #[bench]
    fn bench_fn_semiglobal_score(b: &mut Bencher) {
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        b.iter(|| { score_semiglobal(STR_1, STR_2, -5, -1, &score) });
    }

    #[bench]
    fn bench_fn_local_score(b: &mut Bencher) {
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        b.iter(|| { score_local(STR_1, STR_2, -5, -1, &score) });
    }

    #[bench]
    fn bench_fn_global_score(b: &mut Bencher) {
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        b.iter(|| { score_global(STR_1, STR_2, -5, -1, &score) });
    }
}
