// Copyright 2020 Tianyi Shi
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(non_snake_case)]
use crate::alignment::pairwise::MIN_SCORE;
use crate::alignment::pairwise::{MatchFunc, Scoring};
use crate::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use crate::utils::TextSlice;
use std::cmp::max;

pub struct Aligner<F: MatchFunc> {
    scoring: Scoring<F>,
}

impl<F: MatchFunc> Aligner<F> {
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, gap_extend, match_fn),
        }
    }

    /// Fast global alignment with $O(nm)$ space. Used when `y.len()` is small.
    ///
    /// # Implementation Details
    ///
    /// ## Traceback Matrix `T: Vec<u8>`
    ///
    /// ```ignore
    /// 0b00   0b01    0b10    0b11
    /// start  insert  delete  match_or_subst
    /// ```
    ///
    /// ```ignore
    /// 0b00101101
    ///     / |  \
    ///    S  D   I
    /// ```
    pub fn global(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut S = vec![0; m]; //                                            ! 32 * m bits
        let mut D = vec![MIN_SCORE; m]; //                                    ! 32 * m bits
        let mut T: Vec<TracebackCell> = vec![TracebackCell::new(); n * m]; // ! 8 * n * m bits
        let mut s: i32; // S[i - 1][j - 1] or S[i][j - 1]
        let mut c: i32; // S[i - 1][j] or S[i][j]
        let mut e: i32; // I[i - 1] or I [i]
        let mut idx: usize = 0; // j * m + i
        let mut p: &u8; // x[i] (x[i - 1]) // TODO: is not using p/q faster?
        let mut q: &u8; // y[j] (y[j - 1])
        let mut S_i: &mut i32; // S[i]
        let mut D_i: &mut i32; // D[i]
        let mut score_1: i32; // used when determining D[i] and I[i]
        let mut score_2: i32; // used when determining D[i] and I[i]

        // SAFETY: unchecked indexing is used here. x, y, S, D, T all have a fixed size related to n and/or m;
        //         it should have been implied by the for loops that all indexing operations are in-bound but
        //         the compiler wasn't smart enough to notice this as of October 2020.
        unsafe {
            // T[0] = TracebackCell::new() // origin at T[0 * n + 0]
            let mut t = self.scoring.gap_open;
            for i in 1..m {
                t += self.scoring.gap_extend;
                // I[0][j] = t will not be read
                *S.get_unchecked_mut(i) = t;
                let mut tb = TracebackCell::new();
                tb.set_s_bits(TB_UP);
                idx += 1;
                *T.get_unchecked_mut(idx) = tb; // T[0 * m + i]
            }

            t = self.scoring.gap_open;
            for j in 1..n {
                s = *S.get_unchecked(0);
                t += self.scoring.gap_extend;
                c = t;
                *S.get_unchecked_mut(0) = c;
                e = MIN_SCORE;
                // D[0] = t will not be read
                let mut tb = TracebackCell::new();
                tb.set_s_bits(TB_LEFT);
                idx += 1;
                *T.get_unchecked_mut(idx) = tb; // T[j * m + 0]

                q = y.get_unchecked(j - 1);
                for i in 1..m {
                    S_i = S.get_unchecked_mut(i);
                    D_i = D.get_unchecked_mut(i);
                    p = x.get_unchecked(i - 1);
                    let mut tb = TracebackCell::new();

                    score_1 = e + self.scoring.gap_extend;
                    score_2 = c + self.scoring.gap_open + self.scoring.gap_extend;
                    e = if score_1 > score_2 {
                        tb.set_i_bits(TB_UP);
                        score_1
                    } else {
                        tb.set_i_bits(T.get_unchecked(idx).get_s_bits()); // T[i-1][j]
                        score_2
                    };

                    idx += 1; // ! update idx

                    score_1 = *D_i + self.scoring.gap_extend;
                    score_2 = *S_i + self.scoring.gap_open + self.scoring.gap_extend;
                    *D_i = if score_1 > score_2 {
                        tb.set_d_bits(TB_LEFT);
                        score_1
                    } else {
                        tb.set_d_bits(T.get_unchecked(idx - m).get_s_bits()); //T[i][j-1]
                        score_2
                    };

                    // c is becoming S[i][j]
                    c = s + self.scoring.match_fn.score(*p, *q);
                    tb.set_s_bits(TB_DIAG); // no need to be exact at this stage

                    if e > c {
                        c = e;
                        tb.set_s_bits(TB_UP);
                    }

                    if *D_i > c {
                        c = *D_i;
                        tb.set_s_bits(TB_LEFT);
                    }

                    s = *S_i;
                    *S_i = c;
                    *T.get_unchecked_mut(idx) = tb;
                }
            }

            let (mut operations, i, j) = Self::traceback(&T, x, y, m - 1, n - 1, m);
            operations.resize(operations.len() + i, AlignmentOperation::Ins); // reaching at (i, 0)
            operations.resize(operations.len() + j, AlignmentOperation::Del); // reaching at (0, j)

            operations.reverse();
            Alignment {
                score: *S.get_unchecked(m - 1),
                xstart: 0,
                ystart: 0,
                xend: m - 1,
                yend: n - 1,
                xlen: m - 1,
                ylen: n - 1,
                operations,
                mode: AlignmentMode::Global,
            }
        }
    }

    pub fn local(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut S = vec![0; m]; //                                            ! 32 * m bits
        let mut D = vec![MIN_SCORE; m]; //                                    ! 32 * m bits
        let mut T: Vec<TracebackCell> = vec![TracebackCell::new(); n * m]; // ! 8 * n * m bits
        let mut s: i32;
        let mut c: i32;
        let mut e: i32;
        let mut idx: usize = m - 1;
        let mut p: &u8;
        let mut q: &u8;
        let mut S_i: &mut i32;
        let mut D_i: &mut i32;
        let mut score_1: i32;
        let mut score_2: i32;
        let mut max_score: i32 = MIN_SCORE;
        let mut max_coords: [usize; 2] = [0, 0];

        unsafe {
            // all cells in the first column can be the origin
            for j in 1..n {
                s = 0;
                c = 0;
                *S.get_unchecked_mut(0) = 0;
                e = MIN_SCORE;
                idx += 1;
                *T.get_unchecked_mut(idx) = TracebackCell::new(); // T[j * m + 0] // all cells in the first row can be the origin
                q = y.get_unchecked(j - 1);
                for i in 1..m {
                    S_i = S.get_unchecked_mut(i);
                    D_i = D.get_unchecked_mut(i);
                    p = x.get_unchecked(i - 1);
                    let mut tb = TracebackCell::new();

                    score_1 = e + self.scoring.gap_extend;
                    score_2 = c + self.scoring.gap_open + self.scoring.gap_extend;
                    e = if score_1 > score_2 {
                        tb.set_i_bits(TB_UP);
                        score_1
                    } else {
                        tb.set_i_bits(T.get_unchecked(idx).get_s_bits()); // T[i-1][j]
                        score_2
                    };

                    idx += 1; // ! update idx

                    score_1 = *D_i + self.scoring.gap_extend;
                    score_2 = *S_i + self.scoring.gap_open + self.scoring.gap_extend;
                    *D_i = if score_1 > score_2 {
                        tb.set_d_bits(TB_LEFT);
                        score_1
                    } else {
                        tb.set_d_bits(T.get_unchecked(idx - m).get_s_bits()); //T[i][j-1]
                        score_2
                    };

                    c = s + self.scoring.match_fn.score(*p, *q);
                    tb.set_s_bits(TB_DIAG); // no need to be exact at this stage

                    if e > c {
                        c = e;
                        tb.set_s_bits(TB_UP);
                    }

                    if *D_i > c {
                        c = *D_i;
                        tb.set_s_bits(TB_LEFT);
                    }

                    if c < 0 {
                        c = 0;
                        tb = TracebackCell::new(); // reset origin
                    }

                    if c >= max_score {
                        max_score = c;
                        max_coords = [i, j];
                    }

                    s = *S_i;
                    *S_i = c;
                    *T.get_unchecked_mut(idx) = tb;
                }
            }

            let (xend, yend) = (max_coords[0], max_coords[1]);

            let (mut operations, xstart, ystart) = Self::traceback(&T, x, y, xend, yend, m);

            operations.reverse();
            Alignment {
                score: max_score,
                xstart,
                ystart,
                xend,
                yend,
                xlen: xend - xstart,
                ylen: yend - ystart,
                operations,
                mode: AlignmentMode::Local,
            }
        }
    }

    pub fn traceback(
        T: &Vec<TracebackCell>,
        x: TextSlice,
        y: TextSlice,
        mut i: usize,
        mut j: usize,
        m: usize,
    ) -> (Vec<AlignmentOperation>, usize, usize) {
        let mut operations = Vec::with_capacity(m);
        unsafe {
            let mut next_layer = T.get_unchecked(j * m + i).get_s_bits(); // start from the last tb cell
            loop {
                match next_layer {
                    TB_START => break,
                    TB_UP => {
                        operations.push(AlignmentOperation::Ins);
                        next_layer = T.get_unchecked(j * m + i).get_i_bits();
                        i -= 1;
                    }
                    TB_LEFT => {
                        operations.push(AlignmentOperation::Del);
                        next_layer = T.get_unchecked(j * m + i).get_d_bits();
                        j -= 1;
                    }
                    TB_DIAG => {
                        i -= 1;
                        j -= 1;
                        next_layer = T.get_unchecked(j * m + i).get_s_bits(); // T[i - 1][j - 1]
                        operations.push(if *y.get_unchecked(j) == *x.get_unchecked(i) {
                            AlignmentOperation::Match
                        } else {
                            AlignmentOperation::Subst
                        });
                    }
                    _ => unreachable!(),
                }
            }
        }
        (operations, i, j)
    }
}

#[derive(Copy, Clone)]
pub struct TracebackCell(u8);

const TB_START: u8 = 0b00;
const TB_UP: u8 = 0b01;
const TB_LEFT: u8 = 0b10;
const TB_DIAG: u8 = 0b11;

// Traceback bit positions (LSB)
const I_POS: u8 = 0; // Meaning bits 0,1 corresponds to I and so on
const D_POS: u8 = 2;
const S_POS: u8 = 4;

impl TracebackCell {
    /// Initialize a blank traceback cell
    #[inline(always)]
    pub fn new() -> TracebackCell {
        TracebackCell(0u8)
    }

    /// Sets 2 bits [pos, pos+2) with the 2 LSBs of value
    #[inline(always)]
    fn set_bits(&mut self, pos: u8, value: u8) {
        let bits: u8 = (0b11) << pos;
        self.0 = (self.0 & !bits) // First clear the bits
            | (value << pos) // And set the bits
    }

    #[inline(always)]
    pub fn set_i_bits(&mut self, value: u8) {
        // Traceback corresponding to matrix I
        self.set_bits(I_POS, value);
    }

    #[inline(always)]
    pub fn set_d_bits(&mut self, value: u8) {
        // Traceback corresponding to matrix D
        self.set_bits(D_POS, value);
    }

    #[inline(always)]
    pub fn set_s_bits(&mut self, value: u8) {
        // Traceback corresponding to matrix S
        self.set_bits(S_POS, value);
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_bits(self, pos: u8) -> u8 {
        (self.0 >> pos) & (0b11)
    }

    #[inline(always)]
    pub fn get_i_bits(self) -> u8 {
        self.get_bits(I_POS)
    }

    #[inline(always)]
    pub fn get_d_bits(self) -> u8 {
        self.get_bits(D_POS)
    }

    #[inline(always)]
    pub fn get_s_bits(self) -> u8 {
        self.get_bits(S_POS)
    }
}

impl std::fmt::Debug for TracebackCell {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct(&format!("{:06b}", self.0)).finish()
    }
}

// adapted from pariwise/mod.rs
#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::AlignmentOperation::*;
    use crate::scores::blosum62;
    use std::iter::repeat;

    fn equivalent_operations(o1: &[AlignmentOperation], o2: &[AlignmentOperation]) -> bool {
        if o1.len() != o2.len() {
            return false;
        }
        let (mut i, mut d, mut s, mut m, mut x, mut y) =
            (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
        for o in o1 {
            match o {
                Ins => i += 1,
                Del => d += 1,
                Subst => s += 1,
                Match => m += 1,
                Xclip(n) => x += n,
                Yclip(n) => y += n,
            }
        }
        for o in o2 {
            match o {
                Ins => i -= 1,
                Del => d -= 1,
                Subst => s -= 1,
                Match => m -= 1,
                Xclip(n) => x -= n,
                Yclip(n) => y -= n,
            }
        }
        i == 0 && d == 0 && s == 0 && m == 0 && x == 0 && y == 0
    }

    #[test]
    fn test_global_affine_ins() {
        let x = b"ACGAGAACA";
        let y = b"ACGACA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);

        println!("aln:\n{}", alignment.pretty(x, y));
        assert!(equivalent_operations(
            &alignment.operations,
            &[Match, Match, Match, Ins, Ins, Ins, Match, Match, Match]
        ));
    }

    #[test]
    fn test_global_affine_ins2() {
        let x = b"AGATAGATAGATAGGGAGTTGTGTAGATGATCCACAGT";
        let y = b"AGATAGATAGATGTAGATGATCCACAGT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);

        println!("aln:\n{}", alignment.pretty(x, y));

        let mut correct = Vec::new();
        correct.extend(repeat(Match).take(11));
        correct.extend(repeat(Ins).take(10));
        correct.extend(repeat(Match).take(17));

        assert!(equivalent_operations(&alignment.operations, &correct));
    }

    #[test]
    fn test_global() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);

        println!("\naln:\n{}", alignment.pretty(x, y));
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        // assert_eq!(
        //     alignment.operations,
        //     [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        // );
        assert!(equivalent_operations(
            &alignment.operations,
            &[Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        ))
    }

    #[test]
    fn test_local_affine_ins2() {
        let x = b"ACGTATCATAGATAGATAGGGTTGTGTAGATGATCCACAG";
        let y = b"CGTATCATAGATAGATGTAGATGATCCACAGT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.xstart, 1);
        assert_eq!(alignment.ystart, 0);
    }

    #[test]
    fn test_local() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert!(equivalent_operations(
            &alignment.operations,
            &[Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        ));
    }
}
