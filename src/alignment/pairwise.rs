//! Calculate alignments with a generalized variant of the Smith Waterman algorithm.
//! Complexity: O(n * m) for strings of length m and n.
//!
//! # Example
//!
//! ```
//! use bio::alignment::pairwise::Aligner;
//! use bio::alignment::AlignmentOperation::{Match, Subst, Ins, Del};
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! let score = |&: a: u8, b: u8| if a == b {1i32} else {-1i32};
//! let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
//! let alignment = aligner.semiglobal(x, y);
//! assert_eq!(alignment.i, 4);
//! assert_eq!(alignment.j, 0);
//! assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
//! ```


use std::i32;
use std::iter::repeat;
use std::cmp::max;

use alignment::{Alignment, AlignmentOperation};
use bitencoding::BitEnc;


struct AlignmentState {
    m: usize,
    n: usize,
    best: i32,
}


macro_rules! align {
    (
        $obj:ident, $state:ident,
        $init:block, $inner:block, $outer:block, $ret:block
    ) => (
        {
            let (mut best, mut best_i, mut best_j) = (0, 0, 0);
            let mut col = 0;

            for i in 1..$n+1 {
                col = i % 2;
                let prev = 1 - col;

                // init code
                $init

                let b = y[i - 1];
                let mut score = 0i32;
                for j in 1..$m+1 {
                    let a = x[j - 1];

                    // score for deletion
                    let d_score = max(
                        $obj.S[prev][j] + $obj.gap_open,
                        $obj.D[prev][j] + $obj.gap_extend
                    );
                    // score for insertion
                    let i_score = max(
                        $obj.S[col][j-1] + $obj.gap_open,
                        $obj.I[col][j-1] + $obj.gap_extend
                    );
                    // score for substitution
                    score = $obj.S[prev][j-1] + ($obj.score)(a, b);

                    if d_score > score {
                        score = d_score;
                        $obj.traceback.del(i, j);
                    }
                    else if i_score > score {
                        score = i_score;
                        $obj.traceback.ins(i, j);
                    }
                    else {
                        $obj.traceback.subst(i, j);
                    }

                    // inner code
                    $inner

                    $obj.S[col][j] = score;
                    $obj.D[col][j] = d_score;
                    $obj.I[col][j] = i_score;
                }

                // outer code
                $outer
            }

            // return code
            $ret
        }
    );
}


#[allow(non_snake_case)]
pub struct Aligner<F> where F: Fn(u8, u8) -> i32 {
    S: [Vec<i32>; 2],
    I: [Vec<i32>; 2],
    D: [Vec<i32>; 2],
    traceback: Traceback,
    gap_open: i32,
    gap_extend: i32,
    score: F
}


impl<F> Aligner<F> where F: Fn(u8, u8) -> i32 {
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
    pub fn with_capacity(m: usize, n: usize, gap_open: i32, gap_extend: i32, score: F) -> Self {
        let get_vec = |&:| Vec::with_capacity(m + 1);
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

    fn init(&mut self, m: usize, global: bool) {
        // set minimum score to -inf, and allow to add gap_extend
        // without overflow
        let min_score = i32::MIN - self.gap_extend;
        for k in 0..2us {
            self.S[k].clear();
            self.I[k].clear();
            self.D[k].clear();
            self.I[k].extend(repeat(min_score).take(m + 1));
            self.D[k].extend(repeat(min_score).take(m + 1));
            if global {
                let ref mut s = self.S[k];
                let mut score = self.gap_open;
                for _ in 0..m+1 {
                    s.push(score);
                    score += self.gap_extend;
                }
            }
            else {
                self.S[k].extend(repeat(0).take(m + 1))
            }
        }
    }

    /// Calculate global alignment.
    pub fn global(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, true);
        self.traceback.init(m, n, true);

        align!(
            self, m, n,
            {
                self.S[col][0] = self.gap_open + (i as i32 - 1) * self.gap_extend;
                self.traceback.del(i, 0);
            }, {}, {},
            {
                let score = self.S[col][m];
                self.traceback.get_alignment(n, m, x, y, score)
            }
        );
    }

    /// Calculate semiglobal alignment.
    pub fn semiglobal(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, false);
        self.traceback.init(m, n, false);

        align!(
            self, m, n,
            { self.S[col][0] = 0; },
            {},
            {
                if score > best {
                    best = score;
                    best_i = i;
                    best_j = m;
                }
            },
            { self.traceback.get_alignment(best_i, best_j, x, y, best) }
        );
    }

    /// Calculate local alignment.
    pub fn local(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        let (m, n) = (x.len(), y.len());
        self.init(m, false);
        self.traceback.init(m, n, false);

        align!(
            self, m, n,
            { self.S[col][0] = 0; },
            {
                if score < 0 {
                    self.traceback.start(i, j);
                    score = 0;
                }
                else if score > best {
                    best = score;
                    best_i = i;
                    best_j = j;
                }
            },
            {},
            { self.traceback.get_alignment(best_i, best_j, x, y, best) }
        );
    }
}


struct Traceback {
    matrix: Vec<BitEnc>
}


static TBSTART: u8 = 0b00;
static TBSUBST: u8 = 0b01;
static TBINS: u8   = 0b10;
static TBDEL: u8   = 0b11;


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

    fn init(&mut self, m: usize, n: usize, global: bool) {
        if global {
            // set the first cell to start, the rest to deletions
            for i in 0..n+1 {
                self.matrix[i].clear();
                self.matrix[i].push_values(m + 1, TBDEL);
            }
            self.matrix[0].set(0, TBSTART);
        }
        else {
            for i in 0..n+1 {
                self.matrix[i].clear();
                self.matrix[i].push_values(m + 1, TBSTART);
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

    fn is_subst(&self, i: usize, j: usize) -> bool {
        self.matrix[i].get(j).unwrap() == TBSUBST
    }

    fn is_del(&self, i: usize, j: usize) -> bool {
        self.matrix[i].get(j).unwrap() == TBDEL
    }

    fn is_ins(&self, i: usize, j: usize) -> bool {
        self.matrix[i].get(j).unwrap() == TBINS
    }

    fn get_alignment(&self, mut i: usize, mut j: usize, x: &[u8], y: &[u8], score: i32) -> Alignment {
        let mut ops = Vec::with_capacity(x.len());

        loop {
            let (ii, jj, op) = if self.is_subst(i, j) {
                let op = if y[i-1] == x[j-1] {
                    AlignmentOperation::Match
                }
                else {
                    AlignmentOperation::Subst
                };
                (i - 1, j - 1, op)
            }
            else if self.is_del(i, j) {
                (i - 1, j, AlignmentOperation::Del)
            }
            else if self.is_ins(i, j) {
                (i, j - 1, AlignmentOperation::Ins)
            } else {
                // reached alignment start
                break;
            };
            ops.push(op);
            i = ii;
            j = jj;
        }

        ops.reverse();
        Alignment { i: i, j: j, operations: ops, score: score}
    }
}

#[cfg(test)]
mod tests {
    use super::Aligner;
    use alignment::AlignmentOperation::{Match, Subst, Del};

    #[test]
    fn test_semiglobal() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |&: a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.i, 4);
        assert_eq!(alignment.j, 0);
        assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_local() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |&: a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.i, 4);
        assert_eq!(alignment.j, 0);
        assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_global() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |&: a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.i, 0);
        assert_eq!(alignment.j, 0);
        assert_eq!(alignment.operations, [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }
}
