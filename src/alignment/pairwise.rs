use std::i32;
use std::iter::repeat;
use std::cmp::max;
use std::collections::Bitv;

use alignment::{Alignment, AlignmentOperation};

#[derive(Copy)]
enum AlignmentType {
    Global,
    Semiglobal,
    Local
}


/// Calculate alignments with a generalized variant of the Smith Waterman algorithm.
/// Complexity: O(n * m) for strings of length m and n.
///
/// # Example
///
/// ```
/// use bio::alignment::pairwise::Aligner;
/// use bio::alignment::AlignmentOperation::{Match, Subst, Ins, Del};
/// let x = b"ACCGTGGAT";
/// let y = b"AAAAACCGTTGAT";
/// let score = |&: a: u8, b: u8| if a == b {1i32} else {-1i32};
/// let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
/// let alignment = aligner.semiglobal(x, y);
/// assert_eq!(alignment.i, 4);
/// assert_eq!(alignment.j, 0);
/// assert_eq!(alignment.operations, [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
/// ```
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

    fn init(&mut self, m: usize, alignment_type: AlignmentType) {
        // set minimum score to -inf, and allow to add gap_extend
        // without overflow
        let min_score = i32::MIN - self.gap_extend;
        for k in 0..2us {
            self.S[k].clear();
            self.I[k].clear();
            self.D[k].clear();
            self.I[k].extend(repeat(min_score).take(m + 1));
            self.D[k].extend(repeat(min_score).take(m + 1));
            match alignment_type {
                AlignmentType::Global => {
                    let ref mut s = self.S[k];
                    let mut score = self.gap_open;
                    for _ in 0..m+1 {
                        s.push(score);
                        score += self.gap_extend;
                    }
                },
                _ => self.S[k].extend(repeat(0).take(m + 1))
            }
        }
    }

    /// Calculate global alignment.
    pub fn global(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        self.align(x, y, AlignmentType::Global)
    }

    /// Calculate semiglobal alignment.
    pub fn semiglobal(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        self.align(x, y, AlignmentType::Semiglobal)
    }

    /// Calculate local alignment.
    pub fn local(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        self.align(x, y, AlignmentType::Local)
    }

    fn align(&mut self, x: &[u8], y: &[u8], alignment_type: AlignmentType) -> Alignment {
        let (m, n) = (x.len(), y.len());

        self.init(m, alignment_type);
        self.traceback.init(n, alignment_type);

        let (mut best, mut best_i, mut best_j) = (0, 0, 0);
        let mut col = 0;

        for i in 1..n+1 {
            col = i % 2;
            let prev = 1 - col;

            match alignment_type {
                AlignmentType::Global => {
                    self.S[col][0] = self.gap_open + (i as i32 - 1) * self.gap_extend;
                    self.traceback.del(i, 0);
                },
                _ => {
                    // with local and semiglobal, allow to begin anywhere in y
                    self.S[col][0] = 0;
                }
            }

            let b = y[i - 1];
            let mut score = 0i32;
            for j in 1..m+1 {
                let a = x[j - 1];

                // score for deletion
                let d_score = max(
                    self.S[prev][j] + self.gap_open,
                    self.D[prev][j] + self.gap_extend
                );
                // score for insertion
                let i_score = max(
                    self.S[col][j-1] + self.gap_open,
                    self.I[col][j-1] + self.gap_extend
                );
                // score for substitution
                score = self.S[prev][j-1] + (self.score)(a, b);

                if d_score > score {
                    score = d_score;
                    self.traceback.del(i, j);
                }
                else if i_score > score {
                    score = i_score;
                    self.traceback.ins(i, j);
                }
                else {
                    self.traceback.subst(i, j);
                }

                match alignment_type {
                    AlignmentType::Local => {
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
                    _ => ()
                }

                self.S[col][j] = score;
                self.D[col][j] = d_score;
                self.I[col][j] = i_score;
            }

            match alignment_type {
                AlignmentType::Semiglobal => {
                    if score > best {
                        best = score;
                        best_i = i;
                        best_j = m;
                    }
                },
                _ => ()
            }
        }
        match alignment_type {
            AlignmentType::Global => {
                let score = self.S[col][m];
                self.traceback.get_alignment(n, m, x, y, score)
            },
            _ => self.traceback.get_alignment(best_i, best_j, x, y, best)
        }

        
    }
}


struct Traceback {
    subst: Vec<Bitv>,
    del: Vec<Bitv>,
    ins: Vec<Bitv>
}


impl Traceback {
    fn with_capacity(m: usize, n: usize) -> Self {
        let get_vec = |&:| repeat(Bitv::from_elem(m + 1, false)).take(n + 1).collect::<Vec<Bitv>>();
        Traceback {
            subst: get_vec(),
            del: get_vec(),
            ins: get_vec()
        }
    }

    fn init(&mut self, n: usize, alignment_type: AlignmentType) {
        match alignment_type {
            AlignmentType::Global => {
                for i in 0..n+1 {
                    self.subst[i].set_all();
                    self.subst[i].negate();
                    self.del[i].set_all();
                    self.ins[i].set_all();
                    self.ins[i].negate();
                }
                // set the first cell to start, the rest to deletions
                self.del[0].set(0, false);
            },
            _ => {
                for i in 0..n+1 {
                    self.subst[i].set_all();
                    self.subst[i].negate();
                    self.del[i].set_all();
                    self.del[i].negate();
                    self.ins[i].set_all();
                    self.ins[i].negate();
                }
            }
        }
    }

    fn start(&mut self, i: usize, j: usize) {
        self.subst[i].set(j, false);
        self.del[i].set(j, false);
        self.ins[i].set(j, false);
    }

    fn subst(&mut self, i: usize, j: usize) {
        self.subst[i].set(j, true);
    }

    fn del(&mut self, i: usize, j: usize) {
        self.del[i].set(j, true);
    }

    fn ins(&mut self, i: usize, j: usize) {
        self.ins[i].set(j, true);
    }

    fn is_subst(&self, i: usize, j: usize) -> bool {
        self.subst[i].get(j).unwrap()
    }

    fn is_del(&self, i: usize, j: usize) -> bool {
        self.del[i].get(j).unwrap()
    }

    fn is_ins(&self, i: usize, j: usize) -> bool {
        self.ins[i].get(j).unwrap()
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
