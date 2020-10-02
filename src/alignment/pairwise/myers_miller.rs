// Copyright 2020 Tianyi Shi
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Alignment with affine gap penalty in linear space, by combining Gotoh's (1982) and
//! Hirschberg's (1975) ideas, which was first implemented in C (Myers & Miller 1988).
//!
//! Myers & Miller originally used their technique to implement global alignment only,
//! but alignments of other modes can be achieved by first finding the termini of
//! the non-global alignment and then global-aligning the corresponding substrings.
//!
//! # Time Complexity
//!
//! $O(nm)$ for strings of length $m$ and $n$.
//!
//! # Space Complexity
//!
//! $O(n)$. For exact number of bits, see `global`, `semiglobal`, `local`, and `custom` methods
//!
//! # References
//!
//! - [Eugene W. Myers and Webb Miller (1988) Optimal alignments in linear space. _Bioinformatics_ **4**: 11-17.](https://doi.org/10.1093/bioinformatics/4.1.11)
//! - [Hirschberg, D. S. (1975) A linear space algorithm for computing maximal common subsequences. _Commun. Assoc. Comput. Mach._ **18**: 341-343.](https://doi.org/10.1145/360825.360861)
//! - [Gotoh, O. (1982) An improved algorithm for matching biological sequences. _J. Molec. Biol._ **162**: 705-708.](https://doi.org/10.1016/0022-2836(82)90398-9)

#![allow(non_snake_case)]
use crate::alignment::pairwise::MIN_SCORE;
use crate::alignment::pairwise::{MatchFunc, Scoring};
use crate::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use crate::utils::TextSlice;
use std::cmp::max;

pub struct Aligner<F: MatchFunc + Sync> {
    scoring: Scoring<F>,
    parallel: bool,
}

impl<F: MatchFunc + Sync> Aligner<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions (also see bio::scores)
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, gap_extend, match_fn),
            parallel: true, //cfg!(feature = "parallel"),
        }
    }

    /// Run without parallelization
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use bio::alignment::pairwise::myers_miller::Aligner;
    ///
    /// let aligner = Aligner.new(-5, -1, |a,b|{a==b{1i32}else{-1i32}}).no_parallel();
    /// let result = aligner.global()
    /// ```
    #[cfg(feature = "parallel")]
    pub fn no_parallel(mut self) -> Self {
        self.parallel = false;
        self
    }

    /// Calculate global alignment of `x` against `y`.
    ///
    /// - Time complexity: $O(mn)$, where $m$ and $n$ are the lengths of the first and the second
    ///   sequence
    /// - Space complexity: $O(n)$; specifically, about $64n$ bits. Note that `compute_recursive()`
    ///   uses less and less space as recursion proceeds.
    pub fn global(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        let (m, n) = (x.len(), y.len());
        let operations =
            self.compute_recursive(x, y, m, n, self.scoring.gap_open, self.scoring.gap_open);
        let score = self.cost_only(x, y, false, self.scoring.gap_open).0[x.len()];
        return Alignment {
            score,
            xstart: 0,
            ystart: 0,
            xend: m,
            yend: n,
            xlen: m,
            ylen: n,
            operations,
            mode: AlignmentMode::Global,
        };
    }

    /// Calculate local alignment of x against y.
    ///
    /// This algorithm first find the pair of coordinates corresponding to the termini of one optimal
    /// local alignment path. Then, a global alignment on the corresponding substrings calculates the
    /// operations.
    ///
    /// - Space complexity: $O(n)$
    /// - Time complexity: $(mn + m'n')$ where $m$ and $n$ are lengths of the input sequences; $m'$
    ///   and $n'$ are the lengths of the substrings of the input sequences that correspond to the
    ///   optimal alignment.
    ///
    /// Termini can be determined by two approaches. Huang's method (1991) is faster, but uses about
    /// $446m$ bits. Shamir's method is slightly slower, but only uses $64m$ bits.
    ///
    /// # References
    ///
    /// - [Huang, X. and Miller, W. 1991. A time-efficient linear-space local similarity algorithm. Adv. Appl. Math. 12, 3 (Sep. 1991), 337-357](https://doi.org/10.1016/0196-8858(91)90017-D)
    pub fn local(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        let (score, xstart, ystart, xend, yend) = self.find_local_score_and_termini_huang(x, y);
        let xlen = xend - xstart;
        let ylen = yend - ystart;
        let operations = self.compute_recursive(
            &x[xstart..xend],
            &y[ystart..yend],
            xlen,
            ylen,
            self.scoring.gap_open,
            self.scoring.gap_open,
        );
        Alignment {
            score,
            xstart,
            xend,
            ystart,
            yend,
            xlen,
            ylen,
            operations,
            mode: AlignmentMode::Local,
        }
    }
    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    ///
    /// `xstart` is always 0 and `xend` is always `x.len()`. So this algorithm first finds `ystart`
    /// and `yend`, thus determining the coordinates of the termini of the optimal alignment. Then,
    /// a global alignment on the corresponding substrings calculates the operations.
    ///
    /// - Space complexity: $O(n)$
    /// - Time complexity: $(mn + mn')$ where $m$ and $n$ are lengths of the input sequences `x` and
    ///   `y`; $n'$ is the substrings of `y` that correspond to the optimal alignment.
    pub fn semiglobal(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let (score, ystart, yend) = self.find_semiglobal_score_and_termini(x, y);
        let ylen = yend - ystart;
        let xlen = x.len();
        let operations = self.compute_recursive(
            x,
            &y[ystart..yend],
            xlen,
            ylen,
            self.scoring.gap_open,
            self.scoring.gap_open,
        );
        // (doesn't need to insert Yclip)

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        Alignment {
            score,
            xstart: 0,
            xend: xlen,
            ystart,
            yend,
            xlen: x.len(),
            ylen: y.len(),
            operations,
            mode: AlignmentMode::Semiglobal,
        }
    }

    /// Recursively compute alignments of sub-sequences and concatenating them.
    ///
    /// - Space complexity: $O(n)$. Precisely, about $128n + log\_2{m}$, where $m =$ `x.len()` and $n =$ `y.len()`.
    ///   $128n$ is used for four `Vec<i32>` in `find_mid()` and $log\_2{m}$ is due to the activation stack (the
    ///   latter can be negligible).
    ///
    /// # Panics
    ///
    /// Rust has a [`recursion_limit` which defaults to 128](https://doc.rust-lang.org/reference/attributes/limits.html#the-recursion_limit-attribute),
    /// which means the length of the sequence `x` should not exceed $2^128 = 3.4\times10^38$ (well, you'll
    /// never encounter such gigantic biological sequences)
    pub fn compute_recursive(
        &self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> Vec<AlignmentOperation> {
        if self.parallel {
            self.compute_recursive_parallel(x, y, m, n, tb, te)
        } else {
            self.compute_recursive_noparallel(x, y, m, n, tb, te)
        }
    }

    fn compute_recursive_noparallel(
        &self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> Vec<AlignmentOperation> {
        // * m = x.len(); n = y.len()
        if n == 0 {
            return vec![AlignmentOperation::Ins; m];
        }
        if m == 0 {
            return vec![AlignmentOperation::Del; n];
        }
        if n == 1 {
            return self.nw_onecol(x, y[0], m, max(tb, te));
        }
        let (imid, jmid, join_by_deletion) = self.find_mid_noparallel(x, y, m, n, tb, te); // `find_mid()` uses 128m bits
        return if join_by_deletion {
            let a = self.compute_recursive_noparallel(
                &x[..imid],
                &y[..jmid - 1],
                imid,
                jmid - 1,
                tb,
                0,
            );
            let b = self.compute_recursive_noparallel(
                &x[imid..],
                &y[jmid + 1..],
                m - imid,
                n - jmid - 1,
                0,
                te,
            );
            [a, vec![AlignmentOperation::Del; 2], b].concat()
        } else {
            let a = self.compute_recursive_noparallel(
                &x[..imid],
                &y[..jmid],
                imid,
                jmid,
                tb,
                self.scoring.gap_open,
            );
            let b = self.compute_recursive_noparallel(
                &x[imid..],
                &y[jmid..],
                m - imid,
                n - jmid,
                self.scoring.gap_open,
                te,
            );
            [a, b].concat()
        };
    }

    fn compute_recursive_parallel(
        &self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> Vec<AlignmentOperation> {
        // * m = x.len(); n = y.len()
        if n == 0 {
            return vec![AlignmentOperation::Ins; m];
        }
        if m == 0 {
            return vec![AlignmentOperation::Del; n];
        }
        if n == 1 {
            return self.nw_onecol(x, y[0], m, max(tb, te));
        }
        let (imid, jmid, join_by_deletion) = self.find_mid(x, y, m, n, tb, te); // `find_mid()` uses 128n bits
        return if join_by_deletion {
            // x.len() === (&x[..imid]).len() + (&x[imid..]).len()
            // so the sum of the space used by the two subtasks is still 128m (m= x.len())
            let (a, b) = rayon::join(
                || {
                    self.compute_recursive_parallel(
                        &x[..imid],
                        &y[..jmid - 1],
                        imid,
                        jmid - 1,
                        tb,
                        0,
                    )
                },
                || {
                    self.compute_recursive_parallel(
                        &x[imid..],
                        &y[jmid + 1..],
                        m - imid,
                        n - jmid - 1,
                        0,
                        te,
                    )
                },
            );
            [a, vec![AlignmentOperation::Del; 2], b].concat()
        } else {
            let (a, b) = rayon::join(
                || {
                    self.compute_recursive_parallel(
                        &x[..imid],
                        &y[..jmid],
                        imid,
                        jmid,
                        tb,
                        self.scoring.gap_open,
                    )
                },
                || {
                    self.compute_recursive_parallel(
                        &x[imid..],
                        &y[jmid..],
                        m - imid,
                        n - jmid,
                        self.scoring.gap_open,
                        te,
                    )
                },
            );
            [a, b].concat()
        };
    }

    /// Find the "midpoint" (see module-level documentation)
    ///
    /// - Space complexity: $O(m)$, where $m =$ `x.len()`. Specifically, about $128m$ bits, which are occupied
    ///   by the four `Vec<i32>` of length $n$: `cc_left`, `dd_left`, `cc_right` and `dd_right`.
    /// - Time complexity: $O(nm)$
    fn find_mid(
        &self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> (usize, usize, bool) {
        let jmid = n / 2;
        let (a, b) = rayon::join(
            || self.cost_only(x, &y[..jmid], false, tb),
            || self.cost_only(x, &y[jmid..], true, te),
        );
        let (cc_left, dd_left) = a;
        let (cc_right, dd_right) = b;
        let mut max = MIN_SCORE;
        let mut imid = 0;
        let mut join_by_deletion = false;
        for i in 0..=m {
            let c = cc_left[i] + cc_right[m - i];
            if c > max {
                max = c;
                imid = i;
                join_by_deletion = false;
            }
            let d = dd_left[i] + dd_right[m - i] - self.scoring.gap_open; // subtract duplicating open!
            if d > max {
                max = d;
                imid = i;
                join_by_deletion = true;
            }
        }
        (imid, jmid, join_by_deletion)
    }

    /// The non-parallel version of Self::find_mid
    fn find_mid_noparallel(
        &self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> (usize, usize, bool) {
        let jmid = n / 2;
        let (cc_left, dd_left) = self.cost_only(x, &y[..jmid], false, tb);
        let (cc_right, dd_right) = self.cost_only(x, &y[jmid..], true, te);
        let mut max = MIN_SCORE;
        let mut imid = 0;
        let mut join_by_deletion = false;
        for i in 0..=m {
            let c = cc_left[i] + cc_right[m - i];
            if c > max {
                max = c;
                imid = i;
                join_by_deletion = false;
            }
            let d = dd_left[i] + dd_right[m - i] - self.scoring.gap_open; // subtract duplicating open!
            if d > max {
                max = d;
                imid = i;
                join_by_deletion = true;
            }
        }
        (imid, jmid, join_by_deletion)
    }

    pub fn cost_only(
        &self,
        x: TextSlice,
        y: TextSlice,
        rev: bool,
        tx: i32,
    ) -> (Vec<i32>, Vec<i32>) {
        self.cost_only_no_pq_by_col(x, y, rev, tx)
    }

    /// Cost-only (score-only) Gotoh's algorithm in linear space
    ///
    /// - Space Complexity: $O(m)$; specifically, about $64m$ bits, where $m =$ `x.len() + 1`, which
    ///   are used by two `Vec<i32>`.
    /// - Time complexity: $O(mn)$
    pub fn cost_only_no_pq_by_col(
        &self,
        x: TextSlice,
        y: TextSlice,
        rev: bool,
        tx: i32,
    ) -> (Vec<i32>, Vec<i32>) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; m]; // match/mismatch    32 * n bits
        let mut dd: Vec<i32> = vec![i32::MIN; m]; // deletion          32 * n bits
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut s: i32; // C(i-1, j-1)
        let mut t: i32;
        t = self.scoring.gap_open;
        for i in 1..m {
            t += self.scoring.gap_extend;
            cc[i] = t;
        }
        t = tx; // originally self.scoring.gap_open;
        for j in 1..n {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = t;
            cc[0] = c;
            // dd[0] = c;
            e = i32::MIN;
            for i in 1..m {
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend; // update e to I[i,j]
                dd[i] = max(dd[i], cc[i] + self.scoring.gap_open) + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = if rev {
                    max(
                        max(dd[i], e),
                        s + self.scoring.match_fn.score(x[m - i - 1], y[n - j - 1]),
                    )
                } else {
                    max(
                        max(dd[i], e),
                        s + self.scoring.match_fn.score(x[i - 1], y[j - 1]),
                    )
                };
                s = cc[i];
                cc[i] = c;
            }
        }
        dd[0] = cc[0]; // only need to fill cost of deletion in dd once
        (cc, dd)
    }

    pub fn cost_only_no_pq_by_row(
        &self,
        x: TextSlice,
        y: TextSlice,
        rev: bool,
        tx: i32,
    ) -> (Vec<i32>, Vec<i32>) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; n]; // match/mismatch    32 * n bits
        let mut dd: Vec<i32> = vec![i32::MIN; n]; // deletion          32 * n bits
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut s: i32; // C(i-1, j-1)
        let mut t: i32;
        t = self.scoring.gap_open;
        for j in 1..n {
            t += self.scoring.gap_extend;
            cc[j] = t;
        }
        t = tx; // originally self.scoring.gap_open;
        for i in 1..m {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = t;
            cc[0] = c;
            // dd[0] = c;
            e = i32::MIN;
            for j in 1..n {
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = max(dd[j], cc[j] + self.scoring.gap_open) + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = if rev {
                    max(
                        max(dd[j], e),
                        s + self.scoring.match_fn.score(x[m - i - 1], y[n - j - 1]),
                    )
                } else {
                    max(
                        max(dd[j], e),
                        s + self.scoring.match_fn.score(x[i - 1], y[j - 1]),
                    )
                };
                s = cc[j];
                cc[j] = c;
            }
        }
        dd[0] = cc[0]; // otherwise indels at start/end will be free
        (cc, dd)
    }

    pub fn cost_only_pq(
        &self,
        x: TextSlice,
        y: TextSlice,
        rev: bool,
        tx: i32,
    ) -> (Vec<i32>, Vec<i32>) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; m]; // match/mismatch    32 * n bits
        let mut dd: Vec<i32> = vec![i32::MIN; m]; //    deletion          32 * n bits
        let mut e: i32; // I(i - 1, j)
        let mut c: i32; // C(i - 1, j)
        let mut s: i32; // C(i - 1, j - 1)
        let mut t: i32;
        let mut p: u8;
        let mut q: u8;
        t = self.scoring.gap_open;
        for i in 1..m {
            t += self.scoring.gap_extend;
            cc[i] = t;
        }
        t = tx;
        for j in 1..n {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = t;
            cc[0] = c;
            // dd[0] = c;
            e = i32::MIN;
            q = if rev { y[n - j - 1] } else { y[j - 1] };
            for i in 1..m {
                p = if rev { x[m - i - 1] } else { x[i - 1] };
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend; // update e to I[i,j]
                dd[i] = max(dd[i], cc[i] + self.scoring.gap_open) + self.scoring.gap_extend; // cc[i] = C[i, j - 1]
                c = max(max(dd[i], e), s + self.scoring.match_fn.score(p, q));
                s = cc[i];
                cc[i] = c;
            }
        }
        dd[0] = cc[0]; // only need to fill cost of deletion in dd once
        (cc, dd)
    }

    /// Find the maximum score in a local alignment between `x` and `y` and the coordinates
    /// of the alignment path termini
    ///
    /// - Space complexity: $O(m)$; where $m=$ `x.len() + 1`. This function uses 2 `Vec<i32>`s and 2
    ///   `Vec<[usize; 2]>`s, thus taking up about $320m$ bits (x64 architecture) or $192n$ bits
    ///   (x32 architecture)
    /// - Time complexity: $O(nm)$
    ///
    /// # References
    ///
    /// - [Huang, X. and Miller, W. 1991. A time-efficient linear-space local similarity algorithm. Adv. Appl. Math. 12, 3 (Sep. 1991), 337-357](https://doi.org/10.1016/0196-8858(91)90017-D)
    fn find_local_score_and_termini_huang(
        &self,
        x: TextSlice,
        y: TextSlice,
    ) -> (i32, usize, usize, usize, usize) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; m]; // match/mismatch           32 * n bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; m]; // deletion          32 * n bits
        let mut origin_cc: Vec<[usize; 2]> = Vec::with_capacity(m); // usize * 2 * n bits
        let mut origin_dd: Vec<[usize; 2]> = vec![[0, 0]; m]; //       usize * 2 * n bits
        for i in 0..m {
            origin_cc.push([i, 0]);
            // origin_dd.push([i, 0]); values won't be used;
            // dd[0] = [MIN, MIN, ...], dd[i][0] garanteed to be less than c[i] + self.scoring.gap_open
            // origin_ii.push([i, 0]); values won't be used
        }
        let mut e: i32; // I(i - 1, j)
        let mut e_origin: [usize; 2];
        let mut c: i32; // C(i - 1, j)
        let mut c_origin;
        let mut max_c = MIN_SCORE;
        let mut c_start = [0usize, 0usize];
        let mut c_end = [m, n];
        let mut s: i32; // C(i-1, j-1)
        let mut s_origin: [usize; 2]; // origin of C(i-1, j-1)
        let mut xi: u8;
        let mut yj: u8;
        for j in 1..n {
            // s and cc[0] = 0; cc[0] always equals to 0
            s = 0;
            s_origin = [0, j - 1];
            c = 0;
            c_origin = [0, j];
            origin_cc[0] = [0, j];
            e = MIN_SCORE; // I[0, j] garanteed to be less than c + self.scoring.gap_open
                           // origin_dd[0] = [0, j];  value won't be used
            e_origin = [0, j]; // value won't be used; may be optimised
            yj = y[j - 1];
            for i in 1..m {
                xi = x[i - 1];
                e = if e > c + self.scoring.gap_open {
                    //e_origin = e_origin;
                    e
                } else {
                    e_origin = c_origin;
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend; // update e to I[i,j]
                dd[i] = if dd[i] > cc[i] + self.scoring.gap_open {
                    // origin_dd[j] = origin_dd[j];
                    dd[i]
                } else {
                    origin_dd[i] = origin_cc[i];
                    cc[i]
                } + self.scoring.gap_extend; // cc[i] = C[i, j-1]
                c = s + self.scoring.match_fn.score(xi, yj); // substitution score
                c_origin = s_origin;
                s = cc[i];
                s_origin = origin_cc[i];
                if dd[i] > c {
                    c = dd[i];
                    c_origin = origin_dd[i];
                }
                if e > c {
                    c = e;
                    c_origin = e_origin;
                }
                if c < 0 {
                    // the critical step in local alignment
                    c = 0;
                    c_origin = [i, j]
                }
                if c >= max_c {
                    max_c = c;
                    c_start = c_origin;
                    c_end = [i, j];
                }
                origin_cc[i] = c_origin;
                cc[i] = c;
            }
        }
        (max_c, c_start[0], c_start[1], c_end[0], c_end[1])
    }

    /// Find the maximum score in a local alignment between `x` and `y` and the coordinates
    /// of the alignment path termini, based on Shamir's approach
    ///
    /// - Space complexity: $O(m)$; specifically, about $64m$ bits, where $m=$ `x.len() + 1)
    /// - Time complexity: $O(nm)$ (slightly slower than Huang's approach)
    fn find_local_score_and_termini_shamir(
        &self,
        x: TextSlice,
        y: TextSlice,
    ) -> (i32, usize, usize, usize, usize) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; m]; // match/mismatch           32 * n bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; m]; // deletion          32 * n bits
        let mut e: i32; // I(i - 1, j)
        let mut c: i32; // C(i - 1, j)
        let mut max_c = MIN_SCORE;
        let mut xstart = 0usize;
        let mut ystart = 0usize;
        let mut xend = m;
        let mut yend = n;
        let mut s: i32; // C(i-1, j-1)
        let mut p: u8;
        let mut q: u8;
        for j in 1..n {
            s = 0;
            c = 0;
            e = MIN_SCORE;
            q = y[j - 1];
            for i in 1..m {
                p = x[i - 1];
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend;
                dd[i] = max(dd[i], cc[i] + self.scoring.gap_open) + self.scoring.gap_extend;
                c = max(s + self.scoring.match_fn.score(p, q), max(dd[i], max(e, 0)));
                s = cc[i];
                if c >= max_c {
                    max_c = c;
                    xend = i;
                    yend = j;
                }
                cc[i] = c;
            }
        }
        // run in reverse
        let x = &x[..xend];
        let y = &y[..yend];
        cc = vec![0; xend + 1];
        dd = vec![MIN_SCORE; xend + 1];
        max_c = MIN_SCORE;
        for j in 1..=yend {
            s = 0;
            c = 0;
            e = MIN_SCORE;
            q = y[yend - j];
            for i in 1..=xend {
                p = x[xend - i];
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend;
                dd[i] = max(dd[i], cc[i] + self.scoring.gap_open) + self.scoring.gap_extend;
                c = max(s + self.scoring.match_fn.score(p, q), max(dd[i], max(e, 0)));
                s = cc[i];
                if c >= max_c {
                    max_c = c;
                    xstart = i; // subtract from xend later
                    ystart = j;
                }
                cc[i] = c;
            }
        }
        xstart = xend - xstart;
        ystart = yend - ystart;
        (max_c, xstart, ystart, xend, yend)
    }

    /// Find the maximum score in a semiglobal alignment of `x` against `y` (x is global, y is local),
    /// and the coordinates of the alignment path termini.
    ///
    /// In semiglobal mode, `xstart = 0` and `xend = m`, so only `ystart` and `yend` are computed
    /// and returned.
    ///
    /// - Space complexity: $O(m)$, where $m = $ `x.len()`. This function uses 2 `Vec<i32>` and 3
    ///   `Vec<usize>`, thus using `256m` bits (x64 architecture) or `160m` bits (x32 atchitecture).
    /// - Time complexity: $O(mn)$
    fn find_semiglobal_score_and_termini(&self, x: TextSlice, y: TextSlice) -> (i32, usize, usize) {
        let (m, n) = (x.len(), y.len());
        let mut cc: Vec<i32> = Vec::with_capacity(m + 1); //                   32 * m bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; m + 1]; //                      32 * m bits
        let mut y_origin_cc: Vec<usize> = vec![0; m + 1]; //                usize * m bits
        let mut y_origin_dd: Vec<usize> = vec![0; m + 1]; //                usize * m bits

        let mut e: i32; // I(i, j-1)
        let mut e_y_origin: usize;
        let mut c: i32; // C(i, j-1)
        let mut c_y_origin;
        let mut s: i32; // C(i-1, j-1)
        let mut s_y_origin: usize; // y_origin of C(i-1, j-1)
        let mut p: u8;
        let mut q: u8;
        let mut tmp: i32;
        let mut ystart = 0usize;
        let mut yend = n;
        let mut max_last_row = MIN_SCORE;
        cc.push(0i32);
        let mut t = self.scoring.gap_open;
        for _i in 1..=m {
            t += self.scoring.gap_extend;
            cc.push(t);
        }
        for j in 1..=n {
            s = 0;
            c = 0;
            // cc[0] = 0; Never read
            e = MIN_SCORE;

            s_y_origin = j - 1;
            c_y_origin = j;
            e_y_origin = j;
            // y_origin_cc[0] = j;
            q = y[j - 1];
            for i in 1..=m {
                p = x[i - 1];
                tmp = c + self.scoring.gap_open;
                e = if e > tmp {
                    e
                } else {
                    e_y_origin = c_y_origin;
                    tmp
                } + self.scoring.gap_extend;

                tmp = cc[i] + self.scoring.gap_open;
                dd[i] = if dd[i] > tmp {
                    // y_origin_dd[i] = y_origin_dd[i];
                    dd[i]
                } else {
                    y_origin_dd[i] = y_origin_cc[i];
                    tmp
                } + self.scoring.gap_extend;

                c = s + self.scoring.match_fn.score(p, q);
                c_y_origin = s_y_origin;
                s = cc[i];
                s_y_origin = y_origin_cc[i];
                if dd[i] > c {
                    c = dd[i];
                    c_y_origin = y_origin_dd[i];
                }
                if e > c {
                    c = e;
                    c_y_origin = e_y_origin;
                }
                cc[i] = c;
                y_origin_cc[i] = c_y_origin;

                if i == m && c > max_last_row {
                    // c + self.scoring.yclip_suffix
                    max_last_row = c;
                    ystart = c_y_origin;
                    yend = j;
                }
            }
        }
        (max_last_row, ystart, yend)
    }

    // ! Different from the original custom aligner, this one does not allow a gap to align to another gap!
    /// if all 4 clip scores are set to zero, then in all these alignments, terminal gaps are not penalized
    /// ```ignore
    /// ---TTGGCC  AAATTGG--  --TTGG---  AATTGGCCC
    /// AAATTGG--  ---TTGGCC  AATTGGCCC  --TTGG---
    /// ```
    /// But a gap is not allowed to align to another gap:
    /// ```ignore
    /// --AATTG  ATGAT--
    /// ---ATTG  ATGAT---
    /// ```
    fn find_custom_score_and_termini_no_gap_overlap(
        &self,
        x: TextSlice,
        y: TextSlice,
    ) -> (i32, usize, usize, usize, usize) {
        let (m, n) = (x.len() + 1, y.len() + 1); // TODO: to be standardised
        let mut cc: Vec<i32> = Vec::with_capacity(n); //          32 * n bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; n]; //              32 * n bits
        let mut origin_cc: Vec<usize> = vec![0; n]; //            usize * n bits
        let mut clip_x_cc: Vec<bool> = vec![false; n]; //         8 * n bits
        let mut origin_dd: Vec<usize> = vec![0; n]; //            usize * n bits
        let mut clip_x_dd: Vec<bool> = vec![false; n]; //         8 * n bits
        let mut origin_ii: Vec<usize> = vec![0; n]; //            usize * n bits
        let mut clip_x_ii: Vec<bool> = vec![false; n]; //         8 * n bits

        // About clip_x: false: clip y (start = [0, j]); true: clip x (start = [i, 0]), or start = [0,0]
        cc.push(0);
        let mut t = self.scoring.gap_open;
        for _j in 1..n {
            t += self.scoring.gap_extend; // TODO: may be optimised: no need to add after reaching yclip_prefix; same later for xclip_prefix
            cc.push(if t > self.scoring.yclip_prefix {
                t
            } else {
                self.scoring.yclip_prefix
            });
        } // end of 0th row initialisation
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut max_last_column_or_row = MIN_SCORE; // tracks the maximum of the last column
        let mut c_origin;
        let mut c_clip_x: bool;
        let mut xstart = 0usize;
        let mut ystart = 0usize;
        let mut xend = m;
        let mut yend = n;
        let mut s: i32; // C(i-1, j-1)
        let mut s_origin: usize; // origin of C(i-1, j-1)
        let mut s_clip_x: bool;
        t = self.scoring.gap_open;
        for i in 1..m {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = if t > self.scoring.xclip_prefix {
                t
            } else {
                self.scoring.xclip_prefix
            };
            cc[0] = c;
            e = MIN_SCORE;

            s_origin = i - 1; // origin_cc[0] (prev)
            s_clip_x = true; //  clip_x_cc[0] (prev)
            c_origin = i;
            c_clip_x = true;
            origin_cc[0] = i;
            clip_x_cc[0] = true;
            for j in 1..n {
                e = if e > c + self.scoring.gap_open {
                    origin_ii[j] = origin_ii[j - 1];
                    clip_x_ii[j] = clip_x_ii[j - 1];
                    e
                } else {
                    origin_ii[j] = c_origin;
                    clip_x_ii[j] = c_clip_x;
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    // origin_dd[j] = origin_dd[j];
                    // clip_x_dd[j] = clip_x_dd[j];
                    dd[j]
                } else {
                    origin_dd[j] = origin_cc[j];
                    clip_x_dd[j] = clip_x_cc[j];
                    cc[j]
                } + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = s + self.scoring.match_fn.score(x[i - 1], y[j - 1]); // substitution score
                c_origin = s_origin;
                c_clip_x = s_clip_x;
                s = cc[j];
                s_origin = origin_cc[j];
                s_clip_x = clip_x_cc[j];
                if dd[j] > c {
                    c = dd[j];
                    c_origin = origin_dd[j];
                    c_clip_x = clip_x_dd[j];
                }
                if e > c {
                    c = e;
                    c_origin = origin_ii[j];
                    c_clip_x = clip_x_ii[j];
                }
                cc[j] = c;
                origin_cc[j] = c_origin;
                clip_x_cc[j] = c_clip_x;

                if j == n && c + self.scoring.xclip_suffix > max_last_column_or_row {
                    max_last_column_or_row = c + self.scoring.xclip_suffix;
                    xend = i;
                    // yend = n unchanged;
                    if c_clip_x {
                        xstart = c_origin;
                        ystart = 0;
                    } else {
                        xstart = 0;
                        ystart = c_origin;
                    }
                }
            }
        }
        // last (m-th) row
        max_last_column_or_row = if max_last_column_or_row > cc[n] {
            max_last_column_or_row
        } else {
            if clip_x_cc[n] {
                xstart = origin_cc[n];
                ystart = 0;
            } else {
                xstart = 0;
                ystart = origin_cc[n];
            }
            cc[n]
        };
        for j in 0..n {
            if cc[j] + self.scoring.yclip_suffix > max_last_column_or_row {
                max_last_column_or_row = cc[j] + self.scoring.yclip_suffix;
                xend = m;
                yend = j;
                if clip_x_cc[j] {
                    xstart = origin_cc[j];
                    ystart = 0;
                } else {
                    xstart = 0;
                    ystart = origin_cc[j];
                }
            }
        }
        (max_last_column_or_row, xstart, ystart, xend, yend)
    }

    /// Compute the (global) alignment operations between a single letter sequence `y` and a
    /// second sequence `x`. The second sequence can be empty, i.e. `b""`
    fn nw_onecol(&self, x: TextSlice<'_>, y: u8, m: usize, tx: i32) -> Vec<AlignmentOperation> {
        let score_by_indels_only =
            tx + self.scoring.gap_extend * (m as i32 + 1) + self.scoring.gap_open;
        let mut max = score_by_indels_only;
        let score_with_one_substitution_base =
            (m as i32 - 1) * self.scoring.gap_extend + self.scoring.gap_open; // plus substitution score and possibly one more gap_open
        let mut maxi_ = 0usize;
        for i_ in 0..m {
            // index of sequence instead of matrix; y[j] instead of j[j-1] is the jth character
            let score = score_with_one_substitution_base
                + self.scoring.match_fn.score(x[i_], y)
                + if i_ == 0 || i_ == m - 1 {
                    0
                } else {
                    self.scoring.gap_open
                };
            if score > max {
                max = score;
                maxi_ = i_;
            }
        }
        return if max == score_by_indels_only {
            let mut res = Vec::with_capacity(m + 1);
            res.push(AlignmentOperation::Del);
            for _i in 0..m {
                res.push(AlignmentOperation::Ins)
            }
            res
        } else {
            let mut res = Vec::with_capacity(m);
            for _j in 0..maxi_ {
                res.push(AlignmentOperation::Ins)
            }
            if x[maxi_] == y {
                res.push(AlignmentOperation::Match);
            } else {
                res.push(AlignmentOperation::Subst);
            }
            for _i in 0..(m - maxi_ - 1) {
                res.push(AlignmentOperation::Ins)
            }
            res
        };
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
    fn test_blosum62() {
        let x = b"AAAA";
        let y = b"AAAA";
        let score = &blosum62;
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.score, 16);
        assert!(equivalent_operations(
            &alignment.operations,
            &[Match, Match, Match, Match]
        ));
    }

    #[test]
    fn test_issue11() {
        let y = b"TACC"; //GTGGAC";
        let x = b"AAAAACC"; //GTTGACGCAA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert!(equivalent_operations(
            &alignment.operations,
            &[Ins, Ins, Ins, Subst, Match, Match, Match]
        ));
    }

    #[test]
    fn test_left_aligned_del() {
        let x = b"GTGCATCATGTG";
        let y = b"GTGCATCATCATGTG";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);
        println!("\naln:\n{}", alignment.pretty(x, y));

        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        // assert_eq!(
        //     alignment.operations,
        //     [
        //         Match, Match, Match, Del, Del, Del, Match, Match, Match, Match, Match, Match,
        //         Match, Match, Match,
        //     ]
        // );
        assert!(equivalent_operations(
            &alignment.operations,
            &[
                Match, Match, Match, Del, Del, Del, Match, Match, Match, Match, Match, Match,
                Match, Match, Match,
            ]
        ))
    }

    // Test that trailing deletions are correctly handled
    // in global mode
    #[test]
    fn test_global_right_del() {
        let x = b"AACCACGTACGTGGGGGGA";
        let y = b"CCACGTACGT";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);

        println!("\naln:\n{}", alignment.pretty(x, y));

        println!("score:{}", alignment.score);
        assert_eq!(alignment.score, -9);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert!(equivalent_operations(
            &alignment.operations,
            &[
                Ins, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Ins, Ins, Ins, Ins, Ins, Ins, Ins,
            ]
        ));
    }

    #[test]
    fn test_left_aligned_ins() {
        let x = b"GTGCATCATCATGTG";
        let y = b"GTGCATCATGTG";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.global(x, y);
        println!("\naln:\n{}", alignment.pretty(x, y));

        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert!(equivalent_operations(
            &alignment.operations,
            &[
                Match, Match, Match, Ins, Ins, Ins, Match, Match, Match, Match, Match, Match,
                Match, Match, Match,
            ]
        ));
    }

    // semiglobal
    #[test]
    fn test_semiglobal() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = Aligner::new(-5, -1, score);
        let alignment = aligner.semiglobal(x, y);
        println!("{:?}", alignment);
        println!("{}", alignment.pretty(x, y));
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert!(equivalent_operations(
            &alignment.operations,
            &[Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        ));
    }

    // local

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
