// Copyright 2020 Tianyi Shi
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Alignment with affine gap penalty in linear space, by combining Gotoh's (1982) and
//! Hirschberg's (1975) ideas, which was first implemented in C (Myers & Miller 1988).
//!
//! Myers & Miller originally used their technique to implement global alignment only,
//! but alignments of other modes can be achieved by 1) finding the termini of
//! the "partial" alignment 2) global-aligning the resulting substrings.
//!
//! # Time Complexity
//!
//! $O(nm)$ for strings of length $m$ and $n$.
//!
//! # Space Complexity
//!
//! $O(n)$, i.e. linear. For exact number of bits, see `global`, `semiglobal`, `local`, and `custom` methods
//!
//! # References
//!
//! - [Eugene W. Myers and Webb Miller (1988) Optimal alignments in linear space. _Bioinformatics_ **4**: 11-17.](https://doi.org/10.1093/bioinformatics/4.1.11)
//! - [Hirschberg, D. S. (1975) A linear space algorithm for computing maximal common subsequences. _Commun. Assoc. Comput. Mach._ **18**: 341-343.](https://doi.org/10.1145/360825.360861)
//! - [Gotoh, O. (1982) An improved algorithm for matching biological sequences. _J. Molec. Biol._ **162**: 705-708.](https://doi.org/10.1016/0022-2836(82)90398-9)

use crate::alignment::pairwise::{MatchFunc, Scoring};
use crate::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use crate::utils::TextSlice;
use std::cmp::max;

pub struct Aligner<F: MatchFunc + Sync> {
    scoring: Scoring<F>,
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
        }
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
        let score = self.cost_only(x, y, false, self.scoring.gap_open).0[y.len()];
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
    pub fn semiglobal(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = i32::MIN;
        self.scoring.xclip_suffix = i32::MIN;
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
    /// Recursively compute alignments of sub-sequences and concatenating them
    fn compute_recursive(
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
        if m == 1 {
            return self.nw_onerow(x[0], y, n, tb, te);
        }
        let (imid, jmid, join_by_deletion) = self.find_mid(x, y, m, n, tb, te);
        return if join_by_deletion {
            [
                self.compute_recursive(&x[..imid - 1], &y[..jmid], imid - 1, jmid, tb, 0),
                vec![AlignmentOperation::Ins; 2],
                self.compute_recursive(&x[imid + 1..], &y[jmid..], m - imid - 1, n - jmid, 0, te),
            ]
            .concat()
        } else {
            let (a, b) = rayon::join(
                || {
                    (&self).compute_recursive(
                        &x[..imid],
                        &y[..jmid],
                        imid,
                        jmid,
                        tb,
                        self.scoring.gap_open,
                    )
                },
                || {
                    self.compute_recursive(
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

    /// Find the "midpoint"
    fn find_mid(
        &self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> (usize, usize, bool) {
        let imid = m / 2;
        let (a, b) = rayon::join(
            || self.cost_only(&x[..imid], y, false, tb),
            || self.cost_only(&x[imid..], y, true, te),
        );
        let (cc_upper, dd_upper) = a;
        let (cc_lower, dd_lower) = b;
        let mut max = i32::MIN;
        let mut jmid = 0;
        let mut join_by_deletion = false;
        for j in 0..=n {
            let c = cc_upper[j] + cc_lower[n - j];
            if c > max {
                max = c;
                jmid = j;
                join_by_deletion = false;
            }
            let d = dd_upper[j] + dd_lower[n - j] - self.scoring.gap_open; // subtract duplicating open!
            if d > max {
                max = d;
                jmid = j;
                join_by_deletion = true;
            }
        }
        (imid, jmid, join_by_deletion)
    }

    /// Cost-only (score-only) Gotoh's algorithm in linear space
    ///
    /// - Space Complexity: $O(n)$; specifically, about $64n$ bits, where $n =$ `y.len() + 1`
    /// Use six scalars and two vectors of length (N + 1), where N is the length
    /// of the shorter sequence.
    /// -Time complexity: $O(nm)$
    fn cost_only(&self, x: TextSlice, y: TextSlice, rev: bool, tx: i32) -> (Vec<i32>, Vec<i32>) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; n]; // match/mismatch    32 * n bits
        let mut dd: Vec<i32> = vec![0; n]; // deletion          32 * n bits
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut s: i32; // C(i-1, j-1)
        let mut t: i32;
        t = self.scoring.gap_open;
        for j in 1..n {
            t += self.scoring.gap_extend;
            cc[j] = t;
            dd[j] = i32::MIN;
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

    /// Find the maximum score in a local alignment between `x` and `y` and the coordinates
    /// of the alignment path termini
    ///
    /// - Space complexity: $O(n)$; specifically, about $448n$ bits (x64 architecture)
    ///   or $256n$ bits (x32 architecture), where $n=$ `y.len() + 1)
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
        let mut cc: Vec<i32> = vec![0; n]; // match/mismatch           32 * n bits
        let mut dd: Vec<i32> = vec![i32::MIN; n]; // deletion          32 * n bits
        let mut origin_cc: Vec<[usize; 2]> = Vec::with_capacity(n); // usize * 2 * n bits
        let mut origin_dd: Vec<[usize; 2]> = vec![[0, 0]; n]; // usize * 2 * n bits
        let mut origin_ii: Vec<[usize; 2]> = vec![[0, 0]; n]; // usize * 2 * n bits
        for j in 0..n {
            origin_cc.push([0, j]);
            // origin_dd.push([0, j]); values won't be used;
            // dd[0] = [MIN, MIN, ...], dd[0][j] garanteed to be less than c[j] + self.scoring.gap_open
            // origin_ii.push([0, j]); values won't be used
        }
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut c_origin;
        let mut max_c = i32::MIN;
        let mut c_start = [0usize, 0usize];
        let mut c_end = [m, n];
        let mut s: i32; // C(i-1, j-1)
        let mut s_origin: [usize; 2]; // origin of C(i-1, j-1)
        for i in 1..m {
            // s and cc[0] = 0; cc[0] always equals to 0
            s = 0;
            s_origin = [i - 1, 0];
            c = 0;
            c_origin = [i, 0];
            origin_cc[0] = [i, 0];
            e = i32::MIN; // I[i, 0] garanteed to be less than c + self.scoring.gap_open
                          // origin_dd[0] = [i, 0];  value won't be used
                          // origin_ii[0] = [i, 0];  value won't be used
            for j in 1..n {
                e = if e > c + self.scoring.gap_open {
                    origin_ii[j] = origin_ii[j - 1];
                    e
                } else {
                    origin_ii[j] = c_origin;
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    // origin_dd[j] = origin_dd[j];
                    dd[j]
                } else {
                    origin_dd[j] = origin_cc[j];
                    cc[j]
                } + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = s + self.scoring.match_fn.score(x[i - 1], y[j - 1]); // substitution score
                c_origin = s_origin;
                s = cc[j];
                s_origin = origin_cc[j];
                if dd[j] > c {
                    c = dd[j];
                    c_origin = origin_dd[j];
                }
                if e > c {
                    c = e;
                    c_origin = origin_ii[j];
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
                origin_cc[j] = c_origin;
                cc[j] = c;
            }
        }
        (max_c, c_start[0], c_start[1], c_end[0], c_end[1])
    }

    /// Find the maximum score in a local alignment between `x` and `y` and the coordinates
    /// of the alignment path termini, based on Shamir's approach
    ///
    /// - Space complexity: $O(n)$; specifically, about $64n$ bits, where $n=$ `y.len() + 1)
    /// - Time complexity: $O(nm)$ (slightly slower than the original `find_local_score_and_termini`)
    fn find_local_score_and_termini_shamir(
        &self,
        x: TextSlice,
        y: TextSlice,
    ) -> (i32, usize, usize, usize, usize) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; n]; // match/mismatch           32 * n bits
        let mut dd: Vec<i32> = vec![i32::MIN; n]; // deletion          32 * n bits
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut max_c = i32::MIN;
        let mut xstart = 0usize;
        let mut ystart = 0usize;
        let mut xend = m;
        let mut yend = n;
        let mut s: i32; // C(i-1, j-1)
        for i in 1..m {
            s = 0;
            c = 0;
            e = i32::MIN;
            for j in 1..n {
                e = if e > c + self.scoring.gap_open {
                    e
                } else {
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend;
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    dd[j]
                } else {
                    cc[j]
                } + self.scoring.gap_extend;
                c = s + self.scoring.match_fn.score(x[i - 1], y[j - 1]);
                s = cc[j];
                if dd[j] > c {
                    c = dd[j];
                }
                if e > c {
                    c = e;
                }
                if c < 0 {
                    c = 0;
                }
                if c >= max_c {
                    max_c = c;
                    xend = i;
                    yend = j;
                }
                cc[j] = c;
            }
        }
        // run in reverse
        let x = &x[..xend];
        let y = &y[..yend];
        cc = vec![0; n];
        dd = vec![i32::MIN; n];
        max_c = i32::MIN;
        for i in 1..=xend {
            s = 0;
            c = 0;
            e = i32::MIN;
            for j in 1..=yend {
                e = if e > c + self.scoring.gap_open {
                    e
                } else {
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend;
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    dd[j]
                } else {
                    cc[j]
                } + self.scoring.gap_extend;
                c = s + self.scoring.match_fn.score(x[xend - i], y[yend - j]);
                s = cc[j];
                if dd[j] > c {
                    c = dd[j];
                }
                if e > c {
                    c = e;
                }
                if c < 0 {
                    c = 0;
                }
                if c >= max_c {
                    max_c = c;
                    xstart = i; // subtract from xend later
                    ystart = j;
                }
                cc[j] = c;
            }
        }
        xstart = xend - xstart;
        ystart = yend - ystart;
        (max_c, xstart, ystart, xend, yend)
    }

    /// Find the maximum score in a semiglobal alignment of `x` against `y` (x is global, y is local),
    /// and the coordinates of the alignment path termini
    ///
    /// In semiglobal mode, `xstart === 0` and `xend === m`, so only `ystart` and `yend` are computed
    /// and returned
    fn find_semiglobal_score_and_termini(&self, x: TextSlice, y: TextSlice) -> (i32, usize, usize) {
        let (m, n) = (x.len(), y.len());
        let mut cc: Vec<i32> = vec![0; n + 1]; //                           32 * n bits
        let mut dd: Vec<i32> = vec![i32::MIN; n + 1]; //                    32 * n bits
        let mut y_origin_cc: Vec<usize> = (0..=n).into_iter().collect(); // usize * n bits
        let mut y_origin_dd: Vec<usize> = vec![0; n + 1]; //                usize * n bits
        let mut y_origin_ii: Vec<usize> = vec![0; n + 1]; //                usize * n bits
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut c_y_origin;
        let mut s: i32; // C(i-1, j-1)
        let mut s_y_origin: usize; // y_origin of C(i-1, j-1)
        let mut t = self.scoring.gap_open;
        for i in 1..=m {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = t;
            cc[0] = c;
            e = i32::MIN;

            s_y_origin = 0;
            c_y_origin = 0;
            // y_origin_cc[0] = 0;WON'T CHANGE
            for j in 1..=n {
                e = if e > c + self.scoring.gap_open {
                    y_origin_ii[j] = y_origin_ii[j - 1];
                    e
                } else {
                    y_origin_ii[j] = c_y_origin;
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    // y_origin_dd[j] = y_origin_dd[j];
                    dd[j]
                } else {
                    y_origin_dd[j] = y_origin_cc[j];
                    cc[j]
                } + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = s + self.scoring.match_fn.score(x[i - 1], y[j - 1]); // substitution score
                c_y_origin = s_y_origin;
                s = cc[j];
                s_y_origin = y_origin_cc[j];
                if dd[j] > c {
                    c = dd[j];
                    c_y_origin = y_origin_dd[j];
                }
                if e > c {
                    c = e;
                    c_y_origin = y_origin_ii[j];
                }
                cc[j] = c;
                y_origin_cc[j] = c_y_origin;
            }
        }
        // last (m-th) row
        let mut max_last_row = cc[n];
        let mut ystart = y_origin_cc[n];
        let mut yend = n;
        // in semiglobal mode, xstart === 0; xend === m
        for j in 0..=n {
            let clip_score = cc[j] + self.scoring.yclip_suffix;
            if clip_score > max_last_row {
                max_last_row = clip_score;
                ystart = y_origin_cc[j];
                yend = j;
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
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = Vec::with_capacity(n); //          32 * n bits
        let mut dd: Vec<i32> = vec![i32::MIN; n]; //              32 * n bits
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
        let mut max_last_column_or_row = i32::MIN; // tracks the maximum of the last column
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
            e = i32::MIN;

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

    /// Compute the (global) alignment operations between a single letter sequence `x` and a
    /// second sequence `y`. The second sequence can be empty, i.e. `b""`
    fn nw_onerow(
        &self,
        x: u8,
        y: TextSlice,
        n: usize,
        tb: i32,
        te: i32,
    ) -> Vec<AlignmentOperation> {
        let score_by_indels_only =
            max(tb, te) + self.scoring.gap_extend * (n as i32 + 1) + self.scoring.gap_open;
        let mut max = score_by_indels_only;
        let score_with_one_substitution_base =
            (n as i32 - 1) * self.scoring.gap_extend + self.scoring.gap_open; // plus substitution score and possibly one more gap_open
        let mut maxj_ = 0usize;
        for j_ in 0..n {
            // index of sequence instead of matrix; y[j] instead of j[j-1] is the jth character
            let score = score_with_one_substitution_base
                + self.scoring.match_fn.score(x, y[j_])
                + if j_ == 0 || j_ == n - 1 {
                    0
                } else {
                    self.scoring.gap_open
                };
            if score > max {
                max = score;
                maxj_ = j_;
            }
        }
        return if max == score_by_indels_only {
            let mut res = Vec::with_capacity(n + 1);
            res.push(AlignmentOperation::Ins);
            for _j in 0..n {
                res.push(AlignmentOperation::Del)
            }
            res
        } else {
            let mut res = Vec::with_capacity(n);
            for _j in 0..maxj_ {
                res.push(AlignmentOperation::Del)
            }
            if x == y[maxj_] {
                res.push(AlignmentOperation::Match);
            } else {
                res.push(AlignmentOperation::Subst);
            }
            for _j in 0..(n - maxj_ - 1) {
                res.push(AlignmentOperation::Del)
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
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Ins, Ins, Ins, Match, Match, Match]
        );
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

        assert_eq!(alignment.operations, correct);
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
        assert_eq!(alignment.operations, [Match, Match, Match, Match]);
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
        assert_eq!(
            alignment.operations,
            [Ins, Ins, Ins, Subst, Match, Match, Match]
        );
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
        assert_eq!(
            alignment.operations,
            [
                Ins, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Ins, Ins, Ins, Ins, Ins, Ins, Ins,
            ]
        );
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
        assert_eq!(
            alignment.operations,
            [
                Match, Match, Match, Ins, Ins, Ins, Match, Match, Match, Match, Match, Match,
                Match, Match, Match,
            ]
        );
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
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );
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
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );
    }
}
