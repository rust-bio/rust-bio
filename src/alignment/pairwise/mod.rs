// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Copyright 2020 Tianyi Shi (PR #354: space-efficiency and speed improvement)
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
//! use bio::alignment::AlignmentOperation::*;
//!
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//! let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
//! let alignment = aligner.semiglobal(x, y);
//! // x is global (target sequence) and y is local (reference sequence)
//! assert_eq!(alignment.ystart, 4);
//! assert_eq!(alignment.xstart, 0);
//! assert_eq!(
//!     alignment.operations,
//!     [Match, Match, Match, Match, Match, Subst, Match, Match, Match]
//! );
//!
//! // If you don't know sizes of future sequences, you could
//! // use Aligner::new().
//! // Global alignment:
//! let aligner = Aligner::new(-5, -1, &score);
//! let x = b"ACCGTGGAT";
//! let y = b"AAAAACCGTTGAT";
//! let alignment = aligner.global(x, y);
//! assert_eq!(alignment.ystart, 0);
//! assert_eq!(alignment.xstart, 0);
//! assert_eq!(aligner.local(x, y).score, 7);
//!
//! // In addition to the standard modes (Global, Semiglobal and Local), a custom alignment
//! // mode is supported which supports a user-specified clipping penalty. Clipping is a
//! // special boundary condition where you are allowed to clip off the beginning/end of
//! // the sequence for a fixed penalty. As a starting example, we can use the custom mode
//! // for achieving the three standard modes as follows.
//!
//! // scoring for semiglobal mode
//! let scoring = Scoring::new(-5, -1, &score) // Gap open, gap extend and match score function
//!     .xclip(MIN_SCORE) // Clipping penalty for x set to 'negative infinity', hence global in x
//!     .yclip(0); // Clipping penalty for y set to 0, hence local in y
//! let aligner = Aligner::with_scoring(scoring);
//! let alignment = aligner.custom(x, y); // The custom aligner invocation
//! assert_eq!(alignment.ystart, 4);
//! assert_eq!(alignment.xstart, 0);
//! // Note that in the custom mode, the clips are explicitly mentioned in the operations
//! assert_eq!(
//!     alignment.operations,
//!     [
//!         Yclip(4),
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Subst,
//!         Match,
//!         Match,
//!         Match
//!     ]
//! );
//!
//! // scoring for global mode
//! // scoring can also be created using from_scores if the match and mismatch scores are constants
//! let scoring = Scoring::from_scores(-5, -1, 1, -1) // Gap open, extend, match, mismatch score
//!     .xclip(MIN_SCORE) // Clipping penalty for x set to 'negative infinity', hence global in x
//!     .yclip(MIN_SCORE); // Clipping penalty for y set to 'negative infinity', hence global in y
//! let aligner = Aligner::with_scoring(scoring);
//! let alignment = aligner.custom(x, y); // The custom aligner invocation
//! assert_eq!(alignment.ystart, 0);
//! assert_eq!(alignment.xstart, 0);
//! // Note that in the custom mode, the clips are explicitly mentioned in the operations
//! assert_eq!(
//!     alignment.operations,
//!     [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match]
//! );
//!
//! // Similarly if the clip penalties are both set to 0, we have local alignment mode. The scoring
//! // struct also lets users set different penalties for prefix/suffix clipping, thereby letting
//! // users have the flexibility to create a wide variety of boundary conditions. The xclip() and
//! // yclip() methods sets the prefix and suffix penalties to be equal. The scoring struct can be
//! // explicitly constructed for full flexibility.
//!
//! // The following example considers a modification of the semiglobal mode where you are allowed
//! // to skip a prefix of the target sequence x, for a penalty of -10, but you have to consume
//! // the rest of the string in the alignment
//!
//! let scoring = Scoring {
//!     gap_open: -5,
//!     gap_extend: -1,
//!     match_fn: |a: u8, b: u8| if a == b { 1i32 } else { -3i32 },
//!     match_scores: Some((1, -3)),
//!     xclip_prefix: -10,
//!     xclip_suffix: MIN_SCORE,
//!     yclip_prefix: 0,
//!     yclip_suffix: 0,
//! };
//! let x = b"GGGGGGACGTACGTACGT";
//! let y = b"AAAAACGTACGTACGTAAAA";
//! let aligner = Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring);
//! let alignment = aligner.custom(x, y);
//! println!("{}", alignment.pretty(x, y));
//! assert_eq!(alignment.score, 2);
//! assert_eq!(
//!     alignment.operations,
//!     [
//!         Yclip(4),
//!         Xclip(6),
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Match,
//!         Yclip(4)
//!     ]
//! );
//! ```

#![allow(non_snake_case)]

use std::cmp::max;
use std::i32;

use crate::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use crate::utils::TextSlice;

pub mod banded;

/// Value to use as a 'negative infinity' score. Should be close to `i32::MIN`,
/// but avoid underflow when used with reasonable scoring parameters or even
/// adding two negative infinities. Use ~ `0.4 * i32::MIN`
pub const MIN_SCORE: i32 = -858_993_459;

/// Trait required to instantiate a Scoring instance
pub trait MatchFunc {
    fn score(&self, a: u8, b: u8) -> i32;
}

/// A concrete data structure which implements trait MatchFunc with constant
/// match and mismatch scores
#[derive(Debug, Clone)]
pub struct MatchParams {
    pub match_score: i32,
    pub mismatch_score: i32,
}

impl MatchParams {
    /// Create new MatchParams instance with given match and mismatch scores
    ///
    /// # Arguments
    ///
    /// * `match_score` - the score for a match (should not be negative)
    /// * `mismatch_score` - the score for a mismatch (should not be positive)
    pub fn new(match_score: i32, mismatch_score: i32) -> Self {
        assert!(match_score >= 0, "match_score can't be negative");
        assert!(mismatch_score <= 0, "mismatch_score can't be positive");
        MatchParams {
            match_score,
            mismatch_score,
        }
    }
}

impl MatchFunc for MatchParams {
    #[inline]
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The trait Matchfunc is also implemented for Fn(u8, u8) -> i32 so that Scoring
/// can be instantiated using closures and custom user defined functions
impl<F> MatchFunc for F
where
    F: Fn(u8, u8) -> i32,
{
    fn score(&self, a: u8, b: u8) -> i32 {
        (self)(a, b)
    }
}

/// Details of scoring are encapsulated in this structure.
///
/// An [affine gap score model](https://en.wikipedia.org/wiki/Gap_penalty#Affine)
/// is used so that the gap score for a length `k` is:
/// `GapScore(k) = gap_open + gap_extend * k`
#[derive(Debug, Clone)]
pub struct Scoring<F: MatchFunc> {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub match_fn: F,
    pub match_scores: Option<(i32, i32)>,
    pub xclip_prefix: i32,
    pub xclip_suffix: i32,
    pub yclip_prefix: i32,
    pub yclip_suffix: i32,
}

impl Scoring<MatchParams> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// match and mismatch scores. The clip penalties are set to `MIN_SCORE` by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_score` - the score for a match
    /// * `mismatch_score` - the score for a mismatch
    pub fn from_scores(
        gap_open: i32,
        gap_extend: i32,
        match_score: i32,
        mismatch_score: i32,
    ) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn: MatchParams::new(match_score, mismatch_score),
            match_scores: Some((match_score, mismatch_score)),
            xclip_prefix: MIN_SCORE,
            xclip_suffix: MIN_SCORE,
            yclip_prefix: MIN_SCORE,
            yclip_suffix: MIN_SCORE,
        }
    }
}

impl<F: MatchFunc> Scoring<F> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// and the score function. The clip penalties are set to [`MIN_SCORE`](constant.MIN_SCORE.html) by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn,
            match_scores: None,
            xclip_prefix: MIN_SCORE,
            xclip_suffix: MIN_SCORE,
            yclip_prefix: MIN_SCORE,
            yclip_suffix: MIN_SCORE,
        }
    }

    /// Sets the prefix and suffix clipping penalties for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Clipping penalty for x (both prefix and suffix, should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip(-5);
    /// assert!(scoring.xclip_prefix == -5);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == -5);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn xclip(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_prefix = penalty;
        self.xclip_suffix = penalty;
        self
    }

    /// Sets the prefix clipping penalty for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Prefix clipping penalty for x (should not be positive)
    ///
    /// # Example
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip_prefix(-5);
    /// assert!(scoring.xclip_prefix == -5);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn xclip_prefix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_prefix = penalty;
        self
    }

    /// Sets the suffix clipping penalty for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Suffix clipping penalty for x (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip_suffix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == -5);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn xclip_suffix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_suffix = penalty;
        self
    }

    /// Sets the prefix and suffix clipping penalties for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Clipping penalty for y (both prefix and suffix, should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == -5);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == -5);
    /// ```
    pub fn yclip(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self.yclip_suffix = penalty;
        self
    }

    /// Sets the prefix clipping penalty for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Prefix clipping penalty for y (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip_prefix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == -5);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn yclip_prefix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self
    }

    /// Sets the suffix clipping penalty for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Suffix clipping penalty for y (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip_suffix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == -5);
    /// ```
    pub fn yclip_suffix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_suffix = penalty;
        self
    }
}

/// A generalized Smith-Waterman aligner.
///
/// `M(i,j)` is the best score such that `x[i]` and `y[j]` ends in a match (or substitution)
/// ```ignore
///              .... A   G  x_i
///              .... C   G  y_j
/// ```
/// `I(i,j)` is the best score such that `x[i]` is aligned with a gap
/// ```ignore
///              .... A   G  x_i
///              .... G  y_j  -
/// ```
/// This is interpreted as an insertion into `x` w.r.t reference `y`
///
/// `D(i,j)` is the best score such that `y[j]` is aligned with a gap
/// ```ignore
///              .... A  x_i  -
///              .... G   G  y_j
/// ```
/// This is interpreted as a deletion from `x` w.r.t reference `y`
///
/// `S(i,j)` is the best score for prefixes `x[0..i]`, `y[0..j]`
///
/// To save space, only one column of these matrices is stored at any
/// point — the column `j` is obtained by overwriting column `j-1` (Myers & Miller 1988).
/// Moreover, `M(i,j)` is not explicitly stored.
///
/// `Lx` is the optimal x suffix clipping lengths from each position of the
/// sequence y
///
/// `Ly` is the optimal y suffix clipping lengths from each position of the
/// sequence x
///
/// `Sn` is the last column of the matrix. This is needed to keep track of
/// suffix clipping scores
///
/// `traceback` - see [`bio::alignment::pairwise::TracebackCell`](struct.TracebackCell.html)
///
/// `scoring` - see [`bio::alignment::pairwise::Scoring`](struct.Scoring.html)
///
/// # References
/// - [Eugene W. Myers and Webb Miller (1988) Optimal alignments in linear space. _Bioinformatics_ **4**: 11-17.](https://doi.org/10.1093/bioinformatics/4.1.11)
pub struct Aligner<F: MatchFunc> {
    scoring: Scoring<F>,
}

const DEFAULT_ALIGNER_CAPACITY: usize = 200;

impl<F: MatchFunc> Aligner<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        Aligner::with_capacity(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            gap_open,
            gap_extend,
            match_fn,
        )
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
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn with_capacity(m: usize, n: usize, gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Aligner {
            scoring: Scoring::new(gap_open, gap_extend, match_fn),
        }
    }

    /// Create new aligner instance with given the scoring struct
    ///
    /// # Arguments
    ///
    /// * `scoring` - the scoring struct (see bio::alignment::pairwise::Scoring)
    pub fn with_scoring(scoring: Scoring<F>) -> Self {
        Aligner::with_capacity_and_scoring(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            scoring,
        )
    }

    /// Create new aligner instance with scoring and size hint. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `scoring` - the scoring struct
    pub fn with_capacity_and_scoring(m: usize, n: usize, scoring: Scoring<F>) -> Self {
        assert!(scoring.gap_open <= 0, "gap_open can't be positive");
        assert!(scoring.gap_extend <= 0, "gap_extend can't be positive");
        assert!(
            scoring.xclip_prefix <= 0,
            "Clipping penalty (x prefix) can't be positive"
        );
        assert!(
            scoring.xclip_suffix <= 0,
            "Clipping penalty (x suffix) can't be positive"
        );
        assert!(
            scoring.yclip_prefix <= 0,
            "Clipping penalty (y prefix) can't be positive"
        );
        assert!(
            scoring.yclip_suffix <= 0,
            "Clipping penalty (y suffix) can't be positive"
        );

        Aligner { scoring }
    }

    /// The core function to compute the alignment
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    pub fn custom(self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        use traceback_old::*;
        let (m, n) = (x.len(), y.len());
        let mut I = vec![MIN_SCORE; m + 1];
        let mut D = vec![MIN_SCORE; m + 1];
        let mut S = vec![MIN_SCORE; m + 1];
        let mut Lx = vec![0usize; n + 1];
        let mut Ly = vec![0usize; m + 1];
        let mut Sn = vec![MIN_SCORE; m + 1];
        let mut traceback = Traceback::with_capacity(m + 1, n + 1);
        traceback.init(m + 1, n + 1);

        {
            // k = 0 (first row)
            let mut tb = TracebackCell::new();
            tb.set_all(TB_START);
            traceback.set(0, 0, tb);
            Sn[0] = self.scoring.yclip_suffix;
            Ly[0] = n;
        }
        {
            // j = 0 (manipulation of S, I, D)
            S[0] = 0;

            for i in 1..=m {
                let mut tb = TracebackCell::new();
                tb.set_all(TB_START);
                if i == 1 {
                    I[i] = self.scoring.gap_open + self.scoring.gap_extend;
                    tb.set_i_bits(TB_START);
                } else {
                    // Insert all i characters
                    let i_score = self.scoring.gap_open + self.scoring.gap_extend * (i as i32);
                    let c_score =
                        self.scoring.xclip_prefix + self.scoring.gap_open + self.scoring.gap_extend; // Clip then insert
                    I[i] = if i_score > c_score {
                        tb.set_i_bits(TB_INS);
                        i_score
                    } else {
                        tb.set_i_bits(TB_XCLIP_PREFIX);
                        c_score
                    };
                }

                if i == m {
                    tb.set_s_bits(TB_XCLIP_SUFFIX);
                }

                if I[i] > S[i] {
                    S[i] = I[i];
                    tb.set_s_bits(TB_INS);
                }

                if self.scoring.xclip_prefix > S[i] {
                    S[i] = self.scoring.xclip_prefix;
                    tb.set_s_bits(TB_XCLIP_PREFIX);
                }

                // Track the score if we do a suffix clip (x) after this character
                if i != m && S[i] + self.scoring.xclip_suffix > S[m] {
                    S[m] = S[i] + self.scoring.xclip_suffix;
                    Lx[0] = m - i;
                }

                traceback.set(i, 0, tb);

                // Track the score if we do suffix clip (y) from here
                if S[i] + self.scoring.yclip_suffix > Sn[i] {
                    Sn[i] = S[i] + self.scoring.yclip_suffix;
                    Ly[i] = n;
                }
            }
        }

        let mut d0 = self.scoring.gap_open; // gap_open + j * gap_extend
        let mut s: i32; // becomes prev S[i-1] for the next round of iteration; update s before updating current S[i]
        let mut c: i32; // becomes curr S[i-1]; update c after updating current S[i]
        for j in 1..=n {
            d0 += self.scoring.gap_extend;
            {
                // Handle i = 0 case
                let mut tb = TracebackCell::new();
                I[0] = MIN_SCORE;

                if j == 1 {
                    D[0] = d0;
                    tb.set_d_bits(TB_START);
                } else {
                    // Delete all j characters
                    let c_score =
                        self.scoring.yclip_prefix + self.scoring.gap_open + self.scoring.gap_extend;
                    D[0] = if d0 > c_score {
                        tb.set_d_bits(TB_DEL);
                        d0
                    } else {
                        tb.set_d_bits(TB_YCLIP_PREFIX);
                        c_score
                    };
                }

                s = S[0];
                // ! S update start
                S[0] = if D[0] > self.scoring.yclip_prefix {
                    tb.set_s_bits(TB_DEL);
                    D[0]
                } else {
                    tb.set_s_bits(TB_YCLIP_PREFIX);
                    self.scoring.yclip_prefix
                };

                if j == n && Sn[0] > S[0] {
                    // Check if the suffix clip score is better
                    S[0] = Sn[0];
                    tb.set_s_bits(TB_YCLIP_SUFFIX);
                // Track the score if we do suffix clip (y) from here
                } else if S[0] + self.scoring.yclip_suffix > Sn[0] {
                    Sn[0] = S[0] + self.scoring.yclip_suffix;
                    Ly[0] = n - j;
                }
                // ! S update end
                c = S[0];

                traceback.set(0, j, tb);
            }

            let q = y[j - 1];
            let xclip_score = self.scoring.xclip_prefix + max(self.scoring.yclip_prefix, d0);
            let prev_s_m = S[m];
            S[m] = MIN_SCORE;
            for i in 1..=m {
                let p = x[i - 1];
                let mut tb = TracebackCell::new();

                let m_score = s + self.scoring.match_fn.score(p, q);
                let i_score = I[i - 1] + self.scoring.gap_extend;
                let s_score = c + self.scoring.gap_open + self.scoring.gap_extend;
                let best_i_score = if i_score > s_score {
                    tb.set_i_bits(TB_INS);
                    i_score
                } else {
                    tb.set_i_bits(traceback.get(i - 1, j).get_s_bits());
                    s_score
                };

                let d_score = D[i] + self.scoring.gap_extend; // D prev i
                let s_score = self.scoring.gap_open
                    + self.scoring.gap_extend
                    + if i != m { S[i] } else { prev_s_m }; // S prev i
                let best_d_score = if d_score > s_score {
                    tb.set_d_bits(TB_DEL);
                    d_score
                } else {
                    tb.set_d_bits(traceback.get(i, j - 1).get_s_bits());
                    s_score
                };

                tb.set_s_bits(TB_XCLIP_SUFFIX);
                let mut best_s_score = if i != m { MIN_SCORE } else { S[m] };

                if m_score > best_s_score {
                    best_s_score = m_score;
                    tb.set_s_bits(if p == q { TB_MATCH } else { TB_SUBST });
                }

                if best_i_score > best_s_score {
                    best_s_score = best_i_score;
                    tb.set_s_bits(TB_INS);
                }

                if best_d_score > best_s_score {
                    best_s_score = best_d_score;
                    tb.set_s_bits(TB_DEL);
                }

                if xclip_score > best_s_score {
                    best_s_score = xclip_score;
                    tb.set_s_bits(TB_XCLIP_PREFIX);
                }

                let yclip_score = self.scoring.yclip_prefix
                    + self.scoring.gap_open
                    + self.scoring.gap_extend * (i as i32);
                if yclip_score > best_s_score {
                    best_s_score = yclip_score;
                    tb.set_s_bits(TB_YCLIP_PREFIX);
                }

                s = S[i];
                S[i] = best_s_score; // ! S[i] updated
                c = best_s_score;
                I[i] = best_i_score;
                D[i] = best_d_score;

                // Track the score if we do suffix clip (x) from here
                if i != m && c + self.scoring.xclip_suffix > S[m] {
                    S[m] = c + self.scoring.xclip_suffix;
                    Lx[j] = m - i;
                }

                // Track the score if we do suffix clip (y) from here
                if c + self.scoring.yclip_suffix > Sn[i] {
                    Sn[i] = c + self.scoring.yclip_suffix;
                    Ly[i] = n - j;
                }

                traceback.set(i, j, tb);
            }
        }

        // Handle suffix clipping in the j=n case
        for i in 0..=m {
            if Sn[i] > S[i] {
                S[i] = Sn[i];
                traceback.get_mut(i, n).set_s_bits(TB_YCLIP_SUFFIX);
            }
            if S[i] + self.scoring.xclip_suffix > S[m] {
                S[m] = S[i] + self.scoring.xclip_suffix;
                Lx[n] = m - i;
                traceback.get_mut(m, n).set_s_bits(TB_XCLIP_SUFFIX);
            }
        }

        // Since there could be a change in the last column of S,
        // recompute the last column of I as this could also change
        for i in 1..=m {
            let s_score = S[i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            if s_score > I[i] {
                I[i] = s_score;
                let s_bit = traceback.get(i - 1, n).get_s_bits();
                traceback.get_mut(i, n).set_i_bits(s_bit);
            }
            if s_score > S[i] {
                S[i] = s_score;
                traceback.get_mut(i, n).set_s_bits(TB_INS);
                if S[i] + self.scoring.xclip_suffix > S[m] {
                    S[m] = S[i] + self.scoring.xclip_suffix;
                    Lx[n] = m - i;
                    traceback.get_mut(m, n).set_s_bits(TB_XCLIP_SUFFIX);
                }
            }
        }

        let mut i = m;
        let mut j = n;
        let mut operations = Vec::with_capacity(x.len());
        let mut xstart: usize = 0usize;
        let mut ystart: usize = 0usize;
        let mut xend = m;
        let mut yend = n;

        let mut last_layer = traceback.get(i, j).get_s_bits();

        loop {
            let next_layer: u16;
            match last_layer {
                TB_START => break,
                TB_INS => {
                    operations.push(AlignmentOperation::Ins);
                    next_layer = traceback.get(i, j).get_i_bits();
                    i -= 1;
                }
                TB_DEL => {
                    operations.push(AlignmentOperation::Del);
                    next_layer = traceback.get(i, j).get_d_bits();
                    j -= 1;
                }
                TB_MATCH => {
                    operations.push(AlignmentOperation::Match);
                    next_layer = traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                TB_SUBST => {
                    operations.push(AlignmentOperation::Subst);
                    next_layer = traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                TB_XCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Xclip(i));
                    xstart = i;
                    i = 0;
                    next_layer = traceback.get(0, j).get_s_bits();
                }
                TB_XCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Xclip(Lx[j]));
                    i -= Lx[j];
                    xend = i;
                    next_layer = traceback.get(i, j).get_s_bits();
                }
                TB_YCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Yclip(j));
                    ystart = j;
                    j = 0;
                    next_layer = traceback.get(i, 0).get_s_bits();
                }
                TB_YCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Yclip(Ly[i]));
                    j -= Ly[i];
                    yend = j;
                    next_layer = traceback.get(i, j).get_s_bits();
                }
                _ => panic!("Dint expect this!"),
            }
            last_layer = next_layer;
        }

        operations.reverse();
        Alignment {
            score: S[m],
            ystart,
            xstart,
            yend,
            xend,
            ylen: n,
            xlen: m,
            operations,
            mode: AlignmentMode::Custom,
        }
    }

    /// Calculate global alignment of x against y.
    pub fn global(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        use traceback_new::*;
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

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        use traceback_new::*;
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut S = vec![0; m]; //                                            ! 32 * m bits
        let mut D = vec![MIN_SCORE; m]; //                                    ! 32 * m bits
        let mut T: Vec<TracebackCell> = vec![TracebackCell::new(); n * m]; // ! 8 * n * m bits
        let mut s: i32;
        let mut c: i32;
        let mut e: i32;
        let mut idx: usize = 0;
        let mut p: &u8;
        let mut q: &u8;
        let mut S_i: &mut i32;
        let mut D_i: &mut i32;
        let mut score_1: i32;
        let mut score_2: i32;
        let mut max_last_row: i32 = MIN_SCORE;
        let mut yend: usize = n;

        unsafe {
            let mut t = self.scoring.gap_open;
            for i in 1..m {
                t += self.scoring.gap_extend;
                *S.get_unchecked_mut(i) = t;
                let mut tb = TracebackCell::new();
                tb.set_s_bits(TB_UP);
                idx += 1;
                *T.get_unchecked_mut(idx) = tb; // T[0 * m + i]
            }

            t = self.scoring.gap_open;
            for j in 1..n {
                s = 0;
                c = 0;
                *S.get_unchecked_mut(0) = 0;
                e = MIN_SCORE;
                idx += 1;
                *T.get_unchecked_mut(idx) = TracebackCell::new(); // T[j * m + 0]

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
                    tb.set_s_bits(TB_DIAG);

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

                    if i == m - 1 && c > max_last_row {
                        max_last_row = c;
                        yend = j;
                    }
                }
            }

            let (mut operations, i, ystart) = Self::traceback(&T, x, y, m - 1, yend, m);
            operations.resize(operations.len() + i, AlignmentOperation::Ins); // reaching at (i, 0)

            operations.reverse();
            Alignment {
                score: max_last_row,
                xstart: 0,
                ystart,
                xend: m - 1,
                yend,
                xlen: m - 1,
                ylen: n - 1,
                operations,
                mode: AlignmentMode::Semiglobal,
            }
        }
    }

    /// Calculate local alignment of x against y.
    pub fn local(&self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        use traceback_new::*;
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
                // *S.get_unchecked_mut(0) = 0; unchanged
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
                xlen: m - 1,
                ylen: n - 1,
                operations,
                mode: AlignmentMode::Local,
            }
        }
    }

    pub fn traceback(
        T: &Vec<traceback_new::TracebackCell>,
        x: TextSlice,
        y: TextSlice,
        mut i: usize,
        mut j: usize,
        m: usize,
    ) -> (Vec<AlignmentOperation>, usize, usize) {
        use traceback_new::*;
        let mut operations = Vec::with_capacity(m);
        unsafe {
            let mut next_layer = T.get_unchecked(j * m + i).get_s_bits(); // start from the last tb cell
            loop {
                match next_layer {
                    TB_ORIGIN => break,
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

mod traceback_new {

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
    #[derive(Copy, Clone)]
    pub struct TracebackCell(u8);

    pub const TB_ORIGIN: u8 = 0b00;
    pub const TB_UP: u8 = 0b01;
    pub const TB_LEFT: u8 = 0b10;
    pub const TB_DIAG: u8 = 0b11;

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
}

pub mod traceback_old {

    /// Packed representation of one cell of a Smith-Waterman traceback matrix.
    /// Stores the I, D and S traceback matrix values in two bytes.
    /// Possible traceback moves include : start, insert, delete, match, substitute,
    /// prefix clip and suffix clip for x & y. So we need 4 bits each for matrices I, D, S
    /// to keep track of these 9 moves.
    #[derive(Copy, Clone)]
    pub struct TracebackCell {
        v: u16,
    }

    impl Default for TracebackCell {
        fn default() -> Self {
            TracebackCell { v: 0 }
        }
    }

    // Traceback bit positions (LSB)
    const I_POS: u8 = 0; // Meaning bits 0,1,2,3 corresponds to I and so on
    const D_POS: u8 = 4;
    const S_POS: u8 = 8;

    // Traceback moves
    pub const TB_START: u16 = 0b0000;
    pub const TB_INS: u16 = 0b0001;
    pub const TB_DEL: u16 = 0b0010;
    pub const TB_SUBST: u16 = 0b0011;
    pub const TB_MATCH: u16 = 0b0100;

    pub const TB_XCLIP_PREFIX: u16 = 0b0101; // prefix clip of x
    pub const TB_XCLIP_SUFFIX: u16 = 0b0110; // suffix clip of x
    pub const TB_YCLIP_PREFIX: u16 = 0b0111; // prefix clip of y
    pub const TB_YCLIP_SUFFIX: u16 = 0b1000; // suffix clip of y

    const TB_MAX: u16 = 0b1000; // Useful in checking that the
                                // TB value we got is a valid one

    impl TracebackCell {
        /// Initialize a blank traceback cell
        #[inline(always)]
        pub fn new() -> TracebackCell {
            Default::default()
        }

        /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
        #[inline(always)]
        fn set_bits(&mut self, pos: u8, value: u16) {
            let bits: u16 = (0b1111) << pos;
            assert!(
                value <= TB_MAX,
                "Expected a value <= TB_MAX while setting traceback bits"
            );
            self.v = (self.v & !bits) // First clear the bits
        | (value << pos) // And set the bits
        }

        #[inline(always)]
        pub fn set_i_bits(&mut self, value: u16) {
            // Traceback corresponding to matrix I
            self.set_bits(I_POS, value);
        }

        #[inline(always)]
        pub fn set_d_bits(&mut self, value: u16) {
            // Traceback corresponding to matrix D
            self.set_bits(D_POS, value);
        }

        #[inline(always)]
        pub fn set_s_bits(&mut self, value: u16) {
            // Traceback corresponding to matrix S
            self.set_bits(S_POS, value);
        }

        // Gets 4 bits [pos, pos+4) of v
        #[inline(always)]
        fn get_bits(self, pos: u8) -> u16 {
            (self.v >> pos) & (0b1111)
        }

        #[inline(always)]
        pub fn get_i_bits(self) -> u16 {
            self.get_bits(I_POS)
        }

        #[inline(always)]
        pub fn get_d_bits(self) -> u16 {
            self.get_bits(D_POS)
        }

        #[inline(always)]
        pub fn get_s_bits(self) -> u16 {
            self.get_bits(S_POS)
        }

        /// Set all matrices to the same value.
        pub fn set_all(&mut self, value: u16) {
            self.set_i_bits(value);
            self.set_d_bits(value);
            self.set_s_bits(value);
        }
    }

    /// Internal traceback.
    pub struct Traceback {
        rows: usize,
        cols: usize,
        matrix: Vec<TracebackCell>,
    }

    impl Traceback {
        pub fn with_capacity(m: usize, n: usize) -> Self {
            let rows = m + 1;
            let cols = n + 1;
            Traceback {
                rows,
                cols,
                matrix: Vec::with_capacity(rows * cols),
            }
        }

        pub fn init(&mut self, m: usize, n: usize) {
            self.matrix.clear();
            let mut start = TracebackCell::new();
            start.set_all(TB_START);
            // set every cell to start
            self.resize(m, n, start);
        }

        #[inline(always)]
        pub fn set(&mut self, i: usize, j: usize, v: TracebackCell) {
            debug_assert!(i < self.rows);
            debug_assert!(j < self.cols);
            self.matrix[i * self.cols + j] = v;
        }

        #[inline(always)]
        pub fn get(&self, i: usize, j: usize) -> &TracebackCell {
            debug_assert!(i < self.rows);
            debug_assert!(j < self.cols);
            &self.matrix[i * self.cols + j]
        }

        pub fn get_mut(&mut self, i: usize, j: usize) -> &mut TracebackCell {
            debug_assert!(i < self.rows);
            debug_assert!(j < self.cols);
            &mut self.matrix[i * self.cols + j]
        }

        pub fn resize(&mut self, m: usize, n: usize, v: TracebackCell) {
            self.rows = m + 1;
            self.cols = n + 1;
            self.matrix.resize(self.rows * self.cols, v);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::AlignmentOperation::*;
    use crate::scores::blosum62;
    use std::iter::repeat;

    #[test]
    fn traceback_cell_old() {
        use super::traceback_old::*;
        let mut tb = TracebackCell::new();

        tb.set_all(TB_SUBST);
        assert_eq!(tb.get_i_bits(), TB_SUBST);
        assert_eq!(tb.get_d_bits(), TB_SUBST);
        assert_eq!(tb.get_s_bits(), TB_SUBST);

        tb.set_d_bits(TB_INS);
        assert_eq!(tb.get_d_bits(), TB_INS);

        tb.set_i_bits(TB_XCLIP_PREFIX);
        assert_eq!(tb.get_d_bits(), TB_INS);
        assert_eq!(tb.get_i_bits(), TB_XCLIP_PREFIX);

        tb.set_d_bits(TB_DEL);
        assert_eq!(tb.get_d_bits(), TB_DEL);
        assert_eq!(tb.get_i_bits(), TB_XCLIP_PREFIX);

        tb.set_s_bits(TB_YCLIP_SUFFIX);
        assert_eq!(tb.get_d_bits(), TB_DEL);
        assert_eq!(tb.get_i_bits(), TB_XCLIP_PREFIX);
        assert_eq!(tb.get_s_bits(), TB_YCLIP_SUFFIX);
    }

    #[test]
    fn test_semiglobal() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );
    }

    // Test case for underflow of the SW score.
    #[test]
    fn test_semiglobal_gap_open_lt_mismatch() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -5i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -1, -1, score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Match, Del, Match, Ins, Match, Match, Match,]
        );
    }

    #[test]
    fn test_global_affine_ins() {
        let x = b"ACGAGAACA";
        let y = b"ACGACA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
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
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
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
        let y = b"CGTATCATAGATAGATGTAGATGATCCACAGT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.xstart, 1);
        assert_eq!(alignment.ystart, 0);
    }

    #[test]
    fn test_local() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );
    }

    #[test]
    fn test_global() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
        let alignment = aligner.global(x, y);

        println!("\naln:\n{}", alignment.pretty(x, y));
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );
    }

    #[test]
    fn test_blosum62() {
        let x = b"AAAA";
        let y = b"AAAA";
        let score = &blosum62;
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, score);
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
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Ins, Ins, Ins, Subst, Match, Match, Match]
        );
    }

    #[test]
    fn test_issue12_1() {
        let x = b"CCGGCA";
        let y = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 1);
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Subst, Subst, Subst]
        );
    }

    #[test]
    fn test_issue12_2() {
        let y = b"CCGGCA";
        let x = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 0);

        assert_eq!(
            alignment.operations,
            [Subst, Match, Ins, Ins, Ins, Ins, Ins, Ins, Subst, Match, Match,]
        );
    }

    #[test]
    fn test_issue12_3() {
        let y = b"CCGTCCGGCAA";
        let x = b"AAAAACCGTTGACGCAA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [
                Ins, Ins, Ins, Ins, Ins, Ins, Match, Subst, Subst, Match, Subst, Subst, Subst,
                Match, Match, Match, Match,
            ]
        );

        let aligner = Aligner::with_capacity(y.len(), x.len(), -5, -1, &score);
        let alignment = aligner.semiglobal(y, x);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Match, Subst, Subst, Match, Subst, Subst, Subst, Match, Match, Match, Match,]
        );
    }

    #[test]
    fn test_left_aligned_del() {
        let x = b"GTGCATCATGTG";
        let y = b"GTGCATCATCATGTG";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let alignment = aligner.global(x, y);
        println!("\naln:\n{}", alignment.pretty(x, y));

        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [
                Match, Match, Match, Del, Del, Del, Match, Match, Match, Match, Match, Match,
                Match, Match, Match,
            ]
        );
    }

    // Test that trailing deletions are correctly handled
    // in global mode
    #[test]
    fn test_global_right_del() {
        let x = b"AACCACGTACGTGGGGGGA";
        let y = b"CCACGTACGT";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
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
        let aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
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

    #[test]
    fn test_aligner_new() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, &score);

        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );

        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );

        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(
            alignment.operations,
            [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );
    }

    #[test]
    fn test_semiglobal_simple() {
        let x = b"GAAAACCGTTGAT";
        let y = b"ACCGTGGATGGG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let aligner = Aligner::new(-5, -1, &score);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(
            alignment.operations,
            [Ins, Ins, Ins, Ins, Match, Match, Match, Match, Match, Subst, Match, Match, Match,]
        );
    }

    #[test]
    fn test_insert_only_semiglobal() {
        let x = b"TTTT";
        let y = b"AAAA";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let aligner = Aligner::new(-5, -1, &score);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(alignment.operations, [Ins, Ins, Ins, Ins]);
    }

    #[test]
    fn test_insert_in_between_semiglobal() {
        let x = b"GGGGG";
        let y = b"GGTAGGG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let aligner = Aligner::new(-5, -1, &score);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(
            alignment.operations,
            [Match, Match, Del, Del, Match, Match, Match]
        );
    }

    #[test]
    fn test_xclip_prefix_custom() {
        let x = b"GGGGGGATG";
        let y = b"ATG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(-5);

        let aligner = Aligner::with_scoring(scoring);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Xclip(6), Match, Match, Match]);
    }

    #[test]
    fn test_yclip_prefix_custom() {
        let y = b"GGGGGGATG";
        let x = b"ATG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).yclip(-5);

        let aligner = Aligner::with_scoring(scoring);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Yclip(6), Match, Match, Match]);
    }

    #[test]
    fn test_xclip_suffix_custom() {
        let x = b"GAAAA";
        let y = b"CG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(-5).yclip(0);

        let aligner = Aligner::with_scoring(scoring);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Yclip(1), Match, Xclip(4)]);
    }

    #[test]
    fn test_yclip_suffix_custom() {
        let y = b"GAAAA";
        let x = b"CG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let scoring = Scoring::new(-5, -1, &score).yclip(-5).xclip(0);

        let aligner = Aligner::with_scoring(scoring);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Xclip(1), Match, Yclip(4)]);
    }

    #[test]
    fn test_longer_string_all_operations() {
        let x = b"TTTTTGGGGGGATGGCCCCCCTTTTTTTTTTGGGAAAAAAAAAGGGGGG";
        let y = b"GGGGGGATTTCCCCCCCCCTTTTTTTTTTAAAAAAAAA";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(-5).yclip(0);

        let aligner = Aligner::with_scoring(scoring);
        let alignment = aligner.custom(x, y);

        println!("{}", alignment.pretty(x, y));
        assert_eq!(alignment.score, 7);
    }

    #[test]
    fn test_scoring_from_scores() {
        let y = b"GGGGGGATG";
        let x = b"ATG";

        let scoring = Scoring::from_scores(-5, -1, 1, -1).yclip(-5);

        let aligner = Aligner::with_scoring(scoring);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Yclip(6), Match, Match, Match]);
    }

    #[test]
    fn test_only_clips() {
        let x = b"GGAAAAAAAAAAAAA";
        let y = b"TTTTAATTTGTGTAAAAAATAATA";
        let base_score = Scoring::from_scores(-4, -4, 4, -7);
        let scoring = Scoring {
            xclip_prefix: 0,
            xclip_suffix: 0,
            yclip_suffix: 0,
            ..base_score
        };
        let mut al = Aligner::with_scoring(scoring);
        let alignment = al.custom(x, y);
        assert_eq!(alignment.score, 0);
    }

    #[test]
    fn test_zero_score_clips() {
        let x = b"AA";
        let y = b"CC";
        let base_score = Scoring::from_scores(-1, -1, 1, -1);
        {
            let scoring = Scoring {
                xclip_prefix: 0,
                yclip_prefix: 0,
                ..base_score.clone()
            };
            let mut al = Aligner::with_scoring(scoring);
            let alignment = al.custom(x, y);
            assert_eq!(alignment.score, 0);
        }

        {
            let scoring = Scoring {
                xclip_prefix: 0,
                yclip_suffix: 0,
                ..base_score.clone()
            };
            let mut al = Aligner::with_scoring(scoring);
            let alignment = al.custom(x, y);
            assert_eq!(alignment.score, 0);
        }

        {
            let scoring = Scoring {
                xclip_suffix: 0,
                yclip_prefix: 0,
                ..base_score.clone()
            };
            let mut al = Aligner::with_scoring(scoring);
            let alignment = al.custom(x, y);
            assert_eq!(alignment.score, 0);
        }

        {
            let scoring = Scoring {
                xclip_suffix: 0,
                yclip_suffix: 0,
                ..base_score.clone()
            };
            let mut al = Aligner::with_scoring(scoring);
            let alignment = al.custom(x, y);
            assert_eq!(alignment.score, 0);
        }
    }
}
