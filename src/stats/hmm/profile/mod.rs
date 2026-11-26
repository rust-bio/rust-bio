// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Profile hidden Markov models for multiple sequence alignments.
//!
//! A [`ProfileHmm`] is built from a block of equal-length, pre-aligned rows in which the byte `-`
//! marks a gap. Every column that holds a residue in at least half of the (weighted) rows becomes a
//! position of the model. Each position owns a match state, an insert state and a delete state
//! (the "Plan 7" architecture of Eddy, 1998), so that a query can be aligned with runs of inserted
//! residues and runs of skipped positions. Emission and transition probabilities are estimated
//! from the counts of the block, smoothed with a configurable pseudocount.
//!
//! The model offers
//!
//! * the most probable alignment of a query ([`ProfileHmm::best_alignment`], Viterbi),
//! * the total probability of a query over all alignments ([`ProfileHmm::forward`],
//!   [`ProfileHmm::backward`]) and per-position match posteriors
//!   ([`ProfileHmm::posterior_match`]),
//! * the alignment with maximum expected accuracy ([`ProfileHmm::optimal_accuracy`]),
//! * null-corrected bit scores and a ranked search over many queries
//!   ([`ProfileHmm::bit_score`], [`ProfileHmm::search`]),
//! * local alignment and scanning for repeated hits ([`ProfileHmm::local_alignment`],
//!   [`ProfileHmm::scan`]),
//! * Baum-Welch and Viterbi re-estimation ([`ProfileHmm::reestimate`],
//!   [`ProfileHmm::reestimate_viterbi`]),
//! * sampling of alignments and sequences ([`ProfileHmm::sample_alignment`],
//!   [`ProfileHmm::generate`]),
//! * a plain-text form that round-trips ([`ProfileHmm::to_text`], [`ProfileHmm::from_text`]).
//!
//! All probabilities are handled as natural logarithms. Scores reported as bit scores are
//! log-odds against a background model that emits residues independently with the frequencies
//! observed in the input block.
//!
//! # Example
//!
//! ```
//! use approx::assert_relative_eq;
//! use bio::alphabets::Alphabet;
//! use bio::stats::hmm::profile::{ProfileHmm, ProfileOp};
//!
//! let mut builder = ProfileHmm::builder(&Alphabet::new(b"ACGT"));
//! builder.add_row(b"ACGT");
//! builder.add_row(b"ACGT");
//! builder.add_row(b"AC-T");
//! let profile = builder.build().unwrap();
//! assert_eq!(profile.num_columns(), 4);
//! assert_eq!(profile.consensus(), b"ACGT".to_vec());
//!
//! // The most probable alignment of a query that equals the consensus matches every position.
//! let alignment = profile.best_alignment(b"ACGT").unwrap();
//! assert_eq!(
//!     alignment.ops,
//!     vec![
//!         ProfileOp::Match(b'A'),
//!         ProfileOp::Match(b'C'),
//!         ProfileOp::Match(b'G'),
//!         ProfileOp::Match(b'T'),
//!     ]
//! );
//!
//! // Forward and backward agree on the total probability of a query.
//! let forward = profile.forward(b"ACGT").unwrap();
//! let backward = profile.backward(b"ACGT").unwrap();
//! assert_relative_eq!(forward, backward, epsilon = 1e-9);
//!
//! // Queries resembling the profile score above queries that do not.
//! assert!(profile.bit_score(b"ACGT").unwrap() > profile.bit_score(b"TTTT").unwrap());
//!
//! // The text form can be written out and read back without changing any score.
//! let restored = ProfileHmm::from_text(&profile.to_text()).unwrap();
//! assert_eq!(restored.forward(b"ACGT").unwrap(), forward);
//! ```
//!
//! # References
//!
//! * Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
//!   Cambridge University Press. <http://doi.org/10.1017/CBO9780511790492>
//! * Eddy, S. R. (1998). Profile hidden Markov models. Bioinformatics, 14(9), 755-763.
//!   <https://doi.org/10.1093/bioinformatics/14.9.755>

use thiserror::Error;

use crate::alphabets::Alphabet;

mod alignment;
mod analyze;
mod builder;
mod decode;
mod forward_backward;
mod local;
mod sample;
mod score;
mod text;
mod train;
mod viterbi;

pub use self::alignment::{ProfileAlignment, ProfileOp};
pub use self::analyze::PAD;
pub use self::builder::{ProfileBuilder, GAP};
pub use self::local::LocalHit;
pub use self::score::SearchHit;

pub(crate) const NEG_INF: f64 = f64::NEG_INFINITY;

/// Add two probabilities given as natural logarithms.
#[inline]
pub(crate) fn ln_add(a: f64, b: f64) -> f64 {
    if a == NEG_INF {
        return b;
    }
    if b == NEG_INF {
        return a;
    }
    if a >= b {
        a + (b - a).exp().ln_1p()
    } else {
        b + (a - b).exp().ln_1p()
    }
}

/// Errors that can occur while building, parsing or applying a [`ProfileHmm`].
#[derive(Error, Clone, Debug, PartialEq)]
pub enum ProfileError {
    /// The builder holds no rows.
    #[error("alignment has no rows")]
    NoRows,
    /// A row has a different length than the first row.
    #[error("alignment rows differ in length: expected {expected}, found {found}")]
    RaggedRows { expected: usize, found: usize },
    /// The rows have no columns.
    #[error("alignment columns are empty")]
    EmptyAlignment,
    /// No column holds a residue in at least half of the weighted rows.
    #[error("no column qualifies as a match column")]
    NoMatchColumns,
    /// The pseudocount is zero, negative or not finite.
    #[error("pseudocount must be positive and finite, got {0}")]
    NonPositivePseudocount(f64),
    /// A residue is not part of the model's alphabet.
    #[error("symbol {0} is not in the alphabet")]
    UnknownSymbol(u8),
    /// The text form passed to [`ProfileHmm::from_text`] is malformed.
    #[error("cannot parse profile: {0}")]
    Parse(String),
}

/// A profile hidden Markov model with match, insert and delete states per position.
///
/// Positions are numbered `1..=num_columns()`. Position `0` stands for the begin state; the end is
/// reached by the transitions leaving position `num_columns()`. Use [`ProfileHmm::builder`] to
/// create a model from a block of aligned rows.
#[derive(Clone, Debug)]
pub struct ProfileHmm {
    pub(crate) num_cols: usize,
    pub(crate) symbols: Vec<u8>,
    pub(crate) rank_of: Vec<i32>,
    pub(crate) emission: Vec<Vec<f64>>,
    pub(crate) null: Vec<f64>,
    pub(crate) t_mm: Vec<f64>,
    pub(crate) t_mi: Vec<f64>,
    pub(crate) t_md: Vec<f64>,
    pub(crate) t_im: Vec<f64>,
    pub(crate) t_ii: Vec<f64>,
    pub(crate) t_id: Vec<f64>,
    pub(crate) t_dm: Vec<f64>,
    pub(crate) t_di: Vec<f64>,
    pub(crate) t_dd: Vec<f64>,
}

impl ProfileHmm {
    /// Start building a model over the residues of `alphabet`. The gap byte `-` is never part of
    /// the model's alphabet, even if `alphabet` contains it.
    pub fn builder(alphabet: &Alphabet) -> ProfileBuilder {
        ProfileBuilder::new(alphabet)
    }

    /// The number of positions of the model, that is the number of aligned columns of the block
    /// it was built from.
    pub fn num_columns(&self) -> usize {
        self.num_cols
    }

    pub(crate) fn rank(&self, symbol: u8) -> Result<usize, ProfileError> {
        let idx = symbol as usize;
        if idx < self.rank_of.len() && self.rank_of[idx] >= 0 {
            Ok(self.rank_of[idx] as usize)
        } else {
            Err(ProfileError::UnknownSymbol(symbol))
        }
    }

    pub(crate) fn ranks(&self, query: &[u8]) -> Result<Vec<usize>, ProfileError> {
        query.iter().map(|&c| self.rank(c)).collect()
    }

    pub(crate) fn null_log_prob_ranks(&self, ranks: &[usize]) -> f64 {
        ranks.iter().map(|&r| self.null[r]).sum()
    }

    /// The natural-log probability of `query` under the background model, the sum of the
    /// background log-probabilities of its residues.
    pub fn null_log_prob(&self, query: &[u8]) -> Result<f64, ProfileError> {
        let ranks = self.ranks(query)?;
        Ok(self.null_log_prob_ranks(&ranks))
    }

    /// The most probable alignment of `query` against the model (Viterbi). Ties are broken in
    /// favour of match, then delete, then insert states.
    pub fn best_alignment(&self, query: &[u8]) -> Result<ProfileAlignment, ProfileError> {
        viterbi::viterbi(self, query)
    }

    /// The natural logarithm of the total probability of `query`, summed over all alignments
    /// (forward algorithm).
    pub fn forward(&self, query: &[u8]) -> Result<f64, ProfileError> {
        forward_backward::forward(self, query)
    }

    /// The natural logarithm of the total probability of `query`, summed over all alignments
    /// (backward algorithm). Agrees with [`ProfileHmm::forward`].
    pub fn backward(&self, query: &[u8]) -> Result<f64, ProfileError> {
        forward_backward::backward(self, query)
    }

    /// For every position of the model, the posterior probability in `[0, 1]` that some residue
    /// of `query` is aligned to the match state of that position.
    pub fn posterior_match(&self, query: &[u8]) -> Result<Vec<f64>, ProfileError> {
        forward_backward::posterior_match(self, query)
    }

    /// The alignment of `query` that maximises the expected number of correctly matched
    /// positions. Its `log_prob` never exceeds that of [`ProfileHmm::best_alignment`].
    pub fn optimal_accuracy(&self, query: &[u8]) -> Result<ProfileAlignment, ProfileError> {
        decode::optimal_accuracy(self, query)
    }

    /// The expected number of correctly matched positions of the alignment returned by
    /// [`ProfileHmm::optimal_accuracy`]; a value between `0` and `num_columns()`.
    pub fn expected_accuracy(&self, query: &[u8]) -> Result<f64, ProfileError> {
        decode::expected_accuracy(self, query)
    }

    /// The posterior match probabilities of those positions that the most probable alignment of
    /// `query` matches, in position order.
    pub fn posterior_path(&self, query: &[u8]) -> Result<Vec<f64>, ProfileError> {
        decode::posterior_path(self, query)
    }

    /// The natural-log probability of the given alignment. The operations must consume exactly
    /// `num_columns()` positions, otherwise the result is negative infinity.
    pub fn score_path(&self, ops: &[ProfileOp]) -> Result<f64, ProfileError> {
        decode::path_log_prob(self, ops)
    }

    /// The log-odds score in bits of `query`: `(forward(query) - null_log_prob(query)) / ln 2`.
    pub fn bit_score(&self, query: &[u8]) -> Result<f64, ProfileError> {
        score::bit_score(self, query)
    }

    /// Score every query and return those whose bit score is at least `threshold`, by descending
    /// bit score and by ascending index among equal scores.
    pub fn search(
        &self,
        queries: &[&[u8]],
        threshold: f64,
    ) -> Result<Vec<SearchHit>, ProfileError> {
        score::search(self, queries, threshold)
    }

    /// The best local alignment of the model to a region of `query`. The ends are free, and the
    /// score is never negative; if no region scores above zero, an all-zero hit is returned.
    pub fn local_alignment(&self, query: &[u8]) -> Result<LocalHit, ProfileError> {
        local::local_alignment(self, query)
    }

    /// Find non-overlapping local hits scoring at least `min_bits` in `query`, ordered by
    /// `query_start`. The best hit is taken first and the flanking regions are searched again.
    pub fn scan(&self, query: &[u8], min_bits: f64) -> Result<Vec<LocalHit>, ProfileError> {
        local::scan(self, query, min_bits)
    }

    /// Sample an alignment of `query` from the posterior distribution over alignments. The
    /// result only depends on `query` and `seed`.
    pub fn sample_alignment(
        &self,
        query: &[u8],
        seed: u64,
    ) -> Result<ProfileAlignment, ProfileError> {
        sample::sample_alignment(self, query, seed)
    }

    /// Emit a sequence from the model. The result only depends on `seed`.
    pub fn generate(&self, seed: u64) -> Vec<u8> {
        sample::generate(self, seed)
    }

    /// The most probable residue of every position.
    pub fn consensus(&self) -> Vec<u8> {
        analyze::consensus(self)
    }

    /// The information content in bits of every position, the relative entropy of its emission
    /// distribution against the background model.
    pub fn position_information(&self) -> Vec<f64> {
        analyze::position_information(self)
    }

    /// The sum of [`ProfileHmm::position_information`] over all positions.
    pub fn relative_entropy(&self) -> f64 {
        analyze::relative_entropy(self)
    }

    /// Render the consensus followed by the best alignment of every query as rows of
    /// [`ProfileHmm::align_all`], one per line.
    pub fn to_alignment_text(&self, queries: &[&[u8]]) -> Result<String, ProfileError> {
        analyze::to_alignment_text(self, queries)
    }

    /// Align all queries and lay them out as a rectangular block. A skipped position is written
    /// as [`GAP`]; residues inserted after a position are followed by [`PAD`] bytes so that every
    /// row has the width of the widest insertion at that position.
    pub fn align_all(&self, queries: &[&[u8]]) -> Result<Vec<Vec<u8>>, ProfileError> {
        analyze::align_all(self, queries)
    }

    /// The total natural-log likelihood of all `sequences`.
    pub fn log_likelihood(&self, sequences: &[&[u8]]) -> Result<f64, ProfileError> {
        train::total_log_likelihood(self, sequences)
    }

    /// Re-estimate the parameters with `iterations` rounds of Baum-Welch. Returns the retrained
    /// model and the log-likelihood of `sequences` before the first and after every round, so the
    /// trace has `iterations + 1` entries and never decreases.
    pub fn reestimate(
        &self,
        sequences: &[&[u8]],
        iterations: usize,
    ) -> Result<(ProfileHmm, Vec<f64>), ProfileError> {
        train::reestimate(self, sequences, iterations)
    }

    /// Like [`ProfileHmm::reestimate`], but counts are taken from the most probable alignment of
    /// every sequence (Viterbi training). The trace holds the summed log-probability of those
    /// alignments.
    pub fn reestimate_viterbi(
        &self,
        sequences: &[&[u8]],
        iterations: usize,
    ) -> Result<(ProfileHmm, Vec<f64>), ProfileError> {
        train::reestimate_viterbi(self, sequences, iterations)
    }

    /// Write the model in a plain-text form that [`ProfileHmm::from_text`] reads back without
    /// losing precision.
    pub fn to_text(&self) -> String {
        text::to_text(self)
    }

    /// Read a model written by [`ProfileHmm::to_text`].
    pub fn from_text(text: &str) -> Result<ProfileHmm, ProfileError> {
        text::from_text(text)
    }
}

#[cfg(test)]
pub(crate) fn example_profile() -> ProfileHmm {
    let mut builder = ProfileHmm::builder(&Alphabet::new(b"ACGT"));
    for row in [b"ACGTAC", b"ACGTAC", b"AC--AC", b"AC--AC", b"ACGTAC"] {
        builder.add_row(row);
    }
    builder.build().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_null_log_prob_is_additive() {
        let profile = example_profile();
        let ac = profile.null_log_prob(b"AC").unwrap();
        let a = profile.null_log_prob(b"A").unwrap();
        let c = profile.null_log_prob(b"C").unwrap();
        assert!((ac - (a + c)).abs() < 1e-9);
        assert!(profile.null_log_prob(b"").unwrap().abs() < 1e-12);
    }

    #[test]
    fn test_unknown_symbol_is_rejected() {
        let profile = example_profile();
        assert_eq!(
            profile.forward(b"ACX").unwrap_err(),
            ProfileError::UnknownSymbol(b'X')
        );
        assert_eq!(
            profile.best_alignment(b"N").unwrap_err(),
            ProfileError::UnknownSymbol(b'N')
        );
    }

    #[test]
    fn test_ln_add() {
        assert_eq!(ln_add(NEG_INF, NEG_INF), NEG_INF);
        assert_eq!(ln_add(NEG_INF, -1.5), -1.5);
        assert!((ln_add(0.5f64.ln(), 0.25f64.ln()) - 0.75f64.ln()).abs() < 1e-12);
    }
}
