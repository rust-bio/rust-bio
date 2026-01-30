//! This module contains the implementation of a classic `PairHMM` as described in
//! Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
//! Current Topics in Genome Analysis 2008. http://doi.org/10.1017/CBO9780511790492.
//! It also contains a modified variant `HomopolyPairHMM` with additional homopolymer states suited
//! for dealing with homopolymer runs in sequencing as often encountered in Oxford Nanopore
//! sequencing data.
//!
//! Traits defined in this module apply to both `PairHMM` and `HomopolyPairHMM`.
//!
//! # Examples
//! ```
//! use approx::assert_relative_eq;
//! use bio::stats::pairhmm::{
//!     EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission,
//! };
//! use bio::stats::{LogProb, Prob};
//! use num_traits::Zero;
//!
//! // Two sequences for which we'd like to know if they are likely related.
//! let x = b"AAAA";
//! let y = b"AAAT";
//!
//! // For this example, we disallow gaps, so all probabilities are zero here.
//! struct GapParams;
//! impl GapParameters for GapParams {
//!     fn prob_gap_x(&self) -> LogProb {
//!         LogProb::zero()
//!     }
//!     fn prob_gap_y(&self) -> LogProb {
//!         LogProb::zero()
//!     }
//!     fn prob_gap_x_extend(&self) -> LogProb {
//!         LogProb::zero()
//!     }
//!     fn prob_gap_y_extend(&self) -> LogProb {
//!         LogProb::zero()
//!     }
//! }
//! let gap_params = GapParams;
//!
//! // The PairHMM instance stores the gap params, since these are constant.
//! let mut pairhmm = PairHMM::new(&gap_params);
//!
//! // However, emission parameters depend on the actual sequences
//! struct EmissionParams {
//!     x: &'static [u8],
//!     y: &'static [u8],
//! }
//!
//! const PROB_SUBSTITUTION: f64 = 0.1;
//! const PROB_NO_SUBSTITUION: f64 = 1. - PROB_SUBSTITUTION;
//! impl EmissionParameters for EmissionParams {
//!     fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
//!         if self.x[i] == self.y[j] {
//!             // if two bases match, emit a Match!
//!             XYEmission::Match(LogProb::from(Prob(PROB_NO_SUBSTITUION)))
//!         } else {
//!             // otherwise emit a Mismatch!
//!             // Note that the probability here is `mismatch / 3`, since probabilities should sum
//!             // to 1 and there are 3 possible mismatch configurations
//!             XYEmission::Mismatch(LogProb::from(Prob(PROB_SUBSTITUTION / 3.)))
//!         }
//!     }
//!
//!     // In this example, emitting x[i] is as likely as not observing a mismatch.
//!     // In more complex cases, this might e.g. depend on base qualities reported by the sequencer
//!     fn prob_emit_x(&self, i: usize) -> LogProb {
//!         LogProb::from(Prob(PROB_NO_SUBSTITUION))
//!     }
//!     fn prob_emit_y(&self, j: usize) -> LogProb {
//!         LogProb::from(Prob(PROB_NO_SUBSTITUION))
//!     }
//!
//!     fn len_x(&self) -> usize {
//!         self.x.len()
//!     }
//!     fn len_y(&self) -> usize {
//!         self.y.len()
//!     }
//! }
//!
//! // Since we want to do global alignment here, disallow free start and end gaps in x.
//! struct GlobalAlignmentMode;
//! impl StartEndGapParameters for GlobalAlignmentMode {
//!     fn free_start_gap_x(&self) -> bool {
//!         false
//!     }
//!     fn free_end_gap_x(&self) -> bool {
//!         false
//!     }
//! }
//!
//! // Finally calculate the probability of relatedness between x and y!
//! let prob_related = pairhmm.prob_related(&EmissionParams { x, y }, &GlobalAlignmentMode, None);
//!
//! // … and compare it to a rough estimation
//! let prob_expected = LogProb::from(Prob(PROB_NO_SUBSTITUION.powi(3) * PROB_SUBSTITUTION / 3.));
//! assert_relative_eq!(*prob_related, *prob_expected, epsilon = 1e-5);
//! ```
pub use core::PairHMM;
pub use homopolypairhmm::{BaseSpecificHopParameters, HomopolyPairHMM, HopParameters};

use crate::stats::LogProb;

mod core;
mod homopolypairhmm;

// traits common to pairhmm implementations

/// Trait for parametrization of `PairHMM` emission behavior.
pub trait EmissionParameters {
    /// Emission probability for `(x[i], y[j])`.
    /// Returns a tuple with probability and a boolean indicating whether emissions match
    /// (e.g., are the same DNA alphabet letter).
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission;

    /// Emission probability for `(x[i], -)`.
    fn prob_emit_x(&self, i: usize) -> LogProb;

    /// Emission probability for `(-, y[j])`.
    fn prob_emit_y(&self, j: usize) -> LogProb;

    fn len_x(&self) -> usize;

    fn len_y(&self) -> usize;
}
/// Trait needed for the `HomopolyPairHMM`, because its implementation details
/// depend on the actual bases to distinguish between Match states.
pub trait Emission {
    /// Base emitted at `i` in sequence `x`.
    /// Should be one of b'A', b'C', b'G' or b'T'.
    fn emission_x(&self, i: usize) -> u8;
    /// Base emitted at `i` in sequence `y`.
    /// Should be one of b'A', b'C', b'G' or b'T'.
    fn emission_y(&self, j: usize) -> u8;
}

/// Trait for parametrization of `PairHMM` gap behavior.
pub trait GapParameters {
    /// Probability to open gap in x.
    fn prob_gap_x(&self) -> LogProb;

    /// Probability to open gap in y.
    fn prob_gap_y(&self) -> LogProb;

    /// Probability to extend gap in x.
    fn prob_gap_x_extend(&self) -> LogProb;

    /// Probability to extend gap in y.
    fn prob_gap_y_extend(&self) -> LogProb;
}

/// Trait for parametrization of `PairHMM` start and end gap behavior.
/// This trait can be used to implement global and semiglobal alignments.
///
/// * global: methods return `false` and `LogProb::ln_zero()`.
/// * semiglobal: methods return `true` and `LogProb::ln_one()`.
pub trait StartEndGapParameters {
    /// Probability to start at `x[i]`. This can be left unchanged if you use `free_start_gap_x` and
    /// `free_end_gap_x`.
    #[inline]
    #[allow(unused_variables)]
    fn prob_start_gap_x(&self, i: usize) -> LogProb {
        if self.free_start_gap_x() {
            LogProb::ln_one()
        } else {
            // For global alignment, this has to return 0.0.
            LogProb::ln_zero()
        }
    }

    /// Allow free start gap in x.
    fn free_start_gap_x(&self) -> bool;

    /// Allow free end gap in x.
    fn free_end_gap_x(&self) -> bool;
}

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub enum XYEmission {
    Match(LogProb),
    Mismatch(LogProb),
}

impl XYEmission {
    pub fn prob(&self) -> LogProb {
        match *self {
            XYEmission::Match(p) => p,
            XYEmission::Mismatch(p) => p,
        }
    }

    pub fn is_match(&self) -> bool {
        match *self {
            XYEmission::Match(_) => true,
            XYEmission::Mismatch(_) => false,
        }
    }
}
