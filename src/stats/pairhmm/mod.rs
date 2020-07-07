pub use homopolypairhmm::{HomopolyPairHMM, HopParameters};
pub use pairhmm::PairHMM;

use crate::stats::LogProb;

mod homopolypairhmm;
mod pairhmm;

// traits common to pairhmm implementations

/// Trait for parametrization of `PairHMM` emission behavior.
pub trait EmissionParameters {
    /// Emission probability for (x[i], y[j]).
    /// Returns a tuple with probability and a boolean indicating whether emissions match
    /// (e.g., are the same DNA alphabet letter).
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission;

    /// Emission probability for (x[i], -).
    fn prob_emit_x(&self, i: usize) -> LogProb;

    /// Emission probability for (-, y[j]).
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
    /// Probability to start at x[i]. This can be left unchanged if you use `free_start_gap_x` and
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

#[derive(Debug)]
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
