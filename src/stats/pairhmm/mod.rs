use crate::stats::LogProb;
pub use homopolypairhmm::EmissionParameters as ExtendedEmissionParameters;
pub use homopolypairhmm::{HomopolyPairHMM, HopParameters};
pub use pairhmm::{EmissionParameters, PairHMM};

mod homopolypairhmm;
mod pairhmm;

// traits common to pairhmm implementations

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

pub enum XYEmission {
    Match(LogProb),
    Mismatch(LogProb),
}

impl XYEmission {
    pub fn prob(&self) -> LogProb {
        match self {
            &XYEmission::Match(p) => p,
            &XYEmission::Mismatch(p) => p,
        }
    }

    pub fn is_match(&self) -> bool {
        match self {
            &XYEmission::Match(_) => true,
            &XYEmission::Mismatch(_) => false,
        }
    }
}
