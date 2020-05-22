pub use homopolypairhmm::EmissionParameters as ExtendedEmissionParameters;
pub use homopolypairhmm::{HomopolyPairHMM, HopParameters};
pub use pairhmm::{EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission};

mod homopolypairhmm;
mod pairhmm;
