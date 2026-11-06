// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Utilities for Bayesian statistics.

pub mod bayes_factors;
pub mod model;
pub use self::bayes_factors::BayesFactor;
pub use self::model::Model;

use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::stats::LogProb;

/// For each of the hypothesis tests given as posterior error probabilities
/// (PEPs, i.e. the posterior probability of the null hypothesis), estimate the FDR
/// for the case that all null hypotheses with at most this PEP are rejected.
/// FDR is calculated as presented by Müller, Parmigiani, and Rice,
/// "FDR and Bayesian Multiple Comparisons Rules" (July 2006).
/// Johns Hopkin's University, Dept. of Biostatistics Working Papers. Working Paper 115.
///
/// # Returns
///
/// A vector of expected FDRs in the same order as the given PEPs.
///
/// # Example
///
/// ```
/// use approx::assert_relative_eq;
/// use bio::stats::bayesian::expected_fdr;
/// use bio::stats::{LogProb, Prob};
///
/// // Three hypothesis tests with posterior error probabilities (PEPs) of
/// // 0.1, 0.0 (i.e. certainly not null) and 0.25.
/// let peps = [
///     LogProb::from(Prob(0.1)),
///     LogProb::ln_zero(),
///     LogProb::from(Prob(0.25)),
/// ];
/// let fdrs = expected_fdr(&peps);
///
/// // The hypothesis with PEP = 0.0 is never expected to be a false positive.
/// assert_relative_eq!(*fdrs[1], *LogProb::ln_zero());
/// // Rejecting just the two most confident hypotheses (indices 1 and 0) yields
/// // an expected FDR of (0.0 + 0.1) / 2 = 0.05.
/// assert_relative_eq!(*fdrs[0], *LogProb::from(Prob(0.05)), epsilon = 1e-6);
/// // Rejecting all three yields an expected FDR of (0.0 + 0.1 + 0.25) / 3.
/// assert_relative_eq!(*fdrs[2], *LogProb::from(Prob(0.35 / 3.0)), epsilon = 1e-6);
/// ```
///
/// # Complexity
///
/// Runs in _O(n log n)_ time and requires _O(n)_ additional space, where
/// `n = peps.len()`, dominated by sorting the PEPs by ascending probability.
pub fn expected_fdr(peps: &[LogProb]) -> Vec<LogProb> {
    // sort indices
    let sorted_idx =
        (0..peps.len()).sorted_by(|&i, &j| OrderedFloat(*peps[i]).cmp(&OrderedFloat(*peps[j])));
    // estimate FDR
    let mut expected_fdr = vec![LogProb::ln_zero(); peps.len()];
    for (j, (expected_fp, i)) in LogProb::ln_cumsum_exp(sorted_idx.clone().map(|i| peps[i]))
        .zip(sorted_idx)
        .enumerate()
    {
        let fdr = LogProb(*expected_fp - ((j + 1) as f64).ln());
        expected_fdr[i] = if fdr <= LogProb::ln_one() {
            fdr
        } else {
            LogProb::ln_one()
        };
    }

    expected_fdr
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::LogProb;

    #[test]
    fn test_expected_fdr() {
        let peps = [
            LogProb(0.1f64.ln()),
            LogProb::ln_zero(),
            LogProb(0.25f64.ln()),
        ];
        let fdrs = expected_fdr(&peps);
        println!("{:?}", fdrs);

        assert_relative_eq!(*fdrs[1], *LogProb::ln_zero());
        assert_relative_eq!(*fdrs[0], *LogProb(0.05f64.ln()));
        assert_relative_eq!(*fdrs[2], *LogProb((0.35 / 3.0f64).ln()), epsilon = 0.000001);
    }
}
