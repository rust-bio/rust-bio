// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Utilities for Bayesian statistics.

use itertools::Itertools;
use ordered_float::OrderedFloat;

use stats::LogProb;

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
pub fn expected_fdr(peps: &[LogProb]) -> Vec<LogProb> {
    // sort indices
    let sorted_idx =
        (0..peps.len()).sorted_by(|&i, &j| OrderedFloat(*peps[i]).cmp(&OrderedFloat(*peps[j])));
    // estimate FDR
    let mut expected_fdr = vec![LogProb::ln_zero(); peps.len()];
    for (i, expected_fp) in LogProb::ln_cumsum_exp(sorted_idx.iter().map(|&i| peps[i])).enumerate()
    {
        let fdr = LogProb(*expected_fp - ((i + 1) as f64).ln());
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
    use stats::LogProb;

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
        assert_relative_eq!(*fdrs[2], *LogProb((0.35 / 3.0f64).ln()));
    }
}
