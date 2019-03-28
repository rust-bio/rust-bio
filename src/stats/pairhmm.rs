// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A pair Hidden Markov Model for calculating the probability that two sequences are related to
//! each other. Depending on the used parameters, this can, e.g., be used to calculate the
//! probability that a certain sequencing read comes from a given position in a reference genome.

use std::cmp;
use std::mem;
use std::usize;

use crate::stats::LogProb;

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

/// Fast approximation of sum over the three given proabilities. If the largest is sufficiently
/// large compared to the others, we just return that instead of computing the full (expensive)
/// sum.
#[inline]
fn ln_sum3_exp_approx(mut p0: LogProb, mut p1: LogProb, mut p2: LogProb) -> LogProb {
    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    if p1 > p0 {
        mem::swap(&mut p1, &mut p0);
    }
    if *(p0 - p1) > 10.0 {
        // if p0 is strong enough compared to second, just return the maximum
        p0
    } else {
        // calculate accurate sum (slower)
        LogProb::ln_sum_exp(&[p0, p1, p2])
    }
}

/// A pair Hidden Markov Model for comparing sequences x and y as described by
/// Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
/// Current Topics in Genome Analysis 2008. http://doi.org/10.1017/CBO9780511790492.
pub struct PairHMM {
    fm: [Vec<LogProb>; 2],
    fx: [Vec<LogProb>; 2],
    fy: [Vec<LogProb>; 2],
    min_edit_dist: [Vec<usize>; 2],
    prob_cols: Vec<LogProb>,
}

impl Default for PairHMM {
    fn default() -> Self {
        PairHMM {
            fm: [Vec::new(), Vec::new()],
            fx: [Vec::new(), Vec::new()],
            fy: [Vec::new(), Vec::new()],
            min_edit_dist: [Vec::new(), Vec::new()],
            prob_cols: Vec::new(),
        }
    }
}

impl PairHMM {
    pub fn new() -> Self {
        Default::default()
    }

    /// Calculate the probability of sequence x being related to y via any alignment.
    ///
    /// # Arguments
    ///
    /// * `gap_params` - parameters for opening or extending gaps
    /// * `emission_params` - parameters for emission
    /// * `max_edit_dist` - maximum edit distance to consider; if not `None`, perform banded alignment
    pub fn prob_related<G, E>(
        &mut self,
        gap_params: &G,
        emission_params: &E,
        max_edit_dist: Option<usize>,
    ) -> LogProb
    where
        G: GapParameters + StartEndGapParameters,
        E: EmissionParameters,
    {
        for k in 0..2 {
            self.fm[k].clear();
            self.fx[k].clear();
            self.fy[k].clear();
            self.min_edit_dist[k].clear();
            self.prob_cols.clear();

            self.fm[k].resize(emission_params.len_y() + 1, LogProb::ln_zero());
            self.fx[k].resize(emission_params.len_y() + 1, LogProb::ln_zero());
            self.fy[k].resize(emission_params.len_y() + 1, LogProb::ln_zero());
            self.min_edit_dist[k].resize(emission_params.len_y() + 1, usize::MAX);

            if gap_params.free_end_gap_x() {
                let c = (emission_params.len_x() * 3).saturating_sub(self.prob_cols.capacity());
                self.prob_cols.reserve_exact(c);
            }
        }

        // cache probs
        let prob_no_gap = gap_params
            .prob_gap_x()
            .ln_add_exp(gap_params.prob_gap_y())
            .ln_one_minus_exp();
        let prob_no_gap_x_extend = gap_params.prob_gap_x_extend().ln_one_minus_exp();
        let prob_no_gap_y_extend = gap_params.prob_gap_y_extend().ln_one_minus_exp();
        let prob_gap_x = gap_params.prob_gap_x();
        let prob_gap_y = gap_params.prob_gap_y();
        let prob_gap_x_extend = gap_params.prob_gap_x_extend();
        let prob_gap_y_extend = gap_params.prob_gap_y_extend();
        let do_gap_y_extend = prob_gap_y_extend != LogProb::ln_zero();
        let do_gap_x_extend = prob_gap_x_extend != LogProb::ln_zero();

        let mut prev = 0;
        let mut curr = 1;
        self.fm[prev][0] = LogProb::ln_one();

        // iterate over x
        for i in 0..emission_params.len_x() {
            // allow alignment to start from offset in x (if prob_start_gap_x is set accordingly)
            self.fm[prev][0] = self.fm[prev][0].ln_add_exp(gap_params.prob_start_gap_x(i));
            if gap_params.free_start_gap_x() {
                self.min_edit_dist[prev][0] = 0;
            }

            let prob_emit_x = emission_params.prob_emit_x(i);

            let j_min = if do_gap_x_extend {
                0
            } else {
                // The final upper triangle in the dynamic programming matrix
                // will be only zeros if gap extension is not allowed on x.
                // Hence we can omit it.
                emission_params
                    .len_y()
                    .saturating_sub(emission_params.len_x() - i + 1)
            };

            let j_max = if do_gap_x_extend {
                emission_params.len_y()
            } else {
                // The initial lower triangle in the dynamic programming matrix
                // will be only zeros if gap extension is not allowed on x.
                // Hence we can omit it.
                cmp::min(i + 1, emission_params.len_y())
            };

            // iterate over y
            for j in j_min..j_max {
                let j_ = j + 1;
                let j_minus_one = j_ - 1;

                let min_edit_dist_topleft = self.min_edit_dist[prev][j_minus_one];
                let min_edit_dist_top = self.min_edit_dist[curr][j_minus_one];
                let min_edit_dist_left = self.min_edit_dist[prev][j_];

                if let Some(max_edit_dist) = max_edit_dist {
                    if cmp::min(
                        min_edit_dist_topleft,
                        cmp::min(min_edit_dist_top, min_edit_dist_left),
                    ) > max_edit_dist
                    {
                        // skip this cell if best edit dist is already larger than given maximum
                        continue;
                    }
                }

                let (prob_match_mismatch, prob_gap_x, prob_gap_y, min_edit_dist) = {
                    let fm_curr = &self.fm[curr];
                    let fm_prev = &self.fm[prev];
                    let fx_prev = &self.fx[prev];
                    let fy_curr = &self.fy[curr];
                    let fy_prev = &self.fy[prev];

                    // match or mismatch
                    let emit_xy = emission_params.prob_emit_xy(i, j);
                    let prob_match_mismatch = emit_xy.prob()
                        + ln_sum3_exp_approx(
                            prob_no_gap + fm_prev[j_minus_one],
                            // coming from state X
                            prob_no_gap_x_extend + fx_prev[j_minus_one],
                            // coming from state Y
                            prob_no_gap_y_extend + fy_prev[j_minus_one],
                        );

                    // gap in y
                    let mut prob_gap_y = prob_emit_x
                        + (
                            // open gap
                            prob_gap_y + fm_prev[j_]
                        );
                    if do_gap_y_extend {
                        prob_gap_y = prob_gap_y.ln_add_exp(
                            // extend gap
                            prob_gap_y_extend + fx_prev[j_],
                        );
                    }

                    // gap in x
                    let mut prob_gap_x = emission_params.prob_emit_y(j)
                        + (
                            // open gap
                            prob_gap_x + fm_curr[j_minus_one]
                        );
                    if do_gap_x_extend {
                        prob_gap_x = prob_gap_x.ln_add_exp(
                            // extend gap
                            prob_gap_x_extend + fy_curr[j_minus_one],
                        );
                    }

                    // calculate minimal number of mismatches
                    let min_edit_dist = if max_edit_dist.is_some() {
                        cmp::min(
                            if emit_xy.is_match() {
                                // a match, so nothing changes
                                min_edit_dist_topleft
                            } else {
                                // one new mismatch
                                min_edit_dist_topleft.saturating_add(1)
                            },
                            cmp::min(
                                // gap in y (no new mismatch)
                                min_edit_dist_left.saturating_add(1),
                                // gap in x (no new mismatch)
                                min_edit_dist_top.saturating_add(1),
                            ),
                        )
                    } else {
                        0
                    };

                    (prob_match_mismatch, prob_gap_x, prob_gap_y, min_edit_dist)
                };

                self.fm[curr][j_] = prob_match_mismatch;
                self.fx[curr][j_] = prob_gap_y;
                self.fy[curr][j_] = prob_gap_x;
                if max_edit_dist.is_some() {
                    self.min_edit_dist[curr][j_] = min_edit_dist;
                }
            }

            if gap_params.free_end_gap_x() {
                // Cache column probabilities or simply record the last probability.
                // We can put all of them in one array since we simply have to sum in the end.
                // This is also good for numerical stability.
                self.prob_cols.push(self.fm[curr].last().unwrap().clone());
                self.prob_cols.push(self.fx[curr].last().unwrap().clone());
                // TODO check removing this (we don't want open gaps in x):
                self.prob_cols.push(self.fy[curr].last().unwrap().clone());
            }

            // next column
            mem::swap(&mut curr, &mut prev);
            // reset next column to zeros
            for v in &mut self.fm[curr] {
                *v = LogProb::ln_zero();
            }
        }

        let p = if gap_params.free_end_gap_x() {
            LogProb::ln_sum_exp(&self.prob_cols)
        } else {
            LogProb::ln_sum_exp(&[
                *self.fm[prev].last().unwrap(),
                *self.fx[prev].last().unwrap(),
                *self.fy[prev].last().unwrap(),
            ])
        };
        // take the minimum with 1.0, because sum of paths can exceed probability 1.0
        // especially in case of repeats
        assert!(!p.is_nan());
        if p > LogProb::ln_one() {
            LogProb::ln_one()
        } else {
            p
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::{LogProb, Prob};

    // Single base insertion and deletion rates for R1 according to Schirmer et al.
    // BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
    static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
    static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
    static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

    fn prob_emit_x_or_y() -> LogProb {
        LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
    }

    struct TestEmissionParams {
        x: &'static [u8],
        y: &'static [u8],
    }

    impl EmissionParameters for TestEmissionParams {
        fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
            if self.x[i] == self.y[j] {
                XYEmission::Match(LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST))
            } else {
                XYEmission::Mismatch(LogProb::from(PROB_ILLUMINA_SUBST / Prob(3.0)))
            }
        }

        fn prob_emit_x(&self, _: usize) -> LogProb {
            prob_emit_x_or_y()
        }

        fn prob_emit_y(&self, _: usize) -> LogProb {
            prob_emit_x_or_y()
        }

        fn len_x(&self) -> usize {
            self.x.len()
        }

        fn len_y(&self) -> usize {
            self.y.len()
        }
    }

    struct TestGapParams;

    impl GapParameters for TestGapParams {
        fn prob_gap_x(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_INS)
        }

        fn prob_gap_y(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_DEL)
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }
    }

    impl StartEndGapParameters for TestGapParams {
        fn free_start_gap_x(&self) -> bool {
            false
        }

        fn free_end_gap_x(&self) -> bool {
            false
        }
    }

    pub struct SemiglobalGapParams;

    impl GapParameters for SemiglobalGapParams {
        fn prob_gap_x(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_INS)
        }

        fn prob_gap_y(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_DEL)
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }
    }

    impl StartEndGapParameters for SemiglobalGapParams {
        fn free_start_gap_x(&self) -> bool {
            true
        }

        fn free_end_gap_x(&self) -> bool {
            true
        }
    }

    #[test]
    fn test_same() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p, 0.0, epsilon = 0.1);
    }

    #[test]
    fn test_insertion() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCTGATCGATCT";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(p.exp(), PROB_ILLUMINA_INS.powi(2), epsilon = 1e-11);
    }

    #[test]
    fn test_deletion() {
        let x = b"AGCTCGATCTGATCGATCT";
        let y = b"AGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(p.exp(), PROB_ILLUMINA_DEL.powi(2), epsilon = 1e-10);
    }

    #[test]
    fn test_mismatch() {
        let x = b"AGCTCGAGCGATCGATC";
        let y = b"TGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(
            p.exp(),
            (PROB_ILLUMINA_SUBST / Prob(3.0)).powi(2),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_large() {
        let x = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC\
ATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTAT\
CTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGT";

        let y = b"GGGTATGCACGCGATAGCATTGCGAGATGCTGGAGCTGGAGCACCCTATGTCGC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = SemiglobalGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let p_banded = pair_hmm.prob_related(&gap_params, &emission_params, Some(2));

        assert_relative_eq!(*p, *p_banded, epsilon = 1e-7);
    }
}
