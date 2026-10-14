// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! The pair HMM forward algorithm of [`PairHMM`](super::PairHMM) computed in linear
//! probability space.
//!
//! `PairHMM` keeps its three state matrices as log probabilities, which costs three log-space
//! additions (each an `exp` and a `log1p`) per cell. For the short sequences the banded
//! algorithm is typically used on, plain `f64` probabilities are precise enough, and a cell
//! then costs a few multiplications and additions. Every cell inside the band is reached by a
//! path with at most `max_edit_dist` edits, so its value is bounded below by the product of
//! that path; only when even that product is below the smallest normal `f64` does the result
//! underflow, in which case [`LinearPairHMM::prob_related`] returns `None` and the caller can
//! fall back to `PairHMM`. Terms far below the dominant path may underflow earlier, which is
//! harmless: they are negligible in the sum anyway.
//!
//! The recurrences, the banding by minimal edit distance and the handling of the alignment
//! mode are the same as in `PairHMM`, so both compute the same probability up to the
//! approximations `PairHMM` uses in log space.

use std::cmp;
use std::mem;

use serde::{Deserialize, Serialize};

use crate::stats::pairhmm::band::BandScan;
use crate::stats::pairhmm::{EmissionParameters, GapParameters, StartEndGapParameters};
use crate::stats::LogProb;

/// Emission parameters of a pair HMM in linear probability space.
///
/// Implement this directly when the probabilities can be precomputed (e.g. from base
/// qualities), or wrap any [`EmissionParameters`] in [`LogSpaceEmission`] to derive them from
/// the log-space values at the cost of an `exp` per cell.
pub trait LinearEmissionParameters {
    /// Probability of emitting `x[i]` and `y[j]` together, and whether this is a match.
    fn prob_emit_xy(&self, i: usize, j: usize) -> (f64, bool);

    /// Probability of emitting `x[i]` alone.
    fn prob_emit_x(&self, i: usize) -> f64;

    /// Probability of emitting `y[j]` alone.
    fn prob_emit_y(&self, j: usize) -> f64;

    /// Length of x.
    fn len_x(&self) -> usize;

    /// Length of y.
    fn len_y(&self) -> usize;
}

/// [`LinearEmissionParameters`] derived from any [`EmissionParameters`] by exponentiation.
#[derive(Clone, Copy, Debug)]
pub struct LogSpaceEmission<'a, E>(pub &'a E);

impl<E: EmissionParameters> LinearEmissionParameters for LogSpaceEmission<'_, E> {
    #[inline]
    fn prob_emit_xy(&self, i: usize, j: usize) -> (f64, bool) {
        let emission = self.0.prob_emit_xy(i, j);
        (emission.prob().exp(), emission.is_match())
    }

    #[inline]
    fn prob_emit_x(&self, i: usize) -> f64 {
        self.0.prob_emit_x(i).exp()
    }

    #[inline]
    fn prob_emit_y(&self, j: usize) -> f64 {
        self.0.prob_emit_y(j).exp()
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.0.len_x()
    }

    #[inline]
    fn len_y(&self) -> usize {
        self.0.len_y()
    }
}

#[derive(Default, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
struct LinearGapParamCache {
    no_gap: f64,
    no_gap_x_extend: f64,
    no_gap_y_extend: f64,
    gap_x: f64,
    gap_y: f64,
    gap_x_extend: f64,
    gap_y_extend: f64,
    do_gap_x_extend: bool,
    do_gap_y_extend: bool,
}

/// A pair HMM whose forward algorithm runs in linear probability space.
///
/// It computes the same quantity as [`PairHMM`](super::PairHMM) from plain `f64`
/// probabilities, which is several times faster, and reports an underflow instead of a
/// result (see [`LinearPairHMM::prob_related`]).
#[derive(Default, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub struct LinearPairHMM {
    // state M (match or mismatch), X (x[i] emitted alone) and Y (y[j] emitted alone), as two
    // alternating columns indexed by position in y + 1
    fm: [Vec<f64>; 2],
    fx: [Vec<f64>; 2],
    fy: [Vec<f64>; 2],
    min_edit_dist: [Vec<usize>; 2],
    scan: BandScan,
    // probability of the alignments ending in each column (free end gap in x)
    prob_cols: Vec<f64>,
    gap_params: LinearGapParamCache,
}

impl LinearPairHMM {
    /// Create a new instance with the given gap parameters.
    pub fn new<G>(gap_params: &G) -> Self
    where
        G: GapParameters,
    {
        let lin = |p: LogProb| p.exp();
        let gap_params = LinearGapParamCache {
            no_gap: lin(gap_params
                .prob_gap_x()
                .ln_add_exp(gap_params.prob_gap_y())
                .ln_one_minus_exp()),
            no_gap_x_extend: lin(gap_params.prob_gap_x_extend().ln_one_minus_exp()),
            no_gap_y_extend: lin(gap_params.prob_gap_y_extend().ln_one_minus_exp()),
            gap_x: lin(gap_params.prob_gap_x()),
            gap_y: lin(gap_params.prob_gap_y()),
            gap_x_extend: lin(gap_params.prob_gap_x_extend()),
            gap_y_extend: lin(gap_params.prob_gap_y_extend()),
            do_gap_y_extend: gap_params.prob_gap_y_extend() != LogProb::ln_zero(),
            do_gap_x_extend: gap_params.prob_gap_x_extend() != LogProb::ln_zero(),
        };
        Self {
            gap_params,
            ..Default::default()
        }
    }

    /// Calculate the probability of sequence x being related to y via any alignment, like
    /// [`PairHMM::prob_related`](super::PairHMM::prob_related), or `None` if the result is not
    /// a normal `f64`: either the computation underflowed, or no alignment exists at all
    /// (e.g. within the band, in global mode). The caller is expected to fall back to the
    /// log-space implementation, which tells the two cases apart.
    ///
    /// # Arguments
    ///
    /// * `emission_params` - parameters for emission, in linear space
    /// * `alignment_mode` - parameters for free end/start gaps
    /// * `max_edit_dist` - maximum edit distance to consider; if not `None`, perform banded
    ///   alignment
    pub fn prob_related<E, A>(
        &mut self,
        emission_params: &E,
        alignment_mode: &A,
        max_edit_dist: Option<usize>,
    ) -> Option<LogProb>
    where
        E: LinearEmissionParameters,
        A: StartEndGapParameters,
    {
        let len_x = emission_params.len_x();
        let len_y = emission_params.len_y();

        for k in 0..2 {
            self.fm[k].clear();
            self.fx[k].clear();
            self.fy[k].clear();
            self.min_edit_dist[k].clear();
            self.fm[k].resize(len_y + 1, 0.0);
            self.fx[k].resize(len_y + 1, 0.0);
            self.fy[k].resize(len_y + 1, 0.0);
            self.min_edit_dist[k].resize(len_y + 1, usize::MAX);
        }
        self.prob_cols.clear();
        if alignment_mode.free_end_gap_x() {
            self.prob_cols.reserve(len_x);
        }
        self.scan.reset();

        let mut prev = 0;
        let mut curr = 1;
        self.fm[prev][0] = 1.0;
        // origin cell, see PairHMM
        self.min_edit_dist[prev][0] = 0;

        for i in 0..len_x {
            // Clear what this column buffer holds from two columns ago, so that every cell
            // outside of the ranges written below is empty.
            for &(a, b) in self.scan.stale(curr) {
                self.fm[curr][a..=b].fill(0.0);
                self.fx[curr][a..=b].fill(0.0);
                self.fy[curr][a..=b].fill(0.0);
                self.min_edit_dist[curr][a..=b].fill(usize::MAX);
            }
            self.fm[curr][0] = 0.0;
            self.min_edit_dist[curr][0] = usize::MAX;

            // allow alignment to start from offset in x (if prob_start_gap_x is set accordingly)
            self.fm[prev][0] += alignment_mode.prob_start_gap_x(i).exp();
            // With a band, an alignment starting at column i has to consume the whole of y
            // within the remaining columns, which needs at least len_y - (len_x - i)
            // insertions: beyond that it cannot stay inside the band and would only be
            // computed to be pruned.
            let may_start = match max_edit_dist {
                Some(d) => i <= (len_x + d).saturating_sub(len_y),
                None => true,
            };
            if alignment_mode.free_start_gap_x() {
                if may_start {
                    self.min_edit_dist[prev][0] = 0;
                } else {
                    self.fm[prev][0] = 0.0;
                    self.min_edit_dist[prev][0] = usize::MAX;
                }
            }

            let prob_emit_x = emission_params.prob_emit_x(i);

            self.scan
                .begin_column(curr, prev, len_y, self.min_edit_dist[prev][0] != usize::MAX);
            let mut inside = false;
            while let Some(j_) = self.scan.next_cell(inside) {
                let j = j_ - 1;

                let min_edit_dist_topleft = self.min_edit_dist[prev][j];
                let min_edit_dist_top = self.min_edit_dist[curr][j];
                let min_edit_dist_left = self.min_edit_dist[prev][j_];

                inside = match max_edit_dist {
                    Some(max_edit_dist) => {
                        cmp::min(
                            min_edit_dist_topleft,
                            cmp::min(min_edit_dist_top, min_edit_dist_left),
                        ) <= max_edit_dist
                    }
                    None => true,
                };
                if !inside {
                    continue;
                }

                let (emit_xy, is_match) = emission_params.prob_emit_xy(i, j);

                // match or mismatch, coming from M, X (extended with gap_y_extend) or Y
                // (extended with gap_x_extend)
                let prob_match_mismatch = emit_xy
                    * (self.gap_params.no_gap * self.fm[prev][j]
                        + self.gap_params.no_gap_y_extend * self.fx[prev][j]
                        + self.gap_params.no_gap_x_extend * self.fy[prev][j]);

                // gap in y: x[i] emitted alone, opened or extended
                let mut prob_gap_y = self.gap_params.gap_y * self.fm[prev][j_];
                if self.gap_params.do_gap_y_extend {
                    prob_gap_y += self.gap_params.gap_y_extend * self.fx[prev][j_];
                }
                prob_gap_y *= prob_emit_x;

                // gap in x: y[j] emitted alone, opened or extended
                let mut prob_gap_x = self.gap_params.gap_x * self.fm[curr][j];
                if self.gap_params.do_gap_x_extend {
                    prob_gap_x += self.gap_params.gap_x_extend * self.fy[curr][j];
                }
                prob_gap_x *= emission_params.prob_emit_y(j);

                self.fm[curr][j_] = prob_match_mismatch;
                self.fx[curr][j_] = prob_gap_y;
                self.fy[curr][j_] = prob_gap_x;
                if max_edit_dist.is_some() {
                    self.min_edit_dist[curr][j_] = cmp::min(
                        if is_match {
                            min_edit_dist_topleft
                        } else {
                            min_edit_dist_topleft.saturating_add(1)
                        },
                        cmp::min(
                            min_edit_dist_left.saturating_add(1),
                            min_edit_dist_top.saturating_add(1),
                        ),
                    );
                }
                self.scan.mark(curr, j_);
            }

            if alignment_mode.free_end_gap_x() {
                self.prob_cols
                    .push(self.fm[curr][len_y] + self.fx[curr][len_y] + self.fy[curr][len_y]);
            }
            mem::swap(&mut curr, &mut prev);
        }

        let p = if alignment_mode.free_end_gap_x() {
            self.prob_cols.iter().sum::<f64>()
        } else {
            self.fm[prev][len_y] + self.fx[prev][len_y] + self.fy[prev][len_y]
        };
        if p.is_nan() || p < f64::MIN_POSITIVE {
            // underflow (or no alignment at all): leave it to the log-space implementation
            return None;
        }
        // the sum over all paths can exceed 1.0, especially in case of repeats
        Some(if p > 1.0 {
            LogProb::ln_one()
        } else {
            LogProb(p.ln())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::pairhmm::{PairHMM, XYEmission};
    use crate::stats::Prob;
    use rand::prelude::*;

    struct Gaps {
        extend: bool,
        free: bool,
    }

    impl GapParameters for Gaps {
        fn prob_gap_x(&self) -> LogProb {
            LogProb::from(Prob(1e-3))
        }

        fn prob_gap_y(&self) -> LogProb {
            LogProb::from(Prob(2e-3))
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            if self.extend {
                LogProb::from(Prob(0.3))
            } else {
                LogProb::ln_zero()
            }
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            if self.extend {
                LogProb::from(Prob(0.4))
            } else {
                LogProb::ln_zero()
            }
        }
    }

    impl StartEndGapParameters for Gaps {
        fn free_start_gap_x(&self) -> bool {
            self.free
        }

        fn free_end_gap_x(&self) -> bool {
            self.free
        }
    }

    /// Emissions with per-position base qualities of y and per-position emission
    /// probabilities of x (so that a dropped emission in the X state is visible).
    struct Emission {
        x: Vec<u8>,
        y: Vec<u8>,
        miscall: Vec<f64>,
        emit_x: Vec<f64>,
    }

    impl EmissionParameters for Emission {
        fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
            if self.x[i] == self.y[j] {
                XYEmission::Match(LogProb::from(Prob(1.0 - self.miscall[j])))
            } else {
                XYEmission::Mismatch(LogProb::from(Prob(self.miscall[j] / 3.0)))
            }
        }

        fn prob_emit_x(&self, i: usize) -> LogProb {
            LogProb::from(Prob(self.emit_x[i]))
        }

        fn prob_emit_y(&self, j: usize) -> LogProb {
            LogProb::from(Prob(self.miscall[j]))
        }

        fn len_x(&self) -> usize {
            self.x.len()
        }

        fn len_y(&self) -> usize {
            self.y.len()
        }
    }

    impl LinearEmissionParameters for Emission {
        fn prob_emit_xy(&self, i: usize, j: usize) -> (f64, bool) {
            if self.x[i] == self.y[j] {
                (1.0 - self.miscall[j], true)
            } else {
                (self.miscall[j] / 3.0, false)
            }
        }

        fn prob_emit_x(&self, i: usize) -> f64 {
            self.emit_x[i]
        }

        fn prob_emit_y(&self, j: usize) -> f64 {
            self.miscall[j]
        }

        fn len_x(&self) -> usize {
            self.x.len()
        }

        fn len_y(&self) -> usize {
            self.y.len()
        }
    }

    fn random_emission(rng: &mut StdRng, max_edits: usize) -> (Emission, usize) {
        let len = rng.random_range(20..120);
        let x: Vec<u8> = (0..len).map(|_| *b"ACGT".choose(rng).unwrap()).collect();
        let start = rng.random_range(0..len / 3);
        let end = rng.random_range(start + 10..=len);
        let mut y = x[start..end].to_vec();
        let edits = rng.random_range(0..=max_edits);
        for _ in 0..edits {
            let pos = rng.random_range(0..y.len());
            match rng.random_range(0..3) {
                0 => y[pos] = *b"ACGT".choose(rng).unwrap(),
                1 => y.insert(pos, *b"ACGT".choose(rng).unwrap()),
                _ => {
                    y.remove(pos);
                }
            }
        }
        let miscall = (0..y.len())
            .map(|_| 10f64.powf(-(rng.random_range(2..=40) as f64) / 10.0))
            .collect();
        let emit_x = (0..x.len()).map(|_| rng.random_range(0.05..1.0)).collect();
        (
            Emission {
                x,
                y,
                miscall,
                emit_x,
            },
            edits,
        )
    }

    /// Both forward algorithms compute the same probability, in global and semiglobal mode,
    /// banded and unbanded, with and without gap extension. The log-space one is the
    /// approximate one (fastexp and a three-way log-sum that drops small terms).
    #[test]
    fn test_agrees_with_log_space() {
        let mut rng = StdRng::seed_from_u64(3);
        for extend in [false, true] {
            for free in [false, true] {
                let gaps = Gaps { extend, free };
                let mut linear = LinearPairHMM::new(&gaps);
                let mut log_space = PairHMM::new(&gaps);
                for _ in 0..150 {
                    let (emission, edits) = random_emission(&mut rng, 6);
                    for max_edit_dist in [None, Some(edits + 4)] {
                        let expected = log_space.prob_related(&emission, &gaps, max_edit_dist);
                        if !expected.is_finite() {
                            // no alignment within the band (global mode): reported as None
                            assert_eq!(linear.prob_related(&emission, &gaps, max_edit_dist), None);
                            continue;
                        }
                        let p = linear
                            .prob_related(&emission, &gaps, max_edit_dist)
                            .expect("no underflow expected on these sizes");
                        assert!(
                            (*p - *expected).abs() < 1e-3,
                            "linear {} vs log-space {} (extend {}, free {}, band {:?})",
                            *p,
                            *expected,
                            extend,
                            free,
                            max_edit_dist
                        );
                        let via_adapter = linear
                            .prob_related(&LogSpaceEmission(&emission), &gaps, max_edit_dist)
                            .unwrap();
                        assert!((*p - *via_adapter).abs() < 1e-9);
                    }
                }
            }
        }
    }

    /// A long window of mismatches is below the smallest normal f64: the log-space forward
    /// still has a value, the linear one reports the underflow.
    #[test]
    fn test_underflow_is_reported() {
        let gaps = Gaps {
            extend: false,
            free: false,
        };
        let emission = Emission {
            x: vec![b'A'; 300],
            y: vec![b'C'; 300],
            miscall: vec![1e-4; 300],
            emit_x: vec![1.0; 300],
        };
        assert!(PairHMM::new(&gaps)
            .prob_related(&emission, &gaps, None)
            .is_finite());
        assert_eq!(
            LinearPairHMM::new(&gaps).prob_related(&emission, &gaps, None),
            None
        );
    }

    /// `None` starts where the result leaves the normal range of `f64`, not at zero: a
    /// subnormal result is an underflow too.
    #[test]
    fn test_underflow_threshold_is_the_smallest_normal() {
        let gaps = Gaps {
            extend: false,
            free: false,
        };
        let mut log_space = PairHMM::new(&gaps);
        let mut linear = LinearPairHMM::new(&gaps);
        let ln_min_positive = f64::MIN_POSITIVE.ln();
        let (mut seen_some, mut seen_none) = (false, false);
        for n in (950..1080).step_by(10) {
            let emission = Emission {
                x: vec![b'A'; n],
                y: vec![b'A'; n],
                miscall: vec![0.5; n],
                emit_x: vec![1.0; n],
            };
            let expected = *log_space.prob_related(&emission, &gaps, None);
            let p = linear.prob_related(&emission, &gaps, None);
            if expected > ln_min_positive + 0.05 {
                let p = p.expect("a normal result");
                assert!((*p - expected).abs() < 1e-2, "{} vs {}", *p, expected);
                seen_some = true;
            } else if expected < ln_min_positive - 0.05 {
                assert_eq!(p, None, "subnormal result {} must be reported", expected);
                seen_none = true;
            }
        }
        assert!(seen_some && seen_none);
    }

    /// An empty y has no alignment: zero in log space, `None` here.
    #[test]
    fn test_empty_y_is_none() {
        for free in [false, true] {
            let gaps = Gaps {
                extend: false,
                free,
            };
            let emission = Emission {
                x: b"ACGT".to_vec(),
                y: vec![],
                miscall: vec![],
                emit_x: vec![1.0; 4],
            };
            for max_edit_dist in [None, Some(0), Some(2)] {
                assert_eq!(
                    PairHMM::new(&gaps).prob_related(&emission, &gaps, max_edit_dist),
                    LogProb::ln_zero()
                );
                assert_eq!(
                    LinearPairHMM::new(&gaps).prob_related(&emission, &gaps, max_edit_dist),
                    None
                );
            }
        }
    }
}
