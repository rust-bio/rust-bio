// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A pair Hidden Markov Model for calculating the probability that two sequences are related to
//! each other. Depending on the used parameters, this can, e.g., be used to calculate the
//! probability that a certain sequencing read comes from a given position in a reference genome.

use std::cmp;
use std::fmt::Debug;
use std::mem;
use std::ops::{Add, AddAssign, Div, Mul, Shr, Sub};
use std::usize;

use enum_map::{Enum, EnumMap};
use itertools::Itertools;
use num_traits::{One, Zero};

use crate::stats::pairhhmm::State::*;
use crate::stats::probs::TLogProb;

/// Trait for parametrization of `PairHMM` gap behavior.
pub trait GapParameters {
    /// Probability to open gap in x.
    fn prob_gap_x(&self) -> TLogProb;

    /// Probability to open gap in y.
    fn prob_gap_y(&self) -> TLogProb;

    /// Probability to extend gap in x.
    fn prob_gap_x_extend(&self) -> TLogProb;

    /// Probability to extend gap in y.
    fn prob_gap_y_extend(&self) -> TLogProb;
}

/// Trait for parametrization of `PairHMM` gap behavior.
pub trait HopParameters {
    /// Probability to open gap in x.
    fn prob_hop_x(&self) -> TLogProb;

    /// Probability to open gap in y.
    fn prob_hop_y(&self) -> TLogProb;

    /// Probability to extend gap in x.
    fn prob_hop_x_extend(&self) -> TLogProb;

    /// Probability to extend gap in y.
    fn prob_hop_y_extend(&self) -> TLogProb;
}

/// Trait for parametrization of `PairHMM` start and end gap behavior.
/// This trait can be used to implement global and semiglobal alignments.
///
/// * global: methods return `false` and `TLogProb::ln_zero()`.
/// * semiglobal: methods return `true` and `TLogProb::ln_one()`.
pub trait StartEndGapParameters {
    /// Probability to start at x[i]. This can be left unchanged if you use `free_start_gap_x` and
    /// `free_end_gap_x`.
    #[inline]
    #[allow(unused_variables)]
    fn prob_start_gap_x(&self, i: usize) -> TLogProb {
        if self.free_start_gap_x() {
            TLogProb::one()
        } else {
            // For global alignment, this has to return 0.0.
            TLogProb::zero()
        }
    }

    /// Allow free start gap in x.
    fn free_start_gap_x(&self) -> bool;

    /// Allow free end gap in x.
    fn free_end_gap_x(&self) -> bool;
}

pub enum XYEmission {
    Match(TLogProb),
    Mismatch(TLogProb),
}

impl XYEmission {
    pub fn prob(&self) -> TLogProb {
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

#[derive(Eq, PartialEq, Debug, Enum, Clone, Copy)]
#[repr(usize)]
pub enum State {
    MatchA = 0,
    MatchC = 1,
    MatchG = 2,
    MatchT = 3,
    GapX = 4,
    GapY = 5,
    HopAX = 6,
    HopAY = 7,
    HopCX = 8,
    HopCY = 9,
    HopGX = 10,
    HopGY = 11,
    HopTX = 12,
    HopTY = 13,
}

const STATES: [State; 14] = [
    MatchA, MatchC, MatchG, MatchT, GapX, GapY, HopAX, HopAY, HopCX, HopCY, HopGX, HopGY, HopTX,
    HopTY,
];

const MATCH_STATES: [State; 4] = [MatchA, MatchC, MatchG, MatchT];
const HOP_X_STATES: [State; 4] = [HopAX, HopCX, HopGX, HopTX];
const HOP_Y_STATES: [State; 4] = [HopAY, HopCY, HopGY, HopTY];

fn space_bits(a: u32) -> u64 {
    let mut x = a as u64 & 0x0000_0000_FFFF_FFFF;
    x = (x | (x << 16)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x << 8)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x << 2)) & 0x3333_3333_3333_3333;
    x = (x | (x << 1)) & 0x5555_5555_5555_5555;
    x
}

pub fn interleave_bits(a: u32, b: u32) -> u64 {
    space_bits(a) << 1 | space_bits(b)
}

impl Shr for State {
    type Output = usize;

    fn shr(self, rhs: State) -> Self::Output {
        let a = self as u32;
        let b = rhs as u32;
        interleave_bits(a, b) as usize
    }
}

/// Trait for parametrization of `PairHMM` emission behavior.
pub trait EmissionParameters {
    /// Emission probability for (x[i], y[j]).
    /// Returns a tuple with probability and a boolean indicating whether emissions match
    /// (e.g., are the same DNA alphabet letter).
    fn prob_emit_x_and_y(&self, s: State, i: usize, j: usize) -> XYEmission;

    /// Emission probability for (x[i], -).
    fn prob_emit_x_and_gap(&self, s: State, i: usize) -> TLogProb;

    /// Emission probability for (-, y[j]).
    fn prob_emit_gap_and_y(&self, s: State, j: usize) -> TLogProb;

    /// Emission probability for (x[i], -).
    fn prob_emit_x_and_hop(&self, s: State, i: usize) -> TLogProb;

    /// Emission probability for (-, y[j]).
    fn prob_emit_hop_and_y(&self, s: State, j: usize) -> TLogProb;

    fn len_x(&self) -> usize;

    fn len_y(&self) -> usize;
}

trait Pr:
    Add<Output = Self>
    + AddAssign<Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Div<f64, Output = Self>
    + One
    + Zero
    + Clone
    + Copy
    + Debug
{
}

pub trait Reset<T: Copy> {
    fn reset(&mut self, value: T);
}

impl<T: Copy> Reset<T> for [T] {
    fn reset(&mut self, value: T) {
        for v in self {
            *v = value;
        }
    }
}

/// A pair Hidden Markov Model for comparing sequences x and y as described by
/// Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
/// Current Topics in Genome Analysis 2008. http://doi.org/10.1017/CBO9780511790492.
#[derive(Debug, Clone)]
pub struct PairHHMM {}

lazy_static! {
    static ref MATCH_HOP_X: Vec<(State, State)> = MATCH_STATES
        .iter()
        .copied()
        .zip(HOP_X_STATES.iter().copied())
        .collect_vec();
    static ref MATCH_HOP_Y: Vec<(State, State)> = MATCH_STATES
        .iter()
        .copied()
        .zip(HOP_Y_STATES.iter().copied())
        .collect_vec();
    static ref HOP_X_HOP_X: Vec<(State, State)> = HOP_X_STATES
        .iter()
        .copied()
        .zip(HOP_X_STATES.iter().copied())
        .collect_vec();
    static ref HOP_Y_HOP_Y: Vec<(State, State)> = HOP_Y_STATES
        .iter()
        .copied()
        .zip(HOP_Y_STATES.iter().copied())
        .collect_vec();
    static ref HOP_X_MATCH: Vec<(State, State)> = HOP_X_STATES
        .iter()
        .copied()
        .cartesian_product(MATCH_STATES.iter().copied())
        .collect_vec();
    static ref HOP_Y_MATCH: Vec<(State, State)> = HOP_Y_STATES
        .iter()
        .copied()
        .cartesian_product(MATCH_STATES.iter().copied())
        .collect_vec();
    static ref MATCH_SAME_: Vec<(State, State)> = MATCH_STATES
        .iter()
        .copied()
        .zip(MATCH_STATES.iter().copied())
        .collect_vec();
    static ref MATCH_OTHER: Vec<(State, State)> = MATCH_STATES
        .iter()
        .copied()
        .cartesian_product(MATCH_STATES.iter().copied())
        .filter(|(a, b)| a != b)
        .collect_vec();
}

fn build_transition_table<G: GapParameters, H: HopParameters>(
    gap_params: &G,
    hop_params: &H,
) -> Vec<TLogProb> {
    let n = STATES.len();
    let max_pair_label = (1 << (2 * n)) - 1;
    let mut transition_probs = vec![TLogProb::zero(); max_pair_label];

    let prob_hop_x = hop_params.prob_hop_x();
    let prob_hop_y = hop_params.prob_hop_y();
    let prob_hop_x_extend = hop_params.prob_hop_x_extend();
    let prob_hop_y_extend = hop_params.prob_hop_y_extend();

    let prob_gap_x = gap_params.prob_gap_x();
    let prob_gap_y = gap_params.prob_gap_y();
    let prob_gap_x_extend = gap_params.prob_gap_x_extend();
    let prob_gap_y_extend = gap_params.prob_gap_y_extend();

    MATCH_HOP_X
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = prob_hop_x);
    MATCH_HOP_Y
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = prob_hop_y);
    HOP_X_HOP_X
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = prob_hop_x_extend);
    HOP_Y_HOP_Y
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = prob_hop_y_extend);
    HOP_X_MATCH
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = (TLogProb::one() - prob_hop_x_extend) / 1.);
    HOP_Y_MATCH
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = (TLogProb::one() - prob_hop_y_extend) / 1.);

    let match_same = (TLogProb::one() - (prob_gap_y + prob_gap_x + prob_hop_x + prob_hop_y)) / 1.;
    let match_other = (TLogProb::one() - (prob_gap_y + prob_gap_x + prob_hop_x + prob_hop_y)) / 1.;
    MATCH_SAME_
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = match_same);
    MATCH_OTHER
        .iter()
        .for_each(|(a, b)| transition_probs[*a >> *b] = match_other);

    MATCH_STATES
        .iter()
        .for_each(|&a| transition_probs[a >> GapX] = prob_gap_y);
    MATCH_STATES
        .iter()
        .for_each(|&a| transition_probs[a >> GapY] = prob_gap_x);
    MATCH_STATES
        .iter()
        .for_each(|&b| transition_probs[GapX >> b] = (TLogProb::one() - prob_gap_y_extend) / 1.);
    MATCH_STATES
        .iter()
        .for_each(|&b| transition_probs[GapY >> b] = (TLogProb::one() - prob_gap_x_extend) / 1.);
    transition_probs[GapX >> GapX] = prob_gap_y_extend;
    transition_probs[GapY >> GapY] = prob_gap_x_extend;
    transition_probs
}

pub enum AlignmentMode {
    Global,
    Semiglobal,
}

impl StartEndGapParameters for AlignmentMode {
    fn free_start_gap_x(&self) -> bool {
        match self {
            AlignmentMode::Semiglobal => true,
            AlignmentMode::Global => false,
        }
    }

    fn free_end_gap_x(&self) -> bool {
        match self {
            AlignmentMode::Semiglobal => true,
            AlignmentMode::Global => false,
        }
    }
}

impl PairHHMM {
    pub fn new() -> Self {
        Self {}
    }

    /// Calculate the probability of sequence x being related to y via any alignment.
    ///
    /// # Arguments
    ///
    /// * `gap_params` - parameters for opening or extending gaps
    /// * `hop_params` - parameters for opening or extending hops
    /// * `emission_params` - parameters for emission
    /// * `max_edit_dist` - maximum edit distance to consider; if not `None`, perform banded alignment
    pub fn prob_related<G, H, E, A>(
        &mut self,
        gap_params: &G,
        hop_params: &H,
        emission_params: &E,
        alignment_mode: &A,
        max_edit_dist: Option<usize>,
    ) -> TLogProb
    where
        G: GapParameters,
        H: HopParameters,
        E: EmissionParameters,
        A: StartEndGapParameters,
    {
        let mut prev = 0;
        let mut curr = 1;
        let mut v: [EnumMap<State, Vec<TLogProb>>; 2] = [EnumMap::new(), EnumMap::new()];
        let transition_probs = build_transition_table(gap_params, hop_params);

        let len_y = emission_params.len_y();
        let len_x = emission_params.len_x();
        let mut min_edit_dist: [Vec<usize>; 2] =
            [vec![usize::MAX; len_y + 1], vec![usize::MAX; len_y + 1]];
        let free_end_gap_x = alignment_mode.free_end_gap_x();
        let free_start_gap_x = alignment_mode.free_start_gap_x();
        let mut prob_cols = Vec::with_capacity(len_x * STATES.len());

        for state in &STATES {
            v[prev][*state] = vec![TLogProb::zero(); len_y + 1];
        }

        v[curr] = v[prev].clone();

        for &m in &MATCH_STATES {
            v[prev][m][0] = TLogProb::one() / 4.;
        }

        for i in 0..len_x {
            if free_start_gap_x {
                for &m in &MATCH_STATES {
                    v[prev][m][0] += TLogProb::one() / 4.;
                }
                min_edit_dist[prev][0] = 0;
            }

            // cache probs for x[i]
            let prob_emit_x_and_gap = emission_params.prob_emit_x_and_gap(GapY, i);
            // let probs_emit_x_and_hop: SmallVec<[TLogProb; 4]> = HOP_Y_STATES
            //     .iter()
            //     .map(|&h| emission_params.prob_emit_x_and_hop(h, i))
            //     .collect();

            for j in 0..len_y {
                let j_ = j + 1;
                let j_minus_one = j_ - 1;

                let min_edit_dist_topleft = min_edit_dist[prev][j_minus_one];
                let min_edit_dist_top = min_edit_dist[curr][j_minus_one];
                let min_edit_dist_left = min_edit_dist[prev][j_];

                if let Some(max_edit_dist) = max_edit_dist {
                    if min3(min_edit_dist_topleft, min_edit_dist_top, min_edit_dist_left)
                        > max_edit_dist
                    {
                        // skip this cell if best edit dist is already larger than given maximum
                        continue;
                    }
                }

                let mut any_match = false;
                for &m in &MATCH_STATES {
                    let emission = emission_params.prob_emit_x_and_y(m, i, j);
                    any_match |= emission.is_match();
                    v[curr][m][j_] = emission.prob()
                        * STATES
                            .iter()
                            .map(|&s| transition_probs[s >> m] * v[prev][s][j_minus_one])
                            .sum::<TLogProb>();
                }

                v[curr][GapY][j_] = prob_emit_x_and_gap
                    * (MATCH_STATES
                        .iter()
                        .map(|&s| transition_probs[s >> GapY] * v[prev][s][j_])
                        .sum::<TLogProb>())
                    + transition_probs[GapY >> GapY] * v[prev][GapY][j_];

                MATCH_HOP_Y.iter().for_each(|&(m, h)| {
                    v[curr][h][j_] = emission_params.prob_emit_x_and_hop(h, i)
                        * (transition_probs[m >> h] * v[prev][m][j_]
                            + transition_probs[h >> h] * v[prev][h][j_])
                });

                v[curr][GapX][j_] = emission_params.prob_emit_gap_and_y(GapX, j)
                    * (MATCH_STATES
                        .iter()
                        .map(|&s| transition_probs[s >> GapX] * v[curr][s][j_minus_one])
                        .sum::<TLogProb>())
                    + transition_probs[GapX >> GapX] * v[curr][GapX][j_minus_one];

                MATCH_HOP_X.iter().for_each(|&(m, h)| {
                    v[curr][h][j_] = emission_params.prob_emit_hop_and_y(h, j)
                        * (transition_probs[m >> h] * v[curr][m][j_minus_one]
                            + transition_probs[h >> h] * v[curr][h][j_minus_one])
                });

                // calculate minimal number of mismatches
                if max_edit_dist.is_some() {
                    min_edit_dist[curr][j_] = min3(
                        if any_match {
                            // a match, so nothing changes
                            min_edit_dist_topleft
                        } else {
                            // one new mismatch
                            min_edit_dist_topleft.saturating_add(1)
                        },
                        // gap or hop in y (no new mismatch)
                        min_edit_dist_left.saturating_add(1),
                        // gap or hop in x (no new mismatch)
                        min_edit_dist_top.saturating_add(1),
                    )
                };

                if free_end_gap_x {
                    // Cache column probabilities or simply record the last probability.
                    // We can put all of them in one array since we simply have to sum in the end.
                    // This is also good for numerical stability.
                    prob_cols.extend(MATCH_STATES.iter().map(|&s| v[curr][s][len_y]));
                    prob_cols.extend(HOP_Y_STATES.iter().map(|&s| v[curr][s][len_y]));
                    prob_cols.extend(HOP_X_STATES.iter().map(|&s| v[curr][s][len_y]));
                    prob_cols.push(v[curr][GapY][len_y]);
                    // TODO check removing this (we don't want open gaps in x):
                    prob_cols.push(v[curr][GapX][len_y]);
                }
            }
            mem::swap(&mut prev, &mut curr);
            for &s in &MATCH_STATES {
                v[curr][s].reset(TLogProb::zero());
            }
        }
        if free_end_gap_x {
            prob_cols.iter().cloned().sum::<TLogProb>()
        } else {
            STATES
                .iter()
                .map(|&state| v[prev][state][len_y])
                .sum::<TLogProb>()
        }
    }
}

fn min3<T: Ord>(a: T, b: T, c: T) -> T {
    cmp::min(a, cmp::min(b, c))
}

#[cfg(test)]
mod tests {
    use crate::stats::pairhhmm::AlignmentMode::{Global, Semiglobal};
    use crate::stats::pairhmm::PairHMM;
    use crate::stats::probs::TLogProb;
    use crate::stats::{LogProb, Prob};

    use super::*;

    // Single base insertion and deletion rates for R1 according to Schirmer et al.
    // BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
    static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
    static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
    static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

    struct TestEmissionParams {
        x: &'static [u8],
        y: &'static [u8],
    }
    // log(0.0021)
    const PROB_SUBSTITUTION: TLogProb = TLogProb(-6.165_817_934_252_76);
    // log(2.8e-6)
    const PROB_OPEN_GAP_Y: TLogProb = TLogProb(-12.785_891_140_783_116);
    // log(5.1e-6)
    const PROB_OPEN_GAP_X: TLogProb = TLogProb(-12.186_270_018_233_994);

    const EMIT_MATCH: TLogProb = TLogProb(-0.0021022080918701985);
    const EMIT_GAP_AND_Y: TLogProb = TLogProb(-0.0021022080918701985);
    const EMIT_X_AND_GAP: TLogProb = TLogProb(-0.0021022080918701985);

    const T_MATCH_TO_HOP_X: TLogProb = TLogProb(-11.512925464970229);
    const T_MATCH_TO_HOP_Y: TLogProb = TLogProb(-11.512925464970229);
    const T_HOP_X_TO_HOP_X: TLogProb = TLogProb(-2.3025850929940455);
    const T_HOP_Y_TO_HOP_Y: TLogProb = TLogProb(-2.3025850929940455);

    const T_MATCH_TO_MATCH: TLogProb = TLogProb(-7.900031205113962e-06);
    const T_MATCH_TO_GAP_Y: TLogProb = TLogProb(-12.785_891_140_783_116);
    const T_MATCH_TO_GAP_X: TLogProb = TLogProb(-12.186_270_018_233_994);
    const T_GAP_TO_GAP: TLogProb = TLogProb(-9.210340371976182);

    impl EmissionParameters for TestEmissionParams {
        fn prob_emit_x_and_y(&self, s: State, i: usize, j: usize) -> XYEmission {
            let (a, b) = (self.x[i], self.y[j]);
            let p = match s {
                State::MatchA => match (a, b) {
                    (b'A', b'A') => TLogProb::one() - PROB_SUBSTITUTION,
                    (b'A', _y) => PROB_SUBSTITUTION / 6.,
                    (_x, b'A') => PROB_SUBSTITUTION / 6.,
                    _ => TLogProb::zero(),
                },
                State::MatchC => match (a, b) {
                    (b'C', b'C') => TLogProb::one() - PROB_SUBSTITUTION,
                    (b'C', _y) => PROB_SUBSTITUTION / 6.,
                    (_x, b'C') => PROB_SUBSTITUTION / 6.,
                    _ => TLogProb::zero(),
                },
                State::MatchG => match (a, b) {
                    (b'G', b'G') => TLogProb::one() - PROB_SUBSTITUTION,
                    (b'G', _y) => PROB_SUBSTITUTION / 6.,
                    (_x, b'G') => PROB_SUBSTITUTION / 6.,
                    _ => TLogProb::zero(),
                },
                State::MatchT => match (a, b) {
                    (b'T', b'T') => TLogProb::one() - PROB_SUBSTITUTION,
                    (b'T', _y) => PROB_SUBSTITUTION / 6.,
                    (_x, b'T') => PROB_SUBSTITUTION / 6.,
                    _ => TLogProb::zero(),
                },
                _ => TLogProb::zero(),
            };
            if self.x[i] == self.y[j] {
                XYEmission::Match(p)
            } else {
                XYEmission::Mismatch(p)
            }
        }

        fn prob_emit_x_and_gap(&self, _s: State, _i: usize) -> TLogProb {
            TLogProb::one() - PROB_SUBSTITUTION
        }

        fn prob_emit_gap_and_y(&self, _s: State, _j: usize) -> TLogProb {
            TLogProb::one() - PROB_SUBSTITUTION
        }

        fn prob_emit_x_and_hop(&self, s: State, i: usize) -> TLogProb {
            match s {
                State::HopAY => {
                    if self.x[i] == b'A' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                State::HopCY => {
                    if self.x[i] == b'C' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                State::HopGY => {
                    if self.x[i] == b'G' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                State::HopTY => {
                    if self.x[i] == b'T' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                _ => TLogProb::zero(),
            }
        }

        fn prob_emit_hop_and_y(&self, s: State, j: usize) -> TLogProb {
            match s {
                State::HopAX => {
                    if self.y[j] == b'A' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                State::HopCX => {
                    if self.y[j] == b'C' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                State::HopGX => {
                    if self.y[j] == b'G' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                State::HopTX => {
                    if self.y[j] == b'T' {
                        TLogProb::one() - PROB_SUBSTITUTION
                    } else {
                        PROB_SUBSTITUTION / 3.
                    }
                }
                _ => TLogProb::zero(),
            }
        }

        fn len_x(&self) -> usize {
            self.x.len()
        }

        fn len_y(&self) -> usize {
            self.y.len()
        }
    }

    struct TestSingleGapParams;

    impl GapParameters for TestSingleGapParams {
        fn prob_gap_x(&self) -> TLogProb {
            PROB_OPEN_GAP_Y
        }

        fn prob_gap_y(&self) -> TLogProb {
            PROB_OPEN_GAP_X
        }

        fn prob_gap_x_extend(&self) -> TLogProb {
            TLogProb::zero()
        }

        fn prob_gap_y_extend(&self) -> TLogProb {
            TLogProb::zero()
        }
    }

    struct NoGapParams;

    impl GapParameters for NoGapParams {
        fn prob_gap_x(&self) -> TLogProb {
            TLogProb::zero()
        }

        fn prob_gap_y(&self) -> TLogProb {
            TLogProb::zero()
        }

        fn prob_gap_x_extend(&self) -> TLogProb {
            TLogProb::zero()
        }

        fn prob_gap_y_extend(&self) -> TLogProb {
            TLogProb::zero()
        }
    }

    struct TestExtendGapParams;

    impl GapParameters for TestExtendGapParams {
        fn prob_gap_x(&self) -> TLogProb {
            TLogProb::from(PROB_ILLUMINA_INS)
        }

        fn prob_gap_y(&self) -> TLogProb {
            TLogProb::from(PROB_ILLUMINA_DEL)
        }

        fn prob_gap_x_extend(&self) -> TLogProb {
            T_GAP_TO_GAP
        }

        fn prob_gap_y_extend(&self) -> TLogProb {
            T_GAP_TO_GAP
        }
    }

    struct TestNoHopParams;
    impl HopParameters for TestNoHopParams {
        fn prob_hop_x(&self) -> TLogProb {
            TLogProb::zero()
        }

        fn prob_hop_y(&self) -> TLogProb {
            TLogProb::zero()
        }

        fn prob_hop_x_extend(&self) -> TLogProb {
            TLogProb::zero()
        }

        fn prob_hop_y_extend(&self) -> TLogProb {
            TLogProb::zero()
        }
    }

    struct TestHopParams;

    impl HopParameters for TestHopParams {
        fn prob_hop_x(&self) -> TLogProb {
            T_MATCH_TO_HOP_X
        }

        fn prob_hop_y(&self) -> TLogProb {
            T_MATCH_TO_HOP_Y
        }

        fn prob_hop_x_extend(&self) -> TLogProb {
            T_HOP_X_TO_HOP_X
        }

        fn prob_hop_y_extend(&self) -> TLogProb {
            T_HOP_Y_TO_HOP_Y
        }
    }

    #[test]
    fn impossible_global_alignment() {
        let x = b"AAA";
        let y = b"A";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );
        assert_eq!(p, TLogProb::zero());
    }

    #[test]
    fn test_hompolymer_run_in_y() {
        let x = b"ACGT";
        let y = b"ACCCCGT";
        let emission_params = TestEmissionParams { x, y };
        let hop_params = TestHopParams;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&NO_GAP_PARAMS, &hop_params, &emission_params, &Global, None);
        let p_most_likely_path_with_hops = LogProb(
            *EMIT_MATCH // A A
                + *T_MATCH_TO_MATCH
                + *EMIT_MATCH // C C
                + *T_MATCH_TO_HOP_X
                + *EMIT_MATCH // C CC
                + *T_HOP_X_TO_HOP_X
                + *EMIT_MATCH // C CCC
                + *T_HOP_X_TO_HOP_X
                + *EMIT_MATCH // C CCCC
                + (1. - 0.1f64).ln()
                + *EMIT_MATCH // G G
                + *T_MATCH_TO_MATCH
                + *EMIT_MATCH, // T T
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path_with_hops, *p, epsilon = 0.001);
    }

    #[test]
    fn test_hompolymer_run_in_x() {
        let x = b"ACCCCGT";
        let y = b"ACGT";
        let emission_params = TestEmissionParams { x, y };
        let hop_params = TestHopParams;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&NO_GAP_PARAMS, &hop_params, &emission_params, &Global, None);
        let p_most_likely_path_with_hops = LogProb(
            *EMIT_MATCH // A A
                + *T_MATCH_TO_MATCH
                + *EMIT_MATCH // C C
                + *T_MATCH_TO_HOP_Y
                + *EMIT_MATCH // CC C
                + *T_HOP_Y_TO_HOP_Y
                + *EMIT_MATCH // CCC C
                + *T_HOP_Y_TO_HOP_Y
                + *EMIT_MATCH // CCCC C
                + (1. - 0.1f64).ln()
                + *EMIT_MATCH // G G
                + *T_MATCH_TO_MATCH
                + *EMIT_MATCH, // T T
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path_with_hops, *p, epsilon = 0.001);
    }

    #[test]
    fn test_interleave_gaps_x() {
        let x = b"AGAGAG";
        let y = b"ACGTACGTACGT";

        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );

        let n_matches = 6.;
        let n_insertions = 6.;

        let p_most_likely_path = TLogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_AND_Y * n_insertions
                + *T_MATCH_TO_GAP_X * n_insertions
                + *(TLogProb::one() - PROB_OPEN_GAP_Y) * n_insertions,
        );

        let p_max = TLogProb(*T_MATCH_TO_GAP_X * n_insertions);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_interleave_gaps_y() {
        let x = b"ACGTACGTACGT";
        let y = b"AGAGAG";

        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );

        let n_matches = 6.;
        let n_insertions = 6.;

        let p_most_likely_path = TLogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_X_AND_GAP * n_insertions
                + *T_MATCH_TO_GAP_Y * n_insertions
                + *(TLogProb::one() - PROB_OPEN_GAP_X) * n_insertions,
        );

        let p_max = TLogProb(*T_MATCH_TO_GAP_Y * n_insertions);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }
    static SINGLE_GAP_PARAMS: TestSingleGapParams = TestSingleGapParams;
    static EXTEND_GAP_PARAMS: TestExtendGapParams = TestExtendGapParams;
    static NO_GAP_PARAMS: NoGapParams = NoGapParams;
    static NO_HOP_PARAMS: TestNoHopParams = TestNoHopParams;

    #[test]
    fn test_same() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCGATCGATC";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );
        let n = x.len() as f64;
        let p_most_likely_path = LogProb(*EMIT_MATCH * n + *T_MATCH_TO_MATCH * (n - 1.));
        let p_max = LogProb(*EMIT_MATCH * n);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.001);
        assert_relative_eq!(*p, *p_max, epsilon = 0.001);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_x() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCTGATCGATCT";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );

        let n_matches = 17.;
        let n_insertions = 2.;

        let p_most_likely_path = TLogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_AND_Y * n_insertions
                + *T_MATCH_TO_GAP_X * n_insertions
                + (1. - *PROB_ILLUMINA_INS).ln(),
        );

        let p_max = TLogProb(*T_MATCH_TO_GAP_X * 2.);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_x_2() {
        let x = b"ACAGTA";
        let y = b"ACAGTCA";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );

        let n_matches = 6.;
        let n_insertions = 1.;

        let p_most_likely_path = TLogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_AND_Y * n_insertions
                + *T_MATCH_TO_GAP_X * n_insertions
                + (1. - *PROB_ILLUMINA_INS).ln(),
        );

        let p_max = TLogProb(*T_MATCH_TO_GAP_X * n_insertions);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_y() {
        let x = b"AGCTCGATCTGATCGATCT";
        let y = b"AGCTCGATCGATCGATC";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );

        let n_matches = 17.;
        let n_deletions = 2.;

        let p_most_likely_path = TLogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_deletions)
                + *EMIT_X_AND_GAP * n_deletions
                + *T_MATCH_TO_GAP_Y * n_deletions
                + (1. - *PROB_ILLUMINA_DEL).ln(),
        );

        let p_max = TLogProb(*T_MATCH_TO_GAP_Y * 2.);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_multigap_y() {
        let x = b"AGCTCGATCTGATCGATCT";
        let y = b"AGCTTCTGATCGATCT";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &EXTEND_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );
        let n_matches = 16.;
        let n_consecutive_deletions = 3.;
        let p_most_likely_path = TLogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_consecutive_deletions)
                + *PROB_OPEN_GAP_Y
                + *EMIT_X_AND_GAP * n_consecutive_deletions
                + *T_GAP_TO_GAP * (n_consecutive_deletions - 1.)
                + *(TLogProb::one() - T_GAP_TO_GAP),
        );

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
    }

    #[test]
    fn test_mismatch() {
        let x = b"AGCTCGAGCGATCGATC";
        let y = b"TGCTCGATCGATCGATC";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );

        let n = x.len() as f64;
        let p_most_likely_path = TLogProb(
            *EMIT_MATCH * (n - 2.)
                + *T_MATCH_TO_MATCH * (n - 1.)
                + (*PROB_ILLUMINA_SUBST / 3.).ln() * 2.,
        );
        let p_max = TLogProb((*PROB_ILLUMINA_SUBST / 3.).ln() * 2.);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 1e-2);
        assert_relative_eq!(*p, *p_max, epsilon = 1e-1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_banded() {
        let x = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC\
ATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTAT\
CTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGT";

        let y = b"GGGTATGCACGCGATAGCATTGCGAGATGCTGGAGCTGGAGCACCCTATGTCGC";

        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Semiglobal,
            None,
        );

        let p_banded = pair_hmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Semiglobal,
            Some(2),
        );
        assert_relative_eq!(*p, *p_banded, epsilon = 1e-3);
    }

    #[test]
    fn test_phmm_vs_phhmm() {
        let x = b"AGAGC";
        let y = b"ACGTACGTC";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hhmm = PairHHMM::new();
        let p1 = pair_hhmm.prob_related(
            &SINGLE_GAP_PARAMS,
            &NO_HOP_PARAMS,
            &emission_params,
            &Global,
            None,
        );

        let mut pair_hmm = PairHMM::new();

        struct TestSingleGapParamsPairHMM;
        impl crate::stats::pairhmm::StartEndGapParameters for TestSingleGapParamsPairHMM {
            fn free_start_gap_x(&self) -> bool {
                false
            }

            fn free_end_gap_x(&self) -> bool {
                false
            }
        }
        impl crate::stats::pairhmm::GapParameters for TestSingleGapParamsPairHMM {
            fn prob_gap_x(&self) -> LogProb {
                LogProb::from(PROB_ILLUMINA_DEL)
            }

            fn prob_gap_y(&self) -> LogProb {
                LogProb::from(PROB_ILLUMINA_INS)
            }

            fn prob_gap_x_extend(&self) -> LogProb {
                LogProb::zero()
            }

            fn prob_gap_y_extend(&self) -> LogProb {
                LogProb::zero()
            }
        }

        fn prob_emit_x_or_y() -> LogProb {
            LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
        }

        struct TestEmissionParamsPairHMM {
            x: &'static [u8],
            y: &'static [u8],
        }

        impl crate::stats::pairhmm::EmissionParameters for TestEmissionParamsPairHMM {
            fn prob_emit_xy(&self, i: usize, j: usize) -> crate::stats::pairhmm::XYEmission {
                if self.x[i] == self.y[j] {
                    crate::stats::pairhmm::XYEmission::Match(LogProb::from(
                        Prob(1.0) - PROB_ILLUMINA_SUBST,
                    ))
                } else {
                    crate::stats::pairhmm::XYEmission::Mismatch(LogProb::from(
                        PROB_ILLUMINA_SUBST / Prob(3.0),
                    ))
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

        let p2 = pair_hmm.prob_related(
            &TestSingleGapParamsPairHMM,
            &TestEmissionParamsPairHMM { x, y },
            None,
        );
        assert_relative_eq!(*p1, *p2, epsilon = 1e-4)
    }
}
