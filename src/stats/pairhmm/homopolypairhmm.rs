// Copyright 2014-2016 Johannes Köster.
// Copyright 2020 Till Hartmann.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A pair Hidden Markov Model for calculating the probability that two sequences are related to
//! each other. Depending on the used parameters, this can, e.g., be used to calculate the
//! probability that a certain sequencing read comes from a given position in a reference genome.
//! In contrast to `PairHMM`, this `HomopolyPairHMM` takes into account homopolymer errors as
//! often encountered e.g. in Oxford Nanopore Technologies sequencing.
//!
//! Time complexity: O(n * m) where `n = seq1.len()`, `m = seq2.len()` (or `m = min(seq2.len(), max_edit_dist)` with banding enabled).
//! Memory complexity: O(m) where `m = seq2.len()`.
//! Note that if the number of states weren't fixed in this implementation, we would have to include
//! these in both time and memory complexity above as an additional factor.
//!
//! The `HomopolyPairHMM` introduces the term "hop" for starting and extending homopolymer runs
//! by analogy with "gap". Therefore, the constructor needs an additional parameter `hop_params`
//! implementing `HopParameters`. Also, the emission parameter needs to implement `Emission`,
//! since this HMM model needs to be able to distinguish the four different match states for
//! A, C, G and T (see Details below).
//!
//! # Details
//! The HomopolyPairHMM defined in this module has one Match state for each character from [A, C, G, T],
//! for each of those Match states two corresponding Hop (homopolymer run) states
//! (one for a run in sequence `x`, one for a run in `y`),
//! as well as the usual GapX and GapY states.
//!
//! In states `MatchV` (where `V` ∈ `{A, C, G, T}`), the probability to emit anything other than
//! `(V, V)`, `(V, y != V)`, `(x != V, y)` should be zero.
//!
//! State `HopVZ` (where `V` ∈ `{A, C, G, T}`, `Z` ∈ `{X, Y}`) can only be transitioned to from
//! corresponding state `MatchV`.
//!
//! The transition matrix is given below:
//!     | MA | MC | MG | MT | HAX | HAY | HCX | HCY | HGX | HGY | HTX | HTY | GX | GY
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! MA  |  x |  x |  x |  x |  x  |  x  |     |     |     |     |     |     |  x |  x
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! MC  |  x |  x |  x |  x |     |     |  x  |  x  |     |     |     |     |  x |  x
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! MG  |  x |  x |  x |  x |     |     |     |     |  x  |  x  |     |     |  x |  x
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! MT  |  x |  x |  x |  x |     |     |     |     |     |     |  x  |  x  |  x |  x
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HAX |  x |  x |  x |  x |  x  |     |     |     |     |     |     |     |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HAY |  x |  x |  x |  x |     |  x  |     |     |     |     |     |     |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HCX |  x |  x |  x |  x |     |     |  x  |     |     |     |     |     |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HCY |  x |  x |  x |  x |     |     |     |  x  |     |     |     |     |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HGX |  x |  x |  x |  x |     |     |     |     |  x  |     |     |     |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HGY |  x |  x |  x |  x |     |     |     |     |     |  x  |     |     |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HTX |  x |  x |  x |  x |     |     |     |     |     |     |  x  |     |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! HTY |  x |  x |  x |  x |     |     |     |     |     |     |     |  x  |    |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! GX  |  x |  x |  x |  x |     |     |     |     |     |     |     |     |  x |
//! ----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|----|---
//! GY  |  x |  x |  x |  x |     |     |     |     |     |     |     |     |    |  x

use std::cmp;
use std::fmt::Debug;
use std::mem;
use std::ops::Shr;
use usize;

use enum_map::{Enum, EnumMap};
use itertools::Itertools;
use num_traits::Zero;

use crate::alphabets::dna::iupac_mask;
use crate::stats::pairhmm::homopolypairhmm::State::*;
use crate::stats::pairhmm::{Emission, EmissionParameters, GapParameters, StartEndGapParameters};
use crate::stats::probs::LogProb;
use crate::stats::Prob;

#[repr(usize)]
#[derive(
    Enum, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
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

impl State {
    #[cfg(test)]
    fn supports(&self, x: u8, y: u8) -> bool {
        match self {
            // For each match state, check if the base of the match state is supported by the IUPAC mask of x and y.
            MatchA | MatchC | MatchG | MatchT => {
                iupac_mask(self.base().unwrap()) & (iupac_mask(x) | iupac_mask(y)) != 0
            }
            _ => false,
        }
    }

    fn base(&self) -> Option<u8> {
        match self {
            MatchA | HopAX | HopAY => Some(b'A'),
            MatchC | HopCX | HopCY => Some(b'C'),
            MatchG | HopGX | HopGY => Some(b'G'),
            MatchT | HopTX | HopTY => Some(b'T'),
            _ => None,
        }
    }
}

const NUM_STATES: usize = 14;

const STATES: [State; NUM_STATES] = [
    MatchA, MatchC, MatchG, MatchT, GapX, GapY, HopAX, HopAY, HopCX, HopCY, HopGX, HopGY, HopTX,
    HopTY,
];

const MATCH_STATES: [State; 4] = [MatchA, MatchC, MatchG, MatchT];
const HOP_X_STATES: [State; 4] = [HopAX, HopCX, HopGX, HopTX];
const HOP_Y_STATES: [State; 4] = [HopAY, HopCY, HopGY, HopTY];

// We define Shr (>>) for `State` such that a transition from State `a` to State `b` can be modeled
// as `a >> b`, where `a >> b` is an integer in `0..(1 << (2 * NUM_STATES)) - 1]` used for indexing
// the transition table (see `build_transition_table`).
impl Shr for State {
    type Output = usize;

    fn shr(self, rhs: State) -> Self::Output {
        let a = self as u32;
        let b = rhs as u32;
        interleave_bits(a, b) as usize
    }
}

fn space_bits(a: u32) -> u64 {
    let mut x = a as u64 & 0x0000_0000_FFFF_FFFF;
    x = (x | (x << 16)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x << 8)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x << 2)) & 0x3333_3333_3333_3333;
    x = (x | (x << 1)) & 0x5555_5555_5555_5555;
    x
}

fn interleave_bits(a: u32, b: u32) -> u64 {
    space_bits(a) << 1 | space_bits(b)
}

/// Trait for parametrization of `PairHMM` hop behavior.
pub trait HopParameters {
    /// Probability to start hop in x.
    fn prob_hop_x(&self) -> LogProb;

    /// Probability to start hop in y.
    fn prob_hop_y(&self) -> LogProb;

    /// Probability to extend hop in x.
    fn prob_hop_x_extend(&self) -> LogProb;

    /// Probability to extend hop in y.
    fn prob_hop_y_extend(&self) -> LogProb;
}

/// Trait for parametrization of `PairHMM` hop behavior.
pub trait BaseSpecificHopParameters {
    /// Probability to start hop in x.
    fn prob_hop_x_with_base(&self, base: u8) -> LogProb;

    /// Probability to start hop in y.
    fn prob_hop_y_with_base(&self, base: u8) -> LogProb;

    /// Probability to extend hop in x.
    fn prob_hop_x_extend_with_base(&self, base: u8) -> LogProb;

    /// Probability to extend hop in y.
    fn prob_hop_y_extend_with_base(&self, base: u8) -> LogProb;
}

impl<H: HopParameters> BaseSpecificHopParameters for H {
    fn prob_hop_x_with_base(&self, _base: u8) -> LogProb {
        self.prob_hop_x()
    }

    fn prob_hop_y_with_base(&self, _base: u8) -> LogProb {
        self.prob_hop_y()
    }

    fn prob_hop_x_extend_with_base(&self, _base: u8) -> LogProb {
        self.prob_hop_x_extend()
    }

    fn prob_hop_y_extend_with_base(&self, _base: u8) -> LogProb {
        self.prob_hop_y_extend()
    }
}

/// A pair Hidden Markov Model for comparing sequences x and y as described by
/// Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
/// Current Topics in Genome Analysis 2008. http://doi.org/10.1017/CBO9780511790492.
/// The default model has been extended to consider homopolymer errors, at the cost of more states
/// and transitions.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct HomopolyPairHMM {
    /// Transition probabilities, indexed by `[from as usize][to as usize]`. Transitions that
    /// the model does not have are zero.
    transition_probs: [[LogProb; NUM_STATES]; NUM_STATES],
}

impl Default for HomopolyPairHMM {
    fn default() -> Self {
        HomopolyPairHMM {
            transition_probs: [[LogProb::ln_zero(); NUM_STATES]; NUM_STATES],
        }
    }
}

impl HomopolyPairHMM {
    /// Create a new instance of a HomopolyPairHMM.
    /// # Arguments
    ///
    /// * `gap_params` - parameters for opening or extending gaps
    /// * `hop_params` - parameters for opening or extending hops
    pub fn new<G, H>(gap_params: &G, hop_params: &H) -> Self
    where
        G: GapParameters,
        H: BaseSpecificHopParameters,
    {
        Self {
            transition_probs: build_transition_table(gap_params, hop_params),
        }
    }

    /// Calculate the probability of sequence x being related to y via any alignment.
    ///
    /// # Arguments
    ///
    /// * `emission_params` - parameters for emission
    /// * `alignment_mode` - parameters for free end/start gaps
    /// * `max_edit_dist` - maximum edit distance to consider; if not `None`, perform banded alignment
    pub fn prob_related<E, A>(
        &self,
        emission_params: &E,
        alignment_mode: &A,
        max_edit_dist: Option<usize>,
    ) -> LogProb
    where
        E: EmissionParameters + Emission,
        A: StartEndGapParameters,
    {
        let mut prev = 0;
        let mut curr = 1;
        let mut v: [EnumMap<State, Vec<LogProb>>; 2] = [EnumMap::default(), EnumMap::default()];
        let transition_probs = &self.transition_probs;
        // IUPAC masks of the bases A, C, G and T, in the order of MATCH_STATES (and of the
        // match/hop pairs below).
        let base_masks = MATCH_STATES.map(|m| iupac_mask(m.base().unwrap()));

        let len_y = emission_params.len_y();
        let len_x = emission_params.len_x();
        let mut min_edit_dist: [Vec<usize>; 2] =
            [vec![usize::MAX; len_y + 1], vec![usize::MAX; len_y + 1]];
        min_edit_dist[0][0] = 0;
        let free_end_gap_x = alignment_mode.free_end_gap_x();
        let free_start_gap_x = alignment_mode.free_start_gap_x();

        let mut prob_cols = Vec::with_capacity(len_x * STATES.len());

        for state in &STATES {
            v[prev][*state] = vec![LogProb::zero(); len_y + 1];
        }

        v[curr] = v[prev].clone();

        for &m in &MATCH_STATES {
            v[prev][m][0] = LogProb::from(Prob(1. / 4.));
        }

        for i in 0..len_x {
            if free_start_gap_x {
                let prob_start_gap_x = LogProb(*alignment_mode.prob_start_gap_x(i) - 4f64.ln());
                for &m in &MATCH_STATES {
                    v[prev][m][0] = v[prev][m][0].ln_add_exp(prob_start_gap_x);
                }
                min_edit_dist[prev][0] = 0;
            }

            // cache probs for x[i]
            let prob_emit_x_and_gap = emission_params.prob_emit_x(i);
            let emission_x = emission_params.emission_x(i);
            let mask_x = iupac_mask(emission_x);
            // A hop in y can only continue a run of a base that x supports.
            let hop_y_supported = base_masks.map(|mask| mask & mask_x != 0);

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

                let emission_y = emission_params.emission_y(j);
                let mask_y = iupac_mask(emission_y);
                let mut any_match = false;

                // The match states supporting the pair (x[i], y[j]).
                let mask_xy = mask_x | mask_y;
                let supported = base_masks.map(|mask| mask & mask_xy != 0);
                let num_states = supported.iter().filter(|&&supported| supported).count() as f64;
                let ln_num_states = num_states.ln();

                for (k, &m) in MATCH_STATES.iter().enumerate() {
                    if supported[k] {
                        // The emission is queried per match state, since it can differ between the active states (e.g. for x=T, y=Y, MatchT is a match while MatchC is a mismatch).
                        // Emission prob is the log prob of the emission, adjusted by the number of states. If we have an exact match like ('A', 'A') there is only one active state so we do not need to adjust the emission prob. If there is a mismatch like ('A', 'C') we need to halve the emission prob to account that MatchA and MatchC are active. If there is a IUPAC ambiguity like ('A', 'Y') we need to account for the possible match states MatchA, MatchC, and MatchT.

                        let emission =
                            emission_params.prob_emit_xy_for_base(i, j, m.base().unwrap());
                        let emission_prob = LogProb::from(*emission.prob() - ln_num_states);

                        any_match |= emission.is_match();
                        let mut terms = [LogProb::ln_zero(); NUM_STATES];
                        for (term, &s) in terms.iter_mut().zip(STATES.iter()) {
                            *term =
                                transition_probs[s as usize][m as usize] + v[prev][s][j_minus_one];
                        }
                        v[curr][m][j_] = emission_prob + LogProb::ln_sum_exp(&terms);
                    } else {
                        v[curr][m][j_] = LogProb::zero();
                    }
                }

                let mut gap_y_terms = [LogProb::ln_zero(); 5];
                for (term, &s) in gap_y_terms.iter_mut().zip(MATCH_STATES.iter()) {
                    *term = transition_probs[s as usize][GapY as usize] + v[prev][s][j_];
                }
                gap_y_terms[4] = transition_probs[GapY as usize][GapY as usize] + v[prev][GapY][j_];
                v[curr][GapY][j_] = prob_emit_x_and_gap + LogProb::ln_sum_exp(&gap_y_terms);

                for (k, &(m, h)) in MATCH_HOP_Y.iter().enumerate() {
                    v[curr][h][j_] = if hop_y_supported[k] {
                        (transition_probs[m as usize][h as usize] + v[prev][m][j_])
                            .ln_add_exp(transition_probs[h as usize][h as usize] + v[prev][h][j_])
                    } else {
                        LogProb::zero()
                    }
                }

                let mut gap_x_terms = [LogProb::ln_zero(); 5];
                for (term, &s) in gap_x_terms.iter_mut().zip(MATCH_STATES.iter()) {
                    *term = transition_probs[s as usize][GapX as usize] + v[curr][s][j_minus_one];
                }
                gap_x_terms[4] =
                    transition_probs[GapX as usize][GapX as usize] + v[curr][GapX][j_minus_one];
                v[curr][GapX][j_] =
                    emission_params.prob_emit_y(j) + LogProb::ln_sum_exp(&gap_x_terms);

                for (k, &(m, h)) in MATCH_HOP_X.iter().enumerate() {
                    v[curr][h][j_] = if base_masks[k] & mask_y != 0 {
                        (transition_probs[m as usize][h as usize] + v[curr][m][j_minus_one])
                            .ln_add_exp(
                                transition_probs[h as usize][h as usize] + v[curr][h][j_minus_one],
                            )
                    } else {
                        LogProb::zero()
                    };
                }

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
            }
            if free_end_gap_x {
                // Record the probability of ending in this column, once per column and once the
                // last row has been computed: the hop and gap rows are not reset between
                // columns, so reading them earlier would pick up the previous column's values.
                // All of them go in one array since we simply have to sum in the end, which is
                // also good for numerical stability.
                prob_cols.extend(MATCH_STATES.iter().map(|&s| v[curr][s][len_y]));
                prob_cols.extend(HOP_Y_STATES.iter().map(|&s| v[curr][s][len_y]));
                prob_cols.extend(HOP_X_STATES.iter().map(|&s| v[curr][s][len_y]));
                prob_cols.push(v[curr][GapY][len_y]);
                // TODO check removing this (we don't want open gaps in x):
                prob_cols.push(v[curr][GapX][len_y]);
            }
            mem::swap(&mut prev, &mut curr);
            for &s in &MATCH_STATES {
                v[curr][s].reset(LogProb::zero());
            }
            // In global mode only the origin is allowed to be 0. After the swap we have to overwrite the min_edit_dist in order to not allow free start gaps.
            if max_edit_dist.is_some() {
                min_edit_dist[curr].reset(usize::MAX);
            }
        }
        let p = if free_end_gap_x {
            LogProb::ln_sum_exp(&prob_cols.iter().cloned().collect_vec())
        } else {
            LogProb::ln_sum_exp(
                &STATES
                    .iter()
                    .map(|&state| v[prev][state][len_y])
                    .collect_vec(),
            )
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

// explicitly defined groups of transitions between states
const MATCH_HOP_X: [(State, State); 4] = [
    (MatchA, HopAX),
    (MatchC, HopCX),
    (MatchG, HopGX),
    (MatchT, HopTX),
];
const MATCH_HOP_Y: [(State, State); 4] = [
    (MatchA, HopAY),
    (MatchC, HopCY),
    (MatchG, HopGY),
    (MatchT, HopTY),
];
const HOP_X_HOP_X: [(State, State); 4] = [
    (HopAX, HopAX),
    (HopCX, HopCX),
    (HopGX, HopGX),
    (HopTX, HopTX),
];
const HOP_Y_HOP_Y: [(State, State); 4] = [
    (HopAY, HopAY),
    (HopCY, HopCY),
    (HopGY, HopGY),
    (HopTY, HopTY),
];
const HOP_X_MATCH: [(State, State); 16] = [
    (HopAX, MatchA),
    (HopAX, MatchC),
    (HopAX, MatchG),
    (HopAX, MatchT),
    (HopCX, MatchA),
    (HopCX, MatchC),
    (HopCX, MatchG),
    (HopCX, MatchT),
    (HopGX, MatchA),
    (HopGX, MatchC),
    (HopGX, MatchG),
    (HopGX, MatchT),
    (HopTX, MatchA),
    (HopTX, MatchC),
    (HopTX, MatchG),
    (HopTX, MatchT),
];
const HOP_Y_MATCH: [(State, State); 16] = [
    (HopAY, MatchA),
    (HopAY, MatchC),
    (HopAY, MatchG),
    (HopAY, MatchT),
    (HopCY, MatchA),
    (HopCY, MatchC),
    (HopCY, MatchG),
    (HopCY, MatchT),
    (HopGY, MatchA),
    (HopGY, MatchC),
    (HopGY, MatchG),
    (HopGY, MatchT),
    (HopTY, MatchA),
    (HopTY, MatchC),
    (HopTY, MatchG),
    (HopTY, MatchT),
];
const MATCH_SAME_: [(State, State); 4] = [
    (MatchA, MatchA),
    (MatchC, MatchC),
    (MatchG, MatchG),
    (MatchT, MatchT),
];
const MATCH_OTHER: [(State, State); 12] = [
    (MatchA, MatchC),
    (MatchA, MatchG),
    (MatchA, MatchT),
    (MatchC, MatchA),
    (MatchC, MatchG),
    (MatchC, MatchT),
    (MatchG, MatchC),
    (MatchG, MatchA),
    (MatchG, MatchT),
    (MatchT, MatchC),
    (MatchT, MatchG),
    (MatchT, MatchA),
];

fn build_transition_table<G: GapParameters, H: BaseSpecificHopParameters>(
    gap_params: &G,
    hop_params: &H,
) -> [[LogProb; NUM_STATES]; NUM_STATES] {
    let mut transition_probs = [[LogProb::ln_zero(); NUM_STATES]; NUM_STATES];
    let mut set = |from: State, to: State, p: LogProb| {
        transition_probs[from as usize][to as usize] = p;
    };

    let prob_gap_x = gap_params.prob_gap_x();
    let prob_gap_y = gap_params.prob_gap_y();
    let prob_gap_x_extend = gap_params.prob_gap_x_extend();
    let prob_gap_y_extend = gap_params.prob_gap_y_extend();

    MATCH_HOP_X.iter().for_each(|(a, b)| {
        set(
            *a,
            *b,
            hop_params.prob_hop_x_with_base(b.base().expect("Unsupported base")),
        );
    });
    MATCH_HOP_Y.iter().for_each(|(a, b)| {
        set(
            *a,
            *b,
            hop_params.prob_hop_y_with_base(b.base().expect("Unsupported base")),
        );
    });
    HOP_X_HOP_X.iter().for_each(|(a, b)| {
        assert_eq!(a.base(), b.base());
        set(
            *a,
            *b,
            hop_params.prob_hop_x_extend_with_base(b.base().expect("Unsupported base")),
        );
    });
    HOP_Y_HOP_Y.iter().for_each(|(a, b)| {
        assert_eq!(a.base(), b.base());
        set(
            *a,
            *b,
            hop_params.prob_hop_y_extend_with_base(b.base().expect("Unsupported base")),
        );
    });
    HOP_X_MATCH.iter().for_each(|(a, b)| {
        set(
            *a,
            *b,
            hop_params
                .prob_hop_x_with_base(a.base().expect("Unsupported base"))
                .ln_one_minus_exp(),
        );
    });
    HOP_Y_MATCH.iter().for_each(|(a, b)| {
        set(
            *a,
            *b,
            hop_params
                .prob_hop_y_with_base(a.base().expect("Unsupported base"))
                .ln_one_minus_exp(),
        );
    });

    let prob_hop_x = LogProb::ln_sum_exp(&[
        hop_params.prob_hop_x_with_base(b'A'),
        hop_params.prob_hop_x_with_base(b'C'),
        hop_params.prob_hop_x_with_base(b'G'),
        hop_params.prob_hop_x_with_base(b'T'),
    ]) - LogProb(4.0);
    let prob_hop_y = LogProb::ln_sum_exp(&[
        hop_params.prob_hop_y_with_base(b'A'),
        hop_params.prob_hop_y_with_base(b'C'),
        hop_params.prob_hop_y_with_base(b'G'),
        hop_params.prob_hop_y_with_base(b'T'),
    ]) - LogProb(4.0);
    let match_same =
        LogProb::ln_sum_exp(&[prob_gap_y, prob_gap_x, prob_hop_x, prob_hop_y]).ln_one_minus_exp();
    let match_other =
        LogProb::ln_sum_exp(&[prob_gap_y, prob_gap_x, prob_hop_x, prob_hop_y]).ln_one_minus_exp();
    MATCH_SAME_.iter().for_each(|(a, b)| {
        set(*a, *b, match_same);
    });
    MATCH_OTHER.iter().for_each(|(a, b)| {
        set(*a, *b, match_other);
    });

    // GapX consumes a base of y only (a gap in x), GapY a base of x only (a gap in y)
    MATCH_STATES.iter().for_each(|&a| {
        set(a, GapX, prob_gap_x);
    });
    MATCH_STATES.iter().for_each(|&a| {
        set(a, GapY, prob_gap_y);
    });
    MATCH_STATES.iter().for_each(|&b| {
        set(GapX, b, prob_gap_x_extend.ln_one_minus_exp());
    });
    MATCH_STATES.iter().for_each(|&b| {
        set(GapY, b, prob_gap_y_extend.ln_one_minus_exp());
    });
    set(GapX, GapX, prob_gap_x_extend);
    set(GapY, GapY, prob_gap_y_extend);
    transition_probs
}

trait Reset<T: Copy> {
    fn reset(&mut self, value: T);
}

impl<T: Copy> Reset<T> for [T] {
    fn reset(&mut self, value: T) {
        for v in self {
            *v = value;
        }
    }
}

fn min3<T: Ord>(a: T, b: T, c: T) -> T {
    cmp::min(a, cmp::min(b, c))
}

#[cfg(test)]
mod tests {
    use crate::stats::pairhmm::homopolypairhmm::tests::AlignmentMode::{Global, Semiglobal};
    use crate::stats::pairhmm::{EmissionParameters, PairHMM, XYEmission};
    use crate::stats::{LogProb, Prob};
    use std::iter::repeat;
    use std::sync::LazyLock;

    use super::*;

    // Single base insertion and deletion rates for R1 according to Schirmer et al.
    // BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
    static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
    static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
    static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

    // log(0.0021)
    const PROB_SUBSTITUTION: LogProb = LogProb(-6.165_817_934_252_76);
    // log(2.8e-6): a gap in x is an insertion in y
    const PROB_OPEN_GAP_X: LogProb = LogProb(-12.785_891_140_783_116);
    // log(5.1e-6): a gap in y is a deletion in y
    const PROB_OPEN_GAP_Y: LogProb = LogProb(-12.186_270_018_233_994);

    const EMIT_MATCH: LogProb = LogProb(-0.0021022080918701985);
    const EMIT_GAP_AND_Y: LogProb = LogProb(-0.0021022080918701985);
    const EMIT_X_AND_GAP: LogProb = LogProb(-0.0021022080918701985);

    const T_MATCH_TO_HOP_X: LogProb = LogProb(-11.512925464970229);
    const T_MATCH_TO_HOP_Y: LogProb = LogProb(-11.512925464970229);
    const T_HOP_X_TO_HOP_X: LogProb = LogProb(-2.3025850929940455);
    const T_HOP_Y_TO_HOP_Y: LogProb = LogProb(-2.3025850929940455);

    const T_MATCH_TO_MATCH: LogProb = LogProb(-7.900_031_205_113_962e-6);
    const T_MATCH_TO_GAP_X: LogProb = LogProb(-12.785_891_140_783_116);
    const T_MATCH_TO_GAP_Y: LogProb = LogProb(-12.186_270_018_233_994);
    const T_GAP_TO_GAP: LogProb = LogProb(-9.210340371976182);

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

    struct TestEmissionParams {
        x: Vec<u8>,
        y: Vec<u8>,
    }

    impl EmissionParameters for TestEmissionParams {
        fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
            if self.x[i] == self.y[j] {
                XYEmission::Match(PROB_SUBSTITUTION.ln_one_minus_exp())
            } else {
                XYEmission::Mismatch(LogProb::from(PROB_ILLUMINA_SUBST / Prob(3.)))
            }
        }

        fn prob_emit_x(&self, _i: usize) -> LogProb {
            PROB_SUBSTITUTION.ln_one_minus_exp()
        }

        fn prob_emit_y(&self, _j: usize) -> LogProb {
            PROB_SUBSTITUTION.ln_one_minus_exp()
        }

        fn len_x(&self) -> usize {
            self.x.len()
        }

        fn len_y(&self) -> usize {
            self.y.len()
        }
    }

    impl Emission for TestEmissionParams {
        fn emission_x(&self, i: usize) -> u8 {
            self.x[i]
        }

        fn emission_y(&self, j: usize) -> u8 {
            self.y[j]
        }
    }

    struct TestSingleGapParams;

    impl GapParameters for TestSingleGapParams {
        fn prob_gap_x(&self) -> LogProb {
            PROB_OPEN_GAP_X
        }

        fn prob_gap_y(&self) -> LogProb {
            PROB_OPEN_GAP_Y
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            LogProb::zero()
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            LogProb::zero()
        }
    }

    struct NoGapParams;

    impl GapParameters for NoGapParams {
        fn prob_gap_x(&self) -> LogProb {
            LogProb::zero()
        }

        fn prob_gap_y(&self) -> LogProb {
            LogProb::zero()
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            LogProb::zero()
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            LogProb::zero()
        }
    }

    struct TestExtendGapParams;

    impl GapParameters for TestExtendGapParams {
        fn prob_gap_x(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_INS)
        }

        fn prob_gap_y(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_DEL)
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            T_GAP_TO_GAP
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            T_GAP_TO_GAP
        }
    }

    struct TestNoHopParams;

    impl HopParameters for TestNoHopParams {
        fn prob_hop_x(&self) -> LogProb {
            LogProb::zero()
        }

        fn prob_hop_y(&self) -> LogProb {
            LogProb::zero()
        }

        fn prob_hop_x_extend(&self) -> LogProb {
            LogProb::zero()
        }

        fn prob_hop_y_extend(&self) -> LogProb {
            LogProb::zero()
        }
    }

    struct TestHopParams;

    impl HopParameters for TestHopParams {
        fn prob_hop_x(&self) -> LogProb {
            T_MATCH_TO_HOP_X
        }

        fn prob_hop_y(&self) -> LogProb {
            T_MATCH_TO_HOP_Y
        }

        fn prob_hop_x_extend(&self) -> LogProb {
            T_HOP_X_TO_HOP_X
        }

        fn prob_hop_y_extend(&self) -> LogProb {
            T_HOP_Y_TO_HOP_Y
        }
    }

    static SINGLE_GAP_PARAMS: TestSingleGapParams = TestSingleGapParams;
    static EXTEND_GAP_PARAMS: TestExtendGapParams = TestExtendGapParams;
    static NO_GAP_PARAMS: NoGapParams = NoGapParams;
    static NO_HOP_PARAMS: TestNoHopParams = TestNoHopParams;
    static SINGLE_GAPS_NO_HOPS_PHMM: LazyLock<HomopolyPairHMM> =
        LazyLock::new(|| HomopolyPairHMM::new(&SINGLE_GAP_PARAMS, &NO_HOP_PARAMS));
    static NO_GAPS_WITH_HOPS_PHMM: LazyLock<HomopolyPairHMM> =
        LazyLock::new(|| HomopolyPairHMM::new(&NO_GAP_PARAMS, &TestHopParams));
    static EXTEND_GAPS_NO_HOPS_PHMM: LazyLock<HomopolyPairHMM> =
        LazyLock::new(|| HomopolyPairHMM::new(&EXTEND_GAP_PARAMS, &NO_HOP_PARAMS));

    /// Match states supporting the pair `(x, y)`
    fn supporting(x: u8, y: u8) -> Vec<State> {
        MATCH_STATES
            .iter()
            .copied()
            .filter(|m| m.supports(x, y))
            .collect()
    }

    #[test]
    fn supports_unambiguous_bases() {
        assert_eq!(supporting(b'A', b'A'), [MatchA]);
        assert_eq!(supporting(b'A', b'G'), [MatchA, MatchG]);
        assert_eq!(supporting(b'T', b'C'), [MatchC, MatchT]);
    }

    #[test]
    fn supports_ambiguous_codes() {
        assert_eq!(supporting(b'R', b'C'), [MatchA, MatchC, MatchG]);
        assert_eq!(supporting(b'Y', b'C'), [MatchC, MatchT]);
        assert_eq!(supporting(b'S', b'S'), [MatchC, MatchG]);
        assert_eq!(supporting(b'N', b'A'), [MatchA, MatchC, MatchG, MatchT]);
        assert_eq!(supporting(b'R', b'Y'), [MatchA, MatchC, MatchG, MatchT]);
        assert_eq!(supporting(b'r', b'c'), supporting(b'R', b'C'));
    }

    #[test]
    fn impossible_global_alignment() {
        let x = b"AAA".to_vec();
        let y = b"A".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);
        assert_eq!(p, LogProb::zero());
    }

    #[test]
    fn test_hompolymer_run_in_y() {
        let pair_hmm = &NO_GAPS_WITH_HOPS_PHMM;
        for i in 1..5 {
            let x = b"ACGT".to_vec();
            let y = format!("AC{}GT", repeat("C").take(i).join(""))
                .as_bytes()
                .to_vec();
            let emission_params = TestEmissionParams { x, y };

            let p = pair_hmm.prob_related(&emission_params, &Global, None);
            let p_most_likely_path_with_hops = LogProb(
                *EMIT_MATCH // A A
                    + *T_MATCH_TO_MATCH
                    + *EMIT_MATCH // C C
                    + *T_MATCH_TO_HOP_X // C CC
                    + *T_HOP_X_TO_HOP_X * ((i - 1) as f64)
                    + (1. - 0.1f64).ln()
                    + *EMIT_MATCH // G G
                    + *T_MATCH_TO_MATCH
                    + *EMIT_MATCH, // T T
            );
            assert!(*p <= 0.0);
            assert!(*p >= *p_most_likely_path_with_hops);
            assert!(*p < *p_most_likely_path_with_hops + 1.);
        }
    }

    #[test]
    fn test_hompolymer_run_in_x() {
        let pair_hmm = &NO_GAPS_WITH_HOPS_PHMM;
        for i in 1..5 {
            let x = format!("AC{}GT", repeat("C").take(i).join(""))
                .as_bytes()
                .to_vec();

            let y = b"ACGT".to_vec();

            let emission_params = TestEmissionParams { x, y };

            let p = pair_hmm.prob_related(&emission_params, &Global, None);
            let p_most_likely_path_with_hops = LogProb(
                *EMIT_MATCH // A A
                    + *T_MATCH_TO_MATCH
                    + *EMIT_MATCH // C C
                    + *T_MATCH_TO_HOP_Y // CC C
                    + *T_HOP_Y_TO_HOP_Y * ((i - 1) as f64)
                    + (1. - 0.1f64).ln()
                    + *EMIT_MATCH // G G
                    + *T_MATCH_TO_MATCH
                    + *EMIT_MATCH, // T T
            );
            assert!(*p <= 0.0);
            assert!(*p >= *p_most_likely_path_with_hops);
            assert!(*p < *p_most_likely_path_with_hops + 1.);
        }
    }

    #[test]
    fn test_interleave_gaps_x() {
        let x = b"AGAGAG".to_vec();
        let y = b"ACGTACGTACGT".to_vec();

        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);

        let n_matches = 6.;
        let n_insertions = 6.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_AND_Y * n_insertions
                + *T_MATCH_TO_GAP_X * n_insertions
                + *(PROB_OPEN_GAP_Y.ln_one_minus_exp()) * n_insertions,
        );

        let p_max = LogProb(*T_MATCH_TO_GAP_X * n_insertions);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_interleave_gaps_y() {
        let x = b"ACGTACGTACGT".to_vec();
        let y = b"AGAGAG".to_vec();

        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);

        let n_matches = 6.;
        let n_insertions = 6.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_X_AND_GAP * n_insertions
                + *T_MATCH_TO_GAP_Y * n_insertions
                + *PROB_OPEN_GAP_X.ln_one_minus_exp() * n_insertions,
        );

        let p_max = LogProb(*T_MATCH_TO_GAP_Y * n_insertions);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_same() {
        let x = b"AGCTCGATCGATCGATC".to_vec();
        let y = b"AGCTCGATCGATCGATC".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);
        let n = 17.;
        let p_most_likely_path = LogProb(*EMIT_MATCH * n + *T_MATCH_TO_MATCH * (n - 1.));
        let p_max = LogProb(*EMIT_MATCH * n);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.001);
        assert_relative_eq!(*p, *p_max, epsilon = 0.001);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_x() {
        let x = b"AGCTCGATCGATCGATC".to_vec();
        let y = b"AGCTCGATCTGATCGATCT".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);

        let n_matches = 17.;
        let n_insertions = 2.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_AND_Y * n_insertions
                + *T_MATCH_TO_GAP_X * n_insertions
                + (1. - *PROB_ILLUMINA_INS).ln(),
        );

        let p_max = LogProb(*T_MATCH_TO_GAP_X * 2.);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_x_2() {
        let x = b"ACAGTA".to_vec();
        let y = b"ACAGTCA".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);

        let n_matches = 6.;
        let n_insertions = 1.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_AND_Y * n_insertions
                + *T_MATCH_TO_GAP_X * n_insertions
                + (1. - *PROB_ILLUMINA_INS).ln(),
        );

        let p_max = LogProb(*T_MATCH_TO_GAP_X * n_insertions);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_y() {
        let x = b"AGCTCGATCTGATCGATCT".to_vec();
        let y = b"AGCTCGATCGATCGATC".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);

        let n_matches = 17.;
        let n_deletions = 2.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_deletions)
                + *EMIT_X_AND_GAP * n_deletions
                + *T_MATCH_TO_GAP_Y * n_deletions
                + (1. - *PROB_ILLUMINA_DEL).ln(),
        );

        let p_max = LogProb(*T_MATCH_TO_GAP_Y * 2.);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_multigap_y() {
        let x = b"AGCTCGATCTGATCGATCT".to_vec();
        let y = b"AGCTTCTGATCGATCT".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &EXTEND_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);
        let n_matches = 16.;
        let n_consecutive_deletions = 3.;
        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH_TO_MATCH * (n_matches - n_consecutive_deletions)
                + *PROB_OPEN_GAP_Y
                + *EMIT_X_AND_GAP * n_consecutive_deletions
                + *T_GAP_TO_GAP * (n_consecutive_deletions - 1.)
                + *T_GAP_TO_GAP.ln_one_minus_exp(),
        );

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
    }

    #[test]
    fn test_mismatch() {
        let x = b"AGCTCGAGCGATCGATC".to_vec();
        let y = b"TGCTCGATCGATCGATC".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Global, None);

        let n = 17.;
        let p_most_likely_path = LogProb(
            *EMIT_MATCH * (n - 2.)
                + *T_MATCH_TO_MATCH * (n - 1.)
                + (*PROB_ILLUMINA_SUBST / 3.).ln() * 2.,
        );
        let p_max = LogProb((*PROB_ILLUMINA_SUBST / 3.).ln() * 2.);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 1e-2);
        assert_relative_eq!(*p, *p_max, epsilon = 1e-1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_banded() {
        let x = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC\
ATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTAT\
CTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGT"
            .to_vec();

        let y = b"GGGTATGCACGCGATAGCATTGCGAGATGCTGGAGCTGGAGCACCCTATGTCGC".to_vec();

        let emission_params = TestEmissionParams { x, y };

        let pair_hmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&emission_params, &Semiglobal, None);

        let p_banded = pair_hmm.prob_related(&emission_params, &Semiglobal, Some(2));
        assert_relative_eq!(*p, *p_banded, epsilon = 1e-3);
    }

    /// In global mode only the origin starts with an edit distance of 0, so leading bases of x
    /// that y lacks count against the band
    #[test]
    fn test_global_band_counts_leading_bases_of_x() {
        let e = TestEmissionParams {
            x: b"GGGGACGTACGT".to_vec(),
            y: b"ACGTACGT".to_vec(),
        };
        let pair_hmm = &EXTEND_GAPS_NO_HOPS_PHMM;
        let p = pair_hmm.prob_related(&e, &Global, None);
        assert_ne!(p, LogProb::ln_zero());
        assert_eq!(
            pair_hmm.prob_related(&e, &Global, Some(1)),
            LogProb::ln_zero()
        );
        let p_wide = pair_hmm.prob_related(&e, &Global, Some(e.x.len()));
        assert_relative_eq!(*p, *p_wide, epsilon = 1e-12);
    }

    #[test]
    fn test_phmm_vs_phhmm() {
        let x = b"AGAGAGC".to_vec();
        let y = b"ATACGTACGTC".to_vec();
        let emission_params = TestEmissionParams { x, y };

        let pair_hhmm = &SINGLE_GAPS_NO_HOPS_PHMM;
        let p1 = pair_hhmm.prob_related(&emission_params, &Global, None);

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
                LogProb::from(PROB_ILLUMINA_INS)
            }

            fn prob_gap_y(&self) -> LogProb {
                LogProb::from(PROB_ILLUMINA_DEL)
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

        let mut pair_hmm = PairHMM::new(&TestSingleGapParamsPairHMM);

        let x = b"AGAGAGC";
        let y = b"ATACGTACGTC";
        let p2 = pair_hmm.prob_related(
            &TestEmissionParamsPairHMM { x, y },
            &AlignmentMode::Global,
            None,
        );
        assert_relative_eq!(*p1, *p2, epsilon = 1e-4)
    }

    /// With no hops, `HomopolyPairHMM` is a `PairHMM`, in semiglobal mode too: the sum over the
    /// column ends has to contain each column once, with the values of that column.
    #[test]
    fn test_phmm_vs_phhmm_semiglobal() {
        struct FreeEnds;
        impl crate::stats::pairhmm::StartEndGapParameters for FreeEnds {
            fn free_start_gap_x(&self) -> bool {
                true
            }

            fn free_end_gap_x(&self) -> bool {
                true
            }
        }
        impl crate::stats::pairhmm::GapParameters for FreeEnds {
            fn prob_gap_x(&self) -> LogProb {
                LogProb::from(PROB_ILLUMINA_INS)
            }

            fn prob_gap_y(&self) -> LogProb {
                LogProb::from(PROB_ILLUMINA_DEL)
            }

            fn prob_gap_x_extend(&self) -> LogProb {
                LogProb::zero()
            }

            fn prob_gap_y_extend(&self) -> LogProb {
                LogProb::zero()
            }
        }
        struct Emission {
            x: &'static [u8],
            y: &'static [u8],
        }
        impl crate::stats::pairhmm::EmissionParameters for Emission {
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
                LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
            }

            fn prob_emit_y(&self, _: usize) -> LogProb {
                LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
            }

            fn len_x(&self) -> usize {
                self.x.len()
            }

            fn len_y(&self) -> usize {
                self.y.len()
            }
        }

        // y is a window of x with one substitution, one insertion and one deletion.
        let x = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCAC";
        let y = b"TCTATCACCCTATTAACCTCTCACGGGAGCTCTCCCATGCATTTGGTATTTTCG";
        let p_homopoly = SINGLE_GAPS_NO_HOPS_PHMM.prob_related(
            &TestEmissionParams {
                x: x.to_vec(),
                y: y.to_vec(),
            },
            &Semiglobal,
            None,
        );
        let p_pair = PairHMM::new(&FreeEnds).prob_related(&Emission { x, y }, &FreeEnds, None);
        assert_relative_eq!(*p_homopoly, *p_pair, epsilon = 1e-4);
    }

    /// With base-independent parameters and emissions, renaming the bases consistently in x and y
    /// must not change the probability.
    #[test]
    fn test_base_renaming_does_not_change_probability() {
        let phmm = HomopolyPairHMM::new(&EXTEND_GAP_PARAMS, &TestHopParams);
        let x = b"ACCCAGGGTTTACGAAATCCC".to_vec();
        let y = b"ACCAGGGGTTACGAAAATCC".to_vec();
        let reference = phmm.prob_related(
            &TestEmissionParams {
                x: x.clone(),
                y: y.clone(),
            },
            &Global,
            None,
        );
        for perm in [*b"CAGT", *b"GTAC", *b"TGCA", *b"CGTA", *b"ATGC"] {
            let rename = |s: &[u8]| -> Vec<u8> {
                s.iter()
                    .map(|b| perm[b"ACGT".iter().position(|c| c == b).unwrap()])
                    .collect()
            };
            let p = phmm.prob_related(
                &TestEmissionParams {
                    x: rename(&x),
                    y: rename(&y),
                },
                &Global,
                None,
            );
            assert_relative_eq!(*p, *reference, epsilon = 1e-9);
        }
    }

    /// `prob_gap_x` opens a gap in x, i.e. a base of y emitted alone (an insertion in y), and
    /// `prob_gap_y` a gap in y, as for `PairHMM`: with insertions far more likely than deletions,
    /// an extra base in y has to be more probable than a missing one.
    #[test]
    fn test_gap_parameters_refer_to_the_sequence_with_the_gap() {
        struct InsertionsLikely;
        impl GapParameters for InsertionsLikely {
            fn prob_gap_x(&self) -> LogProb {
                LogProb::from(Prob(1e-2))
            }

            fn prob_gap_y(&self) -> LogProb {
                LogProb::from(Prob(1e-6))
            }

            fn prob_gap_x_extend(&self) -> LogProb {
                LogProb::zero()
            }

            fn prob_gap_y_extend(&self) -> LogProb {
                LogProb::zero()
            }
        }
        let pair_hmm = HomopolyPairHMM::new(&InsertionsLikely, &NO_HOP_PARAMS);
        let x = b"ACGTTGCAGT".to_vec();
        let insertion = TestEmissionParams {
            x: x.clone(),
            y: b"ACGTTGACAGT".to_vec(),
        };
        let deletion = TestEmissionParams {
            x,
            y: b"ACGTTGAGT".to_vec(),
        };
        let p_insertion = pair_hmm.prob_related(&insertion, &Global, None);
        let p_deletion = pair_hmm.prob_related(&deletion, &Global, None);
        assert!(
            p_insertion > p_deletion,
            "insertion {} should be more likely than deletion {}",
            *p_insertion,
            *p_deletion
        );
    }

    /// Reference values computed before the transition table became a dense array.
    #[test]
    fn test_values_are_unchanged() {
        let windows: [(&[u8], &[u8]); 3] = [
            (b"ACCCAGGGTTTACGAAATCCCGATT", b"ACCAGGGGTTACGAAAATCCGATT"),
            (b"GATTTACAGGGGCATTTYRACCCGGNT", b"TTACAGGGCATTTCAACCGGT"),
            (b"TTTTTTACGTACGTAAAAAAGGGCCC", b"TTTTTACGTACGTAAAAAGGGGCCCA"),
        ];
        let phmm = HomopolyPairHMM::new(&EXTEND_GAP_PARAMS, &TestHopParams);
        let mut values = Vec::new();
        for (x, y) in windows {
            for (mode_is_global, band) in [(true, None), (false, None), (false, Some(6))] {
                let e = TestEmissionParams {
                    x: x.to_vec(),
                    y: y.to_vec(),
                };
                let p = if mode_is_global {
                    phmm.prob_related(&e, &Global, band)
                } else {
                    phmm.prob_related(&e, &Semiglobal, band)
                };
                values.push(*p);
            }
        }
        let expected = [
            -32.086378685365,
            -29.001106194398,
            -29.001106194398,
            -75.954685092471,
            -38.906771826337,
            -38.906771826337,
            -29.588020202579,
            -20.095429251557,
            -20.095429251557,
        ];
        for (v, e) in values.iter().zip(expected) {
            assert_relative_eq!(*v, e, epsilon = 1e-9);
        }
    }
}
