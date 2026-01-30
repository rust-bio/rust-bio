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
use std::iter::once;
use std::mem;
use std::ops::Shr;
use usize;

use enum_map::{Enum, EnumMap};
use itertools::Itertools;
use num_traits::Zero;

use crate::stats::pairhmm::homopolypairhmm::State::*;
use crate::stats::pairhmm::{
    Emission, EmissionParameters, GapParameters, StartEndGapParameters, XYEmission,
};
use crate::stats::probs::LogProb;
use crate::stats::Prob;
use std::collections::HashMap;

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
    fn supports(&self, x: u8, y: u8) -> bool {
        match self {
            MatchA if x == b'A' || y == b'A' => true,
            MatchC if x == b'C' || y == b'C' => true,
            MatchG if x == b'G' || y == b'G' => true,
            MatchT if x == b'T' || y == b'T' => true,
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

const STATES: [State; 14] = [
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
#[derive(Default, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct HomopolyPairHMM {
    transition_probs: HashMap<usize, LogProb>,
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

        let len_y = emission_params.len_y();
        let len_x = emission_params.len_x();
        let mut min_edit_dist: [Vec<usize>; 2] =
            [vec![usize::MAX; len_y + 1], vec![usize::MAX; len_y + 1]];
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
                let mut any_match = false;
                for &m in &MATCH_STATES {
                    if m.supports(emission_x, emission_y) {
                        let emission = emission_params.prob_emit_xy(i, j);
                        let emission_prob = match emission {
                            XYEmission::Match(p) => p,
                            // since we have separate match states, we need to halve mismatch probs
                            // (since e.g. ('A', _) and (_, 'A') are distinct cases)
                            XYEmission::Mismatch(p) => LogProb::from(*p - 2f64.ln()),
                        };
                        any_match |= emission.is_match();
                        v[curr][m][j_] = emission_prob
                            + LogProb::ln_sum_exp(
                                &STATES
                                    .iter()
                                    .map(|&s| {
                                        transition_probs.get(&(s >> m)).unwrap_or(&LogProb::zero())
                                            + v[prev][s][j_minus_one]
                                    })
                                    .collect_vec(),
                            );
                    } else {
                        v[curr][m][j_] = LogProb::zero();
                    }
                }

                v[curr][GapY][j_] = prob_emit_x_and_gap
                    + LogProb::ln_sum_exp(
                        &MATCH_STATES
                            .iter()
                            .map(|&s| transition_probs[&(s >> GapY)] + v[prev][s][j_])
                            .chain(once(transition_probs[&(GapY >> GapY)] + v[prev][GapY][j_]))
                            .collect_vec(),
                    );

                MATCH_HOP_Y.iter().for_each(|&(m, h)| {
                    v[curr][h][j_] = (transition_probs[&(m >> h)] + v[prev][m][j_])
                        .ln_add_exp(transition_probs[&(h >> h)] + v[prev][h][j_])
                });

                v[curr][GapX][j_] = emission_params.prob_emit_y(j)
                    + LogProb::ln_sum_exp(
                        &MATCH_STATES
                            .iter()
                            .map(|&s| transition_probs[&(s >> GapX)] + v[curr][s][j_minus_one])
                            .chain(once(
                                transition_probs[&(GapX >> GapX)] + v[curr][GapX][j_minus_one],
                            ))
                            .collect_vec(),
                    );

                MATCH_HOP_X.iter().for_each(|&(m, h)| {
                    v[curr][h][j_] = (transition_probs[&(m >> h)] + v[curr][m][j_minus_one])
                        .ln_add_exp(transition_probs[&(h >> h)] + v[curr][h][j_minus_one])
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
                v[curr][s].reset(LogProb::zero());
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
    (HopCX, MatchC),
    (HopCX, MatchC),
    (HopCX, MatchG),
    (HopCX, MatchT),
    (HopGX, MatchG),
    (HopGX, MatchC),
    (HopGX, MatchG),
    (HopGX, MatchT),
    (HopTX, MatchT),
    (HopTX, MatchC),
    (HopTX, MatchG),
    (HopTX, MatchT),
];
const HOP_Y_MATCH: [(State, State); 16] = [
    (HopAY, MatchA),
    (HopAY, MatchC),
    (HopAY, MatchG),
    (HopAY, MatchT),
    (HopCY, MatchC),
    (HopCY, MatchC),
    (HopCY, MatchG),
    (HopCY, MatchT),
    (HopGY, MatchG),
    (HopGY, MatchC),
    (HopGY, MatchG),
    (HopGY, MatchT),
    (HopTY, MatchT),
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
) -> HashMap<usize, LogProb> {
    let mut transition_probs = HashMap::new();

    let prob_gap_x = gap_params.prob_gap_x();
    let prob_gap_y = gap_params.prob_gap_y();
    let prob_gap_x_extend = gap_params.prob_gap_x_extend();
    let prob_gap_y_extend = gap_params.prob_gap_y_extend();

    MATCH_HOP_X.iter().for_each(|(a, b)| {
        transition_probs.insert(
            *a >> *b,
            hop_params.prob_hop_x_with_base(b.base().expect("Unsupported base")),
        );
    });
    MATCH_HOP_Y.iter().for_each(|(a, b)| {
        transition_probs.insert(
            *a >> *b,
            hop_params.prob_hop_y_with_base(b.base().expect("Unsupported base")),
        );
    });
    HOP_X_HOP_X.iter().for_each(|(a, b)| {
        assert_eq!(a.base(), b.base());
        transition_probs.insert(
            *a >> *b,
            hop_params.prob_hop_x_extend_with_base(b.base().expect("Unsupported base")),
        );
    });
    HOP_Y_HOP_Y.iter().for_each(|(a, b)| {
        assert_eq!(a.base(), b.base());
        transition_probs.insert(
            *a >> *b,
            hop_params.prob_hop_y_extend_with_base(b.base().expect("Unsupported base")),
        );
    });
    HOP_X_MATCH.iter().for_each(|(a, b)| {
        transition_probs.insert(
            *a >> *b,
            hop_params
                .prob_hop_x_with_base(a.base().expect("Unsupported base"))
                .ln_one_minus_exp(),
        );
    });
    HOP_Y_MATCH.iter().for_each(|(a, b)| {
        transition_probs.insert(
            *a >> *b,
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
        transition_probs.insert(*a >> *b, match_same);
    });
    MATCH_OTHER.iter().for_each(|(a, b)| {
        transition_probs.insert(*a >> *b, match_other);
    });

    MATCH_STATES.iter().for_each(|&a| {
        transition_probs.insert(a >> GapX, prob_gap_y);
    });
    MATCH_STATES.iter().for_each(|&a| {
        transition_probs.insert(a >> GapY, prob_gap_x);
    });
    MATCH_STATES.iter().for_each(|&b| {
        transition_probs.insert(GapX >> b, prob_gap_y_extend.ln_one_minus_exp());
    });
    MATCH_STATES.iter().for_each(|&b| {
        transition_probs.insert(GapY >> b, prob_gap_x_extend.ln_one_minus_exp());
    });
    transition_probs.insert(GapX >> GapX, prob_gap_y_extend);
    transition_probs.insert(GapY >> GapY, prob_gap_x_extend);
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
    use itertools::repeat_n;

    use crate::stats::pairhmm::homopolypairhmm::tests::AlignmentMode::{Global, Semiglobal};
    use crate::stats::pairhmm::PairHMM;
    use crate::stats::{LogProb, Prob};
    use std::sync::LazyLock;

    use super::*;

    // Single base insertion and deletion rates for R1 according to Schirmer et al.
    // BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
    static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
    static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
    static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

    // log(0.0021)
    const PROB_SUBSTITUTION: LogProb = LogProb(-6.165_817_934_252_76);
    // log(2.8e-6)
    const PROB_OPEN_GAP_Y: LogProb = LogProb(-12.785_891_140_783_116);
    // log(5.1e-6)
    const PROB_OPEN_GAP_X: LogProb = LogProb(-12.186_270_018_233_994);

    const EMIT_MATCH: LogProb = LogProb(-0.0021022080918701985);
    const EMIT_GAP_AND_Y: LogProb = LogProb(-0.0021022080918701985);
    const EMIT_X_AND_GAP: LogProb = LogProb(-0.0021022080918701985);

    const T_MATCH_TO_HOP_X: LogProb = LogProb(-11.512925464970229);
    const T_MATCH_TO_HOP_Y: LogProb = LogProb(-11.512925464970229);
    const T_HOP_X_TO_HOP_X: LogProb = LogProb(-2.3025850929940455);
    const T_HOP_Y_TO_HOP_Y: LogProb = LogProb(-2.3025850929940455);

    const T_MATCH_TO_MATCH: LogProb = LogProb(-7.900_031_205_113_962e-6);
    const T_MATCH_TO_GAP_Y: LogProb = LogProb(-12.785_891_140_783_116);
    const T_MATCH_TO_GAP_X: LogProb = LogProb(-12.186_270_018_233_994);
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
            PROB_OPEN_GAP_Y
        }

        fn prob_gap_y(&self) -> LogProb {
            PROB_OPEN_GAP_X
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
            let y = format!("AC{}GT", repeat_n("C", i).collect::<String>())
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
            let x = format!("AC{}GT", repeat_n("C", i).collect::<String>())
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
}
