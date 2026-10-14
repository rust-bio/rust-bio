// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::alignment::{ProfileAlignment, ProfileOp};
use crate::stats::hmm::profile::decode::path_log_prob;
use crate::stats::hmm::profile::forward_backward::forward_matrices;
use crate::stats::hmm::profile::{ProfileError, ProfileHmm, NEG_INF};
use std::f64::consts::LN_2;

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Rng {
        Rng {
            state: seed ^ 0x9e3779b97f4a7c15,
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut z = self.state.wrapping_add(0x9e3779b97f4a7c15);
        self.state = z;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn choose(weights: &[f64], draw: f64) -> usize {
    let total: f64 = weights.iter().sum();
    if total <= 0.0 {
        return weights.len() - 1;
    }
    let mut target = draw * total;
    for (idx, &w) in weights.iter().enumerate() {
        if target < w {
            return idx;
        }
        target -= w;
    }
    weights.len() - 1
}

const STATE_M: u8 = 0;
const STATE_I: u8 = 1;
const STATE_D: u8 = 2;

pub(crate) fn sample_alignment(
    profile: &ProfileHmm,
    query: &[u8],
    seed: u64,
) -> Result<ProfileAlignment, ProfileError> {
    let ranks = profile.ranks(query)?;
    let n = ranks.len();
    let l = profile.num_cols;
    let fwd = forward_matrices(profile, &ranks);
    let mut rng = Rng::new(seed);

    let end_weights = vec![
        (fwd.m[n][l] + profile.t_mm[l]).exp(),
        (fwd.i[n][l] + profile.t_im[l]).exp(),
        (fwd.d[n][l] + profile.t_dm[l]).exp(),
    ];
    let mut state = match choose(&end_weights, rng.unit()) {
        0 => STATE_M,
        1 => STATE_I,
        _ => STATE_D,
    };

    let mut ops = Vec::new();
    let mut i = n;
    let mut k = l;
    loop {
        if state == STATE_M && k == 0 && i == 0 {
            break;
        }
        match state {
            STATE_M => {
                ops.push(ProfileOp::Match(query[i - 1]));
                let weights = vec![
                    (fwd.m[i - 1][k - 1] + profile.t_mm[k - 1]).exp(),
                    (fwd.i[i - 1][k - 1] + profile.t_im[k - 1]).exp(),
                    (fwd.d[i - 1][k - 1] + profile.t_dm[k - 1]).exp(),
                ];
                state = pick_state(&weights, &mut rng);
                i -= 1;
                k -= 1;
            }
            STATE_I => {
                ops.push(ProfileOp::Insert(query[i - 1]));
                let weights = vec![
                    (fwd.m[i - 1][k] + profile.t_mi[k]).exp(),
                    (fwd.i[i - 1][k] + profile.t_ii[k]).exp(),
                    (fwd.d[i - 1][k] + profile.t_di[k]).exp(),
                ];
                state = pick_state(&weights, &mut rng);
                i -= 1;
            }
            _ => {
                ops.push(ProfileOp::Delete);
                let weights = vec![
                    (fwd.m[i][k - 1] + profile.t_md[k - 1]).exp(),
                    (fwd.i[i][k - 1] + profile.t_id[k - 1]).exp(),
                    (fwd.d[i][k - 1] + profile.t_dd[k - 1]).exp(),
                ];
                state = pick_state(&weights, &mut rng);
                k -= 1;
            }
        }
    }
    ops.reverse();

    let log_prob = path_log_prob(profile, &ops)?;
    let null_lp = profile.null_log_prob_ranks(&ranks);
    let bit_score = if log_prob == NEG_INF {
        NEG_INF
    } else {
        (log_prob - null_lp) / LN_2
    };
    Ok(ProfileAlignment {
        ops,
        log_prob,
        bit_score,
    })
}

fn pick_state(weights: &[f64], rng: &mut Rng) -> u8 {
    match choose(weights, rng.unit()) {
        0 => STATE_M,
        1 => STATE_I,
        _ => STATE_D,
    }
}

pub(crate) fn generate(profile: &ProfileHmm, seed: u64) -> Vec<u8> {
    let l = profile.num_cols;
    let mut rng = Rng::new(seed);
    let mut out = Vec::new();
    let mut state = STATE_M;
    let mut node = 0usize;

    loop {
        let (to_m, to_i, to_d) = match state {
            STATE_M => (profile.t_mm[node], profile.t_mi[node], profile.t_md[node]),
            STATE_I => (profile.t_im[node], profile.t_ii[node], profile.t_id[node]),
            _ => (profile.t_dm[node], profile.t_di[node], profile.t_dd[node]),
        };
        let weights = vec![to_m.exp(), to_i.exp(), to_d.exp()];
        match choose(&weights, rng.unit()) {
            0 => {
                if node == l {
                    break;
                }
                node += 1;
                out.push(emit(profile, &profile.emission[node], &mut rng));
                state = STATE_M;
            }
            1 => {
                out.push(emit(profile, &profile.null, &mut rng));
                state = STATE_I;
            }
            _ => {
                if node >= l {
                    break;
                }
                node += 1;
                state = STATE_D;
            }
        }
        if out.len() > 100_000 {
            break;
        }
    }
    out
}

fn emit(profile: &ProfileHmm, log_probs: &[f64], rng: &mut Rng) -> u8 {
    let weights: Vec<f64> = log_probs.iter().map(|&p| p.exp()).collect();
    let idx = choose(&weights, rng.unit());
    profile.symbols[idx]
}

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;

    #[test]
    fn test_sampled_alignment_is_reproducible() {
        let profile = example_profile();
        let query = b"ACGTAC";
        let first = profile.sample_alignment(query, 1234).unwrap();
        let second = profile.sample_alignment(query, 1234).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.consumed_residues(), query.len());
        assert!(first.log_prob <= profile.best_alignment(query).unwrap().log_prob + 1e-9);
    }

    #[test]
    fn test_sampled_alignment_covers_every_position() {
        let profile = example_profile();
        let query = b"ACAC";
        for seed in 0..20 {
            let alignment = profile.sample_alignment(query, seed).unwrap();
            assert_eq!(alignment.consumed_residues(), query.len());
            assert_eq!(
                alignment.matched_columns() + alignment.deleted_columns(),
                profile.num_columns()
            );
            let score = profile.score_path(&alignment.ops).unwrap();
            assert!((score - alignment.log_prob).abs() < 1e-9);
        }
    }

    #[test]
    fn test_generate_is_reproducible() {
        let profile = example_profile();
        assert_eq!(profile.generate(77), profile.generate(77));
        for seed in 0..20 {
            let sequence = profile.generate(seed);
            assert!(profile.forward(&sequence).unwrap().is_finite());
        }
    }
}
