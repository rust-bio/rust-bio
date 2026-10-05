// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::alignment::{ProfileAlignment, ProfileOp};
use crate::stats::hmm::profile::forward_backward::{posterior_emit_matrix, posterior_match};
use crate::stats::hmm::profile::{ProfileError, ProfileHmm, NEG_INF};
use std::f64::consts::LN_2;

const STEP_NONE: u8 = 0;
const STEP_MATCH: u8 = 1;
const STEP_DELETE: u8 = 2;
const STEP_INSERT: u8 = 3;

pub(crate) fn optimal_accuracy(
    profile: &ProfileHmm,
    query: &[u8],
) -> Result<ProfileAlignment, ProfileError> {
    let ranks = profile.ranks(query)?;
    let (ops, _) = mea_dp(profile, &ranks);
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

pub(crate) fn expected_accuracy(profile: &ProfileHmm, query: &[u8]) -> Result<f64, ProfileError> {
    let ranks = profile.ranks(query)?;
    let (_, accuracy) = mea_dp(profile, &ranks);
    Ok(accuracy)
}

pub(crate) fn posterior_path(profile: &ProfileHmm, query: &[u8]) -> Result<Vec<f64>, ProfileError> {
    let post = posterior_match(profile, query)?;
    let alignment = profile.best_alignment(query)?;
    let mut out = Vec::new();
    let mut column = 0usize;
    for op in &alignment.ops {
        match op {
            ProfileOp::Match(_) => {
                column += 1;
                out.push(post[column - 1]);
            }
            ProfileOp::Delete => {
                column += 1;
            }
            ProfileOp::Insert(_) => {}
        }
    }
    Ok(out)
}

fn mea_dp(profile: &ProfileHmm, ranks: &[usize]) -> (Vec<ProfileOp>, f64) {
    let n = ranks.len();
    let l = profile.num_cols;
    let post = posterior_emit_matrix(profile, ranks);

    let mut score = vec![vec![0.0f64; l + 1]; n + 1];
    let mut back = vec![vec![STEP_NONE; l + 1]; n + 1];
    for row in back.iter_mut().skip(1) {
        row[0] = STEP_INSERT;
    }
    for cell in back[0].iter_mut().skip(1) {
        *cell = STEP_DELETE;
    }

    for i in 1..=n {
        for k in 1..=l {
            let diag = score[i - 1][k - 1] + post[i][k];
            let up = score[i - 1][k];
            let left = score[i][k - 1];
            let mut best = diag;
            let mut step = STEP_MATCH;
            if up > best {
                best = up;
                step = STEP_INSERT;
            }
            if left > best {
                best = left;
                step = STEP_DELETE;
            }
            score[i][k] = best;
            back[i][k] = step;
        }
    }

    let mut ops = Vec::new();
    let mut i = n;
    let mut k = l;
    while i > 0 || k > 0 {
        match back[i][k] {
            STEP_MATCH => {
                ops.push(ProfileOp::Match(profile.symbols[ranks[i - 1]]));
                i -= 1;
                k -= 1;
            }
            STEP_INSERT => {
                ops.push(ProfileOp::Insert(profile.symbols[ranks[i - 1]]));
                i -= 1;
            }
            STEP_DELETE => {
                ops.push(ProfileOp::Delete);
                k -= 1;
            }
            _ => break,
        }
    }
    ops.reverse();
    (ops, score[n][l])
}

pub(crate) fn path_log_prob(profile: &ProfileHmm, ops: &[ProfileOp]) -> Result<f64, ProfileError> {
    let l = profile.num_cols;
    let mut total = 0.0f64;
    let mut node = 0usize;
    let mut prev = b'M';
    for op in ops {
        match op {
            ProfileOp::Match(sym) => {
                let trans = match prev {
                    b'M' => profile.t_mm[node],
                    b'I' => profile.t_im[node],
                    _ => profile.t_dm[node],
                };
                let rank = profile.rank(*sym)?;
                total += trans + profile.emission[node + 1][rank];
                node += 1;
                prev = b'M';
            }
            ProfileOp::Insert(sym) => {
                let trans = match prev {
                    b'M' => profile.t_mi[node],
                    b'I' => profile.t_ii[node],
                    _ => profile.t_di[node],
                };
                let rank = profile.rank(*sym)?;
                total += trans + profile.null[rank];
                prev = b'I';
            }
            ProfileOp::Delete => {
                let trans = match prev {
                    b'M' => profile.t_md[node],
                    b'I' => profile.t_id[node],
                    _ => profile.t_dd[node],
                };
                total += trans;
                node += 1;
                prev = b'D';
            }
        }
        if total == NEG_INF {
            return Ok(NEG_INF);
        }
    }
    if node != l {
        return Ok(NEG_INF);
    }
    let end = match prev {
        b'M' => profile.t_mm[l],
        b'I' => profile.t_im[l],
        _ => profile.t_dm[l],
    };
    total += end;
    Ok(total)
}

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;
    use crate::stats::hmm::profile::ProfileOp;

    #[test]
    fn test_optimal_accuracy_is_not_more_probable_than_viterbi() {
        let profile = example_profile();
        for query in [
            &b"CAAGCA"[..],
            &b"ACGTAC"[..],
            &b"ACAC"[..],
            &b"AGTAC"[..],
            &b"GTACCG"[..],
            &b"ACGTTTAC"[..],
        ] {
            let best = profile.best_alignment(query).unwrap();
            let optimal = profile.optimal_accuracy(query).unwrap();
            assert_eq!(optimal.consumed_residues(), query.len());
            assert!(
                optimal.log_prob <= best.log_prob + 1e-9,
                "optimal accuracy {} exceeds viterbi {} for {:?}",
                optimal.log_prob,
                best.log_prob,
                query
            );
        }
    }

    #[test]
    fn test_expected_accuracy_is_bounded() {
        let profile = example_profile();
        for query in [
            &b"ACGTAC"[..],
            &b"ACAC"[..],
            &b"ACGTTTAC"[..],
            &b"CAAGCA"[..],
        ] {
            let accuracy = profile.expected_accuracy(query).unwrap();
            let posterior_sum: f64 = profile.posterior_match(query).unwrap().iter().sum();
            assert!(accuracy >= 0.0 && accuracy <= profile.num_columns() as f64 + 1e-9);
            assert!(
                accuracy <= posterior_sum + 1e-9,
                "expected accuracy {} exceeds summed match probability {}",
                accuracy,
                posterior_sum
            );
        }
    }

    #[test]
    fn test_posterior_path_follows_matched_columns() {
        let profile = example_profile();
        let query = b"ACAC";
        let best = profile.best_alignment(query).unwrap();
        let posterior = profile.posterior_match(query).unwrap();
        let mut expected = Vec::new();
        let mut column = 0;
        for op in &best.ops {
            match op {
                ProfileOp::Match(_) => {
                    column += 1;
                    expected.push(posterior[column - 1]);
                }
                ProfileOp::Delete => column += 1,
                ProfileOp::Insert(_) => {}
            }
        }
        let path = profile.posterior_path(query).unwrap();
        assert_eq!(path.len(), best.matched_columns());
        assert_eq!(path, expected);
    }

    #[test]
    fn test_score_path_requires_every_position() {
        let profile = example_profile();
        let ops = vec![ProfileOp::Match(b'A'), ProfileOp::Match(b'C')];
        assert_eq!(profile.score_path(&ops).unwrap(), f64::NEG_INFINITY);
    }
}
