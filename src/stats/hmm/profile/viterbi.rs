// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::alignment::{ProfileAlignment, ProfileOp};
use crate::stats::hmm::profile::{ProfileError, ProfileHmm, NEG_INF};
use std::f64::consts::LN_2;

const FROM_NONE: u8 = 255;
const FROM_M: u8 = 0;
const FROM_I: u8 = 1;
const FROM_D: u8 = 2;

fn priority(state: u8) -> u8 {
    match state {
        FROM_M => 2,
        FROM_D => 1,
        _ => 0,
    }
}

fn best(candidates: &[(f64, u8)]) -> (f64, u8) {
    let mut best_val = NEG_INF;
    let mut best_state = FROM_NONE;
    for &(val, state) in candidates {
        if val == NEG_INF {
            continue;
        }
        if best_state == FROM_NONE
            || val > best_val
            || (val == best_val && priority(state) > priority(best_state))
        {
            best_val = val;
            best_state = state;
        }
    }
    (best_val, best_state)
}

pub(crate) fn viterbi(
    profile: &ProfileHmm,
    query: &[u8],
) -> Result<ProfileAlignment, ProfileError> {
    let ranks = profile.ranks(query)?;
    let n = ranks.len();
    let l = profile.num_cols;

    let mut vm = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut vi = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut vd = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut bm = vec![vec![FROM_NONE; l + 1]; n + 1];
    let mut bi = vec![vec![FROM_NONE; l + 1]; n + 1];
    let mut bd = vec![vec![FROM_NONE; l + 1]; n + 1];

    vm[0][0] = 0.0;
    for k in 1..=l {
        let (val, from) = best(&[
            (vm[0][k - 1] + profile.t_md[k - 1], FROM_M),
            (vi[0][k - 1] + profile.t_id[k - 1], FROM_I),
            (vd[0][k - 1] + profile.t_dd[k - 1], FROM_D),
        ]);
        vd[0][k] = val;
        bd[0][k] = from;
    }

    for i in 1..=n {
        let x = ranks[i - 1];
        for k in 1..=l {
            let (val, from) = best(&[
                (vm[i - 1][k - 1] + profile.t_mm[k - 1], FROM_M),
                (vi[i - 1][k - 1] + profile.t_im[k - 1], FROM_I),
                (vd[i - 1][k - 1] + profile.t_dm[k - 1], FROM_D),
            ]);
            if val != NEG_INF {
                vm[i][k] = val + profile.emission[k][x];
                bm[i][k] = from;
            }
        }
        for k in 0..=l {
            let (val, from) = best(&[
                (vm[i - 1][k] + profile.t_mi[k], FROM_M),
                (vi[i - 1][k] + profile.t_ii[k], FROM_I),
                (vd[i - 1][k] + profile.t_di[k], FROM_D),
            ]);
            if val != NEG_INF {
                vi[i][k] = val + profile.null[x];
                bi[i][k] = from;
            }
        }
        for k in 1..=l {
            let (val, from) = best(&[
                (vm[i][k - 1] + profile.t_md[k - 1], FROM_M),
                (vi[i][k - 1] + profile.t_id[k - 1], FROM_I),
                (vd[i][k - 1] + profile.t_dd[k - 1], FROM_D),
            ]);
            vd[i][k] = val;
            bd[i][k] = from;
        }
    }

    let (total, end_from) = best(&[
        (vm[n][l] + profile.t_mm[l], FROM_M),
        (vi[n][l] + profile.t_im[l], FROM_I),
        (vd[n][l] + profile.t_dm[l], FROM_D),
    ]);

    let mut ops = Vec::new();
    let mut state = end_from;
    let mut i = n;
    let mut k = l;
    while !(state == FROM_M && k == 0 && i == 0) {
        match state {
            FROM_M => {
                ops.push(ProfileOp::Match(query[i - 1]));
                let from = bm[i][k];
                i -= 1;
                k -= 1;
                state = from;
            }
            FROM_I => {
                ops.push(ProfileOp::Insert(query[i - 1]));
                let from = bi[i][k];
                i -= 1;
                state = from;
            }
            FROM_D => {
                ops.push(ProfileOp::Delete);
                let from = bd[i][k];
                k -= 1;
                state = from;
            }
            _ => break,
        }
    }
    ops.reverse();

    let null_lp = profile.null_log_prob_ranks(&ranks);
    let bit_score = (total - null_lp) / LN_2;
    Ok(ProfileAlignment {
        ops,
        log_prob: total,
        bit_score,
    })
}

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;
    use crate::stats::hmm::profile::ProfileOp;

    #[test]
    fn test_exact_match_aligns_every_position() {
        let profile = example_profile();
        let alignment = profile.best_alignment(b"ACGTAC").unwrap();
        assert_eq!(alignment.matched_columns(), profile.num_columns());
        assert_eq!(alignment.deleted_columns(), 0);
        assert_eq!(alignment.consumed_residues(), 6);
        let score = profile.score_path(&alignment.ops).unwrap();
        assert!((alignment.log_prob - score).abs() < 1e-9);
    }

    #[test]
    fn test_insert_run_is_scored_with_self_loop() {
        let profile = example_profile();
        let query = b"ACGTAGGC";
        let alignment = profile.best_alignment(query).unwrap();
        let consecutive_inserts = vec![
            ProfileOp::Match(b'A'),
            ProfileOp::Match(b'C'),
            ProfileOp::Match(b'G'),
            ProfileOp::Match(b'T'),
            ProfileOp::Match(b'A'),
            ProfileOp::Insert(b'G'),
            ProfileOp::Insert(b'G'),
            ProfileOp::Match(b'C'),
        ];
        let baseline = profile.score_path(&consecutive_inserts).unwrap();
        assert!(baseline.is_finite());
        assert!(alignment.log_prob >= baseline - 1e-9);
        assert_eq!(alignment.consumed_residues(), query.len());
    }

    #[test]
    fn test_delete_run_is_scored_with_self_loop() {
        let profile = example_profile();
        let query = b"ACAC";
        let alignment = profile.best_alignment(query).unwrap();
        let consecutive_deletes = vec![
            ProfileOp::Match(b'A'),
            ProfileOp::Match(b'C'),
            ProfileOp::Delete,
            ProfileOp::Delete,
            ProfileOp::Match(b'A'),
            ProfileOp::Match(b'C'),
        ];
        let baseline = profile.score_path(&consecutive_deletes).unwrap();
        assert!(baseline.is_finite());
        assert!(alignment.log_prob >= baseline - 1e-9);
        assert!(alignment.deleted_columns() >= profile.num_columns() - query.len());
    }

    #[test]
    fn test_empty_query_deletes_every_position() {
        let profile = example_profile();
        let alignment = profile.best_alignment(b"").unwrap();
        assert_eq!(alignment.deleted_columns(), profile.num_columns());
        assert_eq!(alignment.consumed_residues(), 0);
    }
}
