// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::{ln_add, ProfileError, ProfileHmm, NEG_INF};

pub(crate) struct Matrices {
    pub m: Vec<Vec<f64>>,
    pub i: Vec<Vec<f64>>,
    pub d: Vec<Vec<f64>>,
    pub total: f64,
}

fn fold(terms: &[f64]) -> f64 {
    let mut acc = NEG_INF;
    for &t in terms {
        acc = ln_add(acc, t);
    }
    acc
}

pub(crate) fn forward_matrices(profile: &ProfileHmm, ranks: &[usize]) -> Matrices {
    let n = ranks.len();
    let l = profile.num_cols;
    let mut fm = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut fi = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut fd = vec![vec![NEG_INF; l + 1]; n + 1];

    fm[0][0] = 0.0;
    for k in 1..=l {
        fd[0][k] = fold(&[
            fm[0][k - 1] + profile.t_md[k - 1],
            fi[0][k - 1] + profile.t_id[k - 1],
            fd[0][k - 1] + profile.t_dd[k - 1],
        ]);
    }

    for i in 1..=n {
        let x = ranks[i - 1];
        for k in 1..=l {
            let inc = fold(&[
                fm[i - 1][k - 1] + profile.t_mm[k - 1],
                fi[i - 1][k - 1] + profile.t_im[k - 1],
                fd[i - 1][k - 1] + profile.t_dm[k - 1],
            ]);
            fm[i][k] = if inc == NEG_INF {
                NEG_INF
            } else {
                inc + profile.emission[k][x]
            };
        }
        for k in 0..=l {
            let inc = fold(&[
                fm[i - 1][k] + profile.t_mi[k],
                fi[i - 1][k] + profile.t_ii[k],
                fd[i - 1][k] + profile.t_di[k],
            ]);
            fi[i][k] = if inc == NEG_INF {
                NEG_INF
            } else {
                inc + profile.null[x]
            };
        }
        for k in 1..=l {
            fd[i][k] = fold(&[
                fm[i][k - 1] + profile.t_md[k - 1],
                fi[i][k - 1] + profile.t_id[k - 1],
                fd[i][k - 1] + profile.t_dd[k - 1],
            ]);
        }
    }

    let total = fold(&[
        fm[n][l] + profile.t_mm[l],
        fi[n][l] + profile.t_im[l],
        fd[n][l] + profile.t_dm[l],
    ]);
    Matrices {
        m: fm,
        i: fi,
        d: fd,
        total,
    }
}

pub(crate) fn backward_matrices(profile: &ProfileHmm, ranks: &[usize]) -> Matrices {
    let n = ranks.len();
    let l = profile.num_cols;
    let mut bm = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut bi = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut bd = vec![vec![NEG_INF; l + 1]; n + 1];

    for i in (0..=n).rev() {
        let next = if i < n { Some(ranks[i]) } else { None };
        for k in (0..=l).rev() {
            let mut m_terms: Vec<f64> = Vec::new();
            let mut i_terms: Vec<f64> = Vec::new();
            let mut d_terms: Vec<f64> = Vec::new();

            if let Some(y) = next {
                m_terms.push(profile.t_mi[k] + profile.null[y] + bi[i + 1][k]);
                i_terms.push(profile.t_ii[k] + profile.null[y] + bi[i + 1][k]);
            }
            if k < l {
                if let Some(y) = next {
                    m_terms.push(profile.t_mm[k] + profile.emission[k + 1][y] + bm[i + 1][k + 1]);
                    i_terms.push(profile.t_im[k] + profile.emission[k + 1][y] + bm[i + 1][k + 1]);
                    d_terms.push(profile.t_dm[k] + profile.emission[k + 1][y] + bm[i + 1][k + 1]);
                    d_terms.push(profile.t_di[k] + profile.null[y] + bi[i + 1][k]);
                }
                m_terms.push(profile.t_md[k] + bd[i][k + 1]);
                i_terms.push(profile.t_id[k] + bd[i][k + 1]);
                d_terms.push(profile.t_dd[k] + bd[i][k + 1]);
            } else {
                if i == n {
                    m_terms.push(profile.t_mm[l]);
                    i_terms.push(profile.t_im[l]);
                    if k >= 1 {
                        d_terms.push(profile.t_dm[l]);
                    }
                }
                if let Some(y) = next {
                    if k >= 1 {
                        d_terms.push(profile.t_di[l] + profile.null[y] + bi[i + 1][l]);
                    }
                }
            }

            bm[i][k] = fold(&m_terms);
            bi[i][k] = fold(&i_terms);
            if k >= 1 {
                bd[i][k] = fold(&d_terms);
            }
        }
    }

    let total = bm[0][0];
    Matrices {
        m: bm,
        i: bi,
        d: bd,
        total,
    }
}

pub(crate) fn forward(profile: &ProfileHmm, query: &[u8]) -> Result<f64, ProfileError> {
    let ranks = profile.ranks(query)?;
    Ok(forward_matrices(profile, &ranks).total)
}

pub(crate) fn backward(profile: &ProfileHmm, query: &[u8]) -> Result<f64, ProfileError> {
    let ranks = profile.ranks(query)?;
    Ok(backward_matrices(profile, &ranks).total)
}

pub(crate) fn posterior_match(
    profile: &ProfileHmm,
    query: &[u8],
) -> Result<Vec<f64>, ProfileError> {
    let ranks = profile.ranks(query)?;
    let n = ranks.len();
    let l = profile.num_cols;
    let fwd = forward_matrices(profile, &ranks);
    let bwd = backward_matrices(profile, &ranks);
    let total = fwd.total;
    let mut post = vec![0.0f64; l];
    for k in 1..=l {
        let mut acc = NEG_INF;
        for i in 1..=n {
            let v = fwd.m[i][k] + bwd.m[i][k];
            if v != NEG_INF {
                acc = ln_add(acc, v);
            }
        }
        post[k - 1] = if acc == NEG_INF {
            0.0
        } else {
            (acc - total).exp()
        };
    }
    Ok(post)
}

pub(crate) fn posterior_emit_matrix(profile: &ProfileHmm, ranks: &[usize]) -> Vec<Vec<f64>> {
    let n = ranks.len();
    let l = profile.num_cols;
    let fwd = forward_matrices(profile, ranks);
    let bwd = backward_matrices(profile, ranks);
    let total = fwd.total;
    let mut post = vec![vec![0.0f64; l + 1]; n + 1];
    for (i, row) in post.iter_mut().enumerate().skip(1) {
        for (k, cell) in row.iter_mut().enumerate().skip(1) {
            let v = fwd.m[i][k] + bwd.m[i][k];
            *cell = if v == NEG_INF { 0.0 } else { (v - total).exp() };
        }
    }
    post
}

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;

    #[test]
    fn test_forward_equals_backward() {
        let profile = example_profile();
        for query in [
            &b"ACGTAC"[..],
            &b"ACAC"[..],
            &b"ACGTTTAC"[..],
            &b"AGTAC"[..],
            &b"AC"[..],
            &b""[..],
        ] {
            let forward = profile.forward(query).unwrap();
            let backward = profile.backward(query).unwrap();
            assert!(
                (forward - backward).abs() < 1e-9,
                "forward {} backward {} for {:?}",
                forward,
                backward,
                query
            );
        }
    }

    #[test]
    fn test_forward_is_at_least_the_best_path() {
        let profile = example_profile();
        for query in [&b"ACGTAC"[..], &b"ACAC"[..], &b"ACGTTTAC"[..]] {
            let forward = profile.forward(query).unwrap();
            let best = profile.best_alignment(query).unwrap().log_prob;
            assert!(forward >= best - 1e-9);
        }
    }

    #[test]
    fn test_posteriors_are_probabilities() {
        let profile = example_profile();
        for query in [&b"ACGTAC"[..], &b"ACAC"[..]] {
            let posterior = profile.posterior_match(query).unwrap();
            assert_eq!(posterior.len(), profile.num_columns());
            for value in &posterior {
                assert!(
                    *value >= -1e-12 && *value <= 1.0 + 1e-9,
                    "posterior {} out of [0,1]",
                    value
                );
            }
        }
    }

    #[test]
    fn test_deleted_position_has_lower_posterior() {
        let profile = example_profile();
        let posterior = profile.posterior_match(b"ACAC").unwrap();
        assert!(posterior[2] < posterior[0]);
        assert!(posterior[3] < posterior[0]);
    }
}
