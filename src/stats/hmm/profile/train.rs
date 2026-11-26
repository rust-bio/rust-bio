// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::alignment::ProfileOp;
use crate::stats::hmm::profile::forward_backward::{backward_matrices, forward_matrices};
use crate::stats::hmm::profile::{ProfileError, ProfileHmm, NEG_INF};

struct Expected {
    emission: Vec<Vec<f64>>,
    mm: Vec<f64>,
    mi: Vec<f64>,
    md: Vec<f64>,
    im: Vec<f64>,
    ii: Vec<f64>,
    id: Vec<f64>,
    dm: Vec<f64>,
    di: Vec<f64>,
    dd: Vec<f64>,
}

impl Expected {
    fn zeros(l: usize, alphabet_size: usize) -> Expected {
        Expected {
            emission: vec![vec![0.0; alphabet_size]; l + 1],
            mm: vec![0.0; l + 1],
            mi: vec![0.0; l + 1],
            md: vec![0.0; l + 1],
            im: vec![0.0; l + 1],
            ii: vec![0.0; l + 1],
            id: vec![0.0; l + 1],
            dm: vec![0.0; l + 1],
            di: vec![0.0; l + 1],
            dd: vec![0.0; l + 1],
        }
    }
}

pub(crate) fn total_log_likelihood(
    profile: &ProfileHmm,
    sequences: &[&[u8]],
) -> Result<f64, ProfileError> {
    let mut total = 0.0;
    for seq in sequences {
        let ranks = profile.ranks(seq)?;
        total += forward_matrices(profile, &ranks).total;
    }
    Ok(total)
}

pub(crate) fn reestimate(
    profile: &ProfileHmm,
    sequences: &[&[u8]],
    iterations: usize,
) -> Result<(ProfileHmm, Vec<f64>), ProfileError> {
    let mut current = profile.clone();
    let mut trace = vec![total_log_likelihood(&current, sequences)?];
    for _ in 0..iterations {
        let expected = accumulate(&current, sequences)?;
        current = maximize(&current, &expected);
        trace.push(total_log_likelihood(&current, sequences)?);
    }
    Ok((current, trace))
}

pub(crate) fn reestimate_viterbi(
    profile: &ProfileHmm,
    sequences: &[&[u8]],
    iterations: usize,
) -> Result<(ProfileHmm, Vec<f64>), ProfileError> {
    let mut current = profile.clone();
    let mut trace = vec![viterbi_objective(&current, sequences)?];
    for _ in 0..iterations {
        let counts = viterbi_counts(&current, sequences)?;
        current = maximize(&current, &counts);
        trace.push(viterbi_objective(&current, sequences)?);
    }
    Ok((current, trace))
}

fn viterbi_objective(profile: &ProfileHmm, sequences: &[&[u8]]) -> Result<f64, ProfileError> {
    let mut total = 0.0;
    for seq in sequences {
        total += profile.best_alignment(seq)?.log_prob;
    }
    Ok(total)
}

fn viterbi_counts(profile: &ProfileHmm, sequences: &[&[u8]]) -> Result<Expected, ProfileError> {
    let l = profile.num_cols;
    let mut acc = Expected::zeros(l, profile.symbols.len());
    let pc = 1.0;
    for k in 0..=l {
        acc.mm[k] = pc;
        acc.mi[k] = pc;
        acc.im[k] = pc;
        acc.ii[k] = pc;
        acc.dm[k] = pc;
        acc.di[k] = pc;
        if k < l {
            acc.md[k] = pc;
            acc.id[k] = pc;
            acc.dd[k] = pc;
        }
        if k >= 1 {
            for a in 0..profile.symbols.len() {
                acc.emission[k][a] = pc;
            }
        }
    }
    for seq in sequences {
        let aln = profile.best_alignment(seq)?;
        let mut node = 0usize;
        let mut prev = b'M';
        for op in &aln.ops {
            match op {
                ProfileOp::Match(sym) => {
                    match prev {
                        b'M' => acc.mm[node] += 1.0,
                        b'I' => acc.im[node] += 1.0,
                        _ => acc.dm[node] += 1.0,
                    }
                    node += 1;
                    acc.emission[node][profile.rank(*sym)?] += 1.0;
                    prev = b'M';
                }
                ProfileOp::Insert(_) => {
                    match prev {
                        b'M' => acc.mi[node] += 1.0,
                        b'I' => acc.ii[node] += 1.0,
                        _ => acc.di[node] += 1.0,
                    }
                    prev = b'I';
                }
                ProfileOp::Delete => {
                    match prev {
                        b'M' => acc.md[node] += 1.0,
                        b'I' => acc.id[node] += 1.0,
                        _ => acc.dd[node] += 1.0,
                    }
                    node += 1;
                    prev = b'D';
                }
            }
        }
        match prev {
            b'M' => acc.mm[l] += 1.0,
            b'I' => acc.im[l] += 1.0,
            _ => acc.dm[l] += 1.0,
        }
    }
    Ok(acc)
}

fn accumulate(profile: &ProfileHmm, sequences: &[&[u8]]) -> Result<Expected, ProfileError> {
    let l = profile.num_cols;
    let alphabet_size = profile.symbols.len();
    let mut acc = Expected::zeros(l, alphabet_size);

    for seq in sequences {
        let ranks = profile.ranks(seq)?;
        let n = ranks.len();
        let fwd = forward_matrices(profile, &ranks);
        let bwd = backward_matrices(profile, &ranks);
        let z = fwd.total;
        if z == NEG_INF {
            continue;
        }

        for i in 1..=n {
            let residue = ranks[i - 1];
            for k in 1..=l {
                let g = (fwd.m[i][k] + bwd.m[i][k] - z).exp();
                acc.emission[k][residue] += g;
            }
        }

        for k in 0..=l {
            for (i, &next) in ranks.iter().enumerate() {
                if k < l {
                    let to_m = profile.emission[k + 1][next] + bwd.m[i + 1][k + 1];
                    acc.mm[k] += (fwd.m[i][k] + profile.t_mm[k] + to_m - z).exp();
                    acc.im[k] += (fwd.i[i][k] + profile.t_im[k] + to_m - z).exp();
                    acc.dm[k] += (fwd.d[i][k] + profile.t_dm[k] + to_m - z).exp();
                }
                let to_i = profile.null[next] + bwd.i[i + 1][k];
                acc.mi[k] += (fwd.m[i][k] + profile.t_mi[k] + to_i - z).exp();
                acc.ii[k] += (fwd.i[i][k] + profile.t_ii[k] + to_i - z).exp();
                acc.di[k] += (fwd.d[i][k] + profile.t_di[k] + to_i - z).exp();
            }
            if k < l {
                for i in 0..=n {
                    let to_d = bwd.d[i][k + 1];
                    acc.md[k] += (fwd.m[i][k] + profile.t_md[k] + to_d - z).exp();
                    acc.id[k] += (fwd.i[i][k] + profile.t_id[k] + to_d - z).exp();
                    acc.dd[k] += (fwd.d[i][k] + profile.t_dd[k] + to_d - z).exp();
                }
            }
        }
        acc.mm[l] += (fwd.m[n][l] + profile.t_mm[l] - z).exp();
        acc.im[l] += (fwd.i[n][l] + profile.t_im[l] - z).exp();
        acc.dm[l] += (fwd.d[n][l] + profile.t_dm[l] - z).exp();
    }
    Ok(acc)
}

fn maximize(profile: &ProfileHmm, acc: &Expected) -> ProfileHmm {
    let l = profile.num_cols;
    let mut next = profile.clone();

    for k in 1..=l {
        let total: f64 = acc.emission[k].iter().sum();
        if total > 0.0 {
            next.emission[k] = acc.emission[k]
                .iter()
                .map(|&c| if c > 0.0 { (c / total).ln() } else { NEG_INF })
                .collect();
        }
    }

    for k in 0..=l {
        let last = k == l;
        let m_total = if last {
            acc.mm[k] + acc.mi[k]
        } else {
            acc.mm[k] + acc.mi[k] + acc.md[k]
        };
        if m_total > 0.0 {
            next.t_mm[k] = ln_share(acc.mm[k], m_total);
            next.t_mi[k] = ln_share(acc.mi[k], m_total);
            if !last {
                next.t_md[k] = ln_share(acc.md[k], m_total);
            }
        }
        let i_total = if last {
            acc.im[k] + acc.ii[k]
        } else {
            acc.im[k] + acc.ii[k] + acc.id[k]
        };
        if i_total > 0.0 {
            next.t_im[k] = ln_share(acc.im[k], i_total);
            next.t_ii[k] = ln_share(acc.ii[k], i_total);
            if !last {
                next.t_id[k] = ln_share(acc.id[k], i_total);
            }
        }
        if k >= 1 {
            let d_total = if last {
                acc.dm[k] + acc.di[k]
            } else {
                acc.dm[k] + acc.di[k] + acc.dd[k]
            };
            if d_total > 0.0 {
                next.t_dm[k] = ln_share(acc.dm[k], d_total);
                next.t_di[k] = ln_share(acc.di[k], d_total);
                if !last {
                    next.t_dd[k] = ln_share(acc.dd[k], d_total);
                }
            }
        }
    }
    next
}

fn ln_share(count: f64, total: f64) -> f64 {
    if count > 0.0 {
        (count / total).ln()
    } else {
        NEG_INF
    }
}

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;

    #[test]
    fn test_reestimate_never_decreases_likelihood() {
        let profile = example_profile();
        let sequences: Vec<&[u8]> = vec![
            b"ACGTAC", b"ACGTTAC", b"ACAC", b"AGGTAC", b"ACGGTAC", b"ACGTAG",
        ];
        let (trained, trace) = profile.reestimate(&sequences, 10).unwrap();
        assert_eq!(trace.len(), 11);
        for window in trace.windows(2) {
            assert!(
                window[1] >= window[0] - 1e-9,
                "likelihood dropped {} -> {}",
                window[0],
                window[1]
            );
        }
        assert_eq!(trace[0], profile.log_likelihood(&sequences).unwrap());
        assert!((trace[10] - trained.log_likelihood(&sequences).unwrap()).abs() < 1e-9);
    }

    #[test]
    fn test_reestimate_viterbi_never_decreases_objective() {
        let profile = example_profile();
        let sequences: Vec<&[u8]> = vec![b"ACGTAC", b"ACGTTAC", b"ACAC", b"AGGTAC", b"ACGGTAC"];
        let (_, trace) = profile.reestimate_viterbi(&sequences, 8).unwrap();
        assert_eq!(trace.len(), 9);
        for window in trace.windows(2) {
            assert!(
                window[1] >= window[0] - 1e-9,
                "viterbi training dropped {} -> {}",
                window[0],
                window[1]
            );
        }
    }

    #[test]
    fn test_zero_iterations_keep_the_model() {
        let profile = example_profile();
        let sequences: Vec<&[u8]> = vec![b"ACGTAC", b"ACAC"];
        let (same, trace) = profile.reestimate(&sequences, 0).unwrap();
        assert_eq!(trace.len(), 1);
        assert_eq!(same.to_text(), profile.to_text());
    }
}
