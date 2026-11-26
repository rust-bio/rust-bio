// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::alignment::ProfileOp;
use crate::stats::hmm::profile::builder::GAP;
use crate::stats::hmm::profile::{ProfileError, ProfileHmm};
use std::f64::consts::LN_2;

/// The byte that pads insert slots in the rows returned by [`ProfileHmm::align_all`].
pub const PAD: u8 = b'.';

pub(crate) fn consensus(profile: &ProfileHmm) -> Vec<u8> {
    let mut out = Vec::with_capacity(profile.num_cols);
    for k in 1..=profile.num_cols {
        let mut best = 0usize;
        let mut best_lp = profile.emission[k][0];
        for a in 1..profile.symbols.len() {
            if profile.emission[k][a] > best_lp {
                best_lp = profile.emission[k][a];
                best = a;
            }
        }
        out.push(profile.symbols[best]);
    }
    out
}

pub(crate) fn position_information(profile: &ProfileHmm) -> Vec<f64> {
    let mut out = Vec::with_capacity(profile.num_cols);
    for k in 1..=profile.num_cols {
        let mut bits = 0.0;
        for a in 0..profile.symbols.len() {
            let p = profile.emission[k][a].exp();
            if p > 0.0 {
                bits += p * (profile.emission[k][a] - profile.null[a]) / LN_2;
            }
        }
        out.push(bits);
    }
    out
}

pub(crate) fn align_all(
    profile: &ProfileHmm,
    queries: &[&[u8]],
) -> Result<Vec<Vec<u8>>, ProfileError> {
    let l = profile.num_cols;
    let mut matched: Vec<Vec<Option<u8>>> = Vec::with_capacity(queries.len());
    let mut inserts: Vec<Vec<Vec<u8>>> = Vec::with_capacity(queries.len());

    for query in queries {
        let aln = profile.best_alignment(query)?;
        let mut row_match: Vec<Option<u8>> = vec![None; l + 1];
        let mut row_inserts: Vec<Vec<u8>> = vec![Vec::new(); l + 1];
        let mut pos = 0usize;
        for op in &aln.ops {
            match op {
                ProfileOp::Match(r) => {
                    pos += 1;
                    row_match[pos] = Some(*r);
                }
                ProfileOp::Delete => {
                    pos += 1;
                    row_match[pos] = None;
                }
                ProfileOp::Insert(r) => {
                    row_inserts[pos].push(*r);
                }
            }
        }
        matched.push(row_match);
        inserts.push(row_inserts);
    }

    let mut insert_width = vec![0usize; l + 1];
    for row in &inserts {
        for k in 0..=l {
            if row[k].len() > insert_width[k] {
                insert_width[k] = row[k].len();
            }
        }
    }

    let mut out = Vec::with_capacity(queries.len());
    for q in 0..queries.len() {
        let mut line = Vec::new();
        append_inserts(&mut line, &inserts[q][0], insert_width[0]);
        for k in 1..=l {
            line.push(match matched[q][k] {
                Some(r) => r,
                None => GAP,
            });
            append_inserts(&mut line, &inserts[q][k], insert_width[k]);
        }
        out.push(line);
    }
    Ok(out)
}

pub(crate) fn relative_entropy(profile: &ProfileHmm) -> f64 {
    position_information(profile).iter().sum()
}

pub(crate) fn to_alignment_text(
    profile: &ProfileHmm,
    queries: &[&[u8]],
) -> Result<String, ProfileError> {
    let cons = consensus(profile);
    let mut rows: Vec<&[u8]> = Vec::with_capacity(queries.len() + 1);
    rows.push(&cons);
    rows.extend_from_slice(queries);
    let block = align_all(profile, &rows)?;
    let mut text = String::new();
    for row in &block {
        for &byte in row {
            text.push(byte as char);
        }
        text.push('\n');
    }
    Ok(text)
}

fn append_inserts(line: &mut Vec<u8>, residues: &[u8], width: usize) {
    for &r in residues {
        line.push(r);
    }
    for _ in residues.len()..width {
        line.push(PAD);
    }
}

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;
    use crate::stats::hmm::profile::{GAP, PAD};

    #[test]
    fn test_consensus_and_information() {
        let profile = example_profile();
        assert_eq!(profile.consensus(), b"ACGTAC".to_vec());
        let information = profile.position_information();
        assert_eq!(information.len(), profile.num_columns());
        assert!(information.iter().all(|&bits| bits > 0.0));
        let sum: f64 = information.iter().sum();
        assert!((profile.relative_entropy() - sum).abs() < 1e-12);
    }

    #[test]
    fn test_align_all_is_rectangular() {
        let profile = example_profile();
        let queries: Vec<&[u8]> = vec![b"ACGTAC", b"ACAC", b"ACGTTTAC", b"AGGTAC"];
        let block = profile.align_all(&queries).unwrap();
        assert_eq!(block.len(), queries.len());
        let width = block[0].len();
        assert!(width >= profile.num_columns());
        for row in &block {
            assert_eq!(row.len(), width);
        }
    }

    #[test]
    fn test_align_all_marks_gaps_and_padding() {
        let profile = example_profile();
        let block = profile
            .align_all(&[&b"ACAC"[..], &b"ACGTTTAC"[..]])
            .unwrap();
        let flat: Vec<u8> = block.iter().flatten().copied().collect();
        assert!(flat.contains(&GAP));
        assert!(flat.contains(&PAD));
    }

    #[test]
    fn test_alignment_text_starts_with_consensus() {
        let profile = example_profile();
        let text = profile
            .to_alignment_text(&[&b"ACGTAC"[..], &b"ACAC"[..]])
            .unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], "ACGTAC");
        assert_eq!(lines[1], "ACGTAC");
        assert_eq!(lines[2], "AC--AC");
    }
}
