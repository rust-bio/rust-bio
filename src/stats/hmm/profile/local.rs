// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::{ProfileError, ProfileHmm, NEG_INF};
use std::f64::consts::LN_2;

/// A region of a query that aligns locally to a range of positions of a
/// [`ProfileHmm`](super::ProfileHmm).
///
/// Coordinates are 1-based and inclusive. A hit without any aligned region has a score of zero
/// and all coordinates zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalHit {
    /// The log-odds of the region against the background model, in bits.
    pub bit_score: f64,
    /// The first aligned residue of the query.
    pub query_start: usize,
    /// The last aligned residue of the query.
    pub query_end: usize,
    /// The first aligned position of the model.
    pub first_position: usize,
    /// The last aligned position of the model.
    pub last_position: usize,
}

const FROM_START: u8 = 0;
const FROM_M: u8 = 1;
const FROM_I: u8 = 2;
const FROM_D: u8 = 3;

pub(crate) fn local_alignment(
    profile: &ProfileHmm,
    query: &[u8],
) -> Result<LocalHit, ProfileError> {
    let ranks = profile.ranks(query)?;
    let n = ranks.len();
    let l = profile.num_cols;

    let mut m = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut ins = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut del = vec![vec![NEG_INF; l + 1]; n + 1];
    let mut back = vec![vec![FROM_START; l + 1]; n + 1];

    let mut best = 0.0f64;
    let mut best_i = 0usize;
    let mut best_k = 0usize;

    for i in 1..=n {
        let residue = ranks[i - 1];
        for k in 1..=l {
            let lo = profile.emission[k][residue] - profile.null[residue];
            let mut prev = 0.0;
            let mut from = FROM_START;
            if m[i - 1][k - 1] > prev {
                prev = m[i - 1][k - 1];
                from = FROM_M;
            }
            if ins[i - 1][k - 1] > prev {
                prev = ins[i - 1][k - 1];
                from = FROM_I;
            }
            if del[i - 1][k - 1] > prev {
                prev = del[i - 1][k - 1];
                from = FROM_D;
            }
            m[i][k] = lo + prev;
            back[i][k] = from;

            let open = m[i - 1][k] + profile.t_mi[k];
            let extend = ins[i - 1][k] + profile.t_ii[k];
            ins[i][k] = open.max(extend);

            let dopen = m[i][k - 1] + profile.t_md[k - 1];
            let dextend = del[i][k - 1] + profile.t_dd[k - 1];
            del[i][k] = dopen.max(dextend);

            if m[i][k] > best {
                best = m[i][k];
                best_i = i;
                best_k = k;
            }
        }
    }

    if best <= 0.0 {
        return Ok(LocalHit {
            bit_score: 0.0,
            query_start: 0,
            query_end: 0,
            first_position: 0,
            last_position: 0,
        });
    }

    let query_end = best_i;
    let last_position = best_k;
    let mut i = best_i;
    let mut k = best_k;
    let mut state = FROM_M;
    loop {
        match state {
            FROM_M => {
                if i == 0 || k == 0 || back[i][k] == FROM_START {
                    break;
                }
                let from = back[i][k];
                i -= 1;
                k -= 1;
                state = from;
            }
            FROM_I => {
                if i == 0 {
                    break;
                }
                let open = m[i - 1][k] + profile.t_mi[k];
                let extend = ins[i - 1][k] + profile.t_ii[k];
                state = if extend > open { FROM_I } else { FROM_M };
                i -= 1;
            }
            _ => {
                if k == 0 {
                    break;
                }
                let dopen = m[i][k - 1] + profile.t_md[k - 1];
                let dextend = del[i][k - 1] + profile.t_dd[k - 1];
                state = if dextend > dopen { FROM_D } else { FROM_M };
                k -= 1;
            }
        }
    }

    Ok(LocalHit {
        bit_score: best / LN_2,
        query_start: i,
        query_end,
        first_position: k,
        last_position,
    })
}

pub(crate) fn scan(
    profile: &ProfileHmm,
    query: &[u8],
    min_bits: f64,
) -> Result<Vec<LocalHit>, ProfileError> {
    let mut hits = Vec::new();
    scan_range(profile, query, 0, min_bits, &mut hits)?;
    hits.sort_by_key(|hit| hit.query_start);
    Ok(hits)
}

fn scan_range(
    profile: &ProfileHmm,
    query: &[u8],
    offset: usize,
    min_bits: f64,
    hits: &mut Vec<LocalHit>,
) -> Result<(), ProfileError> {
    if query.is_empty() {
        return Ok(());
    }
    let hit = local_alignment(profile, query)?;
    if hit.bit_score < min_bits || hit.query_end == 0 {
        return Ok(());
    }
    let left_len = hit.query_start.saturating_sub(1);
    let right_start = hit.query_end;
    hits.push(LocalHit {
        bit_score: hit.bit_score,
        query_start: hit.query_start + offset,
        query_end: hit.query_end + offset,
        first_position: hit.first_position,
        last_position: hit.last_position,
    });
    if left_len > 0 {
        scan_range(profile, &query[..left_len], offset, min_bits, hits)?;
    }
    if right_start < query.len() {
        scan_range(
            profile,
            &query[right_start..],
            offset + right_start,
            min_bits,
            hits,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::alphabets::Alphabet;
    use crate::stats::hmm::profile::example_profile;
    use crate::stats::hmm::profile::ProfileHmm;

    #[test]
    fn test_local_score_is_never_negative() {
        let profile = example_profile();
        for query in [
            &b"TTTTACGTACTTTT"[..],
            &b"GGGGGG"[..],
            &b"ACGTAC"[..],
            &b""[..],
        ] {
            let hit = profile.local_alignment(query).unwrap();
            assert!(hit.bit_score >= 0.0);
            assert!(hit.query_end >= hit.query_start);
            assert!(hit.last_position >= hit.first_position);
        }
    }

    #[test]
    fn test_local_ends_are_free() {
        let profile = example_profile();
        let motif = profile.local_alignment(b"ACGTAC").unwrap();
        let flanked = profile.local_alignment(b"TTTTACGTACTTTT").unwrap();
        assert!(motif.bit_score > 0.0);
        assert!(flanked.bit_score >= motif.bit_score - 1e-9);
    }

    #[test]
    fn test_local_hit_coordinates_are_one_based() {
        let profile = example_profile();
        let hit = profile.local_alignment(b"ACGTAC").unwrap();
        assert_eq!(hit.query_start, 1);
        assert_eq!(hit.query_end, 6);
        assert_eq!(hit.first_position, 1);
        assert_eq!(hit.last_position, 6);
    }

    #[test]
    fn test_scan_hits_do_not_overlap() {
        let profile = example_profile();
        let hits = profile.scan(b"ACGTACTTTTTTTTACGTAC", 3.0).unwrap();
        assert!(!hits.is_empty());
        for window in hits.windows(2) {
            assert!(window[0].query_start < window[1].query_start);
            assert!(window[0].query_end < window[1].query_start);
        }
        for hit in &hits {
            assert!(hit.bit_score >= 3.0);
        }
    }

    #[test]
    fn test_scan_reports_single_residue_hits() {
        let mut builder = ProfileHmm::builder(&Alphabet::new(b"ACGT"));
        for _ in 0..4 {
            builder.add_row(b"AC");
        }
        let profile = builder.build().unwrap();
        let hits = profile.scan(b"CA", 0.5).unwrap();
        assert_eq!(hits.len(), 2);
        assert_eq!((hits[0].query_start, hits[0].query_end), (1, 1));
        assert_eq!((hits[0].first_position, hits[0].last_position), (2, 2));
        assert_eq!((hits[1].query_start, hits[1].query_end), (2, 2));
        assert_eq!((hits[1].first_position, hits[1].last_position), (1, 1));
    }

    #[test]
    fn test_scan_threshold_filters_hits() {
        let profile = example_profile();
        assert!(profile.scan(b"GGGGGG", 4.0).unwrap().is_empty());
        assert!(profile.scan(b"", 0.0).unwrap().is_empty());
    }
}
