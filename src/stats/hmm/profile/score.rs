// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::forward_backward::forward_matrices;
use crate::stats::hmm::profile::{ProfileError, ProfileHmm};
use std::cmp::Ordering;
use std::f64::consts::LN_2;

/// A query that cleared the threshold of [`ProfileHmm::search`](super::ProfileHmm::search).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SearchHit {
    /// The index of the query in the slice passed to the search.
    pub index: usize,
    /// The bit score of the query.
    pub bit_score: f64,
}

pub(crate) fn bit_score(profile: &ProfileHmm, query: &[u8]) -> Result<f64, ProfileError> {
    let ranks = profile.ranks(query)?;
    let total = forward_matrices(profile, &ranks).total;
    let null_lp = profile.null_log_prob_ranks(&ranks);
    Ok((total - null_lp) / LN_2)
}

pub(crate) fn search(
    profile: &ProfileHmm,
    queries: &[&[u8]],
    threshold: f64,
) -> Result<Vec<SearchHit>, ProfileError> {
    let mut hits: Vec<SearchHit> = Vec::new();
    for (index, query) in queries.iter().enumerate() {
        let score = bit_score(profile, query)?;
        if score >= threshold {
            hits.push(SearchHit {
                index,
                bit_score: score,
            });
        }
    }
    hits.sort_by(|a, b| match b.bit_score.partial_cmp(&a.bit_score) {
        Some(Ordering::Equal) | None => a.index.cmp(&b.index),
        Some(order) => order,
    });
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;
    use std::f64::consts::LN_2;

    #[test]
    fn test_bit_score_is_null_corrected() {
        let profile = example_profile();
        for query in [&b"ACGTAC"[..], &b"ACAC"[..], &b"TTTTTT"[..]] {
            let expected =
                (profile.forward(query).unwrap() - profile.null_log_prob(query).unwrap()) / LN_2;
            assert!((profile.bit_score(query).unwrap() - expected).abs() < 1e-9);
        }
    }

    #[test]
    fn test_search_orders_hits_and_breaks_ties_by_index() {
        let profile = example_profile();
        let queries: Vec<&[u8]> = vec![b"ACGTAC", b"TTTTTT", b"ACGTAC", b"ACAC"];
        let hits = profile.search(&queries, f64::NEG_INFINITY).unwrap();
        assert_eq!(hits.len(), queries.len());
        assert_eq!(profile.search(&queries, f64::NEG_INFINITY).unwrap(), hits);
        for window in hits.windows(2) {
            let ordered = window[0].bit_score > window[1].bit_score
                || (window[0].bit_score == window[1].bit_score
                    && window[0].index < window[1].index);
            assert!(ordered);
        }
        let first = hits.iter().position(|hit| hit.index == 0).unwrap();
        let second = hits.iter().position(|hit| hit.index == 2).unwrap();
        assert!(first < second);
    }

    #[test]
    fn test_search_applies_threshold() {
        let profile = example_profile();
        let queries: Vec<&[u8]> = vec![b"ACGTAC", b"TTTTTT", b"ACGTAC", b"ACAC"];
        let hits = profile.search(&queries, 1.0).unwrap();
        assert!(hits.iter().all(|hit| hit.bit_score >= 1.0));
        for (index, query) in queries.iter().enumerate() {
            let present = hits.iter().any(|hit| hit.index == index);
            let qualifies = profile.bit_score(query).unwrap() >= 1.0;
            assert_eq!(present, qualifies);
        }
    }
}
