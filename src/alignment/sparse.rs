// Copyright 2014-2015 Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Calculate 'sparse' alignments from kmer matches. Can be much faster than 
//! Smith-Waterman for long string, when a large enough k is used. 
//! Complexity: O(n * log(n)) for a pair of strings with n k-kmer matches. This
//! approach is useful for generating an approximate 'backbone' alignments 
//! between two long sequences, for example in long-read alignment or 
//! genome-genome alignment. The backbone alignment can be used as-is, or can serve
//! as a guide for a banded alignment.  By tuning k so that len(query) + len(reference) < 4^k,
//! the number of false positive kmer matches is kept small, resulting in very
//! fast run times for long strings.
//!
//! # Example
//!
//! ```
//! use bio::alignment::sparse::*;
//!
//! let s1 =   "ACGTACGATAGGTA";
//! let s2 = "TTACGTACGATAGGTATT";
//! let k = 8;
//! let matches = find_kmer_matches(&s1, &s2, k);
//! let sparse_al = lcskpp(&matches, k);
//! let match_path: Vec<(u32,u32)> = sparse_al.path.iter().map(|i| matches[*i]).collect();
//! assert_eq!(match_path, vec![(0,2), (1,3), (2,4), (3,5), (4,6), (5,7), (6,8)]);
//! assert_eq!(sparse_al.score, 14);

use std::cmp::max;
use data_structures::bit_tree::MaxBitTree;
use data_structures::qgram_index::QGramIndex;
use alphabets::Alphabet;

/// Result of a sparse alignment
pub struct SparseAlignmentResult {
    /// LCSk++ path, represented as vector of indices into the input matches vector.
    pub path: Vec<usize>,
    // Score of the path, which is the number of bases covered by the matched kmers.
    pub score: u32,
    // Full DP vector, which can generally be ignored. (It may be useful for testing purposes).
    pub dp_vector: Vec<(u32, i32)>,
}

/// Sparse DP routine for Longest Common Subsequence in length k substrings.  Also known of LCSk++
/// From LCSk++: Practical similarity metric for long strings. Filip Pavetić, Goran Žužić, Mile Šikić
/// Paper here :https://arxiv.org/abs/1407.2407.  Original implementation here: https://github.com/fpavetic/lcskpp
///
/// # Arguments
///
/// * `matches` - a vector of tuples indicating the (string1 position, string2 position) kmer matches between the strings  
/// * `k` - the kmer length used for matching
///
/// # Return value
///
/// The method returns a `SparseAlignmentResult` struct with the following fields:
/// * `path` is the LCSk++ path, represented as vector of indices into the input matches vector.
/// * `score` is the score of the path, which is the number of bases covered by the matched kmers.
/// * `dp_vector` is the full DP vector, which can generally be ignored. (It may be useful for testing purposes).
pub fn lcskpp(matches: &Vec<(u32, u32)>, k: usize) -> SparseAlignmentResult {

    let k = k as u32;

    // incoming matches must be sorted to let us find the predecssor kmers by binary search.
    for i in 1 .. matches.len() {
        assert!(matches[i-1] < matches[i]);
    }

    let mut events: Vec<(u32, u32, u32)> = Vec::new();
    let mut n = 0;

    for (idx, &(x,y)) in matches.iter().enumerate() {
        events.push((x,y,(idx+matches.len()) as u32));
        events.push((x+k,y+k,idx as u32));

        n = max(n, x+k);
        n = max(n, y+k);
    }
    events.sort(); 


    let mut max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(n as usize);
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (0,-1);

    for _ in 0 .. events.len() {
        dp.push((0,0));
    }

    for ev in events {
        let p = (ev.2 % matches.len() as u32) as usize;
        let j = ev.1;
        let is_start = ev.2 >= (matches.len() as u32);
        

        if is_start {
            dp[p] = (k, -1);
            let (best_value, best_position) = max_col_dp.get(j as usize);
            if best_value > 0 {
                dp[p] = (k + best_value, best_position as i32);
                best_dp = max(best_dp, (dp[p].0, p as i32));
            }
        } else {
            // See if this kmer continues a diffent kmer
            if ev.0 >= k+1 && ev.1 >= k+1 {
                match matches.binary_search(&(ev.0-k-1, ev.1-k-1)) {
                    Ok(cont_idx) => {
                        let prev_score = dp[cont_idx].0;
                        let candidate = (prev_score + 1, cont_idx as i32);
                        dp[p] = max(dp[p], candidate);
                        best_dp = max(best_dp, (dp[p].0, p as i32));
                    },
                    _ => (),
                }
            }

            max_col_dp.set(ev.1 as usize, (dp[p].0, p as u32));
        }
    }

    let mut traceback = Vec::new();
    let (best_score, mut prev_match) = best_dp;
    while prev_match >= 0 {
        traceback.push(prev_match as usize);
        prev_match = dp[prev_match as usize].1;
    }
    traceback.reverse();
    SparseAlignmentResult { path: traceback, score: best_score, dp_vector: dp }
}


/// Find all matches of length k between two strings, using a q-gram
/// index. For very long reference strings, it may be more efficient to use and
/// FMD index to generate the matches. Note that this method is mainly for 
/// demonstration & testing purposes.  For aligning many query sequences
/// against the same reference, you should reuse the QGramIndex of the reference.
pub fn find_kmer_matches<T: AsRef<[u8]>>(query: &T, reference: &T, k: usize) -> Vec<(u32, u32)> {

    let slc1 = query.as_ref();
    let slc2 = reference.as_ref();

    let mut alphabet = Alphabet::new(slc2);
    if !alphabet.is_word(slc1) {
        for i in slc1 {
            alphabet.insert(*i)
        }
    }

    let qgram_index = QGramIndex::new(k as u32, slc2, &alphabet);
    let mut matches = Vec::new();

    /// Generate all exact matches of length q
    for (i, qgram) in qgram_index.qgrams(slc1).enumerate() {
        for &p in qgram_index.qgram_matches(qgram) {
            matches.push((i as u32, p as u32));
        }
    }

    matches.sort();
    matches
}


#[cfg(test)]
mod sparse_alignment {
    use super::find_kmer_matches;

    #[test]
    pub fn test_find_kmer_matches() {
        let s1 = "ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
        let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
        let k = 8;
        //let s1 = "  ACGTACGATAGATCCGTACGTAACA     GTACAGTATATCAGTTATATGCGATA";
        //let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";

        let hits = find_kmer_matches(&s1, &s2, k);
        assert_eq!(hits.len(), (25-k+1) + (24-k+1));
        //println!("hits: {:?}", hits);
    }


    #[test]
    pub fn test_lcskpp0() {
        let s1 =   "ACGTACGATAGGTA";
        let s2 = "TTACGTACGATAGGTATT";
        let k = 8;
        let matches = super::find_kmer_matches(&s1, &s2, k);
        let res = super::lcskpp(&matches, k);
        let match_path: Vec<(u32,u32)> = res.path.iter().map(|i| matches[*i]).collect();
        assert_eq!(match_path, vec![(0,2), (1,3), (2,4), (3,5), (4,6), (5,7), (6,8)]);
        assert_eq!(res.score, 14);
     }


    #[test]
    pub fn test_lcskpp1() {
        let s1 = "ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
        let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
        let k = 8;

        let matches = super::find_kmer_matches(&s1, &s2, k);
        let res = super::lcskpp(&matches, k);
        
        // For debugging: 
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        // Should have 25bp group of matches plus a 24 bp group of matches
        assert_eq!(res.score, 25 + 24);
     }

    #[test]
    pub fn test_lcskpp2() {
        // Match the same string -- should get a diagonal traceback, despite lots of off-diagonal homology
        let s1 = "ACGTACGATAGATCCGACGTACGTACGTTCAGTTATATGACGTACGTACGTAACATTTTTGTA";
        let k = 5;

        let matches = super::find_kmer_matches(&s1, &s1, k);
        let res = super::lcskpp(&matches, k);

        // For debugging: 
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        assert_eq!(res.score, s1.len() as u32);

        for i in 0..res.path.len() {
            assert_eq!(matches[res.path[i] as usize], (i as u32, i as u32));
        }
     }
 }