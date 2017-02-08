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
use std::collections::HashMap;

/// Result of a sparse alignment
#[derive(Debug, PartialEq, Eq)]
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


#[derive(PartialEq, Eq, Ord, PartialOrd, Default, Copy, Clone)]
struct PrevPtr {
    plane: u32,
    score: u32,
    d: u32,
    id: usize
}

impl PrevPtr {
    pub fn new(score: u32, x: u32, y: u32, id: usize, gap_extend: u32) -> PrevPtr {
        let d = x + y;
        PrevPtr {
            plane: score + (d * gap_extend),
            score: score,
            d: d,
            id: id,
        }
    }
}




/// Sparse DP routine generalizing LCSk++ method above penalize alignment gaps. 
/// A gap is an unknown combination of mismatch, insertion and deletions, and incurs
/// a penalty of gap_open + d * gap_extend, where d is the distance along the diagonal of the gap.
/// # Arguments
///
/// * `matches` - a vector of tuples indicating the (string1 position, string2 position) kmer matches between the strings  
/// * `k` - the kmer length used for matching
/// * `match_score` - reward for each matched base
/// * `gap_open` - score of opening a gap, including a mismatch gap. Must be negative.
/// * `gap_extend` - score for extending a gap along the diagonal. Must be negative.
///
/// # Return value
///
/// The method returns a `SparseAlignmentResult` struct with the following fields:
/// * `path` is the SDP path, represented as vector of indices into the input matches vector.
/// * `score` is the score of the path, which is the number of bases covered by the matched kmers.
/// * `dp_vector` is the full DP vector, which can generally be ignored. (It may be useful for testing purposes).
pub fn sdpkpp(matches: &Vec<(u32, u32)>, k: usize, match_score: u32, gap_open: i32, gap_extend: i32) -> SparseAlignmentResult {

    let k = k as u32;
    if gap_open > 0 || gap_extend > 0 {
        panic!("gap parameters cannot be positive")
    }
    let _gap_open = (-gap_open) as u32;
    let _gap_extend = (-gap_extend) as u32;

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


    let mut max_col_dp: MaxBitTree<PrevPtr> = MaxBitTree::new(n as usize);
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
            // Default case -- chain starts at this node
            dp[p] = (k * match_score, -1);

            // Find best previous chain, and extend. 
            let best_prev = max_col_dp.get(j as usize);
            if best_prev.score > 0 {
                
                let prev_d = best_prev.d;
                let cur_d = ev.0 + ev.1;
                let gap = _gap_open + (cur_d - prev_d) * _gap_extend;

                let reward = k * match_score;
                let new_score = (best_prev.score + reward).saturating_sub(gap);

                dp[p] = max(dp[p], (new_score, best_prev.id as i32));
                best_dp = max(best_dp, (dp[p].0, p as i32));
            }
        } else {
            // See if this kmer continues a diffent kmer
            if ev.0 >= k+1 && ev.1 >= k+1 {
                match matches.binary_search(&(ev.0-k-1, ev.1-k-1)) {
                    Ok(cont_idx) => {
                        let prev_score = dp[cont_idx].0;
                        let candidate = (prev_score + match_score, cont_idx as i32);
                        dp[p] = max(dp[p], candidate);
                        best_dp = max(best_dp, (dp[p].0, p as i32));
                    },
                    _ => (),
                }
            }

            let prev_frag = PrevPtr::new(dp[p].0, ev.0, ev.1, p, _gap_extend);
            max_col_dp.set(ev.1 as usize, prev_frag);
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
pub fn find_kmer_matches<T: AsRef<[u8]>>(seq1: &T, seq2: &T, k: usize) -> Vec<(u32, u32)> {

    let slc1 = seq1.as_ref();
    let slc2 = seq2.as_ref();

    let mut set: HashMap<&[u8], Vec<u32>> = HashMap::new();
    let mut matches = Vec::new();

    for i in 0 .. slc1.len() - k + 1 {
        set.entry(&slc1[i..i+k]).or_insert_with(|| Vec::new()).push(i as u32);
    }

    for i in 0 .. slc2.len() - k + 1 {
        let slc = &slc2[i..i+k];
        match set.get(slc) {
            Some(matches1) => {
                for pos1 in matches1 {
                    matches.push((*pos1, i as u32));
                }
            },
            None => (),
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

     pub fn strict_compare_lcskpp_sdpkpp(s1: &str, s2: &str) {
        let k = 8;
        let matches = super::find_kmer_matches(&s1, &s2, k);
        let res1 = super::lcskpp(&matches, k);
        let res2 = super::sdpkpp(&matches, k, 1, 0, 0);

        assert_eq!(res1, res2);
     }

    #[test]
    pub fn test_sdp() {
        let s1 =   "ACGTACGATAGGTA";
        let s2 = "TTACGTACGATAGGTATT";
        strict_compare_lcskpp_sdpkpp(s1, s2);
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
     pub fn test_sdp1() {
         let s1 = "ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
         let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
         strict_compare_lcskpp_sdpkpp(s1, s2);
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

     #[test]
     pub fn test_sdp2() {
         let s1 = "ACGTACGATAGATCCGACGTACGTACGTTCAGTTATATGACGTACGTACGTAACATTTTTGTA";
         strict_compare_lcskpp_sdpkpp(s1, s1);
     }


     // Test case from local SV caller alignments.
     // The query sequence ends in 1-2 copies of tandem repeat element
     // The target sequence end in >2 copies of the element
     // Without a gap penalty (i.e. LCSk++), the alignment of the
     // TRs is arbitrary, and way the implementation breaks ties may introduce
     // a gap while maintaining the same score.
     // The SDP code with gap open & extend penalties should resolve this.
     const QUERY_REPEAT: &'static str = "CCTCCCATCTCCACCCACCCTATCCAACCCTGGGGTGGCAGGTCATGAGTGACAGCCCCAAGGACACCAAGGGATGAAGCTTCTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGCAGAACAGGAGGCCCCATCTCCCGTTCTTACTCAGAAGCTGTCAGCAGGGCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTCCCCAGCTCTTGCTTCAGGGCCTCACGTTTCGCCCCCTGAGGCCTGGGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTG";

    const TARGET_REPEAT: &'static str = "CCTCCCATCTCCACCCACCCTATCCAACCCTGGGGTGGCAGGTCATGAGTGACAGCCCCAAGGACACCAAGGGATGAAGCTTCTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGCAGAACAGGAGGCCCCATCTCCCGTTCTTACTCAGAAGCTGTCAGCAGGGCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTCCCCAGCTCTTGCTTCAGGGCCTCACGTTTCGCCCCCTGAGGCCTGGGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGCACGGCTCCCAACTCTCTTCCGGCCAAGGATCCCGTGTTCCTGAAATGTCTTTCTACCAAACACAGTTGCTGTGTAACCACTCATTTCATTTTCCTAATTTGTGTTGATCCAGGACACGGGAGGAGACCTGGGCAGCGGCGGACTCATTGCAGGTCGCTCTGCGGTGAGGACGCCACAGGCAC";

    #[test]
    fn test_lcskpp_tandem_repeat() {
        let k = 8;
        let matches = super::find_kmer_matches(&QUERY_REPEAT, &TARGET_REPEAT, k);
        let res = super::lcskpp(&matches, k);

        // For debugging: 
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        assert_eq!(res.score, QUERY_REPEAT.len() as u32);

        // NOTE -- this test will fail, because LCSk++ introduces a gap in the placement of the TR
        // Corrected with gap scoring in SDP
        /*
        for i in 0..res.path.len() {
            assert_eq!(matches[res.path[i] as usize], (i as u32, i as u32));
        }
        */  
    }


    #[test]
    fn test_sdpkpp_tandem_repeat() {
        let k = 8;
        let matches = super::find_kmer_matches(&QUERY_REPEAT, &TARGET_REPEAT, k);
        let res = super::sdpkpp(&matches, k, 1, -1, -1);

        // For debugging: 
        /*
        for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
            println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        }
        println!("tb: {:?}", tb);
        */ 

        assert_eq!(res.score, QUERY_REPEAT.len() as u32);

        for i in 0..res.path.len() {
            assert_eq!(matches[res.path[i] as usize], (i as u32, i as u32));
        }
    }
 }