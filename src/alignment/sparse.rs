// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Calculate 'sparse' alignments from kmer matches. Can be much faster than Smith-Waterman for long string, when a large enough k is used. 
//! Complexity: O(n * log(n)) for a pair of strings with n k-kmer matches
//!
//!
//! # Example
//!
//! ```
//! use bio::alignment::sparse::*;
//!
//! let s1 =   "ACGTACGATAGGTA";
//! let s2 = "TTACGTACGATAGGTATT";
//! let k = 8;
//! let matches = find_kmer_matches(s1, s2, k);
//! let (tb, score,  _) = lcskpp(&matches, k as u32);
//! let match_path: Vec<(u32,u32)> = tb.iter().map(|i| matches[*i]).collect();
//! assert_eq!(match_path, vec![(0,2), (1,3), (2,4), (3,5), (4,6), (5,7), (6,8)]);
//! assert_eq!(score, 14);

use std::cmp::max;
use std::collections::HashMap;

/// Sparse DP routine for Longest Common Subsequence in length k substrings. Also known of LCSk++
///  From LCSk++: Practical similarity metric for long strings. Filip Pavetić, Goran Žužić, Mile Šikić
/// Paper here :https://arxiv.org/abs/1407.2407.  Original implementation here: https://github.com/fpavetic/lcskpp
/// matches is a tuple of the positions of all k-mer matches between two string. Match vector must be sorted by caller.
/// The first element of the return tuple is the LCSk++ path, represented as vector of indices into the matches vector.
/// The second element of the return tuple is the score of the path, which is the number of bases covered by the matched kmers.
/// The third element of the return tuple is full DP vector, which can generally be ignore. (It may be useful for testing purposes).
pub fn lcskpp(matches: &Vec<(u32, u32)>, k: u32) -> (Vec<usize>, u32, Vec<(u32,i32)>) {

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

    (traceback, best_score, dp)
}

pub fn find_kmer_matches(str1: &str, str2: &str, k: usize) -> Vec<(u32, u32)> {

    let mut set: HashMap<&str, Vec<u32>> = HashMap::new();
    let mut matches = Vec::new();

    for i in 0 .. str1.len() - k + 1 {
        set.entry(&str1[i..i+k]).or_insert_with(|| Vec::new()).push(i as u32);
    }

    for i in 0 .. str2.len() - k + 1 {
        let slc = &str2[i..i+k];
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


struct MaxBitTree<T: Default + Ord> {
    tree: Vec<T>,
}

/// BIT-tree (Binary Indexed Trees, aka Fenwick Tree)
/// Maintains a prefix-sum or prefix-max that can be efficiently
/// Queried and updated.
/// Implementation outlined here:
/// https://www.topcoder.com/community/data-science/data-science-tutorials/binary-indexed-trees/
impl<T: Ord + Default + Copy> MaxBitTree<T> {
    pub fn new(len: usize) -> MaxBitTree<T> {
        let mut tree = Vec::new();

        // Pad length by one. The first element is unused.  
        // Done this way to make the tree structure work correctly.
        for _ in 0..(len+2) {
            tree.push(T::default());
        }

        MaxBitTree { tree: tree }
    }

    pub fn get(&self, idx: usize) -> T {
        let mut idx = idx + 1;
        let mut sum = T::default();
        while idx > 0 {
            sum = max(sum, self.tree[idx].clone());
            idx -= (idx as isize & -(idx as isize)) as usize;
        }

        sum
    }

    pub fn set(&mut self, idx: usize, val: T)
    {
        let mut idx = idx + 1;
        while idx < self.tree.len() {
            self.tree[idx] = max(self.tree[idx].clone(), val);
            idx += (idx as isize & -(idx as isize)) as usize;
        }
    }
}

#[cfg(test)]
mod sparse_alignment {
    use super::{MaxBitTree, find_kmer_matches};

    #[test]
    pub fn test_bit_tree() {
        let mut bit = MaxBitTree::new(10);

        bit.set(0, (1,0));
        bit.set(1, (1,1));
        bit.set(2, (2,2));
        bit.set(3, (3,3));
        bit.set(4, (2,4));
        bit.set(5, (2,5));
        bit.set(6, (4,6));
        bit.set(7, (5,7));

        assert_eq!(bit.get(0), (1, 0));
        assert_eq!(bit.get(1), (1, 1));
        assert_eq!(bit.get(2), (2, 2));
        assert_eq!(bit.get(3), (3, 3));
        assert_eq!(bit.get(4), (3, 3));
        assert_eq!(bit.get(5), (3, 3));
        assert_eq!(bit.get(6), (4, 6));
        assert_eq!(bit.get(7), (5, 7));
    }

    #[test]
    pub fn test_find_kmer_matches() {
        let s1 = "ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
        let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
        let k = 8;
        //let s1 = "  ACGTACGATAGATCCGTACGTAACA     GTACAGTATATCAGTTATATGCGATA";
        //let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";

        let hits = find_kmer_matches(s1, s2, k);
        assert_eq!(hits.len(), (25-k+1) + (24-k+1));
        //println!("hits: {:?}", hits);
    }


    #[test]
    pub fn test_lcskpp0() {
        let s1 =   "ACGTACGATAGGTA";
        let s2 = "TTACGTACGATAGGTATT";
        let k = 8;
        let matches = super::find_kmer_matches(s1, s2, k);
        let (tb, score,  _) = super::lcskpp(&matches, k as u32);
        let match_path: Vec<(u32,u32)> = tb.iter().map(|i| matches[*i]).collect();
        assert_eq!(match_path, vec![(0,2), (1,3), (2,4), (3,5), (4,6), (5,7), (6,8)]);
        assert_eq!(score, 14);
     }


    #[test]
    pub fn test_lcskpp1() {
        let s1 = "ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
        let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
        let k = 8;

        let matches = super::find_kmer_matches(s1, s2, k);
        let (_, score, _) = super::lcskpp(&matches, k as u32);
        
        // For debugging: 
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        // Should have 25bp group of matches plus a 24 bp group of matches
        assert_eq!(score, 25 + 24);
     }

    #[test]
    pub fn test_lcskpp2() {
        // Match the same string -- should get a diagonal traceback, despite lots of off-diagonal homology
        let s1 = "ACGTACGATAGATCCGACGTACGTACGTTCAGTTATATGACGTACGTACGTAACATTTTTGTA";
        let k = 5;

        let matches = super::find_kmer_matches(s1, s1, k);
        let (tb, score, _) = super::lcskpp(&matches, k as u32);

        // For debugging: 
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        assert_eq!(score, s1.len() as u32);

        for i in 0..tb.len() {
            assert_eq!(matches[tb[i] as usize], (i as u32, i as u32));
        }
     }
 }