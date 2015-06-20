// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Rank/Select data structure based on Gonzalez, Grabowski, Mäkinen, Navarro (2005).
//! This implementation uses only a single level of blocks, and performs well for large n.
//!
//! Example
//!
//! ```
//! #![feature(bitvec)]
//! use bio::data_structures::rank_select::RankSelect;
//! use std::collections::BitVec;
//! let mut bits = BitVec::from_elem(64, false);
//! bits.set(5, true);
//! bits.set(32, true);
//! let rs = RankSelect::new(bits, 1);
//! assert!(rs.rank(1).unwrap() == 0);
//! assert!(rs.rank(5).unwrap() == 1);
//! assert!(rs.rank(6).unwrap() == 1);
//! assert!(rs.rank(32).unwrap() == 2);
//! assert!(rs.rank(33).unwrap() == 2);
//! assert!(rs.rank(64) == None);
//!
//! assert!(rs.select(0).unwrap() == 0);
//! assert!(rs.select(1).unwrap() == 5);
//! ```


use std::collections::BitVec;


pub struct RankSelect {
    n: usize,
    bits: Vec<u8>,
    superblocks: Vec<u32>,
    s: usize,
}


impl RankSelect {
    /// Create a new instance.
    ///
    /// # Arguments
    ///
    /// * `bits` - A bit vector.
    /// * `k` - Determines the size (k * 32 bits) of the superblocks.
    ///   A small k means faster rank query times at the expense of using more
    ///   space and slower select query times.
    ///   The data structure needs O(n + n log n / (k * 32)) bits with n being the bits of the given bitvector.
    ///   The data structure is succinct if k is chosen as a sublinear function of n
    ///   (e.g. k = (log n)² / 32).
    pub fn new(bits: BitVec, k: usize) -> RankSelect {
        let n = bits.len();
        let raw = bits.to_bytes();
        let s = k * 32;

        RankSelect { n: n, s: s, superblocks: superblocks(n, s, &raw), bits: raw }
    }

    /// Get the rank of a given bit, i.e. the number of 1-bits in the bitvector up to i (inclusive).
    /// Complexity: O(k).
    ///
    /// # Arguments
    ///
    /// * `i` - Position of the bit to determine the rank for.
    pub fn rank(&self, i: usize) -> Option<u32> {
        if i >= self.n {
            None
        }
        else {
            let s = i / self.s; // the superblock
            let b = i / 8; // the block
            // take the superblock rank
            let mut rank = self.superblocks[s];
            // add the rank within the block
            rank += (self.bits[b] >> 7 - i % 8).count_ones();
            // add the popcounts of blocks in between
            rank += self.bits[s * 32 / 8..b].iter()
                .map(|&a| a.count_ones())
                .fold(0, |a, b| a + b);

            Some(rank)
        }
    }

    /// Get the smallest bit with a given rank.
    /// Complexity: O(log (n / k) + k).
    ///
    /// # Arguments
    ///
    /// * `j` - The rank to find the smallest bit for.
    pub fn select(&self, j: u32) -> Option<usize> {
        let mut superblock = match self.superblocks.binary_search(&j) {
            Ok(i)  => i, // superblock with same rank exists
            Err(i) => i
        };
        if superblock > 0 {
            superblock -= 1;
        }
        let mut rank = self.superblocks[superblock];

        let first_block = superblock * self.s / 8;
        for (block, &b) in self.bits[first_block..].iter().enumerate() {
            let p = b.count_ones();
            if rank + p >= j {
                let mut bit = 0b10000000;
                for i in 0..8usize {
                    rank += (b & bit > 0) as u32;
                    if rank == j {
                        return Some((first_block + block) * 8 + i);
                    }
                    bit >>= 1;
                }
            }
            rank += p;
        }

        None
    }
}


fn superblocks(n: usize, s: usize, raw_bits: &Vec<u8>) -> Vec<u32> {
    let mut superblocks = Vec::with_capacity(n / s + 1);
    let mut rank: u32 = 0;
    let mut i = 0;
    for &b in raw_bits.iter() {
        if i % s == 0 {
            superblocks.push(rank);
        }
        rank += b.count_ones();
        i += 8;
    }

    superblocks
}
