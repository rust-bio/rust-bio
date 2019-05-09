// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Rank/Select data structure based on Gonzalez, Grabowski, Mäkinen, Navarro (2005).
//! This implementation uses only a single level of blocks, and performs well for large n.
//!
//! Example
//!
//! ```
//! extern crate bv;
//! # extern crate bio;
//! # fn main() {
//! use bio::data_structures::rank_select::RankSelect;
//! use bv::BitVec;
//! use bv::BitsMut;
//!
//! let mut bits: BitVec<u8> = BitVec::new_fill(false, 64);
//! bits.set_bit(5, true);
//! bits.set_bit(32, true);
//! let rs = RankSelect::new(bits, 1);
//! assert!(rs.rank(6).unwrap() == 1);
//! # }
//! ```

use std::cmp;

use bv::BitVec;
use bv::Bits;

/// A rank/select data structure.
#[derive(Serialize, Deserialize)]
pub struct RankSelect {
    n: usize,
    bits: BitVec<u8>,
    superblocks_1: Vec<u64>,
    superblocks_0: Vec<u64>,
    s: usize,
    k: usize,
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
    pub fn new(bits: BitVec<u8>, k: usize) -> RankSelect {
        let n = bits.len() as usize;
        let s = k * 32;

        RankSelect {
            n,
            s,
            k,
            superblocks_1: superblocks(true, n, s, &bits),
            superblocks_0: superblocks(false, n, s, &bits),
            bits,
        }
    }

    /// Return the used k (see `RankSelect::new()`).
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get internal representation of bit vector.
    pub fn bits(&self) -> &BitVec<u8> {
        &self.bits
    }

    /// Return i-th bit.
    pub fn get(&self, i: u64) -> bool {
        self.bits.get_bit(i)
    }

    /// Get the 1-rank of a given bit, i.e. the number of 1-bits in the bitvector up to i (inclusive).
    /// Complexity: O(k).
    ///
    /// # Arguments
    ///
    /// * `i` - Position of the bit to determine the rank for.
    pub fn rank_1(&self, i: u64) -> Option<u64> {
        if i >= self.n as u64 {
            None
        } else {
            let s = i / self.s as u64; // the superblock
            let b = i / 8; // the block
            let j = i % 8; // the bit in the block
                           // take the superblock rank
            let mut rank = self.superblocks_1[s as usize];
            // add the rank within the block
            let mask = ((2u16 << j) - 1) as u8;
            rank += (self.bits.get_block(b as usize) & mask).count_ones() as u64;
            // add the popcounts of blocks from the beginning of the current superblock
            // up to the current block
            for block in (s * self.s as u64 / 8)..b {
                let b = self.bits.get_block(block as usize);
                rank += b.count_ones() as u64;
            }

            Some(rank)
        }
    }

    /// Get the 0-rank of a given bit, i.e. the number of 0-bits in the bitvector up to i (inclusive).
    /// Complexity: O(k).
    ///
    /// # Arguments
    ///
    /// * `i` - Position of the bit to determine the rank for.
    pub fn rank_0(&self, i: u64) -> Option<u64> {
        self.rank_1(i).map(|r| (i + 1) - r)
    }

    /// Alias for `RankSelect::rank_1`.
    pub fn rank(&self, i: u64) -> Option<u64> {
        self.rank_1(i)
    }

    /// Get the smallest bit with a given 1-rank.
    /// Complexity: O(log (n / k) + k).
    ///
    /// # Arguments
    ///
    /// * `j` - The rank to find the smallest bit for.
    pub fn select_1(&self, j: u64) -> Option<u64> {
        self.select_x(
            j,
            &self.superblocks_1,
            |bit| bit != 0,
            |block| block.count_ones(),
        )
    }

    /// Get the smallest bit with a given 0-rank.
    /// Complexity: O(log (n / k) + k).
    ///
    /// # Arguments
    ///
    /// * `j` - The rank to find the smallest bit for.
    pub fn select_0(&self, j: u64) -> Option<u64> {
        self.select_x(
            j,
            &self.superblocks_0,
            |bit| bit == 0,
            |block| block.count_zeros(),
        )
    }

    fn select_x<F: Fn(u8) -> bool, C: Fn(u8) -> u32>(
        &self,
        j: u64,
        superblocks: &[u64],
        is_match: F,
        count_all: C,
    ) -> Option<u64> {
        let mut superblock = match superblocks.binary_search(&j) {
            Ok(i) | Err(i) => i, // superblock with same rank exists
        };
        if superblock > 0 {
            superblock -= 1;
        }
        let mut rank = superblocks[superblock];

        let first_block = superblock * self.s / 8;
        for block in first_block..cmp::min(first_block + self.s / 8, self.bits.block_len()) {
            let b = self.bits.get_block(block);
            let p = count_all(b) as u64;
            if rank + p >= j {
                let mut bit = 0b1;
                // do not look at unused bits of the last block
                let max_bit = cmp::min(8, self.bits.len() - block as u64 * 8);
                for i in 0..max_bit {
                    rank += is_match(b & bit) as u64;
                    if rank == j {
                        return Some((first_block + block) as u64 * 8 + i);
                    }
                    bit <<= 1;
                }
            }
            rank += p;
        }

        None
    }

    /// Alias for `RankSelect::select_1`.
    pub fn select(&self, j: u64) -> Option<u64> {
        self.select_1(j)
    }
}

/// Create `n` superblocks of size `s` from a given bitvector.
fn superblocks(t: bool, n: usize, s: usize, bits: &BitVec<u8>) -> Vec<u64> {
    let mut superblocks = Vec::with_capacity(n / s + 1);
    let mut rank: u64 = 0;
    let mut i = 0;
    let nblocks = (bits.len() as f64 / 8.0).ceil() as usize;
    for block in 0..nblocks {
        let b = bits.get_block(block);
        if i % s == 0 {
            superblocks.push(rank);
        }
        rank += if t {
            b.count_ones() as u64
        } else {
            b.count_zeros() as u64
        };
        i += 8;
    }

    superblocks
}

#[cfg(test)]
mod tests {
    use super::*;
    use bv::bit_vec;
    use bv::BitVec;
    use bv::BitsMut;

    #[test]
    fn test_rank_select() {
        let mut bits: BitVec<u8> = BitVec::new_fill(false, 64);
        bits.set_bit(5, true);
        bits.set_bit(32, true);
        let rs = RankSelect::new(bits, 1);
        assert_eq!(rs.rank_1(1).unwrap(), 0);
        assert_eq!(rs.rank_1(5).unwrap(), 1);
        assert_eq!(rs.rank_1(6).unwrap(), 1);
        assert_eq!(rs.rank_1(7).unwrap(), 1);
        assert_eq!(rs.rank_1(32).unwrap(), 2);
        assert_eq!(rs.rank_1(33).unwrap(), 2);
        assert_eq!(rs.rank_1(64), None);
        assert_eq!(rs.select_1(0).unwrap(), 0);
        assert_eq!(rs.select_1(1).unwrap(), 5);
        assert_eq!(rs.rank_0(1).unwrap(), 2);
        assert_eq!(rs.rank_0(4).unwrap(), 5);
        assert_eq!(rs.rank_0(5).unwrap(), 5);
        assert_eq!(rs.select_0(0), None);
        assert_eq!(rs.select_0(1).unwrap(), 0);
        assert_eq!(rs.get(5), true);
        assert_eq!(rs.get(1), false);
        assert_eq!(rs.get(32), true);
    }

    #[test]
    fn test_select() {
        let bits: BitVec<u8> = bit_vec![true, false];
        let rs = RankSelect::new(bits, 1);

        assert_eq!(rs.select_0(0), Some(0));
        assert_eq!(rs.select_1(0), None);

        assert_eq!(rs.select_0(1), Some(1));
        assert_eq!(rs.select_1(1), Some(0));

        assert_eq!(rs.select_0(2), None);
        assert_eq!(rs.select_1(2), None);
    }

    #[test]
    fn test_single_select() {
        let bits: BitVec<u8> = bit_vec![true];
        let rs = RankSelect::new(bits, 1);
        assert_eq!(rs.select_1(0), None);
        assert_eq!(rs.select_1(1), Some(0));
        assert_eq!(rs.select_0(0), Some(0));
        assert_eq!(rs.select_0(1), None);

        let bits: BitVec<u8> = bit_vec![false];
        let rs = RankSelect::new(bits, 1);
        assert_eq!(rs.select_1(1), None);
        assert_eq!(rs.select_1(0), Some(0));
        assert_eq!(rs.select_0(0), None);
        assert_eq!(rs.select_0(1), Some(0));
        assert_eq!(rs.rank_0(0), Some(1));
        assert_eq!(rs.rank_1(0), Some(0));
    }

    #[test]
    fn test_rank_k() {
        let mut bits: BitVec<u8> = BitVec::new_fill(false, 72);
        bits.set_bit(63, true);
        let rs = RankSelect::new(bits, 2);
        assert_eq!(rs.rank_1(63), Some(1));
        assert_eq!(rs.rank_1(64), Some(1));
        assert_eq!(rs.rank_1(71), Some(1));
    }
}
