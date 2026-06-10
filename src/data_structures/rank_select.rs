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
use std::ops::Deref;

use bv::BitVec;
use bv::Bits;

/// A rank/select data structure.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct RankSelect {
    n: usize,
    bits: BitVec<u8>,
    superblocks_1: Vec<SuperblockRank>,
    /// superblock size in bits
    s: usize,
    /// superblock size in 32 bits
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
            superblocks_1: superblocks(n, s, &bits),
            bits,
        }
    }

    /// Append a bit to the end of the underlying bit vector, updating the
    /// rank/select index in amortized `O(1)` time.
    ///
    /// Appending a bit cannot change the rank of any earlier position, so the
    /// existing superblock samples stay valid and only an occasional new sample
    /// has to be added (once every `s = k * 32` bits). This makes it possible to
    /// build a `RankSelect` incrementally, e.g. when the bits are produced by a
    /// streaming computation or when using it as the bitmap backing a wavelet
    /// tree.
    ///
    /// The resulting structure is identical to one built with
    /// [`RankSelect::new`] from the same bits and the same `k`.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::data_structures::rank_select::RankSelect;
    /// use bv::BitVec;
    ///
    /// let mut rs = RankSelect::new(BitVec::new(), 1);
    /// rs.push(true);
    /// rs.push(false);
    /// rs.push(true);
    /// assert_eq!(rs.rank_1(2), Some(2));
    /// assert_eq!(rs.select_1(1), Some(0));
    /// ```
    pub fn push(&mut self, bit: bool) {
        // A new superblock sample starts whenever the current length is a
        // multiple of the superblock size `s` (this also covers the very first
        // bit, where `n == 0`). The sample stores the number of 1-bits strictly
        // before that boundary, which is the current total number of 1-bits.
        if self.n.is_multiple_of(self.s) {
            let ones_before = if self.n == 0 {
                0
            } else {
                self.rank_1((self.n - 1) as u64)
                    .expect("position is within bounds before pushing")
            };
            // Mirror the `First`/`Some` encoding used by `superblocks`: a sample
            // is `First` when its rank differs from the previous sample (or it is
            // the first one), and `Some` when the rank is unchanged.
            let sample = match self.superblocks_1.last() {
                Some(last) if **last == ones_before => SuperblockRank::Some(ones_before),
                _ => SuperblockRank::First(ones_before),
            };
            self.superblocks_1.push(sample);
        }

        self.bits.push(bit);
        self.n += 1;
    }

    /// Append all bits yielded by an iterator to the end of the underlying bit
    /// vector, updating the rank/select index as it goes.
    ///
    /// This is equivalent to calling [`RankSelect::push`] once per bit and runs
    /// in amortized `O(1)` per bit. It accepts anything that iterates over
    /// `bool`, so it can extend a `RankSelect` by another `BitVec`, a
    /// `Vec<bool>`, or any lazy bit stream.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::data_structures::rank_select::RankSelect;
    /// use bv::BitVec;
    ///
    /// let mut rs = RankSelect::new(BitVec::new(), 1);
    /// rs.extend([true, false, true, true]);
    /// assert_eq!(rs.rank_1(3), Some(3));
    /// assert_eq!(rs.select_1(2), Some(2));
    /// ```
    pub fn extend<I: IntoIterator<Item = bool>>(&mut self, bits: I) {
        for bit in bits {
            self.push(bit);
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
            let mut rank = *self.superblocks_1[s as usize];
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
            |i| *self.superblocks_1[i],
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
            // Derived from `superblocks_1`: at the start of superblock `i`
            // we have seen `i * s` bits, of which `superblocks_1[i]` are 1s.
            |i| (i * self.s) as u64 - *self.superblocks_1[i],
            |bit| bit == 0,
            |block| block.count_zeros(),
        )
    }

    /// Generic implementation for `select_1` and `select_0`.
    ///
    /// `rank_at_superblock(i)` returns the cumulative rank (1s for `select_1`,
    /// 0s for `select_0`) at the start of superblock `i`. `is_match` and
    /// `count_all` describe how to inspect a single byte for the appropriate
    /// bit-class.
    fn select_x<R: Fn(usize) -> u64, F: Fn(u8) -> bool, C: Fn(u8) -> u32>(
        &self,
        j: u64,
        rank_at_superblock: R,
        is_match: F,
        count_all: C,
    ) -> Option<u64> {
        if j == 0 {
            return None;
        }

        // Find the smallest superblock index whose stored rank is `>= j`.
        // Manual binary search rather than slice::binary_search so we can
        // search a virtual array — `select_0` derives its values from
        // `superblocks_1` rather than storing them.
        let n_super = self.superblocks_1.len();
        let mut lo = 0usize;
        let mut hi = n_super;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if rank_at_superblock(mid) < j {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        // `superblocks[k]` holds the rank *before* superblock k starts, so the
        // j-th bit lives in superblock `lo - 1` (or 0 if j is small enough
        // that the very first superblock already meets the rank).
        let superblock = lo.saturating_sub(1);
        let mut rank = rank_at_superblock(superblock);

        // Scan blocks within the chosen superblock byte-by-byte, popcounting
        // each block. Once the running rank would cross `j` inside a block,
        // walk that block bit-by-bit to find the exact position.
        let first_block = superblock * self.s / 8;
        for block in first_block..cmp::min(first_block + self.s / 8, self.bits.block_len()) {
            let b = self.bits.get_block(block);
            let p = count_all(b) as u64;
            if rank + p >= j {
                let mut bit = 0b1;
                // The final block may extend past the bitvector; clamp to
                // `bits.len()` so trailing zero-padding bits aren't counted.
                let max_bit = cmp::min(8, self.bits.len() - block as u64 * 8);
                for i in 0..max_bit {
                    rank += is_match(b & bit) as u64;
                    if rank == j {
                        return Some(block as u64 * 8 + i);
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

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
pub enum SuperblockRank {
    First(u64),
    Some(u64),
}

impl Deref for SuperblockRank {
    type Target = u64;

    fn deref(&self) -> &u64 {
        match self {
            SuperblockRank::First(rank) => rank,
            SuperblockRank::Some(rank) => rank,
        }
    }
}

impl PartialOrd for SuperblockRank {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SuperblockRank {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        let cmp = (**self).cmp(&**other);
        if cmp == cmp::Ordering::Equal {
            match (self, other) {
                (SuperblockRank::First(_), SuperblockRank::Some(_)) => cmp::Ordering::Less,
                (SuperblockRank::Some(_), SuperblockRank::First(_)) => cmp::Ordering::Greater,
                _ => cmp,
            }
        } else {
            cmp
        }
    }
}

/// Build the 1-rank superblock samples for a bitvector of `n` bits with
/// superblock size `s` (in bits). Each entry stores the cumulative rank of
/// 1s at the start of the corresponding superblock. The 0-rank samples can
/// be derived as `i * s - superblocks[i].rank()` and so are not stored.
fn superblocks(n: usize, s: usize, bits: &BitVec<u8>) -> Vec<SuperblockRank> {
    let mut superblocks = Vec::with_capacity(n / s + 1);
    let mut rank: u64 = 0;
    let mut last_rank = None;
    let mut i = 0;
    let nblocks = (bits.len() as f64 / 8.0).ceil() as usize;
    for block in 0..nblocks {
        let b = bits.get_block(block);
        if i % s == 0 {
            superblocks.push(if Some(rank) != last_rank {
                SuperblockRank::First(rank)
            } else {
                SuperblockRank::Some(rank)
            });
            last_rank = Some(rank);
        }
        rank += b.count_ones() as u64;
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
    fn test_select_start() {
        let mut bits: BitVec<u8> = BitVec::new_fill(false, 900);
        bits.set_bit(64, true);

        let rs = RankSelect::new(bits, 1);

        assert_eq!(rs.select_1(1), Some(64));
    }

    #[test]
    fn test_select_end() {
        let mut bits: BitVec<u8> = BitVec::new_fill(false, 900);
        bits.set_bit(50, true);

        let rs = RankSelect::new(bits, 1);
        assert_eq!(rs.select_1(1).unwrap(), 50);
    }

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
        assert_eq!(rs.select_1(0), None);
        assert_eq!(rs.select_1(1).unwrap(), 5);
        assert_eq!(rs.select_1(2).unwrap(), 32);
        assert_eq!(rs.rank_0(1).unwrap(), 2);
        assert_eq!(rs.rank_0(4).unwrap(), 5);
        assert_eq!(rs.rank_0(5).unwrap(), 5);
        assert_eq!(rs.select_0(0), None);
        assert_eq!(rs.select_0(1).unwrap(), 0);
        assert!(rs.get(5));
        assert!(!rs.get(1));
        assert!(rs.get(32));
    }

    #[test]
    fn test_rank_select2() {
        let mut bits: BitVec<u8> = BitVec::new_fill(false, 64);
        bits.set_bit(5, true);
        bits.set_bit(32, true);
        let rs = RankSelect::new(bits, 1);
        assert_eq!(rs.select_1(2).unwrap(), 32);
    }

    #[test]
    fn test_select() {
        let bits: BitVec<u8> = bit_vec![true, false];
        let rs = RankSelect::new(bits, 1);

        assert_eq!(rs.select_0(0), None);
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
        assert_eq!(rs.select_0(0), None);
        assert_eq!(rs.select_0(1), None);

        let bits: BitVec<u8> = bit_vec![false];
        let rs = RankSelect::new(bits, 1);
        assert_eq!(rs.select_1(1), None);
        assert_eq!(rs.select_1(0), None);
        assert_eq!(rs.select_0(0), None);
        assert_eq!(rs.select_0(1), Some(0));
        assert_eq!(rs.rank_0(0), Some(1));
        assert_eq!(rs.rank_1(0), Some(0));
    }

    // Cross-checks `select_1` and `select_0` against a naive linear scan on
    // a sparse bitvector. Exercises the duplicate-rank path for `select_1`
    // and the derived-rank path for `select_0` (see issue #548).
    #[test]
    fn test_select_against_naive_sparse() {
        let mut bits: BitVec<u8> = BitVec::new_fill(false, 1024);
        let one_positions: &[u64] = &[3, 70, 71, 72, 500, 900, 901, 1023];
        for &p in one_positions {
            bits.set_bit(p, true);
        }
        let zero_positions: Vec<u64> = (0..1024).filter(|i| !one_positions.contains(i)).collect();

        for k in [1usize, 2, 4, 8] {
            let rs = RankSelect::new(bits.clone(), k);
            for (i, &expected) in one_positions.iter().enumerate() {
                assert_eq!(rs.select_1((i + 1) as u64), Some(expected), "k={}", k);
            }
            assert_eq!(rs.select_1(one_positions.len() as u64 + 1), None);
            for (i, &expected) in zero_positions.iter().enumerate() {
                assert_eq!(rs.select_0((i + 1) as u64), Some(expected), "k={}", k);
            }
            assert_eq!(rs.select_0(zero_positions.len() as u64 + 1), None);
        }
    }

    #[test]
    fn test_select_against_naive_randomized() {
        use rand::{RngExt, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(0xdead_beef);
        for _ in 0..50 {
            let n: u64 = 64 + rng.random_range(0..4096);
            let mut bits: BitVec<u8> = BitVec::new_fill(false, n);
            let mut ones: Vec<u64> = Vec::new();
            let mut zeros: Vec<u64> = Vec::new();
            for i in 0..n {
                if rng.random_range(0..16) == 0 {
                    bits.set_bit(i, true);
                    ones.push(i);
                } else {
                    zeros.push(i);
                }
            }
            for &k in &[1usize, 2, 4] {
                let rs = RankSelect::new(bits.clone(), k);
                for (i, &expected) in ones.iter().enumerate() {
                    assert_eq!(rs.select_1((i + 1) as u64), Some(expected));
                }
                for (i, &expected) in zeros.iter().enumerate() {
                    assert_eq!(rs.select_0((i + 1) as u64), Some(expected));
                }
            }
        }
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

    #[test]
    fn test_push_basic_rank_select() {
        let mut rs = RankSelect::new(BitVec::new(), 1);
        // Build the pattern 1 0 1 1 0.
        for &b in &[true, false, true, true, false] {
            rs.push(b);
        }
        assert_eq!(rs.rank_1(0), Some(1));
        assert_eq!(rs.rank_1(4), Some(3));
        assert_eq!(rs.rank_0(4), Some(2));
        assert_eq!(rs.select_1(1), Some(0));
        assert_eq!(rs.select_1(3), Some(3));
        assert_eq!(rs.select_0(2), Some(4));
    }

    /// The central guarantee: a `RankSelect` built incrementally with `push`
    /// must be byte-for-byte identical (superblocks included) to one built with
    /// `new` from the same bits. `RankSelect` derives `PartialEq`, so this
    /// compares the full internal state, across several sampling rates `k` and
    /// lengths that straddle superblock and block boundaries.
    #[test]
    fn test_push_equivalent_to_new() {
        // A small deterministic xorshift keeps this test free of `rand` while
        // still exercising irregular bit patterns. The exact sequence does not
        // matter; only that `push` and `new` see the same bits.
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut next_bit = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state & 1 == 1
        };
        for &k in &[1usize, 2, 4] {
            for n in [0u64, 1, 7, 8, 9, 31, 32, 33, 63, 64, 65, 200, 257] {
                let mut bits: BitVec<u8> = BitVec::new();
                let mut pushed = RankSelect::new(BitVec::new(), k);
                for _ in 0..n {
                    let b = next_bit();
                    bits.push(b);
                    pushed.push(b);
                }
                let built = RankSelect::new(bits, k);
                assert_eq!(
                    pushed, built,
                    "push-built and new-built differ for k={}, n={}",
                    k, n
                );
            }
        }
    }

    #[test]
    fn test_extend_basic() {
        let mut rs = RankSelect::new(BitVec::new(), 1);
        rs.extend([true, false, true, true, false]);
        assert_eq!(rs.rank_1(4), Some(3));
        assert_eq!(rs.rank_0(4), Some(2));
        assert_eq!(rs.select_1(3), Some(3));
        assert_eq!(rs.select_0(2), Some(4));
    }

    /// `extend` must be exactly equivalent to calling `push` once per bit, and
    /// extending an existing structure must match building the concatenation in
    /// one shot with `new`. Both invariants are checked across sampling rates
    /// and split points around block/superblock boundaries.
    #[test]
    fn test_extend_equivalent_to_push_and_new() {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next_bit = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state & 1 == 1
        };
        for &k in &[1usize, 2, 4] {
            for split in [0usize, 1, 8, 31, 32, 33, 64, 100] {
                for tail in [0usize, 1, 7, 32, 65] {
                    // Generate the prefix and the tail bits.
                    let prefix: Vec<bool> = (0..split).map(|_| next_bit()).collect();
                    let extra: Vec<bool> = (0..tail).map(|_| next_bit()).collect();

                    // (a) extend vs push-one-by-one, starting from the same prefix.
                    let mut by_extend = RankSelect::new(BitVec::new(), k);
                    by_extend.extend(prefix.iter().copied());
                    let mut by_push = by_extend.clone();
                    by_extend.extend(extra.iter().copied());
                    for &b in &extra {
                        by_push.push(b);
                    }
                    assert_eq!(by_extend, by_push, "extend != push loop (k={k})");

                    // (b) extended structure vs the concatenation built with `new`.
                    let mut all_bits: BitVec<u8> = BitVec::new();
                    for &b in prefix.iter().chain(extra.iter()) {
                        all_bits.push(b);
                    }
                    let built = RankSelect::new(all_bits, k);
                    assert_eq!(by_extend, built, "extend != new (k={k})");
                }
            }
        }
    }
}
