// Copyright 2014-2016 Johannes Köster, Taylor Cramer.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! The Burrows-Wheeler-Transform and related data structures.
//! The implementation is based on the lecture notes
//! "Algorithmen auf Sequenzen", Kopczynski, Marschall, Martin and Rahmann, 2008 - 2015.

use std::iter::repeat;

use alphabets::Alphabet;
use bytecount;
use data_structures::suffix_array::RawSuffixArray;
use std::ops::Deref;
use utils::prescan;

pub type BWT = Vec<u8>;
pub type BWTSlice = [u8];
pub type Less = Vec<usize>;
pub type BWTFind = Vec<usize>;

pub trait DerefBWT: Deref<Target = BWT> {}
impl<T: Deref<Target = BWT>> DerefBWT for T {}

pub trait DerefOcc: Deref<Target = Occ> {}
impl<T: Deref<Target = Occ>> DerefOcc for T {}

pub trait DerefLess: Deref<Target = Less> {}
impl<T: Deref<Target = Less>> DerefLess for T {}

/// Calculate Burrows-Wheeler-Transform of the given text of length n.
/// Complexity: O(n).
///
/// # Arguments
///
/// * `text` - the text ended by sentinel symbol (being lexicographically smallest)
/// * `pos` - the suffix array for the text
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::suffix_array;
/// use bio::data_structures::bwt::bwt;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array(text);
/// let bwt = bwt(text, &pos);
/// assert_eq!(bwt, b"ATTATTCAGGACCC$CTTTCAA");
/// ```
pub fn bwt(text: &[u8], pos: &RawSuffixArray) -> BWT {
    assert_eq!(text.len(), pos.len());
    let n = text.len();
    let mut bwt: BWT = repeat(0).take(n).collect();
    for r in 0..n {
        let p = pos[r];
        bwt[r] = if p > 0 { text[p - 1] } else { text[n - 1] };
    }

    bwt
}

/// Calculate the inverse of a BWT of length n, which is the original text.
/// Complexity: O(n).
///
/// This only works if the last sentinel in the original text is unique
/// and lexicographically the smallest.
///
/// # Arguments
///
/// * `bwt` - the BWT
pub fn invert_bwt(bwt: &BWTSlice) -> Vec<u8> {
    let alphabet = Alphabet::new(bwt);
    let n = bwt.len();
    let bwtfind = bwtfind(bwt, &alphabet);
    let mut inverse = Vec::with_capacity(n);

    let mut r = bwtfind[0];
    for _ in 0..n {
        r = bwtfind[r];
        inverse.push(bwt[r]);
    }

    inverse
}

/// An occurrence array implementation.
#[derive(Serialize, Deserialize)]
pub struct Occ {
    occ: Vec<Vec<usize>>,
    k: u32,
}

impl Occ {
    /// Calculate occ array with sampling from BWT of length n.
    /// Time complexity: O(n).
    /// Space complexity: O(n / k * A) with A being the alphabet size.
    /// Alphabet size is determined on the fly from the BWT.
    /// For large texts, it is therefore advisable to transform
    /// the text before calculating the BWT (see alphabets::rank_transform).
    ///
    /// # Arguments
    ///
    /// * `bwt` - the BWT
    /// * `k` - the sampling rate: every k-th entry will be stored
    pub fn new(bwt: &BWTSlice, k: u32, alphabet: &Alphabet) -> Self {
        let n = bwt.len();
        let m = alphabet
            .max_symbol()
            .expect("Expecting non-empty alphabet.") as usize + 1;
        let mut occ = Vec::with_capacity(n / k as usize);
        let mut curr_occ: Vec<usize> = repeat(0).take(m).collect();
        for (i, &c) in bwt.iter().enumerate() {
            curr_occ[c as usize] += 1;
            if i % k as usize == 0 {
                occ.push(curr_occ.clone());
            }
        }

        Occ { occ, k }
    }

    /// Get occurrence count of symbol a in BWT[..r+1].
    /// Complexity: O(k).
    pub fn get(&self, bwt: &BWTSlice, r: usize, a: u8) -> usize {
        // NOTE:
        //
        // Retrieving byte match counts in this function is critical to the performance of FM Index.
        //
        // The below manual count code is roughly equivalent to:
        // ```
        // let count = bwt[(i * self.k) + 1..r + 1].iter().filter(|&&c| c == a).count();
        // self.occ[i][a as usize] + count
        // ```
        //
        // But there are a couple of reasons to do this manually:
        // 1) As of 2016, versions of rustc/LLVM vectorize this manual loop more reliably
        //    than the iterator adapter version.
        // 2) Manually accumulating the byte match count in a single chunk can allows
        //    us to use a `u32` for that count, which has faster arithmetic on common arches.
        //    This does necessitate storing `k` as a u32.
        //
        // See the conversation in these issues for some of the history here:
        //
        // https://github.com/rust-bio/rust-bio/pull/74
        // https://github.com/rust-bio/rust-bio/pull/76

        // self.k is our sampling rate, so find our last sampled checkpoint
        let i = r / self.k as usize;
        let checkpoint = self.occ[i][a as usize];

        // find the portion of the BWT past the checkpoint which we need to count
        let start = (i * self.k as usize) + 1;
        let end = r + 1;

        // count all the matching bytes b/t the closest checkpoint and our desired lookup
        let count = bytecount::count(&bwt[start..end], a);

        // return the sampled checkpoint for this character + the manual count we just did
        checkpoint + (count as usize)
    }
}

/// Calculate the less array for a given BWT. Complexity O(n).
pub fn less(bwt: &BWTSlice, alphabet: &Alphabet) -> Less {
    let m = alphabet
        .max_symbol()
        .expect("Expecting non-empty alphabet.") as usize + 2;
    let mut less: Less = repeat(0).take(m).collect();
    for &c in bwt.iter() {
        less[c as usize] += 1;
    }
    // calculate +-prescan
    prescan(&mut less[..], 0, |a, b| a + b);

    less
}

/// Calculate the bwtfind array needed for inverting the BWT. Complexity O(n).
pub fn bwtfind(bwt: &BWTSlice, alphabet: &Alphabet) -> BWTFind {
    let n = bwt.len();
    let mut less = less(bwt, alphabet);

    let mut bwtfind: BWTFind = repeat(0).take(n).collect();
    for (r, &c) in bwt.iter().enumerate() {
        bwtfind[less[c as usize]] = r;
        less[c as usize] += 1;
    }

    bwtfind
}

#[cfg(test)]
mod tests {
    use super::{bwt, bwtfind, invert_bwt, Occ};
    use alphabets::Alphabet;
    use data_structures::suffix_array::suffix_array;

    #[test]
    fn test_bwtfind() {
        let text = b"cabca$";
        let alphabet = Alphabet::new(b"abc$");
        let pos = suffix_array(text);
        let bwt = bwt(text, &pos);
        let bwtfind = bwtfind(&bwt, &alphabet);
        assert_eq!(bwtfind, vec![5, 0, 3, 4, 1, 2]);
    }

    #[test]
    fn test_invert_bwt() {
        let text = b"cabca$";
        let pos = suffix_array(text);
        let bwt = bwt(text, &pos);
        let inverse = invert_bwt(&bwt);
        assert_eq!(inverse, text);
    }

    #[test]
    fn test_occ() {
        let bwt = vec![1u8, 3u8, 3u8, 1u8, 2u8, 0u8];
        let alphabet = Alphabet::new(&[0u8, 1u8, 2u8, 3u8]);
        let occ = Occ::new(&bwt, 3, &alphabet);
        assert_eq!(occ.occ, [[0, 1, 0, 0], [0, 2, 0, 2]]);
        assert_eq!(occ.get(&bwt, 4, 2u8), 1);
        assert_eq!(occ.get(&bwt, 4, 3u8), 2);
    }
}
