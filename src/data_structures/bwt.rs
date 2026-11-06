// Copyright 2014-2016 Johannes Köster, Taylor Cramer.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! The [Burrows-Wheeler-Transform](https://www.semanticscholar.org/paper/A-Block-sorting-Lossless-Data-Compression-Algorithm-Burrows-Wheeler/af56e6d4901dcd0f589bf969e604663d40f1be5d) and related data structures.
//! The implementation is based on the lecture notes
//! "Algorithmen auf Sequenzen", Kopczynski, Marschall, Martin and Rahmann, 2008 - 2015.

use std::iter::repeat_n;

use crate::alphabets::Alphabet;
use crate::data_structures::suffix_array::RawSuffixArraySlice;
use crate::utils::prescan;

pub type BWT = Vec<u8>;
pub type BWTSlice = [u8];
pub type Less = Vec<usize>;
pub type BWTFind = Vec<usize>;

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
/// use bio::data_structures::bwt::bwt;
/// use bio::data_structures::suffix_array::suffix_array;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array(text);
/// let bwt = bwt(text, &pos);
/// assert_eq!(bwt, b"ATTATTCAGGACCC$CTTTCAA");
/// ```
pub fn bwt(text: &[u8], pos: RawSuffixArraySlice) -> BWT {
    assert_eq!(text.len(), pos.len());
    let n = text.len();
    let mut bwt: BWT = repeat_n(0, n).collect();
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
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Occ {
    occ: Vec<Vec<usize>>,
    k: u32,
}

impl Occ {
    /// Calculate occ array with sampling from BWT of length n.
    /// Time complexity: O(n).
    /// Space complexity: O(n / k * A) with A being the alphabet size.
    /// The specified alphabet must match the alphabet of the text and its BWT.
    /// For large texts, it is advisable to transform
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
            .expect("Expecting non-empty alphabet.") as usize
            + 1;
        let mut alpha = alphabet.symbols.iter().collect::<Vec<usize>>();
        // include sentinel '$'
        if (b'$' as usize) < m && !alphabet.is_word(b"$") {
            alpha.push(b'$' as usize);
        }
        let mut occ: Vec<Vec<usize>> = vec![Vec::new(); m];
        let mut curr_occ = vec![0usize; m];

        // characters not in the alphabet won't take up much space
        for &a in &alpha {
            occ[a].reserve(n / k as usize);
        }

        for (i, &c) in bwt.iter().enumerate() {
            curr_occ[c as usize] += 1;

            if i % k as usize == 0 {
                // only visit characters in the alphabet
                for &a in &alpha {
                    occ[a].push(curr_occ[a]);
                }
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
        // let count = bwt[(i * self.k) + 1..r + 1]
        //     .iter()
        //     .filter(|&&c| c == a)
        //     .count();
        // self.occ[a as usize][i] + count
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

        // self.k is our sampling rate, so find the checkpoints either side of r.
        let lo_checkpoint = r / self.k as usize;
        // Get the occurences at the low checkpoint
        let lo_occ = self.occ[a as usize][lo_checkpoint];

        // If the sampling rate is infrequent it is worth checking if there is a closer
        // hi checkpoint.
        if self.k > 64 {
            let hi_checkpoint = lo_checkpoint + 1;
            if let Some(&hi_occ) = self.occ[a as usize].get(hi_checkpoint) {
                // Its possible that there are no occurences between the low and high
                // checkpoint in which case we bail early.
                if lo_occ == hi_occ {
                    return lo_occ;
                }

                // If r is closer to the high checkpoint, count backwards from there.
                let hi_idx = hi_checkpoint * self.k as usize;
                if (hi_idx - r) < (self.k as usize / 2) {
                    return hi_occ - bytecount::count(&bwt[r + 1..=hi_idx], a);
                }
            }
        }

        // Otherwise the default case is to count from the low checkpoint.
        let lo_idx = lo_checkpoint * self.k as usize;
        bytecount::count(&bwt[lo_idx + 1..=r], a) + lo_occ
    }
}

/// Calculate the less array for a given BWT. Complexity O(n).
pub fn less(bwt: &BWTSlice, alphabet: &Alphabet) -> Less {
    let m = alphabet
        .max_symbol()
        .expect("Expecting non-empty alphabet.") as usize
        + 2;
    let mut less: Less = repeat_n(0, m).collect();
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

    let mut bwtfind: BWTFind = repeat_n(0, n).collect();
    for (r, &c) in bwt.iter().enumerate() {
        bwtfind[less[c as usize]] = r;
        less[c as usize] += 1;
    }

    bwtfind
}

#[cfg(test)]
mod tests {
    use super::{bwt, bwtfind, invert_bwt, Occ};
    use crate::alphabets::dna;
    use crate::alphabets::Alphabet;
    use crate::data_structures::suffix_array::suffix_array;
    use crate::data_structures::wavelet_matrix::WaveletMatrix;

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
        let alphabet = Alphabet::new([0u8, 1u8, 2u8, 3u8]);
        let occ = Occ::new(&bwt, 3, &alphabet);
        assert_eq!(occ.occ, [[0, 0], [1, 2], [0, 0], [0, 2]]);
        assert_eq!(occ.get(&bwt, 4, 2u8), 1);
        assert_eq!(occ.get(&bwt, 4, 3u8), 2);
    }

    #[test]
    fn test_occwm() {
        let text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = {
            let mut a = dna::n_alphabet();
            a.insert(b'$');
            a
        };
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let occ = Occ::new(&bwt, 3, &alphabet);
        let wm = WaveletMatrix::new(&bwt);

        for c in [b'A', b'C', b'G', b'T', b'$'] {
            for p in 0..text.len() {
                assert_eq!(occ.get(&bwt, p, c) as u64, wm.rank(c, p as u64));
            }
        }
    }
}
