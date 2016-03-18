// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Suffix arrays and related algorithms.
//! The implementation is based on the lecture notes
//! "Algorithmen auf Sequenzen", Kopczynski, Marschall, Martin and Rahmann, 2008 - 2015.

use std::iter;
use std;
use std::fmt::Debug;

use num::{Integer, Unsigned, NumCast};
use num::traits::cast;

use bit_vec::BitVec;
use vec_map::VecMap;

use alphabets::{Alphabet, RankTransform};
use data_structures::smallints::SmallInts;


pub type SuffixArray = Vec<usize>;
pub type SuffixArraySlice = [usize];
pub type LCPArray = SmallInts<i8, isize>;

/// Construct suffix array for given text of length n.
/// Complexity: O(n).
/// This is an implementation of the induced sorting as presented by
/// Ge Nong, Sen Zhang und Wai Hong Chan (2009), also known as SAIS.
/// The implementation is based on the following lecture notes:
/// http://ls11-www.cs.tu-dortmund.de/people/rahmann/algoseq.pdf
///
/// The idea is to first mark positions as L or S, with L being a position
/// the suffix of which is lexicographically larger than that of the next position.
/// Then, LMS-positions (leftmost S) are S-positions right to an L-position.
/// An LMS substring is the substring from one LMS position to the next (inclusive).
/// The algorithm works as follows:
///
/// 1. Sort LMS positions: first step 2 is applied to the unsorted sequence
///    of positions. Suprisingly, this sorts the LMS substrings. If all substrings
///    are different, LMS positions (and their suffixes) are sorted. Else, a reduced
///    text is build (at most half the size of the original text) and we recurse into
///    suffix array construction on the reduced text, yielding the sorted LMS positions.
/// 2. Derive the order of the other positions/suffixes from the (sorted) LMS positions.
///    For this, the (still empty) suffix array is partitioned into buckets.
///    Each bucket denotes an interval of suffixes with the same first symbol.
///    We know that the L-suffixes have to occur first in the buckets, because they
///    have to be lexicographically smaller than the S-suffixes with the same first letter.
///    The LMS-positions can now be used to insert the L-positions in the correct order
///    into the buckets.
///    Then, the S-positions can be inserted, again using the already existing entries
///    in the array.
///
/// # Arguments
///
/// * `text` - the text, ended by sentinel symbol (being lexicographically smallest). The text may
///   also contain multiple sentinel symbols, used to concatenate multiple sequences without mixing
///   their suffixes together.
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::suffix_array;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array(text);
/// assert_eq!(pos, vec![
///     21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18,
///     2, 16, 0, 19, 4, 13, 10, 3, 12, 9
/// ]);
/// ```
pub fn suffix_array(text: &[u8]) -> SuffixArray {
    let n = text.len();
    let alphabet = Alphabet::new(text);
    let sentinel_count = sentinel_count(text);
    let mut sais = SAIS::new(n);

    match alphabet.len() + sentinel_count {
        a if a <= std::u8::MAX as usize => {
            sais.construct(&transform_text::<u8>(text, &alphabet, sentinel_count))
        }
        a if a <= std::u16::MAX as usize => {
            sais.construct(&transform_text::<u16>(text, &alphabet, sentinel_count))
        }
        a if a <= std::u32::MAX as usize => {
            sais.construct(&transform_text::<u32>(text, &alphabet, sentinel_count))
        }
        _ => sais.construct(&transform_text::<u64>(text, &alphabet, sentinel_count)),
    }

    sais.pos
}


/// Construct lcp array for given text and suffix array of length n.
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
/// use bio::data_structures::suffix_array::{suffix_array,lcp};
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array(text);
///
/// // obtain compressed LCP array
/// let lcp = lcp(text, &pos);
///
/// // get most values in O(1).
/// assert_eq!(lcp.get(6).unwrap(), 4);
///
/// // obtain uncompressed LCP array.
/// let uncompressed = lcp.decompress();
/// assert_eq!(
///     uncompressed,
///     [
///         -1, 0, 1, 1, 2, 1, 4,
///         0, 1, 3, 1, 1, 2, 0,
///         4, 0, 2, 2, 2, 1, 3,
///         3, -1
///     ]
/// )
/// ```
pub fn lcp(text: &[u8], pos: &SuffixArraySlice) -> LCPArray {
    assert!(text.len() == pos.len());
    let n = text.len();

    // provide the lexicographical rank for each suffix
    let mut rank: Vec<usize> = iter::repeat(0).take(n).collect();
    for (r, p) in pos.iter().enumerate() {
        rank[*p] = r;
    }

    let mut lcp = SmallInts::from_elem(-1, n + 1);
    let mut l = 0usize;
    for (p, &r) in rank.iter().enumerate().take(n-1) {
        // since the sentinel has rank 0 and is excluded above,
        // we will never have a negative index below
        let pred = pos[r - 1];
        while pred + l < n && p + l < n && text[p + l] == text[pred + l] {
            l += 1;
        }
        lcp.set(r, l as isize);
        l = if l > 0 {
            l - 1
        } else {
            0
        };
    }

    lcp
}


/// Return last character of the text (expected to be the sentinel).
fn sentinel(text: &[u8]) -> u8 {
    text[text.len() - 1]
}


/// Count the sentinels occuring in the text given that the last character is the sentinel.
fn sentinel_count(text: &[u8]) -> usize {
    let sentinel = sentinel(text);
    assert!(text.iter().all(|&a| a >= sentinel),
            "Expecting extra sentinel symbol being lexicographically smallest at the end of the \
             text.");

    text.iter().fold(0, |count, &a| count + (a == sentinel) as usize)
}


/// Transform the given text into integers for usage in `SAIS`.
fn transform_text<T: Integer + Unsigned + NumCast + Copy>(text: &[u8],
                                                          alphabet: &Alphabet,
                                                          sentinel_count: usize)
                                                          -> Vec<T> {
    let sentinel = sentinel(text);
    let transform = RankTransform::new(alphabet);
    let offset = sentinel_count - 1;

    let mut transformed: Vec<T> = Vec::with_capacity(text.len());
    let mut s = sentinel_count;
    for &a in text.iter() {
        if a == sentinel {
            s -= 1;
            transformed.push(cast(s).unwrap());
        } else {
            transformed.push(cast(*(transform.ranks.get(a as usize)).unwrap() as usize + offset)
                                 .unwrap());
        }
    }

    transformed
}


/// SAIS implementation (see function `suffix_array` for description).
struct SAIS {
    pos: Vec<usize>,
    lms_pos: Vec<usize>,
    reduced_text_pos: Vec<usize>,
    bucket_sizes: VecMap<usize>,
    bucket_start: Vec<usize>,
    bucket_end: Vec<usize>,
}


impl SAIS {
    /// Create a new instance.
    fn new(n: usize) -> Self {
        SAIS {
            pos: Vec::with_capacity(n),
            lms_pos: Vec::with_capacity(n),
            reduced_text_pos: vec![0; n],
            bucket_sizes: VecMap::new(),
            bucket_start: Vec::with_capacity(n),
            bucket_end: Vec::with_capacity(n),
        }
    }

    /// Init buckets.
    fn init_bucket_start<T: Integer + Unsigned + NumCast + Copy>(&mut self, text: &[T]) {
        self.bucket_sizes.clear();
        self.bucket_start.clear();

        for &c in text.iter() {
            if !self.bucket_sizes.contains_key(cast(c).unwrap()) {
                self.bucket_sizes.insert(cast(c).unwrap(), 0);
            }
            *(self.bucket_sizes.get_mut(cast(c).unwrap()).unwrap()) += 1;
        }

        let mut sum = 0;
        for &size in self.bucket_sizes.values() {
            self.bucket_start.push(sum);
            sum += size;
        }
    }

    /// Initialize pointers to the last element of the buckets.
    fn init_bucket_end<T: Integer + Unsigned + NumCast + Copy>(&mut self, text: &[T]) {
        self.bucket_end.clear();
        for &r in self.bucket_start[1..].iter() {
            self.bucket_end.push(r - 1);
        }
        self.bucket_end.push(text.len() - 1);
    }

    /// Check if two LMS substrings are equal.
    fn lms_substring_eq<T: Integer + Unsigned + NumCast + Copy>(&self,
                                                                text: &[T],
                                                                pos_types: &PosTypes,
                                                                i: usize,
                                                                j: usize)
                                                                -> bool {
        for k in 0.. {
            let lmsi = pos_types.is_lms_pos(i + k);
            let lmsj = pos_types.is_lms_pos(j + k);
            if text[i + k] != text[j + k] {
                // different symbols
                return false;
            }
            if lmsi != lmsj {
                // different length
                return false;
            }
            if k > 0 && lmsi && lmsj {
                // same symbols and same length
                return true;
            }
        }
        false
    }

    /// Sort LMS suffixes.
    fn sort_lms_suffixes<T: Integer + Unsigned + NumCast + Copy + Debug,
                         S: Integer + Unsigned + NumCast + Copy + Debug>
        (&mut self,
         text: &[T],
         pos_types: &PosTypes,
         lms_substring_count: usize) {

        // if less than 2 LMS substrings are present, no further sorting is needed
        if lms_substring_count > 1 {
            // sort LMS suffixes by recursively building SA on reduced text
            let mut reduced_text: Vec<S> = vec![cast(0).unwrap(); lms_substring_count];
            let mut label = 0;
            reduced_text[self.reduced_text_pos[self.pos[0]]] = cast(label).unwrap();
            let mut prev = self.pos[0];
            for &p in self.pos[1..].iter() {
                if pos_types.is_lms_pos(p) {
                    // choose same label if substrings are equal
                    if !self.lms_substring_eq(text, pos_types, prev, p) {
                        label += 1;
                    }
                    reduced_text[self.reduced_text_pos[p]] = cast(label).unwrap();
                    prev = p;
                }
            }

            // if we have less labels than substrings, we have to sort by recursion
            // because two or more substrings are equal
            if label + 1 < lms_substring_count {
                // backup lms_pos
                let lms_pos = self.lms_pos.clone();
                // recurse SA construction for reduced text
                self.construct(&reduced_text);
                // obtain sorted lms suffixes
                self.lms_pos.clear();
                for &p in &self.pos {
                    self.lms_pos.push(lms_pos[p]);
                }
            }
        }
    }

    /// Construct the suffix array.
    fn construct<T: Integer + Unsigned + NumCast + Copy + Debug>(&mut self, text: &[T]) {
        let pos_types = PosTypes::new(text);
        self.calc_lms_pos(text, &pos_types);
        self.calc_pos(text, &pos_types);
    }

    /// Step 1 of the SAIS algorithm.
    fn calc_lms_pos<T: Integer + Unsigned + NumCast + Copy + Debug>(&mut self,
                                                                    text: &[T],
                                                                    pos_types: &PosTypes) {
        let n = text.len();

        // collect LMS positions
        self.lms_pos.clear();
        let mut i = 0;
        for r in 0..n {
            if pos_types.is_lms_pos(r) {
                self.lms_pos.push(r);
                self.reduced_text_pos[r] = i;
                i += 1;
            }
        }

        // sort LMS substrings by applying step 2 with unsorted LMS positions
        self.calc_pos(text, pos_types);

        let lms_substring_count = self.lms_pos.len();

        if lms_substring_count <= std::u8::MAX as usize {
            self.sort_lms_suffixes::<T, u8>(text, pos_types, lms_substring_count);
        } else if lms_substring_count <= std::u16::MAX as usize {
            self.sort_lms_suffixes::<T, u16>(text, pos_types, lms_substring_count);
        } else if lms_substring_count <= std::u32::MAX as usize {
            self.sort_lms_suffixes::<T, u32>(text, pos_types, lms_substring_count);
        } else {
            self.sort_lms_suffixes::<T, u64>(text, pos_types, lms_substring_count);
        }
    }

    /// Step 2 of the SAIS algorithm.
    fn calc_pos<T: Integer + Unsigned + NumCast + Copy>(&mut self,
                                                        text: &[T],
                                                        pos_types: &PosTypes) {
        let n = text.len();
        self.pos.clear();

        self.init_bucket_start(text);
        self.init_bucket_end(text);

        // init all positions as unknown (n-1 is max position)
        for _ in text.iter() {
            self.pos.push(n);
        }

        // insert LMS positions to the end of their buckets
        for &p in self.lms_pos.iter().rev() {
            let c: usize = cast(text[p]).unwrap();
            self.pos[self.bucket_end[c]] = p;
            // subtract without overflow: last -1 will cause overflow, but it does not matter
            self.bucket_end[c] = self.bucket_end[c].wrapping_sub(1);
        }

        // reset bucket ends
        self.init_bucket_end(text);

        // insert L-positions into buckets
        for r in 0..n {
            let p = self.pos[r];
            // ignore undefined positions and the zero since it has no predecessor
            if p == n || p == 0 {
                continue;
            }
            let pred = p - 1;
            if pos_types.is_l_pos(pred) {
                let c: usize = cast(text[pred]).unwrap();
                self.pos[self.bucket_start[c]] = pred;
                self.bucket_start[c] += 1;
            }
        }

        // insert S-positions into buckets
        for r in (0..n).rev() {
            let p = self.pos[r];
            if p == 0 {
                continue;
            }
            let pred = p - 1;
            if pos_types.is_s_pos(pred) {
                let c: usize = cast(text[pred]).unwrap();
                self.pos[self.bucket_end[c]] = pred;
                // subtract without overflow: last -1 will cause overflow, but it won't be used
                self.bucket_end[c] = self.bucket_end[c].wrapping_sub(1);
            }
        }
    }
}


/// Position types (L or S).
struct PosTypes {
    pos_types: BitVec,
}


impl PosTypes {
    /// Calculate the text position type.
    /// L-type marks suffixes being lexicographically larger than their successor,
    /// S-type marks the others.
    /// This function fills a BitVec, with 1-bits denoting S-type
    /// and 0-bits denoting L-type.
    ///
    /// # Arguments
    ///
    /// * `text` - the text, ending with a sentinel.
    fn new<T: Integer + Unsigned + NumCast + Copy>(text: &[T]) -> Self {
        let n = text.len();
        let mut pos_types = BitVec::from_elem(n, false);
        pos_types.set(n - 1, true);

        for p in (0..n - 1).rev() {
            if text[p] == text[p + 1] {
                // if the characters are equal, the next position determines
                // the lexicographical order
                let v = pos_types.get(p + 1).unwrap();
                pos_types.set(p, v);
            } else {
                pos_types.set(p, text[p] < text[p + 1]);
            }
        }

        PosTypes { pos_types: pos_types }
    }

    /// Check if p is S-position.
    fn is_s_pos(&self, p: usize) -> bool {
        self.pos_types.get(p).unwrap()
    }

    /// Check if p is L-position.
    fn is_l_pos(&self, p: usize) -> bool {
        !self.pos_types.get(p).unwrap()
    }

    /// Check if p is LMS-position.
    fn is_lms_pos(&self, p: usize) -> bool {
        !(p == 0) && self.is_s_pos(p) && self.is_l_pos(p - 1)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::{PosTypes, SAIS, transform_text};
    use bit_vec::BitVec;
    use alphabets::Alphabet;


    #[test]
    fn test_pos_types() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let pos_types = PosTypes::new(&text);
        let mut test = BitVec::from_bytes(&[0b01100110, 0b10010011, 0b01100100]);
        test.truncate(n);
        assert_eq!(pos_types.pos_types, test);
        let lms_pos: Vec<usize> = (0..n).filter(|&p| pos_types.is_lms_pos(p)).collect();
        assert_eq!(lms_pos, vec![1, 5, 8, 11, 14, 17, 21]);
    }


    #[test]
    fn test_buckets() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let mut sais = SAIS::new(n);
        sais.init_bucket_start(&text);
        assert_eq!(sais.bucket_start, vec![0, 1, 7, 13, 15]);
        sais.init_bucket_end(&text);
        assert_eq!(sais.bucket_end, vec![0, 6, 12, 14, 21]);
    }


    #[test]
    fn test_pos() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let mut sais = SAIS::new(n);
        let pos_types = PosTypes::new(&text);
        sais.lms_pos = vec![21, 5, 14, 8, 11, 17, 1];
        sais.calc_pos(&text, &pos_types);
        assert_eq!(sais.pos,
                   vec![21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18, 2, 16, 0, 19, 4, 13, 10, 3,
                        12, 9]);
    }


    #[test]
    fn test_lms_pos() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let mut sais = SAIS::new(n);
        let pos_types = PosTypes::new(&text);
        sais.calc_lms_pos(&text, &pos_types);
    }


    #[test]
    fn test_issue10_1() {
        let text = b"TGTGTGTGTG$";
        let pos = suffix_array(text);
        assert_eq!(pos, [10, 9, 7, 5, 3, 1, 8, 6, 4, 2, 0]);
    }


    #[test]
    fn test_issue10_2() {
        let text = b"TGTGTGTG$";
        let pos = suffix_array(text);
        assert_eq!(pos, [8, 7, 5, 3, 1, 6, 4, 2, 0]);
    }

    #[test]
    fn test_handles_sentinels_properly() {
        let reads = b"TACTCCGCTAGGGACACCTAAATAGATACTCGCAAAGGCGACTGATATATCCTTAGGTCGAAGAGATACCAGAGAAATAGTAGGTCTTAGGCTAGTCCTT$AAGGACTAGCCTAAGACCTACTATTTCTCTGGTATCTCTTCGACCTAAGGATATATCAGTCGCCTTTGCGAGTATCTATTTAGGTGTCCCTAGCGGAGTA$TAGGGACACCTAAATAGATACTCGCAAAGGCGACTGATATATCCTTAGGTCGAAGAGATACCAGAGAAATAGTAGGTCTTAGGCTAGTCCTTGTCCAGTA$TACTGGACAAGGACTAGCCTAAGACCTACTATTTCTCTGGTATCTCTTCGACCTAAGGATATATCAGTCGCCTTTGCGAGTATCTATTTAGGTGTCCCTA$ACGCACCCCGGCATTCGTCGACTCTACACTTAGTGGAACATACAAATTCGCTCGCAGGAGCGCCTCATACATTCTAACGCAGTGATCTTCGGCTGAGACT$AGTCTCAGCCGAAGATCACTGCGTTAGAATGTATGAGGCGCTCCTGCGAGCGAATTTGTATGTTCCACTAAGTGTAGAGTCGACGAATGCCGGGGTGCGT$";
        suffix_array(reads);
    }
}
