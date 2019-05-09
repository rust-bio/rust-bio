// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Myers bit-parallel approximate pattern matching algorithm.
//! Finds all matches up to a given edit distance. The pattern has to fit into a bitvector,
//! and is thus limited to 64 or (since stable Rust version 1.26) to 128 symbols.
//! Complexity: O(n)
//!
//! # Example
//!
//! Iterating over matches in pairs of `(end, distance)` using `u64` as bitvector type:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let myers = Myers::<u64>::new(pattern);
//! let occ: Vec<_> = myers.find_all_end(text, 2).collect();
//!
//! assert_eq!(occ, [(11, 2), (12, 2)]);
//! # }
//! ```
//!
//! Starting with stable Rust 1.26, it is also possible to use `u128` as bitvector (`Myers128`),
//! which enables longer patterns, but is somewhat slower.
//!
//! # Obtaining the starting position of a match
//!
//! The `Myers::find_all` method provides an iterator over tuples of `(start, end, distance)`.
//! Calculating the starting position requires finding the alignment path, therefore this is
//! slower than `Myers::find_all_end`. Note that the end positions differ from above by one.
//! This is intentional, as the iterator returns a range rather an index, and ranges in Rust
//! do not include the end position by default.
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers::<u64>::new(pattern);
//! let occ: Vec<_> = myers.find_all(text, 2).collect();
//!
//! assert_eq!(occ, [(3, 12, 2), (3, 13, 2)]);
//! # }
//! ```
//!
//! # Obtaining alignments
//!
//! [`FullMatches`](struct.FullMatches.html) returned by `Myers::find_all()` also provide a method
//! for obtaining an alignment path:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers;
//! use bio::alignment::Alignment;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers::<u64>::new(pattern);
//! // create an 'empty' alignment instance, which can be reused
//! let mut aln = Alignment::default();
//!
//! let mut matches = myers.find_all(text, 3);
//! while matches.next_alignment(&mut aln) {
//!     println!("Hit fond in range: {}..{} (distance: {})", aln.ystart, aln.yend, aln.score);
//!     println!("{}", aln.pretty(pattern, text));
//! }
//! # }
//! ```
//! **Output:**
//!
//! <pre>
//! Hit fond in range: 3..10 (distance: 3)
//!    TCCTAGGGC
//!    ||||+|\|+
//! TCCTCCT-GAG-GGATTAGCAC
//!
//! Hit fond in range: 3..11 (distance: 3)
//!    TCCTAGGGC
//!    ||||+|\|\
//! TCCTCCT-GAGGGATTAGCAC
//!
//! Hit fond in range: 3..12 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||+
//! TCCTCCTGAGGG-ATTAGCAC
//!
//! Hit fond in range: 3..13 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||\
//! TCCTCCTGAGGGATTAGCAC
//!
//! ... (truncated)
//!
//! </pre>
//!
//! **Note** that the [`Alignment`](../../alignment/struct.Alignment.html) instance is only created
//! once and then reused. Because the Myers algorithm is very fast, the allocation necessary for
//! `Alignment::operations` can have a non-negligible impact on performance; and thus, recycling
//! makes sense.
//!
//! # Finding the best hit
//!
//! In many cases, only the match with the smallest edit distance is actually of interest.
//! Calculating an alignment for every hit is therefore not necessary.
//! [`LazyMatches`](struct.LazyMatches.html) returned by `Myers::find_all_lazy()`
//! provide an iterator over tuples of `(end, distance)` like `Myers::find_all_end()`, but
//! additionally keep the data necessary for calculating the alignment path later at any desired
//! position. Storing the data itself has a slight performance impact and requires more memory
//! compared to `Myers::find_all_end()` [O(n) as opposed to O(m + k)]. Still the following code
//! is faster than using `FullMatches`:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers;
//! use bio::alignment::Alignment;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers::<u64>::new(pattern);
//! let mut aln = Alignment::default();
//!
//! let mut matches = myers.find_all_lazy(text, 2);
//!
//! // first, find the best hit
//! let (best_end, _) = matches
//!     .by_ref()
//!     .min_by_key(|&(_, dist)| dist)
//!     .unwrap();
//!
//! // now calculate the alignment
//! matches.alignment_at(best_end, &mut aln);
//! println!("Best alignment at {}..{} (distance: {})", aln.ystart, aln.yend, aln.score);
//! println!("{}", aln.pretty(pattern, text));
//! # }
//! ```
//!
//! **Output:**
//!
//! <pre>
//! Best alignment at 3..12 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||+
//! TCCTCCTGAGGG-ATTAGCAC
//! </pre>
//!
//! Actually as seen in the previous chapters, there are two hits with the same distance of 2.
//! It may make sense to consider both of them.
//!
//! # Dealing with ambiguities
//!
//! Matching multiple or all symbols at once can be achieved using `MyersBuilder`. This example
//! allows `N` in the search pattern to match all four DNA bases in the text:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::MyersBuilder;
//!
//! # fn main() {
//! let text = b"GTCTGATCTTACC";
//! let pattern = b"TGATCNT";
//!
//! let myers = MyersBuilder::new()
//!     .ambig(b'N', b"ACGT")
//!     .build_64(pattern);
//! assert_eq!(myers.distance(text), 0);
//! # }
//! ```
//!
//! For more examples see the documentation of [`MyersBuilder`](struct.MyersBuilder.html).

use std::borrow::Borrow;
use std::cmp::min;
use std::collections::HashMap;
use std::iter;
use std::mem::size_of;
use std::ops::Range;
use std::ops::*;
use std::u64;

use num_traits::{Bounded, FromPrimitive, One, PrimInt, ToPrimitive, WrappingAdd, Zero};

use crate::alignment::{Alignment, AlignmentMode, AlignmentOperation};

mod builder;
mod traceback;

pub use self::builder::MyersBuilder;
use self::traceback::Traceback;

/// This trait must be implemented for integer types serving as bit vectors.
/// Only unsigned integers will work correctly.
pub trait BitVec: Copy
    + Default
    + Add
    + Sub
    + BitOr
    + BitOrAssign
    + BitAnd
    + BitXor
    + Not
    + Shl<usize>
    + ShlAssign<usize>
    + ShrAssign<usize>
    // These num_traits traits are required; in addition there are Bounded, Zero and One,
    // which are all required by PrimInt and thus included
    + PrimInt
    + WrappingAdd
{
    /// For all currently implemented BitVec types, the maximum possible distance
    /// can be stored in `u8`. Custom implementations using bigger integers can
    /// adjust `DistType` to hold bigger numbers. Note that due to how the traceback
    /// algorithm currently works, `DistType` should be able to represent numbers larger
    /// than the bit-width of the `BitVec` type. For instance, a hypothetical `BitVec` type
    /// of `u256` should use `u16` as distance, since `u8` cannot represent numbers larger
    /// than 256.
    type DistType: Copy
            + Default
            + AddAssign
            + SubAssign
            + PrimInt // includes Bounded, Num, Zero, One
            + FromPrimitive
            + Bounded;
}

macro_rules! impl_bitvec {
    ($type:ty, $dist:ty) => {
        impl BitVec for $type {
            type DistType = $dist;
        }
    };
}

impl_bitvec!(u8, u8);
impl_bitvec!(u16, u8);
impl_bitvec!(u32, u8);
impl_bitvec!(u64, u8);
#[cfg(has_u128)]
impl_bitvec!(u128, u8);

/// Myers algorithm.
pub struct Myers<T = u64>
where
    T: BitVec,
{
    peq: [T; 256],
    bound: T,
    m: T::DistType,
    tb: Traceback<T>,
}

impl<T: BitVec> Myers<T> {
    /// Create a new instance of Myers algorithm for a given pattern.
    pub fn new<C, P>(pattern: P) -> Self
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        let maxsize = T::DistType::from_usize(size_of::<T>() * 8).unwrap();
        let pattern = pattern.into_iter();
        let m = T::DistType::from_usize(pattern.len()).unwrap();
        assert!(m <= maxsize, "Pattern too long");
        assert!(m > T::DistType::zero(), "Pattern is empty");

        let mut peq = [T::zero(); 256];

        for (i, a) in pattern.enumerate() {
            peq[*a.borrow() as usize] |= T::one() << i;
        }

        Myers {
            peq,
            bound: T::one() << (m.to_usize().unwrap() - 1),
            m,
            tb: Traceback::new(),
        }
    }

    #[inline]
    fn step(&self, state: &mut State<T>, a: u8) {
        let eq = self.peq[a as usize];
        let xv = eq | state.mv;
        let xh = ((eq & state.pv).wrapping_add(&state.pv) ^ state.pv) | eq;

        let mut ph = state.mv | !(xh | state.pv);
        let mut mh = state.pv & xh;

        if ph & self.bound > T::zero() {
            state.dist += T::DistType::one();
        } else if mh & self.bound > T::zero() {
            state.dist -= T::DistType::one();
        }

        ph <<= 1;
        mh <<= 1;
        state.pv = mh | !(xv | ph);
        state.mv = ph & xv;
    }

    // Combining these two steps into one function seems beneficial for performance
    fn step_trace(&mut self, state: &mut State<T>, a: u8) {
        self.step(state, a);
        self.tb.add_state(state.clone());
    }

    /// Calculate the global distance of the pattern to the given text.
    pub fn distance<C, I>(&self, text: I) -> T::DistType
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
    {
        let mut state = State::init(self.m);
        let mut dist = T::DistType::max_value();
        for a in text {
            self.step(&mut state, *a.borrow());
            if state.dist < dist {
                dist = state.dist;
            }
        }
        dist
    }

    /// Finds all matches of pattern in the given text up to a given maximum distance.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn find_all_end<C, I>(
        &self,
        text: I,
        max_dist: T::DistType,
    ) -> Matches<'_, T, C, I::IntoIter>
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
    {
        Matches::new(self, text.into_iter(), max_dist)
    }

    /// Find the best match of the pattern in the given text.
    /// if multiple end positions have the same distance, the first is returned.
    pub fn find_best_end<C, I>(&self, text: I) -> (usize, T::DistType)
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
    {
        self.find_all_end(text, T::DistType::max_value())
            .min_by_key(|&(_, dist)| dist)
            .unwrap()
    }

    /// Finds all matches of pattern in the given text up to a given maximum distance.
    /// In contrast to `find_all_end`, matches are returned as an iterator over ranges
    /// of `(start, end, distance)`. Note that the end coordinate is not included in the
    /// range and thus and thus greater by one compared to the end index returned by
    /// `find_all_end()`.
    pub fn find_all<'a, C, I>(
        &'a mut self,
        text: I,
        max_dist: T::DistType,
    ) -> FullMatches<'a, T, C, I::IntoIter>
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
        I::IntoIter: ExactSizeIterator,
    {
        FullMatches::new(self, text.into_iter(), max_dist)
    }

    /// As `find_all_end`, this function returns an iterator over tuples of `(end, distance)`.
    /// Additionally, it keeps the data necessary for later obtaining the starting positions and/or
    /// the alignment path at *any* position that was already searched.
    pub fn find_all_lazy<'a, C, I>(
        &'a mut self,
        text: I,
        max_dist: T::DistType,
    ) -> LazyMatches<'a, T, C, I::IntoIter>
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
        I::IntoIter: ExactSizeIterator,
    {
        LazyMatches::new(self, text.into_iter(), max_dist)
    }
}

/// The current algorithm state.
#[derive(Clone, Debug, Default)]
struct State<T = u64>
where
    T: BitVec,
{
    pv: T,
    mv: T,
    dist: T::DistType,
}

impl<T> State<T>
where
    T: BitVec,
{
    /// Create new state, initiating it
    pub fn init(m: T::DistType) -> Self {
        State {
            pv: T::max_value(),
            mv: T::zero(),
            dist: m,
        }
    }
}

/// Iterator over pairs of end positions and distance of matches.
pub struct Matches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    myers: &'a Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    max_dist: T::DistType,
}

impl<'a, T, C, I> Matches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    fn new(myers: &'a Myers<T>, text: I, max_dist: T::DistType) -> Self {
        let state = State::init(myers.m);
        Matches {
            myers,
            state,
            text: text.enumerate(),
            max_dist,
        }
    }
}

impl<'a, T, C, I> Iterator for Matches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    type Item = (usize, T::DistType);

    fn next(&mut self) -> Option<(usize, T::DistType)> {
        for (i, a) in self.text.by_ref() {
            self.myers.step(&mut self.state, *a.borrow());
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }
}

/// Iterator over tuples of starting position, end position and distance of matches. In addition,
/// methods for obtaining the hit alignment path are provided.
pub struct FullMatches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    myers: &'a mut Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    text_len: usize,
    m: T::DistType,
    max_dist: T::DistType,
    pos: usize, // current end position, has to be stored for alignment() method
    finished: bool,
}

impl<'a, T, C, I> FullMatches<'a, T, C, I>
where
    T: 'a + BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C> + ExactSizeIterator,
{
    fn new(myers: &'a mut Myers<T>, text_iter: I, max_dist: T::DistType) -> Self {
        let state = State::init(myers.m);
        // Maximum number of traceback columns possibly used by a match
        let num_cols = (myers.m + min(max_dist, myers.m)).to_usize().unwrap();
        myers.tb.init(state.clone(), num_cols, myers.m);
        FullMatches {
            state,
            m: myers.m,
            myers,
            text_len: text_iter.len(),
            text: text_iter.enumerate(),
            max_dist,
            pos: 0,
            finished: false,
        }
    }

    /// Searches the next match and returns a tuple of end position and distance
    /// if found. This involves *no* searching for a starting position and is thus
    /// faster than just iterating over `FullMatches`
    #[inline]
    pub fn next_end(&mut self) -> Option<(usize, T::DistType)> {
        for (i, a) in self.text.by_ref() {
            self.pos = i; // used in alignment()
            self.myers.step_trace(&mut self.state, *a.borrow());
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        self.finished = true;
        None
    }

    /// Searches the next match and returns a tuple of starting position, end position and
    /// distance, or `None` if no match was found. In addition, the alignment path is added to
    /// `ops`.
    #[inline]
    pub fn next_path(
        &mut self,
        ops: &mut Vec<AlignmentOperation>,
    ) -> Option<(usize, usize, T::DistType)> {
        self.next_end()
            .map(|(end, dist)| (self.path(ops).unwrap(), end + 1, dist))
    }

    /// Searches the next match and updates the given `Alignment` with its position
    /// and alignment path if found. The distance is stored in `Alignment::score`.
    /// If no next hit is found, `false` is returned and `aln` remains unchanged.
    #[inline]
    pub fn next_alignment(&mut self, aln: &mut Alignment) -> bool {
        if self.next_end().is_some() {
            self.alignment(aln);
            return true;
        }
        false
    }

    /// Returns the starting position of the current hit. If the search is finished and no hit was
    /// found, `None` is returned.
    #[inline]
    pub fn start(&self) -> Option<usize> {
        if self.finished {
            return None;
        }
        let (len, _) = self.myers.tb.traceback(None);
        Some(self.pos + 1 - len)
    }

    /// Adds the path of the current hit alignment to `ops` and returns the starting position of
    /// the current hit. If the search is finished and no hit was found, `None` is returned.
    #[inline]
    pub fn path(&self, ops: &mut Vec<AlignmentOperation>) -> Option<usize> {
        if self.finished {
            return None;
        }
        let (len, _) = self.myers.tb.traceback(Some(ops));
        Some(self.pos + 1 - len)
    }

    /// Updates the given `Alignment` with its position and alignment path. The edit distance is
    /// stored in `Alignment::score`. If no hit has been found, then `false` will be returned
    /// and nothing is done.
    #[inline]
    pub fn alignment(&mut self, aln: &mut Alignment) -> bool {
        if self.finished {
            return false;
        }
        let (len, _) = self.myers.tb.traceback(Some(&mut aln.operations));
        update_aln(
            self.pos,
            len,
            self.text_len,
            self.state.dist.to_usize().unwrap(),
            self.m.to_usize().unwrap(),
            aln,
        );
        true
    }
}

impl<'a, T, C, I> Iterator for FullMatches<'a, T, C, I>
where
    T: 'a + BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C> + ExactSizeIterator,
{
    type Item = (usize, usize, T::DistType);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.next_end()
            .map(|(end, dist)| (self.start().unwrap(), end + 1, dist))
    }
}

/// Iterator over tuples of end position and distance of matches. In addition,
/// methods for obtaining the hit alignment path are provided.
pub struct LazyMatches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    myers: &'a mut Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    text_len: usize,
    m: T::DistType,
    max_dist: T::DistType,
}

impl<'a, T, C, I> Iterator for LazyMatches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C> + ExactSizeIterator,
{
    type Item = (usize, T::DistType);

    fn next(&mut self) -> Option<(usize, T::DistType)> {
        for (i, a) in self.text.by_ref() {
            self.myers.step_trace(&mut self.state, *a.borrow());
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }
}

impl<'a, T, C, I> LazyMatches<'a, T, C, I>
where
    T: 'a + BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C> + ExactSizeIterator,
{
    fn new(myers: &'a mut Myers<T>, text_iter: I, max_dist: T::DistType) -> Self {
        let state = State::init(myers.m);
        myers.tb.init(state.clone(), text_iter.len(), myers.m);
        LazyMatches {
            state,
            m: myers.m,
            myers,
            text_len: text_iter.len(),
            text: text_iter.enumerate(),
            max_dist,
        }
    }

    /// Takes the end position of a hit and returns a tuple of the corresponding starting position
    /// and the hit distance. If the end position is greater than the end position of the previously
    /// returned hit, `None` is returned.
    /// *Note:* A 0-based end position is expected (as returned by `next_end`).
    pub fn hit_at(&self, end_pos: usize) -> Option<(usize, T::DistType)> {
        self.myers
            .tb
            .traceback_at(end_pos, None)
            .map(|(len, dist)| (end_pos + 1 - len, dist))
    }

    /// Takes the end position of a hit and returns a tuple of the corresponding starting position
    /// and the hit distance. The alignment path is added to `ops`.
    /// As in `hit_at`, the end position has to be searched already, otherwise `None` is returned.
    pub fn path_at(
        &self,
        end_pos: usize,
        ops: &mut Vec<AlignmentOperation>,
    ) -> Option<(usize, T::DistType)> {
        self.myers
            .tb
            .traceback_at(end_pos, Some(ops))
            .map(|(len, dist)| (end_pos + 1 - len, dist))
    }

    /// Takes the end position of a hit and returns a tuple of the corresponding starting position
    /// and the hit distance. The alignment `aln` is updated with the position, alignment path
    /// and distance (stored in `Alignment::score`).
    /// If the end position has not yet been searched, nothing is done and `false` is returned.
    pub fn alignment_at(&self, end_pos: usize, aln: &mut Alignment) -> bool {
        if let Some((aln_len, dist)) = self
            .myers
            .tb
            .traceback_at(end_pos, Some(&mut aln.operations))
        {
            update_aln(
                end_pos,
                aln_len,
                self.text_len,
                dist.to_usize().unwrap(),
                self.m.to_usize().unwrap(),
                aln,
            );
            return true;
        }
        false
    }
}

/// Assumes *0-based* end positions, the coordinates will be converted to 1-based
#[inline(always)]
fn update_aln(
    end_pos: usize,
    aln_len: usize,
    text_len: usize,
    dist: usize,
    m: usize,
    aln: &mut Alignment,
) {
    aln.xstart = 0;
    aln.xend = m;
    aln.xlen = m;
    aln.ylen = text_len;
    aln.yend = end_pos + 1;
    aln.ystart = aln.yend - aln_len;
    aln.mode = AlignmentMode::Semiglobal;
    aln.score = dist as i32;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::AlignmentOperation::*;
    use crate::alignment::{Alignment, AlignmentMode};
    use itertools::Itertools;
    use std::iter::repeat;

    #[test]
    fn test_find_all_end() {
        let text = b"ACCGTGGATGAGCGCCATAG";
        let pattern = b"TGAGCGT";
        let myers = Myers::<u64>::new(pattern);
        let occ = myers.find_all_end(text, 1).collect_vec();
        assert_eq!(occ, [(13, 1), (14, 1)]);
    }

    #[test]
    fn test_find_best_end() {
        let text = b"ACCGTGGATGAGCGCCATAG";
        let pattern = b"TGAGCGT";
        let myers = Myers::<u64>::new(pattern);
        let hit = myers.find_best_end(text);
        assert_eq!(hit, (13, 1));
    }

    #[test]
    fn test_distance() {
        let text = b"TGAGCNTA";
        let patt = b"TGAGCGT";

        let myers = Myers::<u64>::new(patt);
        assert_eq!(myers.distance(text), 1);

        let myers_wildcard = MyersBuilder::new().text_wildcard(b'N').build_64(patt);
        assert_eq!(myers_wildcard.distance(text), 0);
    }

    #[test]
    fn test_full_position() {
        let text = b"CAGACATCTT";
        let pattern = b"AGA";

        let mut myers = Myers::<u64>::new(pattern);
        let matches: Vec<_> = myers.find_all(text, 1).collect();
        assert_eq!(&matches, &[(1, 3, 1), (1, 4, 0), (1, 5, 1), (3, 6, 1)]);
    }

    #[test]
    fn test_traceback_path() {
        let text = "TCAGACAT-CTT".replace('-', "").into_bytes();
        let patt = "TC-GACGTGCT".replace('-', "").into_bytes();

        let mut myers = Myers::<u64>::new(&patt);
        let mut matches = myers.find_all(&text, 3);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 10, 3));
        assert_eq!(
            aln,
            &[Match, Match, Del, Match, Match, Match, Subst, Match, Ins, Match, Match]
        );
    }

    #[test]
    fn test_traceback_path2() {
        let text = "TCAG--CAGATGGAGCTC".replace('-', "").into_bytes();
        let patt = "TCAGAGCAG".replace('-', "").into_bytes();

        let mut myers = Myers::<u64>::new(&patt);
        let mut matches = myers.find_all(&text, 2);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 7, 2));
        assert_eq!(
            aln,
            &[Match, Match, Match, Match, Ins, Ins, Match, Match, Match]
        );
    }

    #[test]
    #[rustfmt::skip]
    fn test_alignment() {
        let text =  "GGTCCTGAGGGATTA".replace('-', "");
        let pattern = "TCCT-AGGGA".replace('-', "");

        let mut myers = Myers::<u64>::new(pattern.as_bytes());
        let expected = Alignment {
            score: 1,
            xstart: 0,
            xend: 9,
            xlen: 9,
            ystart: 2,
            yend: 12,
            ylen: 15,
            operations: vec![
                Match, Match, Match, Match, Del, Match, Match, Match, Match, Match,
            ],
            mode: AlignmentMode::Semiglobal,
        };

        let mut aln = Alignment::default();
        {
            let mut matches = myers.find_all(text.as_bytes(), 1);
            assert!(matches.next_alignment(&mut aln));
            assert_eq!(&aln, &expected);

            aln.score = -1;
            matches.alignment(&mut aln);
            assert_eq!(&aln, &expected);
        }
        // Lazy API
        aln.score = -1;
        let end = expected.yend - 1;
        let mut lazy_matches = myers.find_all_lazy(text.as_bytes(), 1);
        assert!(!lazy_matches.alignment_at(end, &mut aln));
        // now search positions
        aln.score = -1;
        assert_eq!(lazy_matches.next(), Some((end, expected.score as u8)));
        assert!(lazy_matches.alignment_at(end, &mut aln));
        assert_eq!(&aln, &expected);
    }

    #[test]
    fn test_position_cmp() {
        // same as position_at, but 0-based positions from
        let text = b"CAGACATCTT";
        let pattern = b"AGA";
        let starts_exp = [1, 1, 1, 3];
        let end_dist_exp = [(2, 1), (3, 0), (4, 1), (5, 1)];

        let mut myers = Myers::<u64>::new(pattern);

        // standard iterator with 0-based ends
        let end_dist: Vec<_> = myers.find_all_end(text, 1).collect();
        assert_eq!(&end_dist, &end_dist_exp);

        // iterator over full ranges where ends are + 1
        let full_hits: Vec<_> = myers.find_all(text, 1).collect();

        // allows to retrive starting position later
        let mut lazy_matches = myers.find_all_lazy(text, 1);

        // compare with each other and lazy matches
        for ((&start, (end, dist)), (f_start, f_end, f_dist)) in
            starts_exp.into_iter().zip(end_dist).zip(full_hits)
        {
            assert_eq!(start, f_start);
            assert_eq!(dist, f_dist);
            assert_eq!(end + 1, f_end);

            // lazy API
            let (lazy_end, lazy_dist) = lazy_matches.next().unwrap();
            assert_eq!(end, lazy_end);
            assert_eq!(dist, lazy_dist);
            assert_eq!(lazy_matches.hit_at(end), Some((start, dist)));
            // For positions above, information is not (yet) available
            assert_eq!(lazy_matches.hit_at(end + 1), None);
        }
    }

    #[test]
    fn test_path_at() {
        let text = b"CAGACATCTT";
        let pattern = b"AGA";

        let mut myers = Myers::<u64>::new(pattern);
        let mut matches = myers.find_all_lazy(text, 1);

        let expected = &[Match, Match, Ins];
        let mut path = vec![];

        // search first hit
        assert_eq!(matches.next(), Some((2, 1)));

        // retrieve first hit at 0-based end position (2)
        assert_eq!(matches.hit_at(2), Some((1, 1)));
        assert_eq!(matches.path_at(2, &mut path), Some((1, 1)));
        assert_eq!(&path, expected);

        // hit out of range
        path.clear();
        assert!(matches.path_at(3, &mut path).is_none());
        assert!(path.is_empty());

        // now search the next hit
        assert_eq!(matches.next(), Some((3, 0)));
        // position 3 is now searched -> path can be retrieved
        assert_eq!(matches.path_at(3, &mut path), Some((1, 0)));
        assert_eq!(&path, &[Match, Match, Match])
    }

    #[test]
    fn test_shorter() {
        let text = "ATG";
        let pat = "CATGC";

        let mut myers = Myers::<u64>::new(pat.as_bytes());
        let mut matches = myers.find_all(text.as_bytes(), 2);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 3, 2));
        assert_eq!(aln, &[Ins, Match, Match, Match, Ins]);
    }

    #[test]
    fn test_long_shorter() {
        let text = "CCACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT";
        let pattern = "CGCGGTGTCCACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT";

        let mut myers = Myers::<u64>::new(pattern.as_bytes());
        let mut matches = myers.find_all(text.as_bytes(), 8);
        assert_eq!(matches.next().unwrap(), (0, 52, 8));
    }

    #[test]
    fn test_ambig() {
        let text = b"TGABCNTR";
        let patt = b"TGRRCGTR";
        //                x  x
        // Matching is asymmetric here (A matches R and G matches N, but the reverse is not true)

        let myers = MyersBuilder::new().ambig(b'R', b"AG").build_64(patt);
        assert_eq!(myers.distance(text), 2);
    }

    #[test]
    fn test_longest_possible() {
        let text = b"CCACGCGT";

        let mut myers: Myers<u8> = Myers::new(text);
        assert_eq!(myers.find_all(text, 0).next(), Some((0, 8, 0)));
    }

    #[test]
    fn test_large_dist() {
        let pattern: Vec<_> = repeat(b'T').take(64).collect();
        let text: Vec<_> = repeat(b'A').take(64).collect();

        let mut myers = Myers::<u64>::new(&pattern);
        let max_dist = myers
            .find_all_end(&text, 64)
            .max_by_key(|&(_, dist)| dist)
            .unwrap()
            .1;
        assert_eq!(max_dist, 64);

        let max_dist = myers
            .find_all(&text, 64)
            .max_by_key(|&(_, _, dist)| dist)
            .unwrap()
            .1;
        assert_eq!(max_dist, 64);
    }

    #[test]
    #[should_panic(expected = "Pattern too long")]
    fn test_pattern_too_long() {
        let pattern: Vec<_> = repeat(b'T').take(65).collect();
        Myers::<u64>::new(&pattern);
    }

    #[test]
    #[should_panic(expected = "Pattern too long")]
    fn test_pattern_too_long_builder() {
        let pattern: Vec<_> = repeat(b'T').take(65).collect();
        MyersBuilder::new().build_64(&pattern);
    }
}
