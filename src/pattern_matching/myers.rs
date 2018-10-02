// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Myers bit-parallel approximate pattern matching algorithm.
//! Finds all matches up to a given edit distance. The pattern has to fit into a bitvector,
//! and is thus limited to 64 or (since Rust version 1.26) to 128 symbols.
//! Complexity: O(n)
//!
//! # Example
//!
//! Iterating over matches in pairs of `(end, distance)` using `u64` as bitvector type:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers64;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let myers = Myers64::new(pattern);
//! let occ: Vec<_> = myers.find_all_end(text, 2).collect();
//!
//! assert_eq!(occ, [(11, 2), (12, 2)]);
//! # }
//! ```
//!
//! Starting with Rust 1.26, it is also possible to use `u128` as bitvector (`Myers128`), which
//! will is somewhat slower.
//!
//! # Obtaining the starting position of a match
//!
//! The `Myers::find_all` method provides an iterator over tuples of `(start, end, distance)`.
//! This is slower, since the alignment path has to be calculated if the starting position is
//! required. Note that the end positions differ from above by one. This is intentional,
//! representing a range, which does not includde the end position.
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers64;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers64::new(pattern);
//! let occ: Vec<_> = myers.find_all(text, 2).collect();
//!
//! assert_eq!(occ, [(3, 12, 2), (3, 13, 2)]);
//! # }
//! ```
//!
//! # Obtaining alignments
//!
//! [FullMatches](struct.FullMatches.html) returned by `Myers::find_all()` also allow obtaining
//! an alignment path:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::{Myers64, new_alignment};
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers64::new(pattern);
//! // create an 'empty' alignment instance, which can be reused
//! let mut aln = new_alignment(); // TODO: change
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
//! (rest omitted)
//!
//! </pre>
//!
//! **Note** that the [Alignment](../../alignment/struct.Alignment.html) instance is only created
//! once and then reused. Because the Myers algorithm is very fast, the allocation necessary for
//! `Alignment::operations` can have a non-negligible performance impact.
//!
//! # Finding the best hit
//!
//! In many cases, only the match with the smallest edit distance is actually of interest.
//! Calculating an alignment for every hit is therefore not necessary. `Myers::find_all_lazy()`
//! iterates over tuples of `(end, distance)` like `Myers::find_all_end()`, but additionally
//! keeps the data necessary for calculating the alignment path later at any desired position.
//! Storing the data itself has a slight performance impact and requires more memory
//! compared to `Myers::find_all_end()` [O(n + m)].
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::{Myers64, new_alignment};
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers64::new(pattern);
//! let mut aln = new_alignment(); // TODO: change
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
//!     .build(pattern);
//! assert_eq!(myers.distance(text), 0);
//! # }
//! ```
//!
//! For more examples see the documentation of [MyersBuilder](struct.MyersBuilder.html).

use std::borrow::Borrow;
use std::collections::HashMap;
use std::iter;
use std::mem::size_of;
use std::ops::Range;
use std::ops::*;
use std::u64;

use num_traits::{Bounded, FromPrimitive, One, PrimInt, ToPrimitive, WrappingAdd, Zero};

use alignment::{Alignment, AlignmentMode, AlignmentOperation};
use utils::{IntoTextIterator, TextIterator};

/// Integer types serving as bit vectors must implement this trait.
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
    /// adjust `DistType` to hold bigger numbers.
    type DistType: Copy
            + Default
            + AddAssign
            + SubAssign
            + PrimInt // includes Bounded, Num, Zero, One
            + FromPrimitive;
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

/// Myers instance using `u64` as bit vectors (pattern length up to 64)
pub type Myers64 = Myers<u64>;
/// Myers instance using `u128` as bit vectors (pattern length up to 128)
#[cfg(has_u128)]
pub type Myers128 = Myers<u128>;

/// Builds a Myers instance, allowing to specify ambiguities.
///
/// # Example:
///
/// This example shows how recognition of IUPAC ambiguities in patterns
/// can be implemented. Ambiguities in the searched text are not recognized; this
/// would require specifying additional ambiguities (A = MRWVHDN, etc...).
///
/// ```
/// # extern crate bio;
/// use bio::pattern_matching::myers::MyersBuilder;
///
/// # fn main() {
/// let ambigs = [
///     (b'M', &b"ACM"[..]),
///     (b'R', &b"AGR"[..]),
///     (b'W', &b"ATW"[..]),
///     (b'S', &b"CGS"[..]),
///     (b'Y', &b"CTY"[..]),
///     (b'K', &b"GTK"[..]),
///     (b'V', &b"ACGMRSV"[..]),
///     (b'H', &b"ACTMWYH"[..]),
///     (b'D', &b"AGTRWKD"[..]),
///     (b'B', &b"CGTSYKB"[..]),
///     (b'N', &b"ACGTMRWSYKVHDBN"[..])
/// ];
///
/// let mut builder = MyersBuilder::new();
///
/// for &(base, equivalents) in &ambigs {
///     builder.ambig(base, equivalents);
/// }
///
/// let text = b"GGATGNGCGCCATAG";
/// let pattern = b"TRANCGG";
/// //                *   * (mismatch)
///
/// let myers = builder.build(pattern);
/// assert_eq!(myers.distance(text), 2);
/// # }
/// ```
#[derive(Default, Clone, Eq, PartialEq)]
pub struct MyersBuilder {
    ambigs: HashMap<u8, Vec<u8>>,
    wildcards: Vec<u8>,
}

impl MyersBuilder {
    pub fn new() -> MyersBuilder {
        Self::default()
    }

    /// Allows to specify ambiguous characters and their equivalents.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::MyersBuilder;
    ///
    /// # fn main() {
    /// let text = b"ACCGTGGATGAGCGCCATAG";
    /// let pattern =      b"TGAGCGN";
    ///
    /// let myers = MyersBuilder::new()
    ///     .ambig(b'N', b"ACGT")
    ///     .build(pattern);
    ///
    /// assert_eq!(myers.distance(text), 0);
    /// # }
    pub fn ambig<I, B>(&mut self, byte: u8, equivalents: I) -> &mut Self
    where
        I: IntoIterator<Item = B>,
        B: Borrow<u8>,
    {
        let eq = equivalents.into_iter().map(|b| *b.borrow()).collect();
        self.ambigs.insert(byte, eq);
        self
    }

    /// Allows to specify a wildcard character, that upon appearance in the search text
    /// shall be matched by any character of the pattern. Multiple wildcards are possible.
    /// For the inverse, that is, wildcards in the pattern matching any character in search
    /// text, use `ambig(byte, 0..255)`.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::MyersBuilder;
    ///
    /// # fn main() {
    /// let text = b"ACCGTGGATGAGCG*CATAG";
    /// let pattern =      b"TGAGCGT";
    ///
    /// let myers = MyersBuilder::new()
    ///     .text_wildcard(b'*')
    ///     .build(pattern);
    ///
    /// assert_eq!(myers.distance(text), 0);
    /// # }
    pub fn text_wildcard(&mut self, wildcard: u8) -> &mut Self {
        self.wildcards.push(wildcard);
        self
    }

    /// Creates a Myers instance given a pattern, using `u64` as bit vector type
    pub fn build<'a, P>(&self, pattern: P) -> Myers<u64>
    where
        P: IntoTextIterator<'a>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build_other(pattern)
    }

    /// Creates a Myers instance given a pattern, using `u128` as bit vector type
    #[cfg(has_u128)]
    pub fn build128<'a, P>(&self, pattern: P) -> Myers<u128>
    where
        P: IntoTextIterator<'a>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build_other(pattern)
    }

    /// Creates a Myers instance given a pattern, using any desired type for bit vectors
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::{MyersBuilder, Myers};
    ///
    /// # fn main() {
    /// let myers: Myers<u32> = MyersBuilder::new()
    ///     .text_wildcard(b'*')
    ///     .build_other(b"TGAGCG*");
    /// // ...
    /// # }
    pub fn build_other<'a, T, P>(&self, pattern: P) -> Myers<T>
    where
        T: BitVec,
        P: IntoTextIterator<'a>,
        P::IntoIter: ExactSizeIterator,
    {
        let maxsize = T::DistType::from_usize(size_of::<T>() * 8).unwrap();
        let pattern = pattern.into_iter();
        let m = T::DistType::from_usize(pattern.len()).unwrap();
        assert!(m <= maxsize, "Pattern too long");
        assert!(m > T::DistType::zero(), "Pattern is empty");

        let mut peq = [T::zero(); 256];

        for (i, &a) in pattern.enumerate() {
            let mask = T::one() << i;
            // equivalent
            peq[a as usize] |= mask;
            // ambiguities
            if let Some(equivalents) = self.ambigs.get(&a) {
                for &eq in equivalents {
                    peq[eq as usize] |= mask;
                }
            }
        }

        for &w in &self.wildcards {
            peq[w as usize] = T::max_value();
        }

        Myers {
            peq: peq,
            bound: T::one() << (m.to_usize().unwrap() - 1),
            m: m,
            tb: Traceback::new(),
        }
    }
}

/// Myers algorithm.
pub struct Myers<T>
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
    pub fn new<'a, P>(pattern: P) -> Self
    where
        P: IntoTextIterator<'a>,
        P::IntoIter: ExactSizeIterator,
    {
        let maxsize = T::DistType::from_usize(size_of::<T>() * 8).unwrap();
        let pattern = pattern.into_iter();
        let m = T::DistType::from_usize(pattern.len()).unwrap();
        assert!(m <= maxsize, "Pattern too long");
        assert!(m > T::DistType::zero(), "Pattern is empty");

        let mut peq = [T::zero(); 256];

        for (i, &a) in pattern.enumerate() {
            peq[a as usize] |= T::one() << i;
        }

        Myers {
            peq: peq,
            bound: T::one() << (m.to_usize().unwrap() - 1),
            m: m,
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
    pub fn distance<'a, I: IntoTextIterator<'a>>(&self, text: I) -> T::DistType {
        let mut state = State::init(self.m);
        let mut dist = T::DistType::max_value();
        for &a in text {
            self.step(&mut state, a);
            if state.dist < dist {
                dist = state.dist;
            }
        }
        dist
    }

    /// Finds all matches of pattern in the given text up to a given maximum distance.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn find_all_end<'a, I: IntoTextIterator<'a>>(
        &'a self,
        text: I,
        max_dist: T::DistType,
    ) -> Matches<T, I::IntoIter> {
        Matches::new(self, text.into_iter(), max_dist)
    }

    /// Finds all matches of pattern in the given text up to a given maximum distance.
    /// In contrast to `find_all_end`, matches are returned as an iterator over ranges
    /// of `(start, end, distance)`. Note that the end coordinate is not included in the
    /// range and thus and thus greater by one compared to the end index returned by
    /// `find_all_end()`.
    pub fn find_all<'a, I>(
        &'a mut self,
        text: I,
        max_dist: T::DistType,
    ) -> FullMatches<'a, T, I::IntoIter>
    where
        I: IntoTextIterator<'a>,
        I::IntoIter: ExactSizeIterator,
    {
        FullMatches::new(self, text.into_iter(), max_dist)
    }

    /// As `find_all_end`, this function returns an iterator over tuples of `(end, distance)`.
    /// Additionally, it keeps the data necessary for later obtaining the starting positions and/or
    /// the alignment path at *any* position that was already searched.
    pub fn find_all_lazy<'a, I>(
        &'a mut self,
        text: I,
        max_dist: T::DistType,
    ) -> LazyMatches<'a, T, I::IntoIter>
    where
        I: IntoTextIterator<'a>,
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
        let maxsize = T::DistType::from_usize(8 * size_of::<T>()).unwrap();
        let pv = if m == maxsize {
            T::max_value()
        } else {
            (T::one() << m.to_usize().unwrap()) - T::one()
        };

        State {
            pv,
            mv: T::zero(),
            dist: m,
        }
    }
}

/// Iterator over pairs of end positions and distance of matches.
pub struct Matches<'a, T, I>
where
    T: 'a + BitVec,
    I: TextIterator<'a>,
{
    myers: &'a Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    max_dist: T::DistType,
}

impl<'a, T, I> Matches<'a, T, I>
where
    T: BitVec,
    I: TextIterator<'a>,
{
    fn new(myers: &'a Myers<T>, text: I, max_dist: T::DistType) -> Self {
        let state = State::init(myers.m);
        Matches {
            myers: myers,
            state: state,
            text: text.enumerate(),
            max_dist: max_dist,
        }
    }
}

impl<'a, T, I> Iterator for Matches<'a, T, I>
where
    T: BitVec,
    I: TextIterator<'a>,
{
    type Item = (usize, T::DistType);

    fn next(&mut self) -> Option<(usize, T::DistType)> {
        for (i, &a) in self.text.by_ref() {
            self.myers.step(&mut self.state, a);
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }
}

/// Iterator over tuples of starting position, end position and distance of matches. In addition,
/// methods for obtaining the hit alignment path are provided.
pub struct FullMatches<'a, T, I>
where
    T: 'a + BitVec,
    I: TextIterator<'a>,
{
    myers: &'a mut Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    text_len: usize,
    m: T::DistType,
    max_dist: T::DistType,
    pos: usize, // current end position, has to be stored for alignment() method
}

impl<'a, T, I> FullMatches<'a, T, I>
where
    T: 'a + BitVec,
    I: TextIterator<'a> + ExactSizeIterator,
{
    fn new(myers: &'a mut Myers<T>, text_iter: I, max_dist: T::DistType) -> Self {
        let state = State::init(myers.m);
        myers.tb.init(
            state.clone(),
            (myers.m + max_dist).to_usize().unwrap(),
            myers.m,
        );
        FullMatches {
            state: state,
            m: myers.m,
            myers: myers,
            text_len: text_iter.len(),
            text: text_iter.enumerate(),
            max_dist: max_dist,
            pos: 0,
        }
    }
    /// Searches the next match and returns a tuple of end position and distance
    /// if found. This involves *no* searching for a starting position and is thus
    /// faster than just iterating over `FullMatches`
    #[inline]
    pub fn next_end(&mut self) -> Option<(usize, T::DistType)> {
        for (i, &a) in self.text.by_ref() {
            self.pos = i; // used in alignment()
            self.myers.step_trace(&mut self.state, a);
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }

    /// Searches the next match and returns a tuple of starting position, end position and distance.
    /// if found. In addition, the alignment path is added to `ops`
    #[inline]
    pub fn next_path(
        &mut self,
        ops: &mut Vec<AlignmentOperation>,
    ) -> Option<(usize, usize, T::DistType)> {
        self.next_end()
            .map(|(end, dist)| (self.path(ops), end + 1, dist))
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

    /// Returns the starting position of the current hit.
    /// There are two corner cases with maybe 'unexpected' return values:
    /// - if searching was not yet started, `0` will be returned.
    /// - if the search is completed and no hit was found, the starting position
    ///   of the hit at the last position of the text will be returned, regardless
    ///   of how many mismatches there are.
    #[inline]
    pub fn start(&self) -> usize {
        let (len, _) = self.myers.tb.traceback(None);
        self.pos + 1 - len
    }

    /// Adds the path of the current hit alignment to `ops`
    /// and returns the starting position of the current hit.
    /// See `start()` for a description of corner cases.
    #[inline]
    pub fn path(&self, ops: &mut Vec<AlignmentOperation>) -> usize {
        let (len, _) = self.myers.tb.traceback(Some(ops));
        self.pos + 1 - len
    }

    /// Updates the given `Alignment` with its position
    /// and alignment path if found. The distance is stored in `Alignment::score`.
    /// See `start()` for a description of corner cases.
    #[inline]
    pub fn alignment(&mut self, aln: &mut Alignment) {
        let (len, _) = self.myers.tb.traceback(Some(&mut aln.operations));
        update_aln(
            self.pos,
            len,
            self.text_len,
            self.state.dist.to_usize().unwrap(),
            self.m.to_usize().unwrap(),
            aln,
        );
    }
}

impl<'a, T, I> Iterator for FullMatches<'a, T, I>
where
    T: 'a + BitVec,
    I: TextIterator<'a> + ExactSizeIterator,
{
    type Item = (usize, usize, T::DistType);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.next_end()
            .map(|(end, dist)| (self.start(), end + 1, dist))
    }
}

/// Iterator over tuples of end position and distance of matches. In addition,
/// methods for obtaining the hit alignment path are provided.
pub struct LazyMatches<'a, T, I>
where
    T: 'a + BitVec,
    I: TextIterator<'a>,
{
    myers: &'a mut Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    text_len: usize,
    m: T::DistType,
    max_dist: T::DistType,
}

impl<'a, T, I> Iterator for LazyMatches<'a, T, I>
where
    T: BitVec,
    I: TextIterator<'a> + ExactSizeIterator,
{
    type Item = (usize, T::DistType);

    fn next(&mut self) -> Option<(usize, T::DistType)> {
        for (i, &a) in self.text.by_ref() {
            self.myers.step_trace(&mut self.state, a);
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }
}

impl<'a, T, I> LazyMatches<'a, T, I>
where
    T: 'a + BitVec,
    I: TextIterator<'a> + ExactSizeIterator,
{
    fn new(myers: &'a mut Myers<T>, text_iter: I, max_dist: T::DistType) -> Self {
        let state = State::init(myers.m);
        myers.tb.init(
            state.clone(),
            text_iter.len() + max_dist.to_usize().unwrap(),
            myers.m,
        );
        LazyMatches {
            state: state,
            m: myers.m,
            myers: myers,
            text_len: text_iter.len(),
            text: text_iter.enumerate(),
            max_dist: max_dist,
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

// Assumes *0-based* end positions, the coordinates will be converted to 1-based
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

struct Traceback<T: BitVec> {
    states: Vec<State<T>>,
    positions: iter::Cycle<Range<usize>>,
    pos: usize,
    m: T::DistType,
}

impl<T> Traceback<T>
where
    T: BitVec,
{
    fn new() -> Traceback<T> {
        Traceback {
            states: vec![],
            positions: (0..0).cycle(),
            m: T::DistType::zero(),
            pos: 0,
        }
    }

    fn init(&mut self, initial_state: State<T>, num_cols: usize, m: T::DistType) {
        // ensure that empty text does not cause panic
        let num_cols = num_cols + 2;

        self.positions = (0..num_cols).cycle();

        // extend or truncate states vector
        let curr_len = self.states.len();
        if num_cols > curr_len {
            self.states.reserve(num_cols);
            self.states
                .extend((0..num_cols - curr_len).map(|_| State::default()));
        } else {
            self.states.truncate(num_cols);
        }
        // important if using unsafe in add_state(), and also for correct functioning of traceback
        debug_assert!(self.states.len() == num_cols);

        // first column is used to ensure a correct path if the text
        // is shorter than the pattern
        // Actual leftmost column starts at second position
        self.add_state(State {
            dist: T::DistType::max_value(),
            pv: T::max_value(), // all 1s
            mv: T::zero(),
        });

        // initial state
        self.add_state(initial_state);

        self.m = m;
    }

    #[inline]
    fn add_state(&mut self, s: State<T>) {
        self.pos = self.positions.next().unwrap();
        //self.states[self.pos] = s;
        // faster and safe since the positions are cycled.
        // Unless the code in init() is wrong, the index should
        // never be out of bounds.
        *unsafe { self.states.get_unchecked_mut(self.pos) } = s;
    }

    // Returns the length of the current match, optionally adding the
    // alignment path to `ops`
    fn traceback(&self, ops: Option<&mut Vec<AlignmentOperation>>) -> (usize, T::DistType) {
        self._traceback_at(self.pos, ops)
    }

    // Returns the length of a match with a given end position, optionally adding the
    // alignment path to `ops`
    // only to be called if the `states` vec contains all states of the text
    fn traceback_at(
        &self,
        pos: usize,
        ops: Option<&mut Vec<AlignmentOperation>>,
    ) -> Option<(usize, T::DistType)> {
        let pos = pos + 2; // in order to be comparable since self.pos starts at 2, not 0
        if pos <= self.pos {
            return Some(self._traceback_at(pos, ops));
        }
        None
    }

    // returns a tuple of alignment length and hit distance, optionally adding the alignment path
    // to `ops`
    #[inline]
    fn _traceback_at(
        &self,
        pos: usize,
        mut ops: Option<&mut Vec<AlignmentOperation>>,
    ) -> (usize, T::DistType) {
        use self::AlignmentOperation::*;

        // Reverse iterator over states. If remembering all positions,
        // the chain() and cycle() are not actually needed, but there seems
        // to be almost no performance loss.
        let mut states = self.states[..pos + 1]
            .iter()
            .rev()
            .chain(self.states.iter().rev().cycle());
        // // Simpler alternative using skip() is slower in some cases:
        // let mut states = self.states.iter().rev().cycle().skip(self.states.len() - pos - 1);

        // self.print_tb_matrix(pos);

        let ops = &mut ops;
        if let Some(o) = ops.as_mut() {
            o.clear();
        }

        // Mask for testing current position. The mask is moved
        // along mv / pv for testing all positions
        let mut current_pos = T::one() << (self.m.to_usize().unwrap() - 1);

        // horizontal column offset from starting point of traceback matrix (bottom right)
        let mut h_offset = 0;
        // vertical column offset from starting point of traceback matrix (bottom right)
        let mut v_offset = 0;

        // Macro for moving up one cell in traceback matrix, adjusting the distance of the state
        // of the state. Expects a bit mask indicating the current row position in the matrix.
        macro_rules! move_state_up {
            ($state:expr, $pos_mask:expr) => {
                if $state.pv & $pos_mask != T::zero() {
                    $state.dist -= T::DistType::one()
                } else if $state.mv & $pos_mask != T::zero() {
                    $state.dist += T::DistType::one()
                }
            };
        }

        // Macro for moving up 'n' cells in traceback matrix, adjusting the distance of the state.
        macro_rules! move_up_many {
            ($state:expr, $n:expr) => {
                // construct mask covering bits to check
                let m = self.m.to_usize().unwrap();
                let mask = if $n.to_usize().unwrap() == size_of::<T>() * 8 {
                    T::max_value()
                } else {
                    (!(T::max_value() << $n)) << (m - $n)
                };
                // adjust distance
                let p = ($state.pv & mask).count_ones() as i32;
                let m = ($state.mv & mask).count_ones() as i32;
                $state.dist = T::DistType::from_i32($state.dist.to_i32().unwrap() - p + m).unwrap();

                // // A loop seems always slower (not sure about systems without popcnt):
                // let mut p =  T::one() << (self.m.to_usize().unwrap() - 1);
                // for _ in 0..$n {
                //     move_state_up!($state, p);
                //     p >>= 1;
                // }
            };
        }

        // current state
        let mut state = states.next().unwrap().clone();
        // distance of the match
        let dist = state.dist;
        // state left to the current state
        let mut lstate = states.next().unwrap().clone();
        // The distance of the left state is always for diagonal of to the current position
        // in the traceback matrix. This allows checking for a substitution by a simple
        // comparison.
        move_state_up!(lstate, current_pos);

        // Macro for collectively moving up the cursor (adjusting the distance) of both
        // the left and the current state in the traceback matrix
        macro_rules! move_up {
            () => {
                v_offset += 1;
                move_state_up!(state, current_pos);
                current_pos >>= 1;
                move_state_up!(lstate, current_pos);
            };
        }

        // Macro for moving to the left, adjusting h_offset
        macro_rules! move_left {
            () => {
                h_offset += 1;
                state = lstate;
                lstate = states.next().unwrap().clone();
                if v_offset != self.m.to_usize().unwrap() {
                    // move up v_offset + 1 (not only v_offset) because left state
                    // is always in diagonal position
                    move_up_many!(lstate, v_offset + 1);
                }
            };
        }

        // Loop for finding the traceback path
        while v_offset < self.m.to_usize().unwrap() {
            let op = if lstate.dist + T::DistType::one() == state.dist {
                // Diagonal (substitution)
                // Since cursor of left state is always in the diagonal position,
                // this is a simple comparison of distances.
                v_offset += 1;
                current_pos >>= 1;
                move_left!();
                Subst
            } else if state.pv & current_pos != T::zero() {
                // Up
                move_up!();
                Ins
            } else if lstate.mv & current_pos != T::zero() {
                // Left
                move_left!();
                state.dist -= T::DistType::one();
                Del
            } else {
                // Diagonal (match)
                v_offset += 1;
                current_pos >>= 1;
                move_left!();
                Match
            };

            if let Some(o) = ops.as_mut() {
                o.push(op);
            }
        }

        if let Some(o) = ops.as_mut() {
            o.reverse();
        }

        (h_offset, dist)
    }

    // Useful for debugging
    #[allow(dead_code)]
    fn print_tb_matrix(&self, pos: usize) {
        let states = self.states[..pos + 1]
            .iter()
            .rev()
            .chain(self.states.iter().rev().cycle());

        let m = self.m.to_usize().unwrap();

        let mut out: Vec<_> = (0..m + 1).map(|_| vec![]).collect();
        for mut s in states.into_iter().cloned() {
            let mut current_pos = T::one() << (m - 1);
            let mut dist = s.dist;
            for i in 0..m + 1 {
                out[i].push(dist.to_usize().unwrap());
                if s.pv & current_pos != T::zero() {
                    dist -= T::DistType::one();
                } else if s.mv & current_pos != T::zero() {
                    dist += T::DistType::one();
                }
                current_pos >>= 1;
            }
            if s.dist == T::DistType::max_value() {
                break;
            }
        }
        out.reverse();
        for mut line in out {
            line.reverse();
            for dist in line {
                print!("{:<4}", dist);
            }
            print!("\n");
        }
    }
}

// temporary impl to create new 'empty' alignment
// This may be added to bio_types in some form
pub fn new_alignment() -> Alignment {
    Alignment {
        score: 0,
        ystart: 0,
        xstart: 0,
        yend: 0,
        xend: 0,
        ylen: 0,
        xlen: 0,
        operations: vec![],
        mode: AlignmentMode::Custom,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alignment::AlignmentOperation::*;
    use alignment::{Alignment, AlignmentMode};
    use itertools::Itertools;
    use std::iter::repeat;

    #[test]
    fn test_find_all_end() {
        let text = b"ACCGTGGATGAGCGCCATAG";
        let pattern = b"TGAGCGT";
        let myers = Myers64::new(pattern);
        let occ = myers.find_all_end(text, 1).collect_vec();
        assert_eq!(occ, [(13, 1), (14, 1)]);
    }

    #[test]
    fn test_distance() {
        let text = b"TGAGCNTA";
        let patt = b"TGAGCGT";

        let myers = Myers64::new(patt);
        assert_eq!(myers.distance(text), 1);

        let myers_wildcard = MyersBuilder::new().text_wildcard(b'N').build(patt);
        assert_eq!(myers_wildcard.distance(text), 0);
    }

    #[test]
    fn test_full_position() {
        let text = b"CAGACATCTT";
        let pattern = b"AGA";

        let mut myers = Myers64::new(pattern);
        let matches: Vec<_> = myers.find_all(text, 1).collect();
        assert_eq!(&matches, &[(1, 3, 1), (1, 4, 0), (1, 5, 1), (3, 6, 1)]);
    }

    #[test]
    fn test_traceback_path() {
        let text = "TCAGACAT-CTT".replace('-', "").into_bytes();
        let patt = "TC-GACGTGCT".replace('-', "").into_bytes();

        let mut myers = Myers64::new(&patt);
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

        let mut myers = Myers64::new(&patt);
        let mut matches = myers.find_all(&text, 2);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 7, 2));
        assert_eq!(
            aln,
            &[Match, Match, Match, Match, Ins, Ins, Match, Match, Match]
        );
    }

    #[test]
    #[cfg_attr(rustfmt, rustfmt_skip)]
    fn test_alignment() {
        let text =  "GGTCCTGAGGGATTA".replace('-', "");
        let pattern = "TCCT-AGGGA".replace('-', "");

        let mut myers = Myers64::new(pattern.as_bytes());
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
        let mut aln = new_alignment();

        {
            // TODO: a constructor for bio_types::alignment::Alignment
            // would be very convenient
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

        let mut myers = Myers64::new(pattern);

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

        let mut myers = Myers64::new(pattern);
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

        let mut myers = Myers64::new(pat.as_bytes());
        let mut matches = myers.find_all(text.as_bytes(), 2);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 3, 2));
        assert_eq!(aln, &[Ins, Match, Match, Match, Ins]);
    }

    #[test]
    fn test_long_shorter() {
        let text = "CCACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT";
        let pattern = "CGCGGTGTCCACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT";

        let mut myers = Myers64::new(pattern.as_bytes());
        let mut matches = myers.find_all(text.as_bytes(), 8);
        assert_eq!(matches.next().unwrap(), (0, 52, 8));
    }

    #[test]
    fn test_ambig() {
        let text = b"TGABCNT";
        let patt = b"TGRRCGT";
        //                x  x
        // Matching is asymmetric here (A matches R and G matches N, but the reverse is not true)

        let myers = MyersBuilder::new().ambig(b'R', b"AG").build(patt);
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

        let mut myers = Myers64::new(&pattern);
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
}
