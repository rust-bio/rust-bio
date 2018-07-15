// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Myers bit-parallel approximate pattern matching algorithm.
//! Finds all matches up to a given edit distance. The pattern has to fit into a bitvector,
//! and is here limited to 64 symbols.
//! Complexity: O(n)
//!
//! # Example
//!
//! ```
//! # extern crate itertools;
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers64;
//! use itertools::Itertools;
//!
//! # fn main() {
//! let text = b"ACCGTGGATGAGCGCCATAG";
//! let pattern = b"TGAGCGT";
//!
//! let myers = Myers64::new(pattern);
//! let occ = myers.find_all_end(text, 1).collect_vec();
//!
//! assert_eq!(occ, [(13, 1), (14, 1)]);
//! # }
//! ```

use std::iter;
use std::ops::Range;
use std::u64;
use std::borrow::Borrow;
use std::mem::size_of;
use std::ops::*;
use std::marker::PhantomData;

use num_traits::*;

use utils::{IntoTextIterator, TextIterator, TextSlice};
use alignment::{Alignment, AlignmentOperation, AlignmentMode};



/// Integer types serving as bit vectors must implement this trait.
pub trait Num: Add<Output=Self> + Sub<Output=Self> +
    BitOr<Output=Self> + BitOrAssign + BitAnd<Output=Self> + BitXor<Output=Self> + Not<Output=Self> +
    Shl<usize, Output=Self> + ShlAssign<usize> + ShrAssign<usize> +
    Copy + Clone + PartialOrd
{
    fn max_value() -> Self;
    fn zero() -> Self;
    fn one() -> Self;
    fn count_ones(self) -> u32;
    fn wrapping_add(&self, other: &Self) -> Self;
}

// Implementation using num_traits types
impl<T> Num for T
where T: Add<Output=Self> + Sub<Output=Self> +
    BitOr<Output=Self> + BitOrAssign + BitAnd<Output=Self> + BitXor<Output=Self> + Not<Output=Self> +
    Shl<usize, Output=Self> + ShlAssign<usize> + ShrAssign<usize> +
    Copy + PartialOrd +
    Bounded + PrimInt + Zero + One + WrappingAdd {
    fn max_value() -> T {
        <Self as Bounded>::max_value()
    }
    fn zero() -> Self {
        <Self as Zero>::zero()
    }
    fn one() -> Self {
        <Self as One>::one()
    }
    fn count_ones(self) -> u32 {
        <Self as PrimInt>::count_ones(self)
    }
    fn wrapping_add(&self, other: &T) -> T {
        <Self as WrappingAdd>::wrapping_add(&self, other)
    }
}

/// Myers instance based on 32-bit integers (pattern length up to 32)
pub type Myers32 = Myers<u32>;
/// Myers instance based on 64-bit integers (pattern length up to 64)
pub type Myers64 = Myers<u64>;
/// Myers instance based on 128-bit integers (pattern length up to 128)
#[cfg(has_u128)]
pub type Myers128 = Myers<u128>;


/// Myers algorithm.
pub struct Myers<T = u64>
where T: Num
{
    peq: [T; 256],
    bound: T,
    m: usize,
    tb: Traceback<T>,
}



impl<T: Num> Myers<T> {
    /// Create a new instance of Myers algorithm for a given pattern.
    pub fn new<'a, P: IntoTextIterator<'a>>(pattern: P) -> Self {
        Self::with_ambig(pattern.into_iter().map(Some))
    }

    // TODO: eventually move to module docs and use shorter example
    /// Like `Myers64::new()`, but additionally allows for specifying
    /// multiple matching characters in a pattern.
    /// As in `new`, an iterator over pattern positions is supplied. However, a nested iterable
    /// additionally allows to return multiple characters which should be recognized
    /// as matching at a given position. Unambiguous positions can be specified as iterables
    /// with a single element (see below). Even positions that don't match any character are
    /// possible by specifying an empty iterable.
    ///
    /// *Note*: Only a single indirection is possible, thus patterns with 'double borrows' like
    ///  this one: `&[&b"abc"[..], &b"de"[..]]` are not possible. The outer slice could be replaced
    /// by a `vec!` to make it work.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::Myers64;
    ///
    /// # fn main() {
    /// let text =    b"TGACGNTGA";
    /// let pattern = b"TGANGCTGA";
    ///
    /// // 'N' has no special meaning:
    /// let myers = Myers64::new(pattern);
    /// assert_eq!(myers.distance(text), 2);
    ///
    /// // 'N' in a pattern matches all four bases, but 'N' in
    /// // the text does not match any (asymmetric matching):
    /// let myers_ambig_asymm = Myers64::with_ambig(pattern.into_iter().map(|&b| {
    ///     // replacing N with all possible bases.
    ///     // Using vectors. Otherwise, the ref_slice crate could be used
    ///     if b == b'N' { b"ACGT".to_vec() } else { vec![b] }
    /// }));
    /// assert_eq!(myers_ambig_asymm.distance(text), 1);
    ///
    /// // 'N' matches both ways:
    /// let myers_ambig_symm = Myers64::with_ambig(pattern.into_iter().map(|&b| {
    ///     if b == b'N' { b"ACGT".to_vec() } else { vec![b, b'N'] }
    /// }));
    /// assert_eq!(myers_ambig_symm.distance(text), 0);
    /// # }
    pub fn with_ambig<P, I, C>(pattern: P) -> Self
        where P: IntoIterator<Item=I>,
              I: IntoIterator<Item=C>,
              C: Borrow<u8>,
    {
        Self::with_ambig_results::<_, I, C, ()>(
            pattern.into_iter().map(Ok)
        ).unwrap()
    }

    // TODO: eventually move to module docs and use shorter example
    /// Like `Myers64::with_ambig()`, but allows returning a `Result` during iteration
    /// of the pattern. The following example enables matching all ambiguous characters according
    /// to the IUPAC nomenclature.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::Myers64;
    ///
    /// # fn main() {
    ///
    /// use std::collections::HashMap;
    ///
    /// let text =    b"TGACGNTGA";
    /// let pattern = b"TGANGCTGR";
    ///
    /// // HashMap containing all character definitions.
    /// // Note: using the 'maplit' crate could make this less wordy.
    /// let ambig_map: HashMap<_, _> = [
    ///     (b'A', &b"A"[..]),
    ///     (b'C', &b"C"[..]),
    ///     (b'G', &b"G"[..]),
    ///     (b'T', &b"T"[..]),
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
    /// ].into_iter()
    ///  .map(|(b1, b2)| (*b1, b2.to_vec()))
    ///  .collect();
    ///
    /// let myers = Myers64::with_ambig_results(pattern.into_iter().map(|b| {
    ///     ambig_map.get(b).ok_or(format!("Invalid base found: {}", *b as char))
    /// })).unwrap();
    ///
    /// assert_eq!(myers.distance(text), 1);
    /// # }
    /// ```
    pub fn with_ambig_results<P, I, C, E>(pattern: P) -> Result<Self, E>
        where P: IntoIterator<Item=Result<I, E>>,
              I: IntoIterator<Item=C>,
              C: Borrow<u8>,
    {
        let mut peq = [T::zero(); 256];
        let w = size_of::<T>() * 8;
        let mut m = 0;
        for (i, var) in pattern.into_iter().enumerate() {
            m += 1;
            assert!(m <= w, "Pattern too long");
            for a in var?.into_iter() {
                peq[*a.borrow() as usize] |= T::one() << i;
            }
        }

        assert!(m > 0, "Pattern is empty");

        Ok(Myers {
            peq: peq,
            bound: T::one() << (m - 1),
            m: m,
            tb: Traceback::new(),
        })
    }

    /// Create a new instance of Myers algorithm for a given pattern and a wildcard character
    /// that shall match any character.
    pub fn with_wildcard(pattern: TextSlice, wildcard: u8) -> Self {
        let mut myers = Self::new(pattern);
        // wildcard matches all symbols of the pattern.
        myers.peq[wildcard as usize] = T::max_value();

        myers
    }

    #[inline]
    fn step(&self, state: &mut State<T>, a: u8) {
        let eq = self.peq[a as usize];
        let xv = eq | state.mv;
        let xh = ((eq & state.pv).wrapping_add(&state.pv) ^ state.pv) | eq;

        let mut ph = state.mv | !(xh | state.pv);
        let mut mh = state.pv & xh;

        if ph & self.bound > T::zero() {
            state.dist += 1;
        } else if mh & self.bound > T::zero() {
            state.dist -= 1;
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
    pub fn distance<'a, I: IntoTextIterator<'a>>(&self, text: I) -> usize {
        let mut state = State::new(self.m);
        for &a in text {
            self.step(&mut state, a);
        }
        state.dist
    }

    /// Find all matches of pattern in the given text up to a given maximum distance.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn find_all_end<'a, I: IntoTextIterator<'a>>(
        &'a self,
        text: I,
        max_dist: usize
    )
    -> Matches<T, I::IntoIter>
    {
        let state = State::new(self.m);
        Matches {
            myers: self,
            state: state,
            text: text.into_iter().enumerate(),
            max_dist: max_dist,
        }
    }

    /// Find all matches of pattern in the given text up to a given maximum distance.
    /// In contrast to `find_all_end`, matches are returned as an iterator over ranges
    /// of (start, end, distance). Note that the end coordinate is a range coordinate,
    /// not included in the range (end index + 1) and is thus not equivalent to the end
    /// position returned by `find_all_end()`. In order to obtain an alignment, use
    ///
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::Myers64;
    /// use bio::alignment::AlignmentOperation::*;
    ///
    /// # fn main() {
    /// let text = b"ACCGTGGATGAGCGCCATAG";
    /// let pattern =      b"TGAGCGT";
    ///
    /// // only range coordinates required
    /// let mut myers = Myers64::new(pattern);
    /// let occ: Vec<_> = myers.find_all_pos(text, 1).collect();
    /// assert_eq!(occ, [(8, 14, 1), (8, 15, 1)]);
    ///
    /// let mut myers = Myers64::new(pattern);
    /// let mut aln = vec![];
    /// let mut matches = myers.find_all_pos(text, 1);
    /// assert_eq!(matches.next_path(&mut aln).unwrap(), (8, 14, 1));
    /// assert_eq!(aln, &[Match, Match, Match, Match, Match, Match, Ins]);
    /// assert_eq!(matches.next_path(&mut aln).unwrap(), (8, 15, 1));
    /// assert_eq!(aln, &[Match, Match, Match, Match, Match, Match, Subst]);
    /// # }
    pub fn find_all_pos<'a, I: IntoTextIterator<'a>>(
        &'a mut self,
         text: I,
         max_dist: usize
    )
    -> FullMatches<'a, T, I::IntoIter, NoRemember>
    {
        self.tb.init(self.m + max_dist, self.m);
        FullMatches {
            state: State::new(self.m),
            m: self.m,
            myers: self,
            text: text.into_iter().enumerate(),
            max_dist: max_dist,
            pos: 0,
            _remember: PhantomData,
        }
    }

    /// Like `find_all_pos`, but additionally allows for obtaining the starting positions and/or
    /// the alignment at *any* position that was already searched.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::Myers64;
    /// use bio::alignment::AlignmentOperation::*;
    ///
    /// # fn main() {
    /// let text = b"ACCGTGGATGAGCGCCATAG";
    /// let pattern =      b"TGAGCGT";
    ///
    /// // only range coordinates required
    /// let mut myers = Myers64::new(pattern);
    /// let occ: Vec<_> = myers.find_all_pos(text, 1).collect();
    /// assert_eq!(occ, [(8, 14, 1), (8, 15, 1)]);
    ///
    /// let mut myers = Myers64::new(pattern);
    /// let mut aln = vec![];
    /// let mut matches = myers.find_all_pos(text, 1);
    /// assert_eq!(matches.next_path(&mut aln).unwrap(), (8, 14, 1));
    /// assert_eq!(aln, &[Match, Match, Match, Match, Match, Match, Ins]);
    /// assert_eq!(matches.next_path(&mut aln).unwrap(), (8, 15, 1));
    /// assert_eq!(aln, &[Match, Match, Match, Match, Match, Match, Subst]);
    /// # }
    pub fn find_all_pos_remember<'a, I>(
        &'a mut self,
         text: I,
         max_dist: usize
    )
    -> FullMatches<'a, T, I::IntoIter, Remember>
    where I: IntoTextIterator<'a>,
          I::IntoIter: ExactSizeIterator,
    {
        let text_iter = text.into_iter();
        self.tb.init(text_iter.len() + max_dist, self.m);
        FullMatches {
            state: State::new(self.m),
            m: self.m,
            myers: self,
            text: text_iter.enumerate(),
            max_dist: max_dist,
            pos: 0,
            _remember: PhantomData,
        }
    }
}

/// The current algorithm state.
#[derive(Clone, Debug, Default)]
struct State<T = u64> {
    pv: T,
    mv: T,
    dist: usize,
}


impl<T> State<T>
where T: Num
{
    /// Create new state.
    pub fn new(m: usize) -> Self {
        State {
            pv: (T::one() << m) - T::one(),
            mv: T::zero(),
            dist: m,
        }
    }
}

/// Iterator over pairs of end positions and distance of matches.
pub struct Matches<'a, T, I>
where T: 'a + Num,
      I: TextIterator<'a>
{
    myers: &'a Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    max_dist: usize,
}


impl<'a, T, I> Iterator for Matches<'a, T, I>
where T: Num,
      I: TextIterator<'a>
 {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<(usize, usize)> {
        for (i, &a) in self.text.by_ref() {
            self.myers.step(&mut self.state, a);
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }
}

// 'Marker' traits
pub struct Remember;
pub struct NoRemember;


/// Iterator over tuples of starting position, end position and distance of matches. In addition,
/// methods for obtaining the hit alignment are provided.
pub struct FullMatches<'a, T, I, S>
where T: 'a + Num,
      I: TextIterator<'a>
{
    myers: &'a mut Myers<T>,
    state: State<T>,
    text: iter::Enumerate<I>,
    m: usize,
    max_dist: usize,
    pos: usize, // current end position, has to be stored for alignment() method
    _remember: PhantomData<S>,
}

impl<'a, T, I, S> FullMatches<'a, T, I, S>
where T: 'a + Num,
      I: TextIterator<'a>
{
    /// Searches the next match and returns a tuple of end position and distance
    /// if found. This involves *no* searching for a starting position and is thus
    /// faster than just iterating over `FullMatches`
    #[inline]
    pub fn next_end(&mut self) -> Option<(usize, usize)> {
        for (i, &a) in self.text.by_ref() {
            self.pos = i; // used in alignment() TODO: where?
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
    pub fn next_path(&mut self, ops: &mut Vec<AlignmentOperation>) -> Option<(usize, usize, usize)> {
        self.next_end()
            .map(|(end, dist)| (self.path(ops), end + 1, dist))
    }

    /// Searches the next match and updates the given `Alignment` with its position
    /// and alignment path if found. The distance is stored in `Alignment::score`.
    /// If no next hit is found, `false` is returned and `aln` remains unchanged.
    #[inline]
    pub fn next_alignment(&mut self, aln: &mut Alignment) -> bool {
        if let Some(_) = self.next_end() {
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
        update_aln_positions(self.pos, len, self.state.dist, self.m, aln);
    }
}


impl<'a, T, I, S> Iterator for FullMatches<'a, T, I, S>
where T: 'a + Num,
      I: TextIterator<'a>
{
    type Item = (usize, usize, usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.next_end()
            .map(|(end, dist)| (self.start(), end + 1, dist))
    }
}


impl<'a, T, I> FullMatches<'a, T, I, Remember>
where T: 'a + Num,
      I: TextIterator<'a>,
{
    /// Takes the end position of a hit and returns a tuple of the corresponding starting position
    /// and the hit distance. If the end position is greater than the end position of the previously
    /// returned hit, `None` is returned.
    /// *Note:* A 0-based end position is expected (as returned by `next_end`).
    pub fn hit_at(&self, end_pos: usize) -> Option<(usize, usize)> {
        self.myers.tb.traceback_at(end_pos, None)
            .map(|(len, dist)| (end_pos + 1 - len, dist))
    }

    /// Takes the end position of a hit and returns a tuple of the corresponding starting position
    /// and the hit distance. The alignment path is added to `ops`.
    /// As in `hit_at`, the end position has to be searched already, otherwise `None` is returned.
    pub fn path_at(&self, end_pos: usize, ops: &mut Vec<AlignmentOperation>) -> Option<(usize, usize)> {
        self.myers.tb.traceback_at(end_pos, Some(ops))
            .map(|(len, dist)| (end_pos + 1 - len, dist))
    }

    /// Takes the end position of a hit and returns a tuple of the corresponding starting position
    /// and the hit distance. The alignment `aln` is updated with the position, alignment path
    /// and distance (stored in `Alignment::score`).
    /// If the end position has not yet been searched, nothing is done and `false` is returned.
    pub fn alignment_at(&self, end_pos: usize, aln: &mut Alignment) -> bool {
        if let Some((aln_len, dist)) = self.myers.tb.traceback_at(end_pos, Some(&mut aln.operations)) {
            update_aln_positions(end_pos, aln_len, dist, self.m, aln);
            return true;
        }
        false
    }
}

// Assumes *0-based* end positions, the coordinates will be converted to 1-based
#[inline(always)]
fn update_aln_positions(end_pos: usize, aln_len: usize, dist: usize, m: usize, aln: &mut Alignment) {
    aln.xstart = 0;
    aln.xend = m;
    aln.xlen = m;
    aln.ylen = aln_len;
    aln.yend = end_pos + 1;
    aln.ystart = aln.yend - aln.ylen;
    aln.mode = AlignmentMode::Semiglobal;
    aln.score = dist as i32;
}


struct Traceback<T: Num> {
    states: Vec<State<T>>,
    positions: iter::Cycle<Range<usize>>,
    pos: usize,
    m: usize,
}

impl<T> Traceback<T>
where T: Num
{
    fn new() -> Traceback<T> {
        Traceback {
            states: vec![],
            positions: (0..0).cycle(),
            m: 0,
            pos: 0,
        }
    }

    fn init(&mut self, num_cols: usize, m: usize) {

        self.positions = (0..num_cols).cycle();

        // extend or truncate states vector
        let curr_len = self.states.len();
        if num_cols > curr_len {
            self.states.reserve(num_cols);
            self.states.extend((0..num_cols - curr_len).map(|_| State::new(0)));
        } else {
            self.states.truncate(num_cols);
        }
        // important if using unsafe in add_state()
        debug_assert!(self.states.len() == num_cols);

        // first column is used to ensure a correct path if the text
        // is shorter than the pattern
        // Actual leftmost column starts at second position
        let state0 = &mut self.states[0];
        state0.dist = m;
        state0.pv = T::max_value(); // all 1s
        // the first position is 1, not 0 (after call to add_state)
        self.pos = self.positions.next().unwrap();

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
    fn traceback(&self, ops: Option<&mut Vec<AlignmentOperation>>) -> (usize, usize) {
        self._traceback_at(self.pos, ops)
    }

    // Returns the length of a match with a given end position, optionally adding the
    // alignment path to `ops`
    // only to be called if the `states` vec contains all states of the text
    fn traceback_at(&self, pos: usize, ops: Option<&mut Vec<AlignmentOperation>>) -> Option<(usize, usize)> {
        let pos = pos + 1; // in order to be comparable since self.pos starts with 1, not 0
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
        mut ops: Option<&mut Vec<AlignmentOperation>>
    )
    -> (usize, usize)
    {

        use self::AlignmentOperation::*;

        // Reverse iterator over states. If remembering all positions,
        // the chain() and cycle() are not actually needed, but there seems
        // to be almost no performance loss.
        let mut states = self.states[..pos + 1].iter().rev().chain(
            self.states.iter().rev().cycle()
        );
        // // Simpler alternative using skip() is slower in some cases:
        // let mut states = self.states.iter().rev().cycle().skip(self.states.len() - pos - 1);

        let ops = &mut ops;
        if let Some(o) = ops.as_mut() {
            o.clear();
        }

        // Mask with one bit switched on at the current position
        let max_mask = T::one() << (self.m - 1);

        macro_rules! move_up {
            ($state:expr) => {
                if $state.pv & max_mask != T::zero() {
                    $state.dist -= 1
                } else if $state.mv & max_mask != T::zero() {
                    $state.dist += 1
                }
                // Sometimes slower, sometimes faster:
                //$state.dist += ($state.mv & max_mask != T::zero()) as usize;
                //$state.dist -= ($state.pv & max_mask != T::zero()) as usize;
                $state.pv <<= 1;
                $state.mv <<= 1;
            };
        }

        macro_rules! move_up_many {
            ($state:expr, $n:expr) => {
                let mask = ((T::one() << $n) - T::one()) << (self.m - $n);
                $state.dist += (($state.mv & mask)).count_ones() as usize;
                $state.dist -= (($state.pv & mask)).count_ones() as usize;
                $state.mv <<= $n;
                $state.pv <<= $n;

                // A loop seems always slower (not sure about systems without popcnt):
                // for _ in 0..$n {
                //     move_up!($state);
                // }
            };
        }

        // horizontal distance from right end
        let mut h_offset = 0;
        // vertical offset from bottom of table
        let mut v_offset = 0;

        // current state
        let mut state = states.next().unwrap().clone();
        // distance of the match
        let dist = state.dist;
        // state left to the current state
        let mut lstate = states.next().unwrap().clone();

        while v_offset < self.m {
            let op =
                if state.pv & max_mask != T::zero() {
                    // up
                    v_offset += 1;
                    move_up!(state);
                    move_up!(lstate);
                    Ins
                } else {
                    let op =
                        if lstate.dist + 1 == state.dist {
                            // left
                            Del
                        } else {
                            // diagonal
                            v_offset += 1;
                            move_up!(lstate);
                            if lstate.dist == state.dist {
                                Match
                            } else {
                                Subst
                            }
                        };
                    // move left
                    state = lstate;
                    lstate = states.next().unwrap().clone();
                    move_up_many!(lstate, v_offset);
                    h_offset += 1;
                    op
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
}


#[cfg(test)]
mod tests {
    use super::*;
    use alignment::AlignmentOperation::*;
    use alignment::{Alignment, AlignmentMode};
    use itertools::Itertools;

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
        let text =    b"TGAGCNT";
        let pattern = b"TGAGCGT";

        let myers = Myers64::new(pattern);
        assert_eq!(myers.distance(text), 1);

        let myers_wildcard = Myers64::with_wildcard(pattern, b'N');
        assert_eq!(myers_wildcard.distance(text), 0);
    }

    #[test]
    fn test_long() {
        let text = b"ACCGTGGATGAGCGCCATAGACCGTGGATGAGCGCCATAGACCGTGGATGAGCGCCATAGACCGTGGATGAGCG\
CCATAGACCGTGGATGAGCGCCATAG";
        let pattern = b"TGAGCGTTGAGCGTTGAGCGTTGAGCGTTGAGCGTTGAGCGT";
        let myers = Myers64::new(&pattern[..]);
        let occ = myers.find_all_end(&text[..], 1).collect_vec();
        // TODO: what is the correct result?
        println!("{:?}", occ);
    }

    #[test]
    fn test_full_position() {
        let text = b"CAGACATCTT";
        let pattern = b"AGA";

        let mut myers = Myers64::new(pattern);
        let matches: Vec<_> = myers.find_all_pos(text, 1).collect();
        assert_eq!(&matches, &[(1, 3, 1), (1, 4, 0), (1, 5, 1), (3, 6, 1)]);
    }

    #[test]
    fn test_traceback_path() {
        let text =    "TCAGACAT-CTT".replace('-', "").into_bytes();
        let pattern = "TC-GACGTGCT".replace('-', "").into_bytes();

        let mut myers = Myers64::new(&pattern);
        let mut matches = myers.find_all_pos(&text, 3);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 10, 3));
        assert_eq!(aln, &[Match, Match, Del, Match, Match, Match, Subst, Match, Ins, Match, Match]);
    }

    #[test]
    fn test_traceback_path2() {
        let text =    "TCAG--CAGATGGAGCTC".replace('-', "").into_bytes();
        let pattern = "TCAGAGCAG".replace('-', "").into_bytes();

        let mut myers = Myers64::new(&pattern);
        let mut matches = myers.find_all_pos(&text, 2);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 7, 2));
        assert_eq!(aln, &[Match, Match, Match, Match, Ins, Ins, Match, Match, Match]);
    }

    #[test]
    fn test_alignment() {
        let text =    "GGTCCTGAGGGA".replace('-', "");
        let pattern =    "CCT-AGGG".replace('-', "");

        let mut myers = Myers64::new(pattern.as_bytes());
        let mut matches = myers.find_all_pos(text.as_bytes(), 1);
        let expected = Alignment {
            score: 1,
            xstart: 0,
            xend: 7,
            xlen: 7,
            ystart: 3,
            yend: 11,
            ylen: 8,
            operations: vec![Match, Match, Match, Del, Match, Match, Match, Match],
            mode: AlignmentMode::Semiglobal,
        };

        let mut aln = Alignment::new();
        assert!(matches.next_alignment(&mut aln));
        assert_eq!(&aln, &expected);

        aln.score = -1;
        matches.alignment(&mut aln);
        assert_eq!(&aln, &expected);
    }

    #[test]
    fn test_position0_at() {
        // same as position_at, but 0-based positions from
        let text = b"CAGACATCTT";
        let pattern = b"AGA";
        let expected_hits = &[(1, 3, 1), (1, 4, 0), (1, 5, 1), (3, 6, 1)];

        let mut myers = Myers64::new(pattern);

        // first, use the standard iterator with 1-based ends
        let matches: Vec<_> = myers.find_all_pos_remember(text, 1).collect();
        assert_eq!(&matches, expected_hits);

        // then, iterate over 0-based ends
        let mut myers_m = myers.find_all_pos_remember(text, 1);
        let mut ends = vec![];
        while let Some(item) = myers_m.next_end() {
            ends.push(item);
        }

        // retrieve start and distance and compare
        for (&(end0, dist), &(exp_start, exp_end1, exp_dist)) in ends.iter().zip(expected_hits) {
            // note: range end is not 0-based position -> supply end - 1
            assert_eq!(end0 + 1, exp_end1);
            assert_eq!(dist, exp_dist);
            assert_eq!(myers_m.hit_at(end0), Some((exp_start, dist)));
        }
    }

    #[test]
    fn test_position_at() {
        let text =   b"CAGACATCTT";
        let pattern = b"AGA";

        let mut myers = Myers64::new(pattern);
        let mut myers_m = myers.find_all_pos_remember(text, 1);
        let matches: Vec<_> = (&mut myers_m).collect();
        assert_eq!(&matches, &[(1, 3, 1), (1, 4, 0), (1, 5, 1), (3, 6, 1)]);

        for (start, end, dist) in matches {
            // note: range end is not 0-based position -> supply end - 1
            assert_eq!(myers_m.hit_at(end - 1), Some((start, dist)));
        }
    }

    #[test]
    fn test_path_at() {
        let text =   b"CAGACATCTT";
        let pattern = b"AGA";

        let mut myers = Myers64::new(pattern);
        let mut matches = myers.find_all_pos_remember(text, 1);

        let expected = &[Match, Match, Ins];
        // search first hit
        let mut path = vec![];
        assert_eq!(matches.next_path(&mut path).unwrap(), (1, 3, 1));
        assert_eq!(&path, expected);
        path.clear();

        // retrieve first hit at 0-based end position (2)
        assert_eq!(matches.hit_at(2), Some((1, 1)));
        assert_eq!(matches.path_at(2, &mut path), Some((1, 1)));
        assert_eq!(&path, expected);

        // hit out of range
        path.clear();
        assert!(matches.path_at(3, &mut path).is_none());
        assert!(path.is_empty());

        // now search the next hit
        assert_eq!(matches.next_path(&mut path).unwrap(), (1, 4, 0));
        // position 3 is now searched -> path can be retrieved
        assert_eq!(matches.path_at(3, &mut path), Some((1, 0)));
        assert_eq!(&path, &[Match, Match, Match])
    }

    #[test]
    fn test_shorter() {
        let text =     "ATG";
        let pattern = "CATGC";

        let mut myers = Myers64::new(pattern.as_bytes());
        let mut matches = myers.find_all_pos(text.as_bytes(), 2);
        let mut aln = vec![];
        assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 3, 2));
        assert_eq!(aln, &[Ins, Match, Match, Match, Ins]);
    }

    #[test]
    fn test_long_shorter() {
        let text =           "CCACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT";
        let pattern ="CGCGGTGTCCACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT";

        let mut myers = Myers64::new(pattern.as_bytes());
        let mut matches = myers.find_all_pos(text.as_bytes(), 8);
        assert_eq!(matches.next().unwrap(), (0, 52, 8));
    }

    #[test]
    fn test_ambig() {
        let text =    b"TGABCNT";
        let pattern = b"TGRRCGT";
        //                x  x
        // Matching is asymmetric here (A matches R and G matches N, but the reverse is not true)

        let myers = Myers64::with_ambig(pattern.into_iter().map(|&b| {
            if b == b'R' { b"AG".to_vec() } else { vec![b] }
        }));
        assert_eq!(myers.distance(text), 2);
    }
}
