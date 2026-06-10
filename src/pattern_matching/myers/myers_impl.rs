use std::fmt;

use crate::pattern_matching::myers::{word_size, BitVec, DistType};
use num_traits::ToPrimitive;

/// The current algorithm state.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub(crate) struct State<T, D>
where
    T: BitVec,
    D: std::fmt::Debug,
{
    pub pv: T,
    pub mv: T,
    pub dist: D,
}

impl<T, D> State<T, D>
where
    T: BitVec,
    D: DistType,
{
    /// Create and initiate leftmost block given the maximum pattern length
    // it should represent (m)
    #[inline]
    pub fn init(m: D) -> Self {
        State {
            pv: T::max_value(),
            mv: T::zero(),
            dist: m,
        }
    }

    /// Initialize with the maximum distance
    #[inline]
    pub fn init_max_dist() -> Self {
        Self::init(D::max_value())
    }

    /// Returns the known edit distance: this is meant to pair with
    /// `long::States::known_dist()` and should not be `None`, as there
    /// is only one "block" with the simple algorithm, and the edit distance
    /// is thus always known after `Myers::step()`.
    #[inline]
    pub fn known_dist(&self) -> Option<D> {
        Some(self.dist)
    }

    // unmodified State::default() (only used for debugging)
    pub fn is_new(&self) -> bool {
        self.dist == D::zero() && self.pv == T::zero() && self.mv == T::zero()
    }

    // unmodified State::init_max_dist() (only used for debugging)
    pub fn is_max(&self) -> bool {
        self.pv >= (T::max_value() >> 1) && self.mv == T::zero()
    }

    /// Adjust the distance score of the current cell C(j) by n positions
    /// up to cell C(j - n) ('move up' in the DP matrix column).
    /// The provided `range_mask` must have bits activated at the positions j-1 to j-n,
    /// (without the current position j).
    /// See also `adjust_one_up` for 'moving up' only by one position.
    ///
    /// This function is called in the traceback phase when moving to the left
    /// in the DP matrix.
    #[inline]
    pub fn adjust_up_by(&mut self, range_mask: T) {
        let p = (self.pv & range_mask).count_ones();
        let m = (self.mv & range_mask).count_ones();
        let mut dist = self.dist.to_u64().unwrap();
        dist = dist.wrapping_add(m.to_u64().unwrap());
        dist = dist.wrapping_sub(p.to_u64().unwrap());
        self.dist = D::from_u64(dist).unwrap();
    }

    /// This method may be used for performance comparison instead of adjust_up_by(),
    /// but requires a single bit to be set in `pos_mask` (as in `adjust_one_up`,
    /// not a whole range).
    #[inline]
    #[allow(dead_code)]
    pub fn _adjust_up_by(&mut self, mut pos_mask: T, n: usize) {
        for _ in 0..n {
            self.adjust_one_up(pos_mask);
            pos_mask <<= 1;
        }
    }

    /// Modify the distance score of the current cell C(j) to C(j-1)
    /// by applying the vertical delta (= move one position up in the DP matrix).
    /// Expects a mask with a single activated bit at position j - 1.
    /// It is the caller's responsibility to track the current position in the
    /// DP matrix.
    #[inline]
    pub fn adjust_one_up(&mut self, pos_mask: T) {
        //debug_assert!(!self.is_max());
        if self.pv & pos_mask != T::zero() {
            self.dist -= D::one();
        } else if self.mv & pos_mask != T::zero() {
            self.dist += D::one();
        }

        // not faster:
        // let diff = ((self.mv & pos_mask) != T::zero()) as isize - ((self.pv & pos_mask) != T::zero()) as isize;
        // self.dist = D::from_usize(self.dist.to_usize().unwrap().wrapping_add(diff as usize)).unwrap();
    }
}

impl<T, D> fmt::Display for State<T, D>
where
    T: BitVec,
    D: std::fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "State {{\n  + {:0width$b}  (dist = {:?})\n  - {:0width$b}\n}}",
            self.pv,
            self.dist,
            self.mv,
            width = word_size::<T>() // TODO: T::BITS
        )
    }
}

#[rustfmt::skip]
// rustfmt::skip prevents automatic indentation.
// This is not optimal, as no checks are done at all..
macro_rules! impl_myers {
    ($DistType:ty, $max_dist:expr, $Myers:ty, $State:ty, $TbHandler:ty) => {
        mod myers_impl {
// Macro implementing common methods in Myers object. Wrapped in a module
// and then re-exported from there to avoid mixing of namespaces.
// Indented at top level for readability.

use super::Myers;
use crate::alignment::{Alignment, AlignmentOperation};
use crate::pattern_matching::myers::traceback::Traceback;
use crate::pattern_matching::myers::{update_aln, BitVec};
#[allow(unused_imports)] // Bounded is required for <$DistType>::max_value()
use num_traits::{Bounded, ToPrimitive};
use std::borrow::Borrow;
use std::cmp::min;
use std::iter;

impl<T: BitVec> $Myers {
    // Combining these two steps into one function seems beneficial for performance
    #[inline]
    fn step_trace<'a>(
        &mut self,
        mut state: &mut $State,
        a: u8,
        max_dist: $DistType,
        traceback: &mut Traceback<'a, T, $DistType, $TbHandler>,
    ) {
        self.step(&mut state, a, max_dist);
        traceback.add_state(&state, &mut self.states_store);
    }

    /// Calculate the distance of the pattern to the given text.
    pub fn distance<C, I>(&self, text: I) -> $DistType
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
    {
        let max_dist = $max_dist;
        let mut dist = max_dist;
        let m = self.m;
        let mut state = self.initial_state(m, max_dist);
        for a in text {
            self.step(&mut state, *a.borrow(), max_dist);
            if let Some(d) = state.known_dist() {
                if d < dist {
                    dist = d;
                }
            }
        }
        dist
    }

    /// Finds all matches of pattern in the given text up to a given maximum distance.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn find_all_end<C, I>(
        &self,
        text: I,
        max_dist: $DistType,
    ) -> Matches<'_, T, C, I::IntoIter>
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
    {
        Matches::new(self, text.into_iter(), max_dist.min($max_dist))
    }

    /// Find the best match of the pattern in the given text.
    /// if multiple end positions have the same distance, the first is returned.
    pub fn find_best_end<C, I>(&self, text: I) -> (usize, $DistType)
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
    {
        self.find_all_end(text, $max_dist)
            .min_by_key(|&(_, dist)| dist)
            .unwrap()
    }

    /// Finds all matches of pattern in the given text up to a given maximum distance.
    /// In contrast to `find_all_end`, matches are returned as an iterator over ranges
    /// of `(start, end, distance)`. Note that the end coordinate is not included in the
    /// range and thus and thus greater by one compared to the end index returned by
    /// `find_all_end()`.
    pub fn find_all<C, I>(
        &mut self,
        text: I,
        max_dist: $DistType,
    ) -> FullMatches<'_, T, C, I::IntoIter>
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
        I::IntoIter: ExactSizeIterator,
    {
        FullMatches::new(self, text.into_iter(), max_dist.min($max_dist))
    }

    /// As `find_all_end`, this function returns an iterator over tuples of `(end, distance)`.
    /// Additionally, it keeps the data necessary for later obtaining the starting positions and/or
    /// the alignment path at *any* position that was already searched.
    pub fn find_all_lazy<C, I>(
        &mut self,
        text: I,
        max_dist: $DistType,
    ) -> LazyMatches<'_, T, C, I::IntoIter>
    where
        C: Borrow<u8>,
        I: IntoIterator<Item = C>,
        I::IntoIter: ExactSizeIterator,
    {
        LazyMatches::new(self, text.into_iter(), max_dist.min($max_dist))
    }
}

/// Iterator over pairs of end positions and distance of matches.
#[derive(Clone, Debug)]
pub struct Matches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    myers: &'a $Myers,
    state: $State,
    text: iter::Enumerate<I>,
    max_dist: $DistType,
}

impl<'a, T, C, I> Matches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    fn new(myers: &'a Myers<T>, text: I, max_dist: $DistType) -> Self {
        let m = myers.m;
        let state = myers.initial_state(m, max_dist);
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
    type Item = (usize, $DistType);

    fn next(&mut self) -> Option<(usize, $DistType)> {
        for (i, a) in self.text.by_ref() {
            self.myers.step(&mut self.state, *a.borrow(), self.max_dist);
            if let Some(dist) = self.state.known_dist() {
                if dist <= self.max_dist {
                    return Some((i, dist));
                }
            }
        }
        None
    }
}

/// Iterator over tuples of (start, end, distance).
/// Further methods for obtaining the hit alignment path are provided.
#[derive(Debug)]
pub struct FullMatches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    myers: &'a mut $Myers,
    traceback: Traceback<'a, T, $DistType, $TbHandler>,
    state: $State,
    text: iter::Enumerate<I>,
    text_len: usize,
    m: $DistType,
    max_dist: $DistType,
    pos: usize, // current end position, has to be stored for alignment() method
    unsuccessfully_finished: bool,
}

impl<'a, T, C, I> FullMatches<'a, T, C, I>
where
    T: 'a + BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C> + ExactSizeIterator,
{
    fn new(myers: &'a mut $Myers, text_iter: I, max_dist: $DistType) -> Self {
        let state = myers.initial_state(myers.m, max_dist);
        // Calculate maximum number of traceback columns possibly used by a match
        let m = myers.m.to_usize().unwrap();
        let num_cols = m + min(max_dist.to_usize().unwrap(), m);
        let tb = Traceback::new(
            &mut myers.states_store,
            &state,
            num_cols.to_usize().unwrap(),
            myers.m,
            <$TbHandler>::new(),
        );
        FullMatches {
            m: myers.m,
            myers,
            traceback: tb,
            state,
            text_len: text_iter.len(),
            text: text_iter.enumerate(),
            max_dist,
            pos: 0,
            unsuccessfully_finished: false,
        }
    }

    /// Searches the next match and returns a tuple of end position and distance
    /// if found. This involves *no* searching for a starting position and is thus
    /// faster than just iterating over `FullMatches`
    #[inline]
    pub fn next_end(&mut self) -> Option<(usize, $DistType)> {
        for (i, a) in self.text.by_ref() {
            self.pos = i; // used in alignment()
            self.myers.step_trace(
                &mut self.state,
                *a.borrow(),
                self.max_dist,
                &mut self.traceback,
            );
            if let Some(dist) = self.state.known_dist() {
                if dist <= self.max_dist {
                    return Some((i, dist));
                }
            }
        }
        self.unsuccessfully_finished = true;
        None
    }

    /// Searches the next match and returns a tuple of starting position, end position and
    /// distance, or `None` if no match was found. In addition, the alignment path is added to
    /// `ops`. Existing data in the vector will be cleared beforehand.
    #[inline]
    pub fn next_path(
        &mut self,
        ops: &mut Vec<AlignmentOperation>,
    ) -> Option<(usize, usize, $DistType)> {
        self.next_end()
            .map(|(end, dist)| (self.path(ops).unwrap(), end + 1, dist))
    }

    /// Like `FullMatches::path_reverse()`, but the operations will be in reverse order. This
    /// is slightly faster, as the traceback algorithm adds them in reverse order,
    /// and `path()` needs to reverse them. Existing data in the vector will be cleared
    /// beforehand.
    #[inline]
    pub fn next_path_reverse(
        &mut self,
        ops: &mut Vec<AlignmentOperation>,
    ) -> Option<(usize, usize, $DistType)> {
        self.next_end()
            .map(|(end, dist)| (self.path_reverse(ops).unwrap(), end + 1, dist))
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
        if self.unsuccessfully_finished {
            return None;
        }
        let (len, _) = self
            .traceback
            .traceback(None, &self.myers.states_store)
            .unwrap();
        Some(self.pos + 1 - len.to_usize().unwrap())
    }

    /// Adds the path of the current hit alignment to `ops` and returns the starting position of
    /// the current hit. If the search is finished and no hit was found, `None` is returned.
    /// Adds the path of the current hit alignment to `ops` and returns the starting position of
    /// the current hit. If the search is finished and no hit was found, `None` is returned.
    /// Existing data in the vector will be cleared beforehand.
    #[inline]
    pub fn path(&self, ops: &mut Vec<AlignmentOperation>) -> Option<usize> {
        self.path_reverse(ops).map(|pos| {
            ops.reverse();
            pos
        })
    }

    /// Like `FullMatches::path()`, but the operations will be in reverse order. This
    /// is slightly faster, as the traceback algorithm adds them in reverse order,
    /// and `path()` needs to reverse them.
    /// Existing data in the vector will be cleared beforehand.
    #[inline]
    pub fn path_reverse(&self, ops: &mut Vec<AlignmentOperation>) -> Option<usize> {
        if self.unsuccessfully_finished {
            return None;
        }
        ops.clear();
        let (len, _) = self
            .traceback
            .traceback(Some(ops), &self.myers.states_store)
            .unwrap();
        Some(self.pos + 1 - len.to_usize().unwrap())
    }

    /// Updates the given `Alignment` with its position and alignment path. The edit distance is
    /// stored in `Alignment::score`. If no hit has been found yet, then `false` will be returned
    /// and nothing is done.
    #[inline]
    pub fn alignment(&mut self, aln: &mut Alignment) -> bool {
        if self.unsuccessfully_finished {
            return false;
        }
        if let Some(dist) = self.state.known_dist() {
            aln.operations.clear();
            let (len, _) = self
                .traceback
                .traceback(Some(&mut aln.operations), &self.myers.states_store)
                .unwrap();
            aln.operations.reverse();
            update_aln(
                self.pos,
                len.to_usize().unwrap(),
                self.text_len,
                dist.to_usize().unwrap(), //self.state.dist().to_usize().unwrap(),
                self.m.to_usize().unwrap(),
                aln,
            );
            true
        } else {
            false
        }
    }
}

impl<'a, T, C, I> Iterator for FullMatches<'a, T, C, I>
where
    T: 'a + BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C> + ExactSizeIterator,
{
    type Item = (usize, usize, $DistType);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.next_end()
            .map(|(end, dist)| (self.start().unwrap(), end + 1, dist))
    }
}

/// Iterator over tuples of (end, distance).
/// Further methods for obtaining the hit alignment path are provided.
#[derive(Debug)]
pub struct LazyMatches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C>,
{
    myers: &'a mut $Myers,
    traceback: Traceback<'a, T, $DistType, $TbHandler>,
    state: $State,
    text: iter::Enumerate<I>,
    text_len: usize,
    m: $DistType,
    max_dist: $DistType,
}

impl<'a, T, C, I> Iterator for LazyMatches<'a, T, C, I>
where
    T: BitVec,
    C: Borrow<u8>,
    I: Iterator<Item = C> + ExactSizeIterator,
{
    type Item = (usize, $DistType);

    #[inline]
    fn next(&mut self) -> Option<(usize, $DistType)> {
        for (i, a) in self.text.by_ref() {
            self.myers.step_trace(
                &mut self.state,
                *a.borrow(),
                self.max_dist,
                &mut self.traceback,
            );
            // self.traceback
            //     .add_state(&self.state, &mut self.myers.states_store);
            if let Some(dist) = self.state.known_dist() {
                if dist <= self.max_dist {
                    return Some((i, dist));
                }
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
    #[inline]
    fn new(myers: &'a mut Myers<T>, text_iter: I, max_dist: $DistType) -> Self {
        let state = myers.initial_state(myers.m, max_dist);
        let tb = Traceback::new(
            &mut myers.states_store,
            &state,
            text_iter.len(),
            myers.m,
            <$TbHandler>::new(),
        );
        LazyMatches {
            m: myers.m,
            myers,
            traceback: tb,
            state,
            text_len: text_iter.len(),
            text: text_iter.enumerate(),
            max_dist,
        }
    }

    /// Returns the edit distance of a match given its end position, or `None`
    /// if the alignment is not known (see [`alignment_at`](Self::alignment_at)).
    #[inline]
    pub fn dist_at(&self, end_pos: usize) -> Option<$DistType> {
        self.traceback.dist_at(end_pos, &self.myers.states_store)
    }

    /// Returns a tuple of (start position, edit distance) for a match given its
    /// end position. This involves caluculating the alignment path.
    /// The function returns `None` if the alignment is not known
    /// (see [`alignment_at`](Self::alignment_at)).
    #[inline]
    pub fn hit_at(&self, end_pos: usize) -> Option<(usize, $DistType)> {
        self.traceback
            .traceback_at(end_pos, None, &self.myers.states_store)
            .map(|(len, dist)| (end_pos + 1 - len.to_usize().unwrap(), dist))
    }

    /// Calculates the alignment path for a match given its end position.
    /// Returns a tuple of (start position, edit distance) or `None` if the
    /// alignment is not known (see [`alignment_at`](Self::alignment_at)).
    #[inline]
    pub fn path_at(
        &self,
        end_pos: usize,
        ops: &mut Vec<AlignmentOperation>,
    ) -> Option<(usize, $DistType)> {
        self.path_at_reverse(end_pos, ops).map(|rv| {
            ops.reverse();
            rv
        })
    }

    /// Like [`path_at`](Self::path_at), but the operations will be in reverse order.
    /// This is slightly faster, as the traceback algorithm adds them in reverse order,
    /// and `path_at()` needs to reverse them.
    #[inline]
    pub fn path_at_reverse(
        &self,
        end_pos: usize,
        ops: &mut Vec<AlignmentOperation>,
    ) -> Option<(usize, $DistType)> {
        self.traceback
            .traceback_at(end_pos, Some(ops), &self.myers.states_store)
            .map(|(len, dist)| (end_pos + 1 - len.to_usize().unwrap(), dist))
    }

    /// Calculates the alignment for a match given its end position, which can
    /// be obtained with `myers.find_all_lazy(text, max_dist)`.
    /// The supplied [`Alignment`](crate::alignment::Alignment) is updated with the
    /// position, alignment path and edit distance (stored in the `score` field).
    ///
    /// Returns `true` if successful or `false` if the alignment is not known
    /// because the given position has not yet been searched ([`LazyMatches`]
    /// iterator not exhausted).
    ///
    /// # Note on more distant matches
    ///
    /// In the simple case with short patterns, this function succeeds for *all*
    /// end positions, even if hits are not within `max_dist`. In contrast, the
    /// the block-based search ([`long`](crate::pattern_matching::myers::long)
    /// submodule) is only guaranteed to know the alignment for matches <= `max_dist`
    /// and may return `false` otherwise.
    #[inline]
    pub fn alignment_at(&self, end_pos: usize, aln: &mut Alignment) -> bool {
        aln.operations.clear();
        if let Some((aln_len, dist)) = self.traceback.traceback_at(
            end_pos,
            Some(&mut aln.operations),
            &self.myers.states_store,
        ) {
            aln.operations.reverse();
            update_aln(
                end_pos,
                aln_len.to_usize().unwrap(),
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

// module end
}

pub use myers_impl::*;

// macro end
};
}
