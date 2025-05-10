//! Block-based version of the algorithm, which does not restrict pattern length.
//!
//! This module implements the block-based version of the Myers pattern matching algorithm.
//! It can be used for searching patterns of any length and obtaining semiglobal alignments
//! of the hits. Apart from that, the `Myers` object in this module provides exactly the same
//! API as the 'simple' version `bio::pattern_matching::myers::Myers`.
//! For short patterns, the 'simple' version is still to be preferred, as the block-based
//! algorithm is slower.

use std::borrow::Borrow;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::iter;
use std::marker::PhantomData;
use std::mem::replace;
use std::slice;

use itertools::Itertools;
use num_traits::ToPrimitive;

use crate::pattern_matching::myers::traceback::{StatesHandler, TracebackHandler};
use crate::pattern_matching::myers::{ceil_div, State};

use super::word_size;
use super::BitVec;

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
struct Peq<T: BitVec> {
    peq: [T; 256],
    // a bit mask with bits set for all positions covered by the pattern
    bit_mask: T,
}

/// Myers algorithm.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Myers<T = u64>
where
    T: BitVec,
{
    peq: Vec<Peq<T>>,
    pub(crate) m: usize,
    pub(crate) states_store: Vec<State<T, usize>>,
}

impl<T: BitVec> Myers<T> {
    /// Create a new instance of Myers algorithm for a given pattern.
    #[inline]
    pub fn new<P, C>(pattern: P) -> Self
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        Self::new_ambig(pattern, None, None)
    }

    #[inline]
    pub(crate) fn new_ambig<P, C>(
        pattern: P,
        opt_ambigs: Option<&HashMap<u8, Vec<u8>>>,
        opt_wildcards: Option<&[u8]>,
    ) -> Self
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        let w = word_size::<T>();
        let pattern = pattern.into_iter();
        let m = pattern.len();
        assert!(m > 0, "Pattern is empty");
        assert!(m <= usize::MAX / 2, "Pattern too long");

        // build peq
        let mut peq = Vec::with_capacity(ceil_div(m, w));
        for chunk in pattern.chunks(w).into_iter() {
            let mut peq_block = [T::zero(); 256];
            let mut chunk_len = 0;
            for symbol in chunk {
                let symbol = *symbol.borrow();
                let mask = T::one() << chunk_len;
                // equivalent
                peq_block[symbol as usize] |= mask;
                // ambiguities
                if let Some(equivalents) = opt_ambigs.and_then(|ambigs| ambigs.get(&symbol)) {
                    for &eq in equivalents {
                        peq_block[eq as usize] |= mask;
                    }
                }
                chunk_len += 1;
            }
            // wildcards
            if let Some(wildcards) = opt_wildcards {
                for &w in wildcards {
                    peq_block[w as usize] = T::max_value();
                }
            }

            peq.push(Peq {
                peq: peq_block,
                bit_mask: T::one() << (chunk_len - 1),
            });
        }

        Myers {
            peq,
            m,
            states_store: Vec::new(),
        }
    }

    #[inline]
    fn step(&self, state: &mut States<T>, a: u8, max_dist: usize) {
        state.step(a, &self.peq, max_dist)
    }

    #[inline]
    fn initial_state(&self, m: usize, max_dist: usize) -> States<T> {
        States::new(m, max_dist)
    }
}

#[inline]
fn advance_block<T: BitVec>(state: &mut State<T, usize>, p: &Peq<T>, a: u8, hin: i8) -> i8 {
    let mut eq = p.peq[a as usize];
    let xv = eq | state.mv;
    if hin < 0 {
        eq |= T::one();
    }
    let xh = ((eq & state.pv).wrapping_add(&state.pv) ^ state.pv) | eq;

    let mut ph = state.mv | !(xh | state.pv);
    let mut mh = state.pv & xh;

    // let hout = if ph & p.bound > T::zero() {
    //     state.dist += 1;
    //     1
    // } else if mh & p.bound > T::zero() {
    //     state.dist -= 1;
    //     -1
    // } else {
    //     0
    // };

    // apparently faster than commented code above
    let mut hout = ((ph & p.bit_mask) != T::zero()) as i8;
    hout -= ((mh & p.bit_mask) != T::zero()) as i8;
    state.dist = state.dist.wrapping_add(hout as usize);

    ph <<= 1;
    mh <<= 1;

    if hin < 0 {
        mh |= T::one();
    }
    if hin > 0 {
        ph |= T::one();
    }
    // not faster:
    // mh |= T::from_u8((hin < 0) as u8).unwrap();
    // ph |= T::from_u8((hin > 0) as u8).unwrap();

    state.pv = mh | !(xv | ph);
    state.mv = ph & xv;

    hout
}

/// Multi-block state
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub(super) struct States<T: BitVec> {
    states: Vec<State<T, usize>>,
    // index of the last block (max_block_i + 1 blocks are needed to cover the full pattern)
    max_block_i: usize,
    // length of the remaining pattern chunk in the last block
    // (see `add_block()`)
    last_m: usize,
}

impl<T> States<T>
where
    T: BitVec,
{
    fn new(m: usize, max_dist: usize) -> Self {
        let w = word_size::<T>();
        let nblock = ceil_div(m, w);
        let mut s = States {
            states: Vec::with_capacity(nblock),
            max_block_i: nblock - 1,
            last_m: m % w,
        };
        // y = ceil(k/w) in Myers' paper, Fig. 9 where k = max_dist
        let min_blocks = max(1, ceil_div(min(max_dist, m), w));
        for _ in 0..min_blocks {
            s.add_block(0);
        }
        s
    }

    /// Adds a new block to the Ukkonen band
    /// `carry` (`vin` in Myers' paper) can be -1, 0, or 1
    #[inline]
    fn add_block(&mut self, carry: i8) {
        let prev_dist = self.states.last().map(|s| s.dist).unwrap_or(0);
        // delta: number of bits in pv/mv covered by the pattern
        let delta = if self.states.len() == self.max_block_i && self.last_m > 0 {
            // For the last ('bottom') block, we add the length of the remaining pattern
            // chunk = m % w (stored in self.last_m) and ignore the remaining bits.
            // This strategy differs from the solution by Myers (p. 407, note 4 / p. 410, Fig. 9),
            // where the computation starts with distance score = w * b [where b = N blocks],
            // and both the pattern and sequence need to be padded with wild-card symbols.
            self.last_m
        } else {
            word_size::<T>()
        };
        self.states.push(State::init(
            prev_dist
                .wrapping_add(delta)
                .wrapping_add(carry as usize)
                .to_usize()
                .unwrap(),
        ));
    }

    #[inline]
    fn step(&mut self, a: u8, peq: &[Peq<T>], max_dist: usize) {
        let mut carry = 0;
        let mut y = self.states.len() - 1;

        // compute blocks
        for (state, block_peq) in self.states.iter_mut().zip(peq) {
            carry = advance_block(state, block_peq, a, carry);
        }

        // adjust Ukkonen band if necessary
        let w = word_size::<T>();
        let last_dist = self.states[y].dist;
        // note: this will fail if edit distances are > isize::MAX,
        // but in practice such long patterns will never occur
        if (last_dist as isize - carry as isize) as usize <= max_dist
            && y < self.max_block_i
            && (peq[y + 1].peq[a as usize] & T::one() == T::one() || carry < 0)
        {
            y += 1;
            self.add_block(-carry);
            // compute the new block as well
            advance_block(&mut self.states[y], &peq[y], a, carry);
        } else {
            // note: max_dist should be <= usize::MAX - w; this is ensured beforehand
            while y > 0 && self.states[y].dist >= max_dist + w {
                y -= 1;
            }
            self.states.truncate(y + 1);
        }
    }

    /// Returns the edit distance if known (i.e., all blocks were computed),
    #[inline]
    fn known_dist(&self) -> Option<usize> {
        self.states.get(self.max_block_i).map(|s| s.dist)
    }
}

#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub(super) struct LongStatesHandler<'a> {
    n_blocks: usize,
    _a: PhantomData<&'a ()>,
}

impl<'a> LongStatesHandler<'a> {
    #[inline]
    pub fn new() -> Self {
        LongStatesHandler {
            n_blocks: 0,
            _a: PhantomData,
        }
    }
}

impl<'a, T> StatesHandler<'a, T, usize> for LongStatesHandler<'a>
where
    T: BitVec + 'a,
{
    type TracebackHandler = LongTracebackHandler<'a, T>;
    type TracebackColumn = States<T>;

    #[inline]
    fn init(&mut self, n: usize, m: usize) -> usize {
        let w = word_size::<T>();
        self.n_blocks = ceil_div(m.to_usize().unwrap(), w);
        n * self.n_blocks
    }

    #[inline]
    fn n_blocks(&self) -> usize {
        self.n_blocks
    }

    #[inline]
    fn set_max_state(&self, pos: usize, states: &mut [State<T, usize>]) {
        let pos = pos * self.n_blocks;
        for s in states.iter_mut().skip(pos).take(self.n_blocks) {
            *s = State::init_max_dist();
        }
    }

    #[inline]
    fn add_state(
        &self,
        source: &Self::TracebackColumn,
        pos: usize,
        states: &mut [State<T, usize>],
    ) {
        let source = &source.states;
        let pos = pos * self.n_blocks;
        states[pos..pos + source.len()].clone_from_slice(source);

        if source.len() < self.n_blocks {
            // While following the traceback path, it can happen that the block to the
            // left was not computed because it is outside of the band.
            // Here we initialize the block below the last computed block with meaningful
            // defaults that act as a "barrier" preventing the algorithm from moving
            // left to this block.
            states[pos + source.len()] = State {
                dist: usize::MAX,
                pv: T::zero(),
                mv: T::zero(),
            };
            // Furthermore, in case arbitrary matches are requested with the
            // `LazyMatches::..._at()` methods, we need to know if all blocks
            // have been computed in the given column. For this, we set the
            // distance of the last block to `usize::MAX`.
            states[pos + self.n_blocks - 1].dist = usize::MAX;
        }
    }

    #[inline]
    fn dist_at(&self, pos: usize, states: &[State<T, usize>]) -> Option<usize> {
        states.get(self.n_blocks * (pos + 1) - 1).and_then(|s| {
            if s.dist == usize::MAX {
                None
            } else {
                Some(s.dist)
            }
        })
    }

    #[inline]
    fn init_traceback(
        &self,
        m: usize,
        pos: usize,
        states: &'a [State<T, usize>],
    ) -> Option<Self::TracebackHandler> {
        LongTracebackHandler::new(self.n_blocks, m, pos, states)
    }
}

type RevColIter<'a, T> = iter::Rev<slice::Chunks<'a, State<T, usize>>>;

pub(super) struct LongTracebackHandler<'a, T: BitVec> {
    // reverse iterator used to set `col` and `left_col`
    states_iter: iter::Chain<RevColIter<'a, T>, iter::Cycle<RevColIter<'a, T>>>,
    // slice of blocks representing the current DP matrix column
    col: &'a [State<T, usize>],
    // slice of blocks representing the DP matrix column to the left
    left_col: &'a [State<T, usize>],
    // current block
    block: State<T, usize>,
    // current left block
    left_block: State<T, usize>,
    // index of `block` in `col`
    block_idx: usize,
    // index of `left_block` in `left_col`
    left_block_idx: usize,
    // mask with one active bit indicating the position of the last letter
    max_mask: T,
    // mask with one active bit indicating the current position (row) in the current block
    pos_mask: T,
    // Bit mask that "covers" the range of rows *below* the left cell
    // in the current block. This is used used to initialize blocks
    // (adjust the distance score) with `State::adjust_up_by(left_adj_mask)`.
    left_adj_mask: T,
    _a: PhantomData<&'a ()>,
}

impl<'a, T: BitVec> LongTracebackHandler<'a, T> {
    #[inline]
    fn new(n_blocks: usize, m: usize, pos: usize, states: &'a [State<T, usize>]) -> Option<Self> {
        // length of the pattern chunk in the lowermost block
        let mut last_m = m.to_usize().unwrap() % word_size::<T>();
        if last_m == 0 {
            last_m = word_size::<T>();
        }
        // mask with single bit indicating index of last letter in the block
        let mask0 = T::one() << (last_m - 1);

        // reverse iterator over DP matrix columns
        // chain + cycle is needed in `FullMatches` (see comment in simple.rs)
        let pos = n_blocks * (pos + 1);
        let mut states_iter = states[..pos]
            .chunks(n_blocks)
            .rev()
            .chain(states.chunks(n_blocks).rev().cycle());

        let col = states_iter.next().unwrap();
        let left_col = states_iter.next().unwrap();

        // If the last block was not computed, we cannot obtain an alignment
        if col.last().unwrap().dist == usize::MAX {
            return None;
        }

        // initial data for left block
        let (left_block_idx, left_adj_mask, max_mask) = if last_m == 1 && n_blocks > 1 {
            (n_blocks - 2, T::zero(), T::one() << (word_size::<T>() - 1))
        } else {
            (n_blocks - 1, mask0, mask0)
        };

        // the left cell has to be  in diagonal position (move one up)
        let mut left_block = left_col[left_block_idx];
        left_block.adjust_up_by(left_adj_mask);

        Some(LongTracebackHandler {
            block: col[n_blocks - 1],
            left_block,
            col,
            left_col,
            block_idx: n_blocks - 1,
            left_block_idx,
            states_iter,
            pos_mask: mask0,
            left_adj_mask,
            max_mask,
            _a: PhantomData,
        })
    }

    /// Adjust 'left_adj_mask' to cover one more cell above.
    /// This may involve switching to the block on top
    /// *before* the mask would cover the whole block.
    /// Actually, only the masks and block index are adjusted,
    /// not the `State` itself.
    #[inline]
    fn adjust_left_up(&mut self) -> bool {
        let at_boundary = self.left_adj_mask & T::from_usize(0b10).unwrap() != T::zero()
            && self.left_block_idx > 0;
        if !at_boundary {
            self.left_adj_mask = (self.left_adj_mask >> 1) | self.max_mask;
        } else {
            self.max_mask = T::one() << (word_size::<T>() - 1);
            self.left_adj_mask = T::zero();
            self.left_block_idx -= 1;
        }
        // println!("left {:064b}\nmax  {:064b}", self.left_adj_mask, self.left_max_mask);
        at_boundary
    }
}

impl<'a, T: BitVec + 'a> TracebackHandler<'a, T, usize> for LongTracebackHandler<'a, T> {
    #[inline]
    fn dist(&self) -> usize {
        self.block.dist
    }

    #[inline]
    fn left_dist(&self) -> usize {
        self.left_block.dist
    }

    #[inline]
    fn try_move_up(&mut self) -> bool {
        if self.block.pv & self.pos_mask != T::zero() {
            self.move_up();
            return true;
        }
        false
    }

    #[inline]
    fn move_up(&mut self) {
        // (1) move up in current column
        // If the block boundary has not been reached yet, we can shift to the
        // upper position. Otherwise, we move to a new block (if there is one!)
        if self.pos_mask != T::one() || self.block_idx == 0 {
            self.block.adjust_one_up(self.pos_mask);
            self.pos_mask >>= 1;
        } else {
            // move to upper block
            self.pos_mask = T::one() << (word_size::<T>() - 1);
            self.block_idx -= 1;
            self.block = self.col[self.block_idx];
        }
        // (2) move up in left column
        let at_boundary = self.adjust_left_up();
        if !at_boundary {
            self.left_block.adjust_one_up(self.pos_mask);
        } else {
            self.left_block = self.left_col[self.left_block_idx];
        }
    }

    #[inline]
    fn prepare_diagonal(&mut self) {
        self.adjust_left_up();
        if self.pos_mask != T::one() || self.block_idx == 0 {
            // not at block boundary
            self.pos_mask >>= 1;
        } else {
            // at boundary: move to upper block
            self.pos_mask = T::one() << (word_size::<T>() - 1);
            self.block_idx -= 1;
        }
    }

    fn try_prepare_left(&mut self) -> bool {
        if self.left_adj_mask != T::zero() {
            // simple case: not at block boundary
            if self.left_block.mv & self.pos_mask != T::zero() {
                self.left_block.dist -= 1;
                return true;
            }
        } else if let Some(b) = self.left_col.get(self.left_block_idx + 1) {
            // more complicated: at lower block boundary, and there is another block below in the band;
            if b.mv & T::one() == T::one() {
                // move one block down again
                let d = self.left_block.dist - 1;
                self.left_block = *b;
                self.left_block.dist = d;
                return true;
            }
        }
        false
    }

    #[inline]
    fn finish_move_left(&mut self) {
        self.col = self.left_col;
        self.left_col = self.states_iter.next().unwrap();
        self.block = replace(&mut self.left_block, self.left_col[self.left_block_idx]);
        self.left_block.adjust_up_by(self.left_adj_mask);
    }

    #[inline]
    fn done(&self) -> bool {
        self.pos_mask == T::zero() && self.block_idx == 0
    }

    fn print_state(&self) {
        eprintln!(
            "--- TB dist ({:?} <-> {:?})",
            self.left_block.dist, self.block.dist
        );
        eprintln!(
            " lm {:0width$b}\nleft block={} {}\n  m {:0width$b}\ncurrent block={} {}\n",
            self.left_adj_mask,
            self.left_block_idx,
            self.left_block,
            self.pos_mask,
            self.block_idx,
            self.block,
            width = word_size::<T>()
        );
    }
}

impl_myers!(
    usize,
    // maximum allowed value for max_dist (see States::step)
    usize::MAX - crate::pattern_matching::myers::word_size::<T>(),
    Myers<T>,
    crate::pattern_matching::myers::long::States<T>,
    crate::pattern_matching::myers::long::LongStatesHandler<'a>
);

#[cfg(test)]
mod tests {
    impl_common_tests!(true, super, u8, usize, build_64);

    #[test]
    fn test_myers_long_overflow() {
        let pattern = b"AAGACGAGAAAAGAAAGTCTAAAGGACTTTTGTGGCAAGACCATCCCTGTTCCCAACCCGACCCCTGGACCTCCCGCCCCGGGCACTCCCGACCCCCCGACCCCCCGACTCCTGGACCAGGAGACTGA";
        let text = b"GGCAAGGGGGACTGTAGATGGGTGAAAAGAGCAGTCAGGGACCAGGTCCTCAGCCCCCCAGCCCCCCAGCCCTCCAGGTCCCCAGCCCTCCAGGTCCCCAGCCCAACCCTTGTCCTTACCAGAACGTTGTTTTCAGGAAGTCTGAAAGACAAGAGCAGAAAGTCAGTCCCATGGAATTTTCGCTTCCCACAG".to_vec();

        let myers: Myers<u64> = Myers::new(pattern.iter().cloned());

        let hits: Vec<_> = myers.find_all_end(text, usize::MAX).collect();
        dbg!(hits);
    }
}
