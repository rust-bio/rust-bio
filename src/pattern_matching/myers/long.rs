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
use std::u64;

use itertools::Itertools;
use num_traits::ToPrimitive;

use crate::pattern_matching::myers::traceback::{StatesHandler, TracebackHandler};
use crate::pattern_matching::myers::{ceil_div, State};

use super::word_size;
use super::BitVec;

struct Peq<T: BitVec> {
    peq: [T; 256],
    bound: T,
}

/// Myers algorithm.
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
        assert!(m <= usize::max_value() / 2, "Pattern too long");

        // build peq
        let mut peq = vec![];
        for chunk in pattern.chunks(w).into_iter() {
            let mut peq_block = [T::zero(); 256];
            let mut i = 0;
            for symbol in chunk {
                let symbol = *symbol.borrow();
                let mask = T::one() << i;
                // equivalent
                peq_block[symbol as usize] |= mask;
                // ambiguities
                if let Some(equivalents) = opt_ambigs.and_then(|ambigs| ambigs.get(&symbol)) {
                    for &eq in equivalents {
                        peq_block[eq as usize] |= mask;
                    }
                }
                i += 1;
            }
            // wildcards
            if let Some(wildcards) = opt_wildcards {
                for &w in wildcards {
                    peq_block[w as usize] = T::max_value();
                }
            }

            peq.push(Peq {
                peq: peq_block,
                bound: T::one() << (i - 1),
            });
        }

        Myers {
            peq,
            m,
            states_store: vec![],
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

    // apparently faster than uncommented code above
    let mut hout = ((ph & p.bound) != T::zero()) as i8;
    hout -= ((mh & p.bound) != T::zero()) as i8;
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

pub(super) struct States<T: BitVec> {
    states: Vec<State<T, usize>>,
    max_block: usize,
    last_m: usize,
}

impl<T> States<T>
where
    T: BitVec,
{
    fn new(m: usize, max_dist: usize) -> Self {
        let w = word_size::<T>();
        let mut s = States {
            states: vec![],
            max_block: ceil_div(m, w) - 1,
            last_m: m % w,
        };
        let min_blocks = max(1, ceil_div(min(max_dist, m), w));
        for _ in 0..min_blocks {
            s.add_state(0);
        }
        s
    }

    #[inline]
    fn add_state(&mut self, offset: i8) {
        let prev_dist = self.states.last().map(|s| s.dist).unwrap_or(0);
        // delta: number of bits of the new state covered by the pattern
        // For the last block, we add m % w, not all bits are necessarily used.
        // This strategy differs from the solution by Myers (p. 407, note 4).
        // We wanted to avoid having to pad pattern and sequence.
        let delta = if self.states.len() == self.max_block && self.last_m > 0 {
            self.last_m
        } else {
            word_size::<T>()
        };
        self.states.push(State::init(
            (prev_dist)
                .wrapping_add(delta)
                .wrapping_add(offset as usize)
                .to_usize()
                .unwrap(),
        ));
    }

    #[inline]
    fn step(&mut self, a: u8, peq: &[Peq<T>], max_dist: usize) {
        let mut carry = 0;
        let mut last_block = self.states.len() - 1;

        for (state, block_peq) in self.states.iter_mut().zip(peq) {
            carry = advance_block(state, block_peq, a, carry);
        }

        let w = word_size::<T>();
        let last_dist = self.states[last_block].dist;
        if (last_dist as isize - carry as isize) as usize <= max_dist
            && last_block < self.max_block
            && (peq[last_block + 1].peq[a as usize] & T::one() == T::one() || carry < 0)
        {
            last_block += 1;
            self.add_state(-carry as i8);
            advance_block(&mut self.states[last_block], &peq[last_block], a, carry);
        } else {
            while last_block > 0 && self.states[last_block].dist >= max_dist + w {
                last_block -= 1;
            }
            self.states.truncate(last_block + 1);
        }
    }

    /// Returns the last distance score of the traceback column if known
    /// (only if all blocks were computed).
    #[inline]
    fn known_dist(&self) -> Option<usize> {
        self.states.get(self.max_block).map(|s| s.dist)
    }
}

#[derive(Default)]
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
    fn set_max_state(&self, pos: usize, states: &mut [State<T, usize>]) {
        let pos = pos * self.n_blocks;
        for s in states.iter_mut().skip(pos).take(self.n_blocks) {
            *s = State::max();
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
            // When following the traceback path, it can happen that the block to the
            // left was not computed because it is outside of the band.
            // In order to prevent the algorithm from going left in this case,
            // we initialize the block below the last computed block with meaningful
            // defaults.
            states[pos + source.len()] = State {
                dist: usize::max_value(),
                pv: T::zero(),
                mv: T::zero(),
            };
        }
    }

    #[inline]
    fn init_traceback(
        &self,
        m: usize,
        pos: usize,
        states: &'a [State<T, usize>],
    ) -> Self::TracebackHandler {
        LongTracebackHandler::new(self.n_blocks, m, pos, states)
    }
}

type RevColIter<'a, T> = iter::Rev<slice::Chunks<'a, State<T, usize>>>;

pub(super) struct LongTracebackHandler<'a, T: BitVec> {
    states_iter: iter::Chain<RevColIter<'a, T>, iter::Cycle<RevColIter<'a, T>>>,
    block_pos: usize,
    left_block_pos: usize,
    col: &'a [State<T, usize>],
    left_col: &'a [State<T, usize>],
    block: State<T, usize>,
    left_block: State<T, usize>,
    left_max_mask: T,
    pos_bitvec: T,
    left_mask: T,
    _a: PhantomData<&'a ()>,
}

impl<'a, T: BitVec> LongTracebackHandler<'a, T> {
    #[inline]
    fn new(n_blocks: usize, m: usize, pos: usize, states: &'a [State<T, usize>]) -> Self {
        let mut last_m = m.to_usize().unwrap() % word_size::<T>();
        if last_m == 0 {
            last_m = word_size::<T>();
        }
        let mask0 = T::one() << (last_m - 1);

        let pos = n_blocks * (pos + 1);
        let mut states_iter = states[..pos]
            .chunks(n_blocks)
            .rev()
            .chain(states.chunks(n_blocks).rev().cycle());

        let col = states_iter.next().unwrap();
        let left_col = states_iter.next().unwrap();

        // This bit mask is supplied to State::adjust_by_mask() in order to adjust the distance
        // of the left block. It is adjusted with every `move_up_left`
        let left_mask = if last_m != 1 {
            T::zero()
        } else {
            T::from_usize(0b10).unwrap()
        };

        LongTracebackHandler {
            block_pos: n_blocks - 1,
            left_block_pos: n_blocks - 1,
            block: col.last().unwrap().clone(),
            left_block: left_col.last().unwrap().clone(),
            col,
            left_col,
            states_iter,
            pos_bitvec: mask0,
            left_mask,
            left_max_mask: mask0,
            _a: PhantomData,
        }
    }
}

impl<'a, T: BitVec + 'a> TracebackHandler<'a, T, usize> for LongTracebackHandler<'a, T> {
    #[inline]
    fn block(&self) -> &State<T, usize> {
        &self.block
    }

    #[inline]
    fn block_mut(&mut self) -> &mut State<T, usize> {
        &mut self.block
    }

    #[inline]
    fn left_block(&self) -> &State<T, usize> {
        &self.left_block
    }

    #[inline]
    fn left_block_mut(&mut self) -> &mut State<T, usize> {
        &mut self.left_block
    }

    #[inline]
    fn pos_bitvec(&self) -> T {
        self.pos_bitvec
    }

    #[inline]
    fn move_up(&mut self, adjust_dist: bool) {
        // If the block boundary has not been reached yet, we can shift to the
        // upper position. Otherwise, we move to a new block (if there is one!)
        if self.pos_bitvec != T::one() || self.block_pos == 0 {
            if adjust_dist {
                self.block.adjust_dist(self.pos_bitvec);
            }
            self.pos_bitvec >>= 1;
        } else {
            // move to upper block
            self.pos_bitvec = T::one() << (word_size::<T>() - 1);
            self.block_pos -= 1;
            if adjust_dist {
                self.block = self.col[self.block_pos].clone();
            }
        }
    }

    #[inline]
    fn move_up_left(&mut self, adjust_dist: bool) {
        // If the block boundary has not been reached yet, we can extend the range mask by
        // activating a new bit.
        // However, we switch to a new block (if there is one!) before the mask would cover
        // the whole block.
        if self.left_mask & T::from_usize(0b10).unwrap() == T::zero() || self.left_block_pos == 0 {
            self.left_mask = (self.left_mask >> 1) | self.left_max_mask;
            if adjust_dist {
                self.left_block.adjust_dist(self.pos_bitvec);
            }
        } else {
            self.left_max_mask = T::one() << (word_size::<T>() - 1);
            self.left_mask = T::zero();
            self.left_block_pos -= 1;
            if adjust_dist {
                self.left_block = self.left_col[self.left_block_pos].clone();
            }
        }
    }

    #[inline]
    fn move_to_left(&mut self) {
        self.col = self.left_col;
        self.left_col = self.states_iter.next().unwrap();
        self.block = replace(
            &mut self.left_block,
            self.left_col[self.left_block_pos].clone(),
        );
        self.left_block.adjust_by_mask(self.left_mask);
    }

    #[inline]
    fn move_left_down_if_better(&mut self) -> bool {
        if self.left_mask != T::zero() {
            // simple case: not at block boundary
            if self.left_block.mv & self.pos_bitvec != T::zero() {
                self.left_block.dist -= 1;
                return true;
            }
        } else if let Some(b) = self.left_col.get(self.left_block_pos + 1) {
            // more complicated: at lower block boundary, and there is a lower block
            if b.mv & T::one() == T::one() {
                let d = self.left_block().dist - 1;
                self.left_block = b.clone();
                self.left_block.dist = d;
                return true;
            }
        }
        false
    }

    #[inline]
    fn column_slice(&self) -> &'a [State<T, usize>] {
        self.col
    }

    #[inline]
    fn finished(&self) -> bool {
        self.pos_bitvec == T::zero() && self.block_pos == 0
    }
}

impl_myers!(
    usize,
    Myers<T>,
    crate::pattern_matching::myers::long::States<T>,
    crate::pattern_matching::myers::long::LongStatesHandler<'a>
);

#[cfg(test)]
mod tests {
    use super::*;

    impl_tests!(super, u8, usize, build_64);

    #[test]
    fn test_myers_long_overflow() {
        let pattern = b"AAGACGAGAAAAGAAAGTCTAAAGGACTTTTGTGGCAAGACCATCCCTGTTCCCAACCCGACCCCTGGACCTCCCGCCCCGGGCACTCCCGACCCCCCGACCCCCCGACTCCTGGACCAGGAGACTGA";
        let text = b"GGCAAGGGGGACTGTAGATGGGTGAAAAGAGCAGTCAGGGACCAGGTCCTCAGCCCCCCAGCCCCCCAGCCCTCCAGGTCCCCAGCCCTCCAGGTCCCCAGCCCAACCCTTGTCCTTACCAGAACGTTGTTTTCAGGAAGTCTGAAAGACAAGAGCAGAAAGTCAGTCCCATGGAATTTTCGCTTCCCACAG".to_vec();

        let myers: Myers<u64> = Myers::new(pattern.iter().cloned());

        let hits: Vec<_> = myers.find_all_end(text, usize::max_value() - 64).collect();
        dbg!(hits);
    }
}
