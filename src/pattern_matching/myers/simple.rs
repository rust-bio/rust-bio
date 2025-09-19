use std::borrow::Borrow;
use std::collections::HashMap;
use std::iter;
use std::marker::PhantomData;
use std::mem::{replace, size_of};
use std::slice;
use u64;

use num_traits::{FromPrimitive, One, ToPrimitive, Zero};

use crate::pattern_matching::myers::traceback::{StatesHandler, TracebackHandler};
use crate::pattern_matching::myers::{word_size, BitVec, State};

/// Myers algorithm.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Myers<T = u64>
where
    T: BitVec,
{
    pub(crate) peq: [T; 256],
    pub(crate) bound: T,
    pub(crate) m: T::DistType,
    pub(crate) states_store: Vec<State<T, T::DistType>>,
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
        let maxsize = T::DistType::from_usize(size_of::<T>() * 8).unwrap();
        let pattern = pattern.into_iter();
        let m = T::DistType::from_usize(pattern.len()).unwrap();
        assert!(m <= maxsize, "Pattern too long");
        assert!(m > T::DistType::zero(), "Pattern is empty");

        let mut peq = [T::zero(); 256];

        for (i, symbol) in pattern.enumerate() {
            let symbol = *symbol.borrow();
            let mask = T::one() << i;
            // equivalent
            peq[symbol as usize] |= mask;
            // ambiguities
            if let Some(equivalents) = opt_ambigs.and_then(|ambigs| ambigs.get(&symbol)) {
                for &eq in equivalents {
                    peq[eq as usize] |= mask;
                }
            }
        }

        if let Some(wildcards) = opt_wildcards {
            for &w in wildcards {
                peq[w as usize] = T::max_value();
            }
        }

        Myers {
            peq,
            bound: T::one() << (m.to_usize().unwrap() - 1),
            m,
            states_store: vec![],
        }
    }

    #[inline]
    fn initial_state(&self, m: T::DistType, _: T::DistType) -> State<T, T::DistType> {
        State::init(m)
    }

    #[inline]
    fn step(&self, state: &mut State<T, T::DistType>, a: u8, _: T::DistType) {
        self._step(state, a);
    }

    #[inline]
    fn _step(&self, state: &mut State<T, T::DistType>, a: u8) {
        let eq = self.peq[a as usize];
        let xv = eq | state.mv;
        let xh = ((eq & state.pv).wrapping_add(&state.pv) ^ state.pv) | eq;

        let mut ph = state.mv | !(xh | state.pv);
        let mut mh = state.pv & xh;

        // if ph & self.bound > T::zero() {
        //     state.dist += T::DistType::one();
        // } else if mh & self.bound > T::zero() {
        //     state.dist -= T::DistType::one();
        // }
        let diff = ((ph & self.bound) != T::zero()) as i8 - ((mh & self.bound) != T::zero()) as i8;
        state.dist =
            T::DistType::from_usize(state.dist.to_usize().unwrap().wrapping_add(diff as usize))
                .unwrap();

        ph <<= 1;
        mh <<= 1;
        state.pv = mh | !(xv | ph);
        state.mv = ph & xv;
    }

    #[inline]
    pub fn m(&self) -> T::DistType {
        self.m
    }
}

#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub(super) struct ShortStatesHandler<'a>(PhantomData<&'a ()>);

impl<'a> ShortStatesHandler<'a> {
    #[inline]
    pub fn new() -> Self {
        ShortStatesHandler(PhantomData)
    }
}

impl<'a, T: BitVec + 'a> StatesHandler<'a, T, T::DistType> for ShortStatesHandler<'a> {
    type TracebackHandler = ShortTracebackHandler<'a, T>;
    type TracebackColumn = State<T, T::DistType>;

    #[inline]
    fn init(&mut self, n: usize, _: T::DistType) -> usize {
        n
    }

    #[inline]
    fn n_blocks(&self) -> usize {
        1
    }

    #[inline]
    fn set_max_state(&self, pos: usize, states: &mut [State<T, T::DistType>]) {
        states[pos] = State::init_max_dist();
    }

    #[inline]
    fn add_state(
        &self,
        source: &Self::TracebackColumn,
        pos: usize,
        states: &mut [State<T, T::DistType>],
    ) {
        states[pos] = *source;
    }

    #[inline]
    fn dist_at(&self, pos: usize, states: &[State<T, T::DistType>]) -> Option<T::DistType> {
        states.get(pos).map(|s| s.dist)
    }

    #[inline]
    fn init_traceback(
        &self,
        m: T::DistType,
        pos: usize,
        states: &'a [State<T, T::DistType>],
    ) -> Option<Self::TracebackHandler> {
        Some(ShortTracebackHandler::new(m, pos, states))
    }
}

type RevColIter<'a, T> = iter::Rev<slice::Iter<'a, State<T, <T as BitVec>::DistType>>>;

pub(super) struct ShortTracebackHandler<'a, T: BitVec> {
    states_iter: iter::Chain<RevColIter<'a, T>, iter::Cycle<RevColIter<'a, T>>>,
    // block representing the current TB matrix column
    block: State<T, T::DistType>,
    // block representing the TB matrix column to the left
    left_block: State<T, T::DistType>,
    // mask with one active bit indicating the position of the last letter
    // TODO: not needed with padding
    max_mask: T,
    // mask with one active bit indicating the current position (row) in the current block
    pos_mask: T,
    // Bit mask that "covers" the range of positions *below* the left cell
    // in the DP matrix. This is used used to initialize blocks
    // (adjust the distance score) with `State::adjust_up_by(left_adj_mask)`.
    left_adj_mask: T,
    _a: PhantomData<&'a ()>,
}

impl<'a, T: BitVec> ShortTracebackHandler<'a, T> {
    #[inline]
    fn new(m: T::DistType, pos: usize, states: &'a [State<T, T::DistType>]) -> Self {
        let pos_mask = T::one() << (m.to_usize().unwrap() - 1);

        // reverse iterator over DP matrix columns
        // chain + cycle is needed because `find_all` (`FullMatches` iterator)
        // only stores m + max_dist states, and progressively cycles through the `states`.
        // `LazyMatches` stores all states (cycle not actually needed, but performance loss small).
        let mut states_iter = states[..=pos]
            .iter()
            .rev()
            .chain(states.iter().rev().cycle());

        // // Simpler alternative using skip() is slower in some cases:
        // let mut states = states.iter().rev().cycle().skip(states.len() - pos - 1);

        let state = *states_iter.next().unwrap();

        let mut left_state = *states_iter.next().unwrap();
        left_state.adjust_one_up(pos_mask);
        // *note*: this mask is initialized for the left state to be in the diagonal
        let left_mask = pos_mask;

        ShortTracebackHandler {
            block: state,
            left_block: left_state,
            states_iter,
            max_mask: pos_mask,
            pos_mask,
            left_adj_mask: left_mask,
            _a: PhantomData,
        }
    }

    #[inline]
    fn left_mask_up(&mut self) {
        self.left_adj_mask = (self.left_adj_mask >> 1) | self.max_mask;
        // println!("left {:064b}\nmax  {:064b}", self.left_mask, self.max_mask);
    }
}

impl<'a, T> TracebackHandler<'a, T, T::DistType> for ShortTracebackHandler<'a, T>
where
    T: BitVec + 'a,
{
    #[inline]
    fn dist(&self) -> T::DistType {
        self.block.dist
    }

    #[inline]
    fn left_dist(&self) -> T::DistType {
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
        self.block.adjust_one_up(self.pos_mask);
        self.pos_mask >>= 1;
        // (2) move up in left column
        self.left_mask_up();
        self.left_block.adjust_one_up(self.pos_mask);
    }

    #[inline]
    fn try_prepare_left(&mut self) -> bool {
        if self.left_block.mv & self.pos_mask != T::zero() {
            self.left_block.dist -= T::DistType::one();
            return true;
        }
        false
    }

    #[inline]
    fn prepare_diagonal(&mut self) {
        self.left_mask_up();
        self.pos_mask >>= 1;
    }

    #[inline]
    fn finish_move_left(&mut self) {
        self.block = replace(&mut self.left_block, *self.states_iter.next().unwrap());
        self.left_block.adjust_up_by(self.left_adj_mask);
    }

    #[inline]
    fn done(&self) -> bool {
        self.pos_mask == T::zero()
    }

    fn print_state(&self) {
        eprintln!(
            "--- TB dist ({:?} <-> {:?})",
            self.left_block.dist, self.block.dist
        );
        eprintln!(
            " lm {:0width$b}\nleft {}\n  m {:0width$b}\ncurrent {}\n",
            self.left_adj_mask,
            self.left_block,
            self.pos_mask,
            self.block,
            width = word_size::<T>()
        );
    }
}

impl_myers!(
    T::DistType,
    T::DistType::max_value(),
    Myers<T>,
    crate::pattern_matching::myers::State<T, T::DistType>,
    crate::pattern_matching::myers::simple::ShortStatesHandler<'a>
);
