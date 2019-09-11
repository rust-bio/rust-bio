use std::borrow::Borrow;
use std::collections::HashMap;
use std::iter;
use std::marker::PhantomData;
use std::mem::{replace, size_of};
use std::slice;
use std::u64;

use num_traits::{FromPrimitive, One, ToPrimitive, Zero};

use crate::pattern_matching::myers::traceback::{StatesHandler, TracebackHandler};
use crate::pattern_matching::myers::{BitVec, State};

/// Myers algorithm.
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

#[derive(Default)]
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
    fn set_max_state(&self, pos: usize, states: &mut [State<T, T::DistType>]) {
        //states[pos] = State::max();
        *unsafe { states.get_unchecked_mut(pos) } = State::max();
    }

    #[inline]
    fn add_state(
        &self,
        source: &Self::TracebackColumn,
        pos: usize,
        states: &mut [State<T, T::DistType>],
    ) {
        //states[pos] = source.clone();
        *unsafe { states.get_unchecked_mut(pos) } = source.clone();
    }

    #[inline]
    fn init_traceback(
        &self,
        m: T::DistType,
        pos: usize,
        states: &'a [State<T, T::DistType>],
    ) -> Self::TracebackHandler {
        ShortTracebackHandler::new(m, pos, states)
    }
}

type RevColIter<'a, T> = iter::Rev<slice::Iter<'a, State<T, <T as BitVec>::DistType>>>;

pub(super) struct ShortTracebackHandler<'a, T: BitVec> {
    states_iter: iter::Chain<RevColIter<'a, T>, iter::Cycle<RevColIter<'a, T>>>,
    state: State<T, T::DistType>,
    left_state: State<T, T::DistType>,
    max_mask: T,
    pos_bitvec: T,
    left_mask: T,
    _a: PhantomData<&'a ()>,
}

impl<'a, T: BitVec> ShortTracebackHandler<'a, T> {
    #[inline]
    fn new(m: T::DistType, pos: usize, states: &'a [State<T, T::DistType>]) -> Self {
        let mask0 = T::one() << (m.to_usize().unwrap() - 1);

        // Reverse iterator over states. If remembering all positions,
        // the chain() and cycle() are not actually needed, but there seems
        // to be almost no performance loss.
        let mut states_iter = states[..=pos]
            .iter()
            .rev()
            .chain(states.iter().rev().cycle());

        // // Simpler alternative using skip() is slower in some cases:
        // let mut states = states.iter().rev().cycle().skip(states.len() - pos - 1);

        ShortTracebackHandler {
            state: states_iter.next().unwrap().clone(),
            left_state: states_iter.next().unwrap().clone(),
            states_iter,
            max_mask: mask0,
            pos_bitvec: mask0,
            left_mask: T::zero(),
            _a: PhantomData,
        }
    }
}

impl<'a, T> TracebackHandler<'a, T, T::DistType> for ShortTracebackHandler<'a, T>
where
    T: BitVec + 'a,
{
    #[inline]
    fn block(&self) -> &State<T, T::DistType> {
        &self.state
    }

    #[inline]
    fn block_mut(&mut self) -> &mut State<T, T::DistType> {
        &mut self.state
    }

    #[inline]
    fn left_block(&self) -> &State<T, T::DistType> {
        &self.left_state
    }

    #[inline]
    fn left_block_mut(&mut self) -> &mut State<T, T::DistType> {
        &mut self.left_state
    }

    #[inline]
    fn pos_bitvec(&self) -> T {
        self.pos_bitvec
    }

    #[inline]
    fn move_up(&mut self, adjust_dist: bool) {
        if adjust_dist {
            self.state.adjust_dist(self.pos_bitvec);
        }
        self.pos_bitvec >>= 1;
    }

    #[inline]
    fn move_up_left(&mut self, adjust_dist: bool) {
        self.left_mask = (self.left_mask >> 1) | self.max_mask;
        if adjust_dist {
            self.left_state.adjust_dist(self.pos_bitvec);
        }
    }

    #[inline]
    fn move_to_left(&mut self) {
        self.state = replace(
            &mut self.left_state,
            self.states_iter.next().unwrap().clone(),
        );
        self.left_state.adjust_by_mask(self.left_mask);
    }

    #[inline]
    fn move_left_down_if_better(&mut self) -> bool {
        if self.left_state.mv & self.pos_bitvec != T::zero() {
            self.left_state.dist -= T::DistType::one();
            return true;
        }
        false
    }

    #[inline]
    fn column_slice(&self) -> &'a [State<T, T::DistType>] {
        unsafe { std::slice::from_raw_parts(&self.state, 1) }
    }

    #[inline]
    fn finished(&self) -> bool {
        self.pos_bitvec == T::zero()
    }
}

impl_myers!(
    T::DistType,
    Myers<T>,
    crate::pattern_matching::myers::State<T, T::DistType>,
    crate::pattern_matching::myers::simple::ShortStatesHandler<'a>
);
