use super::*;

pub(super) struct Traceback<T: BitVec> {
    states: Vec<State<T>>,
    positions: iter::Cycle<Range<usize>>,
    pos: usize,
    m: T::DistType,
}

impl<T> Traceback<T>
where
    T: BitVec,
{
    pub fn new() -> Traceback<T> {
        Traceback {
            states: vec![],
            positions: (0..0).cycle(),
            m: T::DistType::zero(),
            pos: 0,
        }
    }

    pub fn init(&mut self, initial_state: State<T>, num_cols: usize, m: T::DistType) {
        // Correct traceback needs two additional columns at the left of the matrix (see below).
        // Therefore reserving additional space.
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
    pub fn add_state(&mut self, s: State<T>) {
        self.pos = self.positions.next().unwrap();
        //self.states[self.pos] = s;
        // faster and safe since the positions are cycled.
        // Unless the code in init() is wrong, the index should
        // never be out of bounds.
        *unsafe { self.states.get_unchecked_mut(self.pos) } = s;
    }

    /// Returns the length of the current match, optionally adding the
    /// alignment path to `ops`
    pub fn traceback(&self, ops: Option<&mut Vec<AlignmentOperation>>) -> (usize, T::DistType) {
        self._traceback_at(self.pos, ops)
    }

    /// Returns the length of a match with a given end position, optionally adding the
    /// alignment path to `ops`
    /// only to be called if the `states` vec contains all states of the text
    pub fn traceback_at(
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

    /// returns a tuple of alignment length and hit distance, optionally adding the alignment path
    /// to `ops`
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
        let mut states = self.states[..=pos]
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
        // Note: Substitutions are always preferred over insertions,
        // but the order of the Ins / Del blocks could easily be swapped
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
                v_offset += 1;
                move_state_up!(state, current_pos);
                current_pos >>= 1;
                move_state_up!(lstate, current_pos);
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
        let states = self.states[..=pos]
            .iter()
            .rev()
            .chain(self.states.iter().rev().cycle());

        let m = self.m.to_usize().unwrap();

        let mut out: Vec<_> = (0..=m).map(|_| vec![]).collect();
        for s in states {
            let mut current_pos = T::one() << (m - 1);
            let mut dist = s.dist;
            for i in 0..=m {
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
