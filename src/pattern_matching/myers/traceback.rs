use std::default::Default;
use std::iter;
use std::marker::PhantomData;
use std::ops::Range;

use crate::alignment::AlignmentOperation;

use crate::pattern_matching::myers::{word_size, BitVec, DistType, State};

/// Objects implementing this trait handle the addition of calculated blocks (State<T, D>)
/// to a container, and are responsible for creating the respective `TracebackHandler` object.
pub(super) trait StatesHandler<'a, T, D>
where
    T: BitVec + 'a,
    D: DistType,
{
    /// Object that helps obtaining a single traceback path
    type TracebackHandler: TracebackHandler<'a, T, D>;
    /// Type that represents a column in the traceback matrix
    type TracebackColumn: ?Sized;

    /// Prepare for a new search given n (maximum expected number of traceback columns) and
    /// m (pattern length).
    /// Returns the expected size of the vector storing the calculated blocks given this
    /// information. The vector will then be initialized with the given number of 'empty'
    /// State<T, D> objects and supplied to the other methods as slice.
    fn init(&mut self, n: usize, m: D) -> usize;

    /// Fill the column at `pos` with states initialized with the maximum distance
    /// (`State::max()`).
    fn set_max_state(&self, pos: usize, states: &mut [State<T, D>]);

    /// This method copies over all blocks (or the one block) from a tracback column
    /// into the mutable `states` slice at the given column position.
    fn add_state(&self, source: &Self::TracebackColumn, pos: usize, states: &mut [State<T, D>]);

    /// Initiates a `TracebackHandler` object to assist with a traceback, 'starting'
    /// at the given end position.
    fn init_traceback(&self, m: D, pos: usize, states: &'a [State<T, D>])
        -> Self::TracebackHandler;
}

/// Objects implementing this trait should store states and have methods
/// necessary for obtaining a single traceback path. This allows to use the
/// same traceback code for the simple and the block-based Myers pattern
/// matching approaches. It is designed to be as general as possible
/// to allow different implementations.
///
/// Implementors of `TracebackHandler` keep two `State<T, D>` instances,
/// which store the information from two horizontally adjacent traceback
/// columns, encoded in the PV / MV bit vectors. The columns are accessible
/// using the methods `block()` (current / right column) and `left_block()`
/// (left column). Moving horizontally to the next position can be achieved
/// using `move_left()`.
///
/// Implementors also track the vertical cursor positions within the current
/// traceback columns (two separate cursors for left and right column).
/// `block()` and `left_block()` will always return the block that currently
/// contain the cursors.
/// `pos_bitvec()` returns a bit vector with a single activated bit at the current
/// vertical position within the *right (current)* column.
/// Moving to the next vertical position is achieved by `move_up()` and
/// `move_up_left()`. With the block based implementation, this may involve
/// switching to a new block.
pub(super) trait TracebackHandler<'a, T, D>
where
    T: BitVec + 'a,
    D: DistType,
{
    /// Returns a reference to the current (right) block.
    fn block(&self) -> &State<T, D>;

    /// Returns a mutable reference to the current (right) block.
    fn block_mut(&mut self) -> &mut State<T, D>;

    /// Returns a reference to the left block.
    fn left_block(&self) -> &State<T, D>;

    /// Returns a mutable reference to the left block.
    fn left_block_mut(&mut self) -> &mut State<T, D>;

    /// Bit vector representing the position in the traceback. Only the bit
    /// at the current position should be on.
    /// For a search pattern of length 4, the initial bit vector would be
    /// `0b1000`. A call to `move_up_cursor()` will shift the vector, so another
    /// call to `pos_bitvec()` results in `0b100`.
    /// The bit vector has a width of `T`, meaning that it can store
    /// the same number of positions as the PV and MV vectors. In the
    /// case of the block based algorithm, the vector only stores the
    /// position within the current block.
    fn pos_bitvec(&self) -> T;

    /// Move up cursor by one position in traceback matrix.
    ///
    /// # Arguments
    ///
    /// * adjust_dist: If true, the distance score of the block is adjusted
    ///   based on the current cursor position before moving it up.
    ///  *Note concerning the block based Myers algorithm:*
    ///  The the active bit in bit vector returned by `pos_bitvec()`
    ///  is expected to jump back to the maximum (lowest) position
    ///  when reaching the uppermost position (like `rotate_right()` does).
    fn move_up(&mut self, adjust_dist: bool);

    /// Move up left cursor by one position in traceback matrix.
    ///
    /// # Arguments
    ///
    /// * adjust_dist: If true, the distance score of the block is adjusted
    ///   based on the current cursor position before moving it up.
    ///   However, the current cursor position of the **right** block is used,
    ///   **not** the one of the left block. This is an important oddity, which
    ///   makes only sense because of the design of the traceback algorithm.
    fn move_up_left(&mut self, adjust_dist: bool);

    /// Shift the view by one traceback column / block to the left. The
    /// block that was on the left position previously moves to the right /
    /// current block without changes. The cursor positions have to be
    /// adjusted indepentedently if necessary using `move_up(false)` /
    /// `move_up_left(false)`.
    /// `move_left()` adjusts distance score of the new left block to
    /// be correct for the left vertical cursor position. It is therefore
    /// important that the cursor is moved *before* calling `move_left()`.
    fn move_to_left(&mut self);

    /// Rather specialized method that allows having a simpler code in Traceback::_traceback_at()
    /// Checks if the position below the left cursor has a smaller distance, and if so,
    /// moves the cursor to this block and returns `true`.
    ///
    /// The problem is that the current implementation always keeps the left cursor in the
    /// diagonal position for performance reasons. In this case, checking the actual left
    /// distance score can be complicated with the block-based algorithm since the left cursor
    /// may be at the lower block boundary. If so, the function thus has to check the topmost
    /// position of the lower block and keep this block if the distance is better (lower).
    fn move_left_down_if_better(&mut self) -> bool;

    /// Returns a slice containing all blocks of the current traceback column
    /// from top to bottom. Used for debugging only.
    fn column_slice(&self) -> &[State<T, D>];

    /// Returns true if topmost position in the traceback matrix has been reached,
    /// meaning that the traceback is complete.
    /// Technically this means, that `move_up_cursor()` was called so many times
    /// until the uppermost block was reached and the pos_bitvec() does not contain
    /// any bit, since shifting has removed it from the vector.
    fn finished(&self) -> bool;

    /// For debugging only
    fn print_state(&self) {
        println!(
            "--- TB dist ({:?} <-> {:?})",
            self.left_block().dist,
            self.block().dist
        );
        println!(
            "{:064b} m\n{:064b} + ({:?}) (left) d={:?}\n{:064b} - ({:?})\n \
             {:064b} + ({:?}) (current) d={:?}\n{:064b} - ({:?})\n",
            self.pos_bitvec(),
            self.left_block().pv,
            self.left_block().pv,
            self.left_block().dist,
            self.left_block().mv,
            self.left_block().mv,
            self.block().pv,
            self.block().pv,
            self.block().dist,
            self.block().mv,
            self.block().mv
        );
    }
}

pub(super) struct Traceback<'a, T, D, H>
where
    T: BitVec + 'a,
    D: DistType,
    H: StatesHandler<'a, T, D>,
{
    m: D,
    positions: iter::Cycle<Range<usize>>,
    handler: H,
    pos: usize,
    _t: PhantomData<&'a T>,
}

impl<'a, T, D, H> Traceback<'a, T, D, H>
where
    T: BitVec,
    D: DistType,
    H: StatesHandler<'a, T, D>,
{
    #[inline]
    pub fn new(
        states: &mut Vec<State<T, D>>,
        initial_state: &H::TracebackColumn,
        num_cols: usize,
        m: D,
        mut handler: H,
    ) -> Self {
        // Correct traceback needs two additional columns at the left of the matrix (see below).
        // Therefore reserving additional space.
        let num_cols = num_cols + 2;

        let n_states = handler.init(num_cols, m);

        let mut tb = Traceback {
            m,
            positions: (0..num_cols).cycle(),
            handler,
            pos: 0,
            _t: PhantomData,
        };

        // extend or truncate states vector
        let curr_len = states.len();
        if n_states > curr_len {
            states.reserve(n_states);
            states.extend((0..n_states - curr_len).map(|_| State::default()));
        } else {
            states.truncate(n_states);
            states.shrink_to_fit();
        }
        // important if using unsafe in add_state(), and also for correct functioning of traceback
        debug_assert!(states.len() == n_states);

        // first column is used to ensure a correct path if the text (target)
        // is shorter than the pattern (query)
        tb.pos = tb.positions.next().unwrap();
        tb.handler.set_max_state(tb.pos, states);

        // initial state
        tb.add_state(initial_state, states);

        tb
    }

    #[inline]
    pub fn add_state(&mut self, column: &H::TracebackColumn, states: &mut [State<T, D>]) {
        self.pos = self.positions.next().unwrap();
        self.handler.add_state(column, self.pos, states);
    }

    /// Returns the length of the current match, optionally adding the
    /// alignment path to `ops`
    #[inline]
    pub fn traceback(
        &self,
        ops: Option<&mut Vec<AlignmentOperation>>,
        states: &'a [State<T, D>],
    ) -> (D, D) {
        self._traceback_at(self.pos, ops, states)
    }

    /// Returns the length of a match with a given end position, optionally adding the
    /// alignment path to `ops`
    /// only to be called if the `states` vec contains all states of the text
    #[inline]
    pub fn traceback_at(
        &self,
        pos: usize,
        ops: Option<&mut Vec<AlignmentOperation>>,
        states: &'a [State<T, D>],
    ) -> Option<(D, D)> {
        let pos = pos + 2; // in order to be comparable since self.pos starts at 2, not 0
        if pos <= self.pos {
            return Some(self._traceback_at(pos, ops, states));
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
        state_slice: &'a [State<T, D>],
    ) -> (D, D) {
        use self::AlignmentOperation::*;

        // Generic object that holds the necessary data and methods
        let mut h = self.handler.init_traceback(self.m, pos, state_slice);

        // self.print_tb_matrix(pos, state_slice);

        let ops = &mut ops;

        // horizontal column offset from starting point in traceback matrix (bottom right)
        let mut h_offset = D::zero();

        // distance of the match (will be returned)
        let dist = h.block().dist;

        // The cursor of the left state is always for diagonal position in the traceback matrix.
        // This allows checking for a substitution by a simple comparison.
        h.move_up_left(true);

        // Loop for finding the traceback path
        // If there are several possible solutions, substitutions are preferred over InDels
        // (Subst > Ins > Del)
        while !h.finished() {
            let op;
            // This loop is used to allow skipping `move_left()` using break (kind of similar
            // to 'goto'). This was done to avoid having to inline move_left() three times,
            // which would use more space.
            #[allow(clippy::never_loop)]
            loop {
                // h.print_state();

                if h.left_block().dist.wrapping_add(&D::one()) == h.block().dist {
                    // Diagonal (substitution)
                    // Since the left cursor is always in the upper diagonal position,
                    // a simple comparison of distances is enough to determine substitutions.
                    h.move_up(false);
                    h.move_up_left(false);
                    op = Subst;
                } else if h.block().pv & h.pos_bitvec() != T::zero() {
                    // Up
                    h.move_up(true);
                    h.move_up_left(true);
                    op = Ins;
                    break;
                } else if h.move_left_down_if_better() {
                    // Left
                    op = Del;
                } else {
                    // Diagonal (match)
                    h.move_up(false);
                    h.move_up_left(false);
                    op = Match;
                }

                // Moving one position to the left, adjusting h_offset
                h_offset += D::one();
                h.move_to_left();
                break;
            }

            // println!("{:?}", op);

            if let Some(o) = ops.as_mut() {
                o.push(op);
            }
        }

        (h_offset, dist)
    }

    // Useful for debugging
    #[allow(dead_code)]
    fn print_tb_matrix(&self, pos: usize, state_slice: &'a [State<T, D>]) {
        let mut h = self.handler.init_traceback(self.m, pos, state_slice);
        let m = self.m.to_usize().unwrap();
        let mut out = vec![];
        for _ in 0..state_slice.len() {
            let mut col_out = vec![];
            let mut empty = true;
            for (i, state) in h.column_slice().iter().enumerate().rev() {
                if !(state.is_new() || state.is_max()) {
                    empty = false;
                }
                let w = word_size::<T>();
                let end = (i + 1) * w;
                let n = if end <= m { w } else { m % w };
                state.write_dist_column(n, &mut col_out);
            }
            out.push(col_out);
            h.move_to_left();
            if empty {
                break;
            }
        }

        for j in (0..m).rev() {
            print!("{:>4}: ", m - j + 1);
            for col in out.iter().rev() {
                if let Some(d) = col.get(j) {
                    if *d >= (D::max_value() >> 1) {
                        // missing value
                        print!("    ");
                    } else {
                        print!("{:>4?}", d);
                    }
                } else {
                    print!("   -");
                }
            }
            println!();
        }
    }
}
