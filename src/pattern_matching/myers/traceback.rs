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
    /// Type that represents a column in the DP matrix
    type TracebackColumn: ?Sized;

    /// Prepare for a new search given n (maximum expected number of traceback columns) and
    /// m (pattern length).
    /// Returns the expected size of the vector storing the calculated blocks given this
    /// information. The vector will then be initialized with the given number of
    /// State<T, D> objects and supplied to the other methods as slice.
    fn init(&mut self, n: usize, m: D) -> usize;

    /// The max. number of blocks needed for a DP matrix column
    fn n_blocks(&self) -> usize;

    /// Fill the column at `pos` with states initialized with the maximum distance
    /// (`State::init_max_dist()`). The leftmost column of the traceback matrix needs this.
    fn set_max_state(&self, pos: usize, states: &mut [State<T, D>]);

    /// This method copies over all blocks (or the one block) from a tracback column
    /// into the mutable `states` slice at the given column position.
    fn add_state(&self, source: &Self::TracebackColumn, pos: usize, states: &mut [State<T, D>]);

    /// Returns the edit distance at the given position, or `None` if `pos` is
    /// out of bounds.
    fn dist_at(&self, pos: usize, states: &[State<T, D>]) -> Option<D>;

    /// Initiates a `TracebackHandler` object to assist with a traceback, 'starting'
    /// at the given end position.
    /// Should return `None` if the alignment cannot be obtained
    /// (can happen with block-based algorithm if not all blocks were computed)
    fn init_traceback(
        &self,
        m: D,
        pos: usize,
        states: &'a [State<T, D>],
    ) -> Option<Self::TracebackHandler>;
}

/// `TracebackHandler` objects assist with the computation of a traceback path.
///  The trait is implemented both for the simple and the block-based algorithms.
///
/// Implementors must keep track of two cells in the DP matrix:
/// C = "current" cell
/// L = "left" cell in the diagonal position
///
///   |————————————
///   |   |   |   |
///   |————————————
///   |   | L |   |
///   |————————————
///   |   |   | C |
///   |————————————
///
/// The idea behind this is, that diagonal moves in the DP matrix (= matches /
/// substitutions) should be made as easy/fast as possible, while InDels
/// are usually be less common (especially since there are no affine gap penalties)
/// and upward/leftward moves can thus be more computationally intensive.
///
/// Implementors may store two `State<T, D>` blocks, whose `dist` field reflects
/// the distance score of the current TB matrix cell, and along with that some auxiliary
/// bit vectors indicating cursor position and assisting with adjustments to `dist` based
/// deltas from the PV / MV bit vectors.
pub(super) trait TracebackHandler<'a, T, D>
where
    T: BitVec + 'a,
    D: DistType,
{
    /// Returns the distance score of the current cell C(i, j).
    fn dist(&self) -> D;

    /// Returns the distance score of the left cell, which is in the diagonal
    /// position, C(i-1, j-1).
    fn left_dist(&self) -> D;

    /// Checks if the cell to the left C(i, j-1) has a smaller distance score than
    /// the current cell C(i, j). If so, it calls `move_up()`.
    fn try_move_up(&mut self) -> bool;

    /// Moves up by one position in the DP matrix, adjusting the
    /// distance scores/blocks/bit masks for the current and left columns.
    fn move_up(&mut self);

    /// Checks if the cell to the left C(i-1, j) has a smaller distance score than
    /// the cell in the diagonal C(i-1, j-1). If so, adjusts the score/block
    /// for the left column and returns `true`.
    /// The "current" block needs no adjustment, as it will be discarded in
    /// `finish_move_left()`, which needs to be called afterwards to complete
    /// the move.
    fn try_prepare_left(&mut self) -> bool;

    /// Prepares for a diagonal move (which must be completed with `finish_move_left`).
    /// Essentially, only bit masks and block indices need to be adjusted to
    /// reflect the shift upwards by one position.
    /// The blocks/distances themselves need no modification.
    fn prepare_diagonal(&mut self);

    /// Complete a 'move' to the left by
    /// (1) replacing the current column with the left one (same block and score)
    /// (2) adding a new column to the left and adjusting the distance score
    /// to reflect the diagonal position.
    fn finish_move_left(&mut self);

    /// Returns true if topmost DP matrix row has been reached,
    /// meaning that the complete alignment path has been found.
    fn done(&self) -> bool;

    /// For debugging only
    #[allow(dead_code)]
    fn print_state(&self);
}

#[derive(Clone, Debug)]
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
    /// Creates a new `Traceback` instance.
    /// The `states` vector is potentially reused, so all `State` instances
    /// inside it must be initialized to new values.
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
        states.resize_with(n_states, Default::default);

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

    /// Returns distance for the given end position, or `None` if not stored
    #[inline]
    pub fn dist_at(&self, pos: usize, states: &'a [State<T, D>]) -> Option<D> {
        let pos = pos + 2; // in order to be comparable since self.pos starts at 2, not 0
        if pos <= self.pos {
            return self.handler.dist_at(pos, states).map(|d| d as D);
        }
        None
    }

    /// Returns the length of the current match, optionally adding the
    /// alignment path to `ops`
    #[inline]
    pub fn traceback(
        &self,
        ops: Option<&mut Vec<AlignmentOperation>>,
        states: &'a [State<T, D>],
    ) -> Option<(D, D)> {
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
            return self._traceback_at(pos, ops, states);
        }
        None
    }

    /// Returns a tuple of alignment length and hit distance, optionally adding the alignment path
    /// to `ops`
    #[inline]
    fn _traceback_at(
        &self,
        pos: usize,
        mut ops: Option<&mut Vec<AlignmentOperation>>,
        state_slice: &'a [State<T, D>],
    ) -> Option<(D, D)> {
        use self::AlignmentOperation::*;

        // self.print_tb_matrix(pos, state_slice);

        // Generic object that holds the necessary data and methods
        let mut h = self.handler.init_traceback(self.m, pos, state_slice)?;

        // horizontal column offset from starting point in traceback matrix (bottom right)
        let mut h_offset = D::zero();

        // distance of the match (will be returned)
        let dist = h.dist();

        // Loop for finding the traceback path
        while !h.done() {
            let op;
            // This loop is used to allow skipping `advance_left()` using break (kind of similar
            // to 'goto'). This was done to avoid having to inline advance_left() three times,
            // which would use more space.
            #[allow(clippy::never_loop)]
            loop {
                // h.print_state();

                // If there are several possible solutions, substitutions are preferred over InDels
                // (Subst > Ins > Del)
                if h.left_dist().wrapping_add(&D::one()) == h.dist() {
                    // Diagonal (substitution), move up
                    // Since the left cursor is always in the upper diagonal position,
                    // a simple comparison of distances is enough to determine substitutions.
                    h.prepare_diagonal();
                    op = Subst;
                    // then, move to left (below)
                } else if h.try_move_up() {
                    // Up (insertion to target = deletion in query pattern)
                    op = Ins;
                    break;
                } else if h.try_prepare_left() {
                    // Left (deletion in target text = insertion to query pattern)
                    op = Del;
                    // then, move to left (below)
                } else {
                    // Diagonal (match), move up
                    debug_assert!(h.left_dist() == h.dist());
                    h.prepare_diagonal();
                    op = Match;
                    // then, move to left (below)
                }

                // // This is (probably) equivalent to the Edlib strategy (Ins > Del > Subst)
                // // Tests succeed with this strategy as well, and performance is similar.
                // if h.try_move_up() {
                //     op = Ins;
                //     break;
                // } else if h.left_dist() == h.dist() && h.try_prepare_left() {
                //     op = Del;
                // } else {
                //     h.prepare_diagonal();
                //     op = if h.left_dist().wrapping_add(&D::one()) == h.dist() {
                //         Subst
                //     } else {
                //         Match
                //     };
                // }

                // Moving one column to the left, adjusting h_offset
                h_offset += D::one();
                h.finish_move_left();
                break;
            }

            // dbg!(op);

            if let Some(o) = ops.as_mut() {
                o.push(op);
            }
        }

        Some((h_offset, dist))
    }

    // Useful for debugging
    #[allow(dead_code)]
    fn print_tb_matrix(&self, pos: usize, state_slice: &'a [State<T, D>]) {
        let n_blocks = self.handler.n_blocks();
        let pos = n_blocks * (pos + 1);
        let states_iter = state_slice[..pos]
            .chunks(n_blocks)
            .rev()
            .chain(state_slice.chunks(n_blocks).rev().cycle());
        let m = self.m.to_usize().unwrap();
        let mut out = vec![];
        for col in states_iter {
            let mut col_out = Vec::with_capacity(m);
            let mut empty = true;
            for (i, block) in col.iter().enumerate().rev() {
                if !(block.is_new() || block.is_max()) {
                    empty = false;
                }
                let w = word_size::<T>();
                let end = (i + 1) * w;
                let _m = if end <= m { w } else { m % w };
                let mut _block = *block;
                let mut pos_mask = T::one() << (_m - 1);
                col_out.push(_block.dist);
                for _ in 0.._m {
                    // println!("{}\n{:064b}", _block, pos_mask);
                    _block.adjust_one_up(pos_mask);
                    pos_mask >>= 1;
                    col_out.push(_block.dist);
                }
            }
            out.push(col_out);
            if empty {
                break;
            }
        }

        for j in (0..m).rev() {
            eprint!("{:>4}: ", m - j + 1);
            for col in out.iter().rev() {
                if let Some(d) = col.get(j) {
                    if *d >= (D::max_value() >> 1) {
                        // missing value
                        eprint!("    ");
                    } else {
                        eprint!("{:>4?}", d);
                    }
                } else {
                    eprint!("   -");
                }
            }
            eprintln!();
        }
    }
}
