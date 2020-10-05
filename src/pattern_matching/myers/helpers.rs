use std::fmt::{Binary, Debug};
use std::ops::*;

use num_traits::{AsPrimitive, FromPrimitive, PrimInt, ToPrimitive, WrappingAdd};

/// Trait for types that should be used to store the distance score when using the simple
/// Myers algorithm (not the block-based one, which always uses `usize`).
///
/// For all currently implemented BitVec types, the maximum possible distance
/// can be stored in `u8`. Custom implementations using bigger integers can
/// adjust `DistType` to hold bigger numbers. Note that due to how the traceback
/// algorithm currently works, `DistType` should be able to represent numbers larger
/// than the bit-width of the `BitVec` type. For instance, a hypothetical `BitVec` type
/// of `u256` should use `u16` as distance, since `u8` cannot store numbers larger
/// than 255.
pub trait DistType: Copy
        + Debug
        + Default
        + AddAssign
        + SubAssign
        + PrimInt // includes Bounded, Num, Zero, One
        + FromPrimitive
        + ToPrimitive
        + AsPrimitive<usize> + AsPrimitive<i64>
        + WrappingAdd
        + Sub<Output=Self> {}

impl DistType for u8 {}
impl DistType for u16 {}
impl DistType for u32 {}
impl DistType for u64 {}
impl DistType for usize {}

/// This trait must be implemented for integer types serving as bit vectors.
/// Only unsigned integers will work correctly.
pub trait BitVec: Copy
    + Debug
    + Binary
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
    + ToPrimitive
    + FromPrimitive
    + AsPrimitive<usize> + AsPrimitive<i64>
{
    /// Type that should be used to store the distance score when using the simple
    /// Myers algorithm (not the block-based one, which always uses `usize`).
    type DistType: DistType;
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

use crate::alignment::{Alignment, AlignmentMode};

/// Updates an `Alignment` instance with new data (except of path).
/// Assumes *0-based* range end coordinates, they will be converted to 1-based ones
#[inline(always)]
pub(crate) fn update_aln(
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

#[inline]
pub(crate) fn word_size<T>() -> usize {
    std::mem::size_of::<T>() * 8
}

#[inline]
pub(crate) fn ceil_div(x: usize, y: usize) -> usize {
    if x % y != 0 {
        x / y + 1
    } else {
        x / y
    }
}
