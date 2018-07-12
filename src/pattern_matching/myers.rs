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
//! use bio::pattern_matching::myers::Myers;
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
use std::u64;
use std::borrow::Borrow;
use std::mem::size_of;
use std::ops::*;

use num_traits::*;

use utils::{TextSlice, IntoTextIterator, TextIterator};



/// Integer types serving as bit vectors must implement this trait.
/// This allows the code to be less verbose and should make it simpler
/// to implement for arbitrary types.
/// Functions not part of a trait in the standard library are to be implemented
/// separately. There exists an implementation based on the `num_traits` crate.
/// Requiring num_traits traits to be implemented for custom types would be complicated,
/// for instance PrimInt has loads of functions, while only 'count_ones' is used.
pub trait Num: Add<Output=Self> + Sub<Output=Self> +
    BitOr<Output=Self> + BitOrAssign + BitAnd<Output=Self> + BitXor<Output=Self> + Not<Output=Self> +
    Shl<usize, Output=Self> + ShlAssign<usize> + ShrAssign<usize> +
    Copy + PartialOrd
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
// TODO: add feature switch
/// Myers instance based on 128-bit integers (pattern length up to 128)
//#[cfg(feature="nightly")]
pub type Myers128 = Myers<u128>;



/// Myers algorithm.
pub struct Myers<T: Num> {
    peq: [T; 256],
    bound: T,
    m: usize,
}



impl<T: Num> Myers<T> {
    /// Create a new instance of Myers algorithm for a given pattern.
    pub fn new<'a, P: IntoTextIterator<'a>>(pattern: P) -> Self {
        Self::with_ambig(pattern.into_iter().map(Some))
    }

    // TODO: eventually move to module docs and use shorter example
    /// Like `Myers::new()`, but additionally allows for specifying
    /// multiple matching characters in a pattern.
    ///
    /// Example:
    ///
    /// ```
    /// use pattern_matching::myers_st::Myers64;
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
        Self::with_ambig_results::<_, _, _, ()>(
            pattern.into_iter().map(Ok)
        ).unwrap()
    }

    // TODO: eventually move to module docs and use shorter example
    /// Like `Myers::with_ambig()`, but allows returning a `Result` during iteration
    /// of the pattern. The following example enables matching all ambiguous characters according
    /// to the IUPAC nomenclature.
    ///
    /// Example:
    ///
    /// ```
    /// use pattern_matching::myers_st::Myers64;
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
        Matches {
            myers: self,
            state: State::new(self.m),
            text: text.into_iter().enumerate(),
            max_dist: max_dist,
        }
    }
}


/// The current algorithm state.
struct State<T> {
    pv: T,
    mv: T,
    dist: usize,
}


impl<T> State<T> where T: Num {
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


#[cfg(test)]
mod tests {
    use super::*;
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
        let text = b"TGAGCNT";
        let pattern = b"TGAGCGT";

        let myers = Myers64::new(pattern);
        assert_eq!(myers.distance(text), 1);

        let myers_wildcard = Myers64::with_wildcard(pattern, b'N');
        assert_eq!(myers_wildcard.distance(text), 0);
    }

    #[test]
    fn test_long() {
        let text = b"ACCGTGGATGAGCGCCATAGACCGTGGATGAGCGCCATAGACCGTGGATGAGCGCCATAGACCGTGGATGAGCGCCATAGACCGTGGATGAGCGCCATAG";
        let pattern = b"TGAGCGTTGAGCGTTGAGCGTTGAGCGTTGAGCGTTGAGCGT";
        let myers = Myers64::new(&pattern[..]);
        let occ = myers.find_all_end(&text[..], 1).collect_vec();
        println!("{:?}", occ);
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
