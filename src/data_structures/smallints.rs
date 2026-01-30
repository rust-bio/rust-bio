// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A data structure for a sequence of small integers with a few big integers.
//! Small ints are stored in type S (e.g. a byte), big ints are stored separately (in type `B`) in a BTree.
//! The implementation provides vector-like operations on the data structure (e.g. retrieve a position,
//! add an integer, etc.). Getting and setting (by position) time complexity is `O(1)` for small ints,
//! and `O(b)` for big ints, where `b` is the number of big ints stored.
//!
//! # Space usage
//! SmallInts pay the cost of slower retrieval of big integers with smaller overall memory usage;
//! `O(size_of(S) * (s+b) + size_of(B) * b)` where `S` and `B` are the small and large int types, and
//! `s` and `b` are the number of those stored respectively.
//!
//! # Example
//!
//! ```
//! use bio::data_structures::smallints::SmallInts;
//! let mut smallints: SmallInts<u8, usize> = SmallInts::new();
//! smallints.push(3);
//! smallints.push(4);
//! smallints.push(255);
//! smallints.push(305093);
//! assert_eq!(smallints.get(0).unwrap(), 3);
//! smallints.set(0, 50000);
//! let values: Vec<usize> = smallints.iter().collect();
//! assert_eq!(values, [50000, 4, 255, 305093]);
//! ```

use std::collections::BTreeMap;
use std::iter::{repeat_n, Enumerate};
use std::mem::size_of;
use std::slice;

use num_integer::Integer;
use num_traits::{cast, Bounded, Num, NumCast};

/// Data structure for storing a sequence of small integers with few big ones space efficiently
/// while supporting classical vector operations.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct SmallInts<F: Integer + Bounded + NumCast + Copy, B: Integer + NumCast + Copy> {
    smallints: Vec<F>,
    bigints: BTreeMap<usize, B>,
}

impl<S: Integer + Bounded + NumCast + Copy, B: Integer + NumCast + Copy> Default
    for SmallInts<S, B>
{
    fn default() -> Self {
        assert!(
            size_of::<S>() < size_of::<B>(),
            "S has to be smaller than B"
        );
        SmallInts {
            smallints: Vec::new(),
            bigints: BTreeMap::new(),
        }
    }
}

impl<S: Integer + Bounded + NumCast + Copy, B: Integer + NumCast + Copy> SmallInts<S, B> {
    /// Create a new instance.
    pub fn new() -> Self {
        Default::default()
    }

    /// Create a new instance with a given capacity.
    pub fn with_capacity(n: usize) -> Self {
        assert!(
            size_of::<S>() < size_of::<B>(),
            "S has to be smaller than B"
        );
        SmallInts {
            smallints: Vec::with_capacity(n),
            bigints: BTreeMap::new(),
        }
    }

    /// Create a new instance containing `n` times the integer `v` (and `v` is expected to be small).
    pub fn from_elem(v: S, n: usize) -> Self {
        assert!(
            size_of::<S>() < size_of::<B>(),
            "S has to be smaller than B"
        );
        if v > cast(0).unwrap() {
            assert!(v < S::max_value(), "v has to be smaller than maximum value");
        }

        SmallInts {
            smallints: repeat_n(v, n).collect(),
            bigints: BTreeMap::new(),
        }
    }

    /// Return the integer at position `i`. Time complexity `O(1)` if `i` points to a small int,
    /// `O(log(b))` for a big int, where `b` denotes the number of big ints stored.
    pub fn get(&self, i: usize) -> Option<B> {
        if i < self.smallints.len() {
            self.real_value(i, self.smallints[i])
        } else {
            None
        }
    }

    /// Append `v` to the sequence. This will determine whether `v` is big or small and store it
    /// accordingly. Time complexity `O(1)` for small ints and `O(log(b))` for big ints,
    /// where `b` denotes the number of big ints stored.
    pub fn push(&mut self, v: B) {
        let maxv: S = S::max_value();
        match cast(v) {
            Some(v) if v < maxv => self.smallints.push(v),
            _ => {
                let i = self.smallints.len();
                self.smallints.push(maxv);
                self.bigints.insert(i, v);
            }
        }
    }

    /// Set value of position `i` to `v`. This will determine whether `v` is big or small and store it accordingly.
    /// Time complexity `O(1)` for small ints and `O(log(b))` for big ints,
    /// where `b` denotes the number of big ints stored.
    pub fn set(&mut self, i: usize, v: B) {
        let maxv: S = S::max_value();
        match cast(v) {
            Some(v) if v < maxv => self.smallints[i] = v,
            _ => {
                self.smallints[i] = maxv;
                self.bigints.insert(i, v);
            }
        }
    }

    /// Iterate over sequence. Values will be returned in the big integer type (`B`).
    pub fn iter(&self) -> Iter<'_, S, B> {
        Iter {
            smallints: self,
            items: self.smallints.iter().enumerate(),
        }
    }

    /// Decompress into a normal vector of big integers (type `B`).
    pub fn decompress(&self) -> Vec<B> {
        self.iter().collect()
    }

    /// Length of the sequence.
    pub fn len(&self) -> usize {
        self.smallints.len()
    }

    /// is the sequence empty?
    pub fn is_empty(&self) -> bool {
        self.smallints.is_empty()
    }

    fn real_value(&self, i: usize, v: S) -> Option<B> {
        if v < S::max_value() {
            cast(v)
        } else {
            self.bigints.get(&i).cloned()
        }
    }
}

/// Iterator over the elements of a `SmallInts` sequence.
#[derive(Clone, Debug)]
pub struct Iter<'a, S, B>
where
    S: Integer + Bounded + NumCast + Copy,
    B: Integer + NumCast + Copy,
    <S as Num>::FromStrRadixErr: 'a,
    <B as Num>::FromStrRadixErr: 'a,
{
    smallints: &'a SmallInts<S, B>,
    items: Enumerate<slice::Iter<'a, S>>,
}

impl<'a, S, B> Iterator for Iter<'a, S, B>
where
    S: 'a + Integer + Bounded + NumCast + Copy,
    B: 'a + Integer + NumCast + Copy,
    <S as Num>::FromStrRadixErr: 'a,
    <B as Num>::FromStrRadixErr: 'a,
{
    type Item = B;

    fn next(&mut self) -> Option<B> {
        match self.items.next() {
            Some((i, &v)) => self.smallints.real_value(i, v),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::data_structures::smallints::SmallInts;
    #[test]
    fn test_serde() {
        use serde::{Deserialize, Serialize};
        fn impls_serde_traits<'a, S: Serialize + Deserialize<'a>>() {}

        impls_serde_traits::<SmallInts<i8, isize>>();
    }
}
