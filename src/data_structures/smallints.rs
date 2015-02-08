//! A data structure for a sequence of small integers with a few big integers.
//! Small ints are stored in type S (e.g. a byte), big ints are stored separately (in type B) in a BTree.
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

use std::num::{Int, NumCast, FromPrimitive, cast};
use std::iter::{Enumerate, repeat};
use std::mem::size_of;
use std::collections::BTreeMap;
use std::slice;


pub struct SmallInts<S: Int + NumCast + FromPrimitive, B: Int + NumCast> {
    smallints: Vec<S>,
    bigints: BTreeMap<usize, B>
}


impl<S: Int + NumCast + FromPrimitive, B: Int + NumCast> SmallInts<S, B> {
    pub fn new() -> Self {
        assert!(size_of::<S>() < size_of::<B>(), "S has to be smaller than B");
        SmallInts { smallints: Vec::new(), bigints: BTreeMap::new() }
    }

    pub fn with_capacity(n: usize) -> Self {
        assert!(size_of::<S>() < size_of::<B>(), "S has to be smaller than B");
        SmallInts { smallints: Vec::with_capacity(n), bigints: BTreeMap::new() }
    }

    pub fn from_elem(v: S, n: usize) -> Self {
        assert!(size_of::<S>() < size_of::<B>(), "S has to be smaller than B");
        if v > FromPrimitive::from_int(0).unwrap() {
            assert!(v < Int::max_value(), "v has to be smaller than maximum value");
        }
        
        SmallInts { smallints: repeat(v).take(n).collect(), bigints: BTreeMap::new() }
    }

    pub fn get(&self, i: usize) -> Option<B> {
        if i < self.smallints.len() {
            self.real_value(i, self.smallints[i])
        }
        else {
            None
        }
    }

    pub fn push(&mut self, v: B) {
        let maxv: S = Int::max_value();
        match cast(v) {
            Some(v) if v < maxv => self.smallints.push(v),
            _                   => {
                let i = self.smallints.len();
                self.smallints.push(maxv);
                self.bigints.insert(i, v);
            }
        }
    }

    pub fn set(&mut self, i: usize, v: B) {
        let maxv: S = Int::max_value();
        match cast(v) {
            Some(v) if v < maxv => self.smallints[i] = v,
            _                   => {
                self.smallints[i] = maxv;
                self.bigints.insert(i, v);
            }
        }
    }

    pub fn iter<'a>(&'a self) -> Iter<S, B> {
        Iter { smallints: &self, items: self.smallints.iter().enumerate() }
    }

    pub fn decompress(&self) -> Vec<B> {
        self.iter().collect()
    }

    pub fn len(&self) -> usize {
        self.smallints.len()
    }

    fn real_value(&self, i: usize, v: S) -> Option<B> {
        if v < Int::max_value() {
            cast(v)
        }
        else {
            self.bigints.get(&i).cloned()
        }
    }
}


pub struct Iter<'a, S: 'a + Int + NumCast + FromPrimitive, B: 'a + Int + NumCast> {
    smallints: &'a SmallInts<S, B>,
    items: Enumerate<slice::Iter<'a, S>>
}


impl<'a, S: 'a + Int + NumCast + FromPrimitive, B: 'a + Int + NumCast> Iterator for Iter<'a, S, B> {
    type Item = B;

    fn next(&mut self) -> Option<B> {
        match self.items.next() {
            Some((i, &v)) => self.smallints.real_value(i, v),
            None => None
        }
    }
}
