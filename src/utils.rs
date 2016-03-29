// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Common utilities.


use std::cmp;
use num::Float;

/// Type alias for an owned text, i.e. ``Vec<u8>``.
pub type Text = Vec<u8>;
/// Type alias for a text slice, i.e. ``&[u8]``.
pub type TextSlice<'a> = &'a [u8];

/// Type alias for an iterator over a sequence, i.e. ``Iterator<Item=&u8>``.
pub trait TextIterator<'a>: Iterator<Item=&'a u8> {}
impl<'a, I: Iterator<Item=&'a u8>> TextIterator<'a> for I {}

/// Type alias for a type that can be coerced into a ``TextIterator``.
/// This includes ``&Vec<u8>``, ``&[u8]``, ``Iterator<Item=&u8>``.
pub trait IntoTextIterator<'a>: IntoIterator<Item=&'a u8> {}
impl<'a, T: IntoIterator<Item=&'a u8>> IntoTextIterator<'a> for T {}


/// Remove a trailing newline from the given string in place.
pub fn trim_newline(s: &mut String) {
    if s.ends_with('\n') {
        s.pop();
    }
}


/// In place implementation of scan over a slice.
pub fn scan<T: Copy, F: Fn(T, T) -> T>(a: &mut [T], op: F) {
    let mut s = a[0];
    for v in a.iter_mut().skip(1) {
        s = op(s, *v);
        *v = s;
    }
}


// Inplace implementation of prescan over a slice.
pub fn prescan<T: Copy, F: Fn(T, T) -> T>(a: &mut [T], neutral: T, op: F) {
    let mut s = neutral;
    for v in a.iter_mut() {
        let t = *v;
        *v = s;
        s = op(s, t);
    }
}


#[derive(PartialOrd, PartialEq, Debug, Copy, Clone)]
pub struct NonNaNFloat<F: Float>(F);


impl<F: Float> NonNaNFloat<F> {
    pub fn new(v: F) -> Option<Self> {
        if v.is_nan() {
            Some(NonNaNFloat(v))
        } else {
            None
        }
    }

    pub fn unwrap(&self) -> F {
        let &NonNaNFloat(v) = self;
        v
    }
}


impl<F: Float> Eq for NonNaNFloat<F> {}


impl<F: Float> Ord for NonNaNFloat<F> {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    #[test]
    fn test_scan() {
        let mut a = vec![1, 0, 0, 1];
        scan(&mut a[..], |a, b| a + b);
        assert_eq!(a, vec![1, 1, 1, 2]);
    }

    #[test]
    fn test_prescan() {
        let mut a = vec![1, 0, 0, 1];
        prescan(&mut a[..], 0, |a, b| a + b);
        assert_eq!(a, vec![0, 1, 1, 1]);
    }

    #[test]
    fn test_nonnanfloat() {
        let mut v = [NonNaNFloat(5.1), NonNaNFloat(1.3)];
        v.sort();
        assert_eq!(v, [NonNaNFloat(1.3), NonNaNFloat(5.1)]);
    }

    /// This function demonstrates the use of the IntoSequenceIterator alias, which takes both
    /// slices and iterators.
    fn print_sequence<'a, I: IntoTextIterator<'a>>(sequence: I) {
        for c in sequence {
            println!("{}", c);
        }
    }

    #[test]
    fn test_print_sequence() {
        let s = b"ACGT";
        // use iterator
        print_sequence(s.iter().step(1));
        // use slice
        print_sequence(&s[..]);
        // use vec
        print_sequence(&vec![b'A', b'C']);
        // keep ownership
        println!("{:?}", s);
    }
}
