// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Common utilities.

mod text;
pub use self::text::{Text, TextSlice, TextIterator, IntoTextIterator, trim_newline};

mod nonnanfloat;
pub use self::nonnanfloat::NonNaNFloat;

mod interval;
pub use self::interval::{Interval, IntervalError};


/// Strand information.
#[derive(Debug, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}


impl PartialEq for Strand {
    /// Returns true if both are `Forward` or both are `Reverse`, otherwise returns false.
    fn eq(&self, other: &Strand) -> bool {
        match (self, other) {
            (&Strand::Forward, &Strand::Forward) => true,
            (&Strand::Reverse, &Strand::Reverse) => true,
            _ => false
        }
    }
}

impl Strand {

    /// Returns a `Strand` enum representing the given char.
    ///
    /// The mapping is as follows:
    ///     * '+', 'f', or 'F' becomes `Strand::Forward`
    ///     * '-', 'r', or 'R' becomes `Strand::Reverse`
    ///     * '.', '?' becomes `Strand::Unknown`
    ///     * Any other inputs will return an `Err(StrandError::InvalidChar)`
    pub fn from_char(strand_char: &char) -> Result<Strand, StrandError> {
        match *strand_char {
            '+' | 'f' | 'F' => Ok(Strand::Forward),
            '-' | 'r' | 'R' => Ok(Strand::Reverse),
            '.' | '?'  => Ok(Strand::Unknown),
            invalid => Err(StrandError::InvalidChar(invalid)),
        }
    }

    pub fn is_unknown(&self) -> bool {
        match self {
            &Strand::Unknown => true,
            _ => false,
        }
    }
}

quick_error! {
    #[derive(Debug)]
    pub enum StrandError {
        InvalidChar(invalid_char: char) {
            description("invalid character for strand conversion")
            display("character {:?} can not be converted to a Strand", invalid_char)
        }
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

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_strand() {
        assert_eq!(Strand::from_char(&'+').unwrap(), Strand::Forward);
        assert_eq!(Strand::from_char(&'-').unwrap(), Strand::Reverse);
        assert!(Strand::from_char(&'.').unwrap().is_unknown());
        assert!(Strand::from_char(&'o').is_err());
    }
}
