//! This module defines a newtype `Interval` for `std::ops::Range`, which will panic if `end` < `start`.
//!
//! # Examples
//! Create a new `Interval` given a `Range`.
//! ```
//! use bio::utils::Interval;
//! assert_eq!(Interval::new(3..6).unwrap(), (3..6).into());
//! ```
//!
//! Building an `Interval` from a `Range` with start > end should panic.
//! ```should_panic
//! use bio::utils::Interval;
//! Interval::from(7..1);
//! ```
//!
//! If you want to handle invalid ranges properly, use the `new` constructor
//! ```
//! use bio::utils::Interval;
//! match Interval::new(7..1) {
//!     Ok(interval) => println!("{:?}", interval),
//!     Err(error) => eprintln!("interval start > end"),
//! }
//! ```

pub mod errors;

use std::ops::{Deref, Range};

pub use self::errors::{Error, Result};

/// An `Interval` wraps the `std::ops::Range` from the stdlib and is defined by a start and end field
/// where end should be >= start.
#[derive(Default, Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
pub struct Interval<N: Ord + Clone>(Range<N>);

impl<N: Ord + Clone> Interval<N> {
    /// Construct a new `Interval` from the given Range.
    /// Will return `Err` if end < start.
    pub fn new(r: Range<N>) -> Result<Interval<N>> {
        if r.end >= r.start {
            Ok(Interval(r))
        } else {
            Err(Error::InvalidRange)
        }
    }
}

/// Convert a `Range` into an `Interval`. This conversion will panic if the `Range` has end < start
impl<N: Ord + Clone> From<Range<N>> for Interval<N> {
    fn from(r: Range<N>) -> Self {
        match Interval::new(r) {
            Ok(interval) => interval,
            Err(Error::InvalidRange) => panic!("Cannot convert negative width range to interval"),
        }
    }
}

/// Convert a reference to a `Range` to an interval by cloning. This conversion will panic if the
/// `Range` has end < start
impl<N: Ord + Clone> From<&Range<N>> for Interval<N> {
    fn from(r: &Range<N>) -> Self {
        match Interval::new(r.clone()) {
            Ok(interval) => interval,
            Err(Error::InvalidRange) => panic!("Cannot convert negative width range to interval"),
        }
    }
}

/// Use the `Deref` operator to get a reference to `Range` wrapped by the `Interval` newtype.
impl<N: Ord + Clone> Deref for Interval<N> {
    type Target = Range<N>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::Interval;

    #[test]
    #[should_panic]
    #[allow(clippy::reversed_empty_ranges)]
    fn negative_width_range() {
        let _ = Interval::from(10..5);
    }

    #[test]
    fn range_interval_conversions() {
        assert_eq!(Interval::new(1..10).unwrap(), (1..10).into());
        assert_eq!(Interval::from(1..10), Interval::new(1..10).unwrap());
        //deref access
        let r = Interval::new(1..10).unwrap();
        assert_eq!(*r, (1..10));
        assert_eq!(r.start, 1);
        assert_eq!(r.end, 10);
    }
}
