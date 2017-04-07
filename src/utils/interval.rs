use std::ops::{Range, Deref};

/// An `Interval` wraps the `std::ops::Range` from the stdlib and is defined by a start and end field
/// where end should be >= start.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Interval<N: Ord + Clone>(Range<N>);

impl<N: Ord + Clone> Interval<N> {
    /// Construct a new `Interval` from the given Range.
    /// Will return `Err` if end < start.
    pub fn new(r: Range<N>) -> Result<Interval<N>, IntervalError> {
        if r.end >= r.start {
            Ok(Interval(r))
        } else {
            Err(IntervalError::InvalidRange)
        }
    }
}

/// Convert a `Range` into an `Interval`. This conversion will panic if the `Range` has end < start
impl<N: Ord + Clone> From<Range<N>> for Interval<N> {
    fn from(r: Range<N>) -> Self {
        match Interval::new(r) {
            Ok(interval) => interval,
            Err(IntervalError::InvalidRange) => panic!("Cannot convert negative width range to interval"),
        }
    }
}

/// Convert a reference to a `Range` to an interval by cloning. This conversion will panic if the
/// `Range` has end < start
impl<'a, N: Ord + Clone> From<&'a Range<N>> for Interval<N> {
    fn from(r: &Range<N>) -> Self {
        match Interval::new(r.clone()) {
            Ok(interval) => interval,
            Err(IntervalError::InvalidRange) => panic!("Cannot convert negative width range to interval"),
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

quick_error! {
    #[derive(Debug)]
    pub enum IntervalError {
        InvalidRange {
            description("An Interval must have a Range with a positive width")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Interval;

    #[test]
    #[should_panic]
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
