// Copyright 2016 Brett Bowman
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::fmt;

///
/// Possible errors when performing interval operations
///
#[derive(Debug)]
pub enum IntervalError {
    InvalidInterval,
    NoOverlap,
}

impl fmt::Display for IntervalError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self)
    }
}

///
/// An interval, consisting of start and stop position (the latter exclusive).
///
#[derive(PartialEq, Eq, Debug, Copy, Clone)]
pub struct Interval {
    pub start: usize,
    pub stop: usize,
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "<Interval: {0}-{1}>", self.start, self.stop)
    }
}

impl cmp::Ord for Interval {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        // Primary ordering is by Start position
        if self.start < other.start {
            cmp::Ordering::Less
        } else if self.start > other.start {
            cmp::Ordering::Greater
        } else { 
            // Secondary ordering is by Stop position
            if self.stop < other.stop {
                cmp::Ordering::Less
            } else if self.stop > other.stop {
                cmp::Ordering::Greater
            } else {  // If we got here, Intervals are equal
                cmp::Ordering::Equal
            }
        }
    }
}

impl cmp::PartialOrd for Interval {    
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Interval {
    /// Validate that the length is non-negative on initialization
    pub fn new(start: usize, stop: usize) -> Self {
        if stop < start {
            panic!("{0}: start ({1}) greater than stop ({2})", IntervalError::InvalidInterval, start, stop);
        }
        Interval { start: start, stop: stop }
    }

    /// Get the text within the given interval.
    pub fn get<'a>(&self, text: &'a [u8]) -> &'a [u8] {
        &text[self.start..self.stop]
    }

    /// Test whether two intervals cover some of the same range
    pub fn overlaps(&self, other: &Self) -> bool {
        if other.start <= self.start && self.start <= other.stop {
            true  // Overlaps with left edge
        } else if self.start <= other.start && other.start <= self.stop {
            true  // Overlaps with right edge
        } else {
            false
        }
    }

    /// Test whether a value is within this interval
    pub fn contains(&self, value: usize) -> bool {
        if self.start <= value && self.stop >= value {
            true
        } else {
            false
        }
    }
    
    /// Test whether one interval is completely within another
    pub fn covers(&self, other: &Self) -> Result<bool, IntervalError> {
        let intersection = self.intersect(other);
        match intersection {
            Ok(v)  => Ok(&v == other),
            Err(e) => Err(e),
        }
    }

    /// Return the range covered by both intervals
    pub fn intersect(&self, other: &Self) -> Result<Interval, IntervalError> {
        if !self.overlaps(other) {
            Err(IntervalError::NoOverlap)
        } else {
            Ok(Interval{ start: cmp::max(self.start, other.start), 
                         stop: cmp::min(self.stop, other.stop) })
        }
    }

    /// Return the length of this interval
    pub fn length(&self) -> usize {
        self.stop - self.start
    }

    /// Return the range covered by either interval
    pub fn union(&self, other: &Self) -> Result<Interval, IntervalError> {
        if !self.overlaps(other) {
            Err(IntervalError::NoOverlap)
        } else {
            Ok(Interval{ start: cmp::min(self.start, other.start), 
                         stop: cmp::max(self.stop, other.stop) })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

}
