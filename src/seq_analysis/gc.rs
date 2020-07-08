// Copyright 2014-2016 Johannes Köster, Martin Larralde.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! GC counter over an `IntoTextIterator` object.

//! Complexity: O(n), where n is the length of the sequence.

use std::borrow::Borrow;

/// Base gc content counter
fn gcn_content<C: Borrow<u8>, T: IntoIterator<Item = C>>(sequence: T, step: usize) -> f32 {
    let (l, count) = sequence
        .into_iter()
        .step_by(step)
        .fold((0.0, 0.0), |(l, count), n| match *n.borrow() {
            b'c' | b'g' | b'G' | b'C' => (l + 1.0, count + 1.0),
            _ => (l + 1.0, count),
        });
    count / l
}

/// Returns the ratio of bases which are guanine or cytososine
///
/// # Arguments
///
/// * `sequence` - A sequence of bases
///
/// # Example
///
/// ```
/// use bio::seq_analysis::gc::gc_content;
///
/// const seq: &'static [u8] = b"GATATACA";
/// assert_eq!(gc_content(seq), 2. / 8.);
/// ```
pub fn gc_content<C: Borrow<u8>, T: IntoIterator<Item = C>>(sequence: T) -> f32 {
    gcn_content(sequence, 1usize)
}

/// Returns the ratio of bases in the 3rd position which are guanine
/// or cytososine.
///
/// # Arguments
///
/// * `sequence` - A sequence of bases
///
/// # Example
///
/// ```
/// use bio::seq_analysis::gc::gc3_content;
/// const seq: &'static [u8] = b"GATATACA";
/// //                           ^  ^  ^
/// assert_eq!(gc3_content(seq), 2. / 3.);
/// ```
pub fn gc3_content<C: Borrow<u8>, T: IntoIterator<Item = C>>(sequence: T) -> f32 {
    gcn_content(sequence, 3usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content() {
        let gc0 = b"ATAT";
        assert_eq!(gc_content(gc0), 0.0);
        let gc50 = b"ATGC";
        assert_eq!(gc_content(gc50), 0.5);
        let gc100 = b"GCGC";
        assert_eq!(gc_content(gc100), 1.0);
    }
}
