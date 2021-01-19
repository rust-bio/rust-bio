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
        .fold((0usize, 0usize), |(l, count), n| match *n.borrow() {
            b'c' | b'g' | b'G' | b'C' => (l + 1, count + 1),
            _ => (l + 1, count),
        });
    count as f32 / l as f32
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
/// use approx::assert_relative_eq;
/// assert_relative_eq!(gc_content(seq), 2. / 8., epsilon = f32::EPSILON);
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
/// use approx::assert_relative_eq;
/// use bio::seq_analysis::gc::gc3_content;
/// const seq: &'static [u8] = b"GATATACA";
/// //                           ^  ^  ^
/// assert_relative_eq!(gc3_content(seq), 2. / 3., epsilon = f32::EPSILON);
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
        assert_relative_eq!(gc_content(gc0), 0.0, epsilon = f32::EPSILON);
        let gc50 = b"ATGC";
        assert_relative_eq!(gc_content(gc50), 0.5, epsilon = f32::EPSILON);
        let gc100 = b"GCGC";
        assert_relative_eq!(gc_content(gc100), 1.0, epsilon = f32::EPSILON);
    }

    #[test]
    fn test_gc_content_large() {
        const LENGTH: usize = 10_000_000;
        let mut s = vec![b'G'; LENGTH];
        s.extend_from_slice(&[b'T'; LENGTH]);
        let gc_content = gc_content(s);
        assert_relative_eq!(gc_content, 0.5, epsilon = f32::EPSILON);
    }
}
