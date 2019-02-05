// Copyright 2014-2016 Johannes Köster, Martin Larralde.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! GC counter over an `IntoTextIterator` object.

//! Complexity: o(n)

use std::borrow::Borrow;

/// Base gc content counter
fn gcn_content<C: Borrow<u8>, T: IntoIterator<Item = C>>(sequence: T, step: usize) -> f32 {
    let mut l = 0f32;
    let mut count = 0.0;
    for (i, n) in sequence.into_iter().enumerate() {
        if i % step == 0 {
            count += match *n.borrow() {
                b'c' | b'g' | b'G' | b'C' => 1f32, // G or C
                _ => 0f32,
            };
        }
        l = i as f32;
    }
    count / (l + 1f32)
}

/// gc content counter for every nucleotide
pub fn gc_content<C: Borrow<u8>, T: IntoIterator<Item = C>>(sequence: T) -> f32 {
    gcn_content(sequence, 1usize)
}

/// gc content counter for the nucleotide in 3rd position
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
