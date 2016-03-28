// Copyright 2014-2016 Johannes Köster, Martin Larralde.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

// Sequence analysis algorithms such as orf finder or GC% counter

// Maybe not the most clever way to do it, but at least runs in o(n)

/// Base gc content counter
fn gcn_content<'a>(sequence: &'a[u8], step: usize) -> f32 {
    let mut count = 0.0;
    let len = sequence.len();
    let mut i: usize = 0;
    while i < len {
        count += match sequence[i] {
            b'c'|b'g'|b'G'|b'C' => 1f32, // G or C
                              _ => 0f32,
        };
        i += step
    }
    count / len as f32
}

/// gc content counter for every nucleotide
pub fn gc_content<'a>(sequence: &'a[u8]) -> f32 {
    gcn_content(sequence, 1usize)
}

/// gc content counter for the nucleotide in 3rd position
pub fn gc3_content<'a>(sequence: &'a[u8]) -> f32 {
    gcn_content(sequence, 3usize)
}
