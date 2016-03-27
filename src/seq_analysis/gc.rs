// Copyright 2014-2016 Johannes Köster, Martin Larralde.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

// Sequence analysis algorithms such as orf finder or GC% counter

// Maybe not the most clever way to do it, but at least runs in o(n)

pub fn gc <'a>(sequence: &'a [u8]) -> f32 {
    let mut n: f32 = 0.0;
    let l: usize = sequence.len();
    let mut i: usize = 0;
    while i < l {
        n += match sequence[i] {
            67 | 71 => 1f32, // G or C
            _       => 0f32,
        } ;
        i += 1
    }
    n / l as f32
}

pub fn gc3 <'a>(sequence: &'a [u8]) -> f32 {
    let mut n: f32 = 0.0;
    let l: usize = sequence.len();
    let mut i: usize = 2;
    while i < l {
        n += match sequence[i] {
            67 | 71 => 1f32,
            _       => 0f32,
        };
        i += 3;
    }
    n / l as f32
}
