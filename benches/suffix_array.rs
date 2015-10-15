#![feature(test)]

extern crate test;
extern crate bit_vec;
extern crate bio;

use test::Bencher;

use bio::data_structures::suffix_array::*;

#[bench]
fn bench_suffix_array(b: &mut Bencher) {
    b.iter(|| suffix_array(b"GCCTTAACATTATTACGCCTA$"));
}