#![allow(unstable)]

extern crate bio;
extern crate test;


use test::Bencher;
use bio::pattern_matching::shift_and::ShiftAnd;


#[bench]
fn bench_shift_and(b: &mut Bencher) {
    let pattern = b"AAAA";
    let text = b"ACGGCTAGAAAAGGCTAGGAGTAGGATTCTGCATGCACGACTCGAGCACTAGCACT";
    let shiftand = ShiftAnd::new(pattern);
    b.iter(|| {
        shiftand.find_all(text).collect::<Vec<usize>>()
    });
}
