#![feature(test)]

extern crate test;

use bio::alignment::pairwise::fast::TracebackCell;
use test::Bencher;

#[bench]
fn bench_vec_tb(b: &mut Bencher) {
    b.iter(|| vec![TracebackCell::new(); 1000000]);
}

#[bench]
fn bench_vec_u8(b: &mut Bencher) {
    b.iter(|| vec![0u8; 1000000]);
}
