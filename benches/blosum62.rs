#![feature(test)]

extern crate test;

use bio::scores::{blosum62, blosum62_array};
use test::Bencher;

#[bench]
fn bench_blosum62_lookup(b: &mut Bencher) {
    b.iter(|| {
        for i in 65..91 {
            for j in 65..91 {
                blosum62::blosum62(i, j);
            }
        }
    });
}

#[bench]
fn bench_blosum62_lookup_array(b: &mut Bencher) {
    b.iter(|| {
        for i in 65..91 {
            for j in 65..91 {
                blosum62_array::blosum62(i, j);
            }
        }
    });
}
