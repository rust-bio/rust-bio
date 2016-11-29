// Copyright 2016 Brian Rogoff
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use alphabets::dna;
use utils::TextSlice;
use std::f64;

// We represent DNA kmers of size less than or equal to 32 by a 64 bit wide
// unsigned integer wrapped in a 'newtype' struct marker.

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Hash)]
pub struct DNAKMer(u64);

pub fn dna_to_number(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => panic!("Invalid byte")
    }
}

pub fn number_to_dna(n: u64) -> u8 {
    match n {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => panic!("Invalid number")
    }
}

pub fn pattern_to_kmer(pattern: &[u8]) -> DNAKMer {
    let mut mul = 1;
    let mut result = 0;
    let l = pattern.len();
    for i in 0..l {
        result += mul * dna_to_number(pattern[l - 1 - i]);
        mul *= 4;
    }

    return DNAKMer(result);
}

pub fn pattern_to_canonical(pattern: &[u8]) -> DNAKMer {
    let mut mul = 1;
    let mut result0 = 0;
    let mut result1 = 0;
    let l = pattern.len();
    for i in 0..l {
        result0 += mul * dna_to_number(pattern[l - 1 - i]);
        result1 += mul * dna_to_number(dna::complement(pattern[i]));
        mul *= 4;
    }

    if result0 < result1 {
        return DNAKMer(result0);
    } else {
        return DNAKMer(result1);
    }
}

// This function is useful in printing out a dnakmer, e.g.
// let mut kmer_buf: [u8; 32] = [b'A'; 32];
// dnakmer_into_slice(k, kmer, &mut kmer_buf[0..k]);
// println!("kmer = {}, kmer_string = {}",
//          kmer,
//          String::from_utf8(Vec::from(&kmer_buf[0..k])).unwrap()
// );

pub fn dnakmer_into_slice(k: usize, kmer: DNAKMer, slice: &mut [u8]) {
    let four_power_k = 4u64.pow(k as u32);
    assert!(kmer < DNAKMer(four_power_k));
    let DNAKMer(mut remainder) = kmer;
    for i in 0..k {
        let last_digit = remainder % 4;
        slice[k - 1 - i] = number_to_dna(last_digit);
        remainder /= 4;
    }
}

pub struct DNAKMerSeq<'a> {
    k: usize,
    pos: usize,
    data: TextSlice<'a>
}

// Implement `Iterator` for `KMerSeq`. We use the canonical representation of
// the kmer, which is the one with the smaller numerical value. Factoring out
// the function which returns the kmer would allow this to be parameterizable,
// but I haven't measured the cost yet.

impl<'a> Iterator for DNAKMerSeq<'a> {
    type Item = DNAKMer;

    fn next(&mut self) -> Option<DNAKMer> {
        let pos = self.pos;
        let k = self.k;
        let result = 
            if pos + self.k < self.data.len() {
                Some(pattern_to_canonical(&self.data[pos..pos+k-1]))
            } else {
                None
            };
        self.pos += k;
        result
    }
}

// Returns a DNAKMerSeq sequence generator
pub fn dnakmer_seq<'a>(k: usize, s: &'a [u8]) -> DNAKMerSeq<'a> {
    DNAKMerSeq { k: k, data: s, pos: 0 }
}
