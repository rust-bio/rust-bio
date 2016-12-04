// Copyright 2016 Brian Rogoff
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use alphabets::dna;
use utils::TextSlice;

const KMER_MASKS: [u64; 33] =
    [0x0000000000000000,
     0x0000000000000003,
     0x000000000000000f,
     0x000000000000003f,
     0x00000000000000ff,
     0x00000000000003ff,
     0x0000000000000fff,
     0x0000000000003fff,
     0x000000000000ffff,
     0x000000000003ffff,
     0x00000000000fffff,
     0x00000000003fffff,
     0x0000000000ffffff,
     0x0000000003ffffff,
     0x000000000fffffff,
     0x000000003fffffff,
     0x00000000ffffffff,
     0x00000003ffffffff,
     0x0000000fffffffff,
     0x0000003fffffffff,
     0x000000ffffffffff,
     0x000003ffffffffff,
     0x00000fffffffffff,
     0x00003fffffffffff,
     0x0000ffffffffffff,
     0x0003ffffffffffff,
     0x000fffffffffffff,
     0x003fffffffffffff,
     0x00ffffffffffffff,
     0x03ffffffffffffff,
     0x0fffffffffffffff,
     0x3fffffffffffffff,
     0xffffffffffffffff];

// We represent DNA kmers of size less than or equal to 32 by a 64 bit wide
// unsigned integer wrapped in a 'newtype' struct marker.

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Hash)]
pub struct DNAKMer(u64);

pub fn repr(kmer: DNAKMer) -> u64 {
    let DNAKMer(result) = kmer;
    return result;
}

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
    let l = pattern.len();
    assert!(l > 0 && l <= 32);
    let mut mul = 1;
    let mut result = 0;

    for i in 0..l {
        result += mul * dna_to_number(pattern[l - 1 - i]);
        mul *= 4;
    }

    return DNAKMer(result);
}

pub fn pattern_to_canonical(pattern: &[u8]) -> DNAKMer {
    let l = pattern.len();
    assert!(l <= 32);
    let mut mul = 1;
    let mut result0 = 0;
    let mut result1 = 0;
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

pub fn next_dnakmer(kmer: DNAKMer, k: usize, b: u8) -> DNAKMer {
    let DNAKMer(mut result) = kmer;
    result <<= 2;
    result &= KMER_MASKS[k];
    result |= dna_to_number(b);
    return DNAKMer(result); 
}

// This function is useful in printing out a dnakmer, e.g.
// let mut kmer_buf: [u8; 32] = [b'A'; 32];
// dnakmer_into_slice(k, kmer, &mut kmer_buf[0..k]);
// println!("kmer = {}, kmer_string = {}",
//          kmer,
//          String::from_utf8(Vec::from(&kmer_buf[0..k])).unwrap()
// );

pub fn dnakmer_into_slice(k: usize, kmer: DNAKMer, slice: &mut [u8]) {
    let four_power_k = 4u64.pow(2 * k as u32);
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
    curr_kmer: DNAKMer,
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
        if pos + k < self.data.len() {
            let curr_kmer = self.curr_kmer;
            self.curr_kmer = next_dnakmer(curr_kmer,
                                          k,
                                          self.data[pos+k]);
            self.pos += 1;
            return Some(curr_kmer);
        } else {
            return None;
        }
    }
}

// Returns a DNAKMerSeq sequence generator
pub fn dnakmer_seq<'a>(k: usize, s: &'a [u8]) -> DNAKMerSeq<'a> {
    assert!(k <= 32);
    DNAKMerSeq { k: k,
                 data: s,
                 pos: 0,
                 curr_kmer: pattern_to_kmer(&s[0..k])}
}

#[cfg(test)]
mod test {
    use super::pattern_to_canonical;
    use super::pattern_to_kmer;
    use super::DNAKMer;
    use super::dnakmer_seq;

    const KMER_VALS: [u64;11] =
        [0x000000000000008f,
         0x000000000000003c,
         0x00000000000000f1,
         0x00000000000000c4,
         0x0000000000000013,
         0x000000000000004e,
         0x000000000000003b,
         0x00000000000000ec,
         0x00000000000000b0,
         0x00000000000000c3,
         0x000000000000000d];

    #[test]
    fn num_times_test() {
        let mut count = 0;
        let txt = b"GATTACATGTAATC";
        for kmer in dnakmer_seq(4, &txt[0..txt.len()]) {
            let DNAKMer(kmer_int) = kmer;
            assert!(kmer_int == KMER_VALS[count]);
            count += 1;
        }
        assert!(count == 10);
    }

    #[test]
    fn canonical_test() {
        assert!(pattern_to_canonical(b"GATTACA") ==
                pattern_to_canonical(b"TGTAATC"));
        assert!(pattern_to_canonical(b"ACCT") ==
                pattern_to_canonical(b"AGGT"));
        assert!(pattern_to_canonical(b"ACCT") ==
                DNAKMer(23u64));
    }

    #[test]
    fn kmer_test() {
        assert!(pattern_to_canonical(b"GATTACA") ==
                pattern_to_canonical(b"TGTAATC"));
        assert!(pattern_to_kmer(b"GATTGGCT") ==
                pattern_to_kmer(b"GATTGGCT"));
        assert!(pattern_to_canonical(b"ACCT") !=
                pattern_to_kmer(b"GAAC"));
    }

    #[test]
    #[should_panic]
    fn failing_test() {
        assert!(1i32 == 2i32);
    }
}
