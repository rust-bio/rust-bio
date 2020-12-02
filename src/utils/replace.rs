// Copyright 2020 Christopher Sugai.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! 
//! Traits to replace Ns with pseudorandom nucleotide bases and convert U to/from T.
//! The methods in this trait assumes that:
//! - The sequence contains only ASCII characters.
//! - Lowercase and uppercase characters are not mixed.
//! # Examples
//! ```
//! extern crate bio;
//! use rand::seq::SliceRandom;
//! use rand::Rng;
//! use std::string::String;
//! 
//! let mut seq = *b"ATGCT";
//! seq.to_rna();
//! assert_eq!(seq, *b"AUGCU");
//! 
//! let mut rng = rand::thread_rng();
//! let mut seq = *b"ATGCNNNN";
//! seq.replace_n(&mut rng);
//! for c in seq.iter_mut() {
//!     STANDARD_DNA_NUCLEOTIDES.contains(c);
//! }
//! 
//! ```

extern crate bio;
use rand::seq::SliceRandom;
use rand::Rng;
use std::string::String;

const STANDARD_DNA_NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const STANDARD_RNA_NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'U'];

pub trait DnaSeq: Seq {
    /// Replace `T` with `U`
    fn to_rna(&mut self);
    /// Fill N with pseudorandom nucleotides ACTG and n with actg
    fn replace_n<R: Rng>(&mut self, rng: &mut R);
}

pub trait RnaSeq: Seq {
    /// Replace `U` with `T`
    fn to_dna(&mut self);
    /// Fill N with pseudorandom nucleotides ACTG and n with actg
    fn replace_n<R: Rng>(&mut self, rng: &mut R);
}

impl<T> DnaSeq for T
where
    T: AsMut<[u8]>,
{
    fn to_rna(&mut self) {
        self.as_mut().iter_mut().for_each(|c| match c {
            b'T' => *c = b'U',
            _ => {}
        })
    }
    fn replace_n<R: Rng>(&mut self, rng: &mut R) {
        self.as_mut().iter_mut().for_each(|c| match c {
            b'N' => *c = *STANDARD_DNA_NUCLEOTIDES.choose(rng).unwrap(),
            _ => {}
        })
    }
}

impl<T> RnaSeq for T
where
    T: AsMut<[u8]>,
{
    fn to_dna(&mut self) {
        self.as_mut().iter_mut().for_each(|c| match c {
            b'U' => *c = b'T',
            _ => {}
        })
    }
    fn replace_n<R: Rng>(&mut self, rng: &mut R) {
        self.as_mut().iter_mut().for_each(|c| match c {
            b'N' => *c = *STANDARD_RNA_NUCLEOTIDES.choose(rng).unwrap(),
            _ => {}
        })
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_seq_to_rna() {
        let mut seq = *b"ATGCT";
        seq.to_rna();
        assert_eq!(seq, *b"AUGCU");
    }
    #[test]
    fn test_seq_replace() {
        let mut rng = rand::thread_rng();
        let mut seq = *b"ATGCNNNN";
        seq.replace_n(&mut rng);
        for c in seq.iter_mut() {
            STANDARD_DNA_NUCLEOTIDES.contains(c);
        }
    }
}
