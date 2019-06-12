// Copyright 2014-2015 Johannes Köster, Peer Aramillo Irizar.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Implementation of the DNA alphabet.
//!
//! # Example
//!
//! ```
//! use bio::alphabets;
//! let alphabet = alphabets::dna::alphabet();
//! assert!(alphabet.is_word(b"GATTACA"));
//! assert!(alphabet.is_word(b"gattaca"));
//! assert!(!alphabet.is_word(b"ACGU"));
//! ```

use std::borrow::Borrow;

use crate::alphabets::Alphabet;

/// The DNA alphabet (uppercase and lowercase).
pub fn alphabet() -> Alphabet {
    Alphabet::new(b"ACGTacgt")
}

/// The DNA alphabet including N (uppercase and lowercase).
pub fn n_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTNacgtn")
}

/// The IUPAC DNA alphabet (uppercase and lowercase).
pub fn iupac_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTRYSWKMBDHVNZacgtryswkmbdhvnz")
}

lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

/// Return complement of given DNA alphabet character (IUPAC alphabet supported).
pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

/// Calculate reverse complement of given text (IUPAC alphabet supported).
pub fn revcomp<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}
