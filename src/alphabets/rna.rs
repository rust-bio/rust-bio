
// Copyright 2017 Ryan Hagenson.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Implementation of the RNA alphabet.
//!
//! # Example
//!
//! ```
//! use bio::alphabets;
//! let alphabet = alphabets::rna::alphabet();
//! assert!(alphabet.is_word(b"GAUUACA"));
//! assert!(alphabet.is_word(b"gauuaca"));
//! assert!(!alphabet.is_word(b"ACGT"));
//! ```

use alphabets::Alphabet;
use utils::IntoTextIterator;


/// The RNA alphabet (uppercase and lowercase).
pub fn alphabet() -> Alphabet {
    Alphabet::new(b"ACGUacgu")
}


/// The RNA alphabet including N (uppercase and lowercase).
pub fn n_alphabet() -> Alphabet {
    Alphabet::new(b"ACGUNacgun")
}


/// The IUPAC RNA alphabet (uppercase and lowercase).
pub fn iupac_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTURYSWKMBDHVNZacgturyswkmbdhvnz")
}


lazy_static! {
    static ref COMPLEMENT: Vec<u8> = {
        let mut comp = Vec::new();
        comp.resize(256, 0);
        for (v, mut a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCUYRWSKMDVHBNZ".iter().zip(b"UCGARYSWMKHBDVNZ".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}


/// Return complement of given RNA alphabet character (IUPAC alphabet supported).
pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}


/// Calculate reverse complement of given text (IUPAC alphabet supported).
pub fn revcomp<'a, T: IntoTextIterator<'a>>(text: T) -> Vec<u8>
    where T::IntoIter: DoubleEndedIterator
{
    text.into_iter().rev().map(|&a| complement(a)).collect()
}

#[cfg(test)]
mod tests {
    fn validate_height(node: &Node<i64, String>) {
        let left_height = node.left.as_ref().map_or(0, |n| n.height);
        let right_height = node.right.as_ref().map_or(0, |n| n.height);
        assert!((left_height - right_height).abs() <= 1);
        assert_eq!(node.height, cmp::max(left_height, right_height) + 1)
    }
}
