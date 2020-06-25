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

use std::borrow::Borrow;

use crate::alphabets::Alphabet;

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
    Alphabet::new(b"ACGURYSWKMBDHVNZacguryswkmbdhvnz")
}

lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCUYRWSKMDVHBNZ".iter().zip(b"UCGARYWSMKHBDVNZ".iter()) {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_word() {
        assert!(alphabet().is_word(b"GAUUACA"));
    }

    #[test]
    fn is_no_word() {
        assert!(!alphabet().is_word(b"gaTTaca"));
    }

    #[test]
    fn symbol_is_no_word() {
        assert!(!alphabet().is_word(b"#"));
    }

    #[test]
    fn number_is_no_word() {
        assert!(!alphabet().is_word(b"42"));
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(revcomp(b"GAUUACA"), b"UGUAAUC");
    }
}
