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
use crate::alphabets::Alphabet;
use std::borrow::Borrow;
use std::sync::LazyLock;

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

static COMPLEMENT: LazyLock<[u8; 256]> = LazyLock::new(|| {
    let mut comp = [0; 256];
    comp.iter_mut().enumerate().for_each(|(v, a)| {
        *a = v as u8;
    });
    b"AGCTYRWSKMDVHBN"
        .iter()
        .zip(b"TCGARYWSMKHBDVN".iter())
        .for_each(|(&a, &b)| {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;
        });
    comp
});

/// Return complement of given DNA alphabet character (IUPAC alphabet supported).
///
/// Casing of input character is preserved, e.g. `t` → `a`, but `T` → `A`.
/// All `N`s remain as they are.
///
/// ```
/// use bio::alphabets::dna;
///
/// assert_eq!(dna::complement(65), 84); // A → T
/// assert_eq!(dna::complement(99), 103); // c → g
/// assert_eq!(dna::complement(78), 78); // N → N
/// assert_eq!(dna::complement(89), 82); // Y → R
/// assert_eq!(dna::complement(115), 115); // s → s
/// ```
#[inline]
pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

/// Calculate reverse complement of given text (IUPAC alphabet supported).
///
/// Casing of characters is preserved, e.g. `b"NaCgT"` → `b"aCgTN"`.
/// All `N`s remain as they are.
///
/// ```
/// use bio::alphabets::dna;
///
/// assert_eq!(dna::revcomp(b"ACGTN"), b"NACGT");
/// assert_eq!(dna::revcomp(b"GaTtaCA"), b"TGtaAtC");
/// assert_eq!(dna::revcomp(b"AGCTYRWSKMDVHBN"), b"NVDBHKMSWYRAGCT");
/// ```
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
        assert!(alphabet().is_word(b"GATTACA"));
    }

    #[test]
    fn is_no_word() {
        assert!(!alphabet().is_word(b"gaUUaca"));
    }

    #[test]
    fn symbol_is_no_word() {
        assert!(!alphabet().is_word(b"#"));
    }

    #[test]
    fn number_is_no_word() {
        assert!(!alphabet().is_word(b"42"));
    }
}
