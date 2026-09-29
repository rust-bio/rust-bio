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
//!
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

/// Bit mask of the bases an IUPAC code can denote.
pub fn iupac_mask(a: u8) -> u8 {
    const A: u8 = 0b0001;
    const C: u8 = 0b0010;
    const G: u8 = 0b0100;
    const T: u8 = 0b1000;
    match a.to_ascii_uppercase() {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T,
        b'R' => A | G,
        b'Y' => C | T,
        b'S' => C | G,
        b'W' => A | T,
        b'K' => G | T,
        b'M' => A | C,
        b'B' => C | G | T,
        b'D' => A | G | T,
        b'H' => A | C | T,
        b'V' => A | C | G,
        b'N' => A | C | G | T,
        _ => 0,
    }
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

    const IUPAC_CODES: [(u8, u8); 15] = [
        (b'A', 0b0001),
        (b'C', 0b0010),
        (b'G', 0b0100),
        (b'T', 0b1000),
        (b'R', 0b0101), // A | G = 0001 | 0100
        (b'Y', 0b1010), // C | T = 0010 | 1000
        (b'S', 0b0110), // C | G
        (b'W', 0b1001), // A | T
        (b'K', 0b1100), // G | T
        (b'M', 0b0011), // A | C
        (b'B', 0b1110), // C | G | T
        (b'D', 0b1101), // A | G | T
        (b'H', 0b1011), // A | C | T
        (b'V', 0b0111), // A | C | G
        (b'N', 0b1111), // A | C | G | T
    ];

    #[test]
    fn iupac_mask_is_the_union_of_the_denoted_bases() {
        for (code, expected) in IUPAC_CODES {
            assert_eq!(iupac_mask(code), expected);
        }
    }

    #[test]
    fn iupac_mask_of_unknown_symbol_is_zero() {
        for symbol in [b'Z', b'-', b'*', b' ', 0, 255] {
            assert_eq!(iupac_mask(symbol), 0);
        }
    }
}
