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

// You may generate this complement by the code here, minus the lazy_static part.
// lazy_static! {
//     static ref COMPLEMENT: [u8; 256] = {
//         let mut comp = [0; 256];
//         for (v, a) in comp.iter_mut().enumerate() {
//             *a = v as u8;
//         }
//         for (&a, &b) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
//             comp[a as usize] = b;
//             comp[a as usize + 32] = b + 32;  // lowercase variants
//         }
//         comp
//     };
// }
pub const COMPLEMENT: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 86, 71, 72, 69, 70, 67, 68, 73,
    74, 77, 76, 75, 78, 79, 80, 81, 89, 83, 65, 85, 66, 87, 88, 82, 90, 91, 92, 93, 94, 95, 96,
    116, 118, 103, 104, 101, 102, 99, 100, 105, 106, 109, 108, 107, 110, 111, 112, 113, 121, 115,
    97, 117, 98, 119, 120, 114, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
    154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
    173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
    211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
    230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248,
    249, 250, 251, 252, 253, 254, 255,
];

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
