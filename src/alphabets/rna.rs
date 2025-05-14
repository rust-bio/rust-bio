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

// You may generate this complement by the code here, minus the lazy_static part.
// This is slightly different from the COMPLEMENT IN DNA.
// lazy_static! {
//     static ref COMPLEMENT: [u8; 256] = {
//         let mut comp = [0; 256];
//         for (v, a) in comp.iter_mut().enumerate() {
//             *a = v as u8;
//         }
//         for (&a, &b) in b"AGCUYRWSKMDVHBNZ".iter().zip(b"UCGARYWSMKHBDVNZ".iter()) {
//             comp[a as usize] = b;
//             comp[a as usize + 32] = b + 32;  // lowercase variants
//         }
//         comp
//     };
// }
pub const COMPLEMENT: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 85, 86, 71, 72, 69, 70, 67, 68, 73,
    74, 77, 76, 75, 78, 79, 80, 81, 89, 83, 84, 65, 66, 87, 88, 82, 90, 91, 92, 93, 94, 95, 96,
    117, 118, 103, 104, 101, 102, 99, 100, 105, 106, 109, 108, 107, 110, 111, 112, 113, 121, 115,
    116, 97, 98, 119, 120, 114, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
    154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
    173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
    211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
    230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248,
    249, 250, 251, 252, 253, 254, 255,
];

/// Return complement of given RNA alphabet character (IUPAC alphabet supported).
///
/// Casing of input character is preserved, e.g. `u` → `a`, but `U` → `A`.
/// All `N`s and `Z`s remain as they are.
///
/// ```
/// use bio::alphabets::rna;
///
/// assert_eq!(rna::complement(65), 85); // A → U
/// assert_eq!(rna::complement(103), 99); // g → c
/// assert_eq!(rna::complement(89), 82); // Y → R
/// assert_eq!(rna::complement(115), 115); // s → s
/// assert_eq!(rna::complement(78), 78); // N → N
/// ```
#[inline]
pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

/// Calculate reverse complement of given text (IUPAC alphabet supported).
///
/// Casing of characters is preserved, e.g. `b"uAGg"` → `b"cCUa"`.
/// All `N`s and `Z`s remain as they are.
///
/// ```
/// use bio::alphabets::rna;
///
/// assert_eq!(rna::revcomp(b"ACGUN"), b"NACGU");
/// assert_eq!(rna::revcomp(b"GaUuaCA"), b"UGuaAuC");
/// assert_eq!(rna::revcomp(b"AGCUYRWSKMDVHBNZ"), b"ZNVDBHKMSWYRAGCU");
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
