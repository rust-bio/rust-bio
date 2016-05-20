
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

use alphabets::Alphabet;


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
    Alphabet::new(b"ACGTURYSWKMBDHVNacgturyswkmbdhvn")
}


/// Implementation of transformation into reverse complement.
#[cfg_attr(feature = "serde_macros", derive(Serialize, Deserialize))]
pub struct RevComp {
    comp: Vec<u8>,
}


impl RevComp {
    /// Create a new instance of reverse complement algorithm.
    pub fn new() -> Self {
        let mut comp = Vec::new();
        comp.resize(256, 0);
        for (v, mut a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        RevComp { comp: comp }
    }

    pub fn comp(&self, a: u8) -> u8 {
        self.comp[a as usize]
    }

    /// Calculate the reverse complement of given text.
    /// The text has to be in DNA alphabet containing only ACGTacgt symbols.
    /// Other symbols won't be converted.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets::dna::RevComp;
    /// let revcomp = RevComp::new();
    /// let text = b"AAACCTTW";
    /// let revcomp_text = revcomp.get(text);
    /// assert_eq!(revcomp_text, &b"WAAGGTTT"[..]);
    /// assert_eq!(revcomp.get(&revcomp_text[..]), &text[..]);
    /// ```
    pub fn get(&self, text: &[u8]) -> Vec<u8> {
        text.iter().rev().map(|&a| self.comp(a)).collect()
    }
}
