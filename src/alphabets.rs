
/// A module to handle different alphabets
///
/// # Example
///
/// ```rust
/// use bio::alphabets;
/// let alphabet = alphabets::get_dna_alphabet();
/// assert!(alphabet.is_word(b"AACCTgga"));
/// assert!(!alphabet.is_word(b"AXYZ"));
/// ```

use std::collections::{BitvSet, VecMap};

/// Representation of an alphabet.
pub struct Alphabet {
    symbols: BitvSet
}


impl Alphabet {
    pub fn new(letters: &[u8]) -> Self {
        let mut s = BitvSet::new();
        s.extend(letters.iter().map(|&c| c as usize));

        Alphabet { symbols: s }
    }

    pub fn is_word(&self, text: &[u8]) -> bool {
        text.iter().all(|&c| self.symbols.contains(&(c as usize)))
    }

    pub fn get_ranks(&self) -> VecMap<usize> {
        let mut ranks = VecMap::new();
        for (r, c) in self.symbols.iter().enumerate() {
            ranks.insert(c, r);
        }

        ranks
    }

    pub fn len(&self) -> usize {
        self.symbols.len()
    }
}


/// Obtain the DNA alphabet.
pub fn get_dna_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTacgt")
}


/// Obtain the IUPAC DNA alphabet
pub fn get_iupac_dna_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTURYSWKMBDHVNacgturyswkmbdhvn")
}
