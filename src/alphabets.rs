//! Handling different alphabets.
//!
//! # Example
//!
//! ```rust
//! use bio::alphabets;
//! let alphabet = alphabets::get_dna_alphabet();
//! assert!(alphabet.is_word(b"AACCTgga"));
//! assert!(!alphabet.is_word(b"AXYZ"));
//! ```


use std::collections::{BitvSet, VecMap};


pub type SymbolRanks = VecMap<u8>;


/// Representation of an alphabet.
pub struct Alphabet {
    pub symbols: BitvSet
}


impl Alphabet {
    pub fn new(symbols: &[u8]) -> Self {
        Alphabet::from_iter(symbols.iter().map(|&c| c))
    }

    pub fn from_iter<I: Iterator<Item=u8>>(symbols: I) -> Self {
        let mut s = BitvSet::new();
        s.extend(symbols.map(|c| c as usize));

        Alphabet { symbols: s }
    }

    pub fn is_word(&self, text: &[u8]) -> bool {
        text.iter().all(|&c| self.symbols.contains(&(c as usize)))
    }

    pub fn max_symbol(&self) -> Option<u8> {
        self.symbols.iter().max().map(|c| c as u8)
    }

    pub fn len(&self) -> usize {
        self.symbols.len()
    }
}


pub struct RankTransform {
    pub ranks: SymbolRanks
}


impl RankTransform {
    pub fn new(alphabet: &Alphabet) -> Self {
        let mut ranks = VecMap::new();
        for (r, c) in alphabet.symbols.iter().enumerate() {
            ranks.insert(c, r as u8);
        }

        RankTransform { ranks: ranks }
    }

    pub fn transform(&self, text: &[u8]) -> Vec<u8> {
        text.iter()
            .map(|&c| *self.ranks.get(&(c as usize)).expect("Unexpected character in text."))
            .collect()
    }

    pub fn get_alphabet(&self) -> Alphabet {
        Alphabet::from_iter(self.ranks.keys().map(|c| c as u8))
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
