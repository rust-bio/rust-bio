
/// A module to handle different alphabets
///
/// # Example
///
/// ```rust
/// use bio::alphabets;
/// let alphabet = alphabets::dna_alphabet();
/// assert!(alphabet.is_word(b"AACCT"));
/// assert!(!alphabet.is_word(b"AXYZ"));
/// ```

use std::collections::BitvSet;

/// Representation of an alphabet.
pub struct Alphabet {
    set: BitvSet
}


impl Alphabet {
    pub fn new(letters: &[u8]) -> Self {
        let mut s = BitvSet::new();
        s.extend(letters.iter().map(|&c| c as usize));

        Alphabet { set: s }
    }

    pub fn is_word(&self, text: &[u8]) -> bool {
        text.iter().all(|&c| self.set.contains(&(c as usize)))
    }
}


/// Obtain the DNA alphabet.
pub fn dna_alphabet() -> Alphabet {
    Alphabet::new(b"ACGT")
}


/// Obtain the IUPAC DNA alphabet
pub fn iupac_dna_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTURYSWKMBDHVN")
}
