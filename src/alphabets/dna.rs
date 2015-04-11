
use alphabets::Alphabet;


static SYMBOLS: &'static [u8] = b"ACGT";


/// Obtain the DNA alphabet.
pub fn alphabet() -> Alphabet {
    Alphabet::new(b"ACGTacgt")
}


pub fn n_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTNacgtn")
}


/// Obtain the IUPAC DNA alphabet
pub fn iupac_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTURYSWKMBDHVNacgturyswkmbdhvn")
}


pub struct RevComp {
    comp: [u8; 256]
}


impl RevComp {
    pub fn new() -> Self {
        let mut comp = [0u8; 256];
        for a in 0..256 {
            comp[a] = a as u8;
        }
        for (&a, &b) in SYMBOLS.iter().zip(SYMBOLS.iter().rev()) {
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
    /// #![feature(convert)]
    /// use bio::alphabets::dna::RevComp;
    /// let revcomp = RevComp::new();
    /// let text = b"AAACCTT";
    /// let revcomp_text = revcomp.get(text);
    /// assert_eq!(revcomp_text, &b"AAGGTTT"[..]);
    /// assert_eq!(revcomp.get(revcomp_text.as_slice()), &text[..]);
    /// ```
    pub fn get(&self, text: &[u8]) -> Vec<u8> {
        text.iter().rev().map(|&a| self.comp(a)).collect()
    }
}
