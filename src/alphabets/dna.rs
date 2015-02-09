
use alphabets::Alphabet;


static SYMBOLS: &'static [u8] = b"ACGT";


/// Obtain the DNA alphabet.
pub fn alphabet() -> Alphabet {
    Alphabet::new(b"ACGTacgt")
}


/// Obtain the IUPAC DNA alphabet
pub fn iupac_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTURYSWKMBDHVNacgturyswkmbdhvn")
}


#[derive(Copy)]
pub struct RevComp {
    comp: [u8; 256]
}


impl RevComp {
    pub fn new() -> Self {
        let mut comp = [b'N'; 256];
        for (&a, &b) in SYMBOLS.iter().zip(SYMBOLS.iter().rev()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        RevComp { comp: comp }
    }

    /// Calculate the reverse complement of given text.
    /// The text has to be in DNA alphabet containing only ACGTacgt symbols.
    /// Other symbols will be silently converted into N.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets::dna::RevComp;
    /// let revcomp = RevComp::new();
    /// let text = b"AAACCTT";
    /// let revcomp_text = revcomp.get(text);
    /// assert_eq!(revcomp_text, b"AAGGTTT");
    /// assert_eq!(text, revcomp.get(revcomp_text.as_slice()));
    /// ```
    pub fn get(&self, text: &[u8]) -> Vec<u8> {
        text.iter().rev().map(|&a| self.comp[a as usize]).collect()
    }
}
