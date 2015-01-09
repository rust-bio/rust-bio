pub fn trim_newline(s: &mut String) {
    if s.ends_with("\n") {
        s.pop();
    }
}


pub fn is_alphabet(seq: &String, alphabet: &str) -> bool {
    for n in seq.graphemes(true) {
        if !alphabet.contains(n) {
            return false;
        }
    }
    true
}


pub fn is_dna_alphabet(seq: &String) -> bool {
    is_alphabet(seq, "ACGTacgt")
}


pub fn is_iupac_alphabet(seq: &String) -> bool {
    is_alphabet(seq, "ACGTURYSWKMBDHVNacgturyswkmbdhvn")
}
