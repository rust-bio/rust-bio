use lib3dmol::parser;
use lib3dmol::tools;
use std::env;

// Amino acid sequence structure that will take the sequence as a String, and the 
// length of the string as a usize
pub struct AASequence {
    sequence: String,
}

impl AASequence {
    // Creates a new AASequence from a pdb file from argument
    pub fn new() -> Self {
        let args: Vec<String> = env::args().collect();
            let pdb = &args[1];

            let my_structure = parser::read_pdb(pdb, &pdb[..]);

            let fasta = tools::fasta_seq(&my_structure);

            let aaseq = AASequence {sequence: fasta};
            println!("{}", aaseq.sequence);
            return aaseq;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lib3dmol::parser;
    use lib3dmol::tools;

    const PDB: &str = "src/io/test_pdb/chignolin.pdb";
    const PDB_NAME: &str = "chignolin.pdb";

    #[test]
    fn test_new_aasequence() {
        let my_structure = parser::read_pdb(PDB, PDB_NAME);
        let fasta = tools::fasta_seq(&my_structure);
        let new_sequence = AASequence {sequence: fasta};
        assert_eq!(new_sequence.sequence, "GYDPETGTWG");
    }
}
