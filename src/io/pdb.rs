use lib3dmol::parser;
use lib3dmol::tools;
use std::path::Path;

// Creates a Pdb struct with a field containing the path to the pdb file

pub struct Pdb {
    pdb_file: &'static str,
}

impl Pdb {
    pub fn new(path: &'static str) -> Self {
        Self { 
          pdb_file: Path::new(path).to_str().unwrap()
        }
    }
}

// Amino acid sequence structure that will take the sequence as a String, and the 
// length of the string as a usize
pub struct AASequence {
    sequence: String,
}

impl AASequence {
    // Creates a new AASequence from a pdb file from argument
    pub fn from_pdb(pdb: Pdb) -> Self {
        let my_structure = parser::read_pdb(pdb.pdb_file, &pdb.pdb_file[..]);

        let fasta = tools::fasta_seq(&my_structure);

        let aaseq = AASequence {sequence: fasta};
        println!("{}", aaseq.sequence);
        return aaseq;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PDB: &'static str = "src/io/test_pdb/chignolin.pdb";

    #[test]
    fn test_new_aasequence() {
        let test_pdb = Pdb { pdb_file: PDB };

        let new_sequence = AASequence::from_pdb(test_pdb);
        assert_eq!(new_sequence.sequence, "GYDPETGTWG");
    }
}
