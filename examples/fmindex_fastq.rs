// Import some modules
use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use std::io;
// a given text

fn main() {
    let text = b"ACAGCTCGATCGGTA$";
    let pattern = b"ATCG";
    // Create an FM-Index for the given text.
    // instantiate an alphabet
    let alphabet = alphabets::dna::iupac_alphabet();
    // calculate a suffix array
    let sa = suffix_array(text);
    // calculate the Burrows-Wheeler-transform
    let bwt = bwt(text, &sa);
    // calculate the vectors less and Occ (occurrences)
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    // set up FMIndex
    let fmindex = FMIndex::new(&bwt, &less, &occ);
    // do a backwards search for the pattern
    let interval = fmindex.backward_search(pattern.iter());
    let positions = interval.occ(&sa);
    // Iterate over a FASTQ file, use the alphabet to validate read
    // sequences and search for exact matches in the FM-Index.
    // create FASTQ reader
    let mut reader = fastq::Reader::new(io::stdin());
    let mut record = fastq::Record::new();
    reader.read(&mut record).expect("Failed to parse record");
    while !record.is_empty() {
        let check = record.check();
        if check.is_err() {
            panic!("I got a rubbish record!")
        }
        // obtain sequence
        let seq = record.seq();
        // check, whether seq is in the expected alphabet
        if alphabet.is_word(seq) {
            let interval = fmindex.backward_search(seq.iter());
            let positions = interval.occ(&positions);
        }
        reader.read(&mut record).expect("Failed to parse record");
    }
}
