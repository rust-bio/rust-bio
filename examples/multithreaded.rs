use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use std::sync::Arc;
use std::thread;
fn main() {
    let text = b"ACGGATGCTGGATCGGATCGCGCTAGCTA$";
    let patterns = vec![b"ACCG", b"TGCT"];
    // Create an FM-Index for a given text.
    let alphabet = alphabets::dna::iupac_alphabet();
    let sa = suffix_array(text);
    let bwt = Arc::new(bwt(text, &sa));
    let less = Arc::new(less(bwt.as_ref(), &alphabet));
    let occ = Arc::new(Occ::new(bwt.as_ref(), 3, &alphabet));
    let fmindex = Arc::new(FMIndex::new(bwt, less, occ));
    // Spawn threads to perform backward searches for each interval
    let interval_calculators = patterns
        .into_iter()
        .map(|pattern| {
            let fmindex = fmindex.clone();
            thread::spawn(move || fmindex.backward_search(pattern.iter()))
        })
        .collect::<Vec<_>>();
    // Loop through the results, extracting the positions array for each pattern
    for interval_calculator in interval_calculators {
        let positions = interval_calculator.join().unwrap().occ(&sa);
    }
}
