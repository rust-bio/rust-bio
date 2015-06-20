// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This module contains various useful pattern matching algorithms.
//! The implementations are based on the lecture notes
//! "Algorithmen auf Sequenzen", Kopczynski, Marschall, Martin and Rahmann, 2008 - 2015.
//!
//! * Algorithm of Horspool: fastest for a sufficiently large alphabet
//! * Shift And algorithm: fast for patterns with less than 64 symbols and very small alphabets.
//! * BNDM algorithm: fast for patterns with less than 64 symbols.
//! * BOM algorithm: fast for long patterns and small alphabet.
//! * KMP algorithm: the classical ancestor.


pub mod shift_and;
pub mod kmp;
pub mod bom;
pub mod horspool;
pub mod bndm;
pub mod ukkonen;
pub mod myers;


#[cfg(test)]
mod tests {
    use test::Bencher;
    use super::shift_and::ShiftAnd;
    use super::bndm::BNDM;
    use super::kmp::KMP;
    use super::bom::BOM;
    use super::horspool::Horspool;

    static TEXT: &'static [u8] = b"ACGGCTAGAAAAGGCTAGGAGTAGGATTCTGCATGCACGACT\
CGAGCACTAGCACGGGGGGAGGAGTAGGAGATAGATAGAGGATAGATGATACGGCTAGAAAAGGCTA\
GGAGTAGGATTCTGCATGCACGACTCGAGCACTAGCACGGGGGGAGGAGTAGGAGATAGATAGAGGA\
TAGATGATACGGCTAGAAAAGGCTAGGAGTAGGATTCTGCATGCACGACTCGAGCACTAGCACGGGG\
GGAGGAGTAGGAGATAGATAGAGGATAGATGATACGGCTAGAAAAGGCTAGGAGTAGGATTCTGCAT\
GCACGACTCGAGCACTAGCACGGGGGGAGGAGTAGGAGATAGATAGAGGATAGATGATACGGCTAGA\
AAAGGCTAGGAGTAGGATTCTGCATGCACGACTCGAGCACTAGCACGGGGGGAGGAGTAGGAGATAG\
ATAGAGGATAGATGAT";
    static PATTERN: &'static [u8] = b"GGATTCTGCA";


    #[bench]
    fn bench_shift_and(b: &mut Bencher) {
        let shiftand = ShiftAnd::new(PATTERN);
        b.iter(|| {
            shiftand.find_all(TEXT.iter()).collect::<Vec<usize>>()
        });
    }

    #[bench]
    fn bench_bndm(b: &mut Bencher) {
        let bndm = BNDM::new(PATTERN);
        b.iter(|| {
            bndm.find_all(TEXT).collect::<Vec<usize>>()
        });
    }

    #[bench]
    fn bench_bom(b: &mut Bencher) {
        let bom = BOM::new(PATTERN);
        b.iter(|| {
            bom.find_all(TEXT).collect::<Vec<usize>>()
        });
    }

    #[bench]
    fn bench_horspool(b: &mut Bencher) {
        let horspool = Horspool::new(PATTERN);
        b.iter(|| {
            horspool.find_all(TEXT).collect::<Vec<usize>>()
        });
    }

    #[bench]
    fn bench_kmp(b: &mut Bencher) {
        let kmp = KMP::new(PATTERN);
        b.iter(|| {
            kmp.find_all(TEXT.iter()).collect::<Vec<usize>>()
        });
    }
}
