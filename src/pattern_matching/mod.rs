pub mod shift_and;
pub mod kmp;
pub mod bom;
pub mod horspool;
pub mod bndm;


trait PatternMatching {
    fn find_all(&self, text: &[u8]) -> Iterator<Item=usize>;
}


#[cfg(test)]
mod tests {
    use test::Bencher;
    use super::shift_and::ShiftAnd;
    use super::bndm::BNDM;
    use super::kmp::KMP;
    use super::bom::BOM;
    use super::horspool::Horspool;

    static TEXT: &'static [u8] = b"\
ACGGCTAGAAAAGGCTAGGAGTAGGATTCTGCATGCACGACTCGAGCACTAGCACGGGGGGAGGAGTAGGAGATAGAT\
AGAGGATAGATGAT";
    static PATTERN: &'static [u8] = b"GGATTCTGCA";


    #[bench]
    fn bench_shift_and(b: &mut Bencher) {
        let shiftand = ShiftAnd::new(PATTERN);
        b.iter(|| {
            shiftand.find_all(TEXT).collect::<Vec<usize>>()
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
}
