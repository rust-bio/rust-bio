#![feature(test)]

extern crate test;

use bio::io::fastq::{Reader,Record};
use test::{Bencher,black_box};
use std::io;

#[bench]
fn iterate_over_fastq_file(b: &mut Bencher) {
    b.iter(|| {
        let reader = Reader::new(FASTQ_FILE);
        black_box(reader.records().collect::<Vec<io::Result<Record>>>());
    });
}


const MULTI_FASTQ_FILE: &'static [u8] = b"@read0 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read1 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read2 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read3 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read4 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read5 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read6 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read7 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read8 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read9 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read10 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read11 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read12 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read13 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read14 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read15 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read16 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read17 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read18 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read19 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read20 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read21 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read22 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read23 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read24 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read25 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read26 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read27 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read28 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read29 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read30 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read31 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read32 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read33 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read34 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read35 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read36 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read37 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read38 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read39 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read40 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read41 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read42 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read43 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read44 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read45 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read46 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read47 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read48 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read49 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read50 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read51 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read52 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read53 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read54 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read55 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read56 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read57 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read58 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read59 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read60 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read61 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read62 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read63 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read64 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read65 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read66 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read67 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read68 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read69 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read70 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read71 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read72 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read73 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read74 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read75 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read76 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read77 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read78 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read79 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read80 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read81 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read82 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read83 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read84 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read85 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read86 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read87 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read88 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read89 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read90 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read91 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read92 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read93 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read94 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read95 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read96 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read97 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read98 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read99 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read100 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read101 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read102 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read103 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read104 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read105 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read106 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read107 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read108 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read109 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read110 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read111 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read112 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read113 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read114 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read115 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read116 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read117 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read118 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read119 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read120 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read121 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read122 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read123 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read124 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read125 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read126 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read127 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read128 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read129 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read130 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read131 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read132 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read133 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read134 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read135 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read136 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read137 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read138 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read139 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read140 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read141 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read142 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read143 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read144 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read145 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read146 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read147 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read148 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read149 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read150 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read151 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read152 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read153 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read154 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read155 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read156 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read157 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read158 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read159 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read160 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read161 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read162 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read163 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read164 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read165 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read166 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read167 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read168 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read169 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read170 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read171 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read172 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read173 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read174 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read175 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read176 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read177 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read178 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read179 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read180 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read181 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read182 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read183 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read184 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read185 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read186 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read187 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read188 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read189 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read190 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read191 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read192 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read193 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read194 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read195 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read196 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read197 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read198 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read199 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read200 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read201 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read202 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read203 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read204 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read205 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read206 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read207 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read208 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read209 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read210 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read211 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read212 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read213 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read214 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read215 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read216 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read217 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read218 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read219 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read220 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read221 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read222 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read223 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read224 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read225 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read226 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read227 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read228 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read229 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read230 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read231 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read232 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read233 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read234 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read235 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read236 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read237 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read238 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read239 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read240 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read241 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read242 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read243 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read244 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read245 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read246 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read247 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read248 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read249 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read250 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read251 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read252 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read253 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read254 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read255 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read256 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read257 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read258 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read259 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read260 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read261 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read262 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read263 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read264 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read265 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read266 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read267 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read268 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read269 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read270 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read271 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read272 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read273 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read274 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read275 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read276 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read277 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read278 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read279 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read280 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read281 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read282 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read283 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read284 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read285 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read286 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read287 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read288 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read289 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read290 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read291 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read292 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read293 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read294 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read295 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read296 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read297 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read298 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read299 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read300 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read301 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read302 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read303 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read304 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read305 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read306 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read307 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read308 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read309 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read310 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read311 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read312 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read313 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read314 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read315 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read316 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read317 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read318 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read319 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read320 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read321 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read322 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read323 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read324 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read325 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read326 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read327 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read328 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read329 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read330 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read331 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read332 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read333 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read334 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read335 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read336 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read337 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read338 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read339 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read340 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read341 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read342 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read343 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read344 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read345 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read346 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read347 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read348 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read349 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read350 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read351 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read352 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read353 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read354 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read355 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read356 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read357 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read358 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read359 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read360 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read361 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read362 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read363 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read364 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read365 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read366 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read367 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read368 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read369 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read370 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read371 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read372 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read373 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read374 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read375 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read376 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read377 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read378 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read379 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read380 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read381 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read382 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read383 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read384 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read385 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read386 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read387 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read388 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read389 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read390 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read391 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read392 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read393 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read394 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read395 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read396 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read397 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read398 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read399 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read400 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read401 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read402 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read403 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read404 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read405 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read406 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read407 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read408 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read409 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read410 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read411 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read412 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read413 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read414 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read415 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read416 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read417 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read418 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read419 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read420 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read421 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read422 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read423 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read424 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read425 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read426 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read427 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read428 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read429 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read430 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read431 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read432 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read433 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read434 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read435 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read436 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read437 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read438 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read439 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read440 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read441 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read442 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read443 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read444 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read445 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read446 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read447 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read448 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read449 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read450 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read451 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read452 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read453 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read454 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read455 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read456 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read457 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read458 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read459 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read460 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read461 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read462 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read463 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read464 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read465 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read466 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read467 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read468 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read469 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read470 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read471 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read472 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read473 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read474 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read475 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read476 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read477 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read478 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read479 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read480 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read481 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read482 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read483 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read484 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read485 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read486 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read487 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read488 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read489 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read490 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read491 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read492 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read493 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read494 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read495 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read496 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read497 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read498 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
@read499 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
GGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@$$!@$$!@$$!@$$!@$$!
";

const FASTQ_FILE: &'static [u8] = b"@read0 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read1 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read2 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read3 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read4 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read5 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read6 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read7 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read8 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read9 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read10 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read11 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read12 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read13 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read14 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read15 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read16 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read17 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read18 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read19 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read20 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read21 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read22 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read23 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read24 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read25 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read26 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read27 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read28 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read29 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read30 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read31 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read32 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read33 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read34 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read35 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read36 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read37 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read38 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read39 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read40 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read41 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read42 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read43 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read44 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read45 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read46 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read47 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read48 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read49 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read50 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read51 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read52 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read53 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read54 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read55 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read56 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read57 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read58 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read59 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read60 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read61 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read62 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read63 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read64 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read65 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read66 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read67 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read68 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read69 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read70 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read71 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read72 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read73 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read74 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read75 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read76 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read77 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read78 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read79 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read80 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read81 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read82 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read83 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read84 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read85 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read86 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read87 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read88 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read89 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read90 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read91 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read92 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read93 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read94 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read95 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read96 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read97 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read98 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read99 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read100 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read101 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read102 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read103 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read104 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read105 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read106 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read107 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read108 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read109 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read110 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read111 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read112 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read113 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read114 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read115 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read116 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read117 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read118 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read119 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read120 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read121 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read122 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read123 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read124 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read125 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read126 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read127 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read128 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read129 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read130 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read131 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read132 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read133 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read134 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read135 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read136 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read137 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read138 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read139 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read140 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read141 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read142 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read143 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read144 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read145 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read146 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read147 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read148 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read149 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read150 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read151 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read152 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read153 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read154 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read155 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read156 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read157 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read158 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read159 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read160 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read161 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read162 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read163 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read164 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read165 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read166 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read167 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read168 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read169 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read170 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read171 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read172 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read173 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read174 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read175 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read176 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read177 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read178 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read179 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read180 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read181 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read182 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read183 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read184 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read185 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read186 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read187 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read188 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read189 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read190 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read191 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read192 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read193 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read194 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read195 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read196 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read197 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read198 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read199 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read200 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read201 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read202 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read203 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read204 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read205 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read206 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read207 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read208 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read209 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read210 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read211 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read212 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read213 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read214 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read215 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read216 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read217 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read218 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read219 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read220 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read221 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read222 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read223 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read224 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read225 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read226 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read227 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read228 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read229 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read230 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read231 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read232 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read233 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read234 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read235 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read236 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read237 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read238 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read239 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read240 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read241 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read242 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read243 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read244 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read245 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read246 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read247 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read248 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read249 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read250 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read251 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read252 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read253 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read254 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read255 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read256 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read257 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read258 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read259 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read260 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read261 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read262 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read263 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read264 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read265 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read266 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read267 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read268 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read269 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read270 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read271 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read272 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read273 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read274 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read275 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read276 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read277 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read278 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read279 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read280 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read281 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read282 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read283 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read284 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read285 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read286 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read287 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read288 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read289 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read290 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read291 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read292 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read293 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read294 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read295 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read296 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read297 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read298 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read299 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read300 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read301 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read302 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read303 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read304 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read305 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read306 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read307 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read308 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read309 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read310 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read311 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read312 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read313 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read314 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read315 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read316 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read317 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read318 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read319 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read320 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read321 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read322 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read323 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read324 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read325 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read326 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read327 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read328 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read329 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read330 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read331 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read332 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read333 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read334 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read335 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read336 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read337 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read338 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read339 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read340 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read341 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read342 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read343 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read344 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read345 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read346 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read347 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read348 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read349 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read350 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read351 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read352 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read353 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read354 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read355 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read356 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read357 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read358 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read359 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read360 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read361 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read362 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read363 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read364 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read365 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read366 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read367 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read368 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read369 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read370 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read371 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read372 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read373 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read374 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read375 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read376 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read377 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read378 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read379 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read380 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read381 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read382 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read383 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read384 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read385 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read386 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read387 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read388 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read389 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read390 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read391 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read392 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read393 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read394 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read395 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read396 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read397 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read398 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read399 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read400 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read401 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read402 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read403 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read404 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read405 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read406 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read407 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read408 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read409 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read410 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read411 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read412 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read413 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read414 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read415 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read416 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read417 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read418 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read419 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read420 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read421 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read422 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read423 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read424 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read425 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read426 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read427 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read428 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read429 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read430 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read431 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read432 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read433 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read434 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read435 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read436 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read437 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read438 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read439 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read440 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read441 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read442 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read443 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read444 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read445 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read446 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read447 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read448 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read449 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read450 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read451 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read452 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read453 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read454 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read455 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read456 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read457 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read458 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read459 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read460 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read461 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read462 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read463 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read464 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read465 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read466 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read467 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read468 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read469 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read470 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read471 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read472 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read473 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read474 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read475 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read476 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read477 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read478 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read479 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read480 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read481 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read482 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read483 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read484 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read485 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read486 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read487 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read488 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read489 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read490 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read491 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read492 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read493 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read494 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read495 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read496 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read497 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read498 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
@read499 description
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGGGGGGGGGGGGGGGG
+
@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!@$$!
";
