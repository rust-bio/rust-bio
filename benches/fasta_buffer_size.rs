#![feature(test)]

extern crate test;

use test::Bencher;

use bio::io::fasta;

use std::io::Seek;
use std::io::Write;

use rand::SeedableRng;

const NUCS: [u8; 4] = [b'A', b'C', b'T', b'G'];

fn random_seq<RNG>(length: usize, rng: &mut RNG) -> Vec<u8>
where
    RNG: rand::Rng,
{
    (0..length).map(|_| NUCS[rng.random_range(0..=3)]).collect()
}

fn write_sequence(file: &mut std::fs::File, seed: u64) {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    for i in 0..1000 {
        writeln!(
            file,
            ">{}\n{}",
            i,
            String::from_utf8(random_seq(300, &mut rng)).unwrap()
        )
        .unwrap();
    }
}

macro_rules! bench_func {
    ($name:ident, $read_expr:expr) => {
        #[bench]
        fn $name(b: &mut Bencher) {
            let mut tempfile = tempfile::tempfile().unwrap();

            write_sequence(&mut tempfile, 42);

            b.iter(|| {
                tempfile.seek(io::SeekFrom::Start(0)).unwrap();
                $read_expr(&mut tempfile).records().for_each(|_| ())
            });
        }
    };
}

mod fasta_buffer {
    use std::io;

    use super::*;

    bench_func!(default, fasta::Reader::new);
    bench_func!(wrapped_default, |tempfile| fasta::Reader::new(
        io::BufReader::new(tempfile)
    ));

    bench_func!(capacity_default, |tempfile| fasta::Reader::with_capacity(
        8192, tempfile
    ));

    bench_func!(bufread_default, |tempfile| fasta::Reader::from_bufread(
        io::BufReader::new(tempfile)
    ));

    bench_func!(wrapped_32768, |tempfile| fasta::Reader::new(
        io::BufReader::with_capacity(3278, tempfile)
    ));

    bench_func!(capacity_32768, |tempfile| fasta::Reader::with_capacity(
        32768, tempfile
    ));

    bench_func!(bufread_32768, |tempfile| fasta::Reader::from_bufread(
        io::BufReader::with_capacity(32768, tempfile)
    ));

    bench_func!(wrapped_65536, |tempfile| fasta::Reader::new(
        io::BufReader::with_capacity(65536, tempfile)
    ));

    bench_func!(capacity_65536, |tempfile| fasta::Reader::with_capacity(
        65536, tempfile
    ));

    bench_func!(bufread_65536, |tempfile| fasta::Reader::from_bufread(
        io::BufReader::with_capacity(65536, tempfile)
    ));
}
