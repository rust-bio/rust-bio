#![feature(test)]

extern crate test;

use bio::io::{fasta, fastx};
use rand::distributions::Alphanumeric;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::io;
use test::{black_box, Bencher};

const BASES: &[u8] = b"ACTG";
const ID_LEN: usize = 10;
const DESC_LEN: usize = 20;
const SEQ_LEN: usize = 100;
const FASTA_SIZE: usize = 100;

const ITERS: usize = 4000;

fn gen_random_fasta() -> io::Result<Vec<u8>> {
    let mut rng = StdRng::seed_from_u64(42);
    let mut raw_writer = Vec::new();
    {
        let mut w = fasta::Writer::new(&mut raw_writer);
        for _ in 0..FASTA_SIZE {
            let id: String = (0..ID_LEN)
                .map(|_| char::from(rng.sample(&Alphanumeric)))
                .collect();

            let desc: String = (0..DESC_LEN)
                .map(|_| char::from(rng.sample(&Alphanumeric)))
                .collect();

            let seq: Vec<u8> = (0..SEQ_LEN)
                .map(|_| {
                    let idx = rng.gen_range(0..BASES.len());
                    BASES[idx]
                })
                .collect();

            w.write_record(&fasta::Record::with_attrs(&id, Some(&desc), &seq))?;
        }
    }
    Ok(raw_writer)
}

fn fastx_count_bases<T, E, I>(records: I) -> Result<usize, E>
where
    T: fastx::Record,
    E: std::error::Error,
    I: fastx::Records<T, E>,
{
    let mut nb_bases = 0;
    for result in records {
        let record = result?;
        nb_bases += record.seq().len();
    }
    Ok(nb_bases)
}

fn fasta_count_bases<R>(records: fasta::Records<R>) -> io::Result<usize>
where
    R: io::Read,
{
    let mut nb_bases = 0;
    for result in records {
        let record = result?;
        nb_bases += record.seq().len();
    }
    Ok(nb_bases)
}

fn fastx_check<T, E, I>(records: I) -> Result<(), String>
where
    T: fastx::Record,
    E: std::error::Error,
    I: fastx::Records<T, E>,
{
    for result in records {
        let record = result.map_err(|e| format!("{}", e))?;
        record.check().map_err(|e| e.to_owned())?;
    }
    Ok(())
}

fn fasta_check<R>(records: fasta::Records<R>) -> Result<(), String>
where
    R: io::Read,
{
    for result in records {
        let record = result.map_err(|e| format!("{}", e))?;
        record.check().map_err(|e| e.to_owned())?;
    }
    Ok(())
}

#[bench]
fn bench_fasta_count(b: &mut Bencher) -> io::Result<()> {
    let data = gen_random_fasta()?;
    b.iter(|| {
        for _ in 0..ITERS {
            black_box({
                let data = io::Cursor::new(&data);
                let records = fasta::Reader::new(data).records();
                fasta_count_bases(records).unwrap();
            });
        }
    });
    Ok(())
}

#[bench]
fn bench_fastx_fasta_count(b: &mut Bencher) -> io::Result<()> {
    let data = gen_random_fasta()?;
    b.iter(|| {
        for _ in 0..ITERS {
            black_box({
                let data = io::Cursor::new(&data);
                let records = fasta::Reader::new(data).records();
                fastx_count_bases(records).unwrap();
            });
        }
    });
    Ok(())
}

#[bench]
fn bench_either_fasta_count(b: &mut Bencher) -> io::Result<()> {
    let data = gen_random_fasta()?;
    b.iter(|| {
        for _ in 0..ITERS {
            black_box({
                let data = io::Cursor::new(&data);
                let records = fastx::EitherRecords::new(data);
                fastx_count_bases(records).unwrap();
            });
        }
    });
    Ok(())
}

#[bench]
fn bench_fasta_check(b: &mut Bencher) -> io::Result<()> {
    let data = gen_random_fasta()?;
    b.iter(|| {
        for _ in 0..ITERS {
            black_box({
                let data = io::Cursor::new(&data);
                let records = fasta::Reader::new(data.clone()).records();
                fasta_check(records).unwrap();
            });
        }
    });
    Ok(())
}

#[bench]
fn bench_fastx_fasta_check(b: &mut Bencher) -> io::Result<()> {
    let data = gen_random_fasta()?;
    b.iter(|| {
        for _ in 0..ITERS {
            black_box({
                let data = io::Cursor::new(&data);
                let records = fasta::Reader::new(data.clone()).records();
                fastx_check(records).unwrap();
            });
        }
    });
    Ok(())
}

#[bench]
fn bench_either_fasta_check(b: &mut Bencher) -> io::Result<()> {
    let data = gen_random_fasta()?;
    b.iter(|| {
        for _ in 0..ITERS {
            black_box({
                let data = io::Cursor::new(&data);
                let records = fastx::EitherRecords::new(data);
                fastx_check(records).unwrap();
            });
        }
    });
    Ok(())
}
