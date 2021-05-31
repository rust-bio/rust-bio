// Copyright 2021 Todd Morse.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Traits and utilities to read and write files in the FASTA and FASTQ format interchangably
//!
//! Files in the FASTA and FASTQ format have several common fields: ID, Sequence, and Description.
//! These utilities can be used to implement algorithms that only require these fields for files
//! that are in either format.
//!
//! # Example
//!
//! ## Common Statistics
//!
//! In this example, we implement a count_bases function that for supports both the FASTA and FASTQ
//! format.
//!
//! ```
//! use bio::io::{fasta, fastq, fastx};
//! use std::io;
//!
//! fn count_bases<T, E, I>(mut records: I) -> Result<usize, E>
//! where T: fastx::Record,
//!     E: std::error::Error,
//!     I: Iterator<Item = Result<T, E>> {
//!     let mut nb_bases = 0;
//!     for result in records {
//!         let record = result?;
//!         nb_bases += record.seq().len();
//!     }
//!     Ok(nb_bases)
//! }
//!
//! let mut raw_reader = io::Cursor::new(b">id desc
//! ACTG
//! ");
//!
//! match fastx::get_type(&mut raw_reader) {
//!     Ok(fastx::FastxType::Fasta) => {
//!         let mut reader = fasta::Reader::new(raw_reader);
//!         let nb_bases = count_bases(reader.records()).unwrap();
//!         println!("Number of bases: {}", nb_bases);
//!     },
//!     Ok(fastx::FastxType::Fastq) => {
//!         let mut reader = fastq::Reader::new(raw_reader);
//!         let nb_bases = count_bases(reader.records()).unwrap();
//!         println!("Number of bases: {}", nb_bases);
//!     },
//!     _ => println!("Encountered an error"),
//! }
//! ```
//!
//! ## Filtration
//!
//! In this example, we define an at_least_n_bases function that can filter fasta or fastq
//! records based on their sequence lengths. It works seemlessly with fasta records as if
//! it was implemented just on fasta records. In a realistic scenario this function might be
//! defined in a library so callers could use it with both fasta and fastq files as needed.
//!
//! ```
//! use bio::io::{fasta, fastq, fastx};
//! use std::io;
//!
//! fn at_least_n_bases<T, E, I>(mut records: I, n: usize) -> impl Iterator<Item = Result<T, E>>
//! where T: fastx::Record,
//!     E: std::error::Error,
//!     I: Iterator<Item = Result<T, E>> {
//!     records.filter(move |rr| match rr {
//!         Ok(r) => r.seq().len() > n,
//!         _ => true,
//!     })
//! }
//!
//! let mut reader = fasta::Reader::new(io::stdin());
//! let mut writer = fasta::Writer::new(io::stdout());
//!
//! for record in at_least_n_bases(reader.records(), 10) {
//!     writer.write_record(&record.unwrap());
//! }
//! ```
//!
use std::io;
use std::io::SeekFrom;
use std::io::prelude::*;
use std::mem;

use crate::io::{fasta, fastq};
use crate::utils::TextSlice;

macro_rules! passthrough {
    ($name:ident, $t:ty) => {
        fn $name(&self) -> $t {
            self.$name()
        }
    };
}

macro_rules! matchthrough {
    ($name:ident, $t:ty) => {
        fn $name(&self) -> $t {
            match self {
                EitherRecord::FASTA(f) => f.$name(),
                EitherRecord::FASTQ(f) => f.$name(),
            }
        }
    };
}

pub trait Record {
    fn is_empty(&self) -> bool;
    fn check(&self) -> Result<(), &str>;
    fn id(&self) -> &str;
    fn desc(&self) -> Option<&str>;
    fn seq(&self) -> TextSlice<'_>;
}

impl Record for super::fasta::Record {
    passthrough!(is_empty, bool);
    passthrough!(check, Result<(), &str>);
    passthrough!(id, &str);
    passthrough!(desc, Option<&str>);
    passthrough!(seq, TextSlice<'_>);
}

impl Record for super::fastq::Record {
    passthrough!(is_empty, bool);
    passthrough!(check, Result<(), &str>);
    passthrough!(id, &str);
    passthrough!(desc, Option<&str>);
    passthrough!(seq, TextSlice<'_>);
}

pub enum EitherRecord {
    FASTA(fasta::Record),
    FASTQ(fastq::Record),
}

impl From<fasta::Record> for EitherRecord {
    fn from(record: fasta::Record) -> Self {
       EitherRecord::FASTA(record) 
    }
}

impl From<fastq::Record> for EitherRecord {
    fn from(record: fastq::Record) -> Self {
       EitherRecord::FASTQ(record) 
    }
}

impl Record for EitherRecord {
    matchthrough!(is_empty, bool);
    matchthrough!(check, Result<(), &str>);
    matchthrough!(id, &str);
    matchthrough!(desc, Option<&str>);
    matchthrough!(seq, TextSlice<'_>);
}

enum EitherRecords<R: io::Read> {
    FASTA(fasta::Records<R>),
    FASTQ(fastq::Records<R>),
}

pub struct Records<R: io::Read> {
    // This isn't as type safe as I'd like it to be
    records: Option<EitherRecords<R>>,
    reader: Option<R>,
}

impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<EitherRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(reader) = mem::replace(&mut self.reader, None) {
            let mut reader: io::BufReader<R> = io::BufReader::new(reader);
            let mut line = String::new();
            if let Err(e) = reader.read_line(&mut line) {
                return Some(Err(e));
            }
            match line.chars().next() {
                Some('>') => self.records = Some(EitherRecords::FASTA(fasta::Reader::new_with_line(reader, line).records())),
                Some('@') => self.records = Some(EitherRecords::FASTQ(fastq::Reader::new_with_line_buffer(reader, line).records())),
                Some(c) => return Some(Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("File is not a valid FASTA/FASTQ, illegal start character '{}'", c),
                ))),
                None => return None,
            }
        }
        match &mut self.records {
            Some(EitherRecords::FASTA(r)) => r
                .next()
                .map(|record_res| record_res.map(EitherRecord::from)),
            Some(EitherRecords::FASTQ(r)) => r
                .next()
                .map(|record_res| record_res
                     .map(EitherRecord::from)
                     .map_err(|e| io::Error::new(io::ErrorKind::Other, e)),
                 ),
             None => None,
        }
    }
}

impl<R: io::Read> From<R> for Records<R> {
    fn from(reader: R) -> Self {
        Records { records: None, reader: Some(reader) }
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum FastxType {
    Fastq,
    Fasta,
}

pub fn get_type<R: io::Read + io::Seek>(reader: &mut R) -> Result<FastxType, io::Error> {
    let mut first = [0];
    reader.read_exact(&mut first)?;
    reader.seek(SeekFrom::Current(-1))?;

    match first[0] as char {
        '>' => Ok(FastxType::Fasta),
        '@' => Ok(FastxType::Fastq),
        _ => {
            let msg = format!("first character was '{}' which is not a valid first character for a fastq or fasta file", first[0]);
            Err(io::Error::new(io::ErrorKind::InvalidData, msg))
        }
    }
}

impl std::fmt::Display for FastxType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            FastxType::Fasta => "fasta",
            FastxType::Fastq => "fastq",
        };
        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Cursor;

    const FASTA_FILE: &[u8] = b">id desc
ACCGTAGGCTGA
CCGTAGGCTGAA
CGTAGGCTGAAA
GTAGGCTGAAAA
CCCC
>id2
ATTGTTGTTTTA
ATTGTTGTTTTA
ATTGTTGTTTTA
GGGG
";

    const FASTQ_FILE: &[u8] = b"@id desc
ACCGTAGGCTGA
+
IIIIIIJJJJJJ
";
    #[test]
    fn get_fasta_records() {
        let reader = Cursor::new(FASTA_FILE);
        let mut records = Records::from(reader);
        assert_eq!(records.next().unwrap().unwrap().id(), "id");
        assert_eq!(records.next().unwrap().unwrap().id(), "id2");
        assert!(records.next().is_none());
        assert!(records.next().is_none());
    }

    #[test]
    fn get_empty_records() {
        let reader = Cursor::new("");
        let mut records = Records::from(reader);
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn get_invalid_records() {
        let reader = Cursor::new("(");
        let mut records = Records::from(reader);
        assert!(records.next().unwrap().is_err());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn test_get_type_fasta() {
        let mut read_seeker = Cursor::new(FASTA_FILE);
        let fastx_type = get_type(&mut read_seeker).unwrap();
        assert_eq!(FastxType::Fasta, fastx_type);
        assert_eq!(read_seeker.position(), 0);
    }

    #[test]
    fn test_get_type_fastq() {
        let mut read_seeker = Cursor::new(FASTQ_FILE);
        let fastq_type = get_type(&mut read_seeker).unwrap();
        assert_eq!(FastxType::Fastq, fastq_type);
        assert_eq!(read_seeker.position(), 0);
    }

    #[test]
    fn test_get_type_empty() {
        let mut read_seeker = Cursor::new(b"");
        let e = get_type(&mut read_seeker).unwrap_err();
        assert_eq!(io::ErrorKind::UnexpectedEof, e.kind());
    }

    #[test]
    fn test_get_type_invalid() {
        let mut read_seeker = Cursor::new(b"*");
        let e = get_type(&mut read_seeker).unwrap_err();
        assert_eq!(io::ErrorKind::InvalidData, e.kind());
    }
}
