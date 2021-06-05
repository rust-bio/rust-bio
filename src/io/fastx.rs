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
//! This module serves two use cases:
//!
//! 1. Implementing functions that can be used generically with FASTA/FASTQ records. In this use
//! case the type may be known at compile time by the caller of your function.
//! 2. Processing data that may be either in the FASTA/FASTQ format. In this use case the type
//! cannot be known at compile time and you may or may not want to treat FASTA/FASTQ data
//! differently.
//!
//! # Generic Implementation Examples
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
//! match fastx::get_kind_seek(&mut raw_reader) {
//!     Ok(fastx::FastxKind::FASTA) => {
//!         let mut reader = fasta::Reader::new(raw_reader);
//!         let nb_bases = count_bases(reader.records()).unwrap();
//!         println!("Number of bases: {}", nb_bases);
//!     },
//!     Ok(fastx::FastxKind::FASTQ) => {
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
//! In this example, we define an at_least_n_bases function that can filter FASTA or FASTQ
//! records based on their sequence lengths. It works seemlessly with `fasta::Record`s as if
//! it was implemented just for them. In a realistic scenario this function might be
//! defined in a library so callers could use it with both FASTA and FASTQ files as needed.
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
//! # Unknown Type Examples
//!
//! If the type of a record is not known at compile time the record is a `fastx::EitherRecord`.
//! This type is an enum containing either a `fasta::Record` or a `fastq::Record`. There are also
//! utility functions defined on the enum so you can work with it without detecting the type.
//!
//! ## Parsing data of unknown type
//!
//! ```
//! use bio::io::fastx::{Record, Records};
//! use std::io;
//! use std::str;
//!
//! let mut records = Records::from(io::stdin());
//! while let Some(Ok(record)) = records.next() {
//!     println!("id  : {}", record.id());
//!     println!("desc: {}", record.desc().unwrap_or("none"));
//!     println!("seq : {}", str::from_utf8(record.seq()).unwrap_or("error"));
//!     // Add a default quality in case we have a FASTA record
//!     let default_qual = vec![b'I' ; record.seq().len()];
//!     let qual = record.qual().unwrap_or(&default_qual);
//!     println!("qual: {}", str::from_utf8(qual).unwrap_or("error"));
//!     println!("")
//! }
//! ```
//!
//! # Utility Examples
//!
//! ## Type Detection
//!
//! This module provides utilities for detecting if some data is a FASTA/FASTQ.
//!
//! ```
//! use bio::io::fastx::{get_kind, get_kind_seek, get_kind_file, FastxKind};
//! use std::io;
//! use std::io::Read;
//! use std::fs;
//! use std::fs::File;
//! use std::io::prelude::*;
//!
//! // From a Read
//! 
//! fn from_read() -> io::Result<FastxKind> {
//!     let reader = io::stdin();
//!     let (mut new_reader, kind_res) = get_kind(Box::new(reader));
//!     let kind = kind_res?;
//!     println!("{}", kind);
//!     // Read from start of your old reader
//!     let mut buf = [0u8; 8];
//!     new_reader.read(&mut buf)?;
//!     Ok(kind)
//! }
//!
//!
//! // From a Read + Seek
//!
//! fn from_read_seek() -> io::Result<FastxKind> {
//!     let mut read_seek = io::Cursor::new(b">id desc
//! ACTG
//! ");
//!
//!     get_kind_seek(&mut read_seek)
//! }
//!
//! // From a file path
//!
//! fn from_file_path() -> io::Result<FastxKind> {
//!     get_kind_file("foo.fasta")
//! }
//! ```
//!
use std::convert::AsRef;
use anyhow::Context;
use std::fs;
use std::io::SeekFrom;
use std::io::prelude::*;
use std::io;
use std::mem;
use std::path::Path;

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
    fn qual(&self) -> Option<&[u8]>;
    fn kind(&self) -> FastxKind;
}

impl Record for super::fasta::Record {
    passthrough!(is_empty, bool);
    passthrough!(check, Result<(), &str>);
    passthrough!(id, &str);
    passthrough!(desc, Option<&str>);
    passthrough!(seq, TextSlice<'_>);

    fn qual(&self) -> Option<&[u8]> {
        None
    }

    fn kind(&self) -> FastxKind {
        FastxKind::FASTA
    }

}

impl Record for super::fastq::Record {
    passthrough!(is_empty, bool);
    passthrough!(check, Result<(), &str>);
    passthrough!(id, &str);
    passthrough!(desc, Option<&str>);
    passthrough!(seq, TextSlice<'_>);

    fn qual(&self) -> Option<&[u8]> {
        Some(self.qual())
    }

    fn kind(&self) -> FastxKind {
        FastxKind::FASTQ
    }
}

pub enum EitherRecord {
    FASTA(fasta::Record),
    FASTQ(fastq::Record),
}

impl EitherRecord {
    pub fn to_fasta(self) -> fasta::Record {
        return self.into()
    }

    pub fn to_fastq(self, default_qual: u8) -> fastq::Record {
        match self {
            EitherRecord::FASTQ(f) => f,
            EitherRecord::FASTA(f) => {
                let qual = &vec![default_qual; f.seq().len()];
                fastq::Record::with_attrs(f.id(), f.desc(), f.seq(), qual)
            }
        }
    }
}

impl Into<fasta::Record> for EitherRecord {
    fn into(self) -> fasta::Record {
        match self {
            EitherRecord::FASTA(f) => f,
            EitherRecord::FASTQ(f) => fasta::Record::with_attrs(f.id(), f.desc(), f.seq()),
        }
    }
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

    fn qual(&self) -> Option<&[u8]> {
        match &self {
            EitherRecord::FASTA(f) => Record::qual(f),
            EitherRecord::FASTQ(f) => Record::qual(f),
        }
    }

    matchthrough!(kind, FastxKind);
}

#[derive(Debug)]
enum EitherRecords<R: io::Read> {
    FASTA(fasta::Records<R>),
    FASTQ(fastq::Records<R>),
}

#[derive(Debug)]
pub struct Records<R: io::Read> {
    records: Option<EitherRecords<R>>,
    reader: Option<R>,
}

impl Records<fs::File> {
    /// Read from a given file.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(path.as_ref())
            .map(Records::from)
            .with_context(|| format!("Failed to read fastq from {:#?}", path))
    }
}

impl <R: io::Read> Records<R> {
    pub fn new(reader: R) -> Self {
        Records::from(reader)
    }

    pub fn kind(&mut self) -> io::Result<FastxKind> {
        self.initialize()?;
        match self.records {
            Some(EitherRecords::FASTA(_)) => Ok(FastxKind::FASTA),
            Some(EitherRecords::FASTQ(_)) => Ok(FastxKind::FASTQ),
            None => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Data is empty",
            )),
        }
    }

    fn initialize(&mut self) -> io::Result<()> {
        if let Some(reader) = mem::replace(&mut self.reader, None) {
            let mut reader: io::BufReader<R> = io::BufReader::new(reader);
            let mut line = String::new();
            reader.read_line(&mut line)?;
            match line.chars().next() {
                Some('>') => self.records = Some(EitherRecords::FASTA(fasta::Reader::new_with_line(reader, line).records())),
                Some('@') => self.records = Some(EitherRecords::FASTQ(fastq::Reader::new_with_line_buffer(reader, line).records())),
                Some(c) => return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Data is not a valid FASTA/FASTQ, illegal start character '{}'", c),
                )),
                None => (),
            }
        }
        Ok(())
    }
}

impl<R: io::Read> Iterator for Records<R> {
    type Item = Result<EitherRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Err(e) = self.initialize() {
            return Some(Err(Error::IO(e)));
        }
        match &mut self.records {
            Some(EitherRecords::FASTA(r)) => r
                .next()
                .map(|record_res| record_res
                     .map(EitherRecord::from)
                     .map_err(Error::from)
                ),
            Some(EitherRecords::FASTQ(r)) => r
                .next()
                .map(|record_res| record_res
                     .map(EitherRecord::from)
                     .map_err(Error::from)
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

#[derive(Debug)]
pub enum Error {
    IO(io::Error),
    FASTQ(fastq::Error),
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Self {
       Error::IO(err) 
    }
}

impl From<fastq::Error> for Error {
    fn from(err: fastq::Error) -> Self {
       Error::FASTQ(err) 
    }
}

type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Debug, Eq, PartialEq)]
pub enum FastxKind {
    FASTQ,
    FASTA,
}

/// Determine whether a [`Read`](#Read) is a FastA or FastQ.
///
/// This method takes ownership of the [`Read`](#Read) and returns a new [`Read`](#Read)
/// with position identical to the position of the input [`Read`](#Read).
/// This allows you to pass the returned [`Read`](#Read) directly into
/// [`fasta::Reader::new`](#super::fasta::Reader::new) or [`fastq::Reader::new`](#super::fastq::Reader::new).
///
/// Due to the implementation of the function it is sometimes impossible to return 
///
/// # Example
///
/// ```rust
/// use bio::io::{fasta, fastq};
/// use bio::io::fastx::{FastxKind, get_kind};
/// use std::io;
///
/// fn count_records() -> io::Result<usize> {
///     let (reader, kind) = get_kind(Box::new(io::stdin()));
///     match kind? {
///         FastxKind::FASTA => Ok(fasta::Reader::new(reader).records().count()),
///         FastxKind::FASTQ => Ok(fastq::Reader::new(reader).records().count()),
///     }
/// }
/// ```
pub fn get_kind(mut reader: Box<dyn io::Read>) -> (Box<dyn io::Read>, io::Result<FastxKind>) {
    let mut buf = [0];
    if let Err(e) = reader.read_exact(&mut buf) {
        return (reader, Err(e));
    }
    let first = char::from(buf[0]);
    let new_reader = Box::new(io::Cursor::new(buf).chain(reader));

    match first {
        '>' => (new_reader, Ok(FastxKind::FASTA)),
        '@' => (new_reader, Ok(FastxKind::FASTQ)),
        _ => (new_reader, Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Data is not a valid FASTA/FASTQ, illegal start character '{}'", first),
        ))),
    }
}

/// Determine whether a [`Read`](#Read) + [`Seek`](#Seek) is a FastA or FastQ.
///
/// The benefit of this this function compared to [`get_kind`](#get_kind) is that
/// this function does not take ownership of the [`Read`](#Read) so it can
/// be slightly more convenient to use.
pub fn get_kind_seek<R: io::Read + io::Seek>(reader: &mut R) -> io::Result<FastxKind> {
    let mut buf = [0];
    reader.read_exact(&mut buf)?;
    reader.seek(SeekFrom::Current(-1))?;
    let first = char::from(buf[0]);

    match first {
        '>' => Ok(FastxKind::FASTA),
        '@' => Ok(FastxKind::FASTQ),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Data is not a valid FASTA/FASTQ, illegal start character '{}'", first),
        )),
    }
}

/// Determine whether a file is a FastA or FastQ.
pub fn get_kind_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> io::Result<FastxKind> {
    fs::File::open(&path).and_then(|mut f| get_kind_seek(&mut f))
}

impl std::fmt::Display for FastxKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", match self {
            FastxKind::FASTA => "fasta",
            FastxKind::FASTQ => "fastq",
        })
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
        let mut records = Records::from(FASTA_FILE);
        assert_eq!(records.next().unwrap().unwrap().id(), "id");
        assert_eq!(records.next().unwrap().unwrap().id(), "id2");
        assert!(records.next().is_none());
        assert!(records.next().is_none());
    }

    #[test]
    fn get_empty_records() {
        let mut records = Records::from(b"".as_ref());
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn get_invalid_records() {
        let mut records = Records::from(b"(".as_ref());
        assert!(records.next().unwrap().is_err());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn test_get_kind_fasta() {
        let (mut new_read, fastx_kind) = get_kind(Box::new(FASTA_FILE));
        assert_eq!(FastxKind::FASTA, fastx_kind.unwrap());
        let mut buf = [0u8; 1];
        new_read.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTA_FILE[0]);

    }

    #[test]
    fn test_get_kind_fastq() {
        let (mut new_read, fastx_kind) = get_kind(Box::new(FASTQ_FILE));
        assert_eq!(FastxKind::FASTQ, fastx_kind.unwrap());
        let mut buf = [0u8; 1];
        new_read.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTQ_FILE[0]);
    }

    #[test]
    fn test_get_kind_empty() {
        let (_, kind_res) = get_kind(Box::new(Cursor::new(b"")));
        assert_eq!(kind_res.err().unwrap().kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_get_kind_invalid() {
        let mut read_seeker = Cursor::new(b"*");
        let e = get_kind_seek(&mut read_seeker).unwrap_err();
        assert_eq!(io::ErrorKind::InvalidData, e.kind());
    }

    #[test]
    fn test_get_kind_seek_fasta() {
        let mut read_seeker = Cursor::new(FASTA_FILE);
        let fastx_kind = get_kind_seek(&mut read_seeker).unwrap();
        assert_eq!(FastxKind::FASTA, fastx_kind);
        assert_eq!(read_seeker.position(), 0);
    }

    #[test]
    fn test_get_kind_seek_fastq() {
        let mut read_seeker = Cursor::new(FASTQ_FILE);
        let fastq_kind = get_kind_seek(&mut read_seeker).unwrap();
        assert_eq!(FastxKind::FASTQ, fastq_kind);
        assert_eq!(read_seeker.position(), 0);
    }

    #[test]
    fn test_get_kind_seek_empty() {
        let mut read_seeker = Cursor::new(b"");
        let e = get_kind_seek(&mut read_seeker).unwrap_err();
        assert_eq!(io::ErrorKind::UnexpectedEof, e.kind());
    }

    #[test]
    fn test_get_kind_seek_invalid() {
        let mut read_seeker = Cursor::new(b"*");
        let e = get_kind_seek(&mut read_seeker).unwrap_err();
        assert_eq!(io::ErrorKind::InvalidData, e.kind());
    }
}
