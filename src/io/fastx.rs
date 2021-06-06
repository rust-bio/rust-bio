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
//!     I: fastx::Records<T, E> {
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
//!     Ok(fastx::Kind::FASTA) => {
//!         let mut reader = fasta::Reader::new(raw_reader);
//!         let nb_bases = count_bases(reader.records()).unwrap();
//!         println!("Number of bases: {}", nb_bases);
//!     },
//!     Ok(fastx::Kind::FASTQ) => {
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
//!     I: fastx::Records<T, E> {
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
//! use bio::io::fastx::{Record, EitherRecords};
//! use std::io;
//! use std::str;
//!
//! let mut records = EitherRecords::from(io::stdin());
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
//! use bio::io::fastx::{get_kind, get_kind_seek, get_kind_file, Kind};
//! use std::io;
//! use std::io::Read;
//! use std::fs;
//! use std::fs::File;
//! use std::io::prelude::*;
//!
//! // From a Read
//!
//! fn from_read() -> io::Result<Kind> {
//!     let reader = io::stdin();
//!     let (mut new_reader, kind) = get_kind(reader)?;
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
//! fn from_read_seek() -> io::Result<Kind> {
//!     let mut read_seek = io::Cursor::new(b">id desc
//! ACTG
//! ");
//!
//!     get_kind_seek(&mut read_seek)
//! }
//!
//! // From a file path
//!
//! fn from_file_path() -> io::Result<Kind> {
//!     get_kind_file("foo.fasta")
//! }
//! ```
//!
use anyhow::Context;
use std::convert::AsRef;
use std::fs;
use std::io::SeekFrom;
use std::io::prelude::*;
use std::io;
use std::mem;
use std::path::Path;
use thiserror::Error;

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
                EitherRecord::FASTA(f) => Record::$name(f),
                EitherRecord::FASTQ(f) => Record::$name(f),
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
    fn kind(&self) -> Kind;
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

    fn kind(&self) -> Kind {
        Kind::FASTA
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

    fn kind(&self) -> Kind {
        Kind::FASTQ
    }
}

#[derive(Display, Debug)]
pub enum EitherRecord {
    FASTA(fasta::Record),
    FASTQ(fastq::Record),
}

impl EitherRecord {
    pub fn to_fasta(self) -> fasta::Record {
        return self.into();
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

    matchthrough!(kind, Kind);
}

pub trait Records<R: Record, E>: Iterator<Item = Result<R, E>> {}

impl<T: io::Read> Records<fasta::Record, io::Error> for fasta::Records<T> {}

impl<T: io::Read> Records<fastq::Record, fastq::Error> for fastq::Records<T> {}

#[derive(Debug)]
enum EitherRecordsInner<R: io::Read> {
    FASTA(fasta::Records<R>),
    FASTQ(fastq::Records<R>),
}

#[derive(Debug)]
pub struct EitherRecords<R: io::Read> {
    records: Option<EitherRecordsInner<io::Chain<io::Cursor<[u8; 1]>, R>>>,
    reader: Option<R>,
}

impl EitherRecords<fs::File> {
    /// Read from a given file.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(path.as_ref())
            .map(EitherRecords::from)
            .with_context(|| format!("Failed to read fastq from {:#?}", path))
    }
}

impl<R: Read> EitherRecords<R> {
    pub fn new(reader: R) -> Self {
        EitherRecords::from(reader)
    }

    pub fn kind(&mut self) -> io::Result<Kind> {
        self.initialize()?;
        match self.records {
            Some(EitherRecordsInner::FASTA(_)) => Ok(Kind::FASTA),
            Some(EitherRecordsInner::FASTQ(_)) => Ok(Kind::FASTQ),
            None => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Data is empty",
            )),
        }
    }

    fn initialize(&mut self) -> io::Result<()> {
        if let Some(reader) = mem::replace(&mut self.reader, None) {
            match get_kind(reader) {
                Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => (),
                Err(err) => return Err(err),
                Ok((reader, Kind::FASTA)) => self.records = Some(
                    EitherRecordsInner::FASTA(fasta::Reader::new(reader).records())
                ),
                Ok((reader, Kind::FASTQ)) => self.records = Some(
                    EitherRecordsInner::FASTQ(fastq::Reader::new(reader).records())
                ),
            }
        }
        Ok(())
    }
}

impl<R: Read> Iterator for EitherRecords<R> {
    type Item = Result<EitherRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Err(e) = self.initialize() {
            return Some(Err(Error::IO(e)));
        }
        match &mut self.records {
            Some(EitherRecordsInner::FASTA(r)) => r
                .next()
                .map(|record_res| record_res.map(EitherRecord::from).map_err(Error::from)),
            Some(EitherRecordsInner::FASTQ(r)) => r
                .next()
                .map(|record_res| record_res.map(EitherRecord::from).map_err(Error::from)),
            None => None,
        }
    }
}

impl<R: Read> From<R> for EitherRecords<R> {
    fn from(reader: R) -> Self {
        EitherRecords {
            records: None,
            reader: Some(reader),
        }
    }
}

#[derive(Display, Debug, Error)]
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
pub enum Kind {
    FASTQ,
    FASTA,
}

/// Determine whether a [`Read`](Read) is a FastA or FastQ.
///
/// This method takes ownership of the [`Read`](Read) and returns a new [`Read`](Read)
/// with position identical to the position of the input [`Read`](Read).
/// This allows you to pass the returned [`Read`]#Read) directly into
/// [`fasta::Reader::new`](fasta::Reader::new) or [`fastq::Reader::new`](fastq::Reader::new).
///
/// In general this function is superior to [`get_kind_preserve_read`](get_kind) unless
/// you would like to do something with the [`Read`](Read) in cases where
/// there is an error determining the type.
///
/// Due to the implementation of the function it is sometimes impossible to return
///
/// # Example
///
/// ```rust
/// use bio::io::{fasta, fastq};
/// use bio::io::fastx::{Kind, get_kind};
/// use std::io;
///
/// fn count_records() -> io::Result<usize> {
///     let (reader, kind) = get_kind(io::stdin())?;
///     match kind {
///         Kind::FASTA => Ok(fasta::Reader::new(reader).records().count()),
///         Kind::FASTQ => Ok(fastq::Reader::new(reader).records().count()),
///     }
/// }
/// ```
pub fn get_kind<R: io::Read>(mut reader: R) -> io::Result<(io::Chain<io::Cursor<[u8; 1]>, R>, Kind)> {
    let mut buf = [0];
    reader.read_exact(&mut buf)?;
    let first = char::from(buf[0]);
    let new_reader = io::Cursor::new(buf).chain(reader);

    match first {
        '>' => Ok((new_reader, Kind::FASTA)),
        '@' => Ok((new_reader, Kind::FASTQ)),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Data is not a valid FASTA/FASTQ, illegal start character '{}'",
                first
            ),
        )),
    }
}

/// Determine whether a [`Read`](Read) is a FastA or FastQ.
///
/// You should only use this function if you would like to use your input [`Read`](Read)
/// for something else if there is an error determining the data type. Otherwise,
/// [`get_kind`](get_kind) is preferred.
///
/// This function is very similar to [`get_kind`](get_kind) with the
/// differences being that [`get_kind_preserve_read`](get_kind_preserve_read):
///
/// - Returns the [`Read`](Read) in the correct location if there is an
///   error in determining the data type
/// - Returns a tuple containing a [`Result`](Result) instead of a
///   [`Result`](Result), which can be less convenient to work with
/// - Requires [`Box`](Box) and `dyn`
///
/// # Example
///
/// ```rust
/// use bio::io::{fasta, fastq};
/// use bio::io::fastx::{Kind, get_kind_preserve_read};
/// use std::io;
///
/// fn print_type() -> io::Result<()> {
///     let (mut reader, kind) = get_kind_preserve_read(Box::new(io::stdin()));
///     match kind {
///         Ok(Kind::FASTA) => println!("{}", Kind::FASTA),
///         Ok(Kind::FASTQ) => println!("{}", Kind::FASTQ),
///         Err(e) => {
///             println!("Error determining FastA/FastQ: {}", e);
///             println!("Data:");
///             io::copy(&mut reader, &mut io::stdout());
///         }
///     }
///     Ok(())
/// }
/// ```
pub fn get_kind_preserve_read(
    mut reader: Box<dyn io::Read>,
) -> (Box<dyn io::Read>, io::Result<Kind>) {
    let mut buf = [0];
    if let Err(e) = reader.read_exact(&mut buf) {
        return (reader, Err(e));
    }
    let first = char::from(buf[0]);
    let new_reader = Box::new(io::Cursor::new(buf).chain(reader));

    match first {
        '>' => (new_reader, Ok(Kind::FASTA)),
        '@' => (new_reader, Ok(Kind::FASTQ)),
        _ => (
            new_reader,
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Data is not a valid FASTA/FASTQ, illegal start character '{}'",
                    first
                ),
            )),
        ),
    }
}

/// Determine whether a [`Read`](Read) + [`Seek`](Seek) is a FastA or FastQ.
///
/// The benefit of this this function compared to [`get_kind`](get_kind) is that
/// this function does not take ownership of the [`Read`](Read) so it can
/// be slightly more convenient to use.
pub fn get_kind_seek<R: io::Read + io::Seek>(reader: &mut R) -> io::Result<Kind> {
    let mut buf = [0];
    reader.read_exact(&mut buf)?;
    reader.seek(SeekFrom::Current(-1))?;
    let first = char::from(buf[0]);

    match first {
        '>' => Ok(Kind::FASTA),
        '@' => Ok(Kind::FASTQ),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Data is not a valid FASTA/FASTQ, illegal start character '{}'",
                first
            ),
        )),
    }
}

/// Determine whether a file is a FastA or FastQ.
pub fn get_kind_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> io::Result<Kind> {
    fs::File::open(&path).and_then(|mut f| get_kind_seek(&mut f))
}

impl std::fmt::Display for Kind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Kind::FASTA => "FastA",
                Kind::FASTQ => "FastQ",
            }
        )
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
    const INCOMPLETE_FASTA_FILE: &[u8] = b">id desc
>id2
";

    const FASTQ_FILE: &[u8] = b"@id desc
ACCGTAGGCTGA
+
IIIIIIJJJJJJ
";

    const INCOMPLETE_FASTQ_FILE: &[u8] = b"@id desc
ACCGTAGGCTGA
@id2 desc
CGTAAAATGA
";

    struct MockReader<R: Read> {
        error: Option<io::Error>,
        reader: R,
    }

    impl<R: Read> Read for MockReader<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            match mem::replace(&mut self.error, None) {
                Some(e) => Err(e),
                None => self.reader.read(buf),
           }
        }
    }

    #[test]
    fn records_trait() {
        fn count_records<R: Record, E, I: Records<R, E>>(records: I) -> usize {
            records.count()
        }
        let records = fasta::Reader::new(FASTA_FILE).records();
        let count = count_records(records);
        assert_eq!(count, 2);
    }

    #[test]
    fn kind_display() {
        assert_eq!(format!("{}", Kind::FASTA), "FastA");
        assert_eq!(format!("{}", Kind::FASTQ), "FastQ");
    }

    #[test]
    fn get_fasta_either_records() {
        let mut records = EitherRecords::from(FASTA_FILE);
        assert_eq!(records.next().unwrap().unwrap().id(), "id");
        assert_eq!(records.next().unwrap().unwrap().id(), "id2");
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn get_fasta_either_records_err() {
        let mut records = EitherRecords::from(INCOMPLETE_FASTA_FILE);
        assert!(matches!(records.next().unwrap().unwrap_err(), Error::IO(_)));
    }

    #[test]
    fn get_fastq_either_records() {
        let mut records = EitherRecords::from(FASTQ_FILE);
        assert_eq!(records.next().unwrap().unwrap().id(), "id");
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn get_fastq_either_records_err() {
        let mut records = EitherRecords::from(INCOMPLETE_FASTQ_FILE);
        assert!(matches!(records.next().unwrap().unwrap_err(), Error::FASTQ(_)));
    }

    #[test]
    fn get_empty_either_records() {
        let mut records = EitherRecords::from(b"".as_ref());
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn get_invalid_either_records() {
        let mut records = EitherRecords::from(b"(".as_ref());
        assert!(records.next().unwrap().is_err());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn test_get_kind_preserve_read_fasta() {
        let (mut new_read, fastx_kind) = get_kind_preserve_read(Box::new(FASTA_FILE));
        assert_eq!(Kind::FASTA, fastx_kind.unwrap());
        let mut buf = [0u8; 1];
        new_read.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTA_FILE[0]);
    }

    #[test]
    fn test_get_kind_preserve_read_fastq() {
        let (mut new_read, fastx_kind) = get_kind_preserve_read(Box::new(FASTQ_FILE));
        assert_eq!(Kind::FASTQ, fastx_kind.unwrap());
        let mut buf = [0u8; 1];
        new_read.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTQ_FILE[0]);
    }

    #[test]
    fn test_get_kind_preserve_read_empty() {
        let (_, kind_res) = get_kind_preserve_read(Box::new(Cursor::new(b"")));
        assert_eq!(kind_res.err().unwrap().kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_get_kind_preserve_read_invalid() {
        let read = Cursor::new(b"*");
        let (mut new_read, res) = get_kind_preserve_read(Box::new(read));
        assert_eq!(io::ErrorKind::InvalidData, res.err().unwrap().kind());
        let mut buf = [0u8; 1];
        new_read.read_exact(&mut buf).unwrap();
        assert_eq!(b'*', buf[0]);
    }

    #[test]
    fn test_get_kind_preserve_read_error() {
        let reader = Cursor::new(FASTA_FILE);
        let mock_reader = MockReader { reader, error: Some(io::Error::new(io::ErrorKind::Other, "error")) };
        let (mut new_reader, res) = get_kind_preserve_read(Box::new(mock_reader));
        assert!(matches!(res.unwrap_err().kind(), io::ErrorKind::Other));
        let mut buf = [0u8; 1];
        new_reader.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTA_FILE[0]);
    }

    #[test]
    fn test_get_kind_seek_fasta() {
        let mut read_seeker = Cursor::new(FASTA_FILE);
        let fastx_kind = get_kind_seek(&mut read_seeker).unwrap();
        assert_eq!(Kind::FASTA, fastx_kind);
        assert_eq!(read_seeker.position(), 0);
    }

    #[test]
    fn test_get_kind_seek_fastq() {
        let mut read_seeker = Cursor::new(FASTQ_FILE);
        let fastq_kind = get_kind_seek(&mut read_seeker).unwrap();
        assert_eq!(Kind::FASTQ, fastq_kind);
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
