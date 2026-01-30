// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! General TSV format reading and writing.
//!
//! This module is mainly used as starting point to build TSV readers for specific formats, like
//! BED.
//!
//! # Example - basic usage
//!
//! ```
//! use bio::io::common::{Reader, Record, Writer};
//!
//! let example = b"1\t5\t5000\tname1\t0.5";
//! let mut reader: Reader<&[u8], Record> = Reader::new(&example[..]);
//! let mut writer = Writer::new(vec![]);
//! for record in reader.records() {
//!     let rec = record.expect("Error reading record.");
//!     println!("{}", rec.chrom());
//!     writer.write(&rec).expect("Error writing record.");
//! }
//! ```
//!
//! # Example - reader for BED files
//!
//! ```
//! use bio::io::common;
//! use derefable::Derefable;
//! use serde::{Deserialize, Serialize};
//!
//! // Create a type alias for common::Reader which deserializes into BedRecord when calling
//! // self.records()
//! type Reader<R> = common::Reader<R, BedRecord>;
//! type Writer<W> = common::Writer<W>;
//!
//! // The derefable crate allows to inherit the methods of common::Record.
//! #[derive(
//!     Debug,
//!     Default,
//!     Clone,
//!     Eq,
//!     PartialEq,
//!     Ord,
//!     PartialOrd,
//!     Hash,
//!     Serialize,
//!     Deserialize,
//!     Derefable,
//! )]
//! struct BedRecord(#[deref(mutable)] common::Record);
//!
//! impl BedRecord {
//!     fn new() -> Self {
//!         Self::default()
//!     }
//!     fn name(&self) -> Option<&str> {
//!         self.aux(3)
//!     }
//! }
//!
//! let example = b"1\t5\t5000\tname1";
//! let mut reader = Reader::new(&example[..]);
//! let mut writer = Writer::new(vec![]);
//! for record in reader.records() {
//!     let rec = record.expect("Error reading bed record.");
//!     println!("{}", rec.name().expect("Error reading name"));
//!     writer.write(&rec).expect("Error writing bed record.");
//! }
//! ```
use std::convert::AsRef;
use std::fs;
use std::io;
use std::marker::PhantomData;
use std::path::Path;

use anyhow::Context;
use serde::de::DeserializeOwned;
use serde::Deserialize;
use serde::Serialize;

use getset::{CloneGetters, CopyGetters, Getters, MutGetters, Setters, WithSetters};

/// A TSV reader.
#[derive(Debug)]
pub struct Reader<R: io::Read, T: DeserializeOwned> {
    inner: csv::Reader<R>,
    phantom: PhantomData<T>,
}

impl<T: DeserializeOwned> Reader<fs::File, T> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to read from {:#?}", path))
    }
}

impl<R: io::Read, T: DeserializeOwned> Reader<R, T> {
    /// Read from a given reader.
    pub fn new(reader: R) -> Self {
        Reader {
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_reader(reader),
            phantom: PhantomData,
        }
    }

    /// Iterate over all records.
    pub fn records(&mut self) -> Records<'_, R, T> {
        Records {
            inner: self.inner.deserialize(),
        }
    }
}

/// An iterator over the records of a file.
pub struct Records<'a, R: io::Read, T: DeserializeOwned> {
    inner: csv::DeserializeRecordsIter<'a, R, T>,
}

impl<'a, R: io::Read, T: DeserializeOwned> Iterator for Records<'a, R, T> {
    type Item = csv::Result<T>;

    fn next(&mut self) -> Option<csv::Result<T>> {
        self.inner.next()
    }
}

/// The core writer usable for any record.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
}

impl Writer<fs::File> {
    /// Write to a given file path.
    #[allow(clippy::wrong_self_convention)]
    pub fn to_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::create(path).map(Writer::new)
    }
}

impl<W: io::Write> Writer<W> {
    /// Write to a given writer.
    pub fn new(writer: W) -> Self {
        Writer {
            inner: csv::WriterBuilder::new()
                .delimiter(b'\t')
                .flexible(true)
                .has_headers(false)
                .from_writer(writer),
        }
    }

    /// Write a given record.
    pub fn write(&mut self, record: impl Serialize) -> csv::Result<()> {
        self.inner.serialize(record)
    }
}

/// A base BED record
///
/// It can be used as starting point to construct records for other file types like .narrowPeak,
/// .broadPeak, ecc..
///
/// With the the derefable crate you can write a wrapper around common::Record that maintains its methods,
/// while giving you the possibility of expanding them.
///
/// ```
/// use bio::io::common;
/// use derefable::Derefable;
/// use serde::{Deserialize, Serialize};
///
/// type Reader<R> = common::Reader<R, BedRecord>;
///
/// #[derive(
///     Debug,
///     Default,
///     Clone,
///     Eq,
///     PartialEq,
///     Ord,
///     PartialOrd,
///     Hash,
///     Serialize,
///     Deserialize,
///     Derefable,
/// )]
/// pub struct BedRecord(#[deref(mutable)] common::Record);
///
/// /// Add new methods specific for a BED Record
/// impl BedRecord {
///     /// Create a new Record
///     pub fn new() -> Self {
///         Self::default()
///     }
///
///     /// Name of the feature.
///     pub fn name(&self) -> Option<&str> {
///         self.aux(3)
///     }
///
///     /// Set name.
///     pub fn set_name(&mut self, name: &str) {
///         self.set_aux(3, name)
///     }
/// }
///
/// fn main() {
///     let src = b"chr1\t2\t3\tname1";
///     let mut reader = Reader::new(&src[..]);
///     for record in reader.records() {
///         let rec: BedRecord = record.unwrap();
///         assert_eq!(rec.chrom(), "chr1");
///         assert_eq!(rec.start(), 2);
///         assert_eq!(rec.end(), 3);
///         assert_eq!(rec.name(), Some("name1"));
///     }
/// }
/// ```
#[derive(
    Debug,
    Default,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Serialize,
    Deserialize,
    Getters,
    Setters,
    WithSetters,
    MutGetters,
    CopyGetters,
    CloneGetters,
)]
pub struct Record {
    // Chromosome of feature.
    #[getset(get_mut = "pub", set_with = "pub")]
    chrom: String,

    // Start position of feature (0-based).
    #[getset(get_copy = "pub", set = "pub", get_mut = "pub", set_with = "pub")]
    start: u64,

    // End position of feature (0-based, not included).
    #[getset(get_copy = "pub", set = "pub", get_mut = "pub", set_with = "pub")]
    end: u64,

    #[serde(default)]
    aux: Vec<String>,
}

impl Record {
    /// Create new record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Chromosome of the feature.
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Access auxilliary fields after the end field by index
    /// (counting first field (chromosome) as 0).
    pub fn aux(&self, i: usize) -> Option<&str> {
        let j = i - 3;
        if j < self.aux.len() {
            Some(&self.aux[j])
        } else {
            None
        }
    }
}

impl Record {
    /// Set chromosome.
    pub fn set_chrom(&mut self, chrom: &str) {
        self.chrom = chrom.into()
    }

    /// Set auxiliary field after the end field by index.
    /// If the index is bigger than the current auxiliary vector, empty strings will be pushed into
    /// the vector until the index is reached.
    /// (counting first field (chromosome) as 0).
    pub fn set_aux(&mut self, i: usize, field: &str) {
        let j = i - 3;
        if j < self.aux.len() {
            self.aux[j] = field.into();
        } else {
            for _ in self.aux.len()..j {
                self.aux.push("".to_owned());
            }
            self.aux.push(field.into());
        }
    }

    /// Add auxilliary field.
    pub fn push_aux(&mut self, field: &str) {
        self.aux.push(field.into());
    }
}

#[cfg(test)]
mod tests {
    use derefable::Derefable;

    use super::*;

    const BED_FILE: &[u8] = b"1\t5\t5000\tname1
2\t3\t5005\tname2
";
    const BED_FILE_COMMENT: &[u8] = b"\
# this line should be ignored
1\t5\t5000\tname1
# and this one as well
2\t3\t5005\tname2
";
    const BED_FILE_COMPACT: &[u8] = b"1\t5\t5000\n2\t3\t5005\n";

    type TestReader<R> = Reader<R, TestRecord>;

    // Create test record
    #[derive(
        Debug,
        Default,
        Clone,
        Eq,
        PartialEq,
        Ord,
        PartialOrd,
        Hash,
        Serialize,
        Deserialize,
        Derefable,
    )]
    struct TestRecord(#[deref(mutable)] Record);

    // Implement new methods for TestRecord that were not available for the common::Record
    impl TestRecord {
        pub fn name(&self) -> Option<&str> {
            self.aux(3)
        }

        pub fn set_name(&mut self, name: &str) {
            self.set_aux(3, name);
        }
    }

    #[test]
    fn test_core_reader() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];

        let mut reader = Reader::new(BED_FILE);
        for (i, r) in reader.records().enumerate() {
            let record: Record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.aux(3).expect("Error reading name"), names[i]);
        }
    }

    #[test]
    fn test_core_setters() {
        let mut rec = Record::default();
        rec.set_chrom("chr1");
        rec.set_start(1);
        rec.set_end(2);
        rec.set_aux(4, "fourth");

        assert_eq!(rec.chrom(), "chr1");
        assert_eq!(rec.start(), 1);
        assert_eq!(rec.end(), 2);
        assert_eq!(rec.aux(3), Some(""));
        assert_eq!(rec.aux(4), Some("fourth"));
        assert_eq!(rec.aux(7), None);
    }

    #[test]
    fn test_core_reader_with_comment() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];

        let mut reader = Reader::new(BED_FILE_COMMENT);
        for (i, r) in reader.records().enumerate() {
            let record: Record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
        }
    }

    #[test]
    fn test_core_reader_compact() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];

        let mut reader = Reader::new(BED_FILE_COMPACT);
        for (i, r) in reader.records().enumerate() {
            let record: Record = r.unwrap();
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
        }
    }

    #[test]
    fn test_core_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.bed");
        let error = Reader::<std::fs::File, Record>::from_file(path)
            .unwrap_err()
            .downcast::<String>()
            .unwrap();

        assert_eq!(&error, "Failed to read from \"/I/dont/exist.bed\"")
    }

    #[test]
    fn test_core_writer() {
        let mut reader: Reader<&[u8], Record> = Reader::new(BED_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records() {
            writer
                .write(r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), BED_FILE);
    }

    #[test]
    fn test_core_reader_with_implemented_record() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];

        let mut reader: Reader<&[u8], TestRecord> = Reader::new(BED_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
        }
    }

    #[test]
    fn test_implemented_reader() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];

        let mut reader = TestReader::new(BED_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
        }
    }

    #[test]
    fn test_implemented_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.bed");
        let error = TestReader::from_file(path)
            .unwrap_err()
            .downcast::<String>()
            .unwrap();

        assert_eq!(&error, "Failed to read from \"/I/dont/exist.bed\"")
    }

    #[test]
    fn test_implemented_record_setters() {
        let mut rec = TestRecord::default();
        rec.set_chrom("chr1");
        rec.set_start(1);
        rec.set_end(2);
        rec.set_name("name1");
        rec.set_aux(4, "fourth");
        rec.set_aux(6, "sixth");

        assert_eq!(rec.chrom(), "chr1");
        assert_eq!(rec.start(), 1);
        assert_eq!(rec.end(), 2);
        assert_eq!(rec.name(), Some("name1"));
        assert_eq!(rec.aux(3), Some("name1"));
        assert_eq!(rec.aux(4), Some("fourth"));
        assert_eq!(rec.aux(5), Some(""));
        assert_eq!(rec.aux(6), Some("sixth"));
        assert_eq!(rec.aux(7), None);
    }

    #[test]
    fn test_writer_with_implemented_reader() {
        let mut reader = TestReader::new(BED_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records() {
            writer
                .write(r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), BED_FILE);
    }
}
