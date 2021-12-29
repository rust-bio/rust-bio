// Copyright 2014-2018 Johannes Köster, Henning Timm.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Structs and trait to read and write files in FASTQ format.
//!
//! # Example
//!
//! ## Read
//!
//! In this example, we parse a fastq file from stdin and compute some statistics
//!
//! ```
//! use bio::io::fastq;
//! use std::io;
//! let mut reader = fastq::Reader::new(io::stdin());
//!
//! let mut nb_reads = 0;
//! let mut nb_bases = 0;
//!
//! for result in reader.records() {
//!     let record = result.expect("Error during fastq record parsing");
//!
//!     nb_reads += 1;
//!     nb_bases += record.seq().len();
//! }
//!
//! println!("Number of reads: {}", nb_reads);
//! println!("Number of bases: {}", nb_bases);
//! ```
//!
//! We can also use a `while` loop to iterate over records
//! ```
//! use bio::io::fastq;
//! use std::io;
//! let mut records = fastq::Reader::new(io::stdin()).records();
//!
//! let mut nb_reads = 0;
//! let mut nb_bases = 0;
//!
//! while let Some(Ok(record)) = records.next() {
//!     nb_reads += 1;
//!     nb_bases += record.seq().len();
//! }
//!
//! println!("Number of reads: {}", nb_reads);
//! println!("Number of bases: {}", nb_bases);
//! ```
//!
//! ## Write
//!
//! In this example we generate 10 random sequences with length 100 and write them to stdout.
//!
//! ```
//! use std::io;
//! use bio::io::fastq;
//!
//! let mut seed = 42;
//!
//! let nucleotides = [b'A', b'C', b'G', b'T'];
//!
//! let mut writer = fastq::Writer::new(io::stdout());
//!
//! for _ in 0..10 {
//!     let seq = (0..100).map(|_| {
//!         seed = ((seed ^ seed << 13) ^ seed >> 7) ^ seed << 17; // don't use this random generator
//!         nucleotides[seed % 4]
//!     }).collect::<Vec<u8>>();
//!
//!     let qual = (0..100).map(|_| b'!').collect::<Vec<u8>>();
//!
//!    writer.write("random", None, seq.as_slice(), qual.as_slice());
//! }
//! ```
//!
//! ## Read and Write
//!
//! In this example we filter reads from stdin on mean quality (Phred + 33) and write them to stdout
//!
//! ```
//! use bio::io::fastq;
//! use bio::io::fastq::FastqRead;
//! use std::io;
//!
//! let mut reader = fastq::Reader::new(io::stdin());
//! let mut writer = fastq::Writer::new(io::stdout());
//! let mut record = fastq::Record::new();
//!
//! while let Ok(()) = reader.read(&mut record) {
//!     if record.is_empty() {
//!         let check = record.check();
//!         break;
//!     }
//!
//!     let mut sum_qual = record.qual().iter().sum::<u8>() as f64;
//!
//!     if (sum_qual / record.seq().len() as f64 - 33.0) > 30.0 {
//!         writer.write_record(&record);
//!     }
//! }
//! ```

use anyhow::Context;
use std::convert::AsRef;
use std::fmt;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("expected '@' at record start")]
    MissingAt,

    #[error("can't open {path} file: {source}")]
    FileOpen { path: PathBuf, source: io::Error },

    #[error("can't read input")]
    ReadError(#[from] io::Error),

    #[error("Incomplete record. Each FastQ record has to consist of 4 lines: header, sequence, separator and qualities.")]
    IncompleteRecord,
}
pub type Result<T, E = Error> = std::result::Result<T, E>;

use bio_types::sequence::SequenceRead;

use crate::utils::TextSlice;

/// Trait for FastQ readers.
pub trait FastqRead {
    fn read(&mut self, record: &mut Record) -> Result<()>;
}

/// A FastQ reader.
#[derive(Debug)]
pub struct Reader<B> {
    reader: B,
    line_buffer: String,
}

impl Reader<io::BufReader<fs::File>> {
    /// Read from a given file.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(path.as_ref())
            .map_err(|e| Error::FileOpen {
                path: path.as_ref().to_owned(),
                source: e,
            })
            .map(Reader::new)
            .with_context(|| format!("Failed to read fastq from {:#?}", path))
    }
}

impl<R: io::Read> Reader<io::BufReader<R>> {
    /// Read from a given [`io::Read`](https://doc.rust-lang.org/std/io/trait.Read.html).
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line_buffer: String::new(),
        }
    }

    /// Create a new Fastq reader given a capacity and an instance of `io::Read`.
    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Reader {
            reader: io::BufReader::with_capacity(capacity, reader),
            line_buffer: String::new(),
        }
    }
}

impl<B> Reader<B>
where
    B: io::BufRead,
{
    ///  Create a new Fastq reader with an object that implements `io::BufReader`.
    pub fn from_bufread(bufreader: B) -> Self {
        Reader {
            reader: bufreader,
            line_buffer: String::new(),
        }
    }

    /// Return an iterator over the records of this FastQ file.
    ///
    /// # Errors
    ///
    /// This function will return an error if a record is incomplete
    /// or syntax is violated.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fastq;
    ///
    /// let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
    /// let records = fastq::Reader::new(fq)
    ///     .records()
    ///     .map(|record| record.unwrap());
    /// for record in records {
    ///     assert!(record.check().is_ok())
    /// }
    /// ```
    pub fn records(self) -> Records<B> {
        Records { reader: self }
    }
}

impl<B> FastqRead for Reader<B>
where
    B: io::BufRead,
{
    /// Read the next FastQ entry into the given [`Record`](#Record).
    /// An empty record indicates that no more records can be read.
    ///
    /// This method is useful when you want to read records as fast as
    /// possible because it allows the reuse of a `Record` allocation.
    ///
    /// A more ergonomic approach to reading FastQ records is the
    /// [records](Reader::records) iterator.
    ///
    /// FastQ files with wrapped sequence and quality strings are allowed.
    ///
    /// # Errors
    ///
    /// This function will return an error if the record is incomplete,
    /// syntax is violated or any form of I/O error is encountered.
    /// Additionally, if the FastQ file has line-wrapped records, and the wrapping is not
    /// consistent between the sequence and quality string for a record, parsing will fail.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fastq::Record;
    /// use bio::io::fastq::{FastqRead, Reader};
    /// const FASTQ_FILE: &'static [u8] = b"@id desc
    /// AAAA
    /// +
    /// IIII
    /// ";
    /// let mut reader = Reader::new(FASTQ_FILE);
    /// let mut record = Record::new();
    ///
    /// reader.read(&mut record).unwrap();
    ///
    /// assert_eq!(record.id(), "id");
    /// assert_eq!(record.desc().unwrap(), "desc");
    /// assert_eq!(record.seq().to_vec(), b"AAAA");
    /// assert_eq!(record.qual().to_vec(), b"IIII");
    /// ```
    fn read(&mut self, record: &mut Record) -> Result<()> {
        record.clear();
        self.line_buffer.clear();

        self.reader.read_line(&mut self.line_buffer)?;

        if !self.line_buffer.is_empty() {
            if !self.line_buffer.starts_with('@') {
                return Err(Error::MissingAt);
            }
            let mut header_fields = self.line_buffer[1..].trim_end().splitn(2, ' ');
            record.id = header_fields.next().unwrap_or_default().to_owned();
            record.desc = header_fields.next().map(|s| s.to_owned());
            self.line_buffer.clear();

            self.reader.read_line(&mut self.line_buffer)?;

            let mut lines_read = 0;
            while !self.line_buffer.is_empty() && !self.line_buffer.starts_with('+') {
                record.seq.push_str(self.line_buffer.trim_end());
                self.line_buffer.clear();
                self.reader.read_line(&mut self.line_buffer)?;
                lines_read += 1;
            }

            for _ in 0..lines_read {
                self.line_buffer.clear();
                self.reader
                    .read_line(&mut self.line_buffer)
                    .map_err(Error::ReadError)?;
                record.qual.push_str(self.line_buffer.trim_end());
            }

            if record.qual.is_empty() {
                return Err(Error::IncompleteRecord);
            }
        }

        Ok(())
    }
}

/// A FastQ record.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
    qual: String,
}

impl Record {
    /// Create a new, empty FastQ record.
    pub fn new() -> Self {
        Record {
            id: String::new(),
            desc: None,
            seq: String::new(),
            qual: String::new(),
        }
    }

    /// Create a new FastQ record from given attributes.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fastq::Record;
    ///
    /// let record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", b"QQQQQQQ");
    /// assert_eq!(record.id(), "id_str");
    /// assert_eq!(record.desc(), Some("desc"));
    /// assert_eq!(record.seq(), b"ATGCGGG");
    /// assert_eq!(record.qual(), b"QQQQQQQ");
    /// ```
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: TextSlice<'_>, qual: &[u8]) -> Self {
        let desc = desc.map(|desc| desc.to_owned());
        Record {
            id: id.to_owned(),
            desc,
            seq: String::from_utf8(seq.to_vec()).unwrap(),
            qual: String::from_utf8(qual.to_vec()).unwrap(),
        }
    }

    /// Check if a record is empty.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fastq::Record;
    ///
    /// let mut record = Record::new();
    /// assert!(record.is_empty());
    ///
    /// record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", b"QQQQQQQ");
    /// assert!(!record.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty() && self.qual.is_empty()
    }

    /// Check the validity of a FastQ record.
    ///
    /// # Errors
    /// This function will return an `Err` if one of the following conditions is met:
    /// -   The record identifier is empty.
    /// -   There is a non-ASCII character found in either the sequence or quality strings.
    /// -   The sequence and quality strings do not have the same length.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fastq::Record;
    ///
    /// let mut record = Record::with_attrs("id", None, "Prüfung".as_ref(), b"!!!!!!!");
    /// let actual = record.check().unwrap_err();
    /// let expected = "Non-ascii character found in sequence.";
    /// assert_eq!(actual, expected);
    ///
    /// record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", b"QQQQQQQ");
    /// assert!(record.check().is_ok());
    /// ```
    pub fn check(&self) -> Result<(), &str> {
        if self.id().is_empty() {
            return Err("Expecting id for FastQ record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }
        if !self.qual.is_ascii() {
            return Err("Non-ascii character found in qualities.");
        }
        if self.seq().len() != self.qual().len() {
            return Err("Unequal length of sequence an qualities.");
        }

        Ok(())
    }

    /// Return the id of the record.
    pub fn id(&self) -> &str {
        self.id.as_ref()
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Option<&str> {
        match self.desc.as_ref() {
            Some(desc) => Some(desc),
            None => None,
        }
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> TextSlice<'_> {
        self.seq.trim_end().as_bytes()
    }

    /// Return the base qualities of the record.
    pub fn qual(&self) -> &[u8] {
        self.qual.trim_end().as_bytes()
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
        self.qual.clear();
    }
}

impl fmt::Display for Record {
    /// Allows for using `Record` in a given formatter `f`. In general this is for
    /// creating a `String` representation of a `Record` and, optionally, writing it to
    /// a file.
    ///
    /// # Errors
    /// Returns [`std::fmt::Error`](https://doc.rust-lang.org/std/fmt/struct.Error.html)
    /// if there is an issue formatting to the stream.
    ///
    /// # Examples
    ///
    /// Read in a Fastq `Record` and create a `String` representation of it.
    ///
    /// ```rust
    /// use bio::io::fastq::Reader;
    /// use std::fmt::Write;
    /// // create a "fake" fastq file
    /// let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
    /// let mut records = Reader::new(fq).records().map(|r| r.unwrap());
    /// let record = records.next().unwrap();
    ///
    /// let mut actual = String::new();
    /// // populate `actual` with a string representation of our record
    /// write!(actual, "{}", record).unwrap();
    ///
    /// let expected = std::str::from_utf8(fq).unwrap();
    ///
    /// assert_eq!(actual, expected)
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let header = match self.desc() {
            Some(d) => format!("{} {}", self.id().to_owned(), d),
            None => self.id().to_owned(),
        };
        write!(
            f,
            "@{}\n{}\n+\n{}\n",
            header,
            std::str::from_utf8(self.seq()).unwrap(),
            std::str::from_utf8(self.qual()).unwrap()
        )
    }
}

impl SequenceRead for Record {
    fn name(&self) -> &[u8] {
        self.id.as_bytes()
    }

    fn base(&self, i: usize) -> u8 {
        self.seq()[i]
    }

    fn base_qual(&self, i: usize) -> u8 {
        self.qual()[i]
    }

    fn len(&self) -> usize {
        self.seq().len()
    }
}

/// An iterator over the records of a FastQ file.
#[derive(Debug)]
pub struct Records<R: io::Read> {
    reader: Reader<R>,
}

impl<B> Iterator for Records<B>
where
    B: io::BufRead,
{
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Ok(()) if record.is_empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}

/// A FastQ writer.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    writer: io::BufWriter<W>,
}

impl Writer<fs::File> {
    /// Write to a given file path.
    #[allow(clippy::wrong_self_convention)]
    pub fn to_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::create(path).map(Writer::new)
    }

    /// Write to the given file path and a buffer capacity
    pub fn to_file_with_capacity<P: AsRef<Path>>(capacity: usize, path: P) -> io::Result<Self> {
        fs::File::create(path).map(|file| Writer::with_capacity(capacity, file))
    }
}

impl<W: io::Write> Writer<W> {
    /// Write to a given `io::Write`.
    pub fn new(writer: W) -> Self {
        Writer {
            writer: io::BufWriter::new(writer),
        }
    }

    /// Create a new Fastq writer with a capacity of write buffer
    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            writer: io::BufWriter::with_capacity(capacity, writer),
        }
    }

    /// Create a new Fastq writer with a given BufWriter
    pub fn from_bufwriter(bufwriter: io::BufWriter<W>) -> Self {
        Writer { writer: bufwriter }
    }

    /// Directly write a FastQ record.
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write(record.id(), record.desc(), record.seq(), record.qual())
    }

    /// Write a FastQ record with given id, optional description, sequence and qualities.
    pub fn write(
        &mut self,
        id: &str,
        desc: Option<&str>,
        seq: TextSlice<'_>,
        qual: &[u8],
    ) -> io::Result<()> {
        self.writer.write_all(b"@")?;
        self.writer.write_all(id.as_bytes())?;
        if let Some(desc) = desc {
            self.writer.write_all(b" ")?;
            self.writer.write_all(desc.as_bytes())?;
        }
        self.writer.write_all(b"\n")?;
        self.writer.write_all(seq)?;
        self.writer.write_all(b"\n+\n")?;
        self.writer.write_all(qual)?;
        self.writer.write_all(b"\n")?;

        Ok(())
    }

    /// Flush the writer, ensuring that everything is written.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write as FmtWrite;
    use std::io;

    const FASTQ_FILE: &[u8] = b"@id desc
ACCGTAGGCTGA
+
IIIIIIJJJJJJ
";

    #[test]
    fn test_reader() {
        let reader = Reader::new(FASTQ_FILE);
        let records: Vec<Result<Record>> = reader.records().collect();
        assert_eq!(records.len(), 1);
        for res in records {
            let record = res.unwrap();
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), "id");
            assert_eq!(record.desc(), Some("desc"));
            assert_eq!(record.seq(), b"ACCGTAGGCTGA");
            assert_eq!(record.qual(), b"IIIIIIJJJJJJ");
        }

        let reader = Reader::with_capacity(100, FASTQ_FILE);
        let records: Vec<Result<Record>> = reader.records().collect();
        assert_eq!(records.len(), 1);
        for res in records {
            let record = res.unwrap();
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), "id");
            assert_eq!(record.desc(), Some("desc"));
            assert_eq!(record.seq(), b"ACCGTAGGCTGA");
            assert_eq!(record.qual(), b"IIIIIIJJJJJJ");
        }

        let reader = Reader::from_bufread(io::BufReader::new(FASTQ_FILE));
        let records: Vec<Result<Record>> = reader.records().collect();
        assert_eq!(records.len(), 1);
        for res in records {
            let record = res.unwrap();
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), "id");
            assert_eq!(record.desc(), Some("desc"));
            assert_eq!(record.seq(), b"ACCGTAGGCTGA");
            assert_eq!(record.qual(), b"IIIIIIJJJJJJ");
        }
    }

    #[test]
    fn test_display_record_no_desc_id_without_space_after() {
        let fq: &'static [u8] = b"@id\nACGT\n+\n!!!!\n";
        let mut records = Reader::new(fq).records().map(|r| r.unwrap());
        let record = records.next().unwrap();
        let mut actual = String::new();
        write!(actual, "{}", record).unwrap();

        let expected = std::str::from_utf8(fq).unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_display_record_with_desc_id_has_space_between_id_and_desc() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
        let mut records = Reader::new(fq).records().map(|r| r.unwrap());
        let record = records.next().unwrap();
        let mut actual = String::new();
        write!(actual, "{}", record).unwrap();

        let expected = std::str::from_utf8(fq).unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_fqread_trait() {
        let path = "reads.fq.gz";
        let mut fq_reader: Box<dyn FastqRead> = match path.ends_with(".gz") {
            true => Box::new(Reader::new(io::BufReader::new(FASTQ_FILE))),
            false => Box::new(Reader::new(FASTQ_FILE)),
        };
        // The read method can be called, since it is implemented by
        // `Read`. Right now, the records method would not work.
        let mut record = Record::new();
        fq_reader.read(&mut record).unwrap();
        // Check if the returned result is correct.
        assert_eq!(record.check(), Ok(()));
        assert_eq!(record.id(), "id");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ACCGTAGGCTGA");
        assert_eq!(record.qual(), b"IIIIIIJJJJJJ");
    }

    #[test]
    fn test_record_with_attrs() {
        let record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", b"QQQQQQQ");
        assert_eq!(record.id(), "id_str");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ATGCGGG");
        assert_eq!(record.qual(), b"QQQQQQQ");
    }

    #[test]
    fn test_writer() {
        let mut writer = Writer::new(Vec::new());
        writer
            .write("id", Some("desc"), b"ACCGTAGGCTGA", b"IIIIIIJJJJJJ")
            .expect("Expected successful write");
        writer.flush().expect("Expected successful write");
        assert_eq!(writer.writer.get_ref(), &FASTQ_FILE);

        let mut writer = Writer::with_capacity(100, Vec::new());
        writer
            .write("id", Some("desc"), b"ACCGTAGGCTGA", b"IIIIIIJJJJJJ")
            .expect("Expected successful write");
        writer.flush().expect("Expected successful write");
        assert_eq!(writer.writer.get_ref(), &FASTQ_FILE);

        let mut writer = Writer::from_bufwriter(std::io::BufWriter::with_capacity(100, Vec::new()));
        writer
            .write("id", Some("desc"), b"ACCGTAGGCTGA", b"IIIIIIJJJJJJ")
            .expect("Expected successful write");
        writer.flush().expect("Expected successful write");
        assert_eq!(writer.writer.get_ref(), &FASTQ_FILE);
    }

    #[test]
    fn test_check_record_id_is_empty_raises_err() {
        let record = Record::with_attrs("", None, b"ACGT", b"!!!!");

        let actual = record.check().unwrap_err();
        let expected = "Expecting id for FastQ record.";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_record_seq_is_not_ascii_raises_err() {
        let record = Record::with_attrs("id", None, "Prüfung".as_ref(), b"!!!!");

        let actual = record.check().unwrap_err();
        let expected = "Non-ascii character found in sequence.";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_record_quality_is_not_ascii_raises_err() {
        let record = Record::with_attrs("id", None, b"ACGT", "Qualität".as_ref());

        let actual = record.check().unwrap_err();
        let expected = "Non-ascii character found in qualities.";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_record_quality_and_seq_diff_len_raises_err() {
        let record = Record::with_attrs("id", None, b"ACGT", b"!!!");

        let actual = record.check().unwrap_err();
        let expected = "Unequal length of sequence an qualities.";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_valid_record() {
        let record = Record::with_attrs("id", None, b"ACGT", b"!!!!");

        assert!(record.check().is_ok())
    }

    #[test]
    fn test_read_header_does_not_start_with_correct_char_raises_err() {
        let fq: &'static [u8] = b">id description\nACGT\n+\n!!!!\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();

        let error = reader.read(&mut record).unwrap_err();

        assert!(matches!(error, Error::MissingAt))
    }

    #[test]
    fn test_read_quality_is_empty_raises_err() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();

        let error = reader.read(&mut record).unwrap_err();

        assert!(matches!(error, Error::IncompleteRecord))
    }

    #[test]
    fn test_read_sequence_and_quality_are_wrapped_is_handled_with_one_sequence() {
        let fq: &'static [u8] = b"@id description\nACGT\nGGGG\nC\n+\n@@@@\n!!!!\n$\n";
        let mut reader = Reader::new(fq);

        let mut actual = Record::new();
        reader.read(&mut actual).unwrap();
        let expected = Record::with_attrs("id", Some("description"), b"ACGTGGGGC", b"@@@@!!!!$");

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_read_sequence_and_quality_are_wrapped_is_handled_with_three_sequences() {
        let fq: &'static [u8] = b"@id description\nACGT\nGGGG\nC\n+\n@@@@\n!!!!\n$\n@id2 description\nACGT\nGGGG\nC\n+\n@@@@\n!!!!\n$\n@id3 desc1 desc2\nAAA\nAAA\nAA\n+\n^^^\n^^^\n^^\n";
        let mut reader = Reader::new(fq);

        let mut actual = Record::new();
        reader.read(&mut actual).unwrap();
        let expected = Record::with_attrs("id", Some("description"), b"ACGTGGGGC", b"@@@@!!!!$");

        assert_eq!(actual, expected);

        reader.read(&mut actual).unwrap();
        let expected = Record::with_attrs("id2", Some("description"), b"ACGTGGGGC", b"@@@@!!!!$");

        assert_eq!(actual, expected);

        reader.read(&mut actual).unwrap();
        let expected = Record::with_attrs("id3", Some("desc1 desc2"), b"AAAAAAAA", b"^^^^^^^^");

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_read_wrapped_record_with_inconsistent_wrapping_errors() {
        let fq: &'static [u8] = b"@id description\nACGT\nGGGG\nC\n+\n@@@@\n!!!!$\n@id2 description\nACGT\nGGGG\nC\n+\n@@@@\n!!!!\n$\n@id3 desc1 desc2\nAAA\nAAA\nAA\n+\n^^^\n^^^\n^^\n";
        let mut reader = Reader::new(fq);

        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        let error = reader.read(&mut record).unwrap_err();

        assert!(matches!(error, Error::MissingAt))
    }

    #[test]
    fn test_record_iterator_next_read_returns_err_causes_next_to_return_some_err() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n";
        let mut records = Reader::new(fq).records();

        let error = records.next().unwrap().unwrap_err();

        assert!(matches!(error, Error::IncompleteRecord));
    }

    #[test]
    fn test_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.fq");
        let error = Reader::from_file(path)
            .unwrap_err()
            .downcast::<String>()
            .unwrap();

        assert_eq!(&error, "Failed to read fastq from \"/I/dont/exist.fq\"")
    }

    #[test]
    fn test_reader_from_file_path_exists_returns_ok() {
        let path = Path::new("Cargo.toml");

        assert!(Reader::from_file(path).is_ok())
    }

    #[test]
    fn test_sequence_read_for_record_trait_method_name() {
        let record = Record::with_attrs("id", None, b"ACGT", b"!!!!");

        let actual = record.name();
        let expected = b"id";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_sequence_read_for_record_trait_method_base_idx_in_range() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        let idx = 2;

        let actual = record.base(idx);
        let expected = b'G';

        assert_eq!(actual, expected)
    }

    #[test]
    #[should_panic]
    fn test_sequence_read_for_record_trait_method_base_idx_out_of_range() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        // idx 4 is where the newline character would be - we dont want that included
        let idx = 4;

        record.base(idx);
    }

    #[test]
    fn test_sequence_read_for_record_trait_method_base_qual_idx_in_range() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        let idx = 2;

        let actual = record.base_qual(idx);
        let expected = b'!';

        assert_eq!(actual, expected)
    }

    #[test]
    #[should_panic]
    fn test_sequence_read_for_record_trait_method_base_qual_idx_out_of_range() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        // idx 4 is where the newline character would be - we dont want that included
        let idx = 4;

        record.base_qual(idx);
    }

    #[test]
    fn test_sequence_read_for_record_trait_method_len() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n!!!!\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();

        let actual = record.len();
        let expected = 4;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_read_with_missing_plus() {
        let fq: &'static [u8] = b"@id description\nACGT\n*\n!!!!\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();
        let err = reader.read(&mut record).unwrap_err();

        assert!(matches!(err, Error::IncompleteRecord))
    }

    #[test]
    fn test_writer_to_file_dir_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.fq");

        let actual = Writer::to_file(path).unwrap_err();
        let expected = io::Error::new(io::ErrorKind::NotFound, "foo");

        assert_eq!(actual.kind(), expected.kind());
    }

    #[test]
    fn test_writer_to_file_dir_exists_returns_ok() {
        let file = tempfile::NamedTempFile::new().expect("Could not create temp file");
        let path = file.path();

        assert!(Writer::to_file(path).is_ok())
    }

    #[test]
    fn test_write_record() {
        let path = Path::new("test.fq");
        let file = fs::File::create(path).unwrap();
        {
            let handle = io::BufWriter::new(file);
            let mut writer = Writer { writer: handle };
            let record = Record::with_attrs("id", Some("desc"), b"ACGT", b"!!!!");

            let write_result = writer.write_record(&record);
            assert!(write_result.is_ok());
        }

        let actual = fs::read_to_string(path).unwrap();
        let expected = "@id desc\nACGT\n+\n!!!!\n";

        assert!(fs::remove_file(path).is_ok());
        assert_eq!(actual, expected)
    }
}
