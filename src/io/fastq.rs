// Copyright 2014-2018 Johannes Köster, Henning Timm.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! FastQ reading and writing.
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::fastq;
//! let reader = fastq::Reader::new(io::stdin());
//! ```

use std::convert::AsRef;
use std::fmt;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::path::Path;

use bio_types::sequence::SequenceRead;

use crate::utils::TextSlice;

/// Trait for FASTQ readers.
pub trait FastqRead {
    fn read(&mut self, record: &mut Record) -> io::Result<()>;
}

/// A FastQ reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    reader: io::BufReader<R>,
    line_buf: String,
}

impl Reader<fs::File> {
    /// Read from a given file.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}

impl<R: io::Read> Reader<R> {
    /// Read from a given `io::Read`.
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line_buf: String::new(),
        }
    }

    /// Return an iterator over the records of this FastQ file.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}

impl<R> FastqRead for Reader<R>
where
    R: io::Read,
{
    /// Read the next FASTQ entry into the given `Record`.
    /// An empty record indicates that no more records can be read.
    ///
    /// This method is useful when you want to read records as fast as
    /// possible because it allows the reuse of a `Record` allocation.
    ///
    /// A more ergonomic approach to reading FASTQ records, is the
    /// [records](Reader::records) iterator.
    ///
    /// # Errors
    ///
    /// This function will return an error if the record is incomplete,
    /// syntax is violated or any form of I/O error is encountered.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use std::error::Error;
    /// # use bio::io::fastq::{Reader, FastqRead};
    /// # use bio::io::fastq::Record;
    /// # fn main() -> Result<(), Box<Error>> {
    /// const fastq_file: &'static [u8] = b"@id desc
    /// AAAA
    /// +
    /// IIII
    /// ";
    /// let mut reader = Reader::new(fastq_file);
    /// let mut record = Record::new();
    ///
    /// // Check for errors parsing the record
    /// reader.read(&mut record)?;
    ///
    /// assert_eq!(record.id(), "id");
    /// assert_eq!(record.desc().unwrap(), "desc");
    /// assert_eq!(record.seq().to_vec(), b"AAAA");
    /// assert_eq!(record.qual().to_vec(), b"IIII");
    /// # Ok(())
    /// # }
    /// ```
    fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        self.line_buf.clear();

        self.reader.read_line(&mut self.line_buf)?;

        if !self.line_buf.is_empty() {
            if !self.line_buf.starts_with('@') {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Expected @ at record start.",
                ));
            }
            let mut header_fields = self.line_buf[1..].trim_end().splitn(2, ' ');
            record.id = header_fields.next().unwrap_or_default().to_owned();
            record.desc = header_fields.next().map(|s| s.to_owned());
            self.reader.read_line(&mut record.seq)?;
            self.reader.read_line(&mut self.line_buf)?;
            self.reader.read_line(&mut record.qual)?;
            if record.qual.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Incomplete record. Each FastQ record has to consist \
                     of 4 lines: header, sequence, separator and \
                     qualities.",
                ));
            }
        }

        Ok(())
    }
}

/// A FastQ record.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
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
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: TextSlice<'_>, qual: &[u8]) -> Self {
        let desc = match desc {
            Some(desc) => Some(desc.to_owned()),
            _ => None,
        };
        Record {
            id: id.to_owned(),
            desc,
            seq: String::from_utf8(seq.to_vec()).unwrap(),
            qual: String::from_utf8(qual.to_vec()).unwrap(),
        }
    }

    /// Check if record is empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty() && self.qual.is_empty()
    }

    /// Check validity of FastQ record.
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
            Some(desc) => Some(&desc),
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
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let header = match self.desc() {
            Some(d) => format!("{} {}", self.id, d),
            None => self.id.to_owned(),
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
        self.seq.as_bytes()[i]
    }

    fn base_qual(&self, i: usize) -> u8 {
        self.qual.as_bytes()[i]
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

impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
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
    pub fn to_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::create(path).map(Writer::new)
    }
}

impl<W: io::Write> Writer<W> {
    /// Write to a given `io::Write`.
    pub fn new(writer: W) -> Self {
        Writer {
            writer: io::BufWriter::new(writer),
        }
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
        if desc.is_some() {
            self.writer.write_all(b" ")?;
            self.writer.write_all(desc.unwrap().as_bytes())?;
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

    const FASTQ_FILE: &'static [u8] = b"@id desc
ACCGTAGGCTGA
+
IIIIIIJJJJJJ
";

    #[test]
    fn test_reader() {
        let reader = Reader::new(FASTQ_FILE);
        let records: Vec<io::Result<Record>> = reader.records().collect();
        assert!(records.len() == 1);
        for res in records {
            let record = res.ok().unwrap();
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
            .ok()
            .expect("Expected successful write");
        writer.flush().ok().expect("Expected successful write");
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

        let actual = reader.read(&mut record).unwrap_err();
        let expected = io::Error::new(io::ErrorKind::Other, "Expected @ at record start.");

        assert_eq!(actual.kind(), expected.kind());
        assert_eq!(actual.to_string(), expected.to_string())
    }

    #[test]
    fn test_read_quality_is_empty_raises_err() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n";
        let mut reader = Reader::new(fq);
        let mut record = Record::new();

        let actual = reader.read(&mut record).unwrap_err();
        let expected = io::Error::new(io::ErrorKind::Other, "Incomplete record. Each FastQ record has to consist of 4 lines: header, sequence, separator and qualities.");

        assert_eq!(actual.kind(), expected.kind());
        assert_eq!(actual.to_string(), expected.to_string())
    }

    #[test]
    fn test_record_iterator_next_read_returns_err_causes_next_to_return_some_err() {
        let fq: &'static [u8] = b"@id description\nACGT\n+\n";
        let mut records = Reader::new(fq).records();

        let actual = records.next().unwrap().unwrap_err();
        let expected = io::Error::new(io::ErrorKind::Other, "Incomplete record. Each FastQ record has to consist of 4 lines: header, sequence, separator and qualities.");

        assert_eq!(actual.kind(), expected.kind());
        assert_eq!(actual.to_string(), expected.to_string())
    }

    #[test]
    fn test_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.fq");

        let actual = Reader::from_file(path).unwrap_err();
        let expected = io::Error::new(io::ErrorKind::NotFound, "foo");

        assert_eq!(actual.kind(), expected.kind());
        assert!(actual.to_string().starts_with("No such file or directory"))
    }

    #[test]
    fn test_reader_from_file_path_exists_returns_ok() {
        let path = Path::new("Cargo.toml");

        assert!(Reader::from_file(path).is_ok())
    }
}
