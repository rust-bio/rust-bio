// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! TSV format reading and writing.
//!
//! Implements a trait used for parsing TSV formats like BED

use std::convert::AsRef;
use std::fs;
use std::io;
use std::path::Path;

use anyhow::Context;
use serde::de::DeserializeOwned;
use serde::Deserialize;
use serde::Serialize;

/// A TSV reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to read from {:#?}", path))
    }
}

impl<R: io::Read> Reader<R> {
    /// Read from a given reader.
    pub fn new(reader: R) -> Self {
        Reader {
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_reader(reader),
        }
    }

    /// Iterate over all records.
    pub fn records<T: DeserializeOwned>(&mut self) -> Records<'_, R, T> {
        Records {
            inner: self.inner.deserialize(),
        }
    }
}

/// An iterator over the records of a TSV file.
pub struct Records<'a, R: io::Read, T: DeserializeOwned> {
    inner: csv::DeserializeRecordsIter<'a, R, T>,
}

impl<'a, R: io::Read, T: DeserializeOwned> Iterator for Records<'a, R, T> {
    type Item = csv::Result<T>;

    fn next(&mut self) -> Option<csv::Result<T>> {
        self.inner.next()
    }
}

/// A tsv writer.
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

    /// Write a given recor tsv.
    pub fn write(&mut self, record: impl IsRecord) -> csv::Result<()> {
        self.inner.serialize(record)
    }
}

pub trait IsRecord
where
    Self: Sized + Serialize,
{
}

/// A TSV record
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Record {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    #[serde(default)]
    pub aux: Vec<String>,
}

impl Default for Record {
    fn default() -> Self {
        Self {
            chrom: "".to_owned(),
            start: 0,
            end: 0,
            aux: vec![],
        }
    }
}

impl IsRecord for Record {}
impl IsRecord for &Record {}

impl Record {
    pub fn new() -> Self {
        Self::default()
    }

    /// Chromosome of the feature.
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Start position of feature (0-based).
    pub fn start(&self) -> u64 {
        self.start
    }

    /// End position of feature (0-based, not included).
    pub fn end(&self) -> u64 {
        self.end
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

    /// Set auxilliary field after the end field by index and return the previous value if present.
    /// If the index is bigger than the current aux vector, fill all not set fields with empty
    /// strings.
    /// (counting first field (chromosome) as 0).
    pub fn set_aux<S: Into<String>>(&mut self, i: usize, field: S) -> Option<String> {
        let j = i - 3;
        if j < self.aux.len() {
            let previous_value = self.aux[j].clone();
            self.aux[j] = field.into();
            Some(previous_value)
        } else {
            for _ in self.aux.len()..j {
                self.aux.push("".to_owned());
            }
            self.aux.push(field.into());
            None
        }
    }

    /// Set chromosome.
    pub fn set_chrom(&mut self, chrom: &str) {
        self.chrom = chrom.to_owned()
    }

    /// Set start of feature.
    pub fn set_start(&mut self, start: u64) {
        self.start = start
    }

    /// Set end of feature.
    pub fn set_end(&mut self, end: u64) {
        self.end = end
    }

    /// Add auxilliary field. This has to happen after chr, start, and end have been set.
    pub fn push_aux(&mut self, field: &str) {
        self.aux.push(field.to_owned());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BED_FILE: &[u8] = b"1\t5\t5000\tname1\tup
2\t3\t5005\tname2\tup
";
    const BED_FILE_COMMENT: &[u8] = b"\
# this line should be ignored
1\t5\t5000\tname1\tup
# and this one as well
2\t3\t5005\tname2\tup
";
    const BED_FILE_COMPACT: &[u8] = b"1\t5\t5000\n2\t3\t5005\n";

    struct TestReader<R: io::Read>(Reader<R>);

    impl<R: io::Read> TestReader<R> {
        fn new(source: R) -> Self {
            Self(Reader::new(source))
        }

        fn records(&mut self) -> Records<'_, R, TestRecord> {
            self.0.records()
        }
    }

    impl TestReader<fs::File> {
        /// Read from a given file path.
        fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
            Reader::from_file(path).map(Self)
        }
    }

    type TestRecord = Record;

    impl TestRecord {
        pub fn name(&self) -> Option<&str> {
            self.aux(3)
        }
    }

    #[test]
    fn test_core_reader() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];
        let aux_col_7 = ["test", "test1"];

        let mut reader = Reader::new(BED_FILE);
        for (i, r) in reader.records::<TestRecord>().enumerate() {
            let mut record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.aux(3).expect("Error reading name"), names[i]);
            assert_eq!(record.aux(4).expect("Error reading score"), scores[i]);
            record.set_aux(7, aux_col_7[i]);
            assert_eq!(record.aux(5).expect("Error reading empty field"), "");
            assert_eq!(record.aux(6).expect("Error reading empty field"), "");
            assert_eq!(record.aux(7).expect("Error reading aux col"), aux_col_7[i]);
            let old_name = record.set_aux(3, "new_name");
            assert_eq!(old_name, Some(names[i].to_owned()));
            assert_eq!(record.aux(3).expect("Error reading col 3"), "new_name");
        }
    }

    #[test]
    fn test_core_reader_with_comment() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];

        let mut reader = Reader::new(BED_FILE_COMMENT);
        for (i, r) in reader.records::<TestRecord>().enumerate() {
            let record = r.expect("Error reading record");
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
        for (i, r) in reader.records::<TestRecord>().enumerate() {
            let record = r.unwrap();
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
        }
    }

    #[test]
    fn test_core_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.bed");
        let error = Reader::from_file(path)
            .unwrap_err()
            .downcast::<String>()
            .unwrap();

        assert_eq!(&error, "Failed to read from \"/I/dont/exist.bed\"")
    }

    #[test]
    fn test_core_writer() {
        let mut reader = Reader::new(BED_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records::<TestRecord>() {
            writer
                .write(&r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), BED_FILE);
    }

    #[test]
    fn test_implemented_reader() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];

        let mut reader = TestReader::new(BED_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.aux(3).expect("Error reading name"), names[i]);
            assert_eq!(record.aux(4).expect("Error reading score"), scores[i]);
        }
    }

    #[test]
    fn test_writer_with_implemented_reader() {
        let mut reader = TestReader::new(BED_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records() {
            writer
                .write(&r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), BED_FILE);
    }
}
