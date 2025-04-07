// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! BEDPE format reading and writing.
//!
//! # Example
//!
//! ```
//! use bio::io::bedpe;
//! let example = b"1\t5\t500\t1\t600\t1000\tname1\t0.5";
//! let mut reader = bedpe::Reader::new(&example[..]);
//! let mut writer = bedpe::Writer::new(vec![]);
//! for record in reader.records() {
//!     let rec = record.expect("Error reading record.");
//!     println!("{}", rec.chrom1());
//!     writer.write(&rec).expect("Error writing record.");
//! }
//! ```

use std::convert::AsRef;
use std::fs;
use std::io;
use std::marker::Copy;
use std::ops::Deref;
use std::path::Path;

use anyhow::Context;
use bio_types::annot;
use bio_types::annot::loc::Loc;
use bio_types::strand;

/// A BEDPE reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to read bedpe from {:#?}", path))
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
    pub fn records(&mut self) -> Records<'_, R> {
        Records {
            inner: self.inner.deserialize(),
        }
    }
}

/// An iterator over the records of a BEDPE file.
pub struct Records<'a, R: io::Read> {
    inner: csv::DeserializeRecordsIter<'a, R, Record>,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner.next()
    }
}

/// A BEDPE writer.
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
                .from_writer(writer),
        }
    }

    /// Write a given BEDPE record.
    pub fn write(&mut self, record: &Record) -> csv::Result<()> {
        if record.aux.is_empty() {
            self.inner.serialize((
                &record.chrom1,
                record.start1,
                record.end1,
                &record.chrom2,
                record.start2,
                record.end2,
            ))
        } else {
            self.inner.serialize((
                &record.chrom1,
                record.start1,
                record.end1,
                &record.chrom2,
                record.start2,
                record.end2,
                &record.aux,
            ))
        }
    }
}

/// A BEDPE record as defined by BEDtools
/// (http://bedtools.readthedocs.org/en/latest/content/general-usage.html)
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Record {
    chrom1: String,
    start1: u64,
    end1: u64,
    chrom2: String,
    start2: u64,
    end2: u64,
    #[serde(default)]
    aux: Vec<String>,
}

impl Record {
    /// Create a new BEDPE record.
    pub fn new() -> Self {
        Record {
            chrom1: "".to_owned(),
            start1: 0,
            end1: 0,
            chrom2: "".to_owned(),
            start2: 0,
            end2: 0,
            aux: vec![],
        }
    }

    /// Chromosome of first end of the feature.
    pub fn chrom1(&self) -> &str {
        &self.chrom1
    }

    /// Start position of first end of the feature (0-based).
    pub fn start1(&self) -> u64 {
        self.start1
    }

    /// End position of first end of the feature (0-based, not included).
    pub fn end1(&self) -> u64 {
        self.end1
    }

    /// Chromosome of second end of the feature.
    pub fn chrom2(&self) -> &str {
        &self.chrom2
    }

    /// Start position of second end of the feature (0-based).
    pub fn start2(&self) -> u64 {
        self.start2
    }

    /// End position of second end of the feature (0-based, not included).
    pub fn end2(&self) -> u64 {
        self.end2
    }

    /// Name of the feature.
    pub fn name(&self) -> Option<&str> {
        self.aux(6)
    }

    /// Score of the feature.
    pub fn score(&self) -> Option<&str> {
        self.aux(7)
    }

    /// Strand of first end of the feature.
    pub fn strand1(&self) -> Option<strand::Strand> {
        match self.aux(8) {
            Some("+") => Some(strand::Strand::Forward),
            Some("-") => Some(strand::Strand::Reverse),
            _ => None,
        }
    }

    /// Strand of second end of the feature.
    pub fn strand2(&self) -> Option<strand::Strand> {
        match self.aux(9) {
            Some("+") => Some(strand::Strand::Forward),
            Some("-") => Some(strand::Strand::Reverse),
            _ => None,
        }
    }

    /// Access auxilliary fields after the strand field by index
    /// (counting first field (chromosome) as 0).
    pub fn aux(&self, i: usize) -> Option<&str> {
        let j = i - 6;
        if j < self.aux.len() {
            Some(&self.aux[j])
        } else {
            None
        }
    }

    /// Set chromosome of first end of the feature.
    pub fn set_chrom1(&mut self, chrom: &str) {
        self.chrom1 = chrom.to_owned();
    }

    /// Set start of first end of the feature.
    pub fn set_start1(&mut self, start: u64) {
        self.start1 = start;
    }

    /// Set end of first end of the feature.
    pub fn set_end1(&mut self, end: u64) {
        self.end1 = end;
    }

    /// Set chromosome of second end of the feature.
    pub fn set_chrom2(&mut self, chrom: &str) {
        self.chrom2 = chrom.to_owned();
    }

    /// Set start of second end of the feature.
    pub fn set_start2(&mut self, start: u64) {
        self.start2 = start;
    }

    /// Set end of second end of the feature.
    pub fn set_end2(&mut self, end: u64) {
        self.end2 = end;
    }

    /// Set name.
    pub fn set_name(&mut self, name: &str) {
        if self.aux.is_empty() {
            self.aux.push(name.to_owned());
        } else {
            self.aux[0] = name.to_owned();
        }
    }

    /// Set score.
    pub fn set_score(&mut self, score: &str) {
        if self.aux.is_empty() {
            self.aux.push("".to_owned());
        }
        if self.aux.len() < 2 {
            self.aux.push(score.to_owned());
        } else {
            self.aux[1] = score.to_owned();
        }
    }

    /// Add auxilliary field. This has to happen after name and score have been set.
    pub fn push_aux(&mut self, field: &str) {
        self.aux.push(field.to_owned());
    }
}

/// Generate a BEDPE format `Record` for an annotation position.
///
/// This record will have length 1, and when created it will have an
/// empty name.
impl<R, S> From<(annot::pos::Pos<R, S>, annot::pos::Pos<R, S>)> for Record
where
    R: Deref<Target = str>,
    S: Into<strand::Strand> + Copy,
{
    fn from((pos1, pos2): (annot::pos::Pos<R, S>, annot::pos::Pos<R, S>)) -> Self {
        let mut bedpe = Record::new();
        bedpe.set_chrom1(pos1.refid());
        bedpe.set_start1(pos1.pos() as u64);
        bedpe.set_end1((pos1.pos() + 1) as u64);
        bedpe.set_chrom2(pos2.refid());
        bedpe.set_start2(pos2.pos() as u64);
        bedpe.set_end2((pos2.pos() + 1) as u64);
        bedpe.set_name("");
        bedpe.set_score("0");
        bedpe.push_aux(pos1.strand().into().strand_symbol());
        bedpe.push_aux(pos2.strand().into().strand_symbol());
        bedpe
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bio_types::annot::pos::Pos;
    use bio_types::strand::{ReqStrand, Strand};

    const BEDPE_FILE: &[u8] = b"\
1\t5\t5000\t1\t6000\t10000\tname1\tup
2\t3\t5005\t2\t10\t5500\tname2\tup
";
    const BEDPE_FILE_COMMENT: &[u8] = b"\
# this line should be ignored
1\t5\t5000\t1\t6000\t10000\tname1\tup
# and this one as well
2\t3\t5005\t2\t10\t5500\tname2\tup
";
    const BEDPE_FILE_COMPACT: &[u8] = b"1\t5\t5000\t1\t6000\t10000\n2\t3\t5005\t2\t10\t5500\n";

    #[test]
    fn test_reader() {
        let chroms = [("1", "1"), ("2", "2")];
        let starts = [(5, 6000), (3, 10)];
        let ends = [(5000, 10000), (5005, 5500)];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];

        let mut reader = Reader::new(BEDPE_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!((record.chrom1(), record.chrom2()), chroms[i]);
            assert_eq!((record.start1(), record.start2()), starts[i]);
            assert_eq!((record.end1(), record.end2()), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);
        }
    }

    #[test]
    fn test_reader_with_comment() {
        let chroms = [("1", "1"), ("2", "2")];
        let starts = [(5, 6000), (3, 10)];
        let ends = [(5000, 10000), (5005, 5500)];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];

        let mut reader = Reader::new(BEDPE_FILE_COMMENT);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!((record.chrom1(), record.chrom2()), chroms[i]);
            assert_eq!((record.start1(), record.start2()), starts[i]);
            assert_eq!((record.end1(), record.end2()), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);
        }
    }

    #[test]
    fn test_reader_compact() {
        let chroms = [("1", "1"), ("2", "2")];
        let starts = [(5, 6000), (3, 10)];
        let ends = [(5000, 10000), (5005, 5500)];

        let mut reader = Reader::new(BEDPE_FILE_COMPACT);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!((record.chrom1(), record.chrom2()), chroms[i]);
            assert_eq!((record.start1(), record.start2()), starts[i]);
            assert_eq!((record.end1(), record.end2()), ends[i]);
        }
    }

    #[test]
    fn test_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.bedpe");
        let error = Reader::from_file(path)
            .unwrap_err()
            .downcast::<String>()
            .unwrap();

        assert_eq!(&error, "Failed to read bedpe from \"/I/dont/exist.bedpe\"")
    }

    #[test]
    fn test_writer() {
        let mut reader = Reader::new(BEDPE_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records() {
            writer
                .write(&r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), BEDPE_FILE);
    }

    #[test]
    fn test_bed_from_pos() {
        let pos1 = Pos::new("chrXI".to_owned(), 334412, ReqStrand::Forward);
        let pos2 = Pos::new("chrXI".to_owned(), 300000, ReqStrand::Reverse);

        let record = Record::from((pos1, pos2));

        assert_eq!(record.chrom1(), String::from("chrXI"));
        assert_eq!(record.chrom2(), String::from("chrXI"));
        assert_eq!(record.start1(), 334412);
        assert_eq!(record.end1(), 334412 + 1);
        assert_eq!(record.start2(), 300000);
        assert_eq!(record.end2(), 300000 + 1);
        assert_eq!(record.name(), Some(""));
        assert_eq!(record.score(), Some("0"));
        assert_eq!(record.strand1(), Some(Strand::Forward));
        assert_eq!(record.strand2(), Some(Strand::Reverse));
    }
}
