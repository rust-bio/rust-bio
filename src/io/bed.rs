// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! BED format reading and writing.
//!
//! # Example
//!
//! ```
//! use bio::io::bed;
//! let example = b"1\t5\t5000\tname1\t0.5";
//! let mut reader = bed::Reader::new(&example[..]);
//! let mut writer = bed::Writer::new(vec![]);
//! for record in reader.records() {
//!     let rec = record.ok().expect("Error reading record.");
//!     println!("{}", rec.chrom());
//!     writer.write(&rec).ok().expect("Error writing record.");
//! }
//! ```

use std::convert::AsRef;
use std::fs;
use std::io;
use std::path::Path;

use csv;

use utils::Strand;

/// A BED reader.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}

impl<R: io::Read> Reader<R> {
    /// Read from a given reader.
    pub fn new(reader: R) -> Self {
        Reader {
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(reader),
        }
    }

    /// Iterate over all records.
    pub fn records(&mut self) -> Records<R> {
        Records {
            inner: self.inner.deserialize(),
        }
    }
}

/// A BED record.
pub struct Records<'a, R: 'a + io::Read> {
    inner: csv::DeserializeRecordsIter<'a, R, (String, u64, u64, Option<Vec<String>>)>,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner.next().map(|res| {
            res.map(|(chrom, start, end, aux)| Record {
                chrom,
                start,
                end,
                aux: aux.unwrap_or(Vec::new()),
            })
        })
    }
}

/// A BED writer.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
}

impl Writer<fs::File> {
    /// Write to a given file path.
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

    /// Write a given BED record.
    pub fn write(&mut self, record: &Record) -> csv::Result<()> {
        if record.aux.is_empty() {
            self.inner
                .serialize(&(&record.chrom, record.start, record.end))
        } else {
            self.inner
                .serialize(&(&record.chrom, record.start, record.end, &record.aux))
        }
    }
}

/// A BED record as defined by BEDtools
/// (http://bedtools.readthedocs.org/en/latest/content/general-usage.html)
#[derive(Debug, Serialize, Default)]
pub struct Record {
    chrom: String,
    start: u64,
    end: u64,
    aux: Vec<String>,
}

impl Record {
    /// Create a new BED record.
    pub fn new() -> Self {
        Record {
            chrom: "".to_owned(),
            start: 0,
            end: 0,
            aux: vec![],
        }
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

    /// Name of the feature.
    pub fn name(&self) -> Option<&str> {
        self.aux(3)
    }

    /// Score of the feature.
    pub fn score(&self) -> Option<&str> {
        self.aux(4)
    }

    /// Strand of the feature.
    pub fn strand(&self) -> Option<Strand> {
        match self.aux(5) {
            Some("+") => Some(Strand::Forward),
            Some("-") => Some(Strand::Reverse),
            _ => None,
        }
    }

    /// Access auxilliary fields after the strand field by index
    /// (counting first field (chromosome) as 0).
    pub fn aux(&self, i: usize) -> Option<&str> {
        let j = i - 3;
        if j < self.aux.len() {
            Some(&self.aux[j])
        } else {
            None
        }
    }

    /// Set chromosome.
    pub fn set_chrom(&mut self, chrom: &str) {
        self.chrom = chrom.to_owned();
    }

    /// Set start of feature.
    pub fn set_start(&mut self, start: u64) {
        self.start = start;
    }

    /// Set end of feature.
    pub fn set_end(&mut self, end: u64) {
        self.end = end;
    }

    /// Set name.
    pub fn set_name(&mut self, name: &str) {
        if self.aux.len() < 1 {
            self.aux.push(name.to_owned());
        } else {
            self.aux[0] = name.to_owned();
        }
    }

    /// Set score.
    pub fn set_score(&mut self, score: &str) {
        if self.aux.len() < 1 {
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

#[cfg(test)]
mod tests {
    use super::*;

    const BED_FILE: &'static [u8] = b"1\t5\t5000\tname1\tup
2\t3\t5005\tname2\tup
";
    const BED_FILE_COMPACT: &'static [u8] = b"1\t5\t5000\n2\t3\t5005\n";

    #[test]
    fn test_reader() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];

        let mut reader = Reader::new(BED_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);
        }
    }

    // TODO enable this test case once compact bed file reading has been fixed, see
    // https://github.com/rust-bio/rust-bio/pull/156/files#r157506929
    // #[test]
    // /// Test for 'compact' BED files which only have chrom, start, and stop fields.
    // fn test_reader_compact() {
    //     let chroms = ["1", "2"];
    //     let starts = [5, 3];
    //     let ends = [5000, 5005];
    //
    //     let mut reader = Reader::new(BED_FILE_COMPACT);
    //     for (i, r) in reader.records().enumerate() {
    //         let record = r.unwrap();
    //         assert_eq!(record.chrom(), chroms[i]);
    //         assert_eq!(record.start(), starts[i]);
    //         assert_eq!(record.end(), ends[i]);
    //     }
    // }

    #[test]
    fn test_writer() {
        let mut reader = Reader::new(BED_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records() {
            writer
                .write(&r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), BED_FILE);
    }
}
