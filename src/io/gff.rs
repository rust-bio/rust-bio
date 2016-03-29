// Copyright 2016 Pierre Marijon.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


//! GFF format reading and writing.
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::gff;
//! let reader = gff::Reader::new(io::stdin());
//! ```


use std::io;
use std::fs;
use std::fmt;
use std::path::Path;
use std::convert::AsRef;

use csv;


/// A GFF reader.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}


impl Reader<fs::File> {
    /// Read FASTA from given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}


impl<R: io::Read> Reader<R> {
    /// Create a new Fasta reader given an instance of `io::Read`.
    pub fn new(reader: R) -> Self {
        Reader {
            inner: csv::Reader::from_reader(reader).delimiter(b'\t').has_headers(false)
        }
    }

    /// Iterate over all records.
    pub fn records(&mut self) -> Records<R> {
        Records { inner: self.inner.decode() }
    }
}

/// A GFF record.
pub struct Records<'a, R: 'a + io::Read> {
    inner: csv::DecodedRecords<'a, R, (String, String, String, u64, u64, String, String, String, String)>,
}


impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner.next().map(|res| {
            res.map(|(sequence, source, type_name, start, end, score, strand, frame, attributes)| {
                Record {
                    sequence: sequence,
                    source: source,
                    type_name: type_name,
                    start: start,
                    end: end,
                    score: score,
                    strand: strand,
                    frame: frame,
                    attributes: attributes,
                }
            })
        })
    }
}


/// A GFF writer.
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
        Writer { inner: csv::Writer::from_writer(writer).delimiter(b'\t').flexible(true) }
    }

    /// Write a given GFF record.
    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        self.inner.encode(record)
    }
}


/// A GFF record
#[derive(RustcEncodable)]
pub struct Record {
    sequence: String,
    source: String,
    type_name: String,
    start: u64,
    end: u64,
    score: String,
    strand: String,
    frame: String,
    attributes: String,
}

impl Record {
    /// Create a new GFF record.
    pub fn new() -> Self {
        Record {
            sequence: "".to_owned(),
            source: "".to_owned(),
            type_name: "".to_owned(),
            start: 0,
            end: 0,
            score: ".".to_owned(),
            strand: ".".to_owned(),
            frame: "".to_owned(),
            attributes: "".to_owned(),
        }
    }

    /// Sequence of the feature.
    pub fn sequence(&self) -> &str {
        &self.sequence
    }

    /// Source of the feature.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Type name of the feature.
    pub fn type_name(&self) -> &str {
        &self.type_name
    }

    /// Start position of feature (1-based).
    pub fn start(&self) -> u64 {
        self.start
    }

    /// End position of feature (1-based, not included).
    pub fn end(&self) -> u64 {
        self.end
    }

    /// Score of feature
    pub fn score(&self) -> Option<u64> {
        match self.score.as_ref() {
            "." => None,
            _ => self.score.parse::<u64>().ok(),
        }
    }

    /// Strand of the feature.
    pub fn strand(&self) -> Option<Strand> {
        match self.strand.as_ref() {
            "+" => Some(Strand::Forward),
            "-" => Some(Strand::Reverse),
            _ => None,
        }
    }

    /// Frame of the feature.
    pub fn frame(&self) -> &str {
        &self.frame
    }

    /// Attribute of feature
    pub fn attributes(&self) -> &str {
        &self.attributes
    }
    
    /// Set sequence of feature.
    pub fn set_sequence(&mut self, sequence: &str) {
        self.sequence = sequence.to_owned();
    }

    /// Set source of feature.
    pub fn set_source(&mut self, source: &str) {
        self.source = source.to_owned();
    }

    /// Set type_name of feature.
    pub fn set_type_name(&mut self, type_name: &str) {
        self.type_name = type_name.to_owned();
    }

    /// Set start of feature.
    pub fn set_start(&mut self, start: u64) {
        self.start = start;
    }

    /// Set end of feature.
    pub fn set_end(&mut self, end: u64) {
        self.end = end;
    }

    /// Set score of feature.
    pub fn set_score(&mut self, score: &str) {
        self.score = score.to_owned();
    }

    /// Set strand of feature.
    pub fn set_strand(&mut self, strand: &str) {
        self.strand = strand.to_owned();
    }
}

/// Strand information.
#[derive(Debug, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[cfg(test)]
mod tests {
    use super::*;

    const GFF_FILE: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote=Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID=PRO_0000148105;Note=ATP-dependent protease subunit HslV
";

    #[test]
    fn test_reader() {
        let sequence = ["P0A7B8", "P0A7B8"];
        let source = ["UniProtKB", "UniProtKB"];
        let type_name = ["Initiator methionine", "Chain"];
        let starts = [1, 2];
        let ends = [1, 176];
        let scores = [None, Some(50)];
        let strand = [None, Some(Strand::Forward)];
        let frame = [".", "."];
        let attributes = ["Note=Removed", "ID=PRO_0000148105;Note=ATP-dependent protease subunit HslV"];

        let mut reader = Reader::new(GFF_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.sequence(), sequence[i]);
            assert_eq!(record.source(), source[i]);
            assert_eq!(record.type_name(), type_name[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.score(), scores[i]);
            assert_eq!(record.strand(), strand[i]);
            assert_eq!(record.frame(), frame[i]);
            assert_eq!(record.attributes(), attributes[i]);
        }
    }

    #[test]
    fn test_writer() {
        let mut reader = Reader::new(GFF_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records() {
            writer.write(r.ok().expect("Error reading record")).ok().expect("Error writing record");
        }
        assert_eq!(writer.inner.as_bytes(), GFF_FILE);
    }
}
