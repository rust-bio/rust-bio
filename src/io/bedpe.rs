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
//!
//! let example = b"1\t5\t500\t1\t600\t1000\tname1\t0.5";
//! let mut reader = bedpe::Reader::new(&example[..]);
//! let mut writer = bedpe::Writer::new(vec![]);
//! for record in reader.records() {
//!     let rec = record.expect("Error reading record.");
//!     assert_eq!("1", rec.chrom1());
//!     assert_eq!(5, rec.start1());
//!     writer.write(&rec).expect("Error writing record.");
//! }
//! ```

use std::marker::Copy;
use std::ops::Deref;

use bio_types::annot;
use bio_types::annot::loc::Loc;
use bio_types::strand;
use getset::{CloneGetters, CopyGetters, Getters, MutGetters, Setters, WithSetters};

use crate::io::common;

/// A BEDPE writer.
pub type Writer<W> = common::Writer<W>;

/// A BEDPE reader.
pub type Reader<R> = common::Reader<R, Record>;

/// BEDPE record as defined by bedtools
/// https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
#[derive(
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Default,
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
    /// Chromosome of first end of the feature.
    #[getset(get_mut = "pub", set_with = "pub")]
    chrom1: String,

    // Start position of first end of the feature (0-based).
    #[getset(get_copy = "pub", set = "pub", get_mut = "pub", set_with = "pub")]
    start1: u64,

    // End position of first end of the feature (0-based).
    #[getset(get_copy = "pub", set = "pub", get_mut = "pub", set_with = "pub")]
    end1: u64,

    /// Chromosome of second end of the feature.
    #[getset(get_mut = "pub", set_with = "pub")]
    chrom2: String,

    // Start position of second end of the feature (0-based).
    #[getset(get_copy = "pub", set = "pub", get_mut = "pub", set_with = "pub")]
    start2: u64,

    // End position of second end of the feature (0-based, not-included).
    #[getset(get_copy = "pub", set = "pub", get_mut = "pub", set_with = "pub")]
    end2: u64,

    #[serde(default)]
    aux: Vec<String>,
}

impl Record {
    /// Create new record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Chromosome of first end of the feature.
    pub fn chrom1(&self) -> &str {
        &self.chrom1
    }

    /// Chromosome of second end of the feature.
    pub fn chrom2(&self) -> &str {
        &self.chrom2
    }

    /// Access auxilliary fields after the second end field by index
    /// (counting first field (chromosome1) as 0).
    pub fn aux(&self, i: usize) -> Option<&str> {
        let j = i - 6;
        if j < self.aux.len() {
            Some(&self.aux[j])
        } else {
            None
        }
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
}

impl Record {
    /// Set chromosome of first end of the feature.
    pub fn set_chrom1(&mut self, chrom: &str) {
        self.chrom1 = chrom.into()
    }

    /// Set chromosome of second end of the feature.
    pub fn set_chrom2(&mut self, chrom: &str) {
        self.chrom2 = chrom.into()
    }

    /// Set auxiliary field after the second end field by index.
    /// If the index is bigger than the current auxiliary vector, empty strings will be pushed into
    /// the vector until the index is reached.
    /// (counting first field (chromosome1) as 0).
    pub fn set_aux(&mut self, i: usize, field: &str) {
        let j = i - 6;
        if j < self.aux.len() {
            self.aux[j] = field.into()
        } else {
            for _ in self.aux.len()..j {
                self.aux.push("".to_owned())
            }
            self.aux.push(field.into())
        }
    }

    /// Set name.
    pub fn set_name(&mut self, name: &str) {
        self.set_aux(6, name)
    }

    /// Set score.
    pub fn set_score(&mut self, score: &str) {
        self.set_aux(7, score)
    }

    /// Set strand of the first end of the feature.
    pub fn set_strand1(&mut self, strand1: strand::Strand) {
        self.set_aux(8, strand1.strand_symbol());
    }

    /// Set strand of the second end of the feature.
    pub fn set_strand2(&mut self, strand2: strand::Strand) {
        self.set_aux(9, strand2.strand_symbol());
    }

    /// Add auxilliary field.
    pub fn push_aux(&mut self, field: &str) {
        self.aux.push(field.into())
    }
}

/// Generate a BEDPE format `Record` from two annotation positions.
///
/// Each end of this record will have length 1, and when created it will have an
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
        bedpe.set_strand1(pos1.strand().into());
        bedpe.set_strand2(pos2.strand().into());
        bedpe
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bio_types::annot::pos::Pos;
    use bio_types::strand::{ReqStrand, Strand};
    use std::path::Path;

    const BEDPE_FILE: &[u8] = b"\
1\t5\t5000\t1\t6000\t10000\tname1\tup\t-\t+
2\t3\t5005\t2\t10\t5500\tname2\tup\t.\t.
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
        let strands = [(Some(Strand::Reverse), Some(Strand::Forward)), (None, None)];

        let mut reader = Reader::new(BEDPE_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!((record.chrom1(), record.chrom2()), chroms[i]);
            assert_eq!((record.start1(), record.start2()), starts[i]);
            assert_eq!((record.end1(), record.end2()), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);
            assert_eq!((record.strand1(), record.strand2()), strands[i]);
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

        assert_eq!(&error, "Failed to read from \"/I/dont/exist.bedpe\"")
    }

    #[test]
    fn test_writer() {
        let mut reader = Reader::new(BEDPE_FILE);
        let mut output: Vec<u8> = vec![];
        {
            let mut writer = Writer::new(&mut output);
            for r in reader.records() {
                writer
                    .write(r.expect("Error reading record"))
                    .expect("Error writing record");
            }
        }
        assert_eq!(&output[..], BEDPE_FILE);
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
