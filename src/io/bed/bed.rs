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
//!
//! let example = b"1\t5\t5000\tname1\t0.5";
//! let mut reader = bed::Reader::new(&example[..]);
//! let mut writer = bed::Writer::new(vec![]);
//! for record in reader.records() {
//!     let rec = record.expect("Error reading record.");
//!     println!("{}", rec.chrom());
//!     writer.write(&rec).expect("Error writing record.");
//! }
//! ```

use std::fmt::Write;
use std::marker::Copy;
use std::ops::Deref;

use bio_types::annot;
use bio_types::annot::loc::Loc;
use bio_types::strand;
use derefable::Derefable;

use crate::io::common;

/// BED writer.
pub type Writer<W> = common::Writer<W>;

/// BED reader.
pub type Reader<R> = common::Reader<R, Record>;

/// A BED record as defined by BEDtools
/// (http://bedtools.readthedocs.org/en/latest/content/general-usage.html)
#[derive(
    Debug, Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize, Derefable,
)]
pub struct Record(#[deref(mutable)] common::Record);

/// Add new methods specific for a BED Record
impl Record {
    /// Create a new Record
    pub fn new() -> Self {
        Self::default()
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
    pub fn strand(&self) -> Option<strand::Strand> {
        match self.aux(5) {
            Some("+") => Some(strand::Strand::Forward),
            Some("-") => Some(strand::Strand::Reverse),
            _ => None,
        }
    }

    /// Set name.
    pub fn set_name(&mut self, name: &str) {
        self.set_aux(3, name)
    }

    /// Set score.
    pub fn set_score(&mut self, score: &str) {
        self.set_aux(4, score)
    }

    /// Set strand.
    pub fn set_strand(&mut self, strand: strand::Strand) {
        self.set_aux(5, strand.strand_symbol());
    }
}

impl From<&Record> for annot::contig::Contig<String, strand::Strand> {
    /// Returns a `Contig` annotation for the BED record.
    ///
    /// ```
    /// use bio::io::bed;
    /// use bio_types::annot::contig::Contig;
    /// use bio_types::strand::Strand;
    /// let example = b"chr1\t5\t5000\tname1\t0.5";
    /// let mut reader = bed::Reader::new(&example[..]);
    /// let rec = reader
    ///     .records()
    ///     .next()
    ///     .expect("Found no bed record.")
    ///     .expect("Got a csv::Error");
    /// let loc = Contig::from(&rec);
    /// assert_eq!(loc.to_string(), "chr1:5-5000");
    /// ```
    fn from(rec: &Record) -> Self {
        annot::contig::Contig::new(
            rec.chrom().to_string(),
            rec.start() as isize,
            (rec.end() - rec.start()) as usize,
            rec.strand().unwrap_or(strand::Strand::Unknown),
        )
    }
}

/// Generate a BED format `Record` for an annotation position.
///
/// This record will have length 1, and when created it will have an
/// empty name.
impl<R, S> From<annot::pos::Pos<R, S>> for Record
where
    R: Deref<Target = str>,
    S: Into<strand::Strand> + Copy,
{
    fn from(pos: annot::pos::Pos<R, S>) -> Self {
        let mut bed = Record::new();
        bed.set_chrom(pos.refid());
        bed.set_start(pos.pos() as u64);
        bed.set_end((pos.pos() + 1) as u64);
        bed.set_name("");
        bed.set_score("0");
        bed.set_strand(pos.strand().into());
        bed
    }
}

/// Generate a BED format `Record` for the location.
///
/// As created, it will have an empty name.
impl<R, S> From<annot::contig::Contig<R, S>> for Record
where
    R: Deref<Target = str>,
    S: Into<strand::Strand> + Copy,
{
    fn from(contig: annot::contig::Contig<R, S>) -> Self {
        let mut bed = Record::new();
        bed.set_chrom(contig.refid());
        bed.set_start(contig.start() as u64);
        bed.set_end((contig.start() + contig.length() as isize) as u64);
        bed.set_name("");
        bed.set_score("0");
        bed.set_strand(contig.strand().into());
        bed
    }
}

/// Generate a BED format `Record` for the position.
///
/// Splicing information will be represented with the 12-column
/// BED format, using columns 10 through 12 (blockCount,
/// blockSizes, and blockStarts) for exons.
///
/// As created, it will have an empty name and default to using
/// the overall start & end (columns 1 and 2) for the start and
/// end of the "thick" region (columns 7 and 8).
/// ```
/// # extern crate bio;
/// # extern crate bio_types;
/// use bio::io::bed;
/// use bio_types::annot::spliced::{Spliced, SplicingError};
/// use bio_types::annot::AnnotError;
/// use bio_types::strand::ReqStrand;
/// # fn try_main() -> Result<(), Box<SplicingError>> {
/// let tad3 = Spliced::with_lengths_starts(
///     "chrXII".to_owned(),
///     765265,
///     &vec![808, 52, 109],
///     &vec![0, 864, 984],
///     ReqStrand::Reverse,
/// )
/// .expect("Encountered a bio_types::annot::spliced::SplicingError.");
/// assert_eq!(
///     tad3.to_string(),
///     "chrXII:765265-766073;766129-766181;766249-766358(-)"
/// );
/// let tad3_exons = tad3.exon_contigs();
/// assert_eq!(tad3_exons.len(), 3);
/// assert_eq!(tad3_exons[0].to_string(), "chrXII:766249-766358(-)");
/// assert_eq!(tad3_exons[1].to_string(), "chrXII:766129-766181(-)");
/// assert_eq!(tad3_exons[2].to_string(), "chrXII:765265-766073(-)");
/// let mut buf = Vec::new();
/// {
///     let mut writer = bed::Writer::new(&mut buf);
///     let mut tad3_bed = bed::Record::from(tad3);
///     tad3_bed.set_name("YLR316C");
///     writer.write(&tad3_bed).unwrap();
/// }
/// assert_eq!(
///     "chrXII\t765265\t766358\tYLR316C\t0\t-\t765265\t766358\t0\t3\t808,52,109,\t0,864,984,\n",
///     String::from_utf8(buf)
///         .unwrap_or_else(|_| "???".to_owned())
///         .as_str()
/// );
/// # Ok(())
/// # }
/// # fn main() { try_main().unwrap(); }
/// ```
impl<R, S> From<annot::spliced::Spliced<R, S>> for Record
where
    R: Deref<Target = str>,
    S: Into<strand::Strand> + Copy,
{
    fn from(spliced: annot::spliced::Spliced<R, S>) -> Self {
        let mut bed = Record::new();
        bed.set_chrom(spliced.refid());
        bed.set_start(spliced.start() as u64);
        bed.set_end((spliced.start() + spliced.length() as isize) as u64);
        bed.set_name("");
        bed.set_score("0");
        bed.set_strand(spliced.strand().into());
        bed.push_aux(spliced.start().to_string().as_str()); // thickStart = chromStart
        bed.push_aux(
            (spliced.start() + spliced.length() as isize)
                .to_string()
                .as_str(),
        ); // thickEnd = chromEnd
        bed.push_aux("0"); // RGB color = black

        bed.push_aux(spliced.exon_count().to_string().as_str());

        let mut block_sizes = String::new();
        for block_size in spliced.exon_lengths() {
            write!(block_sizes, "{},", block_size).unwrap();
        }
        bed.push_aux(&block_sizes);

        let mut block_starts = String::new();
        for block_start in spliced.exon_starts() {
            write!(block_starts, "{},", block_start).unwrap();
        }
        bed.push_aux(&block_starts);
        bed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    use bio_types::annot::{contig::Contig, pos::Pos, spliced::Spliced};
    use bio_types::strand::{ReqStrand, Strand};

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

    #[test]
    fn test_reader() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];

        let mut reader = Reader::new(BED_FILE);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);
        }
    }

    #[test]
    fn test_reader_with_comment() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];

        let mut reader = Reader::new(BED_FILE_COMMENT);
        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);
        }
    }

    #[test]
    fn test_reader_compact() {
        let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];

        let mut reader = Reader::new(BED_FILE_COMPACT);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
        }
    }

    #[test]
    fn test_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.bed");
        let error = Reader::from_file(path)
            .unwrap_err()
            .downcast::<String>()
            .unwrap();

        assert_eq!(&error, "Failed to read from \"/I/dont/exist.bed\"")
    }

    #[test]
    fn test_writer() {
        let mut reader = Reader::new(BED_FILE);
        let mut output: Vec<u8> = vec![];
        {
            let mut writer = Writer::new(&mut output);
            for r in reader.records() {
                writer
                    .write(r.expect("Error reading record"))
                    .expect("Error writing record");
            }
        }
        assert_eq!(output, BED_FILE);
    }

    #[test]
    fn spliced_to_bed() {
        //chrV    166236  166885  YER007C-A       0       -       166236  166885  0       2       535,11, 0,638,
        let tma20 = Spliced::with_lengths_starts(
            "chrV".to_owned(),
            166236,
            &[535, 11],
            &[0, 638],
            ReqStrand::Reverse,
        )
        .unwrap();
        let mut buf = Vec::new();
        {
            let mut writer = Writer::new(&mut buf);
            let mut tma20_bed = Record::from(tma20);
            tma20_bed.set_name("YER007C-A");
            writer.write(&tma20_bed).unwrap();
        }
        assert_eq!(
            "chrV\t166236\t166885\tYER007C-A\t0\t-\t166236\t166885\t0\t2\t535,11,\t0,638,\n",
            String::from_utf8(buf).unwrap().as_str()
        );

        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &[11, 94, 630],
            &[0, 420, 921],
            ReqStrand::Forward,
        )
        .unwrap();
        let mut buf = Vec::new();
        {
            let mut writer = Writer::new(&mut buf);
            let mut rpl7b_bed = Record::from(rpl7b);
            rpl7b_bed.set_name("YPL198W");
            writer.write(&rpl7b_bed).unwrap();
        }
        assert_eq!(
            "chrXVI\t173151\t174702\tYPL198W\t0\t+\t173151\t174702\t0\t3\t11,94,630,\t0,420,921,\n",
            String::from_utf8(buf).unwrap().as_str()
        );

        //chrXII  765265  766358  YLR316C 0       -       765265  766358  0       3       808,52,109,     0,864,984,
        let tad3 = Spliced::with_lengths_starts(
            "chrXII".to_owned(),
            765265,
            &[808, 52, 109],
            &[0, 864, 984],
            ReqStrand::Reverse,
        )
        .unwrap();
        let mut buf = Vec::new();
        {
            let mut writer = Writer::new(&mut buf);
            let mut tad3_bed = Record::from(tad3);
            tad3_bed.set_name("YLR316C");
            writer.write(&tad3_bed).unwrap();
        }
        assert_eq!("chrXII\t765265\t766358\tYLR316C\t0\t-\t765265\t766358\t0\t3\t808,52,109,\t0,864,984,\n",
               String::from_utf8(buf).unwrap().as_str());
    }

    #[test]
    fn test_bed_from_contig() {
        let contig = Contig::new(
            "chrXI".to_owned(),
            334412,
            334916 - 334412,
            ReqStrand::Reverse,
        );

        let record = Record::from(contig);

        assert_eq!(record.chrom(), String::from("chrXI"));
        assert_eq!(record.start(), 334412);
        assert_eq!(record.end(), 334412 + (334916 - 334412));
        assert_eq!(record.name(), Some(""));
        assert_eq!(record.score(), Some("0"));
        assert_eq!(record.strand(), Some(Strand::Reverse));
    }

    #[test]
    fn test_bed_from_pos() {
        let pos = Pos::new("chrXI".to_owned(), 334412, ReqStrand::Reverse);

        let record = Record::from(pos);

        assert_eq!(record.chrom(), String::from("chrXI"));
        assert_eq!(record.start(), 334412);
        assert_eq!(record.end(), 334412 + 1);
        assert_eq!(record.name(), Some(""));
        assert_eq!(record.score(), Some("0"));
        assert_eq!(record.strand(), Some(Strand::Reverse));
    }
}
