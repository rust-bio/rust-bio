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
use std::fmt::Write;
use std::fs;
use std::io;
use std::marker::Copy;
use std::ops::Deref;
use std::path::Path;

use csv;

use bio_types::annot;
use bio_types::annot::loc::Loc;
use bio_types::strand;

/// A BED reader.
#[derive(Debug)]
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
    pub fn records(&mut self) -> Records<'_, R> {
        Records {
            inner: self.inner.deserialize(),
        }
    }
}

type BedRecordCsv = (String, u64, u64, Option<Vec<String>>);

/// A BED record.
pub struct Records<'a, R: io::Read> {
    inner: csv::DeserializeRecordsIter<'a, R, BedRecordCsv>,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner.next().map(|res| {
            res.map(|(chrom, start, end, aux)| Record {
                chrom,
                start,
                end,
                aux: aux.unwrap_or_else(Vec::new),
            })
        })
    }
}

/// A BED writer.
#[derive(Debug)]
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
#[derive(Debug, Serialize, Default, Deserialize, Clone)]
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
    pub fn strand(&self) -> Option<strand::Strand> {
        match self.aux(5) {
            Some("+") => Some(strand::Strand::Forward),
            Some("-") => Some(strand::Strand::Reverse),
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

impl<'a> From<&'a Record> for annot::contig::Contig<String, strand::Strand> {
    /// Returns a `Contig` annotation for the BED record.
    ///
    /// ```
    /// # extern crate bio;
    /// # extern crate bio_types;
    /// use bio::io::bed;
    /// use bio_types::strand::Strand;
    /// use bio_types::annot::contig::Contig;
    /// # use std::error::Error;
    /// # fn try_main() -> Result<(), Box<Error>> {
    /// let example = b"chr1\t5\t5000\tname1\t0.5";
    /// let mut reader = bed::Reader::new(&example[..]);
    /// let rec = reader.records().next().ok_or("No record available!")??;
    /// let loc = Contig::from(&rec);
    /// assert_eq!(loc.to_string(), "chr1:5-5000");
    /// # Ok(())
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    fn from(rec: &Record) -> Self {
        annot::contig::Contig::new(
            rec.chrom.to_string(),
            rec.start as isize,
            (rec.end - rec.start) as usize,
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
        bed.push_aux(pos.strand().into().strand_symbol());
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
        bed.push_aux(contig.strand().into().strand_symbol());
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
/// use bio_types::strand::ReqStrand;
/// use bio_types::annot::AnnotError;
/// use bio_types::annot::spliced::{Spliced,SplicingError};
/// use bio::io::bed;
/// # fn try_main() -> Result<(), Box<SplicingError>> {
/// let tad3 = Spliced::with_lengths_starts("chrXII".to_owned(), 765265,
///                                         &vec![808,52,109], &vec![0,864,984],
///                                         ReqStrand::Reverse)?;
/// assert_eq!(tad3.to_string(), "chrXII:765265-766073;766129-766181;766249-766358(-)");
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
///     writer.write(&tad3_bed).ok().unwrap();
/// }
/// assert_eq!("chrXII\t765265\t766358\tYLR316C\t0\t-\t765265\t766358\t0\t3\t808,52,109,\t0,864,984,\n",
///            String::from_utf8(buf).unwrap_or_else(|_| "???".to_owned()).as_str());
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
        bed.push_aux(spliced.strand().into().strand_symbol());
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

    use bio_types::annot::spliced::Spliced;
    use bio_types::strand::ReqStrand;

    const BED_FILE: &'static [u8] = b"1\t5\t5000\tname1\tup
2\t3\t5005\tname2\tup
";
    //const BED_FILE_COMPACT: &'static [u8] = b"1\t5\t5000\n2\t3\t5005\n";

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

    #[test]
    fn spliced_to_bed() {
        //chrV    166236  166885  YER007C-A       0       -       166236  166885  0       2       535,11, 0,638,
        let tma20 = Spliced::with_lengths_starts(
            "chrV".to_owned(),
            166236,
            &vec![535, 11],
            &vec![0, 638],
            ReqStrand::Reverse,
        )
        .unwrap();
        let mut buf = Vec::new();
        {
            let mut writer = Writer::new(&mut buf);
            let mut tma20_bed = Record::from(tma20);
            tma20_bed.set_name("YER007C-A");
            writer.write(&tma20_bed).ok().unwrap();
        }
        assert_eq!(
            "chrV\t166236\t166885\tYER007C-A\t0\t-\t166236\t166885\t0\t2\t535,11,\t0,638,\n",
            String::from_utf8(buf).unwrap().as_str()
        );

        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &vec![11, 94, 630],
            &vec![0, 420, 921],
            ReqStrand::Forward,
        )
        .unwrap();
        let mut buf = Vec::new();
        {
            let mut writer = Writer::new(&mut buf);
            let mut rpl7b_bed = Record::from(rpl7b);
            rpl7b_bed.set_name("YPL198W");
            writer.write(&rpl7b_bed).ok().unwrap();
        }
        assert_eq!(
            "chrXVI\t173151\t174702\tYPL198W\t0\t+\t173151\t174702\t0\t3\t11,94,630,\t0,420,921,\n",
            String::from_utf8(buf).unwrap().as_str()
        );

        //chrXII  765265  766358  YLR316C 0       -       765265  766358  0       3       808,52,109,     0,864,984,
        let tad3 = Spliced::with_lengths_starts(
            "chrXII".to_owned(),
            765265,
            &vec![808, 52, 109],
            &vec![0, 864, 984],
            ReqStrand::Reverse,
        )
        .unwrap();
        let mut buf = Vec::new();
        {
            let mut writer = Writer::new(&mut buf);
            let mut tad3_bed = Record::from(tad3);
            tad3_bed.set_name("YLR316C");
            writer.write(&tad3_bed).ok().unwrap();
        }
        assert_eq!("chrXII\t765265\t766358\tYLR316C\t0\t-\t765265\t766358\t0\t3\t808,52,109,\t0,864,984,\n",
                   String::from_utf8(buf).unwrap().as_str());
    }
}
