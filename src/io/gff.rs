// Copyright 2016 Pierre Marijon.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


//! GFF3 format reading and writing.
//!
//! GFF2 definition : http://gmod.org/wiki/GFF2#The_GFF2_File_Format (not yet support)
//! GTF2 definition : http://mblab.wustl.edu/GTF2.html (not yet support)
//! GFF3 definition : http://gmod.org/wiki/GFF3#GFF3_Format
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::gff;
//! let reader = gff::Reader::new(io::stdin(), gff::GffType::GFF3);
//! ```

use std::io;
use std::fs;
use std::path::Path;
use std::convert::AsRef;
use std::collections::HashMap;
use itertools::Itertools;

use csv;

use io::Strand;

/// GffType
///
/// We have three format in the GFF family.
/// The change is in the last field of GFF.
/// For each type we have key value separator and field separator
#[derive(Clone, Copy)]
pub enum GffType {
    /// Attribute format is key1=value, key2=value
    GFF3,
    /// Attribute format is key1 value; key2 value
    GFF2,
    /// Same as GFF2 just possible keyword and possible value change
    GTF2,
    /// Any, first field of tuple is key value separator, second is field separator
    Any(u8, u8),
}

impl GffType {
    #[inline]
    fn separator(&self) -> (u8, u8) {
        match *self {
            GffType::GFF3 => (b'=', b','),
            GffType::GFF2 => (b' ', b';'),
            GffType::GTF2 => (b' ', b';'),
            GffType::Any(x, y) => (x, y),
        }
    }
}

/// A GFF reader.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
    gff_type: GffType,
}

impl Reader<fs::File> {
    /// Read GFF from given file path in given format.
    pub fn from_file<P: AsRef<Path>>(path: P, fileformat: GffType) -> io::Result<Self> {
        fs::File::open(path).map(|f| Reader::new(f, fileformat))
    }
}


impl<R: io::Read> Reader<R> {
    /// Create a new GFF reader given an instance of `io::Read`, in given format.
    pub fn new(reader: R, fileformat: GffType) -> Self {
        Reader {
            inner: csv::Reader::from_reader(reader).delimiter(b'\t').has_headers(false),
            gff_type: fileformat,
        }
    }

    /// Iterate over all records.
    pub fn records(&mut self) -> Records<R> {
        Records {
            inner: self.inner.decode(),
            gff_type: self.gff_type
        }
    }
}

/// A GFF record.
pub struct Records<'a, R: 'a + io::Read> {
    inner: csv::DecodedRecords<'a, R, (String, String, String, u64, u64, String, String, String, String)>,
    gff_type: GffType,
}


impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        let (delim, termi) = self.gff_type.separator();

        self.inner.next().map(|res| {
            res.map(|(seqname, source, feature_type, start, end, score, strand, frame, attributes)| {
                Record {
                    seqname: seqname,
                    source: source,
                    feature_type: feature_type,
                    start: start,
                    end: end,
                    score: score,
                    strand: strand,
                    frame: frame,
                    attributes: csv::Reader::from_string(attributes)
                        .delimiter(delim)
                        .record_terminator(csv::RecordTerminator::Any(termi))
                        .has_headers(false)
                        .decode()
                        .collect::<csv::Result<_>>()
                        .unwrap(),
                }
            })
        })
    }
}


/// A GFF writer.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    delimiter: char,
    terminator: String,
}


impl Writer<fs::File> {
    /// Write to a given file path in given format.
    pub fn to_file<P: AsRef<Path>>(path: P, fileformat: GffType) -> io::Result<Self> {
        fs::File::create(path).map(|f| Writer::new(f, fileformat))
    }
}


impl<W: io::Write> Writer<W> {
    /// Write to a given writer.
    pub fn new(writer: W, fileformat: GffType) -> Self {
        let (delim, termi) = fileformat.separator();

        Writer {
            inner: csv::Writer::from_writer(writer).delimiter(b'\t').flexible(true),
            delimiter: delim as char,
            terminator: String::from_utf8(vec![termi]).unwrap(),
        }
    }

    /// Write a given GFF record.
    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        let attributes;

        if !record.attributes.is_empty() {
            attributes = record.attributes.iter().map(|(a, b)| format!("{}{}{}", a, self.delimiter, b)).join(&self.terminator);
        } else {
            attributes = "".to_owned();
        }

        self.inner.encode((record.seqname, record.source, record.feature_type, record.start, record.end, record.score, record.strand, record.frame, attributes))
    }
}


/// A GFF record
#[derive(RustcEncodable)]
pub struct Record {
    seqname: String,
    source: String,
    feature_type: String,
    start: u64,
    end: u64,
    score: String,
    strand: String,
    frame: String,
    attributes: HashMap<String, String>,
}

impl Record {
    /// Create a new GFF record.
    pub fn new() -> Self {
        Record {
            seqname: "".to_owned(),
            source: "".to_owned(),
            feature_type: "".to_owned(),
            start: 0,
            end: 0,
            score: ".".to_owned(),
            strand: ".".to_owned(),
            frame: "".to_owned(),
            attributes: HashMap::<String, String>::new(),
        }
    }

    /// Sequence name of the feature.
    pub fn seqname(&self) -> &str {
        &self.seqname
    }

    /// Source of the feature.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Type of the feature.
    pub fn feature_type(&self) -> &str {
        &self.feature_type
    }

    /// Start position of feature (1-based).
    pub fn start(&self) -> &u64 {
        &self.start
    }

    /// End position of feature (1-based, not included).
    pub fn end(&self) -> &u64 {
        &self.end
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
    pub fn attributes(&self) -> &HashMap<String, String> {
        &self.attributes
    }
    
    /// Get mutable reference on seqname of feature.
    pub fn seqname_mut(&mut self) -> &mut String {
        return &mut self.seqname;
    }

    /// Get mutable reference on source of feature.
    pub fn source_mut(&mut self) -> &mut String {
        return &mut self.source;
    }

    /// Get mutable reference on type of feature.
    pub fn feature_type_mut(&mut self) -> &mut String {
        return &mut self.feature_type;
    }

    /// Get mutable reference on start of feature.
    pub fn start_mut(&mut self) -> &mut u64 {
        return &mut self.start;
    }

    /// Get mutable reference on end of feature.
    pub fn end_mut(&mut self) -> &mut u64 {
        return &mut self.end;
    }

    /// Get mutable reference on score of feature.
    pub fn score_mut(&mut self) -> &mut String {
        return &mut self.score;
    }

    /// Get mutable reference on strand of feature.
    pub fn strand_mut(&mut self) -> &mut String {
        return &mut self.strand;
    }

    /// Get mutable reference on attributes of feature.
    pub fn attributes_mut(&mut self) -> &mut HashMap<String, String> {
        return &mut self.attributes;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use io::Strand;
    use std::collections::HashMap;
    
    const GFF_FILE: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote=Removed,ID=test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote=ATP-dependent protease subunit HslV,ID=PRO_0000148105
";
    //required because HashMap iter on element randomly
    const GFF_FILE_ONE_ATTRIB: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote=Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID=PRO_0000148105
";

    const GTF_FILE: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed;ID test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote ATP-dependent;ID PRO_0000148105
";

    //required because HashMap iter on element randomly
    const GTF_FILE_ONE_ATTRIB: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID PRO_0000148105
";

    #[test]
    fn test_reader_gff3() {
        let seqname = ["P0A7B8", "P0A7B8"];
        let source = ["UniProtKB", "UniProtKB"];
        let feature_type = ["Initiator methionine", "Chain"];
        let starts = [1, 2];
        let ends = [1, 176];
        let scores = [None, Some(50)];
        let strand = [None, Some(Strand::Forward)];
        let frame = [".", "."];
        let mut attributes = [HashMap::new(), HashMap::new()];
        attributes[0].insert("ID".to_owned(), "test".to_owned());
        attributes[0].insert("Note".to_owned(), "Removed".to_owned());
        attributes[1].insert("ID".to_owned(), "PRO_0000148105".to_owned());
        attributes[1].insert("Note".to_owned(), "ATP-dependent protease subunit HslV".to_owned());

        let mut reader = Reader::new(GFF_FILE, GffType::GFF3);
        for (i, r) in reader.records().enumerate() {
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.seqname(), seqname[i]);
            assert_eq!(record.source(), source[i]);
            assert_eq!(record.feature_type(), feature_type[i]);
            assert_eq!(*record.start(), starts[i]);
            assert_eq!(*record.end(), ends[i]);
            assert_eq!(record.score(), scores[i]);
            assert_eq!(record.strand(), strand[i]);
            assert_eq!(record.frame(), frame[i]);
            assert_eq!(record.attributes(), &attributes[i]);
        }
    }

    #[test]
    fn test_reader_gtf2() {
        let seqname = ["P0A7B8", "P0A7B8"];
        let source = ["UniProtKB", "UniProtKB"];
        let feature_type = ["Initiator methionine", "Chain"];
        let starts = [1, 2];
        let ends = [1, 176];
        let scores = [None, Some(50)];
        let strand = [None, Some(Strand::Forward)];
        let frame = [".", "."];
        let mut attributes = [HashMap::new(), HashMap::new()];
        attributes[0].insert("ID".to_owned(), "test".to_owned());
        attributes[0].insert("Note".to_owned(), "Removed".to_owned());
        attributes[1].insert("ID".to_owned(), "PRO_0000148105".to_owned());
        attributes[1].insert("Note".to_owned(), "ATP-dependent".to_owned());

        let mut reader = Reader::new(GTF_FILE, GffType::GTF2);
        for (i, r) in reader.records().enumerate() {
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.seqname(), seqname[i]);
            assert_eq!(record.source(), source[i]);
            assert_eq!(record.feature_type(), feature_type[i]);
            assert_eq!(*record.start(), starts[i]);
            assert_eq!(*record.end(), ends[i]);
            assert_eq!(record.score(), scores[i]);
            assert_eq!(record.strand(), strand[i]);
            assert_eq!(record.frame(), frame[i]);
            assert_eq!(record.attributes(), &attributes[i]);
        }
    }

    #[test]
    fn test_writer_gff3() {
        let mut reader = Reader::new(GFF_FILE_ONE_ATTRIB, GffType::GFF3);
        let mut writer = Writer::new(vec![], GffType::GFF3);
        for r in reader.records() {
            writer.write(r.ok().expect("Error reading record")).ok().expect("Error writing record");
        }
        assert_eq!(writer.inner.as_string(), String::from_utf8_lossy(GFF_FILE_ONE_ATTRIB))
    }

    #[test]
    fn test_writer_gtf2() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GffType::GTF2);
        let mut writer = Writer::new(vec![], GffType::GTF2);
        for r in reader.records() {
            writer.write(r.ok().expect("Error reading record")).ok().expect("Error writing record");
        }
        assert_eq!(writer.inner.as_string(), String::from_utf8_lossy(GTF_FILE_ONE_ATTRIB))
    }

    #[test]
    fn test_convert_gtf2_to_gff3() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GffType::GTF2);
        let mut writer = Writer::new(vec![], GffType::GFF3);
        for r in reader.records() {
            writer.write(r.ok().expect("Error reading record")).ok().expect("Error writing record");
        }
        assert_eq!(writer.inner.as_string(), String::from_utf8_lossy(GFF_FILE_ONE_ATTRIB))
    }
}
