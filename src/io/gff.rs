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
//! // import functions (at top of script)
//! use bio::io::gff;
//! use std::io;
//! let mut reader = gff::Reader::new(io::stdin(), gff::GffType::GFF3);
//! let mut writer = gff::Writer::new(vec![], gff::GffType::GFF3);
//! for record in reader.records() {
//!     let rec = record.ok().expect("Error reading record.");
//!     println!("{}", rec.seqname());
//!     writer.write(&rec).ok().expect("Error writing record.");
//! }
//! ```

use itertools::Itertools;
use multimap::MultiMap;
use regex::Regex;
use std::convert::AsRef;
use std::fs;
use std::io;
use std::path::Path;
use std::str::FromStr;

use bio_types::strand::Strand;

/// `GffType`
///
/// We have three format in the GFF family.
/// The change is in the last field of GFF.
/// For each type we have key value separator and field separator
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GffType {
    /// Attribute format is: key1=value; key2=value1,value2
    GFF3,
    /// Attribute format is: key1 value; key2 value1; key2 value2
    GFF2,
    /// Same as GFF2 just possible keyword and possible value change
    GTF2,
    /// Any, first field of tuple separates key from value,
    /// second field separates multiple key value pairs, and
    /// third field separates multiple values for the same key
    Any(u8, u8, u8),
}

impl FromStr for GffType {
    type Err = String;

    /// Create a GffType from a string.
    ///
    /// # Arguments
    ///
    /// * `src_str` - The source string to convert to the GffType.
    fn from_str(src_str: &str) -> Result<Self, Self::Err> {
        match src_str {
            "gff3" => Ok(GffType::GFF3),
            "gff2" => Ok(GffType::GFF2),
            "gtf2" => Ok(GffType::GTF2),
            _ => Err(format!(
                "String '{}' is not a valid GFFType (GFF/GTF format version).",
                src_str
            )),
        }
    }
}

impl GffType {
    #[inline]
    /// First field is key value separator.
    /// Second field terminates a key value pair.
    /// Third field
    fn separator(self) -> (u8, u8, u8) {
        match self {
            GffType::GFF3 => (b'=', b';', b','),
            GffType::GFF2 => (b' ', b';', 0u8),
            GffType::GTF2 => (b' ', b';', 0u8),
            GffType::Any(x, y, z) => (x, y, z),
        }
    }
}

/// A GFF reader.
#[derive(Debug)]
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
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_reader(reader),
            gff_type: fileformat,
        }
    }

    /// Iterate over all records.
    pub fn records(&mut self) -> Records<'_, R> {
        let (delim, term, vdelim) = self.gff_type.separator();
        let r = format!(
            r" *(?P<key>[^{delim}{term}\t]+){delim}(?P<value>[^{delim}{term}\t]+){term}?",
            delim = delim as char,
            term = term as char
        );
        let attribute_re = Regex::new(&r).unwrap();
        Records {
            inner: self.inner.deserialize(),
            attribute_re,
            value_delim: vdelim as char,
        }
    }
}

type GffRecordInner = (
    String,
    String,
    String,
    u64,
    u64,
    String,
    String,
    String,
    String,
);

/// An iterator over the records of a GFF file.
pub struct Records<'a, R: io::Read> {
    inner: csv::DeserializeRecordsIter<'a, R, GffRecordInner>,
    attribute_re: Regex,
    value_delim: char,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner.next().map(|res| {
            res.map(
                |(
                    seqname,
                    source,
                    feature_type,
                    start,
                    end,
                    score,
                    strand,
                    frame,
                    raw_attributes,
                )| {
                    let trim_quotes = |s: &str| s.trim_matches('\'').trim_matches('"').to_owned();
                    let mut attributes = MultiMap::new();
                    for caps in self.attribute_re.captures_iter(&raw_attributes) {
                        for value in caps["value"].split(self.value_delim) {
                            attributes.insert(trim_quotes(&caps["key"]), trim_quotes(value));
                        }
                    }
                    Record {
                        seqname,
                        source,
                        feature_type,
                        start,
                        end,
                        score,
                        strand,
                        frame,
                        attributes,
                    }
                },
            )
        })
    }
}

/// A GFF writer.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    delimiter: char,
    terminator: String,
}

impl Writer<fs::File> {
    /// Write to a given file path in given format.
    #[allow(clippy::wrong_self_convention)]
    pub fn to_file<P: AsRef<Path>>(path: P, fileformat: GffType) -> io::Result<Self> {
        fs::File::create(path).map(|f| Writer::new(f, fileformat))
    }
}

impl<W: io::Write> Writer<W> {
    /// Write to a given writer.
    pub fn new(writer: W, fileformat: GffType) -> Self {
        let (delim, termi, _) = fileformat.separator();

        Writer {
            inner: csv::WriterBuilder::new()
                .delimiter(b'\t')
                .flexible(true)
                .from_writer(writer),
            delimiter: delim as char,
            terminator: String::from_utf8(vec![termi]).unwrap(),
        }
    }

    /// Write a given GFF record.
    pub fn write(&mut self, record: &Record) -> csv::Result<()> {
        let attributes = if !record.attributes.is_empty() {
            record
                .attributes
                .iter()
                .map(|(a, b)| format!("{}{}{}", a, self.delimiter, b))
                .join(&self.terminator)
        } else {
            "".to_owned()
        };

        self.inner.serialize((
            &record.seqname,
            &record.source,
            &record.feature_type,
            record.start,
            record.end,
            &record.score,
            &record.strand,
            &record.frame,
            attributes,
        ))
    }
}

/// A GFF record
///
#[derive(
    Debug,
    Default,
    Serialize,
    Deserialize,
    Clone,
    PartialEq,
    Eq,
    Getters,
    Setters,
    MutGetters,
    Builder,
)]
#[getset(get_mut = "pub", set = "pub")]
#[builder(default)]
pub struct Record {
    #[getset(get = "pub")]
    seqname: String,
    #[getset(get = "pub")]
    source: String,
    #[getset(get = "pub")]
    feature_type: String,
    #[getset(get = "pub")]
    start: u64,
    #[getset(get = "pub")]
    end: u64,
    #[builder(default = "\".\".to_owned()")]
    score: String,
    #[builder(default = "\".\".to_owned()")]
    strand: String,
    #[getset(get = "pub")]
    frame: String,
    #[getset(get = "pub")]
    attributes: MultiMap<String, String>,
}

impl Record {
    /// Create a new GFF record.
    pub fn new() -> Self {
        RecordBuilder::default().build().unwrap()
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use multimap::multimap;

    const GFF_FILE: &[u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\t\
Note=Removed,Obsolete;ID=test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote=ATP-dependent protease subunit HslV;\
ID=PRO_0000148105";
    const GFF_FILE_WITH_COMMENT: &[u8] = b"#comment
P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\t\
Note=Removed,Obsolete;ID=test
#comment
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote=ATP-dependent protease subunit HslV;\
ID=PRO_0000148105";
    //required because MultiMap iter on element randomly
    const GFF_FILE_ONE_ATTRIB: &[u8] =
        b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote=Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID=PRO_0000148105
";

    const GTF_FILE: &[u8] =
        b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed;ID test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote ATP-dependent;ID PRO_0000148105
";

    // Another variant of GTF file, modified from a published GENCODE GTF file.
    const GTF_FILE_2: &[u8] = b"chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\t\
gene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\";\
transcript_id \"ENST00000456328.2\"; gene_type \"transcribed_unprocessed_pseudogene\"";

    // GTF file with duplicate attribute keys, taken from a published GENCODE GTF file.
    const GTF_FILE_DUP_ATTR_KEYS: &[u8] = b"chr1\tENSEMBL\ttranscript\t182393\t\
184158\t.\t+\t.\tgene_id \"ENSG00000279928.1\"; transcript_id \"ENST00000624431.1\";\
gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"FO538757.2\";\
transcript_type \"protein_coding\"; transcript_status \"KNOWN\";\
transcript_name \"FO538757.2-201\"; level 3; protein_id \"ENSP00000485457.1\";\
transcript_support_level \"1\"; tag \"basic\"; tag \"appris_principal_1\";";

    //required because MultiMap iter on element randomly
    const GTF_FILE_ONE_ATTRIB: &[u8] =
        b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID PRO_0000148105
";

    #[test]
    fn test_default() {
        let gff = RecordBuilder::default()
            .seqname("".to_owned())
            .source("".to_owned())
            .feature_type("".to_owned())
            .start(0)
            .end(0)
            .score(".".to_owned())
            .strand(".".to_owned())
            .frame("".to_owned())
            .attributes(multimap!())
            .build()
            .unwrap();

        let gff_def = Record::new();
        assert_eq!(gff, gff_def);
    }

    #[test]
    fn test_reader_gff3() {
        let gffs = vec![
            RecordBuilder::default()
                .seqname("P0A7B8".to_owned())
                .source("UniProtKB".to_owned())
                .feature_type("Initiator methionine".to_owned())
                .start(1)
                .end(1)
                .score(".".to_owned())
                .strand(".".to_owned())
                .frame(".".to_owned())
                .attributes(multimap!(
                    "ID".to_owned() => "test".to_owned(),
                    "Note".to_owned() => "Removed".to_owned(),
                    "Note".to_owned() => "Obsolete".to_owned()
                ))
                .build()
                .unwrap(),
            RecordBuilder::default()
                .seqname("P0A7B8".to_owned())
                .source("UniProtKB".to_owned())
                .feature_type("Chain".to_owned())
                .start(2)
                .end(176)
                .score("50".to_owned())
                .strand("+".to_owned())
                .frame(".".to_owned())
                .attributes(multimap!(
                    "ID".to_owned() => "PRO_0000148105".to_owned(),
                    "Note".to_owned() => "ATP-dependent protease subunit HslV".to_owned()
                ))
                .build()
                .unwrap(),
        ];

        let mut reader = Reader::new(GFF_FILE, GffType::GFF3);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record, gffs[i]);
        }

        let mut reader = Reader::new(GFF_FILE_WITH_COMMENT, GffType::GFF3);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record, gffs[i]);
        }
    }

    #[test]
    fn test_gff_type_from_str() {
        let gff3 = GffType::from_str("gff3").expect("Error parsing");
        assert_eq!(gff3, GffType::GFF3);

        let gff2 = GffType::from_str("gff2").expect("Error parsing");
        assert_eq!(gff2, GffType::GFF2);

        let gtf2 = GffType::from_str("gtf2").expect("Error parsing");
        assert_eq!(gtf2, GffType::GTF2);

        let unk = GffType::from_str("unknown").unwrap_err();
        assert_eq!(
            unk,
            "String 'unknown' is not a valid GFFType (GFF/GTF format version)."
        )
    }

    #[test]
    fn test_reader_gtf2() {
        let gffs = vec![
            RecordBuilder::default()
                .seqname("P0A7B8".to_owned())
                .source("UniProtKB".to_owned())
                .feature_type("Initiator methionine".to_owned())
                .start(1)
                .end(1)
                .score(".".to_owned())
                .strand(".".to_owned())
                .frame(".".to_owned())
                .attributes(multimap!(
                    "ID".to_owned() => "test".to_owned(),
                    "Note".to_owned() => "Removed".to_owned()
                ))
                .build()
                .unwrap(),
            RecordBuilder::default()
                .seqname("P0A7B8".to_owned())
                .source("UniProtKB".to_owned())
                .feature_type("Chain".to_owned())
                .start(2)
                .end(176)
                .score("50".to_owned())
                .strand("+".to_owned())
                .frame(".".to_owned())
                .attributes(multimap!(
                    "ID".to_owned() => "PRO_0000148105".to_owned(),
                    "Note".to_owned() => "ATP-dependent".to_owned()
                ))
                .build()
                .unwrap(),
        ];

        let mut reader = Reader::new(GTF_FILE, GffType::GTF2);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record, gffs[i]);
        }
    }

    #[test]
    fn test_reader_gtf2_2() {
        let gffs = vec![
            RecordBuilder::default()
                .seqname("chr1".to_owned())
                .source("HAVANA".to_owned())
                .feature_type("gene".to_owned())
                .start(11869)
                .end(14409)
                .score(".".to_owned())
                .strand("+".to_owned())
                .frame(".".to_owned())
                .attributes(multimap!(
                    "gene_id".to_owned() => "ENSG00000223972.5".to_owned(),
                    "gene_type".to_owned() => "transcribed_unprocessed_pseudogene".to_owned()
                ))
                .build()
                .unwrap(),
            RecordBuilder::default()
                .seqname("chr1".to_owned())
                .source("HAVANA".to_owned())
                .feature_type("transcript".to_owned())
                .start(11869)
                .end(14409)
                .score(".".to_owned())
                .strand("+".to_owned())
                .frame(".".to_owned())
                .attributes(multimap!(
                    "gene_id".to_owned() => "ENSG00000223972.5".to_owned(),
                    "transcript_id".to_owned() => "ENST00000456328.2".to_owned(),
                    "gene_type".to_owned() => "transcribed_unprocessed_pseudogene".to_owned()
                ))
                .build()
                .unwrap(),
        ];

        let mut reader = Reader::new(GTF_FILE_2, GffType::GTF2);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record, gffs[i]);
        }
    }

    #[test]
    fn test_reader_gtf2_dup_attr_keys() {
        let mut reader = Reader::new(GTF_FILE_DUP_ATTR_KEYS, GffType::GTF2);
        let mut records = reader.records().collect::<Vec<_>>();
        assert_eq!(records.len(), 1);
        let record = records.pop().unwrap().expect("expected one record");
        assert_eq!(record.attributes.get("tag"), Some(&"basic".to_owned()));
        assert_eq!(
            record.attributes.get_vec("tag"),
            Some(&vec!["basic".to_owned(), "appris_principal_1".to_owned()])
        );
    }

    #[test]
    fn test_writer_gff3() {
        let mut reader = Reader::new(GFF_FILE_ONE_ATTRIB, GffType::GFF3);
        let mut writer = Writer::new(vec![], GffType::GFF3);
        for r in reader.records() {
            writer
                .write(&r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), GFF_FILE_ONE_ATTRIB)
    }

    #[test]
    fn test_writer_gtf2() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GffType::GTF2);
        let mut writer = Writer::new(vec![], GffType::GTF2);
        for r in reader.records() {
            writer
                .write(&r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), GTF_FILE_ONE_ATTRIB)
    }

    #[test]
    fn test_convert_gtf2_to_gff3() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GffType::GTF2);
        let mut writer = Writer::new(vec![], GffType::GFF3);
        for r in reader.records() {
            writer
                .write(&r.expect("Error reading record"))
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), GFF_FILE_ONE_ATTRIB)
    }
}
