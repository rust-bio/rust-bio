// Copyright 2019 Magnus Manske
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! GAF2 format reading and writing.
//!
//! GAF2.1 definition : https://geneontology.github.io/docs/go-annotation-file-gaf-format-21/
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::gaf;
//! let reader = gaf::Reader::new(io::stdin(), gaf::GafType::GAF2);
//! ```

use itertools::Itertools;
use multimap::MultiMap;
use regex::Regex;
use std::convert::AsRef;
use std::fs;
use std::io;
use std::path::Path;

use csv;

/// `GafType`
///
/// We have one format in the GAF family.
/// The change is in the last field of GAF.
/// For each type we have key value separator and field separator
#[derive(Debug, Clone, Copy)]
pub enum GafType {
    GAF2,
    Any(u8, u8, u8),
}

impl GafType {
    #[inline]
    fn separator(&self) -> (u8, u8, u8) {
        match *self {
            GafType::GAF2 => (b' ', b'|', b':'),
            GafType::Any(x, y, z) => (x, y, z),
        }
    }
}

/// A GAF reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
    gaf_type: GafType,
}

impl Reader<fs::File> {
    /// Read GAF from given file path in given format.
    pub fn from_file<P: AsRef<Path>>(path: P, fileformat: GafType) -> io::Result<Self> {
        fs::File::open(path).map(|f| Reader::new(f, fileformat))
    }
}

impl<R: io::Read> Reader<R> {
    /// Create a new GAF reader given an instance of `io::Read`, in given format.
    pub fn new(reader: R, fileformat: GafType) -> Self {
        Reader {
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_reader(reader),
            gaf_type: fileformat,
        }
    }

    /// Iterate over all records.
    pub fn records(&mut self) -> Records<R> {
        let (delim, term, vdelim) = self.gaf_type.separator();
        let r1 = format!(
            r" *(?P<key>[^{delim}{term}\t]+){delim}(?P<value>[^{delim}{term}\t]+){term}?",
            delim = delim as char,
            term = term as char
        );
        let db_ref_re = Regex::new(&r1).unwrap();

        let r2 = format!(
            r" *(?P<key>[^{delim}{term}\t]+){delim}(?P<value>[^{delim}{term}\t]+){term}?",
            delim = delim as char,
            term = term as char
        );
        let default_re = Regex::new(&r2).unwrap();

        Records {
            inner: self.inner.deserialize(),
            db_ref_re,
            default_re,
            value_delim: vdelim as char,
        }
    }
}

/// A GAF record.
pub struct Records<'a, R: 'a + io::Read> {
    inner: csv::DeserializeRecordsIter<'a, R, Vec<String>>,
    db_ref_re: Regex,
    default_re: Regex,
    value_delim: char,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner.next().map(|res| {
            res.map(|v| {
                let trim_quotes = |s: &str| s.trim_matches('\'').trim_matches('"').to_owned();

                let db = v[0].clone();
                let db_object_id = v[1].clone();
                let db_object_symbol = v[2].clone();
                let raw_qualifier = v[3].clone();
                let go_id = v[4].clone();
                let raw_db_ref = v[5].clone();
                let evidence_code = v[6].clone();
                let raw_with_from = v[7].clone();
                let aspect = v[8].clone();
                let db_object_name = v[9].clone();
                let raw_db_object_synonym = v[10].clone();
                let db_object_type = v[11].clone();
                let raw_taxon = v[12].clone();
                let date = v[13].clone();
                let assigned_by = v[14].clone();
                let raw_annotation_extension = v[15].clone();
                let gene_product_form_id = v[16].clone();

                let mut db_ref = MultiMap::new();
                for caps in self.db_ref_re.captures_iter(&raw_db_ref.to_string()) {
                    for value in caps["value"].split(self.value_delim) {
                        db_ref.insert(trim_quotes(&caps["key"]), trim_quotes(value));
                    }
                }

                let mut qualifier = Vec::<String>::new();
                for caps in self.default_re.captures_iter(&raw_qualifier.to_string()) {
                    qualifier.push(trim_quotes(&caps["value"]));
                }

                let mut with_from = Vec::<String>::new();
                for caps in self.default_re.captures_iter(&raw_with_from.to_string()) {
                    with_from.push(trim_quotes(&caps["value"]));
                }

                let mut db_object_synonym = Vec::<String>::new();
                for caps in self
                    .default_re
                    .captures_iter(&raw_db_object_synonym.to_string())
                {
                    db_object_synonym.push(trim_quotes(&caps["value"]));
                }

                let mut taxon = Vec::<String>::new();
                for caps in self.default_re.captures_iter(&raw_taxon.to_string()) {
                    taxon.push(trim_quotes(&caps["value"]));
                }

                let mut annotation_extension = Vec::<String>::new();
                for caps in self
                    .default_re
                    .captures_iter(&raw_annotation_extension.to_string())
                {
                    annotation_extension.push(trim_quotes(&caps["value"]));
                }

                Record {
                    db,
                    db_object_id,
                    db_object_symbol,
                    qualifier,
                    go_id,
                    db_ref,
                    evidence_code,
                    with_from,
                    aspect,
                    db_object_name,
                    db_object_synonym,
                    db_object_type,
                    taxon,
                    date,
                    assigned_by,
                    annotation_extension,
                    gene_product_form_id,
                }
            })
        })
    }
}

/// A GAF writer.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    delimiter: char,
    terminator: String,
}

impl Writer<fs::File> {
    /// Write to a given file path in given format.
    pub fn to_file<P: AsRef<Path>>(path: P, fileformat: GafType) -> io::Result<Self> {
        fs::File::create(path).map(|f| Writer::new(f, fileformat))
    }
}

impl<W: io::Write> Writer<W> {
    /// Write to a given writer.
    pub fn new(writer: W, fileformat: GafType) -> Self {
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

    /// Write a given GAF record.
    pub fn write(&mut self, record: &Record) -> csv::Result<()> {
        let db_ref = if !record.db_ref.is_empty() {
            record
                .db_ref
                .iter()
                .map(|(a, b)| format!("{}{}{}", a, self.delimiter, b))
                .join(&self.terminator)
        } else {
            "".to_owned()
        };

        let qualifier = if !record.qualifier.is_empty() {
            record.qualifier.join(&self.terminator)
        } else {
            "".to_owned()
        };

        let with_from = if !record.with_from.is_empty() {
            record.with_from.join(&self.terminator)
        } else {
            "".to_owned()
        };

        let db_object_synonym = if !record.db_object_synonym.is_empty() {
            record.db_object_synonym.join(&self.terminator)
        } else {
            "".to_owned()
        };

        let taxon = if !record.taxon.is_empty() {
            record.taxon.join(&self.terminator)
        } else {
            "".to_owned()
        };

        let annotation_extension = if !record.annotation_extension.is_empty() {
            record.annotation_extension.join(&self.terminator)
        } else {
            "".to_owned()
        };

        let v = //Vec::<String>::new();
        vec![
            record.db.clone(),
            record.db_object_id.clone(),
            record.db_object_symbol.clone(),
            qualifier,
            record.go_id.clone(),
            db_ref,
            record.evidence_code.clone(),
            with_from,
            record.aspect.clone(),
            record.db_object_name.clone(),
            db_object_synonym,
            record.db_object_type.clone(),
            taxon,
            record.date.clone(),
            record.assigned_by.clone(),
            annotation_extension,
            record.gene_product_form_id.clone(),
        ] ;
        self.inner.serialize(v)
    }
}

/// A GAF record
#[derive(Debug, Default, Serialize, Deserialize, Clone)]
pub struct Record {
    db: String,
    db_object_id: String,
    db_object_symbol: String,
    qualifier: Vec<String>,
    go_id: String,
    db_ref: MultiMap<String, String>,
    evidence_code: String,
    with_from: Vec<String>,
    aspect: String,
    db_object_name: String,
    db_object_synonym: Vec<String>,
    db_object_type: String,
    taxon: Vec<String>,
    date: String,
    assigned_by: String,
    annotation_extension: Vec<String>,
    gene_product_form_id: String,
}

impl Record {
    /// Create a new GAF record.
    pub fn new() -> Self {
        Record {
            db: "".to_owned(),
            db_object_id: "".to_owned(),
            db_object_symbol: "".to_owned(),
            qualifier: Vec::<String>::new(),
            go_id: "".to_owned(),
            db_ref: MultiMap::<String, String>::new(),
            evidence_code: "".to_owned(),
            with_from: Vec::<String>::new(),
            aspect: "".to_owned(),
            db_object_name: "".to_owned(),
            db_object_synonym: Vec::<String>::new(),
            db_object_type: "".to_owned(),
            taxon: Vec::<String>::new(),
            date: "".to_owned(),
            assigned_by: "".to_owned(),
            annotation_extension: Vec::<String>::new(),
            gene_product_form_id: "".to_owned(),
        }
    }

    /// Get reference on db
    pub fn db(&self) -> &str {
        &self.db
    }

    /// Get reference on db_object_id
    pub fn db_object_id(&self) -> &str {
        &self.db_object_id
    }

    /// Get reference on db_object_symbol
    pub fn db_object_symbol(&self) -> &str {
        &self.db_object_symbol
    }

    /// Get reference on qualifier
    pub fn qualifier(&self) -> &Vec<String> {
        &self.qualifier
    }

    /// Get reference on go_id
    pub fn go_id(&self) -> &str {
        &self.go_id
    }

    /// Get reference on db_ref
    pub fn db_ref(&self) -> &MultiMap<String, String> {
        &self.db_ref
    }

    /// Get reference on evidence_code
    pub fn evidence_code(&self) -> &str {
        &self.evidence_code
    }

    /// Get reference on with_from
    pub fn with_from(&self) -> &Vec<String> {
        &self.with_from
    }

    /// Get reference on aspect
    pub fn aspect(&self) -> &str {
        &self.aspect
    }

    /// Get reference on db_object_name
    pub fn db_object_name(&self) -> &str {
        &self.db_object_name
    }

    /// Get reference on db_object_synonym
    pub fn db_object_synonym(&self) -> &Vec<String> {
        &self.db_object_synonym
    }

    /// Get reference on db_object_type
    pub fn db_object_type(&self) -> &str {
        &self.db_object_type
    }

    /// Get reference on taxon
    pub fn taxon(&self) -> &Vec<String> {
        &self.taxon
    }

    /// Get reference on date
    pub fn date(&self) -> &str {
        &self.date
    }

    /// Get reference on assigned_by
    pub fn assigned_by(&self) -> &str {
        &self.assigned_by
    }

    /// Get reference on annotation_extension
    pub fn annotation_extension(&self) -> &Vec<String> {
        &self.annotation_extension
    }

    /// Get reference on gene_product_form_id
    pub fn gene_product_form_id(&self) -> &str {
        &self.gene_product_form_id
    }

    /// Get mutable reference on db
    pub fn db_mut(&mut self) -> &mut str {
        &mut self.db
    }

    /// Get mutable reference on db_object_id
    pub fn db_object_id_mut(&mut self) -> &mut str {
        &mut self.db_object_id
    }

    /// Get mutable reference on db_object_symbol
    pub fn db_object_symbol_mut(&mut self) -> &mut str {
        &mut self.db_object_symbol
    }

    /// Get mutable reference on qualifier
    pub fn qualifier_mut(&mut self) -> &mut Vec<String> {
        &mut self.qualifier
    }

    /// Get mutable reference on go_id
    pub fn go_id_mut(&mut self) -> &mut str {
        &mut self.go_id
    }

    /// Get mutable reference on db_ref
    pub fn db_ref_mut(&mut self) -> &mut MultiMap<String, String> {
        &mut self.db_ref
    }

    /// Get mutable reference on evidence_code
    pub fn evidence_code_mut(&mut self) -> &mut str {
        &mut self.evidence_code
    }

    /// Get mutable reference on with_from
    pub fn with_from_mut(&mut self) -> &mut Vec<String> {
        &mut self.with_from
    }

    /// Get mutable reference on aspect
    pub fn aspect_mut(&mut self) -> &mut str {
        &mut self.aspect
    }

    /// Get mutable reference on db_object_name
    pub fn db_object_name_mut(&mut self) -> &mut str {
        &mut self.db_object_name
    }

    /// Get mutable reference on db_object_synonym
    pub fn db_object_synonym_mut(&mut self) -> &mut Vec<String> {
        &mut self.db_object_synonym
    }

    /// Get mutable reference on db_object_type
    pub fn db_object_type_mut(&mut self) -> &mut str {
        &mut self.db_object_type
    }

    /// Get mutable reference on taxon
    pub fn taxon_mut(&mut self) -> &mut Vec<String> {
        &mut self.taxon
    }

    /// Get mutable reference on date
    pub fn date_mut(&mut self) -> &mut str {
        &mut self.date
    }

    /// Get mutable reference on assigned_by
    pub fn assigned_by_mut(&mut self) -> &mut str {
        &mut self.assigned_by
    }

    /// Get mutable reference on annotation_extension
    pub fn annotation_extension_mut(&mut self) -> &mut Vec<String> {
        &mut self.annotation_extension
    }

    /// Get mutable reference on gene_product_form_id
    pub fn gene_product_form_id_mut(&mut self) -> &mut str {
        &mut self.gene_product_form_id
    }
}

/*
#[cfg(test)]
mod tests {
    use super::*;
    use bio_types::strand::Strand;
    use multimap::MultiMap;

    const GAF_FILE: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\t\
Note=Removed,Obsolete;ID=test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote=ATP-dependent protease subunit HslV;\
ID=PRO_0000148105";
    const GAF_FILE_WITH_COMMENT: &'static [u8] = b"#comment
P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\t\
Note=Removed,Obsolete;ID=test
#comment
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote=ATP-dependent protease subunit HslV;\
ID=PRO_0000148105";
    //required because MultiMap iter on element randomly
    const GAF_FILE_ONE_ATTRIB: &'static [u8] =
        b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote=Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID=PRO_0000148105
";

    const GTF_FILE: &'static [u8] =
        b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed;ID test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote ATP-dependent;ID PRO_0000148105
";

    // Another variant of GTF file, modified from a published GENCODE GTF file.
    const GTF_FILE_2: &'static [u8] = b"chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\t\
gene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\";\
transcript_id \"ENST00000456328.2\"; gene_type \"transcribed_unprocessed_pseudogene\"";

    // GTF file with duplicate attribute keys, taken from a published GENCODE GTF file.
    const GTF_FILE_DUP_ATTR_KEYS: &'static [u8] = b"chr1\tENSEMBL\ttranscript\t182393\t\
184158\t.\t+\t.\tgene_id \"ENSG00000279928.1\"; transcript_id \"ENST00000624431.1\";\
gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"FO538757.2\";\
transcript_type \"protein_coding\"; transcript_status \"KNOWN\";\
transcript_name \"FO538757.2-201\"; level 3; protein_id \"ENSP00000485457.1\";\
transcript_support_level \"1\"; tag \"basic\"; tag \"appris_principal_1\";";

    //required because MultiMap iter on element randomly
    const GTF_FILE_ONE_ATTRIB: &'static [u8] =
        b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID PRO_0000148105
";

    #[test]
    fn test_reader_gaf3() {
        let seqname = ["P0A7B8", "P0A7B8"];
        let source = ["UniProtKB", "UniProtKB"];
        let feature_type = ["Initiator methionine", "Chain"];
        let starts = [1, 2];
        let ends = [1, 176];
        let scores = [None, Some(50)];
        let strand = [None, Some(Strand::Forward)];
        let frame = [".", "."];
        let mut attributes = [MultiMap::new(), MultiMap::new()];
        attributes[0].insert("ID".to_owned(), "test".to_owned());
        attributes[0].insert("Note".to_owned(), "Removed".to_owned());
        attributes[0].insert("Note".to_owned(), "Obsolete".to_owned());
        attributes[1].insert("ID".to_owned(), "PRO_0000148105".to_owned());
        attributes[1].insert(
            "Note".to_owned(),
            "ATP-dependent protease subunit HslV".to_owned(),
        );

        let mut reader = Reader::new(GAF_FILE, GafType::GAF3);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
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

        let mut reader = Reader::new(GAF_FILE_WITH_COMMENT, GafType::GAF3);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
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
        let mut attributes = [MultiMap::new(), MultiMap::new()];
        attributes[0].insert("ID".to_owned(), "test".to_owned());
        attributes[0].insert("Note".to_owned(), "Removed".to_owned());
        attributes[1].insert("ID".to_owned(), "PRO_0000148105".to_owned());
        attributes[1].insert("Note".to_owned(), "ATP-dependent".to_owned());

        let mut reader = Reader::new(GTF_FILE, GafType::GTF2);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
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
    fn test_reader_gtf2_2() {
        let seqname = ["chr1", "chr1"];
        let source = ["HAVANA", "HAVANA"];
        let feature_type = ["gene", "transcript"];
        let starts = [11869, 11869];
        let ends = [14409, 14409];
        let scores = [None, None];
        let strand = [Some(Strand::Forward), Some(Strand::Forward)];
        let frame = [".", "."];
        let mut attributes = [MultiMap::new(), MultiMap::new()];
        attributes[0].insert("gene_id".to_owned(), "ENSG00000223972.5".to_owned());
        attributes[0].insert(
            "gene_type".to_owned(),
            "transcribed_unprocessed_pseudogene".to_owned(),
        );
        attributes[1].insert("gene_id".to_owned(), "ENSG00000223972.5".to_owned());
        attributes[1].insert("transcript_id".to_owned(), "ENST00000456328.2".to_owned());
        attributes[1].insert(
            "gene_type".to_owned(),
            "transcribed_unprocessed_pseudogene".to_owned(),
        );

        let mut reader = Reader::new(GTF_FILE_2, GafType::GTF2);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
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
    fn test_reader_gtf2_dup_attr_keys() {
        let mut reader = Reader::new(GTF_FILE_DUP_ATTR_KEYS, GafType::GTF2);
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
    fn test_writer_gaf3() {
        let mut reader = Reader::new(GAF_FILE_ONE_ATTRIB, GafType::GAF3);
        let mut writer = Writer::new(vec![], GafType::GAF3);
        for r in reader.records() {
            writer
                .write(&r.ok().expect("Error reading record"))
                .ok()
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), GAF_FILE_ONE_ATTRIB)
    }

    #[test]
    fn test_writer_gtf2() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GafType::GTF2);
        let mut writer = Writer::new(vec![], GafType::GTF2);
        for r in reader.records() {
            writer
                .write(&r.ok().expect("Error reading record"))
                .ok()
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), GTF_FILE_ONE_ATTRIB)
    }

    #[test]
    fn test_convert_gtf2_to_gaf3() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GafType::GTF2);
        let mut writer = Writer::new(vec![], GafType::GAF3);
        for r in reader.records() {
            writer
                .write(&r.ok().expect("Error reading record"))
                .ok()
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.into_inner().unwrap(), GAF_FILE_ONE_ATTRIB)
    }
}
*/
