//! Readers and writers for common bioinformatics file formats.

pub mod bed;
pub mod fasta;
pub mod fastq;
pub mod gff;
#[cfg(feature = "phylogeny")]
pub mod newick;

/// A generic record.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
    qual: String,
}

impl Record {
    /// Create a new, empty generic record.
    pub fn new() -> Self {
        Record {
            id: String::new(),
            desc: None,
            seq: String::new(),
            qual: String::new(),
        }
    }

    /// Return a fasta::Record object of this generic object
    pub fn fasta_record(&self) -> fasta::Record {
        let desc_str:&str = self.desc.as_ref().unwrap();
        let desc_option:std::option::Option<_> = Some(desc_str);
        let seq:&[u8] = self.seq.as_ref();

        fasta::Record::with_attrs(
            &self.id, desc_option, seq,
        )
    }

    /// Return a fastq::Record object of this generic object
    pub fn fastq_record(&self) -> fastq::Record {
        let desc_str:&str = self.desc.as_ref().unwrap();
        let desc_option:std::option::Option<_> = Some(desc_str);
        let seq:&[u8]  = self.seq.as_ref();
        let qual:&[u8] = self.qual.as_ref();

        fastq::Record::with_attrs(
            &self.id, desc_option, seq, qual
        )
    }
}

/// Create a generic record from a fastq record
impl From<fastq::Record> for Record {
    fn from(fastq: fastq::Record) -> Self {
        Record {
            id:   fastq.id().to_string(),
            desc: Some(fastq.desc().unwrap().to_string()),
            seq:  std::str::from_utf8(fastq.seq()).unwrap().to_string(),
            qual: std::str::from_utf8(fastq.qual()).unwrap().to_string()
        }
    }
}
    

/// Create a generic record from a fasta record
impl From<fasta::Record> for Record {
    fn from(fasta: fasta::Record) -> Self {
        Record {
            id:   fasta.id().to_string(),
            desc: Some(fasta.desc().unwrap().to_string()),
            seq:  std::str::from_utf8(fasta.seq()).unwrap().to_string(),
            qual: String::new(),
        }
    }
}
    

