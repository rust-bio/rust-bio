use std::convert::TryFrom;
use std::io;
use std::io::prelude::*;
use std::mem;

use crate::io::fasta;
use crate::io::fasta::FastaRead;
use crate::io::fastq;
use crate::io::fastq::FastqRead;
use crate::utils::TextSlice;

const BASE_QUAL: u8 = 30;

pub enum Reader<R: io::Read> {
    FASTA(fasta::Reader<R>),
    FASTQ(fastq::Reader<R>),
    Uninitialized(Option<R>),
}

pub enum Record {
    FASTA(fasta::Record),
    FASTQ(fastq::Record),
}

impl From<fasta::Record> for Record {
    fn from(f: fasta::Record) -> Self {
        Record::FASTA(f)
    }
}

impl From<Record> for fasta::Record {
    fn from(value: Record) -> Self {
        match value {
            Record::FASTA(f) => f,
            Record::FASTQ(f) => fasta::Record::with_attrs(f.id(), f.desc(), f.seq()),
        }
    }
}

impl TryFrom<Record> for fastq::Record {
    type Error = &'static str;

    fn try_from(value: Record) -> Result<Self, Self::Error> {
        match value {
            Record::FASTA(_) => Err("Only fastx::Record with qual set can be converted"),
            Record::FASTQ(f) => Ok(f),
        }
    }
}

impl From<fastq::Record> for Record {
    fn from(f: fastq::Record) -> Self {
        Record::FASTQ(f)
    }
}

pub struct Records<R: io::Read> {
    reader: Reader<R>,
    error_has_occured: bool,
}

impl Default for Record {
    fn default() -> Self {
        Record::FASTQ(fastq::Record::new())
    }
}

impl Record {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn with_attrs(
        id: &str,
        desc: Option<&str>,
        seq: TextSlice<'_>,
        qual: Option<&[u8]>,
    ) -> Self {
        match qual {
            Some(q) => Record::FASTQ(fastq::Record::with_attrs(id, desc, seq, q)),
            None => Record::FASTA(fasta::Record::with_attrs(id, desc, seq)),
        }
    }

    pub fn id(&self) -> &str {
        match self {
            Record::FASTA(x) => x.id(),
            Record::FASTQ(x) => x.id(),
        }
    }

    pub fn check(&self) -> Result<(), &str> {
        match self {
            Record::FASTA(x) => x.check(),
            Record::FASTQ(x) => x.check(),
        }
    }

    pub fn desc(&self) -> Option<&str> {
        match self {
            Record::FASTA(x) => x.desc(),
            Record::FASTQ(x) => x.desc(),
        }
    }

    pub fn seq(&self) -> TextSlice<'_> {
        match self {
            Record::FASTA(x) => x.seq(),
            Record::FASTQ(x) => x.seq(),
        }
    }

    pub fn qual(&self) -> Option<&[u8]> {
        match self {
            Record::FASTA(_) => None,
            Record::FASTQ(x) => Some(x.qual()),
        }
    }

    pub fn set_base_qual(&mut self, qual: &[u8]) {
        let new_record = fastq::Record::with_attrs(self.id(), self.desc(), self.seq(), qual);
        mem::replace(self, Record::FASTQ(new_record));
    }
}

impl<R> fasta::FastaRead for Reader<R>
where
    R: io::Read,
{
    fn read(&mut self, record: &mut fasta::Record) -> io::Result<()> {
        if let Reader::Uninitialized(_) = self {
            self.init_reader()?;
        };

        match self {
            Reader::FASTA(x) => x.read(record),
            Reader::FASTQ(x) => {
                let mut fq_record = fastq::Record::new();
                x.read(&mut fq_record)?;
                mem::replace(
                    record,
                    fasta::Record::with_attrs(fq_record.id(), fq_record.desc(), fq_record.seq()),
                );
                Ok(())
            }
            Reader::Uninitialized(_) => unreachable!(),
        }
    }
}

impl<R> fastq::FastqRead for Reader<R>
where
    R: io::Read,
{
    fn read(&mut self, record: &mut fastq::Record) -> io::Result<()> {
        if let Reader::Uninitialized(_) = self {
            self.init_reader()?;
        };

        match self {
            Reader::FASTA(x) => {
                let mut fa_record = fasta::Record::new();
                x.read(&mut fa_record)?;
                mem::replace(
                    record,
                    fastq::Record::with_attrs(
                        fa_record.id(),
                        fa_record.desc(),
                        fa_record.seq(),
                        &vec![BASE_QUAL; fa_record.seq().len()],
                    ),
                );
                Ok(())
            }
            Reader::FASTQ(x) => x.read(record),
            Reader::Uninitialized(_) => unreachable!(),
        }
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader::Uninitialized(Some(reader))
    }

    fn init_reader(&mut self) -> io::Result<()> {
        *self = match mem::replace(self, Reader::Uninitialized(None)) {
            Reader::Uninitialized(r) => {
                if r.is_none() {
                    return Err(io::Error::new(
                        io::ErrorKind::Other,
                        "Empty reader, can't proceed",
                    ));
                };

                let mut reader: io::BufReader<R> = io::BufReader::new(r.unwrap());
                let mut line = String::new();

                reader.read_line(&mut line)?;

                if line.starts_with('>') {
                    Reader::FASTA(fasta::Reader::new_with_line(reader, line))
                } else if line.starts_with('@') {
                    Reader::FASTQ(fastq::Reader::new_with_line(reader, line))
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::Other,
                        "Not a FASTA or FASTQ file",
                    ));
                }
            }
            v => v,
        };

        Ok(())
    }

    pub fn records(self) -> Records<R> {
        Records {
            reader: self,
            error_has_occured: false,
        }
    }
}

impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        if self.error_has_occured {
            return None;
        };

        if let Reader::Uninitialized(_) = self.reader {
            if let Err(x) = self.reader.init_reader() {
                self.error_has_occured = true;
                return Some(Err(x));
            }
        };

        match &mut self.reader {
            Reader::FASTA(reader) => {
                let mut record = fasta::Record::new();
                match reader.read(&mut record) {
                    Ok(()) if record.is_empty() => None,
                    Ok(()) => Some(Ok(Record::FASTA(record))),
                    Err(err) => {
                        self.error_has_occured = true;
                        Some(Err(err))
                    }
                }
            }
            Reader::FASTQ(reader) => {
                let mut record = fastq::Record::new();
                match reader.read(&mut record) {
                    Ok(()) if record.is_empty() => None,
                    Ok(()) => Some(Ok(Record::FASTQ(record))),
                    Err(err) => {
                        self.error_has_occured = true;
                        Some(Err(err))
                    }
                }
            }
            Reader::Uninitialized(_) => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::TryInto;
    use std::io;
    use std::iter;

    use crate::io::fasta;
    use crate::io::fasta::tests::ReaderMock;
    use crate::io::fastq;

    const FASTA_FILE: &'static [u8] = b">id desc
ACCGTAGGCTGA
CCGTAGGCTGAA
CGTAGGCTGAAA
GTAGGCTGAAAA
CCCC
>id2
ATTGTTGTTTTA
ATTGTTGTTTTA
ATTGTTGTTTTA
GGGG
";

    const FASTQ_FILE: &'static [u8] = b"@id desc
ACCGTAGGCTGA
+
IIIIIIJJJJJJ
";

    #[test]
    fn test_fareader() {
        let reader = Reader::new(FASTA_FILE);
        let ids = ["id", "id2"];
        let descs = [Some("desc"), None];
        let seqs: [&[u8]; 2] = [
            b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC",
            b"ATTGTTGTTTTAATTGTTGTTTTAATTGTTGTTTTAGGGG",
        ];

        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
            assert_eq!(record.qual(), None);

            let record: fasta::Record = record.into();
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
        }
    }

    #[test]
    fn test_faread_trait() {
        let path = "genome.fa.gz";
        let mut fa_reader: Box<dyn fasta::FastaRead> = match path.ends_with(".gz") {
            true => Box::new(Reader::new(io::BufReader::new(FASTA_FILE))),
            false => Box::new(Reader::new(FASTA_FILE)),
        };
        let mut record = Record::new().into();
        fa_reader.read(&mut record).unwrap();
        // Check if the returned result is correct.
        assert_eq!(record.check(), Ok(()));
        assert_eq!(record.id(), "id");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(
            record.seq().to_vec(),
            b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC".to_vec()
        );
    }

    #[test]
    fn test_reader_wrong_header() {
        let mut reader = Reader::new(&b"!test\nACGTA\n"[..]);
        let mut record = Record::new().into();
        assert!(
            FastaRead::read(&mut reader, &mut record).is_err(),
            "read() should return Err if FASTA header is malformed"
        );
    }

    #[test]
    fn test_reader_no_id() {
        let mut reader = Reader::new(&b">\nACGTA\n"[..]);
        let mut record = Record::new().into();
        FastaRead::read(&mut reader, &mut record).unwrap();
        assert!(
            record.check().is_err(),
            "check() should return Err if FASTA header is empty"
        );
    }

    #[test]
    fn test_reader_non_ascii_sequence() {
        let mut reader = Reader::new(&b">id\nACGTA\xE2\x98\xB9AT\n"[..]);
        let mut record = Record::new().into();
        FastaRead::read(&mut reader, &mut record).unwrap();
        assert!(
            record.check().is_err(),
            "check() should return Err if FASTA sequence is not ASCII"
        );
    }

    #[test]
    fn test_reader_read_fails() {
        let mut reader = Reader::new(ReaderMock {
            seek_fails: false,
            read_fails: true,
        });
        let mut record = Record::new().into();
        assert!(
            FastaRead::read(&mut reader, &mut record).is_err(),
            "read() should return Err if Read::read fails"
        );
    }

    #[test]
    fn test_fareader_with_fqfile() {
        let mut reader = Reader::new(FASTQ_FILE);
        let mut record = Record::new().into();
        FastaRead::read(&mut reader, &mut record).unwrap();
        assert_eq!(record.id(), "id");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ACCGTAGGCTGA");
    }

    #[test]
    fn test_fqreader_with_fafile() {
        let mut reader = Reader::new(FASTA_FILE);
        let mut record = Record::new().try_into().unwrap();
        FastqRead::read(&mut reader, &mut record).unwrap();
        assert_eq!(record.id(), "id");
        assert_eq!(record.desc(), Some("desc"));
        let seq: Vec<u8> = b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC"
            .into_iter()
            .cloned()
            .collect();
        assert_eq!(record.seq().to_vec(), seq);
        assert_eq!(record.qual().to_vec(), vec![BASE_QUAL; 52]);
    }

    #[test]
    fn test_record_with_attrs() {
        let record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", None);
        assert_eq!(record.id(), "id_str");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ATGCGGG");
    }

    #[test]
    fn test_record_with_attrs_qual() {
        let qual = b"QQQQQQQ";
        let record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", Some(qual));
        assert_eq!(record.id(), "id_str");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ATGCGGG");
        assert_eq!(record.qual(), Some(&qual[..]));
    }

    #[test]
    fn test_reader_read_fails_iter() {
        let reader = Reader::new(ReaderMock {
            seek_fails: false,
            read_fails: true,
        });
        let mut records = reader.records();

        assert!(
            records.next().unwrap().is_err(),
            "next() should return Err if Read::read fails"
        );
        assert!(
            records.next().is_none(),
            "next() should return None after error has occurred"
        );
    }

    #[test]
    fn test_fqreader() {
        let reader = Reader::new(FASTQ_FILE);
        let qual = b"IIIIIIJJJJJJ";
        let records: Vec<io::Result<Record>> = reader.records().collect();
        assert_eq!(records.len(), 1);
        for res in records {
            let record = res.ok().unwrap();
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), "id");
            assert_eq!(record.desc(), Some("desc"));
            assert_eq!(record.seq(), b"ACCGTAGGCTGA");
            assert_eq!(record.qual(), Some(&qual[..]));
        }
    }

    #[test]
    fn test_fqread_trait() {
        let path = "reads.fq.gz";
        let mut fq_reader: Box<dyn fastq::FastqRead> = match path.ends_with(".gz") {
            true => Box::new(Reader::new(io::BufReader::new(FASTQ_FILE))),
            false => Box::new(Reader::new(FASTQ_FILE)),
        };
        // The read method can be called, since it is implemented by
        // `Read`. Right now, the records method would not work.
        let mut record = fastq::Record::new();
        fq_reader.read(&mut record).unwrap();
        // Check if the returned result is correct.
        assert_eq!(record.check(), Ok(()));
        assert_eq!(record.id(), "id");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ACCGTAGGCTGA");
        assert_eq!(record.qual(), b"IIIIIIJJJJJJ");
    }

    #[test]
    fn test_set_basequal() {
        let reader = Reader::new(FASTA_FILE);
        let ids = ["id", "id2"];
        let descs = [Some("desc"), None];
        let seqs: [&[u8]; 2] = [
            b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC",
            b"ATTGTTGTTTTAATTGTTGTTTTAATTGTTGTTTTAGGGG",
        ];
        let quals: [&[u8]; 2] = [&[BASE_QUAL; 52], &[BASE_QUAL; 40]];

        for (i, r) in reader.records().enumerate() {
            let new_qual: Vec<u8> = iter::repeat(40u8).take(quals[i].len()).collect();

            let mut record = r.expect("Error reading record");
            record.set_base_qual(&new_qual);
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
            assert_eq!(record.qual(), Some(&new_qual[..]));

            let record: fastq::Record = record.try_into().unwrap();
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
            assert_eq!(record.qual().to_vec(), new_qual);
        }
    }

    #[test]
    fn test_fa_conversions() {
        let fa_record = fasta::Record::with_attrs("id_str", Some("desc"), b"ATGCGGG");
        let mut record: Record = fa_record.into();

        assert_eq!(record.id(), "id_str");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ATGCGGG");

        record.set_base_qual(&[BASE_QUAL; 7]);

        let fa_record_back: fasta::Record = record.into();

        assert_eq!(fa_record_back.id(), "id_str");
        assert_eq!(fa_record_back.desc(), Some("desc"));
        assert_eq!(fa_record_back.seq(), b"ATGCGGG");
    }

    #[test]
    fn test_fq_conversions() {
        let qual = b"QQQQQQQ";
        let fq_record = fastq::Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", qual);
        let record: Record = fq_record.into();

        assert_eq!(record.id(), "id_str");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ATGCGGG");
        assert_eq!(record.qual(), Some(&qual[..]));

        let fa_record_back: fasta::Record = record.into();

        assert_eq!(fa_record_back.id(), "id_str");
        assert_eq!(fa_record_back.desc(), Some("desc"));
        assert_eq!(fa_record_back.seq(), b"ATGCGGG");
    }
}
