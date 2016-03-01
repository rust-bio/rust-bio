// Copyright 2014 Johannes Köster, Christopher Schröder.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


//! FASTA format reading and writing.
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::fasta;
//! let reader = fasta::Reader::new(io::stdin());
//! ```


use std::io;
use std::io::prelude::*;
use std::ascii::AsciiExt;
use std::collections;
use std::fs;
use std::path::Path;
use std::convert::AsRef;

use csv;


/// A FASTA reader.
pub struct Reader<R: io::Read> {
    reader: io::BufReader<R>,
    line: String,
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
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }

    /// Read next FASTA record into the given `Record`.
    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        if self.line.is_empty() {
            try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() {
                return Ok(());
            }
        }

        if !self.line.starts_with(">") {
            return Err(io::Error::new(io::ErrorKind::Other, "Expected > at record start."));
        }
        record.header.push_str(&self.line);
        loop {
            self.line.clear();
            try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() || self.line.starts_with(">") {
                break;
            }
            record.seq.push_str(&self.line.trim_right());
        }

        Ok(())
    }

    /// Return an iterator over the records of this FastQ file.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}


/// A FASTA index as created by SAMtools (.fai).
pub struct Index {
    inner: collections::HashMap<Vec<u8>, IndexRecord>,
    seqs: Vec<Vec<u8>>,
}


impl Index {
    /// Open a FASTA index from a given `io::Read` instance.
    pub fn new<R: io::Read>(fai: R) -> csv::Result<Self> {
        let mut inner = collections::HashMap::new();
        let mut seqs = vec![];
        let mut fai_reader = csv::Reader::from_reader(fai).delimiter(b'\t').has_headers(false);
        for row in fai_reader.decode() {
            let (name, record): (String, IndexRecord) = try!(row);
            seqs.push(name.clone().into_bytes());
            inner.insert(name.into_bytes(), record);
        }
        Ok(Index {
            inner: inner,
            seqs: seqs,
        })
    }

    /// Open a FASTA index from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: &P) -> csv::Result<Self> {
        match fs::File::open(path) {
            Ok(fai) => Self::new(fai),
            Err(e) => Err(csv::Error::Io(e)),
        }
    }

    /// Open a FASTA index given the corresponding FASTA file path (e.g. for ref.fasta we expect ref.fasta.fai).
    pub fn with_fasta_file<P: AsRef<Path>>(fasta_path: &P) -> csv::Result<Self> {
        let mut ext = fasta_path.as_ref().extension().unwrap().to_str().unwrap().to_owned();
        ext.push_str(".fai");
        let fai_path = fasta_path.as_ref().with_extension(ext);

        Self::from_file(&fai_path)
    }

    /// Return a vector of sequences described in the index.
    pub fn sequences(&self) -> Vec<Sequence> {
        self.seqs
            .iter()
            .map(|name| {
                Sequence {
                    name: name.clone(),
                    len: self.inner.get(name).unwrap().len,
                }
            })
            .collect()
    }
}


/// A FASTA reader with an index as created by SAMtools (.fai).
pub struct IndexedReader<R: io::Read + io::Seek> {
    reader: io::BufReader<R>,
    pub index: Index,
}


impl IndexedReader<fs::File> {
    /// Read from a given file path. This assumes the index ref.fasta.fai to be present for FASTA ref.fasta.
    pub fn from_file<P: AsRef<Path>>(path: &P) -> csv::Result<Self> {
        let index = try!(Index::with_fasta_file(path));

        match fs::File::open(path) {
            Ok(fasta) => Ok(IndexedReader::with_index(fasta, index)),
            Err(e) => Err(csv::Error::Io(e)),
        }
    }
}


impl<R: io::Read + io::Seek> IndexedReader<R> {
    /// Read from a FASTA and its index, both given as `io::Read`. FASTA has to be `io::Seek` in addition.
    pub fn new<I: io::Read>(fasta: R, fai: I) -> csv::Result<Self> {
        let index = try!(Index::new(fai));
        Ok(IndexedReader {
            reader: io::BufReader::new(fasta),
            index: index,
        })
    }

    /// Read from a FASTA and its index, the first given as `io::Read`, the second given as index object.
    pub fn with_index(fasta: R, index: Index) -> Self {
        IndexedReader {
            reader: io::BufReader::new(fasta),
            index: index,
        }
    }

    /// For a given seqname, read the whole sequence into the given vector.
    pub fn read_all(&mut self, seqname: &[u8], seq: &mut Vec<u8>) -> io::Result<()> {
        match self.index.inner.get(seqname) {
            Some(&idx) => self.read(seqname, 0, idx.len, seq),
            None => Err(io::Error::new(io::ErrorKind::Other, "Unknown sequence name.")),
        }
    }

    /// Read the given interval of the given seqname into the given vector (stop position is exclusive).
    pub fn read(&mut self,
                seqname: &[u8],
                start: u64,
                stop: u64,
                seq: &mut Vec<u8>)
                -> io::Result<()> {
        match self.index.inner.get(seqname) {
            Some(idx) => {
                seq.clear();
                // derived from
                // http://www.allenyu.info/item/24-quickly-fetch-sequence-from-samtools-faidx-indexed-fasta-sequences.html
                let line = start / idx.line_bases * idx.line_bytes;
                let line_offset = start % idx.line_bases;
                let offset = idx.offset + line + line_offset;
                let lines = stop / idx.line_bases * idx.line_bytes - line;
                let line_stop = stop % idx.line_bases -
                                if lines == 0 {
                    line_offset
                } else {
                    0
                };

                try!(self.reader.seek(io::SeekFrom::Start(offset)));
                let mut buf = vec![0u8; idx.line_bases as usize];
                for _ in 0..lines {
                    // read full lines
                    try!(self.reader.read(&mut buf));
                    seq.extend(&buf);
                }
                // read last line
                println!("linestop {}", line_stop);
                try!(self.reader.read(&mut buf[..line_stop as usize]));
                seq.extend(&buf[..line_stop as usize]);
                Ok(())
            }
            None => Err(io::Error::new(io::ErrorKind::Other, "Unknown sequence name.")),
        }
    }
}


/// Record of a FASTA index.
#[derive(RustcDecodable, Debug, Copy, Clone)]
struct IndexRecord {
    len: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}


/// A sequence record returned by the FASTA index.
pub struct Sequence {
    pub name: Vec<u8>,
    pub len: u64,
}


/// A Fasta writer.
pub struct Writer<W: io::Write> {
    writer: io::BufWriter<W>,
}


impl Writer<fs::File> {
    /// Write to the given file path.
    pub fn to_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::create(path).map(Writer::new)
    }
}


impl<W: io::Write> Writer<W> {
    /// Create a new Fasta writer.
    pub fn new(writer: W) -> Self {
        Writer { writer: io::BufWriter::new(writer) }
    }

    /// Directly write a Fasta record.
    pub fn write_record(&mut self, record: Record) -> io::Result<()> {
        self.write(record.id().unwrap_or(""), record.desc(), record.seq())
    }

    /// Write a Fasta record with given id, optional description and sequence.
    pub fn write(&mut self, id: &str, desc: Option<&str>, seq: &[u8]) -> io::Result<()> {
        try!(self.writer.write(b">"));
        try!(self.writer.write(id.as_bytes()));
        if desc.is_some() {
            try!(self.writer.write(b" "));
            try!(self.writer.write(desc.unwrap().as_bytes()));
        }
        try!(self.writer.write(b"\n"));
        try!(self.writer.write(seq));
        try!(self.writer.write(b"\n"));

        Ok(())
    }

    /// Flush the writer, ensuring that everything is written.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}


/// A FASTA record.
pub struct Record {
    header: String,
    seq: String,
}


impl Record {
    /// Create a new instance.
    pub fn new() -> Self {
        Record {
            header: String::new(),
            seq: String::new(),
        }
    }

    /// Check if record is empty.
    pub fn is_empty(&self) -> bool {
        self.header.is_empty() && self.seq.is_empty()
    }

    /// Check validity of Fasta record.
    pub fn check(&self) -> Result<(), &str> {
        if self.id().is_none() {
            return Err("Expecting id for FastQ record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }

        Ok(())
    }

    /// Return the id of the record.
    pub fn id(&self) -> Option<&str> {
        self.header[1..].trim_right().splitn(2, " ").next()
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Option<&str> {
        self.header[1..].trim_right().splitn(2, " ").skip(1).next()
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.header.clear();
        self.seq.clear();
    }
}


/// An iterator over the records of a Fasta file.
pub struct Records<R: io::Read> {
    reader: Reader<R>,
}


impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Ok(()) if record.is_empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    const FASTA_FILE: &'static [u8] = b">id desc
ACCGTAGGCTGA
>id2
ATTGTTGTTTTA
";
    const FAI_FILE: &'static [u8] = b"id\t12\t9\t60\t61
";

    #[test]
    fn test_reader() {
        let reader = Reader::new(FASTA_FILE);
        let ids = [Some("id"), Some("id2")];
        let descs = [Some("desc"), None];
        let seqs: [&[u8]; 2] = [b"ACCGTAGGCTGA", b"ATTGTTGTTTTA"];

        for (i, r) in reader.records().enumerate() {
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
        }


        // let record = records.ok().nth(1).unwrap();
    }

    #[test]
    fn test_indexed_reader() {
        let mut reader = IndexedReader::new(io::Cursor::new(FASTA_FILE), FAI_FILE)
                             .ok()
                             .expect("Error reading index");
        let mut seq = Vec::new();
        reader.read(b"id", 1, 5, &mut seq).ok().expect("Error reading sequence.");
        assert_eq!(seq, b"CCGT");
    }

    #[test]
    fn test_writer() {
        let mut writer = Writer::new(Vec::new());
        writer.write("id", Some("desc"), b"ACCGTAGGCTGA").ok().expect("Expected successful write");
        writer.write("id2", None, b"ATTGTTGTTTTA").ok().expect("Expected successful write");
        writer.flush().ok().expect("Expected successful write");
        assert_eq!(writer.writer.get_ref(), &FASTA_FILE);
    }
}
