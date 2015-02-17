// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Fasta reading and writing.


use std::io;
use std::io::prelude::*;
use std::ascii::AsciiExt;


pub struct FastaReader<R: io::Read> {
    reader: io::BufReader<R>,
    line: String
}


impl<R: io::Read> FastaReader<R> {
    /// Create a new FastQ reader.
    pub fn new(reader: R) -> Self {
        FastaReader { reader: io::BufReader::new(reader), line: String::new() }
    }

    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        if self.line.is_empty() {
            try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() {
                return Ok(());
            }
        }

        if !self.line.starts_with(">") {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
                None,
            ));
        }
        record.header.push_str(&self.line);
        loop {
            self.line.clear();
            try!(self.reader.read_line(&mut self.line));
            record.seq.push_str(&self.line.trim_right());
            if self.line.is_empty() || self.line.starts_with(">") {
                break;
            }
        }

        Ok(())
    }

    /// Return an iterator over the records of this FastQ file.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}


impl<R: io::Read + io::Seek> FastaReader<R> {
    /// Seek to a given offset. Intended for internal use by IndexedFastaReader.
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        self.reader.get_mut().seek(pos)
    }
}


pub struct Record {
    header: String,
    seq: String,
}


impl Record {
    pub fn new() -> Self {
        Record { header: String::new(), seq: String::new() }
    }

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
        self.header[1..].words().next()
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Vec<&str> {
        self.header[1..].words().skip(1).collect()
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }

    fn clear(&mut self) {
        self.header.clear();
        self.seq.clear();
    }
}


/// An iterator over the records of a Fasta file.
pub struct Records<R: io::Read> {
    reader: FastaReader<R>,
}


impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Ok(()) if record.is_empty() => None,
            Ok(())   => Some(Ok(record)),
            Err(err) => Some(Err(err))
        }
    }
}


/// A Fasta writer.
pub struct FastaWriter<W: io::Write> {
    writer: io::BufWriter<W>,
}


impl<W: io::Write> FastaWriter<W> {
    /// Create a new Fasta writer.
    pub fn new(writer: W) -> Self {
        FastaWriter { writer: io::BufWriter::new(writer) }
    }

    /// Directly write a Fasta record.
    pub fn write_record(&mut self, record: Record) -> io::Result<()> {
        self.write(record.id().unwrap_or(""), &record.desc(), record.seq())
    }

    /// Write a Fasta record with given values.
    ///
    /// # Arguments
    ///
    /// * `id` - the record id
    /// * `desc` - the optional descriptions
    /// * `seq` - the sequence
    pub fn write(&mut self, id: &str, desc: &[&str], seq: &[u8]) -> io::Result<()> {
        try!(self.writer.write(b">"));
        try!(self.writer.write(id.as_bytes()));
        if !desc.is_empty() {
            for d in desc {
                try!(self.writer.write(b" "));
                try!(self.writer.write(d.as_bytes()));
            }
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


#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    const FASTA_FILE: &'static [u8] = b">id desc
ACCGTAGGCTGA
";

    #[test]
    fn test_reader() {
        let reader = FastaReader::new(FASTA_FILE);
        let records: Vec<io::Result<Record>> = reader.records().collect();
        assert!(records.len() == 1);
        for res in records {
            let record = res.ok().unwrap();
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), Some("id"));
            assert_eq!(record.desc(), ["desc"]);
            assert_eq!(record.seq(), b"ACCGTAGGCTGA");
        }
    }

    #[test]
    fn test_writer() {
        let mut writer = FastaWriter::new(Vec::new());
        writer.write("id", &["desc"], b"ACCGTAGGCTGA").ok().expect("Expected successful write");
        writer.flush().ok().expect("Expected successful write");
        assert_eq!(writer.writer.get_ref(), &FASTA_FILE);
    }
}
