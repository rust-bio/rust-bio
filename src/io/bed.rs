// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


//! BED reading and writing.
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::bed;
//! let example = b"1\t5\t5000\tname1\t0.5".as_slice();
//! let mut reader = bed::Reader::new(example);
//! let mut writer = bed::Writer::new(vec![]);
//! for record in reader.records() {
//!     let rec = record.ok().expect("Error reading record.");
//!     println!("{}", rec.chrom());
//!     writer.write(rec).ok().expect("Error writing record.");
//! }
//! ```


use std::io;
use std::path;
use std::fs;


use csv;


/// A BED reader.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}


impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader { inner: csv::Reader::from_reader(reader).delimiter(b'\t').has_headers(false) }
    }

    pub fn from_file<P: path::AsPath>(path: P) -> io::Result<Reader<fs::File>> {
        fs::File::open(path).map(|f| Reader::new(f))
    }

    pub fn records(&mut self) -> Records<R> {
        Records { inner: self.inner.decode() }
    }
}


pub struct Records<'a, R: 'a +io::Read> {
    inner: csv::DecodedRecords<'a, R, (String, u64, u64, Vec<String>)>
}


impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner.next().map(|res| match res {
            Ok((chrom, start, end, aux)) => Ok(Record {
                chrom: chrom, start: start, end: end, aux: aux
            }),
            Err(e) => Err(e)
        })
    }
}


/// A BED writer.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
}


impl<W: io::Write> Writer<W> {
    pub fn new(writer: W) -> Self {
        Writer { inner: csv::Writer::from_writer(writer).delimiter(b'\t').flexible(true) }
    }

    pub fn from_file<P: path::AsPath>(path: P) -> io::Result<Writer<fs::File>> {
        fs::File::create(path).map(|f| Writer::new(f))
    }

    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        if record.aux.len() == 0 {
            self.inner.encode((record.chrom, record.start, record.end))
        }
        else {
            self.inner.encode(record)
        }
    }
}


/// A BED record as defined by BEDtools (http://bedtools.readthedocs.org/en/latest/content/general-usage.html)
#[derive(RustcEncodable)]
pub struct Record {
    chrom: String,
    start: u64,
    end: u64,
    aux: Vec<String>
}


impl Record {
    pub fn new() -> Self {
        Record { chrom: "".to_string(), start: 0, end: 0, aux: vec![] }
    }

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

    pub fn name(&self) -> Option<&str> {
        self.aux(3)
    }

    pub fn score(&self) -> Option<&str> {
        self.aux(4)
    }

    pub fn strand(&self) -> Option<Strand> {
        match self.aux(5) {
            Some("+") => Some(Strand::Forward),
            Some("-") => Some(Strand::Reverse),
            _         => None
        }
    }

    pub fn aux(&self, i: usize) -> Option<&str> {
        let j = i - 3;
        if j < self.aux.len() {
            Some(&self.aux[j])
        }
        else {
            None
        }
    }

    pub fn set_chrom(&mut self, chrom: &str) {
        self.chrom = chrom.to_string();
    }

    pub fn set_start(&mut self, start: u64) {
        self.start = start;
    }

    pub fn set_end(&mut self, end: u64) {
        self.end = end;
    }

    pub fn set_name(&mut self, name: &str) {
        if self.aux.len() < 1 {
            self.aux.push(name.to_string());
        }
        else {
            self.aux[0] = name.to_string();
        }
    }

    pub fn set_score(&mut self, score: &str) {
        if self.aux.len() < 1 {
            self.aux.push("".to_string());
        }
        if self.aux.len() < 2 {
            self.aux.push(score.to_string());
        }
        else {
            self.aux[1] = score.to_string();
        }
    }

    pub fn push_aux(&mut self, field: &str) {
        self.aux.push(field.to_string());
    }
}


pub enum Strand {
    Forward,
    Reverse,
}


#[cfg(test)]
mod tests {
    use super::*;

    const BED_FILE: &'static [u8] = b"1\t5\t5000\tname1\tup
2\t3\t5005\tname2\tup
";

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

    #[test]
    fn test_writer() {
        let mut reader = Reader::new(BED_FILE);
        let mut writer = Writer::new(vec![]);
        for r in reader.records() {
            writer.write(r.ok().expect("Error reading record")).ok().expect("Error writing record");
        }
        assert_eq!(writer.inner.as_bytes(), BED_FILE);
    }
}
