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
//! let reader = bed::Reader::new(io::stdin());
//! ```


use std::io;
use std::str;

use csv;


pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}


impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader { inner: csv::Reader::from_reader(reader).delimiter(b'\t').has_headers(false) }
    }

    pub fn records(&mut self) -> csv::DecodedRecords<R, Record> {
        self.inner.decode()
    }
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
}


impl<W: io::Write> Writer<W> {
    pub fn new(writer: W) -> Self {
        Writer { inner: csv::Writer::from_writer(writer).delimiter(b'\t').flexible(true) }
    }

    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        self.inner.encode(record)
    }
}


/// A BED record as defined by UCSC (https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
#[derive(RustcEncodable)]
#[derive(Debug)]
pub struct Record {
    chrom: String,
    start: u64,
    end: u64,
    fields: Vec<String>,
}


impl Record {
    pub fn new() -> Self {
        Record { chrom: String::new(), start: 0, end: 0, fields: vec![] }
    }

    pub fn chrom(&self) -> &[u8] {
        &self.chrom.as_bytes()
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn name(&self) -> Option<&str> {
        self.field(3)
    }

    pub fn score(&self) -> Option<&str> {
        self.field(4)
    }

    pub fn is_forward_strand(&self) -> Option<bool> {
        match self.field(5) {
            Some(v) => if v == "+" { Some(true) } else if v == "-" { Some(false) } else { None },
            None    => None
        }
    }

    pub fn field(&self, i: usize) -> Option<&str> {
        let j = i - 3;
        if self.fields.len() > j {
            Some(&self.fields[j])
        }
        else {
            None
        }
    }

    pub fn set_chrom(&mut self, chrom: &[u8]) {
        self.chrom = str::from_utf8(chrom).ok().expect("Error during string conversion.").to_string();
    }

    pub fn set_start(&mut self, start: u64) {
        self.start = start;
    }

    pub fn set_end(&mut self, end: u64) {
        self.end = end;
    }

    pub fn set_name(&mut self, name: &str) {
        if self.fields.len() < 1 {
            self.fields.push(name.to_string());
        }
        else {
            self.fields[0] = name.to_string();
        }
    }

    pub fn set_score(&mut self, score: &str) {
        if self.fields.len() < 1 {
            self.fields.push("".to_string());
        }
        if self.fields.len() < 2 {
            self.fields.push(score.to_string());
        }
        else {
            self.fields[1] = score.to_string();
        }
    }

    pub fn push_field(&mut self, field: &str) {
        self.fields.push(field.to_string());
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    const BED_FILE: &'static [u8] = b"chr1\t5\t500";
//2\t3\t5005\tname2\tup";

    #[test]
    fn test_reader() {
        let chroms = [b"1", b"2"];
        let starts = [5, 2];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];

        let mut reader = Reader::new(BED_FILE);
        for (i, r) in reader.records().enumerate() {
            println!("{:?}", r);
            let record = r.ok().expect("Error reading record");
            assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().unwrap(), names[i]);
            assert_eq!(record.score().unwrap(), scores[i]);
        }
    }
}
