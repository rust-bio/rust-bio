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
pub struct Record<'a> {
    pub chrom: &str,
    pub start: u64,
    pub end: u64,
    fields: Vec<&'a str>
}


impl Record {
    pub fn new() -> Self {
        Record { chrom: "", start: 0, end: 0, fields: vec![] }
    }

    pub fn is_empty() -> bool {
        self.chrom.len() == 0
    }

    pub fn clear() {
        self.chrom = "";
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

    pub fn set_name(&mut self, name: &str) {
        if self.fields.len() < 1 {
            self.fields.push(name);
        }
        else {
            self.fields[0] = name;
        }
    }

    pub fn set_score(&mut self, score: &str) {
        if self.fields.len() < 1 {
            self.fields.push(b"");
        }
        if self.fields.len() < 2 {
            self.fields.push(score);
        }
        else {
            self.fields[1] = score;
        }
    }

    pub fn push_field(&mut self, field: &str) {
        self.fields.push(field);
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
