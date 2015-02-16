// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A FastQ reader.

use std::io;
use std::ascii::AsciiExt;
use std::io::BufRead;


pub struct FastqReader<R: io::Read> {
    reader: io::BufReader<R>,
    sep_line: String
}


impl<R: io::Read> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader { reader: io::BufReader::new(reader), sep_line: String::new() }
    }

    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        try!(self.reader.read_line(&mut record.name));

        if !record.name.is_empty() {
            if !record.name.starts_with("@") {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Expected @ at record start.",
                    None,
                ));
            }
            try!(self.reader.read_line(&mut record.seq));
            try!(self.reader.read_line(&mut self.sep_line));
            try!(self.reader.read_line(&mut record.qual));
        }

        Ok(())
    }

    pub fn records(self, expected_seq_len: usize) -> Records<R> {
        Records { reader: self, expected_seq_len: expected_seq_len }
    }
}


pub struct Record {
    name: String,
    seq: String,
    qual: String,
}


impl Record {
    pub fn new(expected_seq_len: usize) -> Self {
        Record {
            name: String::new(),
            seq: String::with_capacity(expected_seq_len),
            qual: String::with_capacity(expected_seq_len),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.name.is_empty() && self.seq.is_empty() && self.qual.is_empty()
    }

    pub fn check(&self) -> Result<(), &str> {
        let eof = self.qual.ends_with("\n");
        if self.seq.len() != self.qual.len() + eof as usize {
            return Err("Unequal length of sequence an qualities.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }
        if !self.qual.is_ascii() {
            return Err("Non-ascii character found in qualities.");
        }

        Ok(())
    }

    pub fn name(&self) -> &str {
        &self.name[1..self.name.len()-1]
    }

    pub fn seq(&self) -> &[u8] {
        self.seq[..self.seq.len()-1].as_bytes()
    }

    pub fn qual(&self) -> &[u8] {
        self.qual.trim_matches('\n').as_bytes()
    }

    pub fn clear(&mut self) {
        self.name.clear();
        self.seq.clear();
        self.qual.clear();
    }
}


pub struct Records<R: io::Read> {
    reader: FastqReader<R>,
    expected_seq_len: usize,
}


impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        let mut record = Record::new(self.expected_seq_len);
        match self.reader.read(&mut record) {
            Ok(()) if record.is_empty() => None,
            Ok(())   => Some(Ok(record)),
            Err(err) => Some(Err(err))
        }
    }
}
