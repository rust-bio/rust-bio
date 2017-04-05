// Copyright 2014-2016 Johannes Köster, Christopher Schröder.
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
use std::cmp::min;

use csv;

use utils::{TextSlice, Text};


/// Maximum size of temporary buffer used for reading indexed FASTA files.
const MAX_FASTA_BUFFER_SIZE: usize = 512;


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

        if !self.line.starts_with('>') {
            return Err(io::Error::new(io::ErrorKind::Other, "Expected > at record start."));
        }
        record.header.push_str(&self.line);
        loop {
            self.line.clear();
            try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() || self.line.starts_with('>') {
                break;
            }
            record.seq.push_str(self.line.trim_right());
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
    inner: collections::HashMap<String, IndexRecord>,
    seqs: Vec<String>,
}


impl Index {
    /// Open a FASTA index from a given `io::Read` instance.
    pub fn new<R: io::Read>(fai: R) -> csv::Result<Self> {
        let mut inner = collections::HashMap::new();
        let mut seqs = vec![];
        let mut fai_reader = csv::Reader::from_reader(fai)
            .delimiter(b'\t')
            .has_headers(false);
        for row in fai_reader.decode() {
            let (name, record): (String, IndexRecord) = try!(row);
            seqs.push(name.clone());
            inner.insert(name, record);
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

    /// Open a FASTA index given the corresponding FASTA file path.
    /// That is, for ref.fasta we expect ref.fasta.fai.
    pub fn with_fasta_file<P: AsRef<Path>>(fasta_path: &P) -> csv::Result<Self> {
        let mut fai_path = fasta_path.as_ref().as_os_str().to_owned();
        fai_path.push(".fai");

        Self::from_file(&fai_path)
    }

    /// Return a vector of sequences described in the index.
    pub fn sequences(&self) -> Vec<Sequence> {
        self.seqs
            .iter()
            .map(|name| {
                     Sequence {
                         name: name.clone(),
                         len: self.inner[name].len,
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
    /// Read from a given file path. This assumes the index ref.fasta.fai to be
    /// present for FASTA ref.fasta.
    pub fn from_file<P: AsRef<Path>>(path: &P) -> csv::Result<Self> {
        let index = try!(Index::with_fasta_file(path));

        match fs::File::open(path) {
            Ok(fasta) => Ok(IndexedReader::with_index(fasta, index)),
            Err(e) => Err(csv::Error::Io(e)),
        }
    }
}


impl<R: io::Read + io::Seek> IndexedReader<R> {
    /// Read from a FASTA and its index, both given as `io::Read`. FASTA has to
    /// be `io::Seek` in addition.
    pub fn new<I: io::Read>(fasta: R, fai: I) -> csv::Result<Self> {
        let index = try!(Index::new(fai));
        Ok(IndexedReader {
               reader: io::BufReader::new(fasta),
               index: index,
           })
    }

    /// Read from a FASTA and its index, the first given as `io::Read`, the
    /// second given as index object.
    pub fn with_index(fasta: R, index: Index) -> Self {
        IndexedReader {
            reader: io::BufReader::new(fasta),
            index: index,
        }
    }

    /// For a given seqname, read the whole sequence into the given vector.
    pub fn read_all(&mut self, seqname: &str, seq: &mut Text) -> io::Result<()> {
        let idx = self.idx(seqname)?;

        self.read_into_buffer(idx, 0, idx.len, seq)
    }

    /// Read the given interval of the given seqname into the given vector
    /// (stop position is exclusive).
    pub fn read(&mut self, seqname: &str, start: u64, stop: u64, seq: &mut Text) -> io::Result<()> {
        let idx = self.idx(seqname)?;

        self.read_into_buffer(idx, start, stop, seq)
    }


    /// For a given seqname, return an iterator yielding that sequence.
    pub fn read_all_iter(&mut self, seqname: &str) -> io::Result<IndexedReaderIterator<R>> {
        let idx = self.idx(seqname)?;

        self.read_into_iter(idx, 0, idx.len)
    }

    /// For a given seqname and a given range in that sequence, return an
    /// iterator yielding the corresponding sequence.
    pub fn read_iter(&mut self,
                     seqname: &str,
                     start: u64,
                     stop: u64)
                     -> io::Result<IndexedReaderIterator<R>> {
        let idx = self.idx(seqname)?;

        self.read_into_iter(idx, start, stop)
    }

    fn read_into_buffer(&mut self,
                        idx: IndexRecord,
                        start: u64,
                        stop: u64,
                        seq: &mut Text)
                        -> io::Result<()> {
        if stop > idx.len {
            return Err(io::Error::new(io::ErrorKind::Other,
                                      "FASTA read interval was out of bounds"));
        } else if start > stop {
            return Err(io::Error::new(io::ErrorKind::Other, "Invalid query interval"));
        }

        let mut bases_left = stop - start;
        let mut line_offset = self.seek_to(&idx, start)?;

        seq.clear();
        while bases_left > 0 {
            bases_left -= self.read_line(&idx, &mut line_offset, bases_left, seq)?;
        }

        Ok(())
    }

    fn read_into_iter(&mut self,
                      idx: IndexRecord,
                      start: u64,
                      stop: u64)
                      -> io::Result<IndexedReaderIterator<R>> {
        if stop > idx.len {
            return Err(io::Error::new(io::ErrorKind::Other,
                                      "FASTA read interval was out of bounds"));
        } else if start > stop {
            return Err(io::Error::new(io::ErrorKind::Other, "Invalid query interval"));
        }

        let bases_left = stop - start;
        let line_offset = self.seek_to(&idx, start)?;
        let capacity = min(MAX_FASTA_BUFFER_SIZE,
                           min(bases_left, idx.line_bases) as usize);

        Ok(IndexedReaderIterator {
               reader: self,
               record: idx,
               bases_left: bases_left,
               line_offset: line_offset,
               buf: Vec::with_capacity(capacity),
               buf_idx: 0,
           })
    }

    /// Return the IndexRecord for the given sequence name or io::Result::Err
    fn idx(&self, seqname: &str) -> io::Result<IndexRecord> {
        match self.index.inner.get(seqname) {
            Some(idx) => Ok(idx.clone()),
            None => Err(io::Error::new(io::ErrorKind::Other, "Unknown sequence name.")),
        }
    }

    /// Seek to the given position in the specified FASTA record. The position
    /// of the cursor on the line that the seek ended on is returned.
    fn seek_to(&mut self, idx: &IndexRecord, start: u64) -> io::Result<u64> {
        assert!(start <= idx.len);

        let line_offset = start % idx.line_bases;
        let line_start = start / idx.line_bases * idx.line_bytes;
        let offset = idx.offset + line_start + line_offset;
        try!(self.reader.seek(io::SeekFrom::Start(offset)));

        Ok(line_offset)
    }

    /// Tries to read up to `bases_left` bases from the current line into `buf`,
    /// returning the actual number of bases read. Depending on the amount of
    /// whitespace per line, the current `line_offset`, and the amount of bytes
    /// returned from `BufReader::fill_buf`, this function may return Ok(0)
    /// multiple times in a row.
    fn read_line(&mut self,
                 idx: &IndexRecord,
                 line_offset: &mut u64,
                 bases_left: u64,
                 buf: &mut Vec<u8>)
                 -> io::Result<u64> {
        let (bytes_to_read, bytes_to_keep) = {
            let src = self.reader.fill_buf()?;
            if src.is_empty() {
                return Err(io::Error::new(io::ErrorKind::UnexpectedEof,
                                          "FASTA file is truncated."));
            }

            let bases_on_line = idx.line_bases - min(idx.line_bases, *line_offset);
            let bases_in_buffer = min(src.len() as u64, bases_on_line);

            let (bytes_to_read, bytes_to_keep) = if bases_in_buffer < bases_left {
                let bytes_to_read = min(src.len() as u64, idx.line_bytes - *line_offset);

                (bytes_to_read, bases_in_buffer)
            } else {
                (bases_left, bases_left)
            };

            buf.extend_from_slice(&src[..bytes_to_keep as usize]);
            (bytes_to_read, bytes_to_keep)
        };

        self.reader.consume(bytes_to_read as usize);

        assert!(bytes_to_read > 0);
        *line_offset += bytes_to_read;
        if *line_offset >= idx.line_bytes {
            *line_offset = 0;
        }

        Ok(bytes_to_keep)
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
#[derive(Debug, PartialEq)]
pub struct Sequence {
    pub name: String,
    pub len: u64,
}


pub struct IndexedReaderIterator<'a, R: io::Read + io::Seek + 'a> {
    reader: &'a mut IndexedReader<R>,
    record: IndexRecord,
    bases_left: u64,
    line_offset: u64,
    buf: Vec<u8>,
    buf_idx: usize,
}


impl<'a, R: io::Read + io::Seek + 'a> IndexedReaderIterator<'a, R> {
    fn fill_buffer(&mut self) -> io::Result<()> {
        assert!(self.bases_left > 0);

        self.buf.clear();
        let bases_to_read = min(self.buf.capacity() as u64, self.bases_left);

        // May loop one or more times; see IndexedReader::read_line.
        while self.buf.is_empty() {
            self.bases_left -= self.reader
                .read_line(&self.record,
                           &mut self.line_offset,
                           bases_to_read,
                           &mut self.buf)?;
        }

        self.buf_idx = 0;
        Ok(())
    }
}


impl<'a, R: io::Read + io::Seek + 'a> Iterator for IndexedReaderIterator<'a, R> {
    type Item = io::Result<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buf_idx < self.buf.len() {
            let item = Some(Ok(self.buf[self.buf_idx]));
            self.buf_idx += 1;
            item
        } else if self.bases_left > 0 {
            if let Err(e) = self.fill_buffer() {
                return Some(Err(e));
            }

            self.buf_idx = 1;
            Some(Ok(self.buf[0]))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let hint = self.bases_left as usize + (self.buf.len() - self.buf_idx);

        (hint, Some(hint))
    }
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
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write(record.id().unwrap_or(""), record.desc(), record.seq())
    }

    /// Write a Fasta record with given id, optional description and sequence.
    pub fn write(&mut self, id: &str, desc: Option<&str>, seq: TextSlice) -> io::Result<()> {
        try!(self.writer.write_all(b">"));
        try!(self.writer.write_all(id.as_bytes()));
        if desc.is_some() {
            try!(self.writer.write_all(b" "));
            try!(self.writer.write_all(desc.unwrap().as_bytes()));
        }
        try!(self.writer.write_all(b"\n"));
        try!(self.writer.write_all(seq));
        try!(self.writer.write_all(b"\n"));

        Ok(())
    }

    /// Flush the writer, ensuring that everything is written.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}


/// A FASTA record.
#[derive(Default)]
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
        self.header[1..].trim_right().splitn(2, ' ').nth(0)
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Option<&str> {
        self.header[1..].trim_right().splitn(2, ' ').nth(1)
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> TextSlice {
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
    const FAI_FILE: &'static [u8] = b"id\t52\t9\t12\t13
id2\t40\t71\t12\t13
";

    const TRUNCATED_FASTA: &'static [u8] = b">id desc\nACCGTAGGCTGA";

    const FASTA_FILE_CRLF: &'static [u8] = b">id desc\r
ACCGTAGGCTGA\r
CCGTAGGCTGAA\r
CGTAGGCTGAAA\r
GTAGGCTGAAAA\r
CCCC\r
>id2\r
ATTGTTGTTTTA\r
ATTGTTGTTTTA\r
ATTGTTGTTTTA\r
GGGG\r
";
    const FAI_FILE_CRLF: &'static [u8] = b"id\t52\t10\t12\t14\r
id2\t40\t78\t12\t14\r
";

    const FASTA_FILE_NO_TRAILING_LF: &'static [u8] = b">id desc
GTAGGCTGAAAA
CCCC";
    const FAI_FILE_NO_TRAILING_LF: &'static [u8] = b"id\t16\t9\t12\t13";


    const WRITE_FASTA_FILE: &'static [u8] = b">id desc
ACCGTAGGCTGA
>id2
ATTGTTGTTTTA
";

    struct ReaderMock {
        seek_fails: bool,
        read_fails: bool,
    }

    impl Read for ReaderMock {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.read_fails {
                Err(io::Error::new(io::ErrorKind::Other, "Read set to fail"))
            } else {
                Ok(buf.len())
            }
        }
    }

    impl Seek for ReaderMock {
        fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
            if let io::SeekFrom::Start(pos) = pos {
                if self.seek_fails {
                    Err(io::Error::new(io::ErrorKind::Other, "Seek set to fail"))
                } else {
                    Ok(pos)
                }
            } else {
                unimplemented!();
            }
        }
    }

    #[test]
    fn test_reader() {
        let reader = Reader::new(FASTA_FILE);
        let ids = [Some("id"), Some("id2")];
        let descs = [Some("desc"), None];
        let seqs: [&[u8]; 2] = [b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC",
                                b"ATTGTTGTTTTAATTGTTGTTTTAATTGTTGTTTTAGGGG"];

        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
        }
    }

    #[test]
    fn test_index_sequences() {
        let reader = IndexedReader::new(io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let sequences = reader.index.sequences();
        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences[0],
                   Sequence {
                       name: "id".into(),
                       len: 52,
                   });
        assert_eq!(sequences[1],
                   Sequence {
                       name: "id2".into(),
                       len: 40,
                   });
    }

    #[test]
    fn test_indexed_reader() {
        _test_indexed_reader(&FASTA_FILE, &FAI_FILE, _read_buffer);
        _test_indexed_reader_truncated(_read_buffer);
    }

    #[test]
    fn test_indexed_reader_crlf() {
        _test_indexed_reader(&FASTA_FILE_CRLF, &FAI_FILE_CRLF, _read_buffer);
    }

    #[test]
    fn test_indexed_reader_iter() {
        _test_indexed_reader(&FASTA_FILE, &FAI_FILE, _read_iter);
        _test_indexed_reader_truncated(_read_iter);
    }

    #[test]
    fn test_indexed_reader_iter_crlf() {
        _test_indexed_reader(&FASTA_FILE_CRLF, &FAI_FILE_CRLF, _read_iter);
    }

    fn _test_indexed_reader<'a, F>(fasta: &'a [u8], fai: &'a [u8], read: F)
        where F: Fn(&mut IndexedReader<io::Cursor<&'a [u8]>>, &str, u64, u64) -> io::Result<Vec<u8>>
    {
        let mut reader = IndexedReader::new(io::Cursor::new(fasta), fai).unwrap();

        // Test reading various substrings of the sequence
        assert_eq!(read(&mut reader, "id", 1, 5).unwrap(), b"CCGT");
        assert_eq!(read(&mut reader, "id", 1, 31).unwrap(),
                   b"CCGTAGGCTGACCGTAGGCTGAACGTAGGC");
        assert_eq!(read(&mut reader, "id", 13, 23).unwrap(), b"CGTAGGCTGA");
        assert_eq!(read(&mut reader, "id", 36, 52).unwrap(),
                   b"GTAGGCTGAAAACCCC");
        assert_eq!(read(&mut reader, "id2", 12, 40).unwrap(),
                   b"ATTGTTGTTTTAATTGTTGTTTTAGGGG");
        assert_eq!(read(&mut reader, "id2", 12, 12).unwrap(), b"");
        assert_eq!(read(&mut reader, "id2", 12, 13).unwrap(), b"A");
        // Minimal sequence spanning new-line
        assert_eq!(read(&mut reader, "id", 11, 13).unwrap(), b"AC");

        assert!(read(&mut reader, "id2", 12, 11).is_err());
        assert!(read(&mut reader, "id2", 12, 1000).is_err());
        assert!(read(&mut reader, "id3", 0, 1).is_err());
    }

    fn _test_indexed_reader_truncated<'a, F>(read: F)
        where F: Fn(&mut IndexedReader<io::Cursor<&'a [u8]>>, &str, u64, u64) -> io::Result<Vec<u8>>
    {
        let mut reader = IndexedReader::new(io::Cursor::new(TRUNCATED_FASTA), FAI_FILE).unwrap();

        assert_eq!(read(&mut reader, "id", 0, 12).unwrap(), b"ACCGTAGGCTGA");
        assert!(read(&mut reader, "id", 0, 13).is_err()); // read past EOF
        assert!(read(&mut reader, "id", 36, 52).is_err()); // seek and read past EOF
        assert!(read(&mut reader, "id2", 12, 40).is_err()); // seek and read past EOF
    }

    fn _read_buffer<T>(reader: &mut IndexedReader<T>,
                       seqname: &str,
                       start: u64,
                       stop: u64)
                       -> io::Result<Vec<u8>>
        where T: Seek + Read
    {
        let mut seq = vec![];
        reader.read(seqname, start, stop, &mut seq)?;

        Ok(seq)
    }

    fn _read_iter<T>(reader: &mut IndexedReader<T>,
                     seqname: &str,
                     start: u64,
                     stop: u64)
                     -> io::Result<Vec<u8>>
        where T: Seek + Read
    {
        let mut seq = vec![];
        for nuc in reader.read_iter(seqname, start, stop)? {
            seq.push(nuc?);
        }

        Ok(seq)
    }

    #[test]
    fn test_indexed_reader_all() {
        _test_indexed_reader_all(&FASTA_FILE, &FAI_FILE, _read_buffer_all);
    }

    #[test]
    fn test_indexed_reader_crlf_all() {
        _test_indexed_reader_all(&FASTA_FILE_CRLF, &FAI_FILE_CRLF, _read_buffer_all);
    }

    #[test]
    fn test_indexed_reader_iter_all() {
        _test_indexed_reader_all(&FASTA_FILE, &FAI_FILE, _read_iter_all);
    }

    #[test]
    fn test_indexed_reader_iter_crlf_all() {
        _test_indexed_reader_all(&FASTA_FILE_CRLF, &FAI_FILE_CRLF, _read_iter_all);
    }

    fn _test_indexed_reader_all<'a, F>(fasta: &'a [u8], fai: &'a [u8], read: F)
        where F: Fn(&mut IndexedReader<io::Cursor<&'a [u8]>>, &str) -> io::Result<Vec<u8>>
    {
        let mut reader = IndexedReader::new(io::Cursor::new(fasta), fai).unwrap();

        assert_eq!(read(&mut reader, "id").unwrap(),
                   &b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC"[..]);
        assert_eq!(read(&mut reader, "id2").unwrap(),
                   &b"ATTGTTGTTTTAATTGTTGTTTTAATTGTTGTTTTAGGGG"[..]);
    }

    fn _read_buffer_all<T>(reader: &mut IndexedReader<T>, seqname: &str) -> io::Result<Vec<u8>>
        where T: Seek + Read
    {
        let mut seq = vec![];
        reader.read_all(seqname, &mut seq)?;

        Ok(seq)
    }

    fn _read_iter_all<T>(reader: &mut IndexedReader<T>, seqname: &str) -> io::Result<Vec<u8>>
        where T: Seek + Read
    {
        let mut seq = vec![];
        for nuc in reader.read_all_iter(seqname)? {
            seq.push(nuc?);
        }

        Ok(seq)
    }

    #[test]
    fn test_indexed_reader_iter_size_hint() {
        let mut reader = IndexedReader::new(io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
        let mut iterator = reader.read_iter("id", 2, 4).unwrap();

        assert_eq!(iterator.size_hint(), (2, Some(2)));
        assert_eq!(iterator.next().unwrap().unwrap(), b'C');
        assert_eq!(iterator.size_hint(), (1, Some(1)));
        assert_eq!(iterator.next().unwrap().unwrap(), b'G');
        assert_eq!(iterator.size_hint(), (0, Some(0)));
        assert!(iterator.next().is_none());
        assert_eq!(iterator.size_hint(), (0, Some(0)));
    }

    #[test]
    fn test_indexed_reader_reused_buffer() {
        let mut reader = IndexedReader::new(io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
        let mut seq = Vec::new();

        reader.read("id", 1, 5, &mut seq).unwrap();
        assert_eq!(seq, b"CCGT");

        reader.read("id", 13, 23, &mut seq).unwrap();
        assert_eq!(seq, b"CGTAGGCTGA");
    }

    #[test]
    fn test_indexed_reader_no_trailing_lf() {
        let mut reader = IndexedReader::new(io::Cursor::new(FASTA_FILE_NO_TRAILING_LF),
                                            FAI_FILE_NO_TRAILING_LF)
                .unwrap();
        let mut seq = Vec::new();

        reader.read("id", 0, 16, &mut seq).unwrap();
        assert_eq!(seq, b"GTAGGCTGAAAACCCC");
    }

    #[test]
    fn test_indexed_reader_bad_reader() {
        let bad_reader = ReaderMock {
            seek_fails: false,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        let mut seq = Vec::new();
        assert!(reader.read("id", 0, 10, &mut seq).is_ok())
    }

    #[test]
    fn test_indexed_reader_read_seek_fails() {
        let bad_reader = ReaderMock {
            seek_fails: true,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        let mut seq = Vec::new();
        assert!(reader.read("id", 0, 10, &mut seq).is_err());
    }

    #[test]
    fn test_indexed_reader_read_read_fails() {
        let bad_reader = ReaderMock {
            seek_fails: false,
            read_fails: true,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        let mut seq = Vec::new();
        assert!(reader.read("id", 0, 10, &mut seq).is_err());
    }

    #[test]
    fn test_indexed_reader_iter_seek_fails() {
        let bad_reader = ReaderMock {
            seek_fails: true,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        assert!(reader.read_iter("id", 0, 10).is_err());
    }

    #[test]
    fn test_indexed_reader_iter_read_fails() {
        let bad_reader = ReaderMock {
            seek_fails: false,
            read_fails: true,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        let mut iterator = reader.read_iter("id", 0, 10).unwrap();
        assert!(iterator.next().unwrap().is_err());
    }

    #[test]
    fn test_writer() {
        let mut writer = Writer::new(Vec::new());
        writer
            .write("id", Some("desc"), b"ACCGTAGGCTGA")
            .unwrap();
        writer.write("id2", None, b"ATTGTTGTTTTA").unwrap();
        writer.flush().unwrap();
        assert_eq!(writer.writer.get_ref(), &WRITE_FASTA_FILE);
    }
}
