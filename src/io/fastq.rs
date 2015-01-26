use std::io;

use utils::trim_newline;


/// A FASTQ record, returned by the FASTQ parser.
pub struct Record {
    pub name: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>
}


/// An iterable FASTQ parser.
///
/// # Example
///
/// ```rust
/// use bio::io::fastq::FastqFile;
/// use std::io::{BufferedReader, stdin};
///
/// let mut buffer = BufferedReader::new(stdin());
/// let fastq = FastqFile::new(buffer);
/// ```
pub struct FastqFile<T> {
    buffer: T
}


impl<T> FastqFile<T> where T: io::Buffer {
    /// Create a new FASTQ parser object.
    ///
    /// # Arguments
    ///
    /// * `buffer` - a buffer object (e.g. the io::BufferedReader over a file or STDIN)
    pub fn new(buffer: T) -> Self {
        FastqFile { buffer: buffer }
    }
}


/// Iterator over the FASTQ file.
impl<T> Iterator for FastqFile<T> where T: io::Buffer {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        let mut name;
        match self.buffer.read_line() {
            Ok(line) => {
                name = line;
            },
            Err(io::IoError { kind: io::IoErrorKind::EndOfFile, .. } ) => return None,
            Err(e) => panic!(e.desc)
        }

        if !name.starts_with("@") {
            panic!("Expecting record to start with '>'.")
        }
        name.remove(0);

        let mut seq = self.buffer.read_line().ok().expect("Incomplete record");

        let sep = self.buffer.read_line().ok().expect("Incomplete record");
        if !sep.starts_with("+") {
            panic!("Expecting '+' separator between sequence and quals.")
        }

        let mut qual = self.buffer.read_line().ok().expect("Incomplete record");
        trim_newline(&mut name);
        trim_newline(&mut seq);
        trim_newline(&mut qual);

        Some(Record { name: name, seq: seq.into_bytes(), qual: qual.into_bytes() })
    }
}
