use std::io;


/// A FASTQ record, returned by the FASTQ parser.
pub struct Record {
    name: String,
    seq: String,
    qual: String
}


/// An iterable FASTQ parser.
pub struct FastqFile<T> {
    buffer: T
}


/// Create a new FASTQ parser object.
///
/// # Arguments
///
/// * `buffer` - a buffer object (e.g. the io::BufferedReader over a file or STDIN)
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
impl<T> FastqFile<T> where T: io::Buffer {
    pub fn new(buffer: T) -> Self {
        FastqFile { buffer: buffer }
    }
}


/// Iterator over the FASTQ file.
impl<T> Iterator for FastqFile<T> where T: io::Buffer {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        let name;
        match self.buffer.read_line() {
            Ok(line) => {
                name = line.slice(1, line.len()).trim().to_string();
            },
            Err(io::IoError { kind: io::IoErrorKind::EndOfFile, .. } ) => return None,
            Err(e) => panic!(e.desc)
        }

        let seq = self.buffer.read_line().ok().expect("Incomplete record").trim().to_string();
        self.buffer.read_line().ok().expect("Incomplete record");
        let qual = self.buffer.read_line().ok().expect("Incomplete record").trim().to_string();

        Some(Record { name: name, seq: seq, qual: qual })
    }
}
