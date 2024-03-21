//! A struct to read distance matrices in the [Phylip `.dist` format](https://mothur.org/wiki/phylip-formatted_distance_matrix/).
//!
//!  # Example
//!
//!  In this example, we parse a square distance matrix from a string.
//!
//!  ```
//!  use bio::io::phylip;
//!
//!  let distances = phylip::from_string("2
//!  A 0 1
//!  B 1 0
//!  ").unwrap();
//!  assert_eq!(distances.matrix_type, bio_types::distancematrix::MatrixType::Square)
//!  ```
//!
//!  In this example, we parse a lower triangular distance matrix from a string.
//!
//!  ```
//!  use bio::io::phylip;
//!
//!  let distances = phylip::from_string("3
//!  A
//!  B 1
//!  C 2 3
//!  ").unwrap();
//!  assert_eq!(distances.matrix_type, bio_types::distancematrix::MatrixType::Lower)
//!  ```

use bio_types::distancematrix::DistanceMatrix;
use itertools::zip;
use std::{
    fs, io,
    num::ParseIntError,
    path::{Path, PathBuf},
};
use thiserror::Error;

/// A `thiserror` error type gathering all the potential bad outcomes
#[derive(Debug, Error)]
pub enum Error {
    #[error("Error while opening {}: {}", filename.display(), source)]
    OpenFile {
        filename: PathBuf,
        source: std::io::Error,
    },

    #[error("Error while reading distance matrix: {0}")]
    Read(#[from] std::io::Error),

    #[error("Error while writing distance matrix: {0}")]
    Write(std::io::Error),

    #[error("DistanceMatrix contains invalid UTF-8: {0}")]
    InvalidContent(#[from] std::str::Utf8Error),

    #[error("Error while parsing distance matrix: {0}")]
    ParsingError(String),
}
type Result<T, E = Error> = std::result::Result<T, E>;

/// Reads a distance matrix from an `&str`-compatible type
pub fn from_string<S: AsRef<str>>(content: S) -> Result<DistanceMatrix> {
    let mut lines = content.as_ref().lines();
    let n = lines
        .next()
        .ok_or(Error::ParsingError("Expected n, number of taxons".into()))?
        .parse()
        .map_err(|e: ParseIntError| Error::ParsingError(e.to_string()))?;

    let mut names = Vec::with_capacity(n);
    let mut distances = Vec::with_capacity(n);

    for line in lines.into_iter() {
        let mut it = line.split_ascii_whitespace();
        names.push(
            it.next()
                .ok_or(Error::ParsingError(
                    "Line does not start with a name".into(),
                ))?
                .to_string(),
        );
        distances.push(it.map(|chars| chars.parse().unwrap()).collect());
    }

    DistanceMatrix::new(names, distances).map_err(|e| Error::ParsingError(e.to_string()))
}

/// Reads a DistanceMatrix from a file
pub fn from_file<P: AsRef<Path>>(path: P) -> Result<DistanceMatrix> {
    fs::File::open(&path)
        .map(read)
        .map_err(|e| Error::OpenFile {
            filename: path.as_ref().to_owned(),
            source: e,
        })?
}

/// Reads a DistanceMatrix from any type implementing `io::Read`
pub fn read<R: io::Read>(reader: R) -> Result<DistanceMatrix> {
    let content_bytes = reader
        .bytes()
        .collect::<Result<Vec<_>, _>>()
        .map_err(Error::Read)?;
    let content_str = std::str::from_utf8(&content_bytes).map_err(Error::InvalidContent)?;
    from_string(&content_str)
}

/// Convert the DistanceMatrix to the phylip format.
pub fn to_string(m: &DistanceMatrix) -> String {
    let mut s = String::new();
    s += &m.len().to_string();
    s += "\n";
    for (name, ds) in zip(&m.names, &m.distances) {
        s += &name;
        for d in ds {
            s += " ";
            s += &d.to_string();
        }
        s += "\n";
    }
    s
}

/// Writes a distance matrix to a file.
pub fn to_file<P: AsRef<Path>>(path: P, m: &DistanceMatrix) -> Result<()> {
    fs::File::open(&path)
        .map(|w| write(w, m))
        .map_err(|e| Error::OpenFile {
            filename: path.as_ref().to_owned(),
            source: e,
        })?
}

/// Writes a distance matrix to any type implementing `io::Write`.
pub fn write<W: io::Write>(mut writer: W, m: &DistanceMatrix) -> Result<()> {
    let s = to_string(&m);
    writer.write_all(&s.into_bytes()).map_err(Error::Write)
}

#[cfg(test)]
mod tests {
    use bio_types::distancematrix::MatrixType;

    use super::*;

    #[test]
    fn distance_matrix_from_to_string() {
        let square = "3
a 0 1 2
B 1 0 3
XYZ 2 3 0
";
        let lower = "3
a
B 1
XYZ 2 3
";
        let upper = "3
a 1 2
B 3
XYZ
";
        for (s, matrix_type) in zip(
            [square, lower, upper],
            [MatrixType::Square, MatrixType::Lower, MatrixType::Upper],
        ) {
            let t = from_string(s).unwrap();
            assert_eq!(t.matrix_type, matrix_type);
            assert_eq!(to_string(&t), s);
        }
    }

    #[test]
    fn floats_and_names() {
        let s = "1
some_long_name?! 1.0
";
        let t = from_string(s).unwrap();
        assert_eq!(t[(0, 0)], 1.0);
        assert_eq!(t.names[0], "some_long_name?!");

        let s = "1
ABCabc 1.0001
";
        let t = from_string(s).unwrap();
        assert_eq!(t[(0, 0)], 1.0001);
        assert_eq!(t.names[0], "ABCabc");
    }
}
