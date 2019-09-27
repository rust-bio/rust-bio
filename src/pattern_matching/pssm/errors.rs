// Copyright 2019 Johannes Köster, University of Duisburg-Essen.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Error definitions for the `pssm` module.

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display(
        "query length {} is shorter than motif length {}",
        query_len,
        motif_len
    ))]
    QueryTooShort { motif_len: usize, query_len: usize },
    #[snafu(display("attempted to build a motif from sequences with mismatched lengths"))]
    InconsistentLen,
    #[snafu(display("monomer '{}' is invalid", char::from(*mono)))]
    InvalidMonomer { mono: u8 },
    #[snafu(display("motif cannot be created from zero sequences"))]
    EmptyMotif,
    #[snafu(display("information-free motif: a motif in which every monomer is equally likely at every position will result in a divide-by-zero exception"))]
    NullMotif,
    #[snafu(display("expected pseudo-score array of length {}; got {}", expected, received))]
    InvalidPseudos { expected: u8, received: u8 },
}
