// Copyright 2019 Johannes Köster, University of Duisburg-Essen.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Error definitions for the `pssm` module.
use thiserror::Error;

#[derive(
    Error, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub enum Error {
    #[error(
        "query length {} is shorter than motif length {}",
        query_len,
        motif_len
    )]
    QueryTooShort { motif_len: usize, query_len: usize },
    #[error("attempted to build a motif from sequences with mismatched lengths")]
    InconsistentLen,
    #[error("monomer '{}' is invalid", char::from(*mono))]
    InvalidMonomer { mono: u8 },
    #[error("motif cannot be created from zero sequences")]
    EmptyMotif,
    #[error("information-free motif: a motif in which every monomer is equally likely at every position will result in a divide-by-zero exception")]
    NullMotif,
    #[error("expected pseudo-score array of length {}; got {}", expected, received)]
    InvalidPseudos { expected: u8, received: u8 },
}

pub type Result<T, E = Error> = std::result::Result<T, E>;
