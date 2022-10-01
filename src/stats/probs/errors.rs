// Copyright 2019 Johannes Köster, University of Duisburg-Essen.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Error definitions for the `probs` module.
use thiserror::Error;

#[derive(Error, Copy, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub enum Error {
    #[error("probabilty {} not in interval [0,1]", prob)]
    InvalidProb { prob: f64 },
}
pub type Result<T, E = Error> = std::result::Result<T, E>;
