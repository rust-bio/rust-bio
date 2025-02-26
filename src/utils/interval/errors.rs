// Copyright 2019 Johannes Köster, University of Duisburg-Essen.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Error definitions for the `interval` module.
use thiserror::Error;

#[derive(
    Error, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub enum Error {
    #[error("an Interval must have a Range with a positive width")]
    InvalidRange,
}
pub type Result<T, E = Error> = std::result::Result<T, E>;
