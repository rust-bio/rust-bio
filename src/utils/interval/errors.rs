// Copyright 2019 Johannes Köster, University of Duisburg-Essen.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Error definitions for the `interval` module.

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display("an Interval must have a Range with a positive width"))]
    InvalidRange,
}
