// Copyright 2019 Johannes Köster, University of Duisburg-Essen.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Error definitions for the `hmm` module.

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display(
        "inferred from A: N_0={}, N_1={} (must be equal), from B: N={}, M={}, from pi: N={}",
        an0,
        an1,
        bn,
        bm,
        pin
    ))]
    InvalidDimension {
        an0: usize,
        an1: usize,
        bn: usize,
        bm: usize,
        pin: usize,
    },
}
