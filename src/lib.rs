#![feature(core)]
#![feature(collections)]
#![feature(test)]
#![feature(io)]
#![feature(str_words)]

// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! # Rust-bio, a rusty bioinformatics library.
//!
//! This library provides implementations of many algorithms and data structures
//! that are useful for bioinformatics.
//! All provided implementations are rigorously tested via continuous
//! integration (https://travis-ci.org/johanneskoester/rust-bio).
//!
//! Currently, rust-bio provides
//!
//! * pattern matching,
//! * pairwise alignment,
//! * suffix arrays,
//! * BWT and FM-Index,
//! * rank/select data structures.
//!
//! Find out more at https://github.com/johanneskoester/rust-bio.


extern crate test;
extern crate "rustc-serialize" as rustc_serialize;
extern crate csv;

pub mod utils;
pub mod alphabets;
pub mod pattern_matching;
pub mod data_structures;
pub mod alignment;
pub mod io;
