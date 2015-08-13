#![feature(test)]
#![feature(step_by)]
#![feature(convert)]
#![feature(iter_arith)]
#![feature(vec_push_all)]

// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! # Rust-bio, a bioinformatics library for Rust.
//! This library provides implementations of many algorithms and data structures
//! that are useful for bioinformatics.
//! All provided implementations are rigorously tested via continuous
//! integration (https://travis-ci.org/johanneskoester/rust-bio).
//! For installation instructions and a general overview, visit
//! https://github.com/rust-bio/rust-bio.
//!
//! Currently, rust-bio provides
//!
//! * most major pattern matching algorithms,
//! * a convenient alphabet implementation,
//! * pairwise alignment,
//! * suffix arrays,
//! * BWT and FM-Index,
//! * FMD-Index for finding supermaximal exact matches,
//! * a q-gram index,
//! * a rank/select data structure,
//! * FASTQ and FASTA and BED readers and writers,
//! * helper functions for combinatorics and dealing with log probabilities.
//!
//! Documentation and examples for each module can be found in the module descriptions below.


extern crate test;
extern crate rustc_serialize;
extern crate csv;
extern crate num;
extern crate itertools;
extern crate bit_vec;
extern crate vec_map;
extern crate bit_set;

pub mod utils;
pub mod alphabets;
pub mod pattern_matching;
pub mod data_structures;
pub mod alignment;
pub mod io;
pub mod stats;
