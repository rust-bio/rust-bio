#![feature(core)]
#![feature(collections)]
#![feature(test)]

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


pub mod utils;
pub mod alphabets;
pub mod pattern_matching;
pub mod data_structures;
pub mod alignment;
pub mod bitencoding;
