// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! # Rust-bio, a bioinformatics library for Rust.
//! This library provides implementations of many algorithms and data structures
//! that are useful for bioinformatics.
//! All provided implementations are rigorously tested via continuous
//! integration.
//! For installation instructions and a general overview, visit
//! https://rust-bio.github.io.
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
//! * an orf research algorithm,
//! * a rank/select data structure,
//! * FASTQ and FASTA and BED readers and writers,
//! * helper functions for combinatorics and dealing with log probabilities.
//!
//! # Example
//!
//! ```rust
//! use bio::alphabets;
//! use bio::data_structures::suffix_array::{suffix_array, SampleableSuffixArray};
//! use bio::data_structures::bwt::{bwt, less, Occ};
//! use bio::data_structures::fmindex::{fmindex_raw, fmindex_sampled, FMIndexable};
//!
//! let text = b"ACGGATGCTGGATCGGATCGCGCTAGCTA$";
//! let pattern = b"ACCG";
//!
//! // Create an FM-Index for a given text.
//! let alphabet = alphabets::dna::iupac_alphabet();
//! let sa = suffix_array(text);
//! let bw = bwt(text, &sa);
//! let lss = less(&bw, &alphabet);
//! let occ = Occ::new(&bw, 3, &alphabet);
//! let ssa = sa.sample(bw, lss, occ, 32);
//! let fm_sampled = fmindex_sampled(ssa);
//!
//! // Store the suffix array interval from the search, then find positions in the text
//! let interval = fm_sampled.backward_search(pattern.iter());
//! let pos_sampled = interval.occ(&sa);
//!
//! // FMIndexable also provides a convenience method for finding positions
//! let pos_sampled2 = fm_sampled.offsets(pattern.iter());
//! assert_eq!(pos_sampled, pos_sampled2);
//!
//! // You can also trade memory usage for raw speed with the raw suffix array, as well
//! let bw = bwt(text, &sa);
//! let occ = Occ::new(&bw, 3, &alphabet);
//! let lss = less(&bw, &alphabet);
//! let fm_raw = fmindex_raw(sa, bw, occ, lss);
//!
//! // Store the suffix array interval from the search, then find positions in the text
//! let interval = fm_raw.backward_search(pattern.iter());
//! let pos_raw = interval.occ(fm_raw.sa());
//!
//! // FMIndexable also provides a convenience method for finding positions
//! let pos_raw2 = fm_raw.offsets(pattern.iter());
//! assert_eq!(pos_raw, pos_raw2);
//! assert_eq!(pos_raw, pos_sampled);
//! ```
//!
//! Documentation and further examples for each module can be found in the module descriptions below.

#![cfg_attr(feature = "serde_macros", feature(const_fn, custom_derive, plugin))]
#![cfg_attr(feature = "serde_macros", plugin(serde_macros))]

#[cfg(feature = "serde_macros")]
extern crate serde;

extern crate rustc_serialize;
extern crate csv;
extern crate num;
extern crate itertools;
extern crate bit_vec;
extern crate vec_map;
extern crate bit_set;
#[macro_use]
extern crate lazy_static;
extern crate nalgebra;
#[macro_use]
extern crate approx;

pub mod utils;
pub mod alphabets;
pub mod pattern_matching;
pub mod data_structures;
pub mod alignment;
pub mod io;
pub mod seq_analysis;
pub mod stats;
pub mod scores;
