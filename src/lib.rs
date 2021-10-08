// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![doc(
    html_logo_url = "https://raw.githubusercontent.com/rust-bio/rust-bio/master/img/bioferris.svg",
    html_favicon_url = "https://raw.githubusercontent.com/rust-bio/rust-bio/master/img/bioferris.svg"
)]

//! # Rust-bio, a bioinformatics library for Rust.
//! This library provides implementations of many algorithms and data structures
//! that are useful for bioinformatics.
//! All provided implementations are rigorously tested via continuous
//! integration.
//!
//! For **getting started** with using `rust-bio`, see [the `Getting started` section below](#getting-started).
//! For navigating the documentation of the available modules, see [the `Modules` section below](#modules).
//! If you want to contribute to `rust-bio`, see [`CONTRIBUTING.md`](https://github.com/rust-bio/rust-bio/CONTRIBUTING.md).
//!
//! Currently, rust-bio provides
//!
//! * most major pattern matching algorithms,
//! * a convenient alphabet implementation,
//! * pairwise alignment,
//! * suffix arrays,
//! * the [Burrows-Wheeler-transform (BWT)](https://en.wikipedia.org/wiki/Burrows–Wheeler_transform)
//! * the [Full-text index in Minute space index (FM-index)](https://doi.org/10.1109/SFCS.2000.892127),
//! * FMD-Index for finding supermaximal exact matches,
//! * a q-gram index,
//! * utilities to work with [PSSMs](https://en.wikipedia.org/wiki/Position_weight_matrix),
//! * an open reading frame (ORF) search algorithm,
//! * a rank/select data structure,
//! * [serde](https://github.com/serde-rs/serde) support for all data structures when built with `nightly` feature,
//! * readers and writers for FASTQ, FASTA and BED,
//! * helper functions for combinatorics and dealing with log probabilities,
//! * an implementation of the Hidden Markov Model and related algorithms.
//!
//! For reading and writing SAM/BAM/CRAM, VCF/BCF files or tabix indexed files, have a look at [rust-htslib](https://docs.rs/rust-htslib).
//!
//! # Getting started
//!
//! We explain how to use Rust-Bio step-by-step.
//! Users who already have experience with Rust can skip right to [Step 3: Use Rust-Bio in your project](#step-3-use-rust-bio-in-your-project).
//! Users who already know `rust-bio` might want to jump right into the [modules docs](#modules)
//!
//! ## Step 1: Setting up Rust
//!
//! Rust can be installed following the instruction for [rustup](https://rustup.rs/).
//!
//!
//! ## Step 2: Setting up a new Rust project
//!
//! Since Rust-Bio is a library, you need to setup your own new Rust project to use Rust-Bio.
//! With Rust, projects and their dependencies are managed with the builtin package manager [Cargo](https://doc.rust-lang.org/cargo/index.html).
//! To create a new Rust project, issue
//!
//! ```bash
//! cargo new hello_world --bin
//! cd hello_world
//! ```
//! in your terminal. The flag `--bin` tells Cargo to create an executable project instead of a library.
//! In [this section](https://doc.rust-lang.org/nightly/book/hello-cargo.html#a-new-project) of the Rust docs, you find details about what Cargo just created for you.
//!
//! Your new project can be compiled with
//! ```bash
//! cargo build
//! ```
//! If dependencies in your project are out of date, update with
//! ```bash
//! cargo update
//! ```
//! Execute the compiled code with
//! ```bash
//! cargo run
//! ```
//! If you are new to Rust, we suggest to proceed with [learning Rust](https://www.rust-lang.org/learn) via the Rust docs.
//!
//! ## Step 3: Use Rust-Bio in your project
//!
//! To use Rust-Bio in your Rust project, add the following to your `Cargo.toml`
//!
//! ```toml
//! [dependencies]
//! bio = "*"
//! ```
//!
//! and import the crate from your source code (if your `cargo --version` is greater than 1.31, you can skip this step):
//!
//! ```rust
//! extern crate bio;
//! ```
//!
//! # Examples
//!
//! - [FM-index and FASTQ](https://github.com/rust-bio/rust-bio/examples/fmindex_fastq.rs)
//! - [Multithreaded](https://github.com/rust-bio/rust-bio/examples/multithreaded.rs)
//!
//! Documentation and further examples for each module can be found in the module descriptions below.
//!
//! # Benchmarks
//!
//! Since Rust-Bio is based on a compiled language, similar performance to C/C++ based libraries can be expected. Indeed, we find the pattern matching algorithms of Rust-Bio to perform in the range of the C++ library Seqan:
//!
//! | Algorithm | Rust-Bio | Seqan   |
//! | --------- | -------: | ------: |
//! | BNDM      | 77ms     | 80ms    |
//! | Horspool  | 122ms    | 125ms   |
//! | BOM       | 103ms    | 107ms   |
//! | Shift-And | 241ms    | 545ms   |
//!
//! We measured 10000 iterations of searching pattern `GCGCGTACACACCGCCCG` in the sequence of the hg38 MT chromosome.
//! Initialization time of each algorithm for the given pattern was included in each iteration. Benchmarks were conducted with *Cargo bench* for Rust-Bio and *Python timeit* for Seqan on an Intel Core i5-3427U CPU.
//! Benchmarking Seqan from *Python timeit* entails an overhead of 1.46ms for calling a C++ binary. This overhead was subtracted from above Seqan run times.
//! Note that this benchmark only compares the two libraries to exemplify that Rust-Bio has comparable speed to C++ libraries: all used algorithms have their advantages for specific text and pattern structures and lengths (see [the pattern matching section in the documentation](https://docs.rs/bio/0.28.2/bio/pattern_matching/index.html))./!

#[macro_use]
extern crate approx;

#[macro_use]
extern crate custom_derive;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate newtype_derive;

#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate strum_macros;

#[macro_use]
extern crate getset;

#[cfg(feature = "phylogeny")]
#[macro_use]
extern crate pest_derive;

pub mod alignment;
pub mod alphabets;
pub mod data_structures;
pub mod io;
pub mod pattern_matching;
pub mod scores;
pub mod seq_analysis;
pub mod stats;
pub mod utils;
