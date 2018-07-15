// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This module contains various useful pattern matching algorithms.
//! The implementations are based on the lecture notes
//! "Algorithmen auf Sequenzen", Kopczynski, Marschall, Martin and Rahmann, 2008 - 2015.
//!
//! * Algorithm of Horspool: fastest for a sufficiently large alphabet
//! * Shift And algorithm: fast for patterns with less than 64 symbols and very small alphabets.
//! * BNDM algorithm: fast for patterns with less than 64 symbols.
//! * BOM algorithm: fast for long patterns and small alphabet.
//! * KMP algorithm: the classical ancestor.
//! * Ukkonens algorithm: approximate pattern matching with dynamic programming.
//! * Myers algorithm: linear-time approximate pattern matching with edit distance for small patterns
//!
//! Another fast pattern matching algorithm is available in the twoway crate: https://crates.io/crates/twoway

pub mod bndm;
pub mod bom;
pub mod horspool;
pub mod kmp;
pub mod myers;
pub mod pssm;
pub mod shift_and;
pub mod ukkonen;
