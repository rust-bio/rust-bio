// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Mathematical and statistical tools.


pub mod combinatorics;
pub mod probs;
pub mod bayesian;
pub mod pairhmm;
pub mod profilehmm;

pub use stats::probs::{Prob, LogProb, PHREDProb};
