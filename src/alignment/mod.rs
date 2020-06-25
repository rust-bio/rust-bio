// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various alignment and distance computing algorithms.

pub mod distance;
pub mod pairwise;
pub mod poa;
pub mod sparse;

// Re-export the alignment types.
pub use bio_types::alignment::*;
