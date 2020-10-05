// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Myers bit-parallel approximate pattern matching algorithm.
//! Finds all matches up to a given edit distance. The pattern has to fit into a bitvector,
//! and is thus limited to 64 or (since stable Rust version 1.26) to 128 symbols.
//! Complexity: O(n)
//!
//! Traceback allows obtaining the starting position and the alignment path of the hit.
//! Its implementation is somehow similar to the one by Edlib (Šošić and Šikić 2017),
//! although there can be small differences when there is more than one possible alignment
//! path with then same edit distance at a position: Edlib prefers to make insertions
//! and deletions to the pattern (query) over substitutions
//! (Insertion > Deletion > Substitution) while our implementation prefers substitutions
//! (Substitution > Insertion > Deletion).
//!
//! *Myers, G. (1999). A fast bit-vector algorithm for approximate string matching based on dynamic
//!  programming. Journal of the ACM (JACM) 46, 395–415.*
//!
//! *Šošić, M., and Šikić, M. (2017). Edlib: a C/C ++ library for fast, exact sequence alignment
//! using edit distance. Bioinformatics 33, 1394–1395.*
//!
//! # Example
//!
//! Iterating over matches in pairs of `(end, distance)` using `u64` as bitvector type:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let myers = Myers::<u64>::new(pattern);
//! let occ: Vec<_> = myers.find_all_end(text, 2).collect();
//!
//! assert_eq!(occ, [(11, 2), (12, 2)]);
//! # }
//! ```
//!
//! Starting with stable Rust 1.26, it is also possible to use `u128` as bitvector
//! (`Myers::<u128>`), which enables longer patterns, but is somewhat slower.
//!
//! # Long patterns
//!
//! With the default implementation, query (pattern) length is limited by the size of the
//! bit-vector; 64 symbols for `Myers::<u64>`. Patterns longer than 128 symbols (when using
//! `u128` as bit-vector) can only be handled by using the block-based Myers implementation,
//! which lives in the [`long`](long/index.html) submodule. An example:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::{long, Myers};
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let myers_64 = Myers::<u64>::new(pattern);
//! let occ_64: Vec<_> = myers_64.find_all_end(text, 2).collect();
//!
//! // the pattern of length 9 is too long to fit into a single `u8` bit-vector
//! // (panics!)
//! // let myers_8 = Myers::<u8>::new(pattern);
//!
//! // However, we can use the block-based implementation with `u8` bit-vectors
//! let myers_long_8 = long::Myers::<u8>::new(pattern);
//! let occ_long_8: Vec<_> = myers_long_8
//!     .find_all_end(text, 2)
//!     .map(|(end, dist)| (end, dist as u8))
//!     .collect();
//!
//! assert_eq!(occ_64, occ_long_8);
//! # }
//! ```
//! Note that `u8` just used for demonstration, using `u64` is still the best in most cases.
//!
//! # Obtaining the starting position of a match
//!
//! The `Myers::find_all` method provides an iterator over tuples of `(start, end, distance)`.
//! Calculating the starting position requires finding the alignment path, therefore this is
//! slower than `Myers::find_all_end`. Note that the end positions differ from above by one.
//! This is intentional, as the iterator returns a range rather an index, and ranges in Rust
//! do not include the end position by default.
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::Myers;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers::<u64>::new(pattern);
//! let occ: Vec<_> = myers.find_all(text, 2).collect();
//!
//! assert_eq!(occ, [(3, 12, 2), (3, 13, 2)]);
//! # }
//! ```
//!
//! # Obtaining alignments
//!
//! [`FullMatches`](struct.FullMatches.html) returned by `Myers::find_all()` also provide a method
//! for obtaining an alignment path:
//!
//! ```
//! # extern crate bio;
//! use bio::alignment::Alignment;
//! use bio::pattern_matching::myers::Myers;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers::<u64>::new(pattern);
//! // create an 'empty' alignment instance, which can be reused
//! let mut aln = Alignment::default();
//!
//! let mut matches = myers.find_all(text, 3);
//! while matches.next_alignment(&mut aln) {
//!     //println!("Hit fond in range: {}..{} (distance: {})", aln.ystart, aln.yend, aln.score);
//!     //println!("{}", aln.pretty(pattern, text));
//! }
//! # }
//! ```
//! **Output:**
//!
//! <pre>
//! Hit fond in range: 3..10 (distance: 3)
//!    TCCTAGGGC
//!    ||||+|\|+
//! TCCTCCT-GAG-GGATTAGCAC
//!
//! Hit fond in range: 3..11 (distance: 3)
//!    TCCTAGGGC
//!    ||||+|\|\
//! TCCTCCT-GAGGGATTAGCAC
//!
//! Hit fond in range: 3..12 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||+
//! TCCTCCTGAGGG-ATTAGCAC
//!
//! Hit fond in range: 3..13 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||\
//! TCCTCCTGAGGGATTAGCAC
//!
//! ... (truncated)
//!
//! </pre>
//!
//! **Note** that the [`Alignment`](../../alignment/struct.Alignment.html) instance is only created
//! once and then reused. Because the Myers algorithm is very fast, the allocation necessary for
//! `Alignment::operations` can have a non-negligible impact on performance; and thus, recycling
//! makes sense.
//!
//! # Finding the best hit
//!
//! In many cases, only the match with the smallest edit distance is actually of interest.
//! Calculating an alignment for every hit is therefore not necessary.
//! [`LazyMatches`](struct.LazyMatches.html) returned by `Myers::find_all_lazy()`
//! provide an iterator over tuples of `(end, distance)` like `Myers::find_all_end()`, but
//! additionally keep the data necessary for calculating the alignment path later at any desired
//! position. Storing the data itself has a slight performance impact and requires more memory
//! compared to `Myers::find_all_end()` [O(n) as opposed to O(m + k)]. Still the following code
//! is faster than using `FullMatches`:
//!
//! ```
//! # extern crate bio;
//! use bio::alignment::Alignment;
//! use bio::pattern_matching::myers::Myers;
//!
//! # fn main() {
//! let text = b"CGGTCCTGAGGGATTAGCAC";
//! let pattern = b"TCCTAGGGC";
//!
//! let mut myers = Myers::<u64>::new(pattern);
//! let mut aln = Alignment::default();
//!
//! let mut matches = myers.find_all_lazy(text, 2);
//!
//! // first, find the best hit
//! let (best_end, _) = matches.by_ref().min_by_key(|&(_, dist)| dist).unwrap();
//!
//! // now calculate the alignment
//! matches.alignment_at(best_end, &mut aln);
//! println!(
//!     "Best alignment at {}..{} (distance: {})",
//!     aln.ystart, aln.yend, aln.score
//! );
//! println!("{}", aln.pretty(pattern, text));
//! # }
//! ```
//!
//! **Output:**
//!
//! <pre>
//! Best alignment at 3..12 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||+
//! TCCTCCTGAGGG-ATTAGCAC
//! </pre>
//!
//! Actually as seen in the previous chapters, there are two hits with the same distance of 2.
//! It may make sense to consider both of them.
//!
//! # Dealing with ambiguities
//!
//! Matching multiple or all symbols at once can be achieved using `MyersBuilder`. This example
//! allows `N` in the search pattern to match all four DNA bases in the text:
//!
//! ```
//! # extern crate bio;
//! use bio::pattern_matching::myers::MyersBuilder;
//!
//! # fn main() {
//! let text = b"GTCTGATCTTACC";
//! let pattern = b"TGATCNT";
//!
//! let myers = MyersBuilder::new().ambig(b'N', b"ACGT").build_64(pattern);
//! assert_eq!(myers.distance(text), 0);
//! # }
//! ```
//!
//! For more examples see the documentation of [`MyersBuilder`](struct.MyersBuilder.html).

#[macro_use]
mod myers_impl;
mod builder;
mod helpers;
#[cfg(test)]
#[macro_use]
pub(crate) mod common_tests;
pub mod long;
mod simple;
mod traceback;

pub use self::builder::MyersBuilder;
pub use self::helpers::*;
use self::myers_impl::*;
pub use self::simple::*;

#[cfg(test)]
mod tests {
    // from common_tests.rs
    impl_tests!(super, u64, u8, build_64);

    use std::iter::repeat;

    #[test]
    #[should_panic(expected = "Pattern too long")]
    fn test_pattern_too_long() {
        let pattern: Vec<_> = repeat(b'T').take(65).collect();
        super::Myers::<u8>::new(&pattern);
    }

    #[test]
    #[should_panic(expected = "Pattern too long")]
    fn test_pattern_too_long_builder() {
        let pattern: Vec<_> = repeat(b'T').take(65).collect();
        super::MyersBuilder::new().build_64(&pattern);
    }
}
