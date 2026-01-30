// Copyright 2014-2025 Johannes Köster and Markus Schlegel.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Myers bit-parallel approximate pattern matching algorithm.
//! Finds all matches up to a given edit distance k and supports ambiguous matching.
//! Methods for obtaining the alignment path are provided.
//! The simple and fast implementation works with patterns of up to
//! 64 or 128 symbols, depending on the bit vector type used.
//! The block-based implementation ([`long`] submodule) supports patterns of
//! arbitrary length.
//!
//! The time complexity is O(n) where n is the length of the text.
//! With unlimited pattern matching ([`long`]), the time additionally scales with
//! the distance threshold k: O(n * k). The search can be considered linear-time
//! if k is kept constant.
//!
//! While the search only yields the end position and edit distance of matches,
//! additional methods are provided for obtaining the start and the alignment path
//! (edit operations) of a hit, in O(m + k) worst-case time.
//! The implementation is very  similar to Edlib's (Šošić and Šikić, 2017),
//! although minor differences may occur when there are multiple possible
//! alignments with the same edit distance.
//! Edlib prioritizes InDels over substitutions (insertion > deletion > substitution),
//! whereas this implementation prioritizes substitutions (substitution > insertion > deletion).
//!
//! Myers, G. (1999). *A fast bit-vector algorithm for approximate string matching based on dynamic
//!  programming*. Journal of the ACM (JACM) 46, 395–415. <https://doi.org/10.1145/316542.316550>
//!
//! Šošić, M., and Šikić, M. (2017). *Edlib: a C/C ++ library for fast, exact sequence alignment
//! using edit distance.* Bioinformatics 33, 1394–1395. <https://doi.org/10.1093/bioinformatics/btw753>
//!
//! # Example
//!
//! Iterating over matches in pairs of `(end, distance)` using `u64` as bit vector type:
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
//! # Long patterns
//!
//! It is also possible to to construct `Myers::<u128>`, which handles patterns with
//! up to 128 symbols using the standard algorithm.
//!
//! In addition, the [`long`] submodule handles unlimited patterns by splitting
//! them into separate blocks. An example:
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
//! // the pattern of length 9 is too long to fit into a single `u8` bit vector
//! // (panics!)
//! // let myers_8 = Myers::<u8>::new(pattern);
//!
//! // However, we can use `long::Myers` with `u8` blocks
//! let myers_long_8 = long::Myers::<u8>::new(pattern);
//! let occ_long_8: Vec<_> = myers_long_8
//!     .find_all_end(text, 2)
//!     .map(|(end, dist)| (end, dist as u8))
//!     .collect();
//!
//! assert_eq!(occ_64, occ_long_8);
//! # }
//! ```
//!
//! <div class="warning">
//!
//! `u8` is used for demonstration, but `long::Myers::<u64>` is still
//! the best in most cases.
//!
//! </div>
//!
//! # Obtaining the start position of a match
//!
//! The `Myers::find_all` method provides an iterator over tuples of `(start, end, distance)`.
//! As calculating the start position requires finding the alignment path, this is
//! slower than `Myers::find_all_end`.
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
//! <div class="warning">
//!
//! The end positions here are greater by one than in `Myers::find_all_end`.
//! This is intentional, since the iterator returns a range rather than an index, and in
//! Rust, ranges do not include the end position.
//!
//! </div>
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
//!     println!(
//!         "Hit found in range: {}..{} (distance: {})",
//!         aln.ystart, aln.yend, aln.score
//!     );
//!     println!("{}", aln.pretty(pattern, text, 80));
//! }
//! # }
//! ```
//! **Output:**
//!
//! <pre>
//! Hit found in range: 3..10 (distance: 3)
//!    TCCTAGGGC
//!    ||||+|\|+
//! TCCTCCT-GAG-GGATTAGCAC
//!
//! Hit found in range: 3..11 (distance: 3)
//!    TCCTAGGGC
//!    ||||+|\|\
//! TCCTCCT-GAGGGATTAGCAC
//!
//! Hit found in range: 3..12 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||+
//! TCCTCCTGAGGG-ATTAGCAC
//!
//! Hit found in range: 3..13 (distance: 2)
//!    TCCT-AGGGC
//!    ||||x||||\
//! TCCTCCTGAGGGATTAGCAC
//!
//! ... (truncated)
//!
//! </pre>
//!
//! **Note** that the [`Alignment`](../../alignment/struct.Alignment.html) instance is created
//! only once and then reused. Since the Myers algorithm is very fast, the allocation
//! of `Alignment::operations` can have a non-negligible impact on performance.
//!
//! # Finding the best hit
//!
//! In many cases, only the match with the smallest edit distance is the only one of.
//! interest. Therefore, calculating an alignment for every hit is not necessary.
//! [`LazyMatches`](struct.LazyMatches.html) returned by `Myers::find_all_lazy()`
//! provide an iterator over tuples of `(end, distance)` as `Myers::find_all_end()`.
//! However, the information needed to calculate the alignment path later on
//! is also stored. This has a slight performance impact and requires O(n) memory,
//! as opposed to O(m + k) for `Myers::find_all()`.
//! Still, the following code is faster than using `Myers::find_all()`:
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
//! println!("{}", aln.pretty(pattern, text, 80));
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
//! Matching of multiple or all symbols at once can be configured through `MyersBuilder`.
//! In the following example, `N` matches `A`, `C`, `G` or `T`:
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
//! See the documentation of [`MyersBuilder`](struct.MyersBuilder.html) for more examples.

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
    impl_common_tests!(false, super, u64, u8, build_64);

    use std::iter::repeat_n;

    #[test]
    #[should_panic(expected = "Pattern too long")]
    fn test_pattern_too_long() {
        let pattern: Vec<_> = repeat_n(b'T', 65).collect::<Vec<u8>>();
        super::Myers::<u8>::new(pattern);
    }

    #[test]
    #[should_panic(expected = "Pattern too long")]
    fn test_pattern_too_long_builder() {
        let pattern: Vec<_> = repeat_n(b'T', 65).collect::<Vec<u8>>();
        super::MyersBuilder::new().build_64(pattern);
    }
}
