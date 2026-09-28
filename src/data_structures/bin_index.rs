//! UCSC-style interval binning.
//!
//! The [UCSC binning scheme](http://genomewiki.ucsc.edu/index.php/Bin_indexing_system)
//! (Kent et al., *The Human Genome Browser at UCSC*, Genome Res. 2002)
//! accelerates overlap queries over genomic intervals. The genome is tiled by a
//! hierarchy of bins of decreasing size (128 kb, 1 Mb, 8 Mb, 64 Mb, 512 Mb, and,
//! in the extended scheme, 4 Gb). An interval is stored in the smallest bin that
//! fully contains it ([`assign_bin`]); to find every interval overlapping a query
//! range it is then enough to inspect the handful of bins returned by
//! [`overlapping_bins`], instead of scanning all stored intervals.
//!
//! Intervals are 0-based and half-open, i.e. `[start, end)`. Coordinates are
//! `u32` and must satisfy `start < end <= `[`MAX_END`] (2 Gb − 1, the largest
//! coordinate the reference implementation supports). The *standard* scheme is
//! used for `end <= 2^29` (512 Mb) and the *extended* scheme above it; the two
//! share a single bin-number space (extended bins are offset by 4681) so a data
//! set may freely mix small and large intervals, exactly as in the UCSC genome
//! browser database tables.
//!
//! # Example
//!
//! ```
//! use bio::data_structures::bin_index::{assign_bin, overlapping_bins};
//!
//! // The smallest bin fully containing the interval [0, 100_000).
//! let bin = assign_bin(0, 100_000);
//! assert_eq!(bin, 585);
//!
//! // Any interval overlapping [50_000, 60_000) must live in one of these bins,
//! // and the containing bin of our feature is among them.
//! let bins = overlapping_bins(50_000, 60_000);
//! assert!(bins.contains(&bin));
//! ```

/// Bin offsets of the standard scheme, one per level from finest (128 kb) to
/// coarsest (512 Mb): `512+64+8+1, 64+8+1, 8+1, 1, 0`.
const BIN_OFFSETS: [u32; 5] = [512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0];

/// Bin offsets of the extended scheme, adding a 4 Gb level on top:
/// `4096+512+64+8+1, ...`.
const BIN_OFFSETS_EXTENDED: [u32; 6] = [
    4096 + 512 + 64 + 8 + 1,
    512 + 64 + 8 + 1,
    64 + 8 + 1,
    8 + 1,
    1,
    0,
];

/// Shift to reach the finest (128 kb) bin level.
const SHIFT_FIRST: u32 = 17;

/// Shift to move from one bin level to the next coarser one.
const SHIFT_NEXT: u32 = 3;

/// Offset separating standard bin numbers (`0..=4680`) from extended ones, so
/// that both schemes coexist in a single bin-number space.
const BIN_OFFSET_OLD_TO_EXTENDED: u32 = 4681;

/// Largest `end` handled by the standard scheme: 512 Mb (2^29).
const MAX_END_512M: u32 = 512 * 1024 * 1024;

/// Largest coordinate the scheme supports: 2 Gb − 1 (2^31 − 1), the limit
/// imposed by the signed 32-bit coordinates of the reference implementation.
pub const MAX_END: u32 = i32::MAX as u32;

/// Validate that `[start, end)` is a non-empty interval within the supported
/// coordinate range, before any arithmetic relies on it.
#[inline]
fn check_range(start: u32, end: u32) {
    assert!(
        start < end,
        "binning requires a non-empty interval start < end (got start={}, end={})",
        start,
        end
    );
    assert!(
        end <= MAX_END,
        "binning requires end <= {} (2 Gb - 1); got end={}",
        MAX_END,
        end
    );
}

/// Return the smallest bin that fully contains the 0-based, half-open interval
/// `[start, end)`.
///
/// Uses the standard scheme for `end <= 2^29` and the extended scheme above it,
/// matching the UCSC `binFromRange` function. For `end <= 2^29` the returned bin
/// number is identical to the one computed by the UCSC genome browser and by the
/// `interval-binning` reference implementations.
///
/// # Panics
///
/// Panics unless `start < end <= `[`MAX_END`] (the interval must be non-empty
/// and within the supported coordinate range).
///
/// # Example
///
/// ```
/// use bio::data_structures::bin_index::assign_bin;
///
/// assert_eq!(assign_bin(0, 1), 585);
/// assert_eq!(assign_bin(0, 131_073), 73);
/// ```
pub fn assign_bin(start: u32, end: u32) -> u32 {
    check_range(start, end);
    if end <= MAX_END_512M {
        smallest_bin(start, end, &BIN_OFFSETS, 0)
    } else {
        smallest_bin(
            start,
            end,
            &BIN_OFFSETS_EXTENDED,
            BIN_OFFSET_OLD_TO_EXTENDED,
        )
    }
}

/// Walk the bin levels from finest to coarsest and return the first bin that
/// contains the whole interval, i.e. where the start and end fall in the same
/// bin at that level. `base` offsets the result into the extended number space.
///
/// The shift is applied iteratively (`>>= SHIFT_NEXT`) rather than as a single
/// `>> (SHIFT_FIRST + SHIFT_NEXT * level)`: the latter would shift a `u32` by 32
/// at the coarsest extended level, which is undefined.
fn smallest_bin(start: u32, end: u32, offsets: &[u32], base: u32) -> u32 {
    // `end > start >= 0`, so `end - 1` does not underflow.
    let mut start_bin = start >> SHIFT_FIRST;
    let mut end_bin = (end - 1) >> SHIFT_FIRST;
    for &offset in offsets {
        if start_bin == end_bin {
            return base + offset + start_bin;
        }
        start_bin >>= SHIFT_NEXT;
        end_bin >>= SHIFT_NEXT;
    }
    // Unreachable: a validated interval always collapses into a single bin at
    // the coarsest level (where both bin indices are 0).
    unreachable!("interval within supported bounds always fits in a bin")
}

/// Return every bin that may hold an interval overlapping the 0-based, half-open
/// query `[start, end)`.
///
/// The bins are returned finest level first, matching the order of the UCSC
/// `hAddBinToQueryGeneral` query builder. To find all stored intervals
/// overlapping the query, inspect exactly these bins: every interval whose
/// [`assign_bin`] falls in this list and which truly overlaps `[start, end)` is
/// covered. In particular, `overlapping_bins(s, e)` always contains
/// `assign_bin(s, e)`.
///
/// # Panics
///
/// Panics unless `start < end <= `[`MAX_END`].
///
/// # Example
///
/// ```
/// use bio::data_structures::bin_index::overlapping_bins;
///
/// assert_eq!(overlapping_bins(0, 1), vec![585, 73, 9, 1, 0, 4681]);
/// ```
pub fn overlapping_bins(start: u32, end: u32) -> Vec<u32> {
    check_range(start, end);
    let mut bins = Vec::new();
    if end <= MAX_END_512M {
        push_standard(start, end, &mut bins, true);
    } else {
        // A query that reaches into the extended range must still inspect the
        // standard bins covering its lower part, then the extended bins.
        if start < MAX_END_512M {
            push_standard(start, MAX_END_512M, &mut bins, false);
        }
        push_levels(
            start,
            end,
            &BIN_OFFSETS_EXTENDED,
            BIN_OFFSET_OLD_TO_EXTENDED,
            &mut bins,
        );
    }
    bins
}

/// Append the standard-scheme bins overlapping `[start, end)`. When
/// `self_contained` is set, also append the single extended top-level bin
/// (`4681`), the catch-all that may hold a very large interval — matching the
/// `or bin=4681` clause of the UCSC standard query.
fn push_standard(start: u32, end: u32, bins: &mut Vec<u32>, self_contained: bool) {
    push_levels(start, end, &BIN_OFFSETS, 0, bins);
    if self_contained {
        bins.push(BIN_OFFSET_OLD_TO_EXTENDED);
    }
}

/// Append, level by level from finest to coarsest, the bins spanned by
/// `[start, end)`. `base` offsets the bins into the extended number space.
fn push_levels(start: u32, end: u32, offsets: &[u32], base: u32, bins: &mut Vec<u32>) {
    // `end > start >= 0`, so `end - 1` does not underflow.
    let mut start_bin = start >> SHIFT_FIRST;
    let mut end_bin = (end - 1) >> SHIFT_FIRST;
    for &offset in offsets {
        for bin in (base + offset + start_bin)..=(base + offset + end_bin) {
            bins.push(bin);
        }
        start_bin >>= SHIFT_NEXT;
        end_bin >>= SHIFT_NEXT;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values obtained from the UCSC reference implementation
    // (kent binRange.c / hdb.c) and cross-checked against the Python
    // `interval-binning` module. See https://github.com/rust-bio/rust-bio/issues/647.
    #[test]
    fn assign_bin_reference_values() {
        // Standard scheme (end <= 2^29).
        assert_eq!(assign_bin(0, 1), 585);
        assert_eq!(assign_bin(0, 1 << 17), 585);
        assert_eq!(assign_bin(0, (1 << 17) + 1), 73);
        assert_eq!(assign_bin(1 << 17, 1 << 18), 586);
        assert_eq!(assign_bin(1_048_575, 1_048_577), 9);
        assert_eq!(assign_bin((1 << 29) - 1, 1 << 29), 4680); // last standard bin
                                                              // Extended scheme (end > 2^29).
        assert_eq!(assign_bin(1 << 29, (1 << 29) + 1), 13458); // first extended
        assert_eq!(assign_bin(1 << 30, (1 << 30) + 1), 17554);
        assert_eq!(assign_bin(MAX_END - 1, MAX_END), 25745); // largest bin
        assert_eq!(assign_bin(0, MAX_END), 4681); // spans everything
    }

    #[test]
    fn overlapping_bins_reference_values() {
        assert_eq!(overlapping_bins(0, 1), vec![585, 73, 9, 1, 0, 4681]);
        assert_eq!(
            overlapping_bins(0, 131_073),
            vec![585, 586, 73, 9, 1, 0, 4681]
        );
        assert_eq!(
            overlapping_bins(1_048_575, 1_048_577),
            vec![592, 593, 73, 74, 9, 1, 0, 4681]
        );
        // Last standard interval.
        assert_eq!(
            overlapping_bins((1 << 29) - 1, 1 << 29),
            vec![4680, 584, 72, 8, 0, 4681]
        );
        // Purely extended query (start == 2^29, no standard block).
        assert_eq!(
            overlapping_bins(1 << 29, (1 << 29) + 1),
            vec![13458, 5778, 4818, 4698, 4683, 4681]
        );
        assert_eq!(
            overlapping_bins(MAX_END - 1, MAX_END),
            vec![25745, 7313, 5009, 4721, 4685, 4681]
        );
    }

    #[test]
    fn overlapping_bins_never_has_duplicates() {
        for &(s, e) in &[
            (0u32, 1u32),
            (0, 1 << 29),
            (1 << 29, (1 << 29) + 1),
            ((1 << 29) - 1000, (1 << 29) + 1000),
            (0, MAX_END),
        ] {
            let bins = overlapping_bins(s, e);
            let mut sorted = bins.clone();
            sorted.sort_unstable();
            sorted.dedup();
            assert_eq!(sorted.len(), bins.len(), "duplicate bin for [{}, {})", s, e);
        }
    }

    #[test]
    fn self_consistency_assign_is_in_overlap() {
        // A feature is always found by a query over its own extent.
        for &(s, e) in &[
            (0u32, 1u32),
            (100, 200),
            (0, 1 << 17),
            ((1 << 29) - 1, 1 << 29),
            (1 << 29, (1 << 29) + 1),
            ((1 << 29) - 10, (1 << 29) + (1 << 20)), // straddles 2^29
            (MAX_END - 1, MAX_END),
        ] {
            assert!(
                overlapping_bins(s, e).contains(&assign_bin(s, e)),
                "assign_bin not in overlapping_bins for [{}, {})",
                s,
                e
            );
        }
    }

    #[test]
    #[should_panic(expected = "non-empty interval")]
    fn empty_interval_panics() {
        assign_bin(100, 100);
    }

    #[test]
    #[should_panic(expected = "non-empty interval")]
    fn zero_zero_panics() {
        // Guards against a `end - 1` underflow at (0, 0).
        overlapping_bins(0, 0);
    }

    #[test]
    #[should_panic(expected = "2 Gb")]
    fn out_of_range_panics() {
        assign_bin(0, MAX_END + 1);
    }

    proptest::proptest! {
        // The fundamental invariant: whenever a feature and a query overlap, the
        // feature's containing bin is among the query's overlapping bins. Ranges
        // are built around a shared pivot (concentrated near and above 2^29) so
        // they actually overlap and exercise the standard/extended boundary.
        #[test]
        fn overlap_implies_bin_in_query(
            pivot in proptest::sample::select(vec![
                1u32, 1 << 16, 1 << 17, (1 << 20) - 1, 1 << 23, 1 << 26,
                (1 << 29) - 2, 1 << 29, (1 << 29) + 2, 1 << 30, MAX_END - 2,
            ]),
            fa in 0u32..50_000,
            fb in 0u32..50_000,
            qa in 0u32..50_000,
            qb in 0u32..50_000,
        ) {
            // Force both half-open intervals to strictly contain `pivot`
            // (start <= pivot < end), so they are guaranteed to overlap there.
            // `pivot + 1` cannot overflow: the largest pivot is MAX_END - 2.
            let clamp = |x: u32| x.min(MAX_END - 1);
            let fs = clamp(pivot.saturating_sub(fa));
            let fe = clamp(pivot.saturating_add(fb)).max(pivot + 1).min(MAX_END);
            let qs = clamp(pivot.saturating_sub(qa));
            let qe = clamp(pivot.saturating_add(qb)).max(pivot + 1).min(MAX_END);
            let feature_bin = assign_bin(fs, fe);
            let query_bins = overlapping_bins(qs, qe);
            proptest::prop_assert!(
                query_bins.contains(&feature_bin),
                "feature [{},{}) bin {} not found by query [{},{})",
                fs,
                fe,
                feature_bin,
                qs,
                qe
            );
        }
    }
}
