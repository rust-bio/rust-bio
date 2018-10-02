#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate bio;

use std::cmp::{min, max};
use bio::pattern_matching::myers::{MyersBuilder, Myers64, new_alignment};
use bio::alignment::AlignmentOperation;
use bio::alignment::AlignmentOperation::*;


fuzz_target!(|data: &[u8]| {
    if data.len() < 3 {
        return;
    }

    let (max_dist, data) = data.split_first().unwrap();
    let (pattern_len, data) = data.split_first().unwrap();
    let max_dist = min(64, *max_dist);
    let pattern_len = max(1, min(64, min(data.len(), *pattern_len as usize)));

    if data.iter().any(|&b| b < 65 || b > 122 || b > 90 && b < 97) {
        return;
    }

    let (pattern, text) = data.split_at(pattern_len);

    let _ = MyersBuilder::new().build(pattern);

    let mut myers = Myers64::new(pattern);

    let end_dist: Vec<_> = myers.find_all_end(text, max_dist).collect();

    // The following code compares the distance from the myers algorithm with the
    // number of substitutions / InDels found in the alignment path. If the distances
    // are equal, then the traceback found a valid alignment path.
    // Additionally, the actual pattern and text are inspected; matches / substitutions
    // that are unexpectedly (un)equal will cause a panic.
    let mut aln = new_alignment();
    {
        let mut matches = myers.find_all(text, max_dist);
        let mut end_dist_iter = end_dist.iter();
        while matches.next_alignment(&mut aln) {
            let d = get_alignment_dist(
                aln.operations.as_slice(),
                pattern,
                &text[aln.ystart..aln.yend]
            );
            assert_eq!(d, aln.score as usize);
            assert!(d as u8 <= max_dist);
            // compare to earlier search
            let (end, dist) = end_dist_iter.next().unwrap();
            assert_eq!(*end + 1, aln.yend);
            assert_eq!(*dist, aln.score as u8);
        }
    }

    {
        // Lazy API
        let mut matches = myers.find_all_lazy(text, max_dist);
        let mut end_dist_iter = end_dist.iter();
        while let Some((end, dist)) = matches.next() {
            assert!(matches.alignment_at(end, &mut aln));
            let d = get_alignment_dist(
                aln.operations.as_slice(),
                pattern,
                &text[aln.ystart..aln.yend]
            );
            assert_eq!(d as u8, dist);
            assert_eq!(d, aln.score as usize);
            assert!(d as u8 <= max_dist);
            // compare to earlier search
            let (end, dist) = end_dist_iter.next().unwrap();
            assert_eq!(*end + 1, aln.yend);
            assert_eq!(*dist, aln.score as u8);
            // larger positions were not yet searched
            assert!(!matches.alignment_at(end + 1, &mut aln));
        }
    }

    // Lazy API with unlimited distance: each position should be found
    let mut matches = myers.find_all_lazy(text, ::std::u8::MAX);
    let mut i = 0;
    while let Some((end, dist)) = matches.next() {
        assert_eq!(end, i);
        assert!(matches.alignment_at(end, &mut aln));
        let d = get_alignment_dist(
            aln.operations.as_slice(),
            pattern,
            &text[aln.ystart..aln.yend]
        );
        assert_eq!(d as u8, dist);
        assert_eq!(d, aln.score as usize);
        // larger positions were not yet searched
        assert!(!matches.alignment_at(end + 1, &mut aln));
        i += 1;
    }
});


// Returns the edit distance between x & y given the alignment operations, but also
// checks if matches and substitutions are really correct given the actual sequences of
// x and y
fn get_alignment_dist(operations: &[AlignmentOperation], x: &[u8], y: &[u8]) -> usize {

    let mut dist = 0;

    let mut ix = 0;
    let mut iy = 0;
    for op in operations {
        match *op {
            Match => {
                assert!(x[ix] == y[iy], "Match operation, but characters are not equal");
                ix += 1;
                iy += 1;
            }
            Subst => {
                assert!(x[ix] != y[iy], "Subst operation, but characters are equal");
                dist += 1;
                ix += 1;
                iy += 1;
            }
            Del => {
                dist += 1;
                iy += 1;
            }
            Ins => {
                dist += 1;
                ix += 1;
            }
            _ => unreachable!()
        }
    }
    dist
}
