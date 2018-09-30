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

    let (dist, data) = data.split_first().unwrap();
    let (pattern_len, data) = data.split_first().unwrap();
    let dist = min(64, *dist);
    let pattern_len = max(1, min(64, min(data.len(), *pattern_len as usize)));

    if data.iter().any(|&b| b < 65 || b > 122 || b > 90 && b < 97) {
        return;
    }

    let (pattern, text) = data.split_at(pattern_len);

    let _ = MyersBuilder::new().build(pattern);

    let mut myers = Myers64::new(pattern);

    for _ in  myers.find_all(text, dist) {}

    // The following code compares the distance from the myers algorithm with the
    // number of substitutions / InDels found in the alignment path. If the distances
    // are equal, then the traceback found a valid alignment path.
    let mut aln = new_alignment();
    let mut matches = myers.find_all_pos_remember(text, ::std::u8::MAX);
    while let Some((end, dist)) = matches.next_end() {
        matches.alignment_at(end, &mut aln);
        // simple counting of non-matching operations
        let nomatch = aln
            .operations
            .iter()
            .filter(|&&o| o != Match)
            .count();
        // Additionally compares actual pattern and text and panics if matches / substitutions
        // are unexpectedly (un)equal.
        let d = get_alignment_dist(
            aln.operations.as_slice(),
            pattern,
            &text[aln.ystart..aln.yend]
        );
        assert!(dist == d as u8 && dist == nomatch as u8);
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
