#![no_main]
#[macro_use] extern crate libfuzzer_sys;

use std::cmp::{min, max};
use bio::pattern_matching::myers::{MyersBuilder, Myers, long};
use bio::alignment::Alignment;
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

    // Test whether builders succeed
    // TODO: Builder testing could be expanded
    let max_dist = max_dist as u8;
    // Test whether builders succeed
    // TODO: Builder testing could be expanded
    let _ = MyersBuilder::new().build_64(pattern);
    let _ = MyersBuilder::new().build_long_64(pattern);

    // Myers objects
    let mut myers = Myers::<u64>::new(pattern);
    let mut myers_long = long::Myers::<u8>::new(pattern);

    // No traceback, just searching
    let end_dist: Vec<_> = myers.find_all_end(text, max_dist).collect();
    let end_dist_long: Vec<_> = myers_long.find_all_end(text, max_dist as usize)
        .map(|(end, dist)| (end, dist as u8))
        .collect();
    assert_eq!(end_dist, end_dist_long);

    // Test traceback algorithm:
    // The following code compares the distance from the myers algorithm with the
    // number of substitutions / InDels found in the alignment path. If the distances
    // are equal, then the traceback found a valid alignment path.
    // Additionally, the actual pattern and text are inspected; matches / substitutions
    // that are unexpectedly (un)equal will cause a panic.

    // 'Default' API
    let mut aln = Alignment::default();
    let mut aln_long = Alignment::default();
    {
        let mut matches = myers.find_all(text, max_dist);
        let mut matches_long = myers_long.find_all(text, max_dist as usize);

        matches.alignment(&mut aln); // all insertions to text
        matches_long.alignment(&mut aln_long);

        let mut end_dist_iter = end_dist.iter();
        while matches.next_alignment(&mut aln) {
            assert!(matches_long.next_alignment(&mut aln_long));
            assert_eq!(aln, aln_long);

            // verify alignment
            validate_alignment(&aln, pattern, text);
            assert!(aln.score as u8 <= max_dist);

            // compare to find_all_end() results
            let (end, dist) = end_dist_iter.next().unwrap();
            assert_eq!(*end + 1, aln.yend);
            assert_eq!(*dist, aln.score as u8);
        }
        assert!(end_dist_iter.next().is_none());
        assert!(!matches_long.next_alignment(&mut aln_long));
    }

    {
        // Lazy API

        let mut matches = myers.find_all_lazy(text, max_dist);
        let mut matches_long = myers_long.find_all_lazy(text, max_dist as usize);
        let mut end_dist_iter = end_dist.iter();

        while let Some((end, dist)) = matches.next() {
            // compare distances
            let (end_long, dist_long) = matches_long.next().unwrap();
            assert_eq!((end, dist), (end_long, dist_long as u8));

            // compare alignments
            assert!(matches.alignment_at(end, &mut aln));
            assert!(matches_long.alignment_at(end, &mut aln_long));
            assert_eq!(aln, aln_long);

            // verify alignment
            validate_alignment(&aln, pattern, text);
            assert_eq!(dist, aln.score as u8);
            assert!(aln.score as u8 <= max_dist);

            // compare to find_all_end() results
            let (end, dist) = end_dist_iter.next().unwrap();
            assert_eq!(*end + 1, aln.yend);
            assert_eq!(*dist, aln.score as u8);

            // larger positions were not yet searched
            assert!(!matches.alignment_at(end + 1, &mut aln));
        }
        assert!(end_dist_iter.next().is_none());
    }

    // Lazy API with unlimited distance: each position should be found

    let mut matches = myers.find_all_lazy(text, u8::max_value());
    let mut matches_long = myers_long.find_all_lazy(text, u8::max_value() as usize);

    let mut i = 0;
    while let Some((end, dist)) = matches.next() {
        assert_eq!(end, i);

        // compare distances
        let (end_long, dist_long) = matches_long.next().unwrap();
        assert_eq!((end, dist), (end_long, dist_long as u8));

        // compare alignments
        assert!(matches.alignment_at(end, &mut aln));
        assert!(matches_long.alignment_at(end, &mut aln_long));
        assert_eq!(aln, aln_long);

        // verify alignment
        validate_alignment(&aln, pattern, text);
        assert_eq!(dist, aln.score as u8);
        assert!(aln.score as u8 <= u8::max_value());

        // larger positions were not yet searched
        assert!(!matches.alignment_at(end + 1, &mut aln));
        assert!(!matches_long.alignment_at(end + 1, &mut aln_long));
        i += 1;
    }
});


// Validates an Alignment based on the sequences that were used to construct it.
// - calculates the score using edit distance (mismatches / gap penalties = 1) and
//   then compares it to the stored score
// - checks if matches and substitutions are really correct given the actual sequences
fn validate_alignment(aln: &Alignment, x: &[u8], y: &[u8]) {

    let y = &y[aln.ystart..aln.yend];

    let mut ix = 0;
    let mut iy = 0;
    let mut calc_dist = 0;
    for op in &aln.operations {
        match *op {
            Match => {
                assert!(x[ix] == y[iy], "Match operation, but characters are not equal");
                ix += 1;
                iy += 1;
            }
            Subst => {
                assert!(x[ix] != y[iy], "Subst operation, but characters are equal");
                calc_dist += 1;
                ix += 1;
                iy += 1;
            }
            Del => {
                calc_dist += 1;
                iy += 1;
            }
            Ins => {
                calc_dist += 1;
                ix += 1;
            }
            _ => unreachable!()
        }
    }

    assert_eq!(calc_dist, aln.score as usize);
}
