#![no_main]
#[macro_use]
extern crate libfuzzer_sys;

use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation::*;
use bio::pattern_matching::myers::{long, Myers, MyersBuilder};

const HEADER_LEN: usize = 11;

fuzz_target!(|data: &[u8]| {
    if data.len() < HEADER_LEN + 1 {
        return;
    }
    let header = &data[..HEADER_LEN];
    let data = &data[HEADER_LEN..];
    if data.iter().any(|b| !b.is_ascii_graphic()) {
        return;
    }

    let mut builder = MyersBuilder::new();

    // 1: max. distance
    let max_dist = header[0] / 32;
    // 2: pattern length (limited to 128 symbols)
    let pattern_len = header[1].rotate_left(max_dist as u32) as usize;
    let pattern_len = pattern_len.max(1).min(128).min(data.len());
    let (pattern, text) = data.split_at(pattern_len);
    let mut text = text.to_vec();

    // 3: wildcard yes/no
    // 4-5: range of wildcard positions in the text
    let mut wildcard = None;
    let w = header[2];
    if w > b'z' && w.is_ascii() {
        let (i0, i1) = (header[3] as usize, header[4] as usize);
        if i1 >= i0 && i1 < text.len() {
            for i in i0..=i1 {
                text[i] = w;
            }
            builder.text_wildcard(w);
            wildcard = Some(w);
        }
    }

    // 6: ambiguous symbol
    // 7-11: equivalent symbols
    let mut ambigs = None;
    let ambig = header[5];
    if ambig.is_ascii_graphic() && ambig > b'z' {
        let equivalents: Vec<u8> = header[6..]
            .iter()
            .cloned()
            .filter(|&b| b.is_ascii_graphic())
            .collect();
        if !equivalents.is_empty() {
            builder.ambig(ambig, &equivalents);
            ambigs = Some((ambig, equivalents));
        }
    }
    let ambigs = ambigs.as_ref().map(|(a, eq)| (*a, eq.as_slice()));

    let text = &text;
    // dbg!(std::str::from_utf8(pattern), std::str::from_utf8(text), max_dist, wildcard.map(|w| w as char), ambigs);

    macro_rules! run_fuzz {
        ($myers:ident, $myers_long:ident) => {
            // No traceback, just searching
            let end_dist: Vec<_> = $myers.find_all_end(text, max_dist).collect();
            let end_dist_long: Vec<_> = $myers_long
                .find_all_end(text, max_dist as usize)
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
                let mut matches = $myers.find_all(text, max_dist);
                let mut matches_long = $myers_long.find_all(text, max_dist as usize);

                matches.alignment(&mut aln); // all insertions to text
                matches_long.alignment(&mut aln_long);

                let mut end_dist_iter = end_dist.iter();
                while matches.next_alignment(&mut aln) {
                    assert!(matches_long.next_alignment(&mut aln_long));
                    assert_eq!(aln, aln_long);

                    // verify alignment
                    validate_alignment(&aln, pattern, text, wildcard, ambigs);
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

                let mut matches = $myers.find_all_lazy(text, max_dist);
                let mut matches_long = $myers_long.find_all_lazy(text, max_dist as usize);
                let mut end_dist_iter = end_dist.iter();

                while let Some((end, dist)) = matches.next() {
                    // compare distances
                    let (end_long, dist_long) = matches_long.next().unwrap();
                    assert_eq!((end, dist), (end_long, dist_long as u8));
                    assert_eq!(matches.dist_at(end), Some(dist));
                    assert_eq!(matches_long.dist_at(end), Some(dist_long));

                    // compare alignments
                    assert!(matches.alignment_at(end, &mut aln));
                    assert!(matches_long.alignment_at(end, &mut aln_long));
                    assert_eq!(aln, aln_long);

                    // verify alignment
                    validate_alignment(&aln, pattern, text, wildcard, ambigs);
                    assert_eq!(dist, aln.score as u8);
                    assert!(aln.score as u8 <= max_dist);

                    // compare to find_all_end() results
                    let (end, dist) = end_dist_iter.next().unwrap();
                    assert_eq!(*end + 1, aln.yend);
                    assert_eq!(*dist, aln.score as u8);

                    // larger positions were not yet searched
                    assert!(!matches.alignment_at(end + 1, &mut aln));
                    assert!(matches.dist_at(end + 1).is_none());
                    assert!(matches_long.dist_at(end + 1).is_none());
                }
                assert!(end_dist_iter.next().is_none());
            }

            // Lazy API with unlimited distance: all positions should be found

            let mut matches = $myers.find_all_lazy(text, u8::max_value());
            let mut matches_long = $myers_long.find_all_lazy(text, u8::max_value() as usize);

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
                validate_alignment(&aln, pattern, text, wildcard, ambigs);
                assert_eq!(dist, aln.score as u8);
                assert!(aln.score as u8 <= u8::max_value());

                // larger positions were not yet searched
                assert!(!matches.alignment_at(end + 1, &mut aln));
                assert!(!matches_long.alignment_at(end + 1, &mut aln_long));
                i += 1;
            }
        };
    }

    let mut myers_long: long::Myers<u8> = builder.build_long(pattern);
    if pattern.len() <= 64 {
        let mut myers: Myers<u64> = builder.build(pattern);
        run_fuzz!(myers, myers_long);
    } else {
        let mut myers: Myers<u128> = builder.build(pattern);
        run_fuzz!(myers, myers_long);
    }

    // compare to edlib (if no wildcard/ambiguities)
    // if pattern.len() <= 64 && wildcard.is_none() && ambigs.is_none() {
    //     use edlib_rs::edlibrs::*;
    //     let mut config = EdlibAlignConfigRs::default();
    //     config.k = max_dist as i32;
    //     config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    //     config.task = EdlibAlignTaskRs::EDLIB_TASK_PATH;

    //     let mut myers: Myers<u64> = builder.build(pattern);
    //     let mut matches = myers.find_all_lazy(text, max_dist);
    //     let mut aln = Alignment::default();

    //     let hits: Vec<_> = matches.by_ref().collect();
    //     let best_hit = hits.iter().min_by_key(|(_, dist)| dist);
    //     if let Some((best_end, best_dist)) = best_hit {
    //         let res = edlibAlignRs(pattern, text, &config);
    //         assert_eq!(res.status, EDLIB_STATUS_OK);
    //         assert_eq!(res.editDistance, *best_dist as i32);
    //     // TODO: due to differences in traceback implementation,
    //     //       so we cannot compare the alignment paths directly 
    // }
});

// Validates an Alignment based on the sequences that were used to construct it.
// - calculates the score using edit distance (mismatches / gap penalties = 1) and
//   then compares it to the stored score
// - checks if matches and substitutions are really correct given the actual sequences
fn validate_alignment(
    aln: &Alignment,
    x: &[u8],
    y: &[u8],
    wildcard: Option<u8>,
    ambigs: Option<(u8, &[u8])>,
) {
    // dbg!(aln, std::str::from_utf8(x), std::str::from_utf8(y), wildcard, ambigs);
    let y = &y[aln.ystart..aln.yend];
    let mut ix = 0;
    let mut iy = 0;
    let mut calc_dist = 0;
    for op in &aln.operations {
        match *op {
            Match | Subst => {
                let matches = x[ix] == y[iy]
                    || wildcard.map(|w| w == y[iy]).unwrap_or(false)
                    || ambigs
                        .map(|(a, equiv)| x[ix] == a && equiv.contains(&y[iy]))
                        .unwrap_or(false);
                assert!(
                    (*op == Subst) ^ matches,
                    "{:?} operation invalid: {} {} {}",
                    op,
                    x[ix] as char,
                    if *op == Match { "!=" } else { "==" },
                    y[iy] as char
                );
                calc_dist += !matches as usize;
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
            _ => unreachable!(),
        }
    }

    assert_eq!(calc_dist, aln.score as usize, "Score mismatch");
}
