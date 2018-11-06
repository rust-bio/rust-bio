#![no_main]
#[macro_use]
extern crate libfuzzer_sys;
extern crate bio;
use std::cmp::{min, max};
use bio::alignment::{Alignment, AlignmentMode, AlignmentOperation};
use bio::alignment::pairwise::{self, banded, Scoring, MatchParams,  MatchFunc};
use bio::utils::TextSlice;

fn validate_alignment_score(al: &Alignment, x: TextSlice, y: TextSlice, scoring: &Scoring<MatchParams>) -> bool {
    use AlignmentOperation::*;
    let path = al.path();
    let mut score = 0;
    if al.mode==AlignmentMode::Custom {
        if al.xstart > 0 {
            score += scoring.xclip_prefix;
        }
        if al.ystart > 0 {
            score += scoring.yclip_prefix;
        }
        if al.xend < al.xlen {
            score += scoring.xclip_suffix;
        }
        if al.yend < al.ylen {
            score += scoring.yclip_suffix;
        }
    }
    let mut last_op = None;
    for (i, j, op) in path {
        score += match op {
            Match | Subst => scoring.match_fn.score(x[i-1], y[j-1]),
            Del => if last_op==Some(Del) { scoring.gap_extend } else { scoring.gap_open + scoring.gap_extend },
            Ins => if last_op==Some(Ins) { scoring.gap_extend } else { scoring.gap_open + scoring.gap_extend },
            _ => 0,
        };
        last_op = Some(op);
    }
    al.score==score
}

fuzz_target!(|data: &[u8]| {
    if data.len() < 50 || data.len() > 300 {
        return;
    }

    let (split_byte, data) = data.split_first().unwrap();
    let (kmer_byte, data) = data.split_first().unwrap();
    let (window_byte, data) = data.split_first().unwrap();
    let (match_score_byte, data) = data.split_first().unwrap();
    let (mismatch_score_byte, data) = data.split_first().unwrap();
    let (gap_open_byte, data) = data.split_first().unwrap();
    let (gap_extend_byte, data) = data.split_first().unwrap();
    let (xclip_prefix_byte, data) = data.split_first().unwrap();
    let (xclip_suffix_byte, data) = data.split_first().unwrap();
    let (yclip_prefix_byte, data) = data.split_first().unwrap();
    let (yclip_suffix_byte, data) = data.split_first().unwrap();

    let alphabets = b"ACGT";
    let v: Vec<_> = data.iter()
        .map(|i| alphabets[(*i as usize) % alphabets.len()])
        .collect();

    let kmer_len: usize = 5 + (*kmer_byte as usize) % 10;
    let window_size: usize = 5 + (*window_byte as usize) % 10;
    let split_pos: usize = min(data.len() - 1, max(*split_byte as usize, 1));

    let match_score = 1 + (*match_score_byte as i32) % 5;
    let mismatch_score = -((*mismatch_score_byte as i32) % 10);
    let gap_open = -((*gap_open_byte as i32) % 20);
    let gap_extend = -((*gap_extend_byte as i32) % 10);

    let (x, y) = v.split_at(split_pos);
    println!("x: {}, y: {}, k: {}, w: {}, scoring ({}, {}, {}, {})",
             String::from_utf8(x.to_vec()).unwrap(),
             String::from_utf8(y.to_vec()).unwrap(),
             kmer_len,
             window_size,
             gap_open,
             gap_extend,
             match_score,
             mismatch_score);
    let base_score = Scoring::from_scores(gap_open, gap_extend, match_score, mismatch_score);

    {
        println!("Clip scores ({}, {}, {}, {})", xclip_prefix_byte, xclip_suffix_byte, yclip_prefix_byte, yclip_suffix_byte);
        let scoring = Scoring {
            xclip_prefix: -(*xclip_prefix_byte as i32),
            xclip_suffix: -(*xclip_suffix_byte as i32),
            yclip_prefix: -(*yclip_prefix_byte as i32),
            yclip_suffix: -(*yclip_suffix_byte as i32),
            ..base_score.clone()
        };
        // Banded
        let mut b_aligner = banded::Aligner::with_scoring(scoring.clone(), kmer_len, window_size);
        let b_alignment = b_aligner.custom(x, y);
        assert!(validate_alignment_score(&b_alignment, x, y, &scoring));

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring.clone());
        let f_alignment = f_aligner.custom(x, y);
        assert!(validate_alignment_score(&f_alignment, x, y, &scoring));

        // Compare
        // Passing an empty match will force the banded alignmer to band the full matrix
        let band_all_alignment = b_aligner.custom_with_matches(x, y, &Vec::new());
        assert_eq!(band_all_alignment.score, f_alignment.score);
    }


    {
        let scoring = Scoring {
            xclip_prefix: 0,
            xclip_suffix: 0,
            yclip_suffix: 0,
            ..base_score.clone()
        };
        // Banded
        let mut b_aligner = banded::Aligner::with_scoring(scoring.clone(), kmer_len, window_size);
        let b_alignment = b_aligner.custom(x, y);
        assert_eq!(b_alignment.ystart, 0);
        assert!(validate_alignment_score(&b_alignment, x, y, &scoring));

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring.clone());
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.ystart, 0);
        assert!(validate_alignment_score(&f_alignment, x, y, &scoring));

        // Compare
        // Passing an empty match will force the banded alignmer to band the full matrix
        let band_all_alignment = b_aligner.custom_with_matches(x, y, &Vec::new());
        assert_eq!(band_all_alignment.score, f_alignment.score);
    }

    {
        let scoring = Scoring {
            xclip_prefix: 0,
            xclip_suffix: 0,
            yclip_prefix: 0,
            ..base_score.clone()
        };
        // Banded
        let mut b_aligner = banded::Aligner::with_scoring(scoring.clone(), kmer_len, window_size);
        let b_alignment = b_aligner.custom(x, y);
        assert_eq!(b_alignment.yend, b_alignment.ylen);
        assert!(validate_alignment_score(&b_alignment, x, y, &scoring));

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring.clone());
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.yend, f_alignment.ylen);
        assert!(validate_alignment_score(&f_alignment, x, y, &scoring));

        // Compare
        // Passing an empty match will force the banded alignmer to band the full matrix
        let band_all_alignment = b_aligner.custom_with_matches(x, y, &Vec::new());
        assert_eq!(band_all_alignment.score, f_alignment.score);
    }

    {
        let scoring = Scoring {
            xclip_suffix: 0,
            yclip_prefix: 0,
            yclip_suffix: 0,
            ..base_score.clone()
        };
        // Banded
        let mut b_aligner = banded::Aligner::with_scoring(scoring.clone(), kmer_len, window_size);
        let b_alignment = b_aligner.custom(x, y);
        assert_eq!(b_alignment.xstart, 0);
        assert!(validate_alignment_score(&b_alignment, x, y, &scoring));

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring.clone());
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.xstart, 0);
        assert!(validate_alignment_score(&f_alignment, x, y, &scoring));

        // Compare
        // Passing an empty match will force the banded alignmer to band the full matrix
        let band_all_alignment = b_aligner.custom_with_matches(x, y, &Vec::new());
        assert_eq!(band_all_alignment.score, f_alignment.score);
    }

    {
        let scoring = Scoring {
            xclip_prefix: 0,
            yclip_prefix: 0,
            yclip_suffix: 0,
            ..base_score.clone()
        };
        // Banded
        let mut b_aligner = banded::Aligner::with_scoring(scoring.clone(), kmer_len, window_size);
        let b_alignment = b_aligner.custom(x, y);
        assert_eq!(b_alignment.xend, b_alignment.xlen);
        assert!(validate_alignment_score(&b_alignment, x, y, &scoring));

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring.clone());
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.xend, f_alignment.xlen);
        assert!(validate_alignment_score(&f_alignment, x, y, &scoring));

        // Compare
        // Passing an empty match will force the banded alignmer to band the full matrix
        let band_all_alignment = b_aligner.custom_with_matches(x, y, &Vec::new());
        assert_eq!(band_all_alignment.score, f_alignment.score);
    }

    {
        // banded
        let mut aligner = banded::Aligner::with_scoring(base_score.clone(), kmer_len, window_size);

        let alignment = aligner.local(x, y);
        assert!(alignment.score >= 0);
        assert!(validate_alignment_score(&alignment, x, y, &base_score));

        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);
        assert!(validate_alignment_score(&alignment, x, y, &base_score));

        let alignment = aligner.global(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, alignment.ylen);
        assert!(validate_alignment_score(&alignment, x, y, &base_score));
    }

    {
        // full
        let mut aligner = pairwise::Aligner::with_scoring(base_score.clone());

        let alignment = aligner.local(x, y);
        assert!(alignment.score >= 0);
        assert!(validate_alignment_score(&alignment, x, y, &base_score));

        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);
        assert!(validate_alignment_score(&alignment, x, y, &base_score));

        let alignment = aligner.global(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, alignment.ylen);
        assert!(validate_alignment_score(&alignment, x, y, &base_score));
    }

});
