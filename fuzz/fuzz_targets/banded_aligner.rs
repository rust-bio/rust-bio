#![no_main]
#[macro_use]
extern crate libfuzzer_sys;
extern crate bio;
use std::cmp::{min, max};
use bio::alignment::pairwise::{self, banded, Scoring};
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

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring);
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.ystart, 0);

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

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring);
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.yend, f_alignment.ylen);

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

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring);
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.xstart, 0);

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

        // Full
        let mut f_aligner = pairwise::Aligner::with_scoring(scoring);
        let f_alignment = f_aligner.custom(x, y);
        assert_eq!(f_alignment.xend, f_alignment.xlen);

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

        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);

        let alignment = aligner.global(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, alignment.ylen);
    }

    {
        // full
        let mut aligner = pairwise::Aligner::with_scoring(base_score.clone());

        let alignment = aligner.local(x, y);
        assert!(alignment.score >= 0);

        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);

        let alignment = aligner.global(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, alignment.xlen);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, alignment.ylen);
    }

});
