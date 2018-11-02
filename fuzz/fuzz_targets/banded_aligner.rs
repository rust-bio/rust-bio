#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate bio;
use std::cmp::min;
use bio::alignment::pairwise::{banded, Scoring};
fuzz_target!(|data: &[u8]| {
    if data.len() < 50 || data.len() > 300 {
        return;
    }

    let (split_byte, data) = data.split_first().unwrap();
    let (kmer_byte, data) = data.split_first().unwrap();
    let (window_byte, data) = data.split_first().unwrap();

    let alphabets = b"ACGT";
    let v: Vec<_> = data.iter().map(|i| alphabets[(*i as usize)%alphabets.len()]).collect();

    let kmer_len: usize = 5 + (*kmer_byte as usize) % 10;
    let window_size: usize = 5 + (*window_byte as usize) % 10;
    let split_pos: usize = min(data.len()-1, *split_byte as usize);    

    let (x, y) = v.split_at(split_pos);
    println!("x: {}, y: {}, k: {}, w: {}", String::from_utf8(x.to_vec()).unwrap(), String::from_utf8(y.to_vec()).unwrap(), kmer_len, window_size);

    {
        let base_score = Scoring::from_scores(0, -1, 1, -1);
        let scoring = Scoring {
            xclip_prefix: 0,
            xclip_suffix: 0,
            yclip_suffix: 0,
            ..base_score
        };
        let mut aligner = banded::Aligner::with_scoring(scoring, kmer_len, window_size);
        let alignment = aligner.custom(x, y);
        assert_eq!(alignment.ystart, 0);
    }

    {
        let base_score = Scoring::from_scores(0, -1, 1, -1);
        let scoring = Scoring {
            xclip_prefix: 0,
            xclip_suffix: 0,
            yclip_prefix: 0,
            ..base_score
        };
        let mut aligner = banded::Aligner::with_scoring(scoring, kmer_len, window_size);
        let alignment = aligner.custom(x, y);
        assert_eq!(alignment.yend, alignment.ylen);
    }

    {
        let base_score = Scoring::from_scores(0, -1, 1, -1);
        let scoring = Scoring {
            xclip_suffix: 0,
            yclip_prefix: 0,
            yclip_suffix: 0,
            ..base_score
        };
        let mut aligner = banded::Aligner::with_scoring(scoring, kmer_len, window_size);
        let alignment = aligner.custom(x, y);
        assert_eq!(alignment.xstart, 0);
    }

    {
        let base_score = Scoring::from_scores(0, -1, 1, -1);
        let scoring = Scoring {
            xclip_prefix: 0,
            yclip_prefix: 0,
            yclip_suffix: 0,
            ..base_score
        };
        let mut aligner = banded::Aligner::with_scoring(scoring, kmer_len, window_size);
        let alignment = aligner.custom(x, y);
        assert_eq!(alignment.xend, alignment.xlen);
    }

});
