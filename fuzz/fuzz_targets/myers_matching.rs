#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate bio;

use std::cmp::{min, max};
use bio::pattern_matching::myers::{MyersBuilder, Myers64};


fuzz_target!(|data: &[u8]| {
    if data.len() < 3 {
        return;
    }

    let (dist, data) = data.split_first().unwrap();
    let (pattern_len, data) = data.split_first().unwrap();

    let dist = min(64, *dist as usize);
    let pattern_len = max(1, min(64, min(data.len(), *pattern_len as usize)));

    let (pattern, text) = data.split_at(pattern_len);

    let _ = MyersBuilder::new().build(pattern);

    let mut myers = Myers64::new(pattern);

    for _ in  myers.find_all_pos(text, dist) {}

    let mut matches = myers.find_all_pos_remember(text, dist);
    for _ in matches.by_ref() {}

    matches.hit_at(1);
});
