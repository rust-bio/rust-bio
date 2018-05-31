// Copyright 2018 Jeff K
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Merge overlapping paired sequences
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::mating::{mate,merge};
//! let read1 = b"tacgattcgat";
//! let read2 = b"acgtaatcgaa";
//! let offset = match(&read1, &read2).unwrap();
//! let contig: TextSlice = merge(read1, read2, offset);
//! assert_eq!(merge(&read1, &read2, offset), b"tacgattcgattacgt");
//! ```

use std::cmp;

fn merge_records(r1: &Record, r2: &Record) -> Option<Record> {
    let r2_rc = dna::revcomp(r2.seq());

    match mate(&r1.seq(), &r2_rc) {
        Ok(overlap) => {
            let seq = merge(&r1.seq(), &r2_rc, overlap);

            // reverse r2 qual in place (this doesn't)
            let qual = merge(&r1.qual(), &r2.qual().reverse(), overlap);
            Record::with_attrs(r1.id(), None, &seq, &qual) 
        },
        None => None,
    }
}

fn mate(r1: &[u8], r2: [&u8]) -> Option<usize> {
    let min_offset = 25;
    let min_score = 0;
    let max_offset = cmp::min(r1.len(), r2.len());
    let mut m: i16 = 0;
    let mut pos: usize = 0;
    for i in min_offset..max_offset {
        let h = score(&r2[0..i], &r1[max_offset-i..max_offset]);
        if h > m {
            m = h;
            pos = i;
        }
    }
    if pos > min_offset && m > min_score {
        Some(pos)
    }
    None
}

fn score(r1: &[u8], r2: &[u8]) -> i16 {
    let mut s: i16 = 0;
    let len = r1.len();
    for i in 0..len {
        if r1[i] == r2[i] {
            s = s + 1;
        } else {
            s = s - 1;
        }
    }
    s
}

fn merge(r1: &[u8], r2: &[u8], overlap: usize) -> Vec<u8> {
    let r1_end = r1.seq().len() - overlap; 
    let r2_end = r2.seq().len() - overlap;
    let len = r1_end + r2_end + overlap;

    let mut seq = vec![0; len];
    seq[0..r1_end].copy_from_slice(&r1[0..r1_end]);
    seq[r1_end..len].copy_from_slice(&r2[0..r2.len()]);
    seq
}

#[cfg(test)]
mod tests {
    use super::{mate, merge};

    #[test]
    fn test_merge_pair() {
        let r1 = b"tacgattcgat";
        let r2 = b"acgtaatcgaa";
//        let offset = mate(&r1, &r2).unwrap();
        let contig: Vec<u8> = merge(r1.to_slice(), r2.to_slice(), 3);
        assert_eq!(contig, b"tacgattcgattacgt");
    }
}
