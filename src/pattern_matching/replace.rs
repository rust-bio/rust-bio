// Copyright 2020 Christopher Sugai.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Utilities to replace IUPAC and unknown characters to pseudorandom nucleotides. 
//!
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::replace;
//! ut_replace("ACTGuyUNn.-".to_string());
//! // ACTGtyTNn.-
//! rand_replace_gap_lowercase("ACTGuyUNn.-".to_string());
//! // ACTGuyUNnag
//! rand_replace_n("ACTGuyUNn.-".to_string());
//! // ACTGuyUAg.-
//! rand_replace_other("ACTGuyUNn.-".to_string());
//! // ACTGucUNn.-
//! //should only be used after other functions if unknown symbols are causing errors:
//! rand_replace("ACTGuyUNn.-@<^>".to_string()); 
//! // ACTGuTUGACTCTGA
//! ```

use rand::seq::SliceRandom;

// convert U to T
pub fn ut_replace(nucl: String) -> String {
    let dna: String = nucl.chars()
    .map(|ch| match ch {
            'U' => 'T',
            'u' => 't',
            _ => ch
        }).collect();
    dna
}


// fill gaps {.,-} with pseudorandom nucleotides ACTG
pub fn rand_replace_gap(nucl: String) -> String {
    let bases = ['A', 'C', 'T', 'G'];
    let mut rng = rand::thread_rng();
    let dna: String = nucl.chars()
    .map(|ch| match ch {
        '.' => *bases.choose(&mut rng).unwrap(),
        '-' => *bases.choose(&mut rng).unwrap(),
        _ => ch
        }).collect();
    dna
}

// fill gaps {.,-} with pseudorandom nucleotides actg
pub fn rand_replace_gap_lowercase(nucl: String) -> String {
    let bases = ['a', 'c', 't', 'g'];
    let mut rng = rand::thread_rng();
    let dna: String = nucl.chars()
    .map(|ch| match ch {
        '.' => *bases.choose(&mut rng).unwrap(),
        '-' => *bases.choose(&mut rng).unwrap(),
        _ => ch
        }).collect();
    dna
}

// fill N with pseudorandom nucleotides ACTG and n with actg
pub fn rand_replace_n(nucl: String) -> String {
    let bases = ['A', 'C', 'T', 'G'];
    let lowercase_bases = ['a', 'c', 't', 'g'];
    let mut rng = rand::thread_rng();
    let dna: String = nucl.chars()
    .map(|ch| match ch {
        'N' => *bases.choose(&mut rng).unwrap(),
        'n' => *lowercase_bases.choose(&mut rng).unwrap(),
        _ => ch
        }).collect();
    dna
}

// pseudorandom nucleotide replacements within IUPAC specifications, e.g. R: either A or G. Does not fill N's or gaps. Refer to other functions for these.
pub fn rand_replace_other(nucl: String) -> String {
    let mut rng = rand::thread_rng();
    let r_bases = ['A', 'G'];
    let r_bases_lowercase = ['a', 'g'];
    let y_bases = ['C', 'T'];
    let y_bases_lowercase = ['c', 't'];
    let s_bases = ['C', 'G'];
    let s_bases_lowercase = ['c', 'g'];
    let w_bases = ['A', 'T'];
    let w_bases_lowercase = ['a', 't'];
    let k_bases = ['T', 'G'];
    let k_bases_lowercase = ['t', 'g'];
    let m_bases = ['A', 'C'];
    let m_bases_lowercase = ['a', 'c'];
    let b_bases = ['C', 'T', 'G'];
    let b_bases_lowercase = ['c', 't', 'g'];
    let d_bases = ['A', 'T', 'G'];
    let d_bases_lowercase = ['a', 't', 'g'];
    let h_bases = ['A', 'C', 'T'];
    let h_bases_lowercase = ['a', 'c', 't'];
    let v_bases = ['A', 'C', 'G'];
    let v_bases_lowercase = ['a', 'c', 'g'];
    let dna: String = nucl.chars()
    .map(|ch| match ch {
        'R' => *r_bases.choose(&mut rng).unwrap(),
        'r' => *r_bases_lowercase.choose(&mut rng).unwrap(),
        'Y' => *y_bases.choose(&mut rng).unwrap(),
        'y' => *y_bases_lowercase.choose(&mut rng).unwrap(),
        'S' => *s_bases.choose(&mut rng).unwrap(),
        's' => *s_bases_lowercase.choose(&mut rng).unwrap(),
        'W' => *w_bases.choose(&mut rng).unwrap(),
        'w' => *w_bases_lowercase.choose(&mut rng).unwrap(),
        'K' => *k_bases.choose(&mut rng).unwrap(),
        'k' => *k_bases_lowercase.choose(&mut rng).unwrap(),
        'M' => *m_bases.choose(&mut rng).unwrap(),
        'm' => *m_bases_lowercase.choose(&mut rng).unwrap(),
        'B' => *b_bases.choose(&mut rng).unwrap(),
        'b' => *b_bases_lowercase.choose(&mut rng).unwrap(),
        'D' => *d_bases.choose(&mut rng).unwrap(),
        'd' => *d_bases_lowercase.choose(&mut rng).unwrap(),
        'H' => *h_bases.choose(&mut rng).unwrap(),
        'h' => *h_bases_lowercase.choose(&mut rng).unwrap(),
        'V' => *v_bases.choose(&mut rng).unwrap(),
        'v' => *v_bases_lowercase.choose(&mut rng).unwrap(),
        _ => ch
        }).collect();
    dna
}

// fill all other than ACGT with pseudorandom nucleotides ACTG. Should be used after other functions for cleanup of unknown characters.
pub fn rand_replace(nucl: String) -> String {
    let bases = ['A', 'C', 'T', 'G'];
    let mut rng = rand::thread_rng();
    let dna: String = nucl.chars()
    .map(|ch| match ch {
        'A' => ch,
        'a' => ch,
        'C' => ch,
        'c' => ch,
        'T' => ch,
        't' => ch,
        'G' => ch,
        'g' => ch,
        'U' => ch,
        'u' => ch,
        _ => *bases.choose(&mut rng).unwrap()
        }).collect();
    dna
}

// fill all other than ACGT with pseudorandom nucleotides actg. Should be used after other functions for cleanup of unknown characters.
pub fn rand_replace_lowercase(nucl: String) -> String {
    let bases = ['a', 'c', 't', 'g'];
    let mut rng = rand::thread_rng();
    let dna: String = nucl.chars()
    .map(|ch| match ch {
        'A' => ch,
        'a' => ch,
        'C' => ch,
        'c' => ch,
        'T' => ch,
        't' => ch,
        'G' => ch,
        'g' => ch,
        'U' => ch,
        'u' => ch,
        _ => *bases.choose(&mut rng).unwrap()
        }).collect();
    dna
}

#[test]
fn test_rand_replace() {
    let test = rand_replace("ACTGuyUNn.-@<^>".to_string());
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".len());
}