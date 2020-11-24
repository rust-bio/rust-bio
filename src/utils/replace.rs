// Copyright 2020 Christopher Sugai.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! 
//! Traits to replace non-ACTG nucleotide characters with pseudorandom bases. Additional functionality to create strings of pseudorandom nucleotides and convert to uppercase.
//! # Example
//! ```
//! 
//! extern crate bio;
//! extern crate rand;
//! 
//! use rand::rngs::ThreadRng;
//! use std::string::String;
//! use rand::seq::SliceRandom;
//! use bio::utils::TextSlice;
//! use bio::utils::Text;
//! 
//! let mut rng = rand::thread_rng(); //create a random number generator
//! let x1 = "ACTGuyUNn.-@<^>".to_string(); //string
//! let x2 = "ACTGuyUNn.-@<^>"; //stringslice
//! let x3: TextSlice = "ACTGuyUNn.-@<^>".as_bytes(); //rust bio textslice
//! let x4: Text = "ACTGuyUNn.-@<^>".as_bytes().to_vec(); //rust bio text
//!
//!
//! let y1 = x1.replace_gap(rng); //replace gap {.-} using random number generator
//!
//! let y2 = x2.to_string(); //convert stringslice to either Text or String but must be owned
//! let y2 = y2.replace_n(rng); //replace ambiguous/uncalled base {nN}
//!
//! let y3 = std::str::from_utf8(x3).expect("Not valid UTF8").to_string(); //convert textslice to either Text or String but must be owned
//! let y3 = y3.replace_u_with_t(); //replace all other text other than ACTGUactgu, does not take rng
//! 
//! let y4 = x4.replace_all(rng); //replace u with t
//! let y4 = std::string::String::from_utf8(y4).expect("Not valid UTF8");
//!
//! assert_eq!(y4.len(), "ACTGuyUNn.-@<^>".to_string().len());
//!
//! ```
//! 
//! 
//! 
//! 

extern crate bio;
use rand::rngs::ThreadRng;
use std::string::String;
use rand::seq::SliceRandom;

use bio::utils::TextSlice;
use bio::utils::Text;

trait ReplaceNucleotide <T: ToOwned> {
    fn replace_u_with_t(&self) -> Self;
    fn replace_gap(&self, rng: ThreadRng) -> Self;
    fn replace_gap_lowercase(&self, rng: ThreadRng) -> Self;
    fn replace_n(&self, rng: ThreadRng) -> Self;
    fn replace_iupac(&self, rng: ThreadRng) -> Self;
    fn replace_lowercase_with_uppercase(&self, rng: ThreadRng) -> Self;
    fn replace_all_other_with_uppercase(&self, rng: ThreadRng) -> Self;
    fn replace_all_other_with_lowercase(&self, rng: ThreadRng) -> Self;
}

impl ReplaceNucleotide<String> for String {
    /// Specifically replace Uu with Tt, lowercase to lowercase, uppercase to uppercase. Leaves all other characters in place.
    fn replace_u_with_t(&self) -> Self {
        self.chars()
            .map(|ch| match ch {
                'U' => 'T',
                'u' => 't',
                _ => ch,
            })
            .collect()
    }

    /// Specifically fill lowercase actgu with ACTGU but leaves all other characters in place.
    fn replace_lowercase_with_uppercase(&self) -> Self {
        self.chars()
            .map(|ch| match ch {
                'A' => ch,
                'a' => 'A',
                'C' => ch,
                'c' => 'C',
                'T' => ch,
                't' => 'T',
                'G' => ch,
                'g' => 'G',
                'U' => ch,
                'u' => 'U',
                _ => ch,
            }
    }    

    /// Fill gaps {.,-} with pseudorandom nucleotides ACTG
    fn replace_gap(&self, mut rng: ThreadRng) -> Self {
        let bases = ['A', 'C', 'T', 'G'];
        self.chars()
            .map(|ch| match ch {
                '.' => *bases.choose(&mut rng).unwrap(),
                '-' => *bases.choose(&mut rng).unwrap(),
                _ => ch,
            })
            .collect()
    }

    /// Fill gaps {.,-} with pseudorandom nucleotides actg
    fn replace_gap_lowercase(&self, mut rng: ThreadRng) -> Self {
        let bases = ['a', 'c', 't', 'g'];
        self.chars()
            .map(|ch| match ch {
                '.' => *bases.choose(&mut rng).unwrap(),
                '-' => *bases.choose(&mut rng).unwrap(),
                _ => ch,
            })
            .collect()
    }

    /// Fill N with pseudorandom nucleotides ACTG and n with actg
    fn replace_n(&self, mut rng: ThreadRng) -> Self {
        let bases = ['A', 'C', 'T', 'G'];
        let lowercase_bases = ['a', 'c', 't', 'g'];
        self.chars()
            .map(|ch| match ch {
                'N' => *bases.choose(&mut rng).unwrap(),
                'n' => *lowercase_bases.choose(&mut rng).unwrap(),
                _ => ch,
            })
            .collect()
    }
    /// Pseudorandom nucleotide replacements within IUPAC specifications, e.g. R: either A or G. Case specific, r: either a or g.
    fn replace_iupac(&self, mut rng: ThreadRng) -> Self {
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
        self.chars()
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
                _ => ch,
            })
            .collect()
    }

    /// fill all other than ACGTUactgu with pseudorandom nucleotides ACTGU. Should be used last after other functions or for cleanup of unknown characters.
    fn replace_all_other_with_uppercase(&self, mut rng: ThreadRng) -> Self {
        let bases = ['A', 'C', 'T', 'G'];
        self.chars()
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
                _ => *bases.choose(&mut rng).unwrap(),
            })
            .collect()
    }
    /// fill all other than ACGTUactgu with pseudorandom nucleotides actgu. Should be used last after other functions or for cleanup of unknown characters.
    fn replace_all_other_with_lowercase(&self, mut rng: ThreadRng) -> Self {
        let bases = ['a', 'c', 't', 'g'];
        self.chars()
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
                _ => *bases.choose(&mut rng).unwrap(),
            })
            .collect()
    }
}

impl ReplaceNucleotide<Text> for Text {
    /// Specifically replace Uu with Tt, lowercase to lowercase, uppercase to uppercase. Leaves all other characters in place.
    fn replace_u_with_t(&self) -> Self {
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
            .map(|ch| match ch {
                'U' => 'T',
                'u' => 't',
                _ => ch,
            })
            .collect::<String>().into_bytes()
    }

    /// Specifically fill lowercase actgu with ACTGU but leaves all other characters in place.
    fn replace_lowercase_with_uppercase(&self) -> Self {
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
            .map(|ch| match ch {
                'A' => ch,
                'a' => 'A',
                'C' => ch,
                'c' => 'C',
                'T' => ch,
                't' => 'T',
                'G' => ch,
                'g' => 'G',
                'U' => ch,
                'u' => 'U',
                _ => ch,
            }
    }

    /// Fill gaps {.,-} with pseudorandom nucleotides ACTG
    fn replace_gap(&self, mut rng: ThreadRng) -> Self {
        let bases = ['A', 'C', 'T', 'G'];
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
            .map(|ch| match ch {
                '.' => *bases.choose(&mut rng).unwrap(),
                '-' => *bases.choose(&mut rng).unwrap(),
                _ => ch,
            })
            .collect::<String>().into_bytes()
    }
    /// Fill gaps {.,-} with pseudorandom nucleotides actg
    fn replace_gap_lowercase(&self, mut rng: ThreadRng) -> Self {
        let bases = ['a', 'c', 't', 'g'];
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
            .map(|ch| match ch {
                '.' => *bases.choose(&mut rng).unwrap(),
                '-' => *bases.choose(&mut rng).unwrap(),
                _ => ch,
            })
            .collect::<String>().into_bytes()
    }
    /// Fill N with pseudorandom nucleotides ACTG and n with actg
    fn replace_n(&self, mut rng: ThreadRng) -> Self {
        let bases = ['A', 'C', 'T', 'G'];
        let lowercase_bases = ['a', 'c', 't', 'g'];
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
            .map(|ch| match ch {
                'N' => *bases.choose(&mut rng).unwrap(),
                'n' => *lowercase_bases.choose(&mut rng).unwrap(),
                _ => ch,
            })
            .collect::<String>().into_bytes()
    }
    /// Pseudorandom nucleotide replacements within IUPAC specifications, e.g. R: either A or G. Case specific, r: either a or g.    
    fn replace_iupac(&self, mut rng: ThreadRng) -> Self {
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
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
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
                _ => ch,
            })
            .collect::<String>().into_bytes()
    }
    /// fill all other than ACGTUactgu with pseudorandom nucleotides ACTGU. Should be used last after other functions or for cleanup of unknown characters.
    fn replace_all_other_with_uppercase(&self, mut rng: ThreadRng) -> Self {
        let bases = ['A', 'C', 'T', 'G'];
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
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
                _ => *bases.choose(&mut rng).unwrap(),
            })
            .collect::<String>().into_bytes()
    }
    /// fill all other than ACGTUactgu with pseudorandom nucleotides actgu. Should be used last after other functions or for cleanup of unknown characters.
    fn replace_all_other_with_lowercase(&self, mut rng: ThreadRng) -> Self {
        let bases = ['a', 'c', 't', 'g'];
        std::string::String::from_utf8(self.to_owned()).expect("Not valid UTF8").chars()
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
                _ => *bases.choose(&mut rng).unwrap(),
            })
            .collect::<String>().into_bytes()
    }
}



#[test]
fn test_replace_u_with_t() {
    let mut rng = rand::thread_rng();
    let test = "ACTGuyUNn.-@<^>".to_string().replace_u_with_t();
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".to_string().len());
}

#[test]
fn test_replace_gap() {
    let mut rng = rand::thread_rng();
    let test = "ACTGuyUNn.-@<^>".to_string().replace_gap(rng);
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".to_string().len());
}

#[test]
fn test_replace_gap_lowercase() {
    let mut rng = rand::thread_rng();
    let test = "ACTGuyUNn.-@<^>".to_string().replace_gap_lowercase(rng);
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".to_string().len());
}

#[test]
fn test_replace_n() {
    let mut rng = rand::thread_rng();
    let test = "ACTGuyUNn.-@<^>".to_string().replace_n(rng);
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".to_string().len());
}

#[test]
fn test_replace_iupac() {
    let mut rng = rand::thread_rng();
    let test = "ACTGuyUNn.-@<^>".to_string().replace_iupac(rng);
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".to_string().len());
}

#[test]
fn test_replace_all_other_with_uppercase() {
    let mut rng = rand::thread_rng();
    let test = "ACTGuyUNn.-@<^>".to_string().replace_all_other_with_uppercase(rng);
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".to_string().len());
}

#[test]
fn test_rreplace_all_other_with_lowercase() {
    let mut rng = rand::thread_rng();
    let test = "ACTGuyUNn.-@<^>".to_string().replace_all_other_with_lowercase(rng);
    assert_eq!(test.len(), "ACTGuyUNn.-@<^>".to_string().len());
}
