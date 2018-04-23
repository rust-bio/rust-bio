// Copyright 2014-2016 Johannes Köster, Martin Larralde.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! One-way orf finder algorithm.
//!
//! Complexity: O(n).
//!
//! # Example
//!
//! ```
//! use bio::seq_analysis::orf::{Finder, Orf};
//! let start_codons = vec!(b"ATG");
//! let stop_codons  = vec!(b"TGA", b"TAG", b"TAA");
//! let min_len = 50;
//! let finder = Finder::new(start_codons, stop_codons, min_len);
//!
//! let sequence = b"ACGGCTAGAAAAGGCTAGAAAA";
//!
//! for Orf{start, end, offset} in finder.find_all(sequence) {
//!    let orf = &sequence[start..end];
//!    //...do something with orf sequence...
//! }
//! ```
//!
//! Right now the only way to check the reverse strand for orf is to use
//! the `alphabet::dna::RevComp` struct and to check for both sequences.
//! But that's not so performance friendly, as the reverse complementation and the orf research
//! could go on at the same time.

use std::collections::VecDeque;
use std::iter;
use utils::{IntoTextIterator, TextIterator};

/// An implementation of a naive algorithm finder
pub struct Finder {
    start_codons: Vec<VecDeque<u8>>,
    stop_codons: Vec<VecDeque<u8>>,
    min_len: usize,
}

impl Finder {
    /// Create a new instance of a finder for the given start and stop codons and a particular length
    pub fn new<'a>(
        start_codons: Vec<&'a [u8; 3]>,
        stop_codons: Vec<&'a [u8; 3]>,
        min_len: usize,
    ) -> Self {
        Finder {
            start_codons: start_codons.into_iter()                          // Convert start_ and
                                      .map(|x| {                            // stop_codons from
                                          x.into_iter()                     // Vec<&[u8;3]> to
                                           .map(|&x| x as u8)               // Vec<VecDeque<u8>>
                                           .collect::<VecDeque<u8>>()       // so they can be
                                      })                                    // easily compared
                                      .collect(), // with codon built
            stop_codons: stop_codons.into_iter()                            // from IntoTextIterator
                                    .map(|x| {                              // object.
                                        x.into_iter()
                                         .map(|&x| x as u8)
                                         .collect::<VecDeque<u8>>()
                                    })
                                    .collect(),
            min_len,
        }
    }

    /// Find all orfs in the given sequence
    pub fn find_all<'a, I: IntoTextIterator<'a>>(&'a self, seq: I) -> Matches<I::IntoIter> {
        Matches {
            finder: self,
            state: State::new(),
            seq: seq.into_iter().enumerate(),
        }
    }
}

/// An orf representation with start and end position of said orf,
/// as well as offset of the reading frame (1,2,3) and strand location
// (current: +, reverse complementary: -).
pub struct Orf {
    pub start: usize,
    pub end: usize,
    pub offset: i8,
}

/// The current algorithm state.
struct State {
    start_pos: [Option<usize>; 3],
    codon: VecDeque<u8>,
}

impl State {
    /// Create new state.
    pub fn new() -> Self {
        State {
            start_pos: [None, None, None],
            codon: VecDeque::new(),
        }
    }
}

/// Iterator over offset, start position, end position and sequence of matched orfs.
pub struct Matches<'a, I: TextIterator<'a>> {
    finder: &'a Finder,
    state: State,
    seq: iter::Enumerate<I>,
}

impl<'a, I: Iterator<Item = &'a u8>> Iterator for Matches<'a, I> {
    type Item = Orf;

    fn next(&mut self) -> Option<Orf> {
        let mut result: Option<Orf> = None;
        let mut offset: usize;

        for (index, &nuc) in self.seq.by_ref() {
            // update the codon
            if self.state.codon.len() >= 3 {
                self.state.codon.pop_front();
            }
            self.state.codon.push_back(nuc);
            offset = (index + 1) % 3;

            // inside orf
            if self.state.start_pos[offset].is_some() {
                // check if leaving orf
                if self.finder.stop_codons.contains(&self.state.codon) {
                    // check if length is sufficient
                    if index + 1 - self.state.start_pos[offset].unwrap() > self.finder.min_len {
                        // build results
                        result = Some(Orf {
                            start: self.state.start_pos[offset].unwrap() - 2,
                            end: index + 1,
                            offset: offset as i8,
                        });
                    }
                    // reinitialize
                    self.state.start_pos[offset] = None;
                }
            // check if entering orf
            } else if self.finder.start_codons.contains(&self.state.codon) {
                self.state.start_pos[offset] = Some(index);
            }
            if result.is_some() {
                return result;
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orf() {
        let start_codons = vec![b"ATG"];
        let stop_codons = vec![b"TGA", b"TAG", b"TAA"];
        let min_len = 50;
        let finder = Finder::new(start_codons, stop_codons, min_len);

        let sequence = b"ACGGCTAGAAAAGGCTAGAAAA";

        for Orf { start, end, .. } in finder.find_all(sequence) {
            let _ = &sequence[start..end];
        }
    }
}
