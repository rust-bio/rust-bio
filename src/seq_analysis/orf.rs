// Copyright 2014-2016 Johannes Köster, Martin Larralde.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! One-way open reading frame (ORF) finder algorithm.
//!
//! Complexity: O(n).
//!
//! # Example
//!
//! ```
//! use bio::seq_analysis::orf::{Finder, Orf};
//! let start_codons = vec![b"ATG"];
//! let stop_codons = vec![b"TGA", b"TAG", b"TAA"];
//! let min_len = 50;
//! let finder = Finder::new(start_codons, stop_codons, min_len);
//!
//! let sequence = b"ACGGCTAGAAAAGGCTAGAAAA";
//!
//! for Orf { start, end, offset } in finder.find_all(sequence) {
//!     let orf = &sequence[start..end];
//!     //...do something with orf sequence...
//! }
//! ```
//!
//! Right now the only way to check the reverse strand for ORF is to use
//! the `alphabet::dna::RevComp` struct and to check for both sequences.
//! But that's not so performance friendly, as the reverse complementation and the ORF research
//! could go on at the same time.

use std::borrow::Borrow;
use std::collections::VecDeque;
use std::iter;

/// An implementation of a naive algorithm finder
// Implementation note:
//
// VecDeque is used rather than the obvious [u8; 3] to represent
// codons because a VecDeque<u8> is used to represent a sliding codon
// (see: State.codon) window which unfortunately, cannot be compared
// to [u8; 3].
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Finder {
    start_codons: Vec<VecDeque<u8>>,
    stop_codons: Vec<VecDeque<u8>>,
    min_len: usize,
}

impl Finder {
    /// Create a new instance of a finder for the given start and stop codons and the minimum
    /// length of an ORF.
    pub fn new<'a>(
        start_codons: Vec<&'a [u8; 3]>,
        stop_codons: Vec<&'a [u8; 3]>,
        min_len: usize,
    ) -> Self {
        Finder {
            start_codons: start_codons
                .iter()
                .map(|x| x.iter().copied().collect::<VecDeque<u8>>())
                .collect(),
            stop_codons: stop_codons
                .iter()
                .map(|x| x.iter().copied().collect::<VecDeque<u8>>())
                .collect(),
            min_len,
        }
    }

    /// Find all ORFs in the given sequence
    pub fn find_all<C, T>(&self, seq: T) -> Matches<'_, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        Matches {
            finder: self,
            state: State::new(),
            seq: seq.into_iter().enumerate(),
        }
    }
}

/// An ORF representation with start and end position of said ORF,
/// as well as offset of the reading frame (1,2,3) and strand location
// (current: +, reverse complementary: -).
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct Orf {
    pub start: usize,
    pub end: usize,
    pub offset: i8,
}

/// The current algorithm state.
#[derive(Clone, Debug)]
struct State {
    start_pos: [Vec<usize>; 3],
    codon: VecDeque<u8>,
    found: VecDeque<Orf>,
}

impl State {
    /// Create new state.
    pub fn new() -> Self {
        State {
            start_pos: [Vec::new(), Vec::new(), Vec::new()],
            codon: VecDeque::new(),
            found: VecDeque::new(),
        }
    }
}

/// Iterator over offset, start position, end position and sequence of matched ORFs.
#[derive(Clone, Debug)]
pub struct Matches<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    finder: &'a Finder,
    state: State,
    seq: iter::Enumerate<T>,
}

impl<'a, C, T> Iterator for Matches<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    type Item = Orf;

    fn next(&mut self) -> Option<Orf> {
        let mut offset: usize;

        // return any orfs already found
        if !self.state.found.is_empty() {
            return self.state.found.pop_front();
        }

        for (index, nuc) in self.seq.by_ref() {
            // update the codon
            if self.state.codon.len() >= 3 {
                self.state.codon.pop_front();
            }
            self.state.codon.push_back(*nuc.borrow());
            offset = (index + 1) % 3;

            // check if entering orf
            if self.finder.start_codons.contains(&self.state.codon) {
                self.state.start_pos[offset].push(index);
            }
            // inside orf
            if !self.state.start_pos[offset].is_empty() {
                // check if leaving orf
                if self.finder.stop_codons.contains(&self.state.codon) {
                    for start_pos in &self.state.start_pos[offset] {
                        // check if length is sufficient
                        if index + 1 - start_pos > self.finder.min_len {
                            // build results
                            self.state.found.push_back(Orf {
                                start: start_pos - 2,
                                end: index + 1,
                                offset: offset as i8,
                            });
                        // if the first orf is too short, so are the others
                        } else {
                            break;
                        }
                    }
                    // reinitialize
                    self.state.start_pos[offset] = Vec::new();
                }
            }
            if !self.state.found.is_empty() {
                return self.state.found.pop_front();
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn basic_finder() -> Finder {
        let start_codons = vec![b"ATG"];
        let stop_codons = vec![b"TGA", b"TAG", b"TAA"];
        let min_len = 5;
        Finder::new(start_codons, stop_codons, min_len)
    }

    #[test]
    fn test_no_orf() {
        let finder = basic_finder();
        let sequence = b"ACGGCTAGAAAAGGCTAGAAAA";
        assert!(finder.find_all(sequence).next().is_none());
    }

    #[test]
    fn test_one_orf_no_offset() {
        let finder = basic_finder();
        let sequence = b"GGGATGGGGTGAGGG";
        let expected = vec![Orf {
            start: 3,
            end: 12,
            offset: 0,
        }];
        assert_eq!(expected, finder.find_all(sequence).collect::<Vec<Orf>>());
    }

    #[test]
    fn test_one_orf_with_offset() {
        let finder = basic_finder();
        let sequence = b"AGGGATGGGGTGAGGG";
        let expected = vec![Orf {
            start: 4,
            end: 13,
            offset: 1,
        }];
        assert_eq!(expected, finder.find_all(sequence).collect::<Vec<Orf>>());
    }

    #[test]
    fn test_two_orfs_different_offsets() {
        let finder = basic_finder();
        let sequence = b"ATGGGGTGAGGGGGATGGAAAAATAAG";
        let expected = vec![
            Orf {
                start: 0,
                end: 9,
                offset: 0,
            },
            Orf {
                start: 14,
                end: 26,
                offset: 2,
            },
        ];
        assert_eq!(expected, finder.find_all(sequence).collect::<Vec<Orf>>());
    }

    #[test]
    fn test_three_nested_and_offset_orfs() {
        let finder = basic_finder();
        let sequence = b"ATGGGGATGGGGGGATGGAAAAATAAGTAG";
        let expected = vec![
            Orf {
                start: 14,
                end: 26,
                offset: 2,
            },
            Orf {
                start: 0,
                end: 30,
                offset: 0,
            },
            Orf {
                start: 6,
                end: 30,
                offset: 0,
            },
        ];
        assert_eq!(expected, finder.find_all(sequence).collect::<Vec<Orf>>());
    }
}
