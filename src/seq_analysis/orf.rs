// Copyright 2014-2016 Johannes Köster, Martin Larralde.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! One-way orf finder algorithm.
//! complexity: O(n).

//! # Example
//!
//! ```
//! use bio::seq_analysis::orf::NaiveFinder;
//! let start_codons = vec!(b"ATG");
//! let stop_codons  = vec!(b"TGA", b"TAG", b"TAA");
//! let min_len:usize = 50;
//! let finder = NaiveFinder::new(start_codons, stop_codons, min_len);
//!
//! let sequence = b"ACGGCTAGAAAAGGCTAGAAAA"
//!
//! for orf in finder.find_all(sequence) {
//!    ...do something...
//! }
//! ```

//! Right now the only way to check the reverse strand for orf is to use
//! the alphabet::dna::RevComp struct and to check for both sequences.
//! But that's not so performance friendly, as the reverse complementation and the orf research
//! could go on at the same time.


/// Finder algorithm
pub struct NaiveFinder <'a> {
    start_codons: Vec<&'a [u8; 3]>,
    stop_codons: Vec<&'a [u8; 3]>,
    min_len: usize,
}

impl<'a> NaiveFinder <'a> {
    /// Create a new instance for given genetic datas
    pub fn new (start_codons: Vec<&'a [u8; 3]>, stop_codons: Vec<&'a [u8; 3]>, min_len: usize) -> Self {
        NaiveFinder {
            start_codons: start_codons,
            stop_codons: stop_codons,
            min_len: min_len,
        }
    }
    /// Find all open reading frames in given sequence. Matches are returned as iterators over vectors
    pub fn find_all(&'a self, sequence: &'a [u8]) -> Matches {
        Matches::new(self, sequence)
    }
}

/// container struct for argument iteration
struct State {
    index: usize,
    length: usize,
    start_pos: [usize; 3],
    end_pos: [usize; 3],
    in_orf: [bool; 3]
}

impl State {
    /// Create new state
    pub fn new(length: usize) -> Self {
        State {
            index: 0,
            length: length,
            start_pos: [0, 0, 0],
            end_pos: [0, 0, 0],
            in_orf: [false, false, false],
        }
    }
}

/// Iterator yielding open reading frames as vector (start_pos, end_pos + 1, sequence)
pub struct Matches<'a> {
    finder: &'a NaiveFinder<'a>,
    sequence: &'a [u8],
    state: State,
}

impl<'a> Matches<'a> {
    pub fn new (finder: &'a NaiveFinder, sequence: &'a [u8]) -> Self {

        //let mut v = Vec::new();
        //v.extend_from_slice(sequence);

        Matches {
            finder: finder,
            sequence: sequence,
            state: State::new(sequence.len()),
        }
    }
}


impl<'a> Iterator for Matches<'a> {
    type Item = (usize, usize, &'a [u8]);

    fn next(&mut self) -> Option<(usize, usize, &'a [u8])> {

        let mut orf: Option<(usize, usize, &'a [u8])> = None;
        let mut codons: [ [u8; 3]; 3]; //

        while self.state.index < self.state.length - 3 {
            //codons array depending on reading frame
            codons = [
                [ self.sequence[self.state.index]  , self.sequence[self.state.index+1], self.sequence[self.state.index+2] ],
                [ self.sequence[self.state.index+1], self.sequence[self.state.index+2], self.sequence[self.state.index+3] ],
                [ self.sequence[self.state.index+2], self.sequence[self.state.index+3], self.sequence[self.state.index+4] ],
            ];
            //x is the shift of the reading frame
            for x in 0..3 as usize {

                if self.state.in_orf[x] {
                    if self.finder.stop_codons.contains(&&codons[x]) {
                        //exiting orf
                        self.state.end_pos[x] = self.state.index + x + 3;
                        self.state.in_orf[x] = false;
                        //slice the sequence to get the frame
                        let slice = &self.sequence[self.state.start_pos[x]..self.state.end_pos[x]];

                        // (!) still don't know if it's better to return slice or vec
                        // if vec is better, do:
                        //   let mut frame = Vec::new();
                        //   frame.extend_from_slice(slice);
                        // and return Some(.., frame) instead of Some(.., slice)

                        //check if orf length is enough
                        if slice.len() > self.finder.min_len {
                            orf = Some( (self.state.start_pos[x], self.state.end_pos[x], slice) );
                        }
                        self.state.start_pos[x] = 0;
                        self.state.end_pos[x] = 0
                    }
                } else {
                    if self.finder.start_codons.contains(&&codons[x]) {
                        //entering orf
                        self.state.start_pos[x] = self.state.index + x;
                        self.state.in_orf[x] = true;
                    }
                }

            }
            self.state.index += 3;
            if orf.is_some() {
                return orf;
            }
        }
        None
    }
}
