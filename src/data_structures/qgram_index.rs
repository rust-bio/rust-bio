// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::num::{Int, UnsignedInt, NumCast, cast, Float};
use std::collections;
use std::slice;
use std;

use alphabets::{Alphabet, RankTransform};
use utils;


/// Iterator over the q-grams of a given text. Q-grams are encoded as integers.
/// The number of bits for encoding a single symbol is chosen as log2(A) with A being the alphabet
/// size.
///
/// The type Q has to be chosen such that the q-gram fits into it.
pub struct QGrams<'a, Q: UnsignedInt + NumCast> {
    text: &'a [u8],
    bits: usize,
    mask: Q,
    ranks: RankTransform,
}


impl<'a, Q: UnsignedInt + NumCast> QGrams<'a, Q> {
    /// Create new instance.
    ///
    /// # Arguments
    ///
    /// * `q` - the length of the q-gram
    /// * `text` - the text
    /// * `alphabet` - the alphabet to use
    pub fn new(q: usize, text: &'a [u8], alphabet: &Alphabet) -> Self {
        let ranks = RankTransform::new(alphabet);
        let bits = (alphabet.len() as f32).log2().ceil() as usize;
        assert!(bits * q <= Q::
        QGrams {
            text: text,
            ranks: ranks,
            bits: bits,
            mask: cast((1 << q * bits) - 1).unwrap(),
        }
    }

    pub fn iter(&self) -> QGramIter {
        let mut iter = QGramIter { qgrams: self, text: text.iter(), qgram: cast(0).unwrap() };
        for _ in 0..q-1 {
            iter.next();
        }

        iter
    }
}


pub struct QGramIter<'a, Q: UnsignedInt + NumCast> {
    qgrams: &'a QGrams,
    text: slice::Iter<'a, u8>,
    qgram: Q,
}


impl<'a, Q: UnsignedInt + NumCast> QGramIter<'a, Q> {
    fn qgram_push(&mut self, a: u8) {
        self.qgram = self.qgram << self.qgrams.bits;
        self.qgram = (self.qgram | cast(a).unwrap()) & self.qgrams.mask;
    }
}


impl<'a, Q: UnsignedInt + NumCast> Iterator for QGramIter<'a, Q> {
    type Item = Q;

    fn next(&mut self) -> Option<Q> {
        match self.text.next() {
            Some(a) => {
                let b = self.qgrams.ranks.get(*a);
                self.qgram_push(b);
                Some(self.qgram)
            },
            None    => None
        }
    }
}


pub struct QGramIndex<'a> {
    q: usize,
    alphabet: &'a Alphabet,
    address: Vec<usize>,
    pos: Vec<usize>,
}


impl<'a> QGramIndex<'a> {
    pub fn new(q: usize, text: &[u8], alphabet: &'a Alphabet) -> Self {
        QGramIndex::with_max_count(q, text, alphabet, std::usize::MAX)
    }

    pub fn with_max_count(q: usize, text: &[u8], alphabet: &'a Alphabet, max_count: usize) -> Self {
        let qgram_count = alphabet.len().pow(q as u32);
        let mut address = vec![0; qgram_count + 1];
        let mut pos = vec![0; text.len()];

        for qgram in QGrams::<u32>::new(q, text, alphabet) {
            address[qgram as usize] += 1;
        }

        for g in 1..address.len() {
            if address[g] > max_count {
                // mask qgram
                address[g] = 0;
            }
        }

        utils::prescan(&mut address, 0, |a, b| a + b);

        {
            let mut offset = vec![0; qgram_count];
            for (i, qgram) in QGrams::<u32>::new(q, text, alphabet).enumerate() {
                let a = address[qgram as usize];
                if address[qgram as usize + 1] - a != 0 {
                    // if not masked, insert positions
                    pos[a + offset[qgram as usize]] = i;
                    offset[qgram as usize] += 1;
                }
            }
        }

        QGramIndex { q: q, alphabet: alphabet, address: address, pos: pos }
    }

    pub fn matches(&self, qgram: u32) -> &[usize] {
        &self.pos[self.address[qgram as usize]..self.address[qgram as usize + 1]]
    }

    pub fn diagonals(&self, pattern: &[u8]) -> Vec<Diagonal> {
        let mut diagonals = collections::HashMap::new();
        for (i, qgram) in QGrams::<u32>::new(self.q, pattern, self.alphabet).enumerate() {
            for p in self.matches(qgram) {
                let diagonal = p - i;
                if !diagonals.contains_key(&diagonal) {
                    diagonals.insert(diagonal, 1);
                }
                else {
                    *diagonals.get_mut(&diagonal).unwrap() += 1;
                }
            }
        }
        diagonals.into_iter().map(|(diagonal, count)| Diagonal { pos: diagonal, count: count }).collect()
    }

    pub fn exact_matches(&self, pattern: &[u8]) -> Vec<ExactMatch> {
        let mut diagonals: collections::HashMap<usize, ExactMatch> = collections::HashMap::new();
        let mut intervals = Vec::new();
        for (i, qgram) in QGrams::<u32>::new(self.q, pattern, self.alphabet).enumerate() {
            for &p in self.matches(qgram) {
                let diagonal = p - i;
                if !diagonals.contains_key(&diagonal) {
                    // nothing yet, start new match
                    diagonals.insert(diagonal, ExactMatch {
                            pattern: Interval{ start: i, stop: i + self.q},
                            text: Interval { start: p, stop: p + self.q },
                    });
                }
                else {
                    let interval = diagonals.get_mut(&diagonal).unwrap();
                    if interval.pattern.stop - self.q + 1 == i {
                        // extend exact match
                        interval.pattern.stop = i + self.q;
                        interval.text.stop = p + self.q;
                    }
                    else {
                        // report previous match
                        intervals.push(interval.clone());
                        // mismatch or indel, start new match
                        interval.pattern.start = i;
                        interval.pattern.stop = i + self.q;
                        interval.text.start = p;
                        interval.text.stop = p + self.q;
                    }

                }
            }
        }
        // report remaining intervals
        for (_, interval) in diagonals.into_iter() {
            intervals.push(interval);
        }
        intervals
    }
}


#[derive(PartialEq)]
#[derive(Debug)]
pub struct Diagonal {
    pub pos: usize,
    pub count: usize,
}


/// An interval, consisting of start and stop position (the latter exclusive).
#[derive(Clone)]
#[derive(PartialEq)]
#[derive(Debug)]
pub struct Interval {
    pub start: usize,
    pub stop: usize
}


impl Interval {
    pub fn get<'a>(&self, text: &'a [u8]) -> &'a [u8] {
        &text[self.start..self.stop]
    }
}


#[derive(Clone)]
#[derive(PartialEq)]
#[derive(Debug)]
pub struct ExactMatch {
    pub pattern: Interval,
    pub text: Interval,
}


#[cfg(test)]
mod tests {
    use super::*;
    use alphabets;

    fn setup() -> (&'static [u8], alphabets::Alphabet) {
        let text = b"ACGGCTGAGATGAT";
        let alphabet = alphabets::dna::alphabet();

        (text, alphabet)
    }

    #[test]
    fn test_matches() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let qgram = QGrams::new(q, b"TGA", &alphabet).next().unwrap();

        let matches = qgram_index.matches(qgram);
        assert_eq!(matches, [5, 10]);
    }

    #[test]
    fn test_diagonals() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let pattern = b"GCTG";
        let diagonals = qgram_index.diagonals(pattern);
        assert_eq!(diagonals, [Diagonal { pos: 3, count: 2 }]);
    }

    #[test]
    fn test_exact_matches() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let pattern = b"GCTGA";
        let exact_matches = qgram_index.exact_matches(pattern);
        assert!(exact_matches.len() == 2);
        for m in exact_matches {
            assert_eq!(m.pattern.get(pattern), m.text.get(text));
        }
    }
}
