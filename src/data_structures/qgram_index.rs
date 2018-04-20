// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A classical, flexible, q-gram index implementation.
//!
//! # Example
//!
//! ```
//! use bio::data_structures::qgram_index;
//! use bio::alphabets;
//!
//! let text = b"ACGGCTGAGATGAT";
//! let alphabet = alphabets::dna::alphabet();
//! let q = 3;
//! let qgram_index = qgram_index::QGramIndex::new(q, text, &alphabet);
//!
//! let pattern = b"GCTG";
//! let matches = qgram_index.matches(pattern, 1);
//! assert_eq!(matches, [
//!     qgram_index::Match {
//!         pattern: qgram_index::Interval { start: 0, stop: 4 },
//!         text: qgram_index::Interval { start: 3, stop: 7 },
//!         count: 2
//!     }
//! ]);
//! ```

use std;
use std::cmp;
use std::collections;
use std::collections::hash_map::Entry;

use alphabets::{Alphabet, RankTransform};
use utils;

/// A classical, flexible, q-gram index implementation.
#[derive(Serialize, Deserialize)]
pub struct QGramIndex {
    q: u32,
    address: Vec<usize>,
    pos: Vec<usize>,
    ranks: RankTransform,
}

impl QGramIndex {
    /// Create a new q-gram index.
    /// The q has to be smaller than b / log2(|A|) with |A| being the alphabet size and b the number
    /// bits with the `usize` data type.
    pub fn new(q: u32, text: &[u8], alphabet: &Alphabet) -> Self {
        QGramIndex::with_max_count(q, text, alphabet, std::usize::MAX)
    }

    /// Create a new q-gram index, only considering q-grams that occur at most `max_count` times.
    /// The q has to be smaller than b / log2(|A|) with |A| being the alphabet size and b the number
    /// bits with the `usize` data type.
    pub fn with_max_count(q: u32, text: &[u8], alphabet: &Alphabet, max_count: usize) -> Self {
        let ranks = RankTransform::new(alphabet);

        let qgram_count = alphabet.len().pow(q as u32);
        let mut address = vec![0; qgram_count + 1];
        let mut pos = vec![0; text.len()];

        for qgram in ranks.qgrams(q, text) {
            address[qgram] += 1;
        }

        for a in address.iter_mut().skip(1) {
            if *a > max_count {
                // mask qgram
                *a = 0;
            }
        }

        utils::prescan(&mut address, 0, |a, b| a + b);

        {
            let mut offset = vec![0; qgram_count];
            for (i, qgram) in ranks.qgrams(q, text).enumerate() {
                let a = address[qgram as usize];
                if address[qgram as usize + 1] - a != 0 {
                    // if not masked, insert positions
                    pos[a + offset[qgram as usize]] = i;
                    offset[qgram as usize] += 1;
                }
            }
        }

        QGramIndex {
            q,
            address,
            pos,
            ranks,
        }
    }

    /// The used q.
    pub fn q(&self) -> u32 {
        self.q
    }

    /// Return text positions with matching q-gram. Complexity O(1).
    pub fn qgram_matches(&self, qgram: usize) -> &[usize] {
        &self.pos[self.address[qgram]..self.address[qgram + 1]]
    }

    /// Return matches of the given pattern.
    /// Complexity O(m + k) for pattern of length m and k being the number of matching q-grams.
    pub fn matches(&self, pattern: &[u8], min_count: usize) -> Vec<Match> {
        let q = self.q as usize;
        let mut diagonals = collections::HashMap::new();
        for (i, qgram) in self.ranks.qgrams(self.q, pattern).enumerate() {
            for &p in self.qgram_matches(qgram) {
                let diagonal = p - i;
                match diagonals.entry(diagonal) {
                    Entry::Vacant(v) => {
                        v.insert(Match {
                            pattern: Interval {
                                start: i,
                                stop: i + q,
                            },
                            text: Interval {
                                start: p,
                                stop: p + q,
                            },
                            count: 1,
                        });
                    }
                    Entry::Occupied(mut o) => {
                        let m = o.get_mut();
                        m.pattern.stop = i + q;
                        m.text.stop = p + q;
                        m.count += 1;
                    }
                }
            }
        }
        diagonals
            .into_iter()
            .filter_map(|(_, m)| if m.count >= min_count { Some(m) } else { None })
            .collect()
    }

    /// Return exact matches (substrings) of the given pattern.
    /// Complexity O(m + k) for pattern of length m and k being the number of matching q-grams.
    pub fn exact_matches(&self, pattern: &[u8]) -> Vec<ExactMatch> {
        let q = self.q as usize;
        let mut diagonals = collections::HashMap::new();
        let mut matches = Vec::new();

        for (i, qgram) in self.ranks.qgrams(self.q, pattern).enumerate() {
            for &p in self.qgram_matches(qgram) {
                let diagonal = p as i32 - i as i32;
                match diagonals.entry(diagonal) {
                    Entry::Vacant(v) => {
                        v.insert(ExactMatch {
                            pattern: Interval {
                                start: i,
                                stop: i + q,
                            },
                            text: Interval {
                                start: p,
                                stop: p + q,
                            },
                        });
                    }
                    Entry::Occupied(mut o) => {
                        let m = o.get_mut();
                        if m.pattern.stop - q + 1 == i {
                            m.pattern.stop = i + q;
                            m.text.stop = p + q;
                        } else {
                            // discontinue match
                            matches.push(*m);
                            // start new match
                            m.pattern.start = i;
                            m.pattern.stop = i + q;
                            m.text.start = p;
                            m.text.stop = p + q;
                        }
                    }
                }
            }
        }
        for (_, m) in diagonals {
            matches.push(m);
        }

        matches
    }
}

/// An interval, consisting of start and stop position (the latter exclusive).
#[derive(PartialEq, Eq, Debug, Copy, Clone)]
pub struct Interval {
    pub start: usize,
    pub stop: usize,
}

impl Interval {
    /// Get the text within the given interval.
    pub fn get<'a>(&self, text: &'a [u8]) -> &'a [u8] {
        &text[self.start..self.stop]
    }
}

/// A match between the pattern and the text.
#[derive(PartialEq, Eq, Debug, Copy, Clone)]
pub struct Match {
    pub pattern: Interval,
    pub text: Interval,
    pub count: usize,
}

impl cmp::Ord for Match {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.count.cmp(&other.count)
    }
}

impl cmp::PartialOrd for Match {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// An exact match between the pattern and the text.
#[derive(PartialEq, Debug, Copy, Clone)]
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
    fn test_qgram_matches() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let ranks = alphabets::RankTransform::new(&alphabet);

        let qgram = ranks.qgrams(q, b"TGA").next().unwrap();

        let matches = qgram_index.qgram_matches(qgram);
        assert_eq!(matches, [5, 10]);
    }

    #[test]
    fn test_matches() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let pattern = b"GCTG";
        let matches = qgram_index.matches(pattern, 1);
        assert_eq!(
            matches,
            [Match {
                pattern: Interval { start: 0, stop: 4 },
                text: Interval { start: 3, stop: 7 },
                count: 2,
            }]
        );
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

    #[test]
    fn test_exact_matches_self() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let exact_matches = qgram_index.exact_matches(text);
        assert!(exact_matches.len() >= 1);
    }

    #[test]
    #[cfg(feature = "nightly")]
    fn test_serde() {
        use serde::{Deserialize, Serialize};
        fn impls_serde_traits<S: Serialize + Deserialize>() {}

        impls_serde_traits::<QGramIndex>();
    }
}
