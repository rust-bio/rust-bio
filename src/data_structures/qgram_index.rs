// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A classical, flexible, q-gram index implementation.
//!
//! # Example
//!
//! ```
//! use bio::alphabets;
//! use bio::data_structures::qgram_index;
//!
//! let text = b"ACGGCTGAGATGAT";
//! let alphabet = alphabets::dna::alphabet();
//! let q = 3;
//! let qgram_index = qgram_index::QGramIndex::new(q, text, &alphabet);
//!
//! let pattern = b"GCTG";
//! let matches = qgram_index.matches(pattern, 1);
//! assert_eq!(
//!     matches,
//!     [qgram_index::Match {
//!         pattern: qgram_index::Interval { start: 0, stop: 4 },
//!         text: qgram_index::Interval { start: 3, stop: 7 },
//!         count: 2
//!     }]
//! );
//! ```

use std::cmp;
use std::collections;
use std::collections::hash_map::Entry;

use crate::alphabets::{Alphabet, RankTransform};
use crate::utils;

/// A classical, flexible, q-gram index implementation.
///
/// Uses |alphabet|^q + k words of memory, where k is the number of q-grams in the text with count at most `max_count` (if specified).
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct QGramIndex {
    q: u32,
    // For each q-gram, the position in `pos` where positions for this q-gram are stored.
    address: Vec<usize>,
    // The positions in `text` where each q-gram occurs.
    pos: Vec<usize>,
    ranks: RankTransform,
}

impl QGramIndex {
    /// Create a new q-gram index.
    /// The q has to be smaller than b / log2(|A|) with |A| being the alphabet size and b the number
    /// bits with the `usize` data type.
    pub fn new<'a, T, I>(q: u32, text: T, alphabet: &Alphabet) -> Self
    where
        I: Iterator<Item = &'a u8> + ExactSizeIterator + Clone,
        T: IntoIterator<Item = &'a u8, IntoIter = I> + Sized,
    {
        QGramIndex::with_max_count(q, text, alphabet, usize::MAX)
    }

    /// Create a new q-gram index, only considering q-grams that occur at most `max_count` times.
    /// The q has to be smaller than b / log2(|A|) with |A| being the alphabet size and b the number
    /// bits with the `usize` data type.
    pub fn with_max_count<'a, T, I>(q: u32, text: T, alphabet: &Alphabet, max_count: usize) -> Self
    where
        I: Iterator<Item = &'a u8> + ExactSizeIterator + Clone,
        T: IntoIterator<Item = &'a u8, IntoIter = I> + Sized,
    {
        let text = text.into_iter();
        let ranks = RankTransform::new(alphabet);

        let qgram_count = alphabet.len().pow(q);
        let mut address = vec![0; qgram_count + 1];

        for qgram in ranks.qgrams(q, text.clone()) {
            address[qgram] += 1;
        }

        for a in address.iter_mut() {
            if *a > max_count {
                // mask qgram
                *a = 0;
            }
        }

        utils::prescan(&mut address, 0, |a, b| a + b);

        // Address has at least size 1, so unwrap is fine.
        let mut pos = vec![0; *address.last().unwrap()];

        {
            let mut offset = vec![0; qgram_count];
            for (i, qgram) in ranks.qgrams(q, text).enumerate() {
                let a = address[qgram];
                if address[qgram + 1] - a != 0 {
                    // if not masked, insert positions
                    pos[a + offset[qgram]] = i;
                    offset[qgram] += 1;
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

    /// Return matches of the given pattern, matching in at least `min_count` q-grams.
    /// Complexity O(m + k) for pattern of length m and k being the number of matching q-grams.
    ///
    /// A match is a substring of `pattern` and a corresponding substring of the text that share at least `min_count` q-grams.
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
    ///
    /// An exact match is a substring of `pattern` occurring in the text of length at least `q`.
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
                        if m.pattern.stop - q + 1 != i {
                            // discontinue match
                            matches.push(*m);
                            // start new match
                            m.pattern.start = i;
                            m.text.start = p;
                        }
                        m.pattern.stop = i + q;
                        m.text.stop = p + q;
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
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
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
#[derive(Default, Copy, Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
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
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct ExactMatch {
    pub pattern: Interval,
    pub text: Interval,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabets;

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
    fn test_qgram_with_max_count() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::with_max_count(q, text, &alphabet, 1);

        let ranks = alphabets::RankTransform::new(&alphabet);

        let qgram = ranks.qgrams(q, b"TGA").next().unwrap();

        // Should be pruned because the count of 2 is larger than the max_count of 1.
        let matches = qgram_index.qgram_matches(qgram);
        assert_eq!(matches, []);
    }

    #[test]
    fn test_qgram_with_max_count_index_0() {
        let (_, alphabet) = setup();
        let text = b"AAAAA";
        let q = 3;
        let qgram_index = QGramIndex::with_max_count(q, text, &alphabet, 1);

        let ranks = alphabets::RankTransform::new(&alphabet);

        let qgram = ranks.qgrams(q, b"AAA").next().unwrap();

        // Should be pruned because the count of 3 is larger than the max_count of 1.
        let matches = qgram_index.qgram_matches(qgram);
        assert_eq!(matches, []);
    }

    #[test]
    fn test_qgram_sizeof_pos() {
        let (_, alphabet) = setup();
        let text = b"AAAAA";
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let ranks = alphabets::RankTransform::new(&alphabet);

        let qgram = ranks.qgrams(q, b"AAA").next().unwrap();

        // Should be pruned because the count of 3 is larger than the max_count of 1.
        let matches = qgram_index.qgram_matches(qgram);
        assert_eq!(matches, [0, 1, 2]);
    }

    #[test]
    fn test_matches() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        // A fully matching pattern.
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

        // A pattern that matches in one position on two disjoint q-grams.
        let pattern = b"GCTAAGA";
        let matches = qgram_index.matches(pattern, 2);
        assert_eq!(
            matches,
            [Match {
                pattern: Interval { start: 0, stop: 7 },
                text: Interval { start: 3, stop: 10 },
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

        // A pattern that matches in one position on two disjoint q-grams.
        let pattern = b"GCTAAGA";
        let matches = qgram_index.exact_matches(pattern);
        assert_eq!(
            matches,
            [
                ExactMatch {
                    pattern: Interval { start: 0, stop: 3 },
                    text: Interval { start: 3, stop: 6 }
                },
                ExactMatch {
                    pattern: Interval { start: 4, stop: 7 },
                    text: Interval { start: 7, stop: 10 }
                }
            ]
        );
    }

    #[test]
    fn test_exact_matches_self() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text, &alphabet);

        let exact_matches = qgram_index.exact_matches(text);
        assert!(!exact_matches.is_empty());
    }

    #[test]
    fn test_iterator() {
        let (text, alphabet) = setup();
        let q = 3;
        let qgram_index = QGramIndex::new(q, text.iter(), &alphabet);

        let exact_matches = qgram_index.exact_matches(text);
        assert!(!exact_matches.is_empty());
    }

    #[test]
    // #[cfg(feature = "nightly")]
    fn test_serde() {
        use serde::{Deserialize, Serialize};
        fn impls_serde_traits<'a, S: Serialize + Deserialize<'a>>() {}

        impls_serde_traits::<QGramIndex>();
    }
}
