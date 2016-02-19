// Copyright 2014-2015 Johannes Köster, Peer Aramillo Irizar.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Implementation of alphabets and useful utilities.
//!
//! # Example
//!
//! ```rust
//! use bio::alphabets;
//! let alphabet = alphabets::dna::alphabet();
//! assert!(alphabet.is_word(b"AACCTgga"));
//! assert!(!alphabet.is_word(b"AXYZ"));
//! ```

use std::iter::FromIterator;
use std::slice;
use std::mem;

use bit_set::BitSet;
use vec_map::VecMap;

pub mod dna;
pub mod protein;

pub type SymbolRanks = VecMap<u8>;


/// Representation of an alphabet.
pub struct Alphabet {
    pub symbols: BitSet
}


impl Alphabet {
    /// Create new alphabet from given symbols.
    pub fn new(symbols: &[u8]) -> Self {
        Alphabet::from_iter(symbols.iter())
    }

    /// Insert symbol into alphabet.
    pub fn insert(&mut self, a: u8) {
        self.symbols.insert(a as usize);
    }

    /// Check if given text is a word over the alphabet.
    pub fn is_word(&self, text: &[u8]) -> bool {
        text.iter().all(|&c| self.symbols.contains(c as usize))
    }

    /// Return lexicographically maximal symbol.
    pub fn max_symbol(&self) -> Option<u8> {
        self.symbols.iter().max().map(|a| a as u8)
    }

    /// Return size of the alphabet.
    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    /// Is this alphabet empty?
    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }
}

impl<'a> FromIterator<&'a u8> for Alphabet {
    /// Create a new alphabet from the given iterator.
    fn from_iter<T: IntoIterator<Item=&'a u8>>(iterator: T) -> Self {
        let mut s = BitSet::new();
        s.extend(iterator.into_iter().map(|&c| c as usize));

        Alphabet { symbols: s }
    }
}

/// Tools based on transforming the alphabet symbols to their lexicographical ranks.
#[cfg_attr(feature = "serde_macros", derive(Serialize, Deserialize))]
pub struct RankTransform {
    pub ranks: SymbolRanks
}


impl RankTransform {
    /// Construct a new `RankTransform`.
    pub fn new(alphabet: &Alphabet) -> Self {
        let mut ranks = VecMap::new();
        for (r, c) in alphabet.symbols.iter().enumerate() {
            ranks.insert(c, r as u8);
        }

        RankTransform { ranks: ranks }
    }

    /// Get the rank of symbol `a`.
    pub fn get(&self, a: u8) -> u8 {
        *self.ranks.get(a as usize).expect("Unexpected character.")
    }

    /// Transform a given `text`.
    pub fn transform(&self, text: &[u8]) -> Vec<u8> {
        text.iter()
            .map(|&c| *self.ranks.get(c as usize).expect("Unexpected character in text."))
            .collect()
    }

    /// Iterate over q-grams (substrings of length q) of given `text`. The q-grams are encoded
    /// as `usize` by storing the symbol ranks in log2(|A|) bits (with |A| being the alphabet size).
    ///
    /// If q is larger than usize::BITS / log2(|A|), this method fails with an assertion.
    pub fn qgrams<'a>(&'a self, q: u32, text: &'a [u8]) -> QGrams {
        let bits = (self.ranks.len() as f32).log2().ceil() as u32;
        assert!((bits * q) as usize <= mem::size_of::<usize>() * 8, "Expecting q to be smaller than usize / log2(|A|)");

        let mut qgrams = QGrams {
            text: text.iter(),
            ranks: self,
            bits: bits,
            mask: (1 << (q * bits)) - 1,
            qgram: 0,
        };

        for _ in 0..q-1 {
            qgrams.next();
        }

        qgrams
    }

    /// Restore alphabet from transform.
    pub fn alphabet<'a>(&self) -> Alphabet {
        let mut symbols = BitSet::with_capacity(self.ranks.len());
        symbols.extend(self.ranks.keys());
        Alphabet { symbols: symbols }
    }
}


/// Iterator over q-grams.
pub struct QGrams<'a> {
    text: slice::Iter<'a, u8>,
    ranks: &'a RankTransform,
    bits: u32,
    mask: usize,
    qgram: usize,
}


impl<'a> QGrams<'a> {
    fn qgram_push(&mut self, a: u8) {
        self.qgram <<= self.bits;
        self.qgram |= a as usize;
        self.qgram &= self.mask;
    }
}


impl<'a> Iterator for QGrams<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a);
                self.qgram_push(b);
                Some(self.qgram)
            },
            None    => None
        }
    }
}

#[cfg(tests)]
mod tests {
    #[test]
    #[cfg(feature = "nightly")]
    fn test_serde() {
        use serde::{Serialize, Deserialize};
        fn impls_serde_traits<S: Serialize + Deserialize>() {}

        impls_serde_traits::<RankTransform>();
    }
}
