// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::num::{Int, UnsignedInt, NumCast, cast};
use std::collections;
use std::slice;

use alphabets::{Alphabet, RankTransform};


fn qgram_push<Q: UnsignedInt + NumCast>(qgram: &mut Q, a: u8) {
    *qgram = *qgram << 2; // TODO generalize this for n-bit encoding!
    *qgram = *qgram | cast(a).unwrap();
    // TODO mask unused bits with q
}


struct QGrams<'a, Q: UnsignedInt + NumCast> {
    text: slice::Iter<'a, u8>,
    qgram: Q,
    q: usize,
    ranks: RankTransform,
}


impl<'a, Q: UnsignedInt + NumCast> QGrams<'a, Q> {
    pub fn new(q: usize, text: &'a [u8], alphabet: &Alphabet) -> Self {
        let ranks = RankTransform::new(alphabet);
        let mut qgrams = QGrams { text: text.iter(), qgram: cast(0).unwrap(), q: q, ranks: ranks };
        for _ in 0..q-1 {
            qgrams.next();
        }

        qgrams
    }
}


impl<'a, Q: UnsignedInt + NumCast> Iterator for QGrams<'a, Q> {
    type Item = Q;

    fn next(&mut self) -> Option<Q> {
        match self.text.next() {
            Some(a) => {
                qgram_push(&mut self.qgram, self.ranks.get(*a));
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
        let qgram_count = alphabet.len().pow(q as u32);
        let mut address = vec![0; qgram_count];
        let mut pos = vec![0; text.len()];

        for qgram in QGrams::<u32>::new(q, text, alphabet) {
            address[qgram as usize] += 1;
        }

        for i in 1..address.len() {
            address[i] += address[i - 1];
        }

        {
            let mut offset = vec![0; qgram_count];
            for (i, qgram) in QGrams::<u32>::new(q, text, alphabet).enumerate() {
                pos[address[qgram as usize] + offset[qgram as usize]] = i;
                offset[qgram as usize] += 1;
            }
        }

        QGramIndex { q: q, alphabet: alphabet, address: address, pos: pos }
    }

    pub fn matches(&self, qgram: u32) -> &[usize] {
        &self.pos[self.address[qgram as usize]..self.address[qgram as usize + 1]]
    }

    pub fn exact_matches(&self, pattern: &[u8]) -> Vec<ExactMatch> {
        let mut diagonals: collections::HashMap<usize, Vec<ExactMatch>> = collections::HashMap::new();
        for (i, qgram) in QGrams::<u32>::new(self.q, pattern, self.alphabet).enumerate() {
            for &p in self.matches(qgram) {
                let diagonal = p - i;
                if !diagonals.contains_key(&diagonal) {
                    // nothing yet, start new match
                    let mut intervals = vec![
                        ExactMatch {
                            pattern_start: i,
                            pattern_stop: i + self.q,
                            text_start: p,
                            text_stop: p + self.q
                        }
                    ];
                    diagonals.insert(diagonal, intervals);
                }
                else {
                    let mut intervals = diagonals.get_mut(&diagonal).unwrap();

                    let exact_match = intervals.iter().last().unwrap().pattern_stop == i;
                    if exact_match {
                        let mut interval = intervals.iter_mut().last().unwrap();
                        // extend exact match
                        interval.pattern_stop = i + self.q;
                        interval.text_stop = p + self.q;
                    }
                    else {
                        // mismatch or indel, start new match
                        intervals.push(
                            ExactMatch {
                                pattern_start: i,
                                pattern_stop: i + self.q,
                                text_start: p,
                                text_stop: p + self.q
                            }
                        )
                    }

                }
            }
        }
        diagonals.into_iter().flat_map(|(usize, intervals)| intervals.into_iter()).collect()
    }
}


pub struct ExactMatch {
    pub pattern_start: usize,
    pub pattern_stop: usize,
    pub text_start: usize,
    pub text_stop: usize,
}
