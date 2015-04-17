// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Myers bit-parallel approximate pattern matching algorithm.
//! Finds all matches up to a given edit distance. The pattern has to fit into a bitvector,
//! and is here limited to 64 symbols.
//! Complexity: O(n)


use std::iter;


pub struct Myers {
    peq: [u64; 256],
    bound: u64,
    m: u8,
}


impl Myers {
    pub fn new(pattern: &[u8]) -> Self {
        assert!(pattern.len() <= 64 && pattern.len() > 0);

        let mut peq = [0; 256];
        for (i, &a) in pattern.iter().enumerate() {
            peq[a as usize] |= 1 << i;
        }

        Myers {
            peq: peq,
            bound: 1 << (pattern.len() -1),
            m: pattern.len() as u8,
        }
    }

    fn step(&self, state: &mut State, a: u8) {
        let eq = self.peq[a as usize];
        let xv = eq | state.mv;
        let xh = (((eq & state.pv) + state.pv) ^ state.pv) | eq;

        let mut ph = state.mv | !( xh | state.pv);
        let mut mh = state.pv & xh;

        if ph & self.bound > 0 {
            state.dist += 1;
        }
        else if mh & self.bound > 0 {
            state.dist -= 1;
        }

        ph <<= 1;
        mh <<= 1;
        state.pv = mh | !(xv | ph);
        state.mv = ph & xv;
    }

    pub fn distance(&self, text: &[u8]) -> u8 {
        let mut state = State::new(self.m);
        for &a in text {
            self.step(&mut state, a);
        }
        state.dist
    }

    // Find all occurences of pattern in the given text.
    pub fn find_all_end<'a, I: Iterator<Item=&'a u8>>(&'a self, text: I, max_dist: u8) -> Matches<I> {
        Matches { myers: self, state: State::new(self.m), text: text.enumerate(), max_dist: max_dist }
    }
}


struct State {
    pv: u64,
    mv: u64,
    dist: u8,
}


impl State {
    pub fn new(m: u8) -> Self {
        State { pv: (1 << m) - 1, mv: 0, dist: m }
    }
}


pub struct Matches<'a, I: Iterator<Item=&'a u8>> {
    myers: &'a Myers,
    state: State,
    text: iter::Enumerate<I>,
    max_dist: u8,
}


impl<'a, I: Iterator<Item=&'a u8>> Iterator for Matches<'a, I> {
    type Item = (usize, u8);

    fn next(&mut self) -> Option<(usize, u8)> {
        for (i, &a) in self.text.by_ref() {
            self.myers.step(&mut self.state, a);
            if self.state.dist <= self.max_dist {
                return Some((i, self.state.dist));
            }
        }
        None
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    #[test]
    fn test_find_all_end() {
        let text = b"ACCGTGGATGAGCGCCATAG";
        let pattern = b"TGAGCGT";
        let myers = Myers::new(pattern);
        let occ = myers.find_all_end(text.iter(), 1).collect_vec();
        assert_eq!(occ, [(13, 1), (14, 1)]);
    }
}
