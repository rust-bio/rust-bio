// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::iter::DoubleEndedIterator;

use data_structures::bwt::{Occ, Less, less, BWT};
use alphabets::Alphabet;


pub struct FMIndex<'a> {
    bwt: &'a BWT,
    less: Less,
    occ: Occ
}


impl<'a> FMIndex<'a> {
    pub fn new(bwt: &'a BWT, k: usize, alphabet: &Alphabet) -> Self {
        FMIndex { bwt: bwt, less: less(bwt, alphabet), occ: Occ::new(bwt, k, alphabet)}
    }

    /// Perform backward search, yielding suffix array
    /// interval denoting positions where the given pattern occurs.
    ///
    /// # Arguments
    ///
    /// * `pattern` - the pattern to search
    ///
    /// # Example
    ///
    /// ```
    /// use bio::data_structures::bwt::bwt;
    /// use bio::data_structures::fmindex::FMIndex;
    /// use bio::data_structures::suffix_array::suffix_array;
    /// use bio::alphabets::dna;
    /// let text = b"GCCTTAACATTATTACGCCTA$";
    /// let alphabet = dna::alphabet();
    /// let pos = suffix_array(text);
    /// let bwt = bwt(text, &pos);
    /// let fm = FMIndex::new(&bwt, 3, &alphabet);
    /// let pattern = b"TTA";
    /// let sai = fm.backward_search(pattern.iter());
    /// assert_eq!(sai, (19, 21));
    /// ```
    pub fn backward_search<'b, P: Iterator<Item=&'b u8> + DoubleEndedIterator>(&self, pattern: P) -> (usize, usize) {
        let (mut l, mut r) = (0, self.bwt.len() - 1);
        for &a in pattern.rev() {
            let less = self.less[a as usize];
            l = less + if l > 0 { self.occ.get(self.bwt, l - 1, a) } else { 0 };
            r = less + self.occ.get(self.bwt, r, a) - 1;
        }

        (l, r)
    }
}
