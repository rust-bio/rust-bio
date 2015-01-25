use std::collections::{Bitv, VecMap};
use std::iter::count;

use alphabets::Alphabet;


/// Construct suffix array for given text.
///
/// # Arguments
///
/// * `text` - the text
///
/// ```rust
/// use bio::data_structures::suffix_array;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array::get_suffix_array(text);
/// assert_eq!(pos, vec![
///     21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18,
///     2, 16, 0, 19, 4, 13, 10, 3, 12, 9
/// ]);
/// ```
pub fn get_suffix_array(text: &[u8]) -> Vec<usize> {
    let n = text.len();
    let transformed_text = transform_text(text);
    if transformed_text[n-1] != 0 {
        panic!("Expecting extra sentinel character being lexicographically smallest at the end of the text.");
    }

    let mut sais = SAIS::new(n);
    sais.construct(transformed_text.as_slice());

    sais.pos
}


fn transform_text(text: &[u8]) -> Vec<usize> {
    let ranks = Alphabet::new(text).get_ranks();

    text.iter()
        .map(|&c| *ranks.get(&(c as usize)).unwrap())
        .collect()
}


struct SAIS {
    pos: Vec<usize>,
    lms_pos: Vec<usize>,
    bucket_sizes: VecMap<usize>,
    bucket_start: Vec<usize>,
    bucket_end: Vec<usize>
}


impl SAIS {
    fn new(n: usize) -> Self {
        SAIS {
            pos: Vec::with_capacity(n),
            lms_pos: Vec::with_capacity(n),
            bucket_sizes: VecMap::new(),
            bucket_start: Vec::with_capacity(n),
            bucket_end: Vec::with_capacity(n)
        }
    }

    fn init_bucket_start(&mut self, text: &[usize]) {
        self.bucket_sizes.clear();
        self.bucket_start.clear();

        for &c in text.iter() {
            if !self.bucket_sizes.contains_key(&c) {
                self.bucket_sizes.insert(c, 0);
            }
            *(self.bucket_sizes.get_mut(&c).unwrap()) += 1;
        }

        let mut sum = 0;
        for &size in self.bucket_sizes.values() {
            self.bucket_start.push(sum);
            sum += size;
        }
    }

    /// initialize pointers to the last element of the buckets
    fn init_bucket_end(&mut self, text: &[usize]) {
        self.bucket_end.clear();
        for &r in self.bucket_start[1..].iter() {
            self.bucket_end.push(r - 1);
        }
        self.bucket_end.push(text.len() - 1);
    }

    fn lms_substring_eq(
        &self,
        text: &[usize],
        pos_types: &PosTypes,
        i: usize,
        j: usize
    ) -> bool {
        for k in count(0us, 1us) {
            let lmsi = pos_types.is_lms_pos(i + k);
            let lmsj = pos_types.is_lms_pos(j + k);
            if text[i + k] != text[j + k] {
                // different symbols
                return false;
            }
            if lmsi != lmsj {
                // different length
                return false;
            }
            if k > 0 && lmsi && lmsj {
                // same symbols and same length
                return true;
            }
        }
        false
    }

    fn construct(&mut self, text: &[usize]) {
        let pos_types = PosTypes::new(text);
        self.calc_lms_pos(text, &pos_types);
        self.calc_pos(text, &pos_types);
    }

    /// Step 1 of the SAIS algorithm.
    fn calc_lms_pos(&mut self, text: &[usize], pos_types: &PosTypes) {
        let n = text.len();

        // collect LMS positions
        self.lms_pos.clear();
        for r in (0..n) {
            if pos_types.is_lms_pos(r) {
                self.lms_pos.push(r);
            }
        }

        // sort LMS substrings by applying step 2 with unsorted LMS positions
        self.calc_pos(text, pos_types);

        // obtain pre-sorted LMS positions from sorted LMS substrings
        self.lms_pos.clear();
        for &p in self.pos.iter() {
            if pos_types.is_lms_pos(p) {
                self.lms_pos.push(p);
            }
        }
        let lms_substring_count = self.lms_pos.len() - 1;

        // if less than 2 LMS substrings are present, no further sorting is needed
        if lms_substring_count > 1 {
            // sort LMS suffixes by recursively building SA on reduced text
            let mut reduced_text = Vec::with_capacity(self.lms_pos.len());
            let mut label = 0;
            reduced_text.push(label);
            let mut prev = self.lms_pos[0];
            for &p in self.lms_pos[1..].iter() {
                // choose same label if substrings are equal
                if !self.lms_substring_eq(
                    text, pos_types, prev, p
                ) {
                    label += 1;
                }
                reduced_text.push(label);
                prev = p;
            }

            // if we have less labels than substrings, we have to sort by recursion
            // because two or more substrings are equal
            if label + 1 < lms_substring_count {
                // backup lms_pos
                let lms_pos = self.lms_pos.clone();
                // recurse SA construction for reduced text
                self.construct(reduced_text.as_slice());
                // obtain sorted lms suffixes
                self.lms_pos.clear();
                for &p in self.pos.iter() {
                    self.lms_pos.push(lms_pos[p]);
                }
            }
        }
    }

    /// Step 2 of the SAIS algorithm.
    fn calc_pos(&mut self, text: &[usize], pos_types: &PosTypes) {
        let n = text.len();
        self.pos.clear();

        self.init_bucket_start(text);
        self.init_bucket_end(text);

        // init all positions as unknown (n-1 is max position)
        for _ in text.iter() {
            self.pos.push(n);
        }

        // insert LMS positions to the end of their buckets
        for &p in self.lms_pos.iter().rev() {
            let c = text[p];
            self.pos[self.bucket_end[c]] = p;
            self.bucket_end[c] -= 1;
        }

        // reset bucket ends
        self.init_bucket_end(text);

        // insert L-positions into buckets
        for r in (0..n) {
            let p = self.pos[r];
            // ignore undefined positions and the zero since it has no predecessor
            if p == n || p == 0 {
                continue;
            }
            let pred = p - 1;
            if pos_types.is_l_pos(pred) {
                let c = text[pred];
                self.pos[self.bucket_start[c]] = pred;
                self.bucket_start[c] += 1;
            }
        }

        // insert S-positions into buckets
        for r in (0..n).rev() {
            let p = self.pos[r];
            if p == 0 {
                continue;
            }
            let pred = p - 1;
            if pos_types.is_s_pos(pred) {
                let c = text[pred];
                self.pos[self.bucket_end[c]] = pred;
                self.bucket_end[c] -= 1;
            }
        }
    }
}


struct PosTypes {
    pos_types: Bitv,
}


impl PosTypes {
    /// Calculate the text position type.
    /// S-type marks suffixes being lexicographically smaller than their successor,
    /// L-type marks those being larger.
    /// This function fills a Bitv, with 1-bits denoting S-type
    /// and 0-bits denoting L-type.
    ///
    /// # Arguments
    ///
    /// * `text` - the text, ending with a sentinel.
    fn new(text: &[usize]) -> Self {
        let n = text.len();
        let mut pos_types = Bitv::from_elem(n, false);
        pos_types.set(n-1, true);

        for p in (0..n-1).rev() {
            if text[p] == text[p + 1] {
                // if the characters are equal, the next position determines
                // the lexicographical order
                let v = pos_types.get(p + 1).unwrap();
                pos_types.set(p, v);
            }
            else {
                pos_types.set(p, text[p] < text[p + 1]);
            }
        }
        
        PosTypes { pos_types: pos_types }
    }

    fn is_s_pos(&self, p: usize) -> bool {
        self.pos_types.get(p).unwrap()
    }

    fn is_l_pos(&self, p: usize) -> bool {
        !self.pos_types.get(p).unwrap()
    }

    fn is_lms_pos(&self, p: usize) -> bool {
        !(p == 0) && self.is_s_pos(p) && self.is_l_pos(p-1)
    }
}


#[test]
fn test_pos_types() {
    let text = transform_text(b"GCCTTAACATTATTACGCCTA$");
    let n = text.len();

    let pos_types = PosTypes::new(text.as_slice());
    let mut test = Bitv::from_bytes(&[0b01100110, 0b10010011,  0b01100100]);
    test.truncate(n);
    assert_eq!(pos_types.pos_types, test);
    let lms_pos: Vec<usize> = (0..n).filter(|&p| pos_types.is_lms_pos(p)).collect();
    assert_eq!(lms_pos, vec![1, 5, 8, 11, 14, 17, 21]);
}


#[test]
fn test_buckets() {
    let text = transform_text(b"GCCTTAACATTATTACGCCTA$");
    let n = text.len();

    let mut sais = SAIS::new(n);
    sais.init_bucket_start(text.as_slice());
    assert_eq!(sais.bucket_start, vec![0, 1, 7, 13, 15]);
    sais.init_bucket_end(text.as_slice());
    assert_eq!(sais.bucket_end, vec![0, 6, 12, 14, 21]);
}


#[test]
fn test_pos() {
    let text = transform_text(b"GCCTTAACATTATTACGCCTA$");
    let n = text.len();

    let mut sais = SAIS::new(n);
    let pos_types = PosTypes::new(text.as_slice());
    sais.lms_pos = vec![21, 5, 14, 8, 11, 17, 1];
    sais.calc_pos(text.as_slice(), &pos_types);
    assert_eq!(sais.pos, vec![
        21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18,
        2, 16, 0, 19, 4, 13, 10, 3, 12, 9
    ]);
}


#[test]
fn test_lms_pos() {
    let text = transform_text(b"GCCTTAACATTATTACGCCTA$");
    let n = text.len();

    let mut sais = SAIS::new(n);
    let pos_types = PosTypes::new(text.as_slice());
    sais.calc_lms_pos(text.as_slice(), &pos_types);
}
