use std::collections::{Bitv, VecMap};
use std::collections::bitv::Iter;


struct SAIS {
    pos: Vec<usize>,
    pos_types: Bitv,
    lms_pos: Vec<usize>,
    bucket_sizes: VecMap<usize>,
    bucket_start: Vec<usize>,
    bucket_end: Vec<usize>
}


impl SAIS {
    fn new(n: usize) -> Self {
        SAIS {
            pos: Vec::with_capacity(n),
            pos_types: Bitv::with_capacity(n),
            lms_pos: Vec::with_capacity(n),
            bucket_sizes: VecMap::new(),
            bucket_start: Vec::with_capacity(n),
            bucket_end: Vec::with_capacity(n)
        }
    }

    fn init_pos_types(&mut self, text: &[usize]) {
        let n = text.len();
        self.pos_types.clear();

        for _ in 0..n-1 {
            self.pos_types.push(false);
        }
        self.pos_types.push(true);

        for p in (0..n-1).rev() {
            if text[p] == text[p + 1] {
                // if the characters are equal, the next position determines
                // the lexicographical order
                let v = self.pos_types.get(p + 1).unwrap();
                self.pos_types.set(p, v);
            }
            else {
                self.pos_types.set(p, text[p] < text[p + 1]);
            }
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

    fn init_bucket_end(&mut self, text: &[usize]) {
        self.bucket_end.clear();
        for &r in self.bucket_start[1..].iter() {
            self.bucket_end.push(r);
        }
        self.bucket_end.push(text.len());
    }

    fn is_s_pos(&self, r: usize) -> bool {
        self.pos_types.get(r).unwrap()
    }

    fn is_l_pos(&self, r: usize) -> bool {
        !self.pos_types.get(r).unwrap()
    }

    fn is_lms_pos(&self, r: usize) -> bool {
        self.is_s_pos(r) && self.is_l_pos(r)
    }

    /// Step 1 of the SAIS algorithm.
    fn calc_lms_pos(&mut self, text: &[usize]) {
        let n = text.len();

        // collect LMS positions
        self.lms_pos.clear();
        for r in (0..n) {
            if self.is_lms_pos(r) {
                self.lms_pos.push(r);
            }
        }

        // sort LMS substrings
        self.calc_pos(text);

        
    }

    /// Step 2 of the SAIS algorithm.
    fn calc_pos(&mut self, text: &[usize]) {
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
            if p == n {
                continue;
            }
            let prev = p - 1;
            if self.is_l_pos(prev) {
                let c = text[prev];
                self.pos[self.bucket_start[c]] = prev;
                self.bucket_start[c] += 1;
            }
        }

        // insert S-positions into buckets
        for r in (0..n).rev() {
            let prev = self.pos[r] - 1;
            if self.is_s_pos(prev) {
                let c = text[prev];
                self.pos[self.bucket_end[c]] = prev;
                self.bucket_end[c] -= 1;
            }
        }
    }
}


/// DEPRECATED: Calculate the text position type.
/// S-type marks suffixes being lexicographically smaller than their successor,
/// L-type marks those being larger.
/// This function returns a Bitv, with 1-bits denoting S-type
/// and 0-bits denoting L-type.
///
/// # Arguments
///
/// * `text` - the text, ending with a sentinel '$'.
///
/// # Example
///
/// ```rust
/// use std::collections::Bitv;
/// use bio::data_structures::suffix_array::get_pos_type;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos_type = get_pos_type(text);
/// let mut test = Bitv::from_bytes(&[0b01100110, 0b10010011,  0b01100100]);
/// test.truncate(text.len());
/// assert_eq!(pos_type, test);
/// ```
pub fn get_pos_types(text: &[usize]) -> Bitv {
    let n = text.len();
    let mut pos_types = Bitv::from_elem(n, false);
    pos_types.set(n - 1, true);
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

    pos_types
}   
