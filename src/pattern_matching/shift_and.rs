use std::iter::Enumerate;
use std::slice;


/// ShiftAnd algorithm for pattern matching.
///
/// # Example
///
/// ```rust
/// use bio::pattern_matching::shift_and;
/// let pattern = b"AAAA";
/// let text = b"ACGGCTAGAAAAGGCTAG";
/// let shiftand = shift_and::ShiftAnd::new(pattern);
/// let occ = shiftand.find_all(text).next().unwrap();
/// assert_eq!(occ, 8);
/// ```
#[derive(Copy)]
pub struct ShiftAnd {
    m: u8,
    masks: [u64; 256],
    accept: u64
}


impl ShiftAnd {
    /// Create new ShiftAnd instance.
    pub fn new(pattern: &[u8]) -> ShiftAnd {
        if pattern.len() > 64 {
            panic!("Expecting pattern of size at most 64");
        }
        let mut masks = [0; 256];

        let mut bit = 1;
        for c in pattern.iter() {
            masks[*c as usize] |= bit;
            bit *= 2;
        }

        ShiftAnd { m: pattern.len() as u8, masks: masks, accept: bit / 2 }

    }

    /// Find all occurences of pattern in the given text.
    pub fn find_all<'a>(&'a self, text: &'a [u8]) -> FindAll {
        FindAll { shiftand: self, active: 0, text: text.iter().enumerate() }
    }
}


pub struct FindAll<'a> {
    shiftand: &'a ShiftAnd,
    active: u64,
    text: Enumerate<slice::Iter<'a, u8>>
}


impl<'a> Iterator for FindAll<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        for (i, &c) in self.text {
            self.active = ((self.active << 1) | 1) & self.shiftand.masks[c as usize];
            if self.active & self.shiftand.accept > 0 {
                return Some(i - self.shiftand.m as usize + 1);
            }
        }

        None
    }
}


#[cfg(test)]
mod tests {
    use test::Bencher;
    use super::ShiftAnd;

    #[bench]
    fn bench_find_all(b: &mut Bencher) {
        let pattern = b"AAAA";
        let text = b"ACGGCTAGAAAAGGCTAGGAGTAGGATTCTGCATGCACGACTCGAGCACTAGCACT";
        let shiftand = ShiftAnd::new(pattern);
        b.iter(|| {
            shiftand.find_all(text).collect::<Vec<usize>>()
        });
    }
}
