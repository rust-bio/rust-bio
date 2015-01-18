use core::iter::Enumerate;
use core::slice;

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
    m: uint,
    masks: [u64; 255],
    accept: u64
}


impl ShiftAnd {
    /// Create new ShiftAnd instance.
    pub fn new(pattern: &[u8]) -> ShiftAnd {
        let mut masks = [0; 255];

        let mut bit = 1;
        for c in pattern.iter() {
            masks[*c as uint] |= bit;
            bit *= 2;
        }

        ShiftAnd { m: pattern.len(), masks: masks, accept: bit / 2 }

    }

    /// Find all occurences of pattern in the given text.
    pub fn find_all(self, text: &[u8]) -> FindAll {
        FindAll { shiftand: self, active: 0, text: text.iter().enumerate() }
    }
}


pub struct FindAll<'a> {
    shiftand: ShiftAnd,
    active: u64,
    text: Enumerate<slice::Iter<'a, u8>>
}


impl<'a> Iterator for FindAll<'a> {
    type Item = uint;

    fn next(&mut self) -> Option<uint> {
        for (i, c) in self.text {
            self.active = ((self.active << 1) | 1) & self.shiftand.masks[*c as uint];
            if self.active & self.shiftand.accept > 0 {
                return Some(i - self.shiftand.m + 1);
            }
        }

        None
    }
}
