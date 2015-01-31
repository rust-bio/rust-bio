use std::iter::range_step;

/// A fixed-width bit encoding implementation.
///
/// # Example
///
/// ```
/// use bio::bitencoding::BitEnc;
/// let mut bitenc = BitEnc::new(2);
/// bitenc.push(0);
/// bitenc.push(2);
/// bitenc.push(1);
/// let values: Vec<u8> = bitenc.iter().collect();
/// assert_eq!(values, [0, 2, 1]);
/// ```
pub struct BitEnc {
    pub storage: Vec<u32>,
    width: usize,
    mask: u32,
    len: usize
}


fn get_mask(width: usize) -> u32 {
    (1 << width) - 1
}


impl BitEnc {
    pub fn new(width: usize) -> Self {
        assert!(width <= 8, "Only encoding widths up to 8 supported");
        BitEnc { storage: Vec::new(), width: width, mask: get_mask(width), len: 0 }
    }

    pub fn with_capacity(width: usize, n: usize) -> Self {
        assert!(width <= 8, "Only encoding widths up to 8 supported");
        BitEnc { storage: Vec::with_capacity(n * width / 32), width: width, mask: get_mask(width), len: 0 }
    }

    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.get_addr(self.len);
        if bit == 0 {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    pub fn push_values(&mut self, mut n: usize, value: u8) {
        {
            // fill the last block
            let (block, bit) = self.get_addr(self.len);
            if bit > 0 {
                for bit in range_step(bit, 32, self.width) {
                    self.set_by_addr(block, bit, value);
                    n -= 1;
                }
            }
        }

        // pack the value into a block
        let mut value_block = 0;
        {
            let mut v = value as u32;
            for _ in 0..32/self.width {
                value_block |= v;
                v <<= self.width;
            }
        }

        // push as many value blocks as needed
        let i = self.len + n;
        let (block, bit) = self.get_addr(i);
        for _ in self.storage.len()..block {
            self.storage.push(value_block);
        }

        if bit > 0 {
            // add the remaining values to a final block
            self.storage.push(value_block >> 32 - bit);
        }

        self.len = i;
    }

    pub fn set(&mut self, i: usize, value: u8) {
        let (block, bit) = self.get_addr(i);
        self.set_by_addr(block, bit, value);
    }

    pub fn get(&self, i: usize) -> Option<u8> {
        if i >= self.len {
            None
        }
        else {
            let (block, bit) = self.get_addr(i);
            Some(self.get_by_addr(block, bit))
        }
    }

    pub fn iter(&self) -> BitEncIter {
        BitEncIter { bitenc: self, i: 0 }
    }

    pub fn clear(&mut self) {
        self.storage.clear();
        self.len = 0;
    }

    fn get_by_addr(&self, block: usize, bit: usize) -> u8 {
        ((self.storage[block] >> bit) & self.mask) as u8
    }

    fn set_by_addr(&mut self, block: usize, bit: usize, value: u8) {
        let mask = self.mask << bit;
        self.storage[block] |= mask;
        self.storage[block] ^= mask;
        self.storage[block] |= (value as u32 & self.mask) << bit;
    }

    fn get_addr(&self, i: usize) -> (usize, usize) {
        let k = i * self.width;
        (k / 32, k % 32)
    }

    pub fn len(&self) -> usize {
        self.len
    }
}


pub struct BitEncIter<'a> {
    bitenc: &'a BitEnc,
    i: usize
}


impl<'a> Iterator for BitEncIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        let value = self.bitenc.get(self.i);
        self.i += 1;
        value
    }
}


#[cfg(test)]
mod tests {
    use super::BitEnc;

    #[test]
    fn test_bitenc() {
        let mut bitenc = BitEnc::new(2);
        bitenc.push(0);
        bitenc.push(2);
        bitenc.push(1);
        let mut values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [0, 2, 1]);
        bitenc.set(1, 3);
        values = bitenc.iter().collect();
        assert_eq!(values, [0, 3, 1]);
    }

    #[test]
    fn test_push_values() {
        let mut bitenc = BitEnc::new(2);
        bitenc.push_values(32, 0);
        assert_eq!(bitenc.storage, [0,0]);
    }
}
