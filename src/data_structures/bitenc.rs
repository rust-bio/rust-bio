// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A fixed-width bit encoding implementation. This allows to store a sequence of values over
//! a reduced alphabet by packing them bit-encoded into a sequence of bytes.
//!
//! # Example
//!
//! ```
//! use bio::data_structures::bitenc::BitEnc;
//! let mut bitenc = BitEnc::new(2);
//! bitenc.push(0);
//! bitenc.push(2);
//! bitenc.push(1);
//! let values: Vec<u8> = bitenc.iter().collect();
//! assert_eq!(values, [0, 2, 1]);
//! ```

/// A sequence of bitencoded values.
#[derive(Serialize, Deserialize, PartialEq, Hash)]
pub struct BitEnc {
    storage: Vec<u32>,
    width: usize,
    mask: u32,
    len: usize,
    bits: usize,
}

impl Eq for BitEnc {}

fn mask(width: usize) -> u32 {
    (1 << width) - 1
}

impl BitEnc {
    /// Create a new instance with a given encoding width (e.g. width=2 for using two bits per value).
    pub fn new(width: usize) -> Self {
        assert!(width <= 8, "Only encoding widths up to 8 supported");
        BitEnc {
            storage: Vec::new(),
            width,
            mask: mask(width),
            len: 0,
            bits: 32 - 32 % width,
        }
    }

    /// Create a new instance with a given capacity and encoding width (e.g. width=2 for using two bits per value).
    pub fn with_capacity(width: usize, n: usize) -> Self {
        assert!(width <= 8, "Only encoding widths up to 8 supported");
        BitEnc {
            storage: Vec::with_capacity(n * width / 32),
            width,
            mask: mask(width),
            len: 0,
            bits: 32 - 32 % width,
        }
    }

    /// Append a value.
    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.addr(self.len);
        if bit == 0 {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    /// Append `n` times the given value.
    pub fn push_values(&mut self, mut n: usize, value: u8) {
        {
            // fill the last block
            let (block, mut bit) = self.addr(self.len);
            if bit > 0 {
                // TODO use step_by once it has been stabilized: for bit in (bit..32).step_by(self.width) {
                while bit <= 32 {
                    self.set_by_addr(block, bit, value);
                    n -= 1;
                    bit += self.width
                }
            }
        }

        // pack the value into a block
        let mut value_block = 0;
        {
            let mut v = u32::from(value);
            for _ in 0..32 / self.width {
                value_block |= v;
                v <<= self.width;
            }
        }

        // push as many value blocks as needed
        let i = self.len + n;
        let (block, bit) = self.addr(i);
        for _ in self.storage.len()..block {
            self.storage.push(value_block);
        }

        if bit > 0 {
            // add the remaining values to a final block
            self.storage.push(value_block >> (32 - bit));
        }

        self.len = i;
    }

    /// Set the value as position `i`.
    pub fn set(&mut self, i: usize, value: u8) {
        let (block, bit) = self.addr(i);
        self.set_by_addr(block, bit, value);
    }

    /// Get the value at position `i`.
    pub fn get(&self, i: usize) -> Option<u8> {
        if i >= self.len {
            None
        } else {
            let (block, bit) = self.addr(i);
            Some(self.get_by_addr(block, bit))
        }
    }

    /// Iterate over stored values (values will be unpacked into bytes).
    pub fn iter(&self) -> BitEncIter<'_> {
        BitEncIter { bitenc: self, i: 0 }
    }

    /// Clear the sequence.
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
        self.storage[block] |= (u32::from(value) & self.mask) << bit;
    }

    fn addr(&self, i: usize) -> (usize, usize) {
        let k = i * self.width;
        (k / self.bits, k % self.bits)
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// Iterator over values of a bitencoded sequence (values will be unpacked into bytes).
pub struct BitEncIter<'a> {
    bitenc: &'a BitEnc,
    i: usize,
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
        assert_eq!(bitenc.storage, [0, 0]);
    }

    #[test]
    fn test_issue29() {
        for w in 2..9 {
            let mut vec = BitEnc::with_capacity(w, 1000);
            for i in 0..1000 {
                println!("Push {}", i);
                vec.push(1);
            }
        }
    }
}
