
use std::iter::{repeat, Enumerate};
use std::slice;
use std::collections::VecMap;

use alphabets::Alphabet;


type LPS = Vec<usize>;


pub struct KMP {
    m: usize,
    table: Vec<VecMap<usize>>
}


impl KMP {
    pub fn new(pattern: &[u8], alphabet: Alphabet) -> Self {
        //assert!(alphabet.is_word(pattern));
        let k = alphabet.max_symbol()
            .expect("Expecting non-empty alphabet.") as usize + 1;
        let m = pattern.len();

        let mut init = VecMap::with_capacity(k);
        for c in alphabet.symbols.iter() {
            init.insert(c, 0);
        }
        *init.get_mut(&(pattern[0] as usize)).unwrap() = 1;

        let lps = get_lps(pattern);

        let mut table = Vec::with_capacity(m + 1);
        table.push(init);
        for q in 1..m+1 {
            let mut dq = VecMap::with_capacity(k);
            for c in alphabet.symbols.iter() {
                dq.insert(c, *table[lps[q - 1]].get(&c).unwrap());
            }
            if q < m {
                *dq.get_mut(&(pattern[q] as usize)).unwrap() = q + 1;
            }
            table.push(dq);
        }

        KMP { table: table, m: m }
    }

    fn delta(&self, q: usize, a: u8) -> usize {
        *self.table[q].get(&(a as usize)).expect("Missing symbol in alphabet (is the text a word of the given alphabet?)")
    }

    pub fn find_all<'a>(&'a self, text: &'a [u8]) -> FindAll {
        FindAll { kmp: self, q: 0, text: text.iter().enumerate() }
    }
}


fn get_lps(pattern: &[u8]) -> LPS {
    let (m, mut q) = (pattern.len(), 0us);
    let mut lps: LPS = repeat(0).take(m).collect();
    for i in 1..m {
        while q > 0 && pattern[q] != pattern[i] {
            q = lps[q];
        }
        if pattern[q] == pattern[i] {
            q += 1;
        }
        lps[i] = q;
    }

    lps
}


pub struct FindAll<'a> {
    kmp: &'a KMP,
    q: usize,
    text: Enumerate<slice::Iter<'a, u8>>
}


impl<'a> Iterator for FindAll<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        for (i, &c) in self.text {
            self.q = self.kmp.delta(self.q, c);
            if self.q == self.kmp.m {
                return Some(i - self.kmp.m + 1);
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::{get_lps, KMP};
    use alphabets::Alphabet;

    #[test]
    fn test_get_lps() {
        let pattern = b"ababaca";
        let lps = get_lps(pattern);
        assert_eq!(lps, [0, 0, 1, 2, 3, 0, 1]);
    }

    #[test]
    fn test_delta() {
        let pattern = b"abbab";
        let alphabet = Alphabet::new(pattern);
        let kmp = KMP::new(pattern, alphabet);
        assert_eq!(kmp.delta(0, b'a'), 1);
        assert_eq!(kmp.delta(0, b'b'), 0);
        assert_eq!(kmp.delta(1, b'a'), 1);
        assert_eq!(kmp.delta(1, b'b'), 2);
        assert_eq!(kmp.delta(2, b'a'), 1);
        assert_eq!(kmp.delta(2, b'b'), 3);
        assert_eq!(kmp.delta(3, b'a'), 4);
        assert_eq!(kmp.delta(3, b'b'), 0);
        assert_eq!(kmp.delta(4, b'a'), 1);
        assert_eq!(kmp.delta(4, b'b'), 5);
        assert_eq!(kmp.delta(5, b'a'), 1);
        assert_eq!(kmp.delta(5, b'b'), 3);
    }

    #[test]
    fn test_find_all() {
        let text = b"aaaaabbabbbbbbbabbab";
        let pattern = b"abbab";
        let alphabet = Alphabet::new(pattern);
        let kmp = KMP::new(pattern, alphabet);
        let occ: Vec<usize> = kmp.find_all(text).collect();
        assert_eq!(occ, [4, 15]);
    }
}
