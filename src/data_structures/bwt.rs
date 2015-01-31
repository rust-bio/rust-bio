//! The Burrows-Wheeler-Transform and related data structures.

use std::iter::repeat;
use std::iter::AdditiveIterator;

use data_structures::suffix_array::SuffixArray;
use utils::prescan;
use alphabets::Alphabet;


pub type BWT = Vec<u8>;
type Less = Vec<usize>;
type BWTFind = Vec<usize>;


/// Calculate Burrows-Wheeler-Transform of the given text of length n.
/// Complexity: O(n).
///
/// # Arguments
///
/// * `text` - the text ended by sentinel symbol (being lexicographically smallest)
/// * `pos` - the suffix array for the text
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::get_suffix_array;
/// use bio::data_structures::bwt::get_bwt;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = get_suffix_array(text);
/// let bwt = get_bwt(text, &pos);
/// assert_eq!(bwt, b"ATTATTCAGGACCC$CTTTCAA");
/// ```
pub fn get_bwt(text: &[u8], pos: &SuffixArray) -> BWT {
    assert!(text.len() == pos.len());
    let n = text.len();
    let mut bwt: BWT = repeat(0).take(n).collect();
    for r in 0..n {
        let p = pos[r];
        bwt[r] = if p > 0 {text[p-1]} else {text[n-1]};
    }

    bwt
}


/// Calculate the inverse of a BWT of length n, which is the original text.
/// Complexity: O(n).
///
/// # Arguments
///
/// * `bwt` - the BWT
pub fn get_inverse(bwt: &BWT) -> Vec<u8> {
    let alphabet = Alphabet::new(bwt.as_slice());
    let n = bwt.len();
    let bwtfind = get_bwtfind(bwt, &alphabet);
    let mut inverse = Vec::with_capacity(n);

    let mut r = bwtfind[0];
    for _ in 0..n {
        r = bwtfind[r];
        inverse.push(bwt[r]);
    }

    inverse
}


pub struct FMIndex<'a> {
    bwt: &'a BWT,
    less: Less,
    occ: Occ
}


impl<'a> FMIndex<'a> {
    pub fn new(bwt: &'a BWT, k: usize, alphabet: &Alphabet) -> Self {
        FMIndex { bwt: bwt, less: get_less(bwt, alphabet), occ: Occ::new(bwt, k, alphabet)}
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
    /// use bio::data_structures::bwt::{get_bwt, FMIndex};
    /// use bio::data_structures::suffix_array::get_suffix_array;
    /// use bio::alphabets::get_dna_alphabet;
    /// let text = b"GCCTTAACATTATTACGCCTA$";
    /// let alphabet = get_dna_alphabet();
    /// let pos = get_suffix_array(text);
    /// let bwt = get_bwt(text, &pos);
    /// let fm = FMIndex::new(&bwt, 3, &alphabet);
    /// let pattern = b"TTA";
    /// let sai = fm.backward_search(pattern);
    /// assert_eq!(sai, (19, 21));
    /// ```
    pub fn backward_search(&self, pattern: &[u8]) -> (usize, usize) {
        let (mut l, mut r) = (0, self.bwt.len() - 1);
        for &a in pattern.iter().rev() {
            let less = self.less[a as usize];
            l = less + if l > 0 { self.occ.get_occ(self.bwt, l - 1, a) } else { 0 };
            r = less + self.occ.get_occ(self.bwt, r, a) - 1;
        }

        (l, r)
    }
}


pub struct Occ {
    occ: Vec<Vec<usize>>,
    k: usize
}


impl Occ {
    /// Calculate occ array with sampling from BWT of length n.
    /// Time complexity: O(n).
    /// Space complexity: O(n / k * A) with A being the alphabet size.
    /// Alphabet size is determined on the fly from the BWT.
    /// For large texts, it is therefore advisable to transform
    /// the text before calculating the BWT (see alphabets::rank_transform).
    ///
    /// # Arguments
    ///
    /// * `bwt` - the BWT
    /// * `k` - the sampling rate: every k-th entry will be stored
    pub fn new(bwt: &BWT, k: usize, alphabet: &Alphabet) -> Self {
        let n = bwt.len();
        let m = alphabet.max_symbol().expect("Expecting non-empty alphabet.") as usize + 1;
        let mut occ = Vec::with_capacity(n / k);
        let mut curr_occ: Vec<usize> = repeat(0).take(m).collect();
        for (i, &c) in bwt.iter().enumerate() {
            curr_occ[c as usize] += 1;
            if i % k == 0 {
                occ.push(curr_occ.clone());
            }
        }

        Occ { occ: occ, k: k }
    }

    pub fn get_occ(&self, bwt: &BWT, r: usize, a: u8) -> usize {
        let i = r / self.k;
        self.occ[i][a as usize] +
        bwt[(i * self.k) + 1 .. r + 1].iter().map(|&c| (c == a) as usize).sum()
    }
}


fn get_less(bwt: &BWT, alphabet: &Alphabet) -> Less {
    let m = alphabet.max_symbol().expect("Expecting non-empty alphabet.") as usize + 1;
    let mut less: Less = repeat(0)
        .take(m).collect();
    for &c in bwt.iter() {
        less[c as usize] += 1;
    }
    // calculate +-prescan
    prescan(less.as_mut_slice(), 0, |a, b| a + b);

    less
}


/// Calculate the bwtfind array needed for inverting the BWT.
fn get_bwtfind(bwt: &BWT, alphabet: &Alphabet) -> BWTFind {
    let n = bwt.len();
    let mut less = get_less(bwt, alphabet);

    let mut bwtfind: BWTFind = repeat(0).take(n).collect();
    for (r, &c) in bwt.iter().enumerate() {
        bwtfind[less[c as usize]] = r;
        less[c as usize] += 1;
    }

    bwtfind
}


#[cfg(test)]
mod tests {
    use super::{get_bwtfind, get_bwt, get_inverse, Occ};
    use data_structures::suffix_array::get_suffix_array;
    use alphabets::Alphabet;

    #[test]
    fn test_get_bwtfind() {
        let text = b"cabca$";
        let alphabet = Alphabet::new(b"abc$");
        let pos = get_suffix_array(text);
        let bwt = get_bwt(text, &pos);
        let bwtfind = get_bwtfind(&bwt, &alphabet);
        assert_eq!(bwtfind, vec![5, 0, 3, 4, 1, 2]);
    }

    #[test]
    fn test_get_inverse() {
        let text = b"cabca$";
        let pos = get_suffix_array(text);
        let bwt = get_bwt(text, &pos);
        let inverse = get_inverse(&bwt);
        assert_eq!(inverse, text);
    }

    #[test]
    fn test_occ () {
        let bwt = vec![1u8, 3u8, 3u8, 1u8, 2u8, 0u8];
        let alphabet = Alphabet::new([0u8, 1u8, 2u8, 3u8].as_slice());
        let occ = Occ::new(&bwt, 3, &alphabet);
        assert_eq!(occ.occ, [
            [0, 1, 0, 0],
            [0, 2, 0, 2]
        ]);
        assert_eq!(occ.get_occ(&bwt, 4, 2u8), 1);
        assert_eq!(occ.get_occ(&bwt, 4, 3u8), 2);
    }
}
