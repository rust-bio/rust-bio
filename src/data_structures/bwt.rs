
use std::iter::repeat;

use data_structures::suffix_array::SuffixArray;

pub type BWT = Vec<u8>;


/// Calculate Burrows-Wheeler-Transform of the given text of length n.
/// Complexity: O(n)
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


//pub fn backward_search(
