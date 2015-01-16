

/// Calculate the masks for the given pattern.
fn get_masks(pattern: &[u8]) -> ([u64; 255], u64) {
    let mut masks = [0; 255];

    let mut bit = 1;
    for c in pattern.iter() {
        masks[*c as uint] |= bit;
        bit *= 2;
    }

    (masks, bit / 2)
}


/// Find the first occurrence of a pattern in a text.
///
/// # Arguments
///
/// * `pattern` - the pattern to search
/// * `text` - the text to search in
///
/// # Example
///
/// ```rust
/// use bio::patternmatching::shift_and;
/// let pattern = b"AAAA";
/// let text = b"ACGGCTAGAAAAGGCTAG";
/// let occ = shift_and::find_first(pattern, text).unwrap();
/// assert_eq!(occ, 8);
/// ```
pub fn find_first(pattern: &[u8], text: &[u8]) -> Option<uint> {
    let m = pattern.len();
    let (masks, accept) = get_masks(pattern);

    // simulate the NFA
    let mut active = 0;
    for (i, c) in text.iter().enumerate() {
        active = ((active << 1) | 1) & masks[*c as uint];
        if active & accept > 0 {
            return Some(i - m + 1);
        }
    }

    None
}
