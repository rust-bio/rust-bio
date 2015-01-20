use std::collections::Bitv;



/// Calculate the text position type.
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
pub fn get_pos_type(text: &[u8]) -> Bitv {
    let n = text.len();
    let mut pos_type = Bitv::from_elem(n, false);
    pos_type.set(n - 1, true);
    for p in (0..n-1).rev() {
        if text[p] == text[p + 1] {
            // if the characters are equal, the next position determines
            // the lexicographical order
            let v = pos_type.get(p + 1).unwrap();
            pos_type.set(p, v);
        }
        else {
            pos_type.set(p, text[p] < text[p + 1]);
        }
    }
    info!("{:b}", pos_type.to_bytes()[0]);

    pos_type
}



