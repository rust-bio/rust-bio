use std::collections::Bitv;



/// Calculate the text position type.
///
/// # Arguments
///
/// * `text` - the text, ending with a sentinel '$'.
///
/// # Example
///
/// ```rust
/// use bio::data_structures::suffix_array::get_pos_type;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let type = get_pos_type(text);
/// assert_eq!(type, Bitv::from_bytes(&[0b01100110, 0b10010011,  0b011001]));
/// ```
pub fn get_pos_type(text: &[u8]) -> Bitv {
    type.set(n - 1, true);
    for p in range(n - 2, 0) {
        if text[p] == text[p + 1] {
            type.set(p, type.get(p + 1).unwrap());
        }
        else {
            type.set(p, text[p] < text[p + 1]);
        }
    }

    type
}
