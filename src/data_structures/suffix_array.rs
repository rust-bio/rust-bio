use std::collections::Bitv;

fn step2(text: &[usize], pos_types: Bitv, lms_pos: &[usize]) -> Vec<usize> {
    let n = text.len();
    let mut pos = Vec::with_capacity(text.len());

    let mut bucket_start = get_bucket_start(text);
    let mut bucket_end = bucket_start.clone()[1..].to_vec();
    bucket_end.push(text.len());
    let mut bucket_end_lms = bucket_end.clone();

    // init all positions as unknown
    for _ in text.iter() {
        pos.push(n);
    }

    // insert LMS positions to the end of their buckets
    for &p in lms_pos.iter().rev() {
        let c = text[p];
        pos[bucket_end_lms[c]] = p;
        bucket_end_lms[c] -= 1;
    }

    // insert L-positions into buckets
    for r in 0..pos.len() {
        let p = pos[r];
        if p == n {
            continue;
        }
        let prev = p - 1;
        if !pos_types.get(prev).unwrap() {
            // L-position, insert into bucket
            let c = text[prev];
            pos[bucket_start[c]] = prev;
            bucket_start[c] += 1;
        }
    }

    for r in (0..pos.len()).rev() {
        let prev = pos[r] - 1;
        if pos_types.get(prev).unwrap() {
            let c = text[prev];
            pos[bucket_end[c]] = prev;
            bucket_end[c] -= 1;
        }
    }

    pos
}


fn get_bucket_start(text: &[usize]) -> Vec<usize> {
    let alphabet_size = *text.iter().max().unwrap();
    let mut bucket_sizes = Vec::with_capacity(alphabet_size);

    for _ in (0..alphabet_size) {
        bucket_sizes.push(0);
    }

    for &c in text.iter() {
        bucket_sizes[c as usize] += 1;
    }

    let mut bucket_start = Vec::with_capacity(alphabet_size as usize);
    let mut sum = 0;
    for &size in bucket_sizes.iter() {
        bucket_start.push(sum);
        sum += size;
    }

    bucket_start
}


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
pub fn get_pos_types(text: &[usize]) -> Bitv {
    let n = text.len();
    let mut pos_types = Bitv::from_elem(n, false);
    pos_types.set(n - 1, true);
    for p in (0..n-1).rev() {
        if text[p] == text[p + 1] {
            // if the characters are equal, the next position determines
            // the lexicographical order
            let v = pos_types.get(p + 1).unwrap();
            pos_types.set(p, v);
        }
        else {
            pos_types.set(p, text[p] < text[p + 1]);
        }
    }

    pos_types
}
