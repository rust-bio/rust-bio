/// Type alias for an owned text, i.e. ``Vec<u8>``.
pub type Text = Vec<u8>;
/// Type alias for a text slice, i.e. ``&[u8]``.
pub type TextSlice<'a> = &'a [u8];

/// Type alias for an iterator over a sequence, i.e. ``Iterator<Item=&u8>``.
pub trait TextIterator<'a>: Iterator<Item = &'a u8> {}
impl<'a, I: Iterator<Item = &'a u8>> TextIterator<'a> for I {}

/// Type alias for a type that can be coerced into a `TextIterator`.
/// This includes ``&Vec<u8>``, ``&[u8]``, ``Iterator<Item=&u8>``.
pub trait IntoTextIterator<'a>: IntoIterator<Item = &'a u8> {}
impl<'a, T: IntoIterator<Item = &'a u8>> IntoTextIterator<'a> for T {}


/// Remove a trailing newline from the given string in place.
pub fn trim_newline(s: &mut String) {
    if s.ends_with('\n') {
        s.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    /// This function demonstrates the use of the IntoSequenceIterator alias, which takes both
    /// slices and iterators.
    fn print_sequence<'a, I: IntoTextIterator<'a>>(sequence: I) {
        for c in sequence {
            println!("{}", c);
        }
    }

    #[test]
    fn test_print_sequence() {
        let s = b"ACGT";
        // use iterator
        print_sequence(s.iter().step(1));
        // use slice
        print_sequence(&s[..]);
        // use vec
        print_sequence(&vec![b'A', b'C']);
        // keep ownership
        println!("{:?}", s);
    }
}
    
