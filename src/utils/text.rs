/// Type alias for an owned text, i.e. ``Vec<u8>``.
pub type Text = Vec<u8>;
/// Type alias for a text slice, i.e. ``&[u8]``.
pub type TextSlice<'a> = &'a [u8];

/// Remove a trailing newline from the given string in place.
pub fn trim_newline(s: &mut String) {
    if s.ends_with('\n') {
        s.pop();
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Deref;

    /// This function demonstrates the use of the IntoSequenceIterator alias, which takes both
    /// slices and iterators.
    fn print_sequence<'a, Item: Deref<Target = u8>, T: IntoIterator<Item = Item>>(sequence: T) {
        for c in sequence {
            println!("{}", *c);
        }
    }

    #[test]
    fn test_print_sequence() {
        let s = b"ACGT";
        // use iterator
        print_sequence(s.iter().step_by(1));
        // use slice
        print_sequence(&s[..]);
        // use vec
        print_sequence(&vec![b'A', b'C']);
        // keep ownership
        println!("{:?}", s);
    }
}
