// Copyright 2015 Peer Aramillo Irizar.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Implementation of the protein alphabet.
//!
//! # Example
//!
//! ```
//! use bio::alphabets;
//! let alphabet = alphabets::protein::alphabet();
//! assert!(alphabet.is_word(b"DEQsga"));
//! assert!(!alphabet.is_word(b"BzJ"));
//! ```

use crate::alphabets::Alphabet;

/// Returns the standard protein alphabet, containing the 20 common amino acids.
pub fn alphabet() -> Alphabet {
    Alphabet::new(&b"ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv"[..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_word() {
        assert!(alphabet().is_word(b"PRSkl"));
    }

    #[test]
    fn is_no_word() {
        assert!(!alphabet().is_word(b"Bb"));
    }

    #[test]
    fn symbol_is_no_word() {
        assert!(!alphabet().is_word(b"#"));
    }

    #[test]
    fn number_is_no_word() {
        assert!(!alphabet().is_word(b"42"));
    }
}
