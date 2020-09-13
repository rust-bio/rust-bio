// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Substitution matrices taken from [SeqAn](https://github.com/seqan/seqan/blob/master/include%2Fseqan%2Fscore%2Fscore_matrix_data.h)
//!
//! Note these special characters in the alphabet:
//!
//! | character | 3-letter code |               Definition                |
//! | :-------: | :-----------: | :-------------------------------------: |
//! |     B     |      Asx      | Asparagine or Aspartic acid (Aspartate) |
//! |     Z     |      Glx      | Glutamine or Glutamic acid (Glutamate)  |
//! |     X     |      Xaa      |        Any amino acid	All codons        |
//! |     *     |      END      |  Termination codon (translation stop)   |
//!
//! # References
//!
//! - https://www.mathworks.com/help/bioinfo/ref/aminolookup.html
pub use self::blosum62::blosum62;
pub use self::pam120::pam120;
pub use self::pam200::pam200;
pub use self::pam250::pam250;
pub use self::pam40::pam40;

pub mod blosum62;
pub mod pam120;
pub mod pam200;
pub mod pam250;
pub mod pam40;

#[inline]
fn lookup(a: u8) -> usize {
    if a == b'*' {
        26 as usize
    } else {
        (a - 65) as usize
    }
}
