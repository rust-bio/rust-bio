// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub use self::blosum30::blosum30;
pub use self::blosum45::blosum45;
pub use self::blosum62::blosum62;
pub use self::pam120::pam120;
pub use self::pam200::pam200;
pub use self::pam250::pam250;
pub use self::pam40::pam40;

pub mod blosum30;
pub mod blosum45;
pub mod blosum62;
pub mod pam120;
pub mod pam200;
pub mod pam250;
pub mod pam40;

#[inline]
fn lookup(a: u8) -> usize {
    if a == b'Y' {
        23
    } else if a == b'Z' {
        24
    } else if a == b'X' {
        25
    } else if a == b'*' {
        26
    } else {
        (a - 65) as usize
    }
}
