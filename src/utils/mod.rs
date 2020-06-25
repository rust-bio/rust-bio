// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Common utilities.

mod fastexp;
pub use self::fastexp::FastExp;

mod text;
pub use self::text::{trim_newline, Text, TextSlice};

mod interval;
pub use self::interval::Interval;

/// In place implementation of scan over a slice.
pub fn scan<T: Copy, F: Fn(T, T) -> T>(a: &mut [T], op: F) {
    let mut s = a[0];
    for v in a.iter_mut().skip(1) {
        s = op(s, *v);
        *v = s;
    }
}

// Inplace implementation of prescan over a slice.
pub fn prescan<T: Copy, F: Fn(T, T) -> T>(a: &mut [T], neutral: T, op: F) {
    let mut s = neutral;
    for v in a.iter_mut() {
        let t = *v;
        *v = s;
        s = op(s, t);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scan() {
        let mut a = vec![1, 0, 0, 1];
        scan(&mut a[..], |a, b| a + b);
        assert_eq!(a, vec![1, 1, 1, 2]);
    }

    #[test]
    fn test_prescan() {
        let mut a = vec![1, 0, 0, 1];
        prescan(&mut a[..], 0, |a, b| a + b);
        assert_eq!(a, vec![0, 1, 1, 1]);
    }
}
