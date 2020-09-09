# Copyright 2020 Tianyi Shi
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms.

"""Convert a amino acid substitution matrice taken from seqAn for use in rust-bio

Amino acid substitution matrices taken from [seqAn](https://github.com/seqan/seqan/blob/master/include%2Fseqan%2Fscore%2Fscore_matrix_data.h) has rows and columns in this order:

0  1  2 ... 22 23 24 25 26
A  B  C ...  W  Y  Z  X  *

where X, Y and Z are not in alphabetical order.

This script converts such a matrix into alphabetical order, i.e.

0  1  2 ... 22 23 24 25 26
A  B  C ...  W  X  Y  Z  *

and adds some boilerplate around it.

Usage:
python substitution-matrix-conv.py '<matrix>' <name-of-matrix>

Example:
python substitution-matrix-conv.py '[3, 0, -3, 0, 0, -4, 1, -3, -1, -2, -2, ...]' pam120
"""

import os
import sys
import json

m, name = sys.argv[1:3]
m = json.loads(m)


for i in range(27):
    a = 27 * i
    y, z, x = m[a + 23 : a + 26]
    m[a + 23 : a + 26] = [x, y, z]


y = m[27 * 23 : 27 * 24]
z = m[27 * 24 : 27 * 25]
x = m[27 * 25 : 27 * 26]

m[27 * 23 : 27 * 24] = x
m[27 * 24 : 27 * 25] = y
m[27 * 25 : 27 * 26] = z

prefix = """
use crate::scores::lookup;

const MAT: [i32; 729] = """

suffix = """;

/// Return the substitution score in {name} matrix between two amino acids (single characters as `u8`)
pub fn {name}(a: u8, b: u8) -> i32 {{
    MAT[lookup(a) * 27 + lookup(b)]
}}

#[cfg(test)]
mod tests {{
    use super::*;

    #[test]
    fn test_{name}() {{
        let score1 = {name}(b'A', b'A');
        assert_eq!(score1, {AA});
        let score2 = {name}(b'O', b'*');
        assert_eq!(score2, {O*});
        let score3 = {name}(b'A', b'*');
        assert_eq!(score3, {A*});
        let score4 = {name}(b'*', b'*');
        assert_eq!(score4, {**});
        let score5 = {name}(b'X', b'X');
        assert_eq!(score5, {XX});
        let score6 = {name}(b'X', b'Z');
        assert_eq!(score6, {XZ});
    }}
}}
"""

print(
    prefix
    + json.dumps(m)
    + suffix.format(
        **{
            "name": name,
            "AA": m[0],
            "O*": m[15 * 27 - 1],
            "A*": m[27],
            "**": m[728],
            "XX": m[23 * 27 + 23],
            "XZ": m[23 * 27 + 25],
        }
    )
)
