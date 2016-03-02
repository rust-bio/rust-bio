// Copyright 2014-2015 Johannes Köster, Vadim Nazarov.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various alignment and distance computing algorithms.


pub mod pairwise;
pub mod distance;


/// Alignment operations (Match, Subst, Del and Ins).
#[derive(PartialEq, Debug, Copy, Clone)]
pub enum AlignmentOperation {
    Match,
    Subst,
    Del,
    Ins,
}


/// An alignment, consisting of a score, a start in sequence y, a start in sequence x, the length
/// of sequence x and its edit operations (see alignment::pairwise for meaning of x and y).
#[derive(Debug)]
pub struct Alignment {
    pub score: i32,
    pub ystart: usize,
    pub xstart: usize,
    pub xlen: usize,
    pub operations: Vec<AlignmentOperation>,
}


impl Alignment {
    /// Calculate the cigar string.
    pub fn cigar(&self, hard_clip: bool) -> String {
        let add_op = |op, k, cigar: &mut String| {
            cigar.push_str(&format!("{}{}",
                                    k,
                                    match op {
                                        AlignmentOperation::Match => "=",
                                        AlignmentOperation::Subst => "X",
                                        AlignmentOperation::Del => "D",
                                        AlignmentOperation::Ins => "I",
                                    }));
        };

        let op_len = |op: AlignmentOperation| {
            (op == AlignmentOperation::Match || op == AlignmentOperation::Subst ||
             op == AlignmentOperation::Ins) as usize
        };

        let clip_str = if hard_clip {
            "H"
        } else {
            "S"
        };

        let mut cigar = String::new();

        if self.operations.is_empty() {
            cigar.push_str(&format!("{}{}", self.xlen, clip_str));
        } else {
            if self.xstart > 0 {
                cigar.push_str(&format!("{}{}", self.xstart, clip_str));
            }

            let mut last = self.operations[0];
            let mut k = 1;
            let mut alen = op_len(last);
            for &op in self.operations[1..].iter() {
                if op == last {
                    k += 1;
                } else {
                    add_op(last, k, &mut cigar);
                    k = 1;
                }
                last = op;
                alen += op_len(op);
            }
            add_op(last, k, &mut cigar);

            let clip = self.xlen - alen;
            if clip > 0 {
                cigar.push_str(&format!("{}{}", clip, clip_str));
            }
        }

        cigar
    }

    /// Return the pretty formatted alignment as a String.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alignment::pairwise::*;
    ///
    /// let x = b"CCGTCCGGCAA";
    /// let y = b"AAAAACCGTTGACGCAA";
    /// let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    ///
    /// // ------CCGTCCGGCAA
    /// //       |  |   ||||
    /// // AAAAACCGTTGACGCAA
    /// let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
    /// let alignment = aligner.semiglobal(x, y);
    /// println!("Semiglobal: \n{}\n", alignment.pretty(x, y));
    ///
    /// // -----CCGTCCGGCAA-
    /// //      ||||
    /// // AAAAACCGTTGACGCAA
    /// let alignment = aligner.local(x, y);
    /// println!("Local: \n{}\n", alignment.pretty(x, y));
    ///
    /// // ------CCGTCCGGCAA
    /// //       |  |   ||||
    /// // AAAAACCGTTGACGCAA
    /// let alignment = aligner.global(x, y);
    /// println!("Global: \n{}\n", alignment.pretty(x, y));
    /// ```
    pub fn pretty(&self, x: &[u8], y: &[u8]) -> String {
        let mut x_pretty = String::new();
        let mut y_pretty = String::new();
        let mut inb_pretty = String::new();

        if !self.operations.is_empty() {
            let mut x_i: usize = self.xstart;
            let mut y_i: usize = self.ystart;

            // Add '-' before aligned subsequences and un-aligned 5' substrings of sequences.
            if x_i > y_i {
                let diff = x_i - y_i;
                x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&x[0..x_i])));
                for _ in 0..diff {
                    y_pretty.push('-');
                }
                y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&x[0..y_i])));
            } else if x_i < y_i {
                let diff = y_i - x_i;
                for _ in 0..diff {
                    x_pretty.push('-');
                }
                x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&x[0..x_i])));
                y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&y[0..y_i])));
            } else {
                x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&x[0..x_i])));
                y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&y[0..y_i])));
            }
            for _ in 0..x_pretty.len() {
                inb_pretty.push(' ');
            }

            // Process the alignment.
            for i in 0..self.operations.len() {
                match self.operations[i] {
                    AlignmentOperation::Match => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push_str("|");

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Subst => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push(' ');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Del => {
                        x_pretty.push('-');

                        inb_pretty.push(' ');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Ins => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push(' ');

                        y_pretty.push('-');
                    }
                }
            }

            // Add un-aligned 3' substrings of sequences.
            for &item in x.iter().skip(x_i) {
                x_pretty.push_str(&String::from_utf8_lossy(&[item]));
            }
            for &item in y.iter().skip(y_i) {
                y_pretty.push_str(&String::from_utf8_lossy(&[item]));
            }

            // Add trailing '-'.
            if x_pretty.len() > y_pretty.len() {
                for _ in y_pretty.len()..x_pretty.len() {
                    y_pretty.push('-');
                }
            } else {
                for _ in x_pretty.len()..y_pretty.len() {
                    x_pretty.push('-');
                }
            }
        }

        format!("{}\n{}\n{}", x_pretty, inb_pretty, y_pretty)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::AlignmentOperation::*;

    #[test]
    fn test_cigar() {
        let alignment = Alignment {
            score: 5,
            xstart: 3,
            ystart: 0,
            xlen: 10,
            operations: vec![Match, Match, Match, Subst, Ins, Ins, Del, Del],
        };
        assert_eq!(alignment.cigar(false), "3S3=1X2I2D4S");
    }
}
