// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various alignment algorithms.

pub mod pairwise;


/// Alignment operations (Match, Subst, Del and Ins).
#[derive(PartialEq, Debug, Copy, Clone)]
pub enum AlignmentOperation {
    Match,
    Subst,
    Del,
    Ins
}


/// An alignment, consisting of a score, a start in sequence y, a start in sequence x, a length
/// and its edit operations (see alignment::pairwise for meaning of x and y).
#[derive(Debug)]
pub struct Alignment {
    pub score: i32,
    pub ystart: usize,
    pub xstart: usize,
    pub xlen: usize,
    pub operations: Vec<AlignmentOperation>
}


impl Alignment {
    /// Calculate the cigar string.
    pub fn cigar(&self, hard_clip: bool) -> String {
        let add_op = |op, k, cigar: &mut String| {
            cigar.push_str(&format!("{}{}", k, match op {
                AlignmentOperation::Match => "=",
                AlignmentOperation::Subst => "X",
                AlignmentOperation::Del => "D",
                AlignmentOperation::Ins => "I",
            }));
        };

        let op_len = |op: AlignmentOperation| {
            (op == AlignmentOperation::Match || op == AlignmentOperation::Subst || op == AlignmentOperation::Ins) as usize
        };

        let clip_str = if hard_clip {"H"} else {"S"};

        let mut cigar = String::new();

        if !self.operations.is_empty() {
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
        else {
            cigar.push_str(&format!("{}{}", self.xlen, clip_str));
        }

        cigar
    }

    /// Return the pretty formatted alignment as a string.
    pub fn pretty(&self, x: &[u8], y: &[u8]) -> String {
        let mut x_pretty = String::new();
        let mut y_pretty = String::new();
        let mut inb_pretty = String::new();

        if !self.operations.is_empty() {
            let mut x_add = &[u8];
            let mut y_add = String::new();
            let mut inb_add = String::new();

            for i in 1..self.operations.len() {
                match self.operations[i] {
                    AlignmentOperation::Match => {
                        x_add = &x[i];
                        inb_add = "|";
                        y_add = y[i];
                    },
                    AlignmentOperation::Subst => {
                        x_add = x[i];
                        inb_add = " ";
                        y_add = y[i];
                    },
                    AlignmentOperation::Del => {
                        x_add = "-";
                        inb_add = " ";
                        y_add = y[i];
                    },
                    AlignmentOperation::Ins => {
                        x_add = x[i];
                        inb_add = " ";
                        y_add = "-";
                    },
                }

                x_pretty.push_str(x_add);
                inb_pretty.push_str(inb_add);
                y_pretty.push_str(y_add);
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
        let alignment = Alignment { score: 5, xstart: 3, ystart: 0, xlen: 10, operations: vec![Match, Match, Match, Subst, Ins, Ins, Del, Del] };
        assert_eq!(alignment.cigar(false), "3S3=1X2I2D4S");
    }
}
