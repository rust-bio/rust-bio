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


/// An alignment, consisting of a score, a start in sequence y, a start in sequence x, the length
/// of sequence x and its edit operations (see alignment::pairwise for meaning of x and y).
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
