// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Alignment algorithms.

pub mod pairwise;


#[derive(Debug)]
#[derive(PartialEq)]
#[derive(Copy)]
pub enum AlignmentOperation {
    Match,
    Subst,
    Del,
    Ins
}


#[derive(Debug)]
pub struct Alignment {
    pub score: i32,
    pub i: usize,
    pub j: usize,
    pub operations: Vec<AlignmentOperation>
}


impl Alignment {
    pub fn get_cigar(&self) -> String {
        let add_op = |&: op, k, cigar: &mut String| {
            cigar.push_str(format!("{}{}", k, match op {
                AlignmentOperation::Match => "=",
                AlignmentOperation::Subst => "X",
                AlignmentOperation::Del => "D",
                AlignmentOperation::Ins => "I",
            }).as_slice());
        };
        let mut cigar = String::new();

        if self.j > 0 {
            cigar.push_str(format!("{}{}", self.j, "S").as_slice());
        }

        if !self.operations.is_empty() {
            let mut last = self.operations[0];
            let mut k = 1;
            for &op in self.operations.iter() {
                if op == last {
                    k += 1;
                } else {
                    add_op(last, k, &mut cigar);
                    k = 1;
                }
                last = op;
            }
        }

        // TODO add end clipping

        cigar
    }
}
