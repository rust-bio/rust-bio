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


#[derive(Show)]
pub struct Alignment {
    pub score: i32,
    pub i: usize,
    pub j: usize,
    pub operations: Vec<AlignmentOperation>
}
