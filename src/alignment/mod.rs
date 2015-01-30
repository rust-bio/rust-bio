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


pub struct Alignment {
    pub i: usize,
    pub j: usize,
    pub operations: Vec<AlignmentOperation>
}
