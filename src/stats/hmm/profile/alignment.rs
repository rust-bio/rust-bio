// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

/// One step of an alignment between a query and a [`ProfileHmm`](super::ProfileHmm).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProfileOp {
    /// The residue is aligned to the next position of the model.
    Match(u8),
    /// The residue is carried between two positions of the model.
    Insert(u8),
    /// The next position of the model is skipped.
    Delete,
}

/// An alignment of a query against a [`ProfileHmm`](super::ProfileHmm).
#[derive(Clone, Debug, PartialEq)]
pub struct ProfileAlignment {
    /// The steps of the alignment, from the first residue and position to the last.
    pub ops: Vec<ProfileOp>,
    /// The natural-log probability of the alignment.
    pub log_prob: f64,
    /// The log-odds of the alignment against the background model, in bits.
    pub bit_score: f64,
}

impl ProfileAlignment {
    /// The number of positions that are matched to a residue.
    pub fn matched_columns(&self) -> usize {
        self.ops
            .iter()
            .filter(|op| matches!(op, ProfileOp::Match(_)))
            .count()
    }

    /// The number of positions that are skipped.
    pub fn deleted_columns(&self) -> usize {
        self.ops
            .iter()
            .filter(|op| matches!(op, ProfileOp::Delete))
            .count()
    }

    /// The number of query residues that the alignment accounts for, matched or inserted.
    pub fn consumed_residues(&self) -> usize {
        self.ops
            .iter()
            .filter(|op| matches!(op, ProfileOp::Match(_) | ProfileOp::Insert(_)))
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_counts() {
        let alignment = ProfileAlignment {
            ops: vec![
                ProfileOp::Match(b'A'),
                ProfileOp::Insert(b'G'),
                ProfileOp::Insert(b'G'),
                ProfileOp::Delete,
                ProfileOp::Match(b'C'),
            ],
            log_prob: 0.0,
            bit_score: 0.0,
        };
        assert_eq!(alignment.matched_columns(), 2);
        assert_eq!(alignment.deleted_columns(), 1);
        assert_eq!(alignment.consumed_residues(), 4);
    }
}
