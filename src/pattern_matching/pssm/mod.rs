// Copyright 2018 Kieran Hervold
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Implementation of position-specific scoring matrix (PSSM), aka position
//! weight matrix (PWM), for both DNA and amino-acid sequences
//!
//! Matrices can either be built from example sequences, as shown below,
//! or directly from a float matrix (type Array2<f32>).
//!
//! # Examples
//!
//! use bio::pattern_matching::pssm::DNAMotif;
//! let pssm = DNAMotif::from(vec![
//!            b"AAAA".to_vec(),
//!            b"AATA".to_vec(),
//!            b"AAGA".to_vec(),
//!            b"AAAA".to_vec(),
//!        ]);
//! let start_pos = pssm.score(b"CCCCCAATA").unwrap().loc;
//! println!("motif found at position {}", start_pos);

use ndarray::prelude::Array2;
use std::f32::NEG_INFINITY;
use utils::IntoTextIterator;

mod dnamotif;
mod protmotif;

pub use self::dnamotif::DNAMotif;
pub use self::protmotif::ProtMotif;

pub const EPSILON: f32 = 1e-5;
pub const INVALID_MONO: u8 = 255;

/// Errors
#[derive(Debug, Clone)]
pub enum PSSMError {
    /// attempted to build a pattern from sequences with mismatched lengths
    InconsistentLen,
    /// target sequence is shorter than pattern
    TargetTooShort,
    /// unknown monomer, eg, DNA base other than ATGC
    InvalidMonomer,
    /// information-free pattern
    NullPattern,
}

/// Represents motif score & location of match
#[derive(Debug, Clone)]
pub struct ScoredPos {
    pub loc: usize,
    pub sum: f32,
    pub scores: Vec<f32>,
}

impl Default for ScoredPos {
    fn default() -> ScoredPos {
        ScoredPos {
            loc: 0,
            sum: NEG_INFINITY,
            scores: Vec::new(),
        }
    }
}

/// Trait containing code shared between DNA and protein implementations.
pub trait Motif {
    /// lookup table mapping monomer -> index
    const LK: [u8; 127] = [INVALID_MONO; 127];
    const MONOS: &'static [u8] = b"";

    /// Use lookup table LK to find index of given monomer in the scores matrix.
    fn lookup(mono: u8) -> Result<usize,PSSMError> {
        if mono >= 127 {
            Err(PSSMError::InvalidMonomer)
        } else {
            let idx = Self::LK[mono as usize];
            if idx == INVALID_MONO {
                Err(PSSMError::InvalidMonomer)
            } else {
                Ok(idx as usize)
            }
        }
    }

    /// Reverse lookup: given index, return monomer.
    fn rev_lk(idx: usize) -> u8;

    /// Length of motif
    fn len(&self) -> usize;

    /// Represent motif using ambiguous codes.
    /// Primarily useful for DNA motifs, where ambiguous codes are
    /// common (eg, 'M' for 'A or C'); less so for proteins, where we
    /// represent any position without a dominant amino acid as an 'X'
    fn degenerate_consensus(&self) -> Result<Vec<u8>,PSSMError>;

    /// Helper method: accessor for scores matrix.
    fn get_scores(&self) -> &Array2<f32>;

    /// Helper method: sum of "worst" base at each position
    fn get_min_score(&self) -> f32;

    /// Helper method: sum of "best" base at each position
    fn get_max_score(&self) -> f32;

    /// Helper method: information content of a single position.  Used by info_content
    /// `info_content` method.
    /// FIXME: this should be replaced with a CTFE ... or maybe just a constant
    fn get_bits() -> f32;

    /// Standard PSSM scoring is calibrated to (1) pattern len and (2) the min and
    /// max possible scores.  This is just a un-normalized sum of matching bases, useful for
    /// comparing matches from motifs of different lengths.
    fn raw_score<'a, T: IntoTextIterator<'a>>(&self, seq_it: T) -> Result<(usize, f32, Vec<f32>), PSSMError> {
        let pssm_len = self.len();

        let mut best_start = 0;
        let mut best_score = -1.0;
        let mut best_m = Vec::new();
        // we have to look at slices, so a simple iterator won't do
        let seq = seq_it.into_iter().cloned().collect::<Vec<u8>>();
        let scores = self.get_scores();
        for start in 0..seq.len() - pssm_len + 1 {
            let m: Vec<f32> = match (0..pssm_len)
                .map(|i| match Self::lookup(seq[start + i]) {
                    Err(e) => Err(e),
                    Ok(pos) => Ok(scores[[i, pos]])
                })
                .collect() {
                    Ok(m) => m,
                    Err(e) => return Err(e),
                };
            let tot = m.iter().sum();
            if tot > best_score {
                best_score = tot;
                best_start = start;
                best_m = m;
            }
        }
        Ok((best_start, best_score, best_m))
    }

    /// Apply PSSM to sequence, finding the offset with the highest score.
    /// Return None if sequence is too short
    /// see:
    ///   MATCHTM: a tool for searching transcription factor binding sites in DNA sequences
    ///   Nucleic Acids Res. 2003 Jul 1; 31(13): 3576–3579
    ///   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169193/
    ///
    /// # Example
    /// let pssm = DNAMotif::from(vec![
    ///            b"AAAA".to_vec(),
    ///            b"AATA".to_vec(),
    ///            b"AAGA".to_vec(),
    ///            b"AAAA".to_vec(),
    ///        ]);
    /// let start_pos = pssm.score(b"CCCCCAATA").unwrap().loc;
    fn score<'a, T: IntoTextIterator<'a>>(&self, seq_it: T) -> Result<ScoredPos,PSSMError> {
        let pssm_len = self.len();
        let seq = seq_it.into_iter().cloned().collect::<Vec<u8>>();
        if seq.len() < pssm_len {
            return Err(PSSMError::TargetTooShort);
        }
        let min_score = self.get_min_score();
        let max_score = self.get_max_score();

        if max_score == min_score {
            return Err(PSSMError::NullPattern);
        }

        let (best_start, best_score, best_m) = self.raw_score(&seq)?;

        Ok(ScoredPos {
            loc: best_start,
            sum: (best_score - min_score) / (max_score - min_score),
            scores: best_m,
        })
    }

    /// Represents the information content of a motif; roughly the inverse of Shannon Entropy.
    /// Adapted from the information content described here:
    ///    https://en.wikipedia.org/wiki/Sequence_logo#Logo_creation
    fn info_content(&self) -> f32 {
        fn ent<'a, I>(probs: I) -> f32
        where
            I: Iterator<Item = &'a f32>,
        {
            probs
                .map(|p| {
                    if *p == 0.0 {
                        0.0
                    } else {
                        -1.0 * *p * p.log(2.0)
                    }
                })
                .sum()
        }
        let bits = Self::get_bits();
        let scores = self.get_scores();
        let mut tot = 0.0;
        for row in scores.genrows() {
            tot += bits - ent(row.iter());
        }
        tot
    }
}
