// Copyright 2018 Kieran Hervold
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Create a weight matrix representing a set of aligned reference sequences
//! that constitute a motif, and use this matrix to scan query sequences for
//! occurrences of this motif.
//! Complexity: O(n*m) for motif length n and query length m
//!
//! The position-specific scoring matrix (PSSM), aka position weight matrix (PWM),
//! algorithm is implemented for both DNA and amino-acid sequences.
//!
//! # Examples
//!
//! use bio::pattern_matching::pssm::DNAMotif;
//! let pssm = DNAMotif::from_seqs(vec![
//!            b"AAAA".to_vec(),
//!            b"AATA".to_vec(),
//!            b"AAGA".to_vec(),
//!            b"AAAA".to_vec(),
//!        ].as_ref(), None).unwrap();
//! let start_pos = pssm.score(b"CCCCCAATA").unwrap().loc;
//! println!("motif found at position {}", start_pos);
//!
//! /* amino acid sequences are supported, too */
//! use pssm::pattern_matching::pssm::ProtMotif;
//! let pssm = ProtMotif::from_seqs(vec![
//!            b"ARNNYM".to_vec(),
//!            b"ARNRYM".to_vec(),
//!            b"ARNNCM".to_vec(),
//!            b"ARNNYM".to_vec(),
//!        ].as_ref(), None).unwrap();

use ndarray::prelude::Array2;
use std::char;
use std::f32::NEG_INFINITY;
use utils::IntoTextIterator;

mod dnamotif;
mod protmotif;

pub use self::dnamotif::DNAMotif;
pub use self::protmotif::ProtMotif;

/// default pseudocount - used to prevent 0 tallies
pub const DEF_PSEUDO: f32 = 0.5;
/// approximately zero
pub const EPSILON: f32 = 1e-5;
/// value representing an invalid monomer in lookup table
pub const INVALID_MONO: u8 = 255;

/// Errors
quick_error! {
    #[derive(Debug, PartialEq)]
    pub enum PSSMError {
        QueryTooShort(motif_len: usize, query_len: usize) {
            description("query cannot be shorter than motif")
            display("query length {} is shorter than motif length {}", query_len, motif_len)
        }
        InconsistentLen {
            description("attempted to build a motif from sequences with mismatched lengths")
            display("mismatched sequence lengths")
        }
        InvalidMonomer(mono: u8) {
            description("unknown monomer, eg, DNA base other than ATGC")
            display("monomer '{}' is invalid", char::from(*mono))
        }
        NullMotif {
            description("a motif in which every monomer is equally likely at every position will result in a divide-by-zero exception")
            display("information-free motif")
        }
        EmptyMotif {
            description("attempted to create a motif from zero sequences")
            display("motif cannot be created from zero sequences")
        }
        InvalidPseudos(expected: u8, received: u8) {
            description("pseudo-score array should have on entry per monomer")
            display("expected pseudo-score array of length {}; got {}", expected, received)
        }
    }
}

/// Represents motif score & location of match
#[derive(Debug, Clone, PartialEq)]
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

/// Trait containing code shared between DNA and protein implementations
/// of the position-specific scoring matrix.
pub trait Motif {
    /// Lookup table mapping monomer -> index
    const LK: [u8; 127] = [INVALID_MONO; 127];
    /// All monomers, in order corresponding to lookup table
    const MONOS: &'static [u8] = b"";
    /// Monomer count - equal to length of `MONOS`
    const MONO_CT: usize = 0;

    /// Returns a weight matrix representing the sequences provided.
    /// This code is shared by implementations of `from_seqs`
    /// # Arguments
    /// * `seqs` - sequences incorporated into motif
    /// * `pseudos` - array slice with a pseudocount for each monomer;
    ///    defaults to DEF_PSEUDO for all if None is supplied
    ///
    /// FIXME: pseudos should be an array of size MONO_CT, but that
    /// is currently unsupported
    fn seqs_to_weights(
        seqs: &Vec<Vec<u8>>,
        _pseudos: Option<&[f32]>,
    ) -> Result<Array2<f32>, PSSMError> {
        let p1 = vec![DEF_PSEUDO; Self::MONO_CT];
        let pseudos = match _pseudos {
            Some(ref p2) => p2,
            None => p1.as_slice(),
        };

        if pseudos.len() != Self::MONO_CT {
            return Err(PSSMError::InvalidPseudos(
                Self::MONO_CT as u8,
                pseudos.len() as u8,
            ));
        }

        if seqs.len() == 0 {
            return Err(PSSMError::EmptyMotif);
        }

        let seqlen = seqs[0].len();
        let mut counts = Array2::zeros((seqlen, Self::MONO_CT));
        for i in 0..seqlen {
            for base in 0..Self::MONO_CT {
                counts[[i, base]] = pseudos[base];
            }
        }

        for seq in seqs.iter() {
            if seq.len() != seqlen {
                return Err(PSSMError::InconsistentLen);
            }

            for (idx, base) in seq.iter().enumerate() {
                match Self::lookup(*base) {
                    Err(e) => return Err(e),
                    Ok(pos) => counts[[idx, pos]] += 1.0,
                }
            }
        }
        Ok(counts)
    }

    /// Returns the index of given monomer in the scores matrix using the lookup table `LK`
    /// # Arguments
    /// * `mono` - monomer, eg, b'A' for DNA or b'R' for protein
    /// # Errors
    /// * `PSSMError::InvalidMonomer(mono)` - `mono` wasn't found in the lookup table
    fn lookup(mono: u8) -> Result<usize, PSSMError> {
        if mono >= 127 {
            Err(PSSMError::InvalidMonomer(mono))
        } else {
            let idx = Self::LK[mono as usize];
            if idx == INVALID_MONO {
                Err(PSSMError::InvalidMonomer(mono))
            } else {
                Ok(idx as usize)
            }
        }
    }

    /// Returns the monomer associated with the given index; the reverse of `lookup`.
    /// Returns INVALID_MONO if the index isn't associated with a monomer.
    /// # Arguments
    /// * `idx` - the index in question
    fn rev_lk(idx: usize) -> u8;

    /// Returns the length of motif
    fn len(&self) -> usize;

    /// Returns a representation of the motif using ambiguous codes.
    /// Primarily useful for DNA motifs, where ambiguous codes are
    /// common (eg, 'M' for 'A or C'); less so for proteins, where we
    /// represent any position without a dominant amino acid as an 'X'
    fn degenerate_consensus(&self) -> Vec<u8>;

    /// Accessor - returns scores matrix
    fn get_scores(&self) -> &Array2<f32>;

    /// Return sum of "worst" base at each position
    fn get_min_score(&self) -> f32;

    /// Return sum of "best" base at each position
    fn get_max_score(&self) -> f32;

    /// Returns information content of a single position.
    /// Used `info_content` method.
    /// FIXME: this should be replaced with a CTFE ... or maybe just a constant
    fn get_bits() -> f32;

    /// Returns the un-normalized sum of matching bases, useful for comparing matches from
    /// motifs of different lengths
    ///
    /// # Arguments
    /// * `seq_it` - iterator representing the query sequence
    ///
    /// # Errors
    /// * `PSSMError::InvalidMonomer(mono)` - sequence `seq_it` contained invalid monomer `mono`
    fn raw_score<'a, T: IntoTextIterator<'a>>(
        &self,
        seq_it: T,
    ) -> Result<(usize, f32, Vec<f32>), PSSMError> {
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
                    Ok(pos) => Ok(scores[[i, pos]]),
                })
                .collect()
            {
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

    /// Returns a `ScoredPos` struct representing the best match within the query sequence
    /// see:
    ///   MATCHTM: a tool for searching transcription factor binding sites in DNA sequences
    ///   Nucleic Acids Res. 2003 Jul 1; 31(13): 3576–3579
    ///   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169193/
    ///
    /// # Arguments
    /// * `seq_it` - iterator representing the query sequence
    ///
    /// # Errors
    /// * `PSSMError::InvalidMonomer(mono)` - sequence `seq_it` contained invalid monomer `mono`
    /// * `PSSMError::QueryTooShort` - sequence `seq_id` was too short
    ///
    /// # Example
    /// let pssm = DNAMotif::from_seqs(vec![
    ///            b"AAAA".to_vec(),
    ///            b"AATA".to_vec(),
    ///            b"AAGA".to_vec(),
    ///            b"AAAA".to_vec(),
    ///        ].as_ref(), None).unwrap();
    /// let start_pos = pssm.score(b"CCCCCAATA").unwrap().loc;
    fn score<'a, T: IntoTextIterator<'a>>(&self, seq_it: T) -> Result<ScoredPos, PSSMError> {
        let pssm_len = self.len();
        let seq = seq_it.into_iter().cloned().collect::<Vec<u8>>();
        if seq.len() < pssm_len {
            return Err(PSSMError::QueryTooShort(pssm_len, seq.len()));
        }
        let min_score = self.get_min_score();
        let max_score = self.get_max_score();

        if max_score == min_score {
            return Err(PSSMError::NullMotif);
        }

        let (best_start, best_score, best_m) = self.raw_score(&seq)?;

        Ok(ScoredPos {
            loc: best_start,
            sum: (best_score - min_score) / (max_score - min_score),
            scores: best_m,
        })
    }

    /// Returns a float representing the information content of a motif; roughly the
    /// inverse of Shannon Entropy.
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
