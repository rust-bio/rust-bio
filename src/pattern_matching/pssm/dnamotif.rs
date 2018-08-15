// Copyright 2018 Kieran Hervold
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use super::*;
use ndarray::prelude::Array2;
use std::f32;
use std::f32::{INFINITY, NEG_INFINITY};

/// Position-specific scoring matrix for DNA sequences
#[derive(Clone, Debug, PartialEq)]
pub struct DNAMotif {
    /// matrix holding weights at each position, indexed by [position, base]
    pub scores: Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score: f32,
    /// sum of "best" base at each position
    pub max_score: f32,
}

impl DNAMotif {
    /// Returns a Motif representing the sequences provided.
    /// # Arguments
    /// * `seqs` - sequences incorporated into motif
    /// * `pseudos` - array slice with a pseudocount for each monomer;
    ///    defaults to pssm::DEF_PSEUDO for all if None is supplied
    ///
    /// FIXME: pseudos should be an array of size MONO_CT, but that
    /// is currently impossible - see
    /// https://github.com/rust-lang/rust/issues/42863
    pub fn from_seqs(seqs: &Vec<Vec<u8>>, pseudos: Option<&[f32]>) -> Result<Self, PSSMError> {
        let w = Self::seqs_to_weights(seqs, pseudos)?;
        let mut m = DNAMotif {
            scores: w,
            min_score: 0.0,
            max_score: 0.0,
        };
        m.normalize();
        m.calc_minmax();
        Ok(m)
    }

    // helper function -- normalize self.scores
    fn normalize(&mut self) {
        for i in 0..self.len() {
            let mut tot: f32 = 0.0;
            // FIXME: slices would be cleaner
            for base_i in 0..4 {
                tot += self.scores[[i, base_i]];
            }
            for base_i in 0..4 {
                self.scores[[i, base_i]] /= tot;
            }
        }
    }

    // helper function: populate min_score and max_score
    fn calc_minmax(&mut self) {
        let pssm_len = self.len();

        // score corresponding to sum of "worst" bases at each position
        self.min_score = 0.0;
        for i in 0..pssm_len {
            // can't use the regular min/max on f32, so we use f32::min
            let min_sc = (0..4).map(|b| self.scores[[i, b]]).fold(INFINITY, f32::min);
            self.min_score += min_sc;
        }

        // score corresponding to "best" base at each position
        self.max_score = 0.0;
        for i in 0..pssm_len {
            let max_sc = (0..4)
                .map(|b| self.scores[[i, b]])
                .fold(NEG_INFINITY, f32::max);
            self.max_score += max_sc;
        }
    }
}

impl Motif for DNAMotif {
    const LK: [u8; 127] = [
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 3, 255, 255, 255, 2, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 1, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 0, 255, 3, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 1, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    ];
    const MONOS: &'static [u8] = b"ATGC";
    const MONO_CT: usize = 4;

    fn rev_lk(idx: usize) -> u8 {
        match idx {
            0 => b'A',
            1 => b'T',
            2 => b'G',
            3 => b'C',
            _ => INVALID_MONO,
        }
    }

    fn len(&self) -> usize {
        self.scores.dim().0
    }

    fn get_scores(&self) -> &Array2<f32> {
        &self.scores
    }
    fn get_min_score(&self) -> f32 {
        self.min_score
    }
    fn get_max_score(&self) -> f32 {
        self.max_score
    }
    fn get_bits() -> f32 {
        2.0
    }

    fn degenerate_consensus(&self) -> Vec<u8> {
        // derived from
        // https://github.com/biopython/biopython/blob/master/Bio/motifs/matrix.py#L205
        fn two(_a: u8, _b: u8) -> u8 {
            let (a, b) = if _b > _a { (_a, _b) } else { (_b, _a) };
            match (a, b) {
                (b'A', b'C') => b'M',
                (b'A', b'G') => b'R',
                (b'A', b'T') => b'W',
                (b'C', b'G') => b'S',
                (b'C', b'T') => b'Y',
                (b'G', b'T') => b'K',
                _ => unreachable!(), // no other combinations exist
            }
        }
        let len = self.len();
        let mut res = Vec::with_capacity(len);
        for pos in 0..len {
            let mut fracs = (0..4)
                .map(|b| (self.scores[[pos, b]], b))
                .collect::<Vec<(f32, usize)>>();
            // note: reverse sort
            fracs.sort_by(|a, b| b.partial_cmp(a).unwrap());

            res.push(if fracs[0].0 > 0.5 && fracs[0].0 > 2.0 * fracs[1].0 {
                Self::MONOS[fracs[0].1]
            } else if 4.0 * (fracs[0].0 + fracs[1].0) > 3.0 {
                two(Self::MONOS[fracs[0].1], Self::MONOS[fracs[1].1])
            } else if fracs[3].0 < EPSILON {
                let base = Self::MONOS[fracs[3].1];
                match base {
                    b'T' => b'V',
                    b'G' => b'H',
                    b'C' => b'D',
                    b'A' => b'B',
                    _ => unreachable!(), // no other bases exist
                }
            } else {
                b'N'
            });
        }
        res
    }
}

/// Return a DNAMotif wrapping an Array2 representing amino acid
/// weights at each position.  The dimensions and contents of this
/// array are unchecked, and it is incumbent on the user to ensure
/// the correct dimensions are used (ie, SEQ_LEN x 4), and no zeros
/// appear in the array.
impl From<Array2<f32>> for DNAMotif {
    fn from(scores: Array2<f32>) -> Self {
        let mut m = DNAMotif {
            scores: scores,
            min_score: 0.0,
            max_score: 0.0,
        };
        m.normalize();
        m.calc_minmax();
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pattern_matching::pssm::ScoredPos;
    #[test]
    fn simple_pssm() {
        let pssm: DNAMotif = DNAMotif::from_seqs(
            vec![
                b"AAAA".to_vec(),
                b"TTTT".to_vec(),
                b"GGGG".to_vec(),
                b"CCCC".to_vec(),
            ].as_ref(),
            None,
        ).unwrap();
        assert_eq!(pssm.scores, Array2::from_elem((4, 4), 0.25));
    }
    #[test]
    fn find_motif() {
        let pssm = DNAMotif::from_seqs(vec![b"ATGC".to_vec()].as_ref(), None).unwrap();
        let seq = b"GGGGATGCGGGG";
        if let Ok(ScoredPos {
            ref loc, ref sum, ..
        }) = pssm.score(seq)
        {
            assert_eq!(*loc, 4);
            assert_eq!(*sum, 1.0);
        } else {
            assert!(false);
        }
    }

    #[test]
    fn test_info_content() {
        // matrix w/ 100% match to A at each position
        let pssm =
            DNAMotif::from_seqs(vec![b"AAAA".to_vec()].as_ref(), Some(&[0.0, 0.0, 0.0, 0.0]))
                .unwrap();
        // 4 bases * 2 bits per base = 8
        assert_eq!(pssm.info_content(), 8.0);
    }

    #[test]
    fn test_mono_err() {
        let pssm = DNAMotif::from_seqs(vec![b"ATGC".to_vec()].as_ref(), None).unwrap();
        assert_eq!(
            pssm.score(b"AAAAXAAAAAAAAA"),
            Err(PSSMError::InvalidMonomer(b'X'))
        );
    }

    #[test]
    fn test_inconsist_err() {
        assert_eq!(
            DNAMotif::from_seqs(
                vec![b"AAAA".to_vec(), b"TTTT".to_vec(), b"C".to_vec()].as_ref(),
                Some(&[0.0; 4])
            ),
            Err(PSSMError::InconsistentLen)
        );
    }
}
