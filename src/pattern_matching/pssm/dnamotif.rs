// Copyright 2018 Kieran Hervold
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use super::*;
use ndarray::prelude::Array2;
use std::f32;
use std::f32::{INFINITY, NEG_INFINITY};

/// monomers, ie, DNA bases
/// default pseudocount - used to prevent 0 tallies
const DEF_PSEUDO: f32 = 0.5;

#[derive(Clone)]
pub struct DNAMotif {
    pub seq_ct: usize,
    pub scores: Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score: f32,
    /// sum of "best" base at each position
    pub max_score: f32,
}

impl DNAMotif {
    pub fn from_seqs_with_pseudocts(seqs: Vec<Vec<u8>>, pseudos: &[f32; 4]) -> Result<DNAMotif,PSSMError> {
        // null case
        if seqs.len() == 0 {
            return Ok(DNAMotif {
                seq_ct: 0,
                scores: Array2::zeros((0, 0)),
                min_score: 0.0,
                max_score: 0.0,
            });
        }
        let seqlen = seqs[0].len();
        let mut counts = Array2::zeros((seqlen, 4));
        for i in 0..seqlen {
            for base in 0..4 {
                counts[[i, base]] = pseudos[base];
            }
        }

        for seq in seqs.iter() {
            if seq.len() != seqlen {
                return Err(PSSMError::InconsistentLen)
            }

            for (idx, base) in seq.iter().enumerate() {
                match Self::lookup(*base) {
                    Err(e) => return Err(e),
                    Ok(pos) => { counts[[idx, pos]] += 1.0 }
                }
            }
        }
        let mut m = DNAMotif {
            seq_ct: seqs.len(),
            scores: counts,
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

    // helper function
    pub fn calc_minmax(&mut self) {
        let pssm_len = self.len();

        // score corresponding to sum of "worst" bases at each position
        // FIXME: iter ...
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

    /// derived from
    /// https://github.com/biopython/biopython/blob/master/Bio/motifs/matrix.py#L205
    fn degenerate_consensus(&self) -> Result<Vec<u8>,PSSMError> {
        fn two(_a: u8, _b: u8) -> Result<u8,PSSMError> {
            let (a, b) = if _b > _a { (_a, _b) } else { (_b, _a) };
            match (a, b) {
                (b'A', b'C') => Ok(b'M'),
                (b'A', b'G') => Ok(b'R'),
                (b'A', b'T') => Ok(b'W'),
                (b'C', b'G') => Ok(b'S'),
                (b'C', b'T') => Ok(b'Y'),
                (b'G', b'T') => Ok(b'K'),
                _ => return Err(PSSMError::InvalidMonomer),
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
                two(Self::MONOS[fracs[0].1], Self::MONOS[fracs[1].1])?
            } else if fracs[3].0 < EPSILON {
                match Self::MONOS[fracs[3].1] {
                    b'T' => b'V',
                    b'G' => b'H',
                    b'C' => b'D',
                    b'A' => b'B',
                    _ => return Err(PSSMError::InvalidMonomer)
                }
            } else {
                b'N'
            });
        }
        Ok(res)
    }
}

/// calculate scores matrix from a list of equal-length sequences
/// use DEF_PSEUDO as default pseudocount
impl From<Vec<Vec<u8>>> for DNAMotif {
    fn from(seqs: Vec<Vec<u8>>) -> Self {
        DNAMotif::from_seqs_with_pseudocts(seqs, &[DEF_PSEUDO, DEF_PSEUDO, DEF_PSEUDO, DEF_PSEUDO]).expect("from_seqs_with_pseudocts failed")
    }
}

impl From<Array2<f32>> for DNAMotif {
    fn from(scores: Array2<f32>) -> Self {
        let mut m = DNAMotif {
            seq_ct: 0,
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
        let pssm = DNAMotif::from(vec![
            b"AAAA".to_vec(),
            b"TTTT".to_vec(),
            b"GGGG".to_vec(),
            b"CCCC".to_vec(),
        ]);
        assert_eq!(pssm.scores, Array2::from_elem((4, 4), 0.25));
    }
    #[test]
    fn find_motif() {
        let pssm = DNAMotif::from(vec![b"ATGC".to_vec()]);
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
            DNAMotif::from_seqs_with_pseudocts(vec![b"AAAA".to_vec()], &[0.0, 0.0, 0.0, 0.0]).unwrap();
        // 4 bases * 2 bits per base = 8
        assert_eq!(pssm.info_content(), 8.0);
    }
}
