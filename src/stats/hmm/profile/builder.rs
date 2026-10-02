// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::alphabets::Alphabet;
use crate::stats::hmm::profile::{ProfileError, ProfileHmm, NEG_INF};

/// The byte that marks a gap in the aligned rows passed to a [`ProfileBuilder`].
pub const GAP: u8 = b'-';

/// Collects aligned rows and turns them into a [`ProfileHmm`].
///
/// A column becomes a position of the model if the rows holding a residue in it carry at least
/// half of the total row weight. Counts of residues and of transitions between the match, insert
/// and delete states are smoothed with a pseudocount before they are turned into probabilities.
/// Residues in the other columns count as insertions after the preceding position.
pub struct ProfileBuilder {
    symbols: Vec<u8>,
    rank_of: Vec<i32>,
    pseudocount: f64,
    rows: Vec<Vec<u8>>,
    weights: Vec<f64>,
}

impl ProfileBuilder {
    /// Create a builder over the residues of `alphabet` with a pseudocount of `1.0`.
    pub fn new(alphabet: &Alphabet) -> ProfileBuilder {
        let mut symbols: Vec<u8> = alphabet
            .symbols
            .iter()
            .filter(|&s| s as u8 != GAP)
            .map(|s| s as u8)
            .collect();
        symbols.sort_unstable();
        let mut rank_of = vec![-1i32; 256];
        for (rank, &sym) in symbols.iter().enumerate() {
            rank_of[sym as usize] = rank as i32;
        }
        ProfileBuilder {
            symbols,
            rank_of,
            pseudocount: 1.0,
            rows: Vec::new(),
            weights: Vec::new(),
        }
    }

    /// Set the pseudocount added to every residue and transition count. It must be positive and
    /// finite, which [`ProfileBuilder::build`] checks.
    pub fn pseudocount(mut self, value: f64) -> ProfileBuilder {
        self.pseudocount = value;
        self
    }

    /// Add one aligned row with weight `1.0`.
    pub fn add_row(&mut self, aligned: &[u8]) -> &mut ProfileBuilder {
        self.add_weighted_row(aligned, 1.0)
    }

    /// Add one aligned row. A row of weight `n` counts like `n` copies of the row; weights should
    /// not be negative.
    pub fn add_weighted_row(&mut self, aligned: &[u8], weight: f64) -> &mut ProfileBuilder {
        self.rows.push(aligned.to_vec());
        self.weights.push(weight);
        self
    }

    /// Build the model from the rows added so far.
    ///
    /// Fails if there are no rows, the rows differ in length or are empty, a row holds a residue
    /// outside the alphabet, no column qualifies as a position, or the pseudocount is not
    /// positive and finite.
    pub fn build(&self) -> Result<ProfileHmm, ProfileError> {
        if self.pseudocount <= 0.0 || !self.pseudocount.is_finite() {
            return Err(ProfileError::NonPositivePseudocount(self.pseudocount));
        }
        let num_rows = self.rows.len();
        if num_rows == 0 {
            return Err(ProfileError::NoRows);
        }
        let num_cols = self.rows[0].len();
        if num_cols == 0 {
            return Err(ProfileError::EmptyAlignment);
        }
        for row in &self.rows {
            if row.len() != num_cols {
                return Err(ProfileError::RaggedRows {
                    expected: num_cols,
                    found: row.len(),
                });
            }
            for &c in row.iter() {
                if c != GAP && (c as usize >= self.rank_of.len() || self.rank_of[c as usize] < 0) {
                    return Err(ProfileError::UnknownSymbol(c));
                }
            }
        }

        let alphabet_size = self.symbols.len();
        let pc = self.pseudocount;

        let total_weight: f64 = self.weights.iter().sum();
        let mut is_match = vec![false; num_cols];
        for c in 0..num_cols {
            let non_gap_weight: f64 = self
                .rows
                .iter()
                .zip(self.weights.iter())
                .filter(|(row, _)| row[c] != GAP)
                .map(|(_, &w)| w)
                .sum();
            is_match[c] = non_gap_weight * 2.0 >= total_weight;
        }
        let col_to_position: Vec<usize> = {
            let mut pos = 0usize;
            is_match
                .iter()
                .map(|&m| {
                    if m {
                        pos += 1;
                        pos
                    } else {
                        pos
                    }
                })
                .collect()
        };
        let num_match: usize = is_match.iter().filter(|&&m| m).count();
        if num_match == 0 {
            return Err(ProfileError::NoMatchColumns);
        }
        let l = num_match;

        let mut emit_counts = vec![vec![0.0f64; alphabet_size]; l + 1];
        let mut background = vec![0.0f64; alphabet_size];
        let mut count_mm = vec![0.0f64; l + 1];
        let mut count_mi = vec![0.0f64; l + 1];
        let mut count_md = vec![0.0f64; l + 1];
        let mut count_im = vec![0.0f64; l + 1];
        let mut count_ii = vec![0.0f64; l + 1];
        let mut count_id = vec![0.0f64; l + 1];
        let mut count_dm = vec![0.0f64; l + 1];
        let mut count_di = vec![0.0f64; l + 1];
        let mut count_dd = vec![0.0f64; l + 1];

        #[derive(Clone, Copy, PartialEq)]
        enum State {
            M,
            I,
            D,
        }

        for (row, &weight) in self.rows.iter().zip(self.weights.iter()) {
            let mut prev_state = State::M;
            let mut prev_node = 0usize;
            for c in 0..num_cols {
                let symbol = row[c];
                if is_match[c] {
                    let node = col_to_position[c];
                    if symbol != GAP {
                        let rank = self.rank_of[symbol as usize] as usize;
                        emit_counts[node][rank] += weight;
                        background[rank] += weight;
                        match prev_state {
                            State::M => count_mm[prev_node] += weight,
                            State::I => count_im[prev_node] += weight,
                            State::D => count_dm[prev_node] += weight,
                        }
                        prev_state = State::M;
                    } else {
                        match prev_state {
                            State::M => count_md[prev_node] += weight,
                            State::I => count_id[prev_node] += weight,
                            State::D => count_dd[prev_node] += weight,
                        }
                        prev_state = State::D;
                    }
                    prev_node = node;
                } else if symbol != GAP {
                    let rank = self.rank_of[symbol as usize] as usize;
                    background[rank] += weight;
                    match prev_state {
                        State::M => count_mi[prev_node] += weight,
                        State::I => count_ii[prev_node] += weight,
                        State::D => count_di[prev_node] += weight,
                    }
                    prev_state = State::I;
                }
            }
            match prev_state {
                State::M => count_mm[l] += weight,
                State::I => count_im[l] += weight,
                State::D => count_dm[l] += weight,
            }
        }

        let null = log_normalize(&add_pseudo(&background, pc));

        let mut emission = vec![vec![NEG_INF; alphabet_size]; l + 1];
        for k in 1..=l {
            emission[k] = log_normalize(&add_pseudo(&emit_counts[k], pc));
        }

        let mut t_mm = vec![NEG_INF; l + 1];
        let mut t_mi = vec![NEG_INF; l + 1];
        let mut t_md = vec![NEG_INF; l + 1];
        let mut t_im = vec![NEG_INF; l + 1];
        let mut t_ii = vec![NEG_INF; l + 1];
        let mut t_id = vec![NEG_INF; l + 1];
        let mut t_dm = vec![NEG_INF; l + 1];
        let mut t_di = vec![NEG_INF; l + 1];
        let mut t_dd = vec![NEG_INF; l + 1];

        for k in 0..=l {
            let last = k == l;
            let m_dist: Vec<f64> = if last {
                vec![count_mm[k] + pc, count_mi[k] + pc]
            } else {
                vec![count_mm[k] + pc, count_mi[k] + pc, count_md[k] + pc]
            };
            let m_logs = log_normalize(&m_dist);
            t_mm[k] = m_logs[0];
            t_mi[k] = m_logs[1];
            if !last {
                t_md[k] = m_logs[2];
            }

            let i_dist: Vec<f64> = if last {
                vec![count_im[k] + pc, count_ii[k] + pc]
            } else {
                vec![count_im[k] + pc, count_ii[k] + pc, count_id[k] + pc]
            };
            let i_logs = log_normalize(&i_dist);
            t_im[k] = i_logs[0];
            t_ii[k] = i_logs[1];
            if !last {
                t_id[k] = i_logs[2];
            }

            if k >= 1 {
                let d_dist: Vec<f64> = if last {
                    vec![count_dm[k] + pc, count_di[k] + pc]
                } else {
                    vec![count_dm[k] + pc, count_di[k] + pc, count_dd[k] + pc]
                };
                let d_logs = log_normalize(&d_dist);
                t_dm[k] = d_logs[0];
                t_di[k] = d_logs[1];
                if !last {
                    t_dd[k] = d_logs[2];
                }
            }
        }

        Ok(ProfileHmm {
            num_cols: l,
            symbols: self.symbols.clone(),
            rank_of: self.rank_of.clone(),
            emission,
            null,
            t_mm,
            t_mi,
            t_md,
            t_im,
            t_ii,
            t_id,
            t_dm,
            t_di,
            t_dd,
        })
    }
}

fn add_pseudo(counts: &[f64], pc: f64) -> Vec<f64> {
    counts.iter().map(|&c| c + pc).collect()
}

fn log_normalize(weights: &[f64]) -> Vec<f64> {
    let total: f64 = weights.iter().sum();
    weights
        .iter()
        .map(|&w| if w <= 0.0 { NEG_INF } else { (w / total).ln() })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn builder() -> ProfileBuilder {
        ProfileHmm::builder(&Alphabet::new(b"ACGT"))
    }

    #[test]
    fn test_pseudocount_keeps_unseen_residues_scorable() {
        let mut builder = builder();
        builder.add_row(b"ACGT");
        let profile = builder.build().unwrap();
        assert!(profile.forward(b"ACGTA").unwrap().is_finite());
        assert!(profile.bit_score(b"ACGTA").unwrap().is_finite());
        assert!(profile.forward(b"ACG").unwrap().is_finite());
        assert!(profile.bit_score(b"ACG").unwrap().is_finite());
        assert!(profile
            .best_alignment(b"ACGTA")
            .unwrap()
            .log_prob
            .is_finite());
    }

    #[test]
    fn test_pseudocount_changes_scores() {
        let mut low = builder().pseudocount(0.5);
        let mut high = builder().pseudocount(10.0);
        for row in [b"ACGTAC", b"ACGTAC", b"ACGTAC"] {
            low.add_row(row);
            high.add_row(row);
        }
        let low = low.build().unwrap().bit_score(b"ACGTAC").unwrap();
        let high = high.build().unwrap().bit_score(b"ACGTAC").unwrap();
        assert!((low - high).abs() > 1e-6);
    }

    #[test]
    fn test_match_column_needs_half_of_the_rows() {
        let mut half = builder();
        half.add_row(b"AA");
        half.add_row(b"A-");
        assert_eq!(half.build().unwrap().num_columns(), 2);

        let mut below = builder();
        below.add_row(b"AAA");
        below.add_row(b"A--");
        below.add_row(b"A--");
        below.add_row(b"A--");
        assert_eq!(below.build().unwrap().num_columns(), 1);
    }

    #[test]
    fn test_weight_equals_repetition() {
        let mut repeated = builder();
        repeated.add_row(b"ACGTAC");
        repeated.add_row(b"ACGTAC");
        repeated.add_row(b"ACGTAC");
        repeated.add_row(b"AC--AC");
        repeated.add_row(b"AGGTAC");
        let repeated = repeated.build().unwrap();

        let mut weighted = builder();
        weighted.add_weighted_row(b"ACGTAC", 3.0);
        weighted.add_row(b"AC--AC");
        weighted.add_row(b"AGGTAC");
        let weighted = weighted.build().unwrap();

        for query in [&b"ACGTAC"[..], &b"ACAC"[..], &b"ACGTTTAC"[..]] {
            let a = repeated.bit_score(query).unwrap();
            let b = weighted.bit_score(query).unwrap();
            assert!((a - b).abs() < 1e-9);
        }
    }

    #[test]
    fn test_weights_decide_match_columns() {
        let mut builder = builder();
        builder.add_weighted_row(b"AA", 1.0);
        builder.add_weighted_row(b"A-", 3.0);
        assert_eq!(builder.build().unwrap().num_columns(), 1);
    }

    #[test]
    fn test_invalid_input_is_rejected() {
        assert_eq!(builder().build().unwrap_err(), ProfileError::NoRows);

        let mut empty = builder();
        empty.add_row(b"");
        assert_eq!(empty.build().unwrap_err(), ProfileError::EmptyAlignment);

        let mut ragged = builder();
        ragged.add_row(b"ACG");
        ragged.add_row(b"AC");
        assert_eq!(
            ragged.build().unwrap_err(),
            ProfileError::RaggedRows {
                expected: 3,
                found: 2
            }
        );

        let mut unknown = builder();
        unknown.add_row(b"ACX");
        assert_eq!(
            unknown.build().unwrap_err(),
            ProfileError::UnknownSymbol(b'X')
        );

        let mut gaps = builder();
        gaps.add_row(b"---");
        assert_eq!(gaps.build().unwrap_err(), ProfileError::NoMatchColumns);

        let mut bad = builder().pseudocount(0.0);
        bad.add_row(b"ACG");
        assert_eq!(
            bad.build().unwrap_err(),
            ProfileError::NonPositivePseudocount(0.0)
        );
    }

    #[test]
    fn test_gap_is_not_part_of_the_alphabet() {
        let mut builder = ProfileHmm::builder(&Alphabet::new(b"AC-"));
        builder.add_row(b"AC");
        let profile = builder.build().unwrap();
        assert_eq!(
            profile.forward(b"-").unwrap_err(),
            ProfileError::UnknownSymbol(b'-')
        );
    }
}
