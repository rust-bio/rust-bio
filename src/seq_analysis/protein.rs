//! The data is based on Gasteiger _et al_ (2005) i.e. in accordance with
//! [ExPASy's ProtParam Tool](https://web.expasy.org/protparam/)
//! # References
//! - Gasteiger E., Hoogland C., Gattiker A., Duvaud S., Wilkins M.R., Appel R.D., Bairoch A.
//!   Protein Identification and Analysis Tools on the ExPASy Server, in _The Proteomics Protocols
//!   Handbook_, Humana Press (2005). pp. 571-607

use crate::utils::TextSlice;
use std::collections::BTreeMap;

pub type AminoAcidCount = BTreeMap<u8, u32>;

use lazy_static::lazy_static;

// -------------- used for calculating isoelectric point -------------------
const N_TERM_PKA_DEFAULT: f32 = 7.5;
const C_TERM_PKA_DEFAULT: f32 = 3.55;

pub enum Charge {
    Positive,
    Negative,
}

lazy_static! {
    pub static ref pKa_table: BTreeMap<u8, (f32, Charge)> = {
        let mut m = BTreeMap::new();
        m.insert(b'K', (10.0, Charge::Positive));
        m.insert(b'R', (12.0, Charge::Positive));
        m.insert(b'H', (5.98, Charge::Positive));
        m.insert(b'D', (4.05, Charge::Negative));
        m.insert(b'E', (4.45, Charge::Negative));
        m.insert(b'C', (9.00, Charge::Negative));
        m.insert(b'Y', (10.0, Charge::Negative));
        m
    };
    pub static ref n_terminal_pKa_table: BTreeMap<u8, f32> = {
        let mut m = BTreeMap::new();
        m.insert(b'A', 7.59);
        m.insert(b'M', 7.00);
        m.insert(b'S', 6.93);
        m.insert(b'P', 8.36);
        m.insert(b'T', 6.82);
        m.insert(b'V', 7.44);
        m.insert(b'E', 7.70);
        m
    };
    pub static ref c_terminal_pKa_table: BTreeMap<u8, f32> = {
        let mut m = BTreeMap::new();
        m.insert(b'D', 4.55);
        m.insert(b'E', 4.75);
        m
    };
}

// --------------------------------

#[derive(Debug, Default)]
pub struct ProteinSeqAnalysis<'a> {
    pub seq: TextSlice<'a>,
    pub aa_count: AminoAcidCount,
    pub isoelectric_point: f32,
    pub molar_extinction_coefficient: (u32, u32),
}

impl<'a> ProteinSeqAnalysis<'a> {
    pub fn new(seq: TextSlice<'a>) -> Self {
        let mut res = Self::default();
        res.seq = seq;
        res.aa_count = Self::count_aa(seq);
        res
    }

    pub fn analyze(seq: TextSlice<'a>) -> Self {
        let mut res = Self::new(seq);
        res.isoelectric_point = res.calc_isoelectric_point();
        res.molar_extinction_coefficient = res.calc_molar_extinction_coefficient();
        res
    }

    pub fn count_aa(seq: TextSlice) -> AminoAcidCount {
        let mut res: BTreeMap<u8, u32> = BTreeMap::new();
        for &aa in seq {
            let count = res.entry(aa).or_insert(0);
            *count += 1;
        }
        res
    }

    /// Calculate the molar extinction coefficient.
    ///
    /// Calculates the molar extinction coefficient assuming cysteines (reduced) and cystines residues (oxidised)
    pub fn calc_molar_extinction_coefficient(&self) -> (u32, u32) {
        let mut mec_reduced = 0;
        if let Some(n) = self.aa_count.get(&b'W') {
            mec_reduced += 5500 * n;
        }
        if let Some(n) = self.aa_count.get(&b'Y') {
            mec_reduced += 1490 * n;
        }
        let mut mec_oxidised = mec_reduced;
        if let Some(n) = self.aa_count.get(&b'C') {
            mec_oxidised += 125 * n / 2;
        }
        (mec_reduced, mec_oxidised)
    }

    /// Estimate the isoelectric point of a polypeptide chain based on its primary structure.
    pub fn calc_isoelectric_point(&self) -> f32 {
        self.pi_recursive(4.05, 12.0, 7.775)
    }

    fn pi_recursive(&self, mut x1: f32, mut x2: f32, xmid: f32) -> f32 {
        if x2 - x1 < 0.0001 {
            return xmid;
        }
        let charge = self.charge_at_pH(xmid);
        if charge > 0.0 {
            x1 = xmid;
        } else {
            x2 = xmid;
        }
        self.pi_recursive(x1, x2, (x1 + x2) / 2.0)
    }
    pub fn charge_at_pH(&self, pH: f32) -> f32 {
        let mut charge = 0f32;
        for (aa, &count) in self.aa_count.iter() {
            if let Some((pKa, positivity)) = pKa_table.get(aa) {
                match positivity {
                    Charge::Negative => {
                        let partial_charge = 1.0 / (10f32.powf(pKa - pH) + 1.0);
                        charge -= partial_charge * count as f32;
                    }
                    Charge::Positive => {
                        let partial_charge = 1.0 / (10f32.powf(pH - pKa) + 1.0);
                        charge += partial_charge * count as f32;
                    }
                }
            }
        }
        let n_term_pKa = if let Some(&pKa) = n_terminal_pKa_table.get(&self.seq[0]) {
            pKa
        } else {
            N_TERM_PKA_DEFAULT
        };
        charge += 1.0 / (10f32.powf(pH - n_term_pKa) + 1.0);
        let c_term_pKa = if let Some(&pKa) = c_terminal_pKa_table.get(&self.seq[self.seq.len() - 1])
        {
            pKa
        } else {
            C_TERM_PKA_DEFAULT
        };
        charge -= 1.0 / (10f32.powf(c_term_pKa - pH) + 1.0);
        charge
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static S1: &[u8;152] = &b"MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV";

    lazy_static! {
        static ref RES1: ProteinSeqAnalysis<'static> = ProteinSeqAnalysis::analyze(S1);
    }

    #[test]
    fn test_isoelectric_point() {
        assert!((RES1.isoelectric_point - 7.72).abs() < 0.01)
    }
}
