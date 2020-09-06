//! Estimate the isoelectric point of a polypeptide chain based on its primary structure.

use crate::seq_analysis::protein::count_aa;
use crate::utils::TextSlice;
use lazy_static::lazy_static;
use std::collections::BTreeMap;
const N_TERM_PKA_DEFAULT: f32 = 7.5; // positive
const C_TERM_PKA_DEFAULT: f32 = 3.55; // negative

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

pub fn isoelectric_point(seq: TextSlice) -> f32 {
    pi_recursive(seq, 4.05, 12.0, 7.775)
}

fn pi_recursive(seq: TextSlice, mut x1: f32, mut x2: f32, xmid: f32) -> f32 {
    if x2 - x1 < 0.0001 {
        return xmid;
    }
    let charge = charge_at_pH(seq, xmid);
    if charge > 0.0 {
        x1 = xmid;
    } else {
        x2 = xmid;
    }
    pi_recursive(seq, x1, x2, (x1 + x2) / 2.0)
}

fn charge_at_pH(seq: TextSlice, pH: f32) -> f32 {
    let mut aa_count = count_aa(seq);
    let mut charge = 0f32;
    for (aa, &count) in aa_count.iter() {
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
    let n_term_pKa = if let Some(&pKa) = n_terminal_pKa_table.get(&seq[0]) {
        pKa
    } else {
        N_TERM_PKA_DEFAULT
    };
    charge += 1.0 / (10f32.powf(pH - n_term_pKa) + 1.0);
    let c_term_pKa = if let Some(&pKa) = c_terminal_pKa_table.get(&seq[seq.len() - 1]) {
        pKa
    } else {
        C_TERM_PKA_DEFAULT
    };
    charge -= 1.0 / (10f32.powf(c_term_pKa - pH) + 1.0);
    charge
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isoelectric_point() {
        let seq = b"MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV";
        let pi = isoelectric_point(seq);
        assert!(pi - 7.72 < 0.01)
    }
}
