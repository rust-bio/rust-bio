#![feature(test)]

extern crate test;

use bio::stats::pairhmm::*;
use bio::stats::{LogProb, Prob};
use test::Bencher;

static TEXT: &'static [u8] = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC\
ATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTAT\
CTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGT";

static PATTERN: &'static [u8] = b"GGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGC";

// Single base insertion and deletion rates for R1 according to Schirmer et al.
// BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

fn prob_emit_x_or_y() -> LogProb {
    LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
}

pub struct TestEmissionParams {
    x: &'static [u8],
    y: &'static [u8],
}

impl EmissionParameters for TestEmissionParams {
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
        if self.x[i] == self.y[j] {
            XYEmission::Match(LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST))
        } else {
            XYEmission::Mismatch(LogProb::from(PROB_ILLUMINA_SUBST / Prob(3.0)))
        }
    }

    fn prob_emit_x(&self, _: usize) -> LogProb {
        prob_emit_x_or_y()
    }

    fn prob_emit_y(&self, _: usize) -> LogProb {
        prob_emit_x_or_y()
    }

    fn len_x(&self) -> usize {
        self.x.len()
    }

    fn len_y(&self) -> usize {
        self.y.len()
    }
}

pub struct SemiglobalGapParams;

impl GapParameters for SemiglobalGapParams {
    fn prob_gap_x(&self) -> LogProb {
        LogProb::from(PROB_ILLUMINA_INS)
    }

    fn prob_gap_y(&self) -> LogProb {
        LogProb::from(PROB_ILLUMINA_DEL)
    }

    fn prob_gap_x_extend(&self) -> LogProb {
        LogProb::ln_zero()
    }

    fn prob_gap_y_extend(&self) -> LogProb {
        LogProb::ln_zero()
    }
}

pub struct SemiglobalAlignment;

impl StartEndGapParameters for SemiglobalAlignment {
    fn free_start_gap_x(&self) -> bool {
        true
    }

    fn free_end_gap_x(&self) -> bool {
        true
    }
}

#[bench]
fn pairhmm_semiglobal(b: &mut Bencher) {
    let emission_params = TestEmissionParams {
        x: TEXT,
        y: PATTERN,
    };
    let gap_params = SemiglobalGapParams;

    let mut pair_hmm = PairHMM::new(&gap_params);
    pair_hmm.prob_related(&emission_params, &SemiglobalAlignment, Some(4));

    b.iter(|| {
        let p = pair_hmm.prob_related(&emission_params, &SemiglobalAlignment, Some(4));
        assert!(*p <= 0.0);
    });
}
