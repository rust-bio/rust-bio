use crate::stats::LogProb;

pub mod evidence {
    /// Scale of evidence as defined by
    /// [Kass and Raftery 1995](http://www.andrew.cmu.edu/user/kk3n/simplicity/KassRaftery1995.pdf).
    #[derive(
        Display,
        Debug,
        Clone,
        Copy,
        PartialEq,
        Eq,
        PartialOrd,
        Ord,
        Serialize,
        Deserialize,
        EnumString,
        EnumIter,
        IntoStaticStr,
        VariantNames,
    )]
    pub enum KassRaftery {
        #[strum(serialize = "none")]
        None,
        #[strum(serialize = "barely")]
        Barely,
        #[strum(serialize = "positive")]
        Positive,
        #[strum(serialize = "strong")]
        Strong,
        #[strum(serialize = "very-strong")]
        VeryStrong,
    }
}

custom_derive! {
    /// A newtype for Bayes factors.
    #[derive(NewtypeFrom, NewtypeDeref, NewtypeAdd(*), NewtypeSub(*), Default, Copy, Clone, PartialEq, PartialOrd, Debug)]
    pub struct BayesFactor(pub f64);
}

impl BayesFactor {
    /// Calculate Bayes factor from given probabilities.
    pub fn new(a: LogProb, b: LogProb) -> Self {
        BayesFactor((a - b).exp())
    }

    /// Calculate strength of evidence as defined by
    /// [Kass and Raftery 1995](http://www.andrew.cmu.edu/user/kk3n/simplicity/KassRaftery1995.pdf).
    pub fn evidence_kass_raftery(&self) -> evidence::KassRaftery {
        let k = **self;
        if k <= 1.0 {
            evidence::KassRaftery::None
        } else if k <= 3.0 {
            evidence::KassRaftery::Barely
        } else if k <= 20.0 {
            evidence::KassRaftery::Positive
        } else if k <= 150.0 {
            evidence::KassRaftery::Strong
        } else {
            evidence::KassRaftery::VeryStrong
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bayes_factor() {
        let bf = BayesFactor::new(LogProb(0.5_f64.ln()), LogProb(0.1_f64.ln()));
        assert_relative_eq!(*bf, 5.0, epsilon = 1e-9);
        assert_eq!(bf.evidence_kass_raftery(), evidence::KassRaftery::Positive);
    }
}
