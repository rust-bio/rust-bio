use stats::LogProb;

pub mod evidence {
    use std::str::FromStr;

    /// Scale of evidence as defined by
    /// [Kass and Raftery 1995](http://www.andrew.cmu.edu/user/kk3n/simplicity/KassRaftery1995.pdf).
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
    pub enum KassRaftery {
        None,
        Barely,
        Positive,
        Strong,
        VeryStrong,
    }

    impl FromStr for KassRaftery {
        type Err = ParseError;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s {
                "none" => Ok(KassRaftery::None),
                "barely" => Ok(KassRaftery::Barely),
                "positive" => Ok(KassRaftery::Positive),
                "strong" => Ok(KassRaftery::Strong),
                "very-strong" => Ok(KassRaftery::VeryStrong),
                _ => Err(ParseError::Some),
            }
        }
    }

    quick_error! {
        #[derive(Debug)]
        pub enum ParseError {
            Some {
                description("unable to parse Kass Raftery score")
                display("valid values are none, barely, positive, strong, very-strong")
            }
        }
    }
}

custom_derive! {
    /// A newtype for Bayes factors.
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        PartialEq,
        PartialOrd,
        Copy,
        Clone,
        Debug,
    )]
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
        if k >= 1.0 && k <= 3.0 {
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
