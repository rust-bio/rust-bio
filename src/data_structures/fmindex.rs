// Cloneright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! FM-Index and FMD-Index for finding suffix array intervals matching a given pattern in linear time.

use std::iter::DoubleEndedIterator;

use data_structures::suffix_array::SuffixArray;
use data_structures::bwt::{Occ, Less, less, BWT};
use alphabets::dna;
use std::fmt;
use std::mem::swap;
use std::ops::Deref;

/// A suffix array interval.
#[derive(Clone, Copy)]
pub struct Interval<SA: Deref<Target = SuffixArray> + Clone> {
    suffix_array: SA,
    lower: usize,
    upper: usize,
}

impl<SA: Deref<Target = SuffixArray> + Clone> Interval<SA> {
    pub fn occ(&self) -> Vec<usize> {
        self.suffix_array.range(self.lower..self.upper).expect("Interval should be in range of suffix array")
    }
}

impl<SA: Deref<Target = SuffixArray> + Clone> fmt::Debug for Interval<SA> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        fmt.debug_struct("Interval")
            .field("fmindex", &"hidden")
            .field("lower", &self.lower)
            .field("upper", &self.upper)
            .finish()
    }
}

pub trait FMIndexAble<SA: Deref<Target = SuffixArray> + Clone> {
    fn suffix_array(&self) -> SA;
    /// Get occurrence count of symbol a in BWT[..r+1].
    fn occ(&self, r: usize, a: u8) -> usize;
    /// Also known as
    fn less(&self, a: u8) -> usize;
    fn bwt(&self) -> &BWT;

    /// Perform backward search, yielding suffix array
    /// interval denoting exact occurences of the given pattern of length m in the text.
    /// Complexity: O(m).
    ///
    /// # Arguments
    ///
    /// * `pattern` - the pattern to search
    ///
    /// # Example
    ///
    /// ```
    /// use bio::data_structures::bwt::bwt;
    /// use bio::data_structures::fmindex::FMIndex;
    /// use bio::data_structures::suffix_array::suffix_array;
    /// use bio::alphabets::dna;
    ///
    /// let text = b"GCCTTAACATTATTACGCCTA$";
    /// let alphabet = dna::alphabet();
    /// let pos = suffix_array(text);
    /// let fm = FMIndex::new(bwt(text, &pos), 3, &alphabet);
    ///
    /// let pattern = b"TTA";
    /// let sai = fm.backward_search(pattern.iter());
    ///
    /// let occ = sai.occ(&pos);
    ///
    /// assert_eq!(occ, [3, 12, 9]);
    /// ```
    fn backward_search<'b, P: Iterator<Item = &'b u8> + DoubleEndedIterator> (&self, pattern: P) -> Interval<SA> where Self: Sized {
        let (mut l, mut r) = (0, self.bwt().len() - 1);
        for &a in pattern.rev() {
            let less = self.less(a);
            l = less +
                if l > 0 {
                self.occ(l - 1, a)
            } else {
                0
            };
            r = less + self.occ(r, a) - 1;
        }

        Interval {
            suffix_array: self.suffix_array(),
            lower: l,
            upper: r + 1,
        }
    }


}

/// The Fast Index in Minute space (FM-Index, Ferragina and Manzini, 2000) for finding suffix array
/// intervals matching a given pattern.

#[cfg_attr(feature = "serde_macros", derive(Serialize, Deserialize))]
pub struct FMIndex<
        SA: Deref<Target = SuffixArray> + Clone,
        DBWT: Deref<Target = BWT> + Clone,
        DLess: Deref<Target = Less> + Clone,
        DOcc: Deref<Target = Occ> + Clone> {
    sa: SA,
    bwt: DBWT,
    less: DLess,
    occ: DOcc,
}

impl<
    SA: Deref<Target = SuffixArray> + Clone,
    DBWT: Deref<Target = BWT> + Clone,
    DLess: Deref<Target = Less> + Clone,
    DOcc: Deref<Target = Occ> + Clone> FMIndexAble<SA> for FMIndex<SA, DBWT, DLess, DOcc> {

    fn suffix_array(&self) -> SA {
        self.sa.clone()
    }
    fn occ(&self, r: usize, a: u8) -> usize {
        self.occ.get(&self.bwt, r, a)
    }
    fn less(&self, a: u8) -> usize {
        self.less[a as usize]
    }
    /// Provide a reference to the underlying BWT.
    fn bwt(&self) -> &BWT {
        &self.bwt
    }
}

impl<
    SA: Deref<Target = SuffixArray> + Clone,
    DBWT: Deref<Target = BWT> + Clone,
    DLess: Deref<Target = Less> + Clone,
    DOcc: Deref<Target = Occ> + Clone> FMIndex<SA, DBWT, DLess, DOcc> {

    /// Construct a new instance of the FM index.
    ///
    /// # Arguments
    ///
    /// * `sa` - the suffix array (or sample)
    /// * `bwt` - the BWT
    /// * `k` - the sampling rate of the occ array: every k-th entry will be stored (higher k means
    ///   less memory usage, but worse performance)
    /// * `alphabet` - the alphabet of the underlying text, omitting the sentinel
    pub fn new(sa: SA, bwt: DBWT, less: DLess, occ: DOcc) -> Self {
        FMIndex {
            sa: sa,
            bwt: bwt,
            less: less,
            occ: occ,
        }
    }

    /// Construct a new instance of the FMD index (see Heng Li (2012) Bioinformatics).
    /// This expects a BWT that was created from a text over the DNA alphabet with N
    /// (`alphabets::dna::n_alphabet()`) consisting of the
    /// concatenation with its reverse complement, separated by the sentinel symbol `$`.
    /// I.e., let T be the original text and R be its reverse complement.
    /// Then, the expected text is T$R$. Further, multiple concatenated texts are allowed, e.g.
    /// T1$R1$T2$R2$T3$R3$.
    ///
    pub fn as_fmdindex(self) -> FMDIndex<SA, DBWT, DLess, DOcc> {
        let mut alphabet = dna::n_alphabet();
        alphabet.insert(b'$');
        assert!(alphabet.is_word(self.bwt()),
                "Expecting BWT over the DNA alphabet (including N) with the sentinel $.");

        FMDIndex {
            fmindex: self,
            revcomp: dna::RevComp::new(),
        }
    }
}

/// A bi-interval on suffix array of the forward and reverse strand of a DNA text.
#[derive(Clone, Copy)]
pub struct BiInterval<SA: Deref<Target = SuffixArray> + Clone> {
    suffix_array: SA,
    lower: usize,
    lower_rev: usize,
    size: usize,
    match_size: usize,
}

impl<SA: Deref<Target = SuffixArray> + Clone> fmt::Debug for BiInterval<SA> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        fmt.debug_struct("BiInterval")
            .field("fmindex", &"hidden")
            .field("lower", &self.lower)
            .field("lower_rev", &self.lower_rev)
            .field("size", &self.size)
            .field("match_size", &self.match_size)
            .finish()
    }
}

impl<SA: Deref<Target = SuffixArray> + Clone> BiInterval<SA> {
    pub fn forward(&self) -> Interval<SA> {
        Interval {
            suffix_array: self.suffix_array.clone(),
            upper: self.lower + self.size,
            lower: self.lower
        }
    }
    pub fn revcomp(&self) -> Interval<SA> {
        Interval {
            suffix_array: self.suffix_array.clone(),
            upper: self.lower_rev + self.size,
            lower: self.lower_rev
        }
    }

    fn swapped(&self) -> BiInterval<SA> {
        BiInterval {
            suffix_array: self.suffix_array.clone(),
            lower: self.lower_rev,
            lower_rev: self.lower,
            size: self.size,
            match_size: self.match_size,
        }
    }
}


/// The FMD-Index for linear time search of supermaximal exact matches on forward and reverse
/// strand of DNA texts (Li, 2012).
#[cfg_attr(feature = "serde_macros", derive(Serialize, Deserialize))]
pub struct FMDIndex<
    SA: Deref<Target = SuffixArray> + Clone,
    DBWT: Deref<Target = BWT> + Clone,
    DLess: Deref<Target = Less> + Clone,
    DOcc: Deref<Target = Occ> + Clone> {

    fmindex: FMIndex<SA, DBWT, DLess, DOcc>,
    revcomp: dna::RevComp,
}

impl<
    SA: Deref<Target = SuffixArray> + Clone,
    DBWT: Deref<Target = BWT> + Clone,
    DLess: Deref<Target = Less> + Clone,
    DOcc: Deref<Target = Occ> + Clone> FMIndexAble<SA> for FMDIndex<SA, DBWT, DLess, DOcc> {

    fn suffix_array(&self) -> SA {
        self.fmindex.suffix_array()
    }

    fn occ(&self, r: usize, a: u8) -> usize {
        self.fmindex.occ(r, a)
    }

    fn less(&self, a: u8) -> usize {
        self.fmindex.less(a)
    }

    /// Provide a reference to the underlying BWT.
    fn bwt(&self) -> &BWT {
        self.fmindex.bwt()
    }
}

impl<
    SA: Deref<Target = SuffixArray> + Clone,
    DBWT: Deref<Target = BWT> + Clone,
    DLess: Deref<Target = Less> + Clone,
    DOcc: Deref<Target = Occ> + Clone>  FMDIndex<SA, DBWT, DLess, DOcc> {

    /// Find supermaximal exact matches of given pattern that overlap position i in the pattern.
    /// Complexity O(m) with pattern of length m.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::data_structures::fmindex::FMDIndex;
    /// use bio::data_structures::suffix_array::suffix_array;
    /// use bio::data_structures::bwt::bwt;
    ///
    /// let text = b"ATTC$GAAT$";
    /// let pos = suffix_array(text);
    /// let fmdindex = FMDIndex::new(bwt(text, &pos), 3);
    ///
    /// let pattern = b"ATT";
    /// let intervals = fmdindex.smems(pattern, 2);
    /// let occ = intervals[0].occ(&pos);
    /// let occ_revcomp = intervals[0].occ_revcomp(&pos);
    ///
    /// assert_eq!(occ, [0]);
    /// assert_eq!(occ_revcomp, [6]);
    /// ```
    pub fn smems(&self, pattern: &[u8], i: usize) -> Vec<BiInterval<SA>> {

        let curr = &mut Vec::new();
        let prev = &mut Vec::new();
        let mut matches = Vec::new();

        let mut interval = self.init_interval(pattern, i);

        for &a in pattern[i + 1..].iter() {
            // forward extend interval
            let forward_interval = self.forward_ext(&interval, a);

            // if size changed, add last interval to list
            if interval.size != forward_interval.size {
                curr.push(interval.clone());
            }
            // if new interval size is zero, stop, as no further forward extension is possible
            if forward_interval.size == 0 {
                break;
            }
            interval = forward_interval;
        }
        // add the last non-zero interval
        curr.push(interval.clone());
        // reverse intervals such that longest comes first
        curr.reverse();

        swap(curr, prev);
        let mut j = pattern.len() as isize;

        for k in (-1..i as isize).rev() {
            let a = if k == -1 {
                b'$'
            } else {
                pattern[k as usize]
            };
            curr.clear();
            // size of the last confirmed interval
            let mut last_size = -1;

            for interval in prev.iter() {
                // backward extend interval
                let forward_interval = self.backward_ext(&interval, a);

                if (forward_interval.size == 0 || k == -1) &&
                        // interval could not be extended further
                        // if no interval has been extended this iteration,
                        // interval is maximal and can be added to the matches
                        curr.is_empty() && k < j {
                    j = k;
                    matches.push((*interval).clone());
                }
                // add _interval to curr (will be further extended next iteration)
                if forward_interval.size != 0 && forward_interval.size as isize != last_size {
                    last_size = forward_interval.size as isize;
                    curr.push(forward_interval);
                }
            }
            if curr.is_empty() {
                break;
            }
            swap(curr, prev);
        }

        matches
    }

    fn init_interval(&self, pattern: &[u8], i: usize) -> BiInterval<SA> {
        let a = pattern[i];
        let comp_a = self.revcomp.comp(a);
        let lower = self.fmindex.less(a);

        BiInterval {
            suffix_array: self.suffix_array(),
            lower: lower,
            lower_rev: self.fmindex.less(comp_a),
            size: self.fmindex.less(a + 1) - lower,
            match_size: 1,
        }
    }

    fn backward_ext(&self, interval: &BiInterval<SA>, a: u8) -> BiInterval<SA> {
        let mut s = 0;
        let mut o = 0;
        let mut l = interval.lower_rev;
        // Interval [l(c(aP)), u(c(aP))] is a subinterval of [l(c(P)), u(c(P))] for each a,
        // starting with the lexicographically smallest ($),
        // then c(T) = A, c(G) = C, c(C) = G, N, c(A) = T, ...
        // Hence, we calculate lower revcomp bounds by iterating over
        // symbols and updating from previous one.
        for &b in b"$TGCNAtgcna".iter() {
            l = l + s;
            o = self.fmindex.occ(interval.lower - 1, b);
            // calculate size
            s = self.fmindex.occ(interval.lower + interval.size - 1, b) - o;
            if b == a {
                break;
            }
        }
        // calculate lower bound
        let k = self.fmindex.less(a) + o;

        BiInterval {
            suffix_array: self.suffix_array(),
            lower: k,
            lower_rev: l,
            size: s,
            match_size: interval.match_size + 1,
        }
    }


    fn forward_ext(&self, interval: &BiInterval<SA>, a: u8) -> BiInterval<SA> {
        let comp_a = self.revcomp.comp(a);

        self.backward_ext(&interval.swapped(), comp_a)
            .swapped()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use alphabets::dna;
    use data_structures::suffix_array::{suffix_array, SuffixArray};
    use data_structures::bwt::{bwt, less, Occ};
    use std::rc::Rc;

    #[test]
    fn test_smems() {
        let revcomp = dna::RevComp::new();
        let orig_text = b"GCCTTAACAT";
        let revcomp_text = revcomp.get(orig_text);
        let text_builder: Vec<&[u8]> = vec![orig_text, b"$", &revcomp_text[..], b"$"];
        let text = text_builder.concat();

        let alphabet = dna::n_alphabet();
        let sa = Rc::new(suffix_array(&text));
        let bwt = Rc::new(bwt(&text, sa.as_ref()));
        let less = Rc::new(less(&bwt, &alphabet));
        let occ = Rc::new(Occ::new(&bwt, 3, &alphabet));

        let fmindex = FMIndex::new(sa as Rc<SuffixArray>, bwt, less, occ);
        let fmdindex = fmindex.as_fmdindex();
        {
            let pattern = b"AA";
            let intervals = fmdindex.smems(pattern, 0);
            assert_eq!(intervals[0].forward().occ(), [5, 16]);
            assert_eq!(intervals[0].revcomp().occ(), [3, 14]);
        }
        {
            let pattern = b"CTTAA";
            let intervals = fmdindex.smems(pattern, 1);
            assert_eq!(intervals[0].forward().occ(), [2]);
            assert_eq!(intervals[0].revcomp().occ(), [14]);
            assert_eq!(intervals[0].match_size, 5)
        }
    }


    #[test]
    fn test_init_interval() {
        let text = b"ACGT$TGCA$";

        let alphabet = dna::n_alphabet();
        let sa = Rc::new(suffix_array(text));
        let bwt = Rc::new(bwt(text, sa.as_ref()));
        let less = Rc::new(less(&bwt, &alphabet));
        let occ = Rc::new(Occ::new(&bwt, 3, &alphabet));

        let fmindex = FMIndex::new(sa as Rc<SuffixArray>, bwt, less, occ);
        let fmdindex = fmindex.as_fmdindex();
        let pattern = b"T";
        let interval = fmdindex.init_interval(pattern, 0);
        assert_eq!(interval.forward().occ(), [3, 5]);
        assert_eq!(interval.revcomp().occ(), [8, 0]);
    }

    #[test]
    #[cfg(feature = "nightly")]
    fn test_serde() {
        use serde::{Serialize, Deserialize};
        fn impls_serde_traits<S: Serialize + Deserialize>() {}

        impls_serde_traits::<FMIndex>();
        impls_serde_traits::<FMDIndex>();
    }

    #[test]
    fn test_issue39() {
        let reads = b"GGCGTGGTGGCTTATGCCTGTAATCCCAGCACTTTGGGAGGTCGAAGTGGGCGG$CCGC\
                       CCACTTCGACCTCCCAAAGTGCTGGGATTACAGGCATAAGCCACCACGCC$CGAAGTGG\
                       GCGGATCACTTGAGGTCAGGAGTTGGAGACTAGCCTGGCCAACACGATGAAACCCCGTC\
                       TCTAATA$TATTAGAGACGGGGTTTCATCGTGTTGGCCAGGCTAGTCTCCAACTCCTGA\
                       CCTCAAGTGATCCGCCCACTTCG$AGCTCGAAAAATGTTTGCTTATTTTGGTAAAATTA\
                       TTCATTGACTATGCTCAGAAATCAAGCAAACTGTCCATATTTCATTTTTTG$CAAAAAA\
                       TGAAATATGGACAGTTTGCTTGATTTCTGAGCATAGTCAATGAATAATTTTACCAAAAT\
                       AAGCAAACATTTTTCGAGCT$AGCTCGAAAAATGTTTGCTTATTTTGGTAAAATTATTC\
                       ATTGACTATGCTCAGAAATCAAGCAAACTGTCCATATTTCATTTTTTGAAATTACATAT\
                       $ATATGTAATTTCAAAAAATGAAATATGGACAGTTTGCTTGATTTCTGAGCATAGTCAA\
                       TGAATAATTTTACCAAAATAAGCAAACATTTTTCGAGCT$TAAAATTTCCTCTGACAGT\
                       GTAAAAGAGATCTTCATACAAAAATCAGAATTTATATAGTCTCTTTCCAAAAGACCATA\
                       AAACCAATCAGTTAATAGTTGAT$ATCAACTATTAACTGATTGGTTTTATGGTCTTTTG\
                       GAAAGAGACTATATAAATTCTGATTTTTGTATGAAGATCTCTTTTACACTGTCAGAGGA\
                       AATTTTA$CACCTATCTACCCTGAATCTAAGTGCTAACAGGAAAGGATGCCAGATTGCA\
                       TGCCTGCTGATAAAGCCACAGTTTGGACTGTCACTCAATCACCATCGTTC$GAACGATG\
                       GTGATTGAGTGACAGTCCAAACTGTGGCTTTATCAGCAGGCATGCAATCTGGCATCCTT\
                       TCCTGTTAGCACTTAGATTCAGGGTAGATAGGTG$CATCGTTCCTCCTGTGACTCAGTA\
                       TAACAAGATTGGGAGAATACTCTACAGTTCCTGATTCCCCCACAG$CTGTGGGGGAATC\
                       AGGAACTGTAGAGTATTCTCCCAATCTTGTTATACTGAGTCACAGGAGGAACGATG$TG\
                       TAAATTCTGAGAAAAATTTGCAGGTCTTTCTTCAGGAGCATGTAATCTCTTGCTCTCTT\
                       TGTTATCTATCTATAGTACTGTAGGTTATCTGGAGTTGCT$AGCAACTCCAGATAACCT\
                       ACAGTACTATAGATAGATAACAAAGAGAGCAAGAGATTACATGCTCCTGAAGAAAGACC\
                       TGCAAATTTTTCTCAGAATTTACA$CACTTCTCCTTGTCTTTACAGACTGGTTTTGCAC\
                       TGGGAAATCCTTTCACCAGTCAGCCCAGTTAGAGATTCTG$CAGAATCTCTAACTGGGC\
                       TGACTGGTGAAAGGATTTCCCAGTGCAAAACCAGTCTGTAAAGACAAGGAGAAGTG$AA\
                       TGGAGGTATATAAATTATCTGGCAAAGTGACATATCCTGACACATTCTCCAGGATAGAT\
                       CAAATGTTAGGTCACAAAGAGAGTCTTAACAAAATT$AATTTTGTTAAGACTCTCTTTG\
                       TGACCTAACATTTGATCTATCCTGGAGAATGTGTCAGGATATGTCACTTTGCCAGATAA\
                       TTTATATACCTCCATT$TTAATTTTGTTAAGACTCTCTTTGTGACCTAACATTTGATCT\
                       ATCCTGGAGAATGTGTCAGGATATGTCACTTTGCCAGATAATTTATATACCTCCATTTT\
                       $AAAATGGAGGTATATAAATTATCTGGCAAAGTGACATATCCTGACACATTCTCCAGGA\
                       TAGATCAAATGTTAGGTCACAAAGAGAGTCTTAACAAAATTAA$TTCTTCTTTGACTCA\
                       TTGGTTGTTCAATAGTATGTTGTTTAATTTCCATATATTTGTAAATGTTTCCGTTTTCC\
                       TTCTACTATTGAATTTTTGCTTCATC$GATGAAGCAAAAATTCAATAGTAGAAGGAAAA\
                       CGGAAACATTTACAAATATATGGAAATTAAACAACATACTATTGAACAACCAATGAGTC\
                       AAAGAAGAA$AGGAAAACGGAAACATTTACAAATATATGGAAATTAAACAACATACTAT\
                       TGAACAACCAATGAGTCAAAGAAGAAATCAAAAAGAATATTAGAAAAC$GTTTTCTAAT\
                       ATTCTTTTTGATTTCTTCTTTGACTCATTGGTTGTTCAATAGTATGTTGTTTAATTTCC\
                       ATATATTTGTAAATGTTTCCGTTTTCCT$TTAGAAAACAAGCTGACAAAAAAATAAAAA\
                       AACACAACATAGCAAAACTTAGAAATGCAGCAAAGGCAGTACTAAAGAGGGAAATTTAT\
                       AGCAATAAATGC$GCATTTATTGCTATAAATTTCCCTCTTTAGTACTGCCTTTGCTGCA\
                       TTTCTAAGTTTTGCTATGTTGTGTTTTTTTATTTTTTTGTCAGCTTGTTTTCTAA$TTT\
                       ATTGCTATAAATTTCCCTCTTTAGTACTGCCTTTGCTGCATTTCTAAGTTTTGCTATGT\
                       TGTGTTTTTTTATTTTTTTGTCAGCTTGTTTTCTA$TAGAAAACAAGCTGACAAAAAAA\
                       TAAAAAAACACAACATAGCAAAACTTAGAAATGCAGCAAAGGCAGTACTAAAGAGGGAA\
                       ATTTATAGCAATAAA$TCTTTCTTCTTTTTTAAGGTAGGCATTTATTGCTATAAATTTC\
                       CCTCTTTAGTACTGCCTTTG$CAAAGGCAGTACTAAAGAGGGAAATTTATAGCAATAAA\
                       TGCCTACCTTAAAAAAGAAGAAAGA$";

        let alphabet = dna::n_alphabet();
        let sa = Rc::new(suffix_array(reads));
        let bwt = Rc::new(bwt(reads, sa.as_ref()));
        let less = Rc::new(less(&bwt, &alphabet));
        let occ = Rc::new(Occ::new(&bwt, 3, &alphabet));

        let fmindex = FMIndex::new(sa as Rc<SuffixArray>, bwt, less, occ);
        let fmdindex = fmindex.as_fmdindex();

        let read = b"GGCGTGGTGGCTTATGCCTGTAATCCCAGCACTTTGGGAGGTCGAAGTGGGCGG";
        let read_pos = 0;

        for i in 0..read.len() {
            println!("i {}", i);
            let intervals = fmdindex.smems(read, i);
            println!("{:?}", intervals);
            let matches = intervals.iter()
                                   .flat_map(|interval| interval.forward().occ() )
                                   .collect::<Vec<usize>>();
            assert_eq!(matches, vec![read_pos]);
        }
    }
}
