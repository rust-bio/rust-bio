// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Banded Smith-Waterman alignment for fast comparison of long strings.
//! Use sparse dynamic programming to find a 'backbone' alignment from exact
//! k-mer matches, then compute the SW alignment in a 'band' surrounding the
//! backbone, with a configurable width w. This method is not guaranteed
//! to recover the Smith-Waterman alignment, but will usually find the same 
//! alignment if a) there is a reasonable density of exact k-mer matches
//! between the sequences, and b) the width parameter w is larger than the
//! excursion of the alignment path from diagonal between successive kmer 
//! matches.  This technique is employed in long-read aligners (e.g. BLASR and BWA)
//! to drastically reduce runtime compared to Smith Waterman. Currently only local
//! alignment is implemented.
//! Complexity roughly O(min(m,n) * w)
//!
//! # Example
//!
//! ```
//! use bio::alignment::pairwise::banded::*;
//!
//! let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
//! let y = b"AGCACACGTGTGCGCTATACAGTACACGTGTCACAGTTGTACTAGCATGAC";
//! let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
//! let k = 8;
//! let w = 6;
//! let mut aligner = Aligner::new(-5, -1, &score, k, w);
//! let alignment = aligner.local(x, y);
//! assert_eq!(alignment.ystart, 0);
//! assert_eq!(alignment.xstart, 0);
//! ```

use std::i32;
use std::iter::repeat;

use std::collections::HashMap;
use alignment::{Alignment, AlignmentOperation};
use utils::TextSlice;
use std::cmp::min;
use std::cmp::max;
use std::ops::Range;

use super::*;
use alignment::sparse;
use std::rc::Rc;


macro_rules! align_banded {
    (
        $aligner:ident, $x:ident, $y:ident, $state:ident,
        $init:block, $inner:block, $outer:block, $ret:block
    ) => (
        {
            let mut $state = AlignmentState {
                m: $x.len(), n: $y.len(),
                best: 0, best_i: 0, best_j: 0,
                best_layer: TBSUBST,
                i: 1, j: 1,
                score: 0, i_score: 0, d_score: 0,
                col: 0
            };

            let ref band = $aligner.band;

            while $state.i <= $state.n {
                $state.col = $state.i % 2;
                let prev = 1 - $state.col;

                // init code
                $init

                // read next y symbol
                let b = $y[$state.i - 1];
                $state.j = max(1, band.ranges[$state.i].start);

                for jj in max(1, band.ranges[$state.i].start.saturating_sub(1))..$state.j {
                    $aligner.I[$state.col][jj] = MIN_SCORE;
                    $aligner.D[$state.col][jj] = MIN_SCORE;
                    $aligner.S[$state.col][jj] = MIN_SCORE;
                } 

                while $state.j <= min($state.m, band.ranges[$state.i].end) {
                    // read next x symbol
                    let a = $x[$state.j - 1];

                    let mut tb = TracebackCell::new();

                    // score for deletion
                    let d_open = $aligner.S[prev][$state.j] + $aligner.gap_open;
                    let d_extend = $aligner.D[prev][$state.j] + $aligner.gap_extend;

                    if d_open > d_extend {
                        tb.set_d(TBSUBST);
                        $aligner.D[$state.col][$state.j] = d_open;
                    } else {
                        tb.set_d(TBDEL);
                        $aligner.D[$state.col][$state.j] = d_extend;
                    }

                    // score for insertion
                    let i_open = $aligner.S[$state.col][$state.j-1] + $aligner.gap_open;
                    let i_extend = $aligner.I[$state.col][$state.j-1] + $aligner.gap_extend;

                    if i_open > i_extend {
                        tb.set_i(TBSUBST);
                        $aligner.I[$state.col][$state.j] = i_open;
                    } else {
                        tb.set_i(TBINS);
                        $aligner.I[$state.col][$state.j] = i_extend;
                    };


                    // score for substitution
                    let match_score = ($aligner.score)(a, b);
                    $state.score = $aligner.S[prev][$state.j-1] + match_score;
                    tb.set_s(TBSUBST);

                    let from_d = $aligner.D[prev][$state.j-1] + match_score;
                    let from_i = $aligner.I[prev][$state.j-1] + match_score;

                    if from_d > $state.score {
                        $state.score = from_d;
                        tb.set_s(TBDEL);
                    }

                    if from_i > $state.score {
                        $state.score = from_i;
                        tb.set_s(TBINS);
                    }

                    $aligner.traceback.set($state.i, $state.j, tb);

                    // inner code
                    $inner

                    $aligner.S[$state.col][$state.j] = $state.score;
                    $state.j += 1;
                }

                for jj in $state.j..min($state.m, band.ranges[min($state.i+1, $state.n)].end) {
                    $aligner.I[$state.col][jj] = MIN_SCORE;
                    $aligner.D[$state.col][jj] = MIN_SCORE;
                    $aligner.S[$state.col][jj] = MIN_SCORE;
                }

                // outer code
                $outer

                $state.i += 1;
            }

            // return code
            $ret
        }
    );
}


/// A generalized Smith-Waterman aligner.
#[allow(non_snake_case)]
pub struct Aligner<'a, F>
    where F: 'a + Fn(u8, u8) -> i32
{
    S: [Vec<i32>; 2],
    I: [Vec<i32>; 2],
    D: [Vec<i32>; 2],
    traceback: Traceback,
    gap_open: i32,
    gap_extend: i32,
    score: &'a F,

    band: Rc<Band>,
    k: usize,
    w: usize,
}


const DEFAULT_ALIGNER_CAPACITY: usize = 200;

impl<'a, F> Aligner<'a, F>
    where F: Fn(u8, u8) -> i32
{
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `score` - function that returns the score for substitutions (also see bio::scores)
    ///
    pub fn new(gap_open: i32, gap_extend: i32, score: &'a F, k: usize, w: usize) -> Self {
        Aligner::with_capacity(DEFAULT_ALIGNER_CAPACITY,
                               DEFAULT_ALIGNER_CAPACITY,
                               gap_open,
                               gap_extend,
                               score, k, w)
    }

    /// Create new aligner instance. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `score` - function that returns the score for substitutions (also see bio::scores)
    ///
    pub fn with_capacity(m: usize, n: usize, gap_open: i32, gap_extend: i32, score: &'a F, k: usize, w: usize) -> Self {
        let get_vec = || Vec::with_capacity(m + 1);

        let band = Rc::new(Band { ranges: vec![(0..1)] });

        let al = Aligner {
            band: band.clone(),
            S: [get_vec(), get_vec()],
            I: [get_vec(), get_vec()],
            D: [get_vec(), get_vec()],
            traceback: Traceback::with_capacity(m,n),
            gap_open: gap_open,
            gap_extend: gap_extend,
            score: score,
            k: k,
            w: w,
        };

        al
    }

    /// Create new aligner instance with unit (equal to '-1') penalties for gap open and gap extend
    /// and unit score function ('1' if two letters are equal, '-1' if not). This is
    /// effectively equal to Levenshtein metric.
    // pub fn with_unit_cost() -> Self {
    //     let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    //     Aligner::new(-1, -1, *&score)
    // }

    #[inline(never)]
    fn init(&mut self, x: TextSlice, y: TextSlice, alignment_type: AlignmentType) {
        let m = x.len();
        let n = y.len();

        let k = min(self.k, min(m, n));
        // Update the band.
        {
            let _band = Band::create_local(x,y, k, self.w, self.gap_open, self.gap_extend);
            let band = Rc::new(_band);
            self.band = band.clone();
            self.traceback = Traceback::with_capacity(m, n);
        }

        self.traceback.init(m, n, alignment_type);

        // set minimum score to -inf, and allow to add gap_extend
        // without overflow
        for k in 0..2 {
            self.S[k].clear();
            self.I[k].clear();
            self.D[k].clear();

            self.D[k].extend(repeat(MIN_SCORE).take(m + 1));

            match alignment_type {
                AlignmentType::Semiglobal |
                AlignmentType::Global => {
                    let mut i = &mut self.I[k];

                    // need one insertion to establish a gap
                    i.push(MIN_SCORE);

                    // other cells are gaps
                    let mut score = self.gap_open;
                    for _ in 1..m + 1 {
                        i.push(score);
                        score += self.gap_extend;
                    }

                    self.S[k].push(0);
                    // Impossible to reach S state after first position in first column
                    self.S[k].extend(repeat(MIN_SCORE).take(m));
                },

                AlignmentType::Local => {
                    self.S[k].extend(repeat(0).take(m + 1));
                    self.I[k].extend(repeat(MIN_SCORE).take(m + 1));
                },
            }
        }
    }


    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        self.init(x, y, AlignmentType::Semiglobal);

        align_banded!(self,
               x,
               y,
               state,
               {
                   self.S[state.col][0] = 0;
               },
               {},
               {
                   // the second condition ensures that score is overwritten if best
                   // does not reflect a full x-column (can happen in first iteration)
                   if state.score > state.best || state.best_j != state.m {
                       state.best = state.score;
                       state.best_i = state.i;
                       state.best_j = state.m;
                   }
               },
               {
                   self.alignment(state.best_i, state.best_j, x, y, state.best)
               })
    }


    /// Calculate local alignment of x against y.
    pub fn local(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        self.init(x, y, AlignmentType::Local);

        align_banded!(self,
               x,
               y,
               state,
               {
                   self.S[state.col][0] = 0;
               },
               {
                   if state.score < 0 {
                       self.traceback.get_mut(state.i, state.j).set_s(TBSTART);
                       state.score = 0;
                   } else if state.score > state.best {
                       state.best = state.score;
                       state.best_i = state.i;
                       state.best_j = state.j;
                   }
               },
               {},
               {
                   self.alignment(state.best_i, state.best_j, x, y, state.best)
               })
    }

    fn alignment(&self, yend: usize, xend: usize, x: TextSlice, y: TextSlice, score: i32) -> Alignment {

        let mut i = yend;
        let mut j = xend;
        //self.print_traceback_matrices(i,j);

        let mut ops = Vec::with_capacity(x.len());

        let get = move |i,j,ty| {
                match ty {
                    TBDEL => self.traceback.get(i,j).get_d(),
                    TBINS => self.traceback.get(i,j).get_i(),
                    _ => self.traceback.get(i,j).get_s(),
                }
        };

        let mut which_mat = TBSUBST;

        loop {
            let tb = get(i, j, which_mat);

            if tb == TBSTART {
                break;
            }

            let (ii, jj, op) = match which_mat {
                TBSUBST => {
                    let op = if y[i - 1] == x[j - 1] {
                        AlignmentOperation::Match
                    } else {
                        AlignmentOperation::Subst
                    };
                    (i - 1, j - 1, op)
                }
                TBDEL => (i - 1, j, AlignmentOperation::Del),
                TBINS => (i, j - 1, AlignmentOperation::Ins),
                _ => {
                    break;
                }
            };

            ops.push(op);
            i = ii;
            j = jj;
            which_mat = tb;
        }

        ops.reverse();
        Alignment {
            ystart: i,
            xstart: j,
            yend: yend,
            xend: xend,
            xlen: x.len(),
            operations: ops,
            score: score,
        }
    }

    // Debugging helper function for visualizing traceback matrices
    #[allow(dead_code)]
    fn print_traceback_matrices(&self, i: usize, j: usize)
    {
        println!("--");
        for tb in &[TBSUBST, TBINS, TBDEL] {
            println!("--");
            for jj in 0..(j+1) {
                let mut s = String::new();
                for ii in 0..(i+1) {
                    match self.traceback.get(ii,jj).get(*tb) {
                        TBSUBST => s.push_str(" M"),
                        TBDEL => s.push_str(" D"),
                        TBINS => s.push_str(" I"),
                        TBSTART => s.push_str(" S"),
                        _ => (),
                    }
                }
                println!("{}", s);
            }
        }
    }
}

#[derive(Clone, Debug)]
struct Band {
    ranges: Vec<Range<usize>>
}

pub fn key_line<I: Iterator<Item=(u32, u32)>>(mut kmer_starts: I, _k: usize) -> Vec<(u32,u32)> {
    let k = _k as u32;
    let mut line = Vec::new();

    let mut state = 
        match kmer_starts.next() {
            Some(s) => {
                line.push(s);
                s
            },
            None => {
                return line;
            },
        };

    loop {
        let next = kmer_starts.next();

        if next.is_some() {
            let next = next.unwrap();
            if state.0 + 1 == next.0 && state.1 + 1 == next.1 {
                state = next;
            } else {
                line.push((state.0 + k, state.1 + k));
                line.push(next);
                state = next;
            }
        } else {
            break;
        }
    }

    line.push((state.0 + k, state.1 + k));
    line
}

impl Band {


    pub fn find_kmer_matches<T: AsRef<[u8]>>(seq1: &T, seq2: &T, k: usize) -> Vec<(u32, u32)> {

        let slc1 = seq1.as_ref();
        let slc2 = seq2.as_ref();

        let mut set: HashMap<&[u8], Vec<u32>> = HashMap::new();
        let mut matches = Vec::new();

        for i in 0 .. slc1.len() - k + 1 {
            set.entry(&slc1[i..i+k]).or_insert_with(|| Vec::new()).push(i as u32);
        }

        for i in 0 .. slc2.len() - k + 1 {
            let slc = &slc2[i..i+k];
            match set.get(slc) {
                Some(matches1) => {
                    for pos1 in matches1 {
                        matches.push((*pos1, i as u32));
                    }
                },
                None => (),
            }
        }

        matches.sort();
        matches
    }

    pub fn create_local(x: TextSlice, y: TextSlice, k: usize, w: usize, gap_open: i32, gap_extend: i32) -> Band {
        let matches = Band::find_kmer_matches(&x,&y,k);
        let res = sparse::sdpkpp(&matches, k, 2, gap_open, gap_extend);
        let ps = res.path[0];
        let pe = res.path[res.path.len()-1];
        println!("sparse: rstart:{} tstart:{} rend:{}, tend:{}, hits:{}", matches[ps].0, matches[ps].1, matches[pe].0, matches[pe].1, res.score);

        // each entry in matches that is included in the returned path, generates a diagonal line k bases long
        // each gap generates a line from the end of the previous k-match to the start of the next one.
        // we want the band to be a dilated version of this line, such it covers and position within p bases of 
        // this line in x or y.


        //let line = key_line(path.map(|idx| matches[idx]));
        let mut _prev : Option<(usize, usize)> = None;

        let mut ranges: Vec<Range<usize>> = Vec::with_capacity(y.len()+1);
        for _ in 0 .. y.len()+1 {
            ranges.push(x.len()+1..0);
        }

        let ranges_len = ranges.len();

        { // borrow window of ranges 
        let mut reg = |r: usize, c:usize| {

            let c1 = min(c + w, ranges_len - 1);
            ranges[c].start = min(ranges[c].start, r.saturating_sub(w));
            ranges[c1].start = min(ranges[c1].start, r.saturating_sub(w));

            let c2 = c.saturating_sub(w);
            ranges[c].end = max(ranges[c].end, min(r + w, x.len() + 1));
            ranges[c2].end = max(ranges[c2].end, min(r + w, x.len() + 1));
        };

        for idx in res.path {
            let _m = matches[idx];
            let m = (_m.0 as usize, _m.1 as usize);

            match _prev {
                // No previous match -- fill diagonal leading to first hit
                None => {
                    for c in 0..m.1 {
                        let r = m.0 - min(m.0, (m.1 - c));
                        reg(r, c);
                    }

                    // band rest of kmer
                    for i in 0 .. k {
                        reg(m.0 + i, m.1 + i)
                    }
                },

                // Prev match exists -- fill diagonal
                Some(prev) => {
                    // Check if prev kmer overlaps w/ current kmer
                    // case 1: they do: must be successive kmers, so do 
                    //         fast-path of single column at end of 2nd kmer
                    // case 2: connect along line between end of prev and
                    //         start of current, then along the current kmer

                    // case 1:
                    if prev.0 == m.0 - 1 && prev.1 == m.1 - 1 {
                        reg(m.0 + k - 1, m.1 + k - 1);

                    // case 2:
                    } else {
                        for c in prev.1 .. m.1 {
                            let r = prev.0 + (m.0 - prev.0) * (c - prev.1) / (m.1 - prev.1);
                            reg(r, c)
                        }
                        
                        // band rest of kmer
                        for i in 0 .. k {
                            reg(m.0 + i, m.1 + i)
                        }
                    }
                }
            }

            _prev = Some(m);
        }

        // Fill out diagonal along final match to end of matrix
        match _prev {
            Some(last) => {
                for c in last.1..ranges_len {
                    let r = last.0 + (c - last.1);
                    reg(r, c);
                }
            },
            None => (),
        }

        } // end borrow of ranges

        // Don't allow negative size ranges
        for i in 0 .. y.len()+1 {
            ranges[i].start = min(ranges[i].start, ranges[i].end);
        }

        Band { ranges: ranges }
    }
}


#[cfg(test)]
mod banded {
    use alignment::pairwise::{self, banded};
    use utils::TextSlice;

    // Check that the banded alignment is equivalent to the exhaustive SW alignment
    fn compare_to_full_alignment(x: TextSlice, y: TextSlice) {

        let score = |a: u8, b: u8| {
            if a == b {
                1i32
            } else {
                -1i32
            }
        };

        
        let mut banded_aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let banded_alignment = banded_aligner.local(x, y);

        let mut full_aligner = pairwise::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let full_alignment = full_aligner.local(x, y);

        assert_eq!(banded_alignment, full_alignment);
    }

    #[test]
    fn test_same() {
        let x = b"ACGTATCATAGACCCTAGATAGGGTTGTGTAGATGATCCACAGACGTATCATAGATTAGATAGGGTTGTGTAGATGATTCCACAG";
        let y = x.clone();
        compare_to_full_alignment(x,y);
    }

    /*
    #[test]
    fn test_speed() {
        let x = b"ACGTATCGATAGCTGAGTTACCAAGATAGGGTTGTGTAGATGAAGAGATAGCTGAGTTACCAGTTCACGTGCAACTAGACCCTAGATAGGGTTGTGTAGATGATCCACAGACGTATCATAGATTATCAAGAGATAGCTGAGTTACCAAGATAGGGTTGTGTAGATGATTCCACAG";
        let y = x.clone();

        for _ in 0..1000000 {
            compare_to_full_alignment(x,y);
        }
    }
    */

    #[test]
    fn test_big() {
        let query = b"CATCTCCACCCACCCTATCCAACCCTGGGGTGGCAGGTCGTGAGTGACAGCCCCAAGGACACCAAGGGATGAAGCTTCTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGCAGAACAGGAGGCCCCATCTCCCGTTCTTACTCAGAAGCTGTCAGCAGGGCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTCCCCAGCTCTTGCTTCAGGGCCTCACGTTTCGCCCCCTGAGGCCTGGGGGCTCCATCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTG";

        let target = b"CCTCCCATCTCCACCCACCCTATCCAACCCTGGGGTGGCAGGTCATGAGTGACAGCCCCAAGGACACCAAGGGATGAAGCTTCTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGCAGAACAGGAGGCCCCATCTCCCGTTCTTACTCAGAAGCTGTCAGCAGGGCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTCCCCAGCTCTTGCTTCAGGGCCTCACGTTTCGCCCCCTGAGGCCTGGGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGCACGGCTCCCAACTCTCTTCCGGCCAAGGATCCCGTGTTCCTGAAATGTCTTTCTACCAAACACAGTTGCTGTGTAACCACTCATTTCATTTTCCTAATTTGTGTTGATCCAGGACACGGGAGGAGACCTGGGCAGCGGCGGACTCATTGCAGGTCGCTCTGCGGTGAGGACGCCACAGGCAC";

        compare_to_full_alignment(query, target);
    }

    #[test]
    fn test_deletion() {
        let x = b"AGCACACGTGTGCGCTATACAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACACGTGTGCGCTATACAGTAAAAAAAACACGTGTCACAGTTGTACTAGCATGAC";
        compare_to_full_alignment(x,y);
    }

    #[test]
    fn test_insertion() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACACGTGTGCGCTATACAGTACACGTGTCACAGTTGTACTAGCATGAC";
        compare_to_full_alignment(x,y);
    }

    #[test]
    fn test_substitutions() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment(x,y);
    }

    #[test]
    fn test_overhangs1() {
        let x =             b"CGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment(x,y);
    }

    #[test]
    fn test_overhangs2() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y =                b"TATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment(x,y);
    }

    #[test]
    fn test_overhangs3() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTG";
        let y = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment(x,y);
    }

    #[test]
    fn test_overhangs4() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCA";
        compare_to_full_alignment(x,y);
    }
}
