// Copyright 2017 Nicholas Ingolia
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Positions on a named sequence, e.g., 123,456 on chromosome 3.

use std::fmt;
use std::str::FromStr;

use io::bed;
use seq_loc::{break_fwd_rev,break_refid,ParseError,ReqStrand};

/// Position on a named sequence (without strand information)
///
/// The display format for a `SeqPos` is _chr:pos_.
///
/// ```
/// # use bio::seq_loc::ParseError;
/// # fn try_main() -> Result<(), Box<ParseError>> {
/// use bio::seq_loc::seq_pos::SeqPos;
/// let start = SeqPos::new("chrIV".to_owned(), 683946);
/// let start_str = start.to_string();
/// assert_eq!(start_str, "chrIV:683946");
/// let start_str_pos: SeqPos = start_str.parse()?;
/// assert_eq!(start, start_str_pos);
/// # Ok(()) 
/// # }
/// # fn main() { try_main().unwrap(); }
/// ```

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct SeqPos {
    refid: String,
    pos: isize,
}

impl SeqPos {
    /// Construct a new sequence position
    ///
    /// ```
    /// use bio::seq_loc::seq_pos::SeqPos;
    /// let start = SeqPos::new("chrIV".to_owned(), 683946);
    /// ```
    pub fn new(refid: String, pos: isize) -> SeqPos {
        SeqPos{ refid: refid, pos: pos }
    }

    /// Name of the reference sequence (chromosome name, etc.)
    pub fn refid(&self) -> &str { &self.refid.as_str() }
    /// Position on the reference sequence (0-based).
    pub fn pos(&self) -> isize { self.pos }

    /// Construct a stranded sequence position on the specified `strand`.
    pub fn on_strand(self, strand: ReqStrand) -> SeqStrandPos {
        SeqStrandPos::new(self.refid, self.pos, strand)
    }

    /// Generate a BED format `Record` for the position.
    ///
    /// This record will have length 1 and no strand, and when created it will have an empty name.
    pub fn to_bed(&self) -> bed::Record {
        let mut bed = bed::Record::new();
        bed.set_chrom(&self.refid);
        bed.set_start(self.pos as u64);
        bed.set_end((self.pos + 1) as u64);
        bed.set_name("");
        bed.set_score("0");
        bed.push_aux(".");
        bed
    }
}

impl fmt::Display for SeqPos {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}", self.refid, self.pos)
    }
}

impl FromStr for SeqPos {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (refid, posstr) = break_refid(s)?;
        let pos = posstr.parse::<isize>()
            .map_err(|e| ParseError::BadInteger(s.to_owned(), e))?;
        Ok( SeqPos::new(refid.to_owned(), pos) )
    }
}

impl From<SeqStrandPos> for SeqPos {
    fn from(ssp: SeqStrandPos) -> Self {
        SeqPos::new(ssp.refid, ssp.pos)
    }
}

/// Position on a named sequence, with a specific strand
///
/// The display format for a `SeqStrandPos` is _chr:pos(+/-)_.
///
/// ```
/// # use bio::seq_loc::ParseError;
/// # fn try_main() -> Result<(), Box<ParseError>> {
/// use bio::seq_loc::seq_pos::SeqStrandPos;
/// use bio::seq_loc::ReqStrand;
/// let start = SeqStrandPos::new("chrIV".to_owned(), 683946, ReqStrand::Reverse);
/// let start_str = start.to_string();
/// assert_eq!(start_str, "chrIV:683946(-)");
/// let start_str_pos: SeqStrandPos = start_str.parse()?;
/// assert_eq!(start, start_str_pos);
/// # Ok(()) 
/// # }
/// # fn main() { try_main().unwrap(); }
/// ```
#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct SeqStrandPos {
    refid: String,
    pos: isize,
    strand: ReqStrand,
}

impl SeqStrandPos {
    /// Construct a new sequence position
    ///
    /// ```
    /// use bio::seq_loc::seq_pos::SeqStrandPos;
    /// use bio::seq_loc::ReqStrand;
    /// let start = SeqStrandPos::new("chrIV".to_owned(), 683946, ReqStrand::Reverse);
    /// ```
    pub fn new(refid: String, pos: isize, strand: ReqStrand) -> Self {
        SeqStrandPos{ refid: refid, pos: pos, strand: strand }
    }

    /// Name of the reference sequence (chromosome name, etc.)
    pub fn refid(&self) -> &str { &self.refid.as_str() }
    /// Position on the reference sequence (0-based).
    pub fn pos(&self) -> isize { self.pos }
    /// `Strand` of the position
    pub fn strand(&self) -> ReqStrand { self.strand }

    /// Generate a BED format `Record` for the position.
    ///
    /// This record will have length 1, and when created it will have an empty name.
    pub fn to_bed(&self) -> bed::Record {
        let mut bed = bed::Record::new();
        bed.set_chrom(&self.refid);
        bed.set_start(self.pos as u64);
        bed.set_end((self.pos + 1) as u64);
        bed.set_name("");
        bed.set_score("0");
        bed.push_aux(match self.strand {
            ReqStrand::Forward => "+",
            ReqStrand::Reverse => "-",
        });
        bed
    }
}

impl fmt::Display for SeqStrandPos {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}{}", self.refid, self.pos, 
               match self.strand {
                   ReqStrand::Forward => "(+)",
                   ReqStrand::Reverse => "(-)",
               })
    }
}

impl FromStr for SeqStrandPos {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (refid, rest) = break_refid(s)?;
        let (posstr, strand) = break_fwd_rev(rest)?;
        let pos = posstr.parse::<isize>()
            .map_err(|e| ParseError::BadInteger(s.to_owned(), e))?;
        Ok( SeqStrandPos::new(refid.to_owned(), pos, strand) )
    }
}

