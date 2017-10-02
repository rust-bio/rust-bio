use std::num::ParseIntError;

pub mod seq_pos;

#[derive(Debug,Clone,Hash,PartialEq,Eq,Ord,PartialOrd,Copy)]
pub enum ReqStrand {
    Forward,
    Reverse
}

#[derive(Debug,Clone)]
pub enum ParseError {
    NoColon(String),
    NoRefid(String),
    NoPosition(String),
    BadInteger(String, ParseIntError),
    BadStrand(String),
}

// Break a position display string into a reference name part and the
// "rest"
fn break_refid(s: &str) -> Result<(&str,&str),ParseError> {
    let breakpt = s.find(':').ok_or_else(|| ParseError::NoColon(s.to_owned()))?;
    let refid = s.get(..breakpt).ok_or_else(|| ParseError::NoRefid(s.to_owned()))?;
    let rest = s.get((breakpt+1)..).ok_or_else(|| ParseError::NoPosition(s.to_owned()))?;
    Ok( (refid, rest) )
}

// Break a trailing strand flag off of a string
fn break_fwd_rev(s: &str) -> Result<(&str,ReqStrand),ParseError> {
    match s.get((s.len()-3)..) {
        Some("(+)") => Ok( (s.get(0..(s.len()-3)).unwrap_or(""), ReqStrand::Forward) ),
        Some("(-)") => Ok( (s.get(0..(s.len()-3)).unwrap_or(""), ReqStrand::Reverse) ),
        _ => Err(ParseError::BadStrand(s.to_owned())),
    }
}
