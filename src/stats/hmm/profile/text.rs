// Copyright 2026 Sahil Rajput
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::stats::hmm::profile::{ProfileError, ProfileHmm, NEG_INF};
use std::fmt::Write as _;

const MAGIC: &str = "PHMM";
const VERSION: &str = "1";

pub(crate) fn to_text(profile: &ProfileHmm) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "{} {}", MAGIC, VERSION);
    let _ = writeln!(out, "COLUMNS {}", profile.num_cols);
    let alphabet: String = profile.symbols.iter().map(|&b| b as char).collect();
    let _ = writeln!(out, "ALPHABET {}", alphabet);
    let _ = write!(out, "NULL");
    for &v in &profile.null {
        let _ = write!(out, " {}", v);
    }
    let _ = writeln!(out);
    for k in 0..=profile.num_cols {
        let _ = write!(out, "NODE {}", k);
        if k >= 1 {
            let _ = write!(out, " E");
            for &v in &profile.emission[k] {
                let _ = write!(out, " {}", v);
            }
        }
        let _ = write!(out, " T");
        for &v in &[
            profile.t_mm[k],
            profile.t_mi[k],
            profile.t_md[k],
            profile.t_im[k],
            profile.t_ii[k],
            profile.t_id[k],
            profile.t_dm[k],
            profile.t_di[k],
            profile.t_dd[k],
        ] {
            let _ = write!(out, " {}", v);
        }
        let _ = writeln!(out);
    }
    out
}

fn parse_floats(tokens: &[&str]) -> Result<Vec<f64>, ProfileError> {
    tokens
        .iter()
        .map(|t| {
            t.parse::<f64>()
                .map_err(|_| ProfileError::Parse(format!("expected number, got '{}'", t)))
        })
        .collect()
}

pub(crate) fn from_text(text: &str) -> Result<ProfileHmm, ProfileError> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());

    let header: Vec<&str> = lines
        .next()
        .ok_or_else(|| ProfileError::Parse("empty input".to_string()))?
        .split_whitespace()
        .collect();
    if header.len() != 2 || header[0] != MAGIC || header[1] != VERSION {
        return Err(ProfileError::Parse("bad header".to_string()));
    }

    let cols_line: Vec<&str> = lines
        .next()
        .ok_or_else(|| ProfileError::Parse("missing COLUMNS".to_string()))?
        .split_whitespace()
        .collect();
    if cols_line.len() != 2 || cols_line[0] != "COLUMNS" {
        return Err(ProfileError::Parse("bad COLUMNS line".to_string()));
    }
    let l: usize = cols_line[1]
        .parse()
        .map_err(|_| ProfileError::Parse("bad column count".to_string()))?;

    let alpha_line: Vec<&str> = lines
        .next()
        .ok_or_else(|| ProfileError::Parse("missing ALPHABET".to_string()))?
        .split_whitespace()
        .collect();
    if alpha_line.len() != 2 || alpha_line[0] != "ALPHABET" {
        return Err(ProfileError::Parse("bad ALPHABET line".to_string()));
    }
    let symbols: Vec<u8> = alpha_line[1].bytes().collect();
    let alphabet_size = symbols.len();
    if alphabet_size == 0 {
        return Err(ProfileError::Parse("empty alphabet".to_string()));
    }
    let mut rank_of = vec![-1i32; 256];
    for (rank, &sym) in symbols.iter().enumerate() {
        rank_of[sym as usize] = rank as i32;
    }

    let null_line: Vec<&str> = lines
        .next()
        .ok_or_else(|| ProfileError::Parse("missing NULL".to_string()))?
        .split_whitespace()
        .collect();
    if null_line.is_empty() || null_line[0] != "NULL" {
        return Err(ProfileError::Parse("bad NULL line".to_string()));
    }
    let null = parse_floats(&null_line[1..])?;
    if null.len() != alphabet_size {
        return Err(ProfileError::Parse("NULL length mismatch".to_string()));
    }

    let mut emission = vec![vec![NEG_INF; alphabet_size]; l + 1];
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
        let tokens: Vec<&str> = lines
            .next()
            .ok_or_else(|| ProfileError::Parse(format!("missing NODE {}", k)))?
            .split_whitespace()
            .collect();
        if tokens.len() < 2 || tokens[0] != "NODE" {
            return Err(ProfileError::Parse(format!("bad NODE {} line", k)));
        }
        let node_index: usize = tokens[1]
            .parse()
            .map_err(|_| ProfileError::Parse("bad node index".to_string()))?;
        if node_index != k {
            return Err(ProfileError::Parse(format!(
                "expected node {}, found {}",
                k, node_index
            )));
        }
        let mut pos = 2;
        if k >= 1 {
            if tokens.get(pos) != Some(&"E") {
                return Err(ProfileError::Parse(format!("node {} missing E", k)));
            }
            pos += 1;
            if pos + alphabet_size > tokens.len() {
                return Err(ProfileError::Parse(format!("node {} short emission", k)));
            }
            let emit = parse_floats(&tokens[pos..pos + alphabet_size])?;
            emission[k] = emit;
            pos += alphabet_size;
        }
        if tokens.get(pos) != Some(&"T") {
            return Err(ProfileError::Parse(format!("node {} missing T", k)));
        }
        pos += 1;
        let trans = parse_floats(&tokens[pos..])?;
        if trans.len() != 9 {
            return Err(ProfileError::Parse(format!(
                "node {} expects 9 transitions",
                k
            )));
        }
        t_mm[k] = trans[0];
        t_mi[k] = trans[1];
        t_md[k] = trans[2];
        t_im[k] = trans[3];
        t_ii[k] = trans[4];
        t_id[k] = trans[5];
        t_dm[k] = trans[6];
        t_di[k] = trans[7];
        t_dd[k] = trans[8];
    }

    Ok(ProfileHmm {
        num_cols: l,
        symbols,
        rank_of,
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

#[cfg(test)]
mod tests {
    use crate::stats::hmm::profile::example_profile;
    use crate::stats::hmm::profile::{ProfileError, ProfileHmm};

    #[test]
    fn test_round_trip_preserves_scores() {
        let profile = example_profile();
        let restored = ProfileHmm::from_text(&profile.to_text()).unwrap();
        assert_eq!(restored.num_columns(), profile.num_columns());
        assert_eq!(restored.to_text(), profile.to_text());
        for query in [
            &b"ACGTAC"[..],
            &b"ACAC"[..],
            &b"TTTT"[..],
            &b"ACGTTTAC"[..],
            &b""[..],
        ] {
            assert_eq!(
                profile.bit_score(query).unwrap(),
                restored.bit_score(query).unwrap()
            );
            assert_eq!(
                profile.forward(query).unwrap(),
                restored.forward(query).unwrap()
            );
        }
    }

    #[test]
    fn test_malformed_text_is_rejected() {
        let text = example_profile().to_text();
        let parse = |input: &str| ProfileHmm::from_text(input).unwrap_err();
        assert!(matches!(parse(""), ProfileError::Parse(_)));
        assert!(matches!(
            parse(&text.replacen("PHMM 1", "PHMM 2", 1)),
            ProfileError::Parse(_)
        ));
        assert!(matches!(
            parse(&text.replacen("COLUMNS 6", "COLUMNS 7", 1)),
            ProfileError::Parse(_)
        ));
        let truncated: Vec<&str> = text.lines().take(4).collect();
        assert!(matches!(
            parse(&truncated.join("\n")),
            ProfileError::Parse(_)
        ));
    }
}
