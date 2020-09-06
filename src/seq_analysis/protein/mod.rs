//! The data is based on Gasteiger _et al_ (2005) i.e. in accordance with
//! [ExPASy's ProtParam Tool](https://web.expasy.org/protparam/)
//! # References
//! - Gasteiger E., Hoogland C., Gattiker A., Duvaud S., Wilkins M.R., Appel R.D., Bairoch A.
//!   Protein Identification and Analysis Tools on the ExPASy Server, in _The Proteomics Protocols
//!   Handbook_, Humana Press (2005). pp. 571-607

use crate::utils::TextSlice;
use std::collections::BTreeMap;

pub mod isoelectric_point;
pub use isoelectric_point::isoelectric_point;

pub fn count_aa(seq: TextSlice) -> BTreeMap<u8, u32> {
    let mut res: BTreeMap<u8, u32> = BTreeMap::new();
    for &aa in seq {
        let count = res.entry(aa).or_insert(0);
        *count += 1;
    }
    res
}
