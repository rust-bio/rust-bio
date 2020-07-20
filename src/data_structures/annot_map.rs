//! Efficient container for locations annotated across a set of named
//! reference sequences.
//!
//! # Example
//!
//! ```
//! extern crate bio_types;
//! use bio::data_structures::annot_map::AnnotMap;
//! use bio_types::annot::contig::Contig;
//! use bio_types::strand::ReqStrand;
//!
//! // Insert a String annotation into the annotation map at a specified location.
//! let mut genes: AnnotMap<String, String> = AnnotMap::new();
//! let tma22 = Contig::new(
//!     "chrX".to_owned(),
//!     461829,
//!     462426 - 461829,
//!     ReqStrand::Forward,
//! );
//! genes.insert_at("TMA22".to_owned(), &tma22);
//!
//! // Find annotations that overlap a specific query
//! let query = Contig::new("chrX".to_owned(), 462400, 100, ReqStrand::Forward);
//! let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
//! assert_eq!(hits, vec!["TMA22"]);
//! ```

use std::collections::HashMap;
use std::hash::Hash;

use crate::data_structures::interval_tree;
use crate::data_structures::interval_tree::{IntervalTree, IntervalTreeIterator};
use crate::utils::Interval;
use bio_types::annot::loc::Loc;

/// Efficient container for querying annotations, using `HashMap` and
/// `IntervalTree`.
///
/// The container is parameterized over the type of the reference
/// sequence names `R` (which is often a `String`) and the type of the
/// contained objects `T`.
///
/// The container finds annotations that overlap a specific query
/// location. Overlaps are identified without regard for strandedness
/// and without regard for e.g. spliced-out introns within the
/// annotation or the query.
///
/// Thus, the overlapping annotations identified by querying a
/// `AnnotMap` may need further filtering.
#[derive(Debug, Clone)]
pub struct AnnotMap<R, T>
where
    R: Hash + Eq,
{
    refid_itrees: HashMap<R, IntervalTree<isize, T>>,
}

impl<R, T> Default for AnnotMap<R, T>
where
    R: Eq + Hash,
{
    fn default() -> Self {
        AnnotMap {
            refid_itrees: HashMap::new(),
        }
    }
}

impl<R, T> AnnotMap<R, T>
where
    R: Eq + Hash,
{
    /// Create a new, empty `AnnotMap`. Used in conjunction with `insert_at`
    /// or `insert_loc`.
    pub fn new() -> Self {
        Default::default()
    }

    /// Insert an object into the container at a specified location (`Loc`).
    ///
    /// # Arguments
    ///
    /// * `data` - any type of data to be inserted at the location / region
    /// * `location` - any object with the `Loc` trait implemented, determining
    ///   the Range at which to insert the `data`
    ///
    /// # Example
    ///
    /// ```
    /// extern crate bio_types;
    /// use bio::data_structures::annot_map::AnnotMap;
    /// use bio_types::annot::contig::Contig;
    /// use bio_types::strand::ReqStrand;
    ///
    /// let mut genes: AnnotMap<String, String> = AnnotMap::new();
    /// let tma22 = Contig::new(
    ///     "chrX".to_owned(),
    ///     461829,
    ///     462426 - 461829,
    ///     ReqStrand::Forward,
    /// );
    /// genes.insert_at("TMA22".to_owned(), &tma22);
    /// ```
    pub fn insert_at<L>(&mut self, data: T, location: &L)
    where
        R: Eq + Hash + Clone,
        L: Loc<RefID = R>,
    {
        let itree = self
            .refid_itrees
            .entry(location.refid().clone())
            .or_insert_with(IntervalTree::new);
        let rng = location.start()..(location.start() + (location.length() as isize));
        itree.insert(rng, data);
    }

    /// Create an `Iterator` that will visit all entries that overlap
    /// a query location.
    pub fn find<'a, L>(&'a self, location: &'a L) -> AnnotMapIterator<'a, R, T>
    where
        L: Loc<RefID = R>,
    {
        if let Some(itree) = self.refid_itrees.get(location.refid()) {
            let interval = location.start()..(location.start() + (location.length() as isize));
            let itree_iter = itree.find(interval);
            AnnotMapIterator {
                itree_iter: Some(itree_iter),
                refid: location.refid(),
            }
        } else {
            AnnotMapIterator {
                itree_iter: None,
                refid: location.refid(),
            }
        }
    }
}

impl<R, T> AnnotMap<R, T>
where
    R: Eq + Hash + Clone,
    T: Loc<RefID = R>,
{
    /// Insert an object with the `Loc` trait into the container at
    /// its location.
    ///
    /// This inserts all of `data` at the Range of length `data.length()`
    /// that starts at `data.start()`.
    ///
    /// # Example
    ///
    /// ```
    /// extern crate bio_types;
    /// use bio::data_structures::annot_map::AnnotMap;
    /// use bio_types::annot::contig::Contig;
    /// use bio_types::strand::ReqStrand;
    ///
    /// let mut gene_locs = AnnotMap::new();
    /// let tma19 = Contig::new(
    ///     String::from("chrXI"),
    ///     334412,
    ///     (334916 - 334412),
    ///     ReqStrand::Reverse,
    /// );
    /// let assert_copy = tma19.clone();
    /// gene_locs.insert_loc(tma19);
    /// // Find annotations that overlap a specific query
    /// let query = Contig::new(String::from("chrXI"), 334400, 100, ReqStrand::Reverse);
    /// let hits: Vec<&Contig<String, ReqStrand>> = gene_locs.find(&query).map(|e| e.data()).collect();
    /// assert_eq!(hits, vec![&assert_copy]);
    /// ```
    pub fn insert_loc(&mut self, data: T) {
        let itree = self
            .refid_itrees
            .entry(data.refid().clone())
            .or_insert_with(IntervalTree::new);
        let rng = data.start()..(data.start() + (data.length() as isize));
        itree.insert(rng, data);
    }
}

/// A view of one annotation in a `AnnotMap` container.
#[derive(Debug, Clone)]
pub struct Entry<'a, R, T>
where
    R: Eq + Hash,
{
    itree_entry: interval_tree::Entry<'a, isize, T>,
    refid: &'a R,
}

impl<'a, R, T> Entry<'a, R, T>
where
    R: Eq + Hash,
{
    /// Return a reference to the data value in the `AnnotMap`.
    pub fn data(&self) -> &'a T {
        self.itree_entry.data()
    }

    /// Return a reference to the interval spanned by the annotation.
    pub fn interval(&self) -> &'a Interval<isize> {
        self.itree_entry.interval()
    }

    /// Return a reference to the identifier of the annotated reference sequence.
    pub fn refid(&self) -> &'a R {
        self.refid
    }
}

/// An iterator over annotation entries (of type `Entry`) in a
/// `AnnotMap`.
///
/// This struct is created by the `find` function on `AnnotMap`.
pub struct AnnotMapIterator<'a, R, T>
where
    R: Eq + Hash,
{
    itree_iter: Option<IntervalTreeIterator<'a, isize, T>>,
    refid: &'a R,
}

impl<'a, R, T> Iterator for AnnotMapIterator<'a, R, T>
where
    R: 'a + Eq + Hash,
    T: 'a,
{
    type Item = Entry<'a, R, T>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.itree_iter {
            Some(ref mut iter) => match iter.next() {
                Some(next_itree) => Some(Entry {
                    itree_entry: next_itree,
                    refid: self.refid,
                }),
                None => None,
            },
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bio_types::annot::contig::Contig;
    use bio_types::strand::ReqStrand;

    #[test]
    fn lookup() {
        let mut genes: AnnotMap<String, String> = AnnotMap::new();
        genes.insert_at(
            "TMA22".to_owned(),
            &Contig::new(
                "chrX".to_owned(),
                461829,
                462426 - 461829,
                ReqStrand::Forward,
            ),
        );
        genes.insert_at(
            "TMA19".to_owned(),
            &Contig::new(
                "chrXI".to_owned(),
                334412,
                334916 - 334412,
                ReqStrand::Reverse,
            ),
        );

        let query = Contig::new("chrX".to_owned(), 462400, 100, ReqStrand::Forward);
        let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        assert_eq!(hits, vec!["TMA22"]);

        let query = Contig::new("chrXI".to_owned(), 334400, 100, ReqStrand::Forward);
        let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        assert_eq!(hits, vec!["TMA19"]);

        let query = Contig::new("chrXI".to_owned(), 334916, 100, ReqStrand::Forward);
        let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        assert!(hits.is_empty());

        let query = Contig::new("chrX".to_owned(), 461729, 100, ReqStrand::Forward);
        let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        assert!(hits.is_empty());

        let query = Contig::new("chrXI".to_owned(), 462400, 100, ReqStrand::Forward);
        let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        assert!(hits.is_empty());

        let query = Contig::new("NotFound".to_owned(), 0, 0, ReqStrand::Forward);
        let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        assert!(hits.is_empty());
    }

    #[test]
    fn overlaps() {
        let mut genes: AnnotMap<String, String> = AnnotMap::new();

        genes.insert_at(
            "a".to_owned(),
            &Contig::new("chr01".to_owned(), 1000, 1000, ReqStrand::Forward),
        );
        genes.insert_at(
            "b".to_owned(),
            &Contig::new("chr01".to_owned(), 1300, 1000, ReqStrand::Forward),
        );
        genes.insert_at(
            "c".to_owned(),
            &Contig::new("chr01".to_owned(), 1700, 1000, ReqStrand::Forward),
        );
        genes.insert_at(
            "d".to_owned(),
            &Contig::new("chr01".to_owned(), 2200, 1000, ReqStrand::Forward),
        );

        let query = Contig::new("chr01".to_owned(), 1050, 100, ReqStrand::Forward);
        let mut hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        hits.sort();
        assert_eq!(hits, vec!["a"]);

        let query = Contig::new("chr01".to_owned(), 1450, 100, ReqStrand::Forward);
        let mut hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        hits.sort();
        assert_eq!(hits, vec!["a", "b"]);

        let query = Contig::new("chr01".to_owned(), 1850, 100, ReqStrand::Forward);
        let mut hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        hits.sort();
        assert_eq!(hits, vec!["a", "b", "c"]);

        let query = Contig::new("chr01".to_owned(), 2250, 100, ReqStrand::Forward);
        let mut hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        hits.sort();
        assert_eq!(hits, vec!["b", "c", "d"]);

        let query = Contig::new("chr01".to_owned(), 2650, 100, ReqStrand::Forward);
        let mut hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
        hits.sort();
        assert_eq!(hits, vec!["c", "d"]);
    }
}
