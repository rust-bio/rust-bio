//! Efficient container for locations annotated across a set of named
//! reference sequences.
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

impl<R, T> AnnotMap<R, T>
where
    R: Eq + Hash,
{
    /// Creates a new, empty `AnnotMap`.
    ///
    /// ```
    /// extern crate bio;
    /// use bio::data_structures::annot_map::AnnotMap;
    /// let mut genes: AnnotMap<String,String> = AnnotMap::new();
    /// ```
    pub fn new() -> Self {
        AnnotMap {
            refid_itrees: HashMap::new(),
        }
    }

    /// Inserts an object into the container at a specified location.
    ///
    /// # Arguments
    ///
    /// `data` is the data item to be inserted into the annotation map
    ///
    /// `location` is the location associated with the data item.
    ///
    /// ```
    /// extern crate bio_types;
    /// extern crate bio;
    /// use bio_types::strand::ReqStrand;
    /// use bio_types::annot::contig::Contig;
    /// use bio::data_structures::annot_map::AnnotMap;
    /// let mut genes: AnnotMap<String,String> = AnnotMap::new();
    /// let tma22 = Contig::new("chrX".to_owned(), 461829, 462426 - 461829, ReqStrand::Forward);
    /// genes.insert_at("TMA22".to_owned(), &tma22);
    /// let tma19 = Contig::new("chrXI".to_owned(), 334412, (334916 - 334412), ReqStrand::Reverse);
    /// genes.insert_at("TMA19".to_owned(), &tma19);
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

    /// Creates an `Iterator` that will visit all entries that overlap
    /// a query location.
    ///
    ///
    /// # Arguments
    ///
    /// `location` is the annotation location to be searched for
    /// overlapping data items.
    ///
    /// ```
    /// extern crate bio_types;
    /// extern crate bio;
    /// use bio_types::strand::ReqStrand;
    /// use bio_types::annot::contig::Contig;
    /// use bio::data_structures::annot_map::AnnotMap;
    /// let mut genes: AnnotMap<String,String> = AnnotMap::new();
    /// let tma22 = Contig::new("chrX".to_owned(), 461829, 462426 - 461829, ReqStrand::Forward);
    /// genes.insert_at("TMA22".to_owned(), &tma22);
    /// let tma19 = Contig::new("chrXI".to_owned(), 334412, (334916 - 334412), ReqStrand::Reverse);
    /// genes.insert_at("TMA19".to_owned(), &tma19);
    /// let query = Contig::new("chrX".to_owned(), 462400, 100, ReqStrand::Forward);
    /// let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
    /// assert_eq!(hits, vec!["TMA22"]);
    /// let query = Contig::new("chrXI".to_owned(), 334400, 100, ReqStrand::Forward);
    /// let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
    /// assert_eq!(hits, vec!["TMA19"]);
    /// ```
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
    /// Inserts an object with the `Loc` trait into the container at
    /// its location.
    ///
    /// # Argument
    ///
    /// `data` is the data to be inserted based on its location.
    ///
    /// Equivalent to inserting `data` at `data.contig()`.
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
    /// Returns a reference to the data value in the `AnnotMap`.
    pub fn data(&self) -> &'a T {
        self.itree_entry.data()
    }

    /// Returns a reference to the interval spanned by the annotation.
    pub fn interval(&self) -> &'a Interval<isize> {
        self.itree_entry.interval()
    }

    /// Returns a reference to the identifier of the annotated reference sequence.
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
