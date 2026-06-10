//! Interval tree, a data structure for efficiently storing and searching intervals.
//!
//! This data structure uses an `Interval` type to define the intervals. It is implemented by
//! wrapping the `std::ops::Range` in a newtype. An interval must have a positive with and the interval bounds
//! may be specified by any type satisfies both the `std::cmp::Ord` and `Clone` trait. Because
//! `Interval` implements `From<Range>` you can also uses normal `Range` arguments in the
//! `insert` and `find` functions. This implicit conversion will panic if a negative-width range is
//! supplied.
//!
//! Upon inserting an interval may be associated with a data value. The intervals are stored in
//! an augmented AVL-tree which allows for efficient inserting and querying.
//!
//! # Example
//! ```
//! use bio::data_structures::interval_tree::IntervalTree;
//! use bio::utils::Interval;
//!
//! let mut tree = IntervalTree::new();
//! tree.insert(11..20, "Range_1");
//! tree.insert(25..30, "Range_2");
//! for r in tree.find(15..25) {
//!     assert_eq!(r.interval().start, 11);
//!     assert_eq!(r.interval().end, 20);
//!     assert_eq!(r.interval(), &(Interval::from(11..20)));
//!     assert_eq!(r.data(), &"Range_1");
//! }
//! ```

use crate::utils::Interval;
use std::cmp;
use std::iter::FromIterator;
use std::mem;

/// An interval tree for storing intervals with data
#[derive(Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
pub struct IntervalTree<N: Ord + Clone, D> {
    root: Option<Node<N, D>>,
}

impl<N: Ord + Clone, D> Default for IntervalTree<N, D> {
    fn default() -> Self {
        Self { root: None }
    }
}

/// A `find` query on the interval tree does not directly return references to the nodes in the tree, but
/// wraps the fields `interval` and `data` in an `Entry`.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Serialize)]
pub struct Entry<'a, N: Ord + Clone, D> {
    data: &'a D,
    interval: &'a Interval<N>,
}

impl<'a, N: Ord + Clone + 'a, D: 'a> Entry<'a, N, D> {
    /// Get a reference to the data for this entry
    pub fn data(&self) -> &'a D {
        self.data
    }

    /// Get a reference to the interval for this entry
    pub fn interval(&self) -> &'a Interval<N> {
        self.interval
    }
}

/// An `IntervalTreeIterator` is returned by `Intervaltree::find` and iterates over the entries
/// overlapping the query
#[derive(Default, Clone, Eq, PartialEq, Hash, Debug, Serialize)]
pub struct IntervalTreeIterator<'a, N: Ord + Clone, D> {
    nodes: Vec<&'a Node<N, D>>,
    interval: Interval<N>,
}

impl<'a, N: Ord + Clone + 'a, D: 'a> Iterator for IntervalTreeIterator<'a, N, D> {
    type Item = Entry<'a, N, D>;

    fn next(&mut self) -> Option<Entry<'a, N, D>> {
        loop {
            let candidate = self.nodes.pop()?;

            // stop traversal if the query interval is beyond the current node and all children
            if self.interval.start < candidate.max {
                if let Some(ref left) = candidate.left {
                    self.nodes.push(left);
                }

                // don't traverse right if the query interval is completely before the current
                // interval
                if self.interval.end > candidate.interval.start {
                    if let Some(ref right) = candidate.right {
                        self.nodes.push(right);
                    }

                    // overlap is only possible if both tests pass
                    if intersect(&self.interval, &candidate.interval) {
                        return Some(Entry {
                            data: &candidate.value,
                            interval: &candidate.interval,
                        });
                    }
                }
            }
        }
    }
}

/// A `find_mut` query on the interval tree does not directly return references to the nodes in the tree, but
/// wraps the fields `interval` and `data` in an `EntryMut`. Only the data part can be mutably accessed
/// using the `data` method
#[derive(Eq, PartialEq, Hash, Debug, Serialize)]
pub struct EntryMut<'a, N: Ord + Clone, D> {
    data: &'a mut D,
    interval: &'a Interval<N>,
}

impl<'a, N: Ord + Clone + 'a, D: 'a> EntryMut<'a, N, D> {
    /// Get a mutable reference to the data for this entry
    pub fn data(&'a mut self) -> &'a mut D {
        self.data
    }

    /// Get a reference to the interval for this entry
    pub fn interval(&self) -> &'a Interval<N> {
        self.interval
    }
}

/// An `IntervalTreeIteratorMut` is returned by `Intervaltree::find_mut` and iterates over the entries
/// overlapping the query allowing mutable access to the data `D`, not the `Interval`.
#[derive(Default, Eq, PartialEq, Hash, Debug, Serialize)]
pub struct IntervalTreeIteratorMut<'a, N: Ord + Clone, D> {
    nodes: Vec<&'a mut Node<N, D>>,
    interval: Interval<N>,
}

impl<'a, N: Ord + Clone + 'a, D: 'a> Iterator for IntervalTreeIteratorMut<'a, N, D> {
    type Item = EntryMut<'a, N, D>;

    fn next(&mut self) -> Option<EntryMut<'a, N, D>> {
        loop {
            let candidate = self.nodes.pop()?;

            // stop traversal if the query interval is beyond the current node and all children
            if self.interval.start < candidate.max {
                if let Some(ref mut left) = candidate.left {
                    self.nodes.push(left);
                }

                // don't traverse right if the query interval is completely before the current interval
                if self.interval.end > candidate.interval.start {
                    if let Some(ref mut right) = candidate.right {
                        self.nodes.push(right);
                    }

                    // overlap is only possible if both tests pass
                    if intersect(&self.interval, &candidate.interval) {
                        return Some(EntryMut {
                            data: &mut candidate.value,
                            interval: &candidate.interval,
                        });
                    }
                }
            }
        }
    }
}

impl<N: Clone + Ord, D> IntervalTree<N, D> {
    /// Creates a new empty `IntervalTree`
    pub fn new() -> Self {
        Default::default()
    }

    /// Inserts an `Interval` into the tree and associates it with `data`
    pub fn insert<I: Into<Interval<N>>>(&mut self, interval: I, data: D) {
        let interval = interval.into();
        match self.root {
            Some(ref mut n) => n.insert(interval, data),
            None => self.root = Some(Node::new(interval, data)),
        };
    }

    /// Uses the provided `Interval` to find overlapping intervals in the tree and returns an
    /// `IntervalTreeIterator`
    pub fn find<I: Into<Interval<N>>>(&self, interval: I) -> IntervalTreeIterator<'_, N, D> {
        let interval = interval.into();
        match self.root {
            Some(ref n) => IntervalTreeIterator {
                nodes: vec![n],
                interval,
            },
            None => {
                let nodes = vec![];
                IntervalTreeIterator { nodes, interval }
            }
        }
    }

    /// Uses the provided `Interval` to find overlapping intervals in the tree and returns an
    /// `IntervalTreeIteratorMut` that allows mutable access to the `data`
    pub fn find_mut<I: Into<Interval<N>>>(
        &mut self,
        interval: I,
    ) -> IntervalTreeIteratorMut<'_, N, D> {
        let interval = interval.into();
        match self.root {
            Some(ref mut n) => IntervalTreeIteratorMut {
                nodes: vec![n],
                interval,
            },
            None => {
                let nodes = vec![];
                IntervalTreeIteratorMut { nodes, interval }
            }
        }
    }
}

impl<N: Clone + Ord, D, R: Into<Interval<N>>> FromIterator<(R, D)> for IntervalTree<N, D> {
    fn from_iter<I: IntoIterator<Item = (R, D)>>(iter: I) -> Self {
        let mut tree = IntervalTree::new();
        for r in iter {
            tree.insert(r.0, r.1);
        }
        tree
    }
}

#[derive(Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
struct Node<N: Ord + Clone, D> {
    // actual interval data
    interval: Interval<N>,
    value: D,
    // tree metadata
    max: N,
    height: i64,
    left: Option<Box<Node<N, D>>>,
    right: Option<Box<Node<N, D>>>,
}

impl<N: Ord + Clone, D> Node<N, D> {
    fn new(interval: Interval<N>, data: D) -> Self {
        let max = interval.end.clone();
        Node {
            interval,
            max,
            height: 1,
            value: data,
            left: None,
            right: None,
        }
    }

    fn insert(&mut self, interval: Interval<N>, data: D) {
        if interval.start <= self.interval.start {
            if let Some(ref mut son) = self.left {
                son.insert(interval, data);
            } else {
                self.left = Some(Box::new(Node::new(interval, data)));
            }
        } else if let Some(ref mut son) = self.right {
            son.insert(interval, data);
        } else {
            self.right = Some(Box::new(Node::new(interval, data)));
        }
        self.repair();
    }

    fn update_height(&mut self) {
        let left_h = self.left.as_ref().map_or(0, |n| n.height);
        let right_h = self.right.as_ref().map_or(0, |n| n.height);
        self.height = 1 + cmp::max(left_h, right_h);
    }

    fn update_max(&mut self) {
        self.max = self.interval.end.clone();
        if let Some(ref n) = self.left {
            if self.max < n.max {
                self.max = n.max.clone();
            }
        }
        if let Some(ref n) = self.right {
            if self.max < n.max {
                self.max = n.max.clone();
            }
        }
    }

    fn repair(&mut self) {
        let left_h = self.left.as_ref().map_or(0, |n| n.height);
        let right_h = self.right.as_ref().map_or(0, |n| n.height);
        // each case - update both height and max
        if (left_h - right_h).abs() <= 1 {
            self.update_height();
            self.update_max();
        } else if right_h > left_h {
            {
                let right = self
                    .right
                    .as_mut()
                    .expect("Invalid tree: leaf is taller than its sibling.");
                let right_left_h = right.left.as_ref().map_or(0, |n| n.height);
                let right_right_h = right.right.as_ref().map_or(0, |n| n.height);
                if right_left_h > right_right_h {
                    right.rotate_right();
                }
            }
            self.rotate_left();
        } else {
            {
                let left = self
                    .left
                    .as_mut()
                    .expect("Invalid tree: leaf is taller than its sibling.");
                let left_right_h = left.right.as_ref().map_or(0, |n| n.height);
                let left_left_h = left.left.as_ref().map_or(0, |n| n.height);
                if left_right_h > left_left_h {
                    left.rotate_left();
                }
            }
            self.rotate_right();
        }
    }

    fn rotate_left(&mut self) {
        let mut new_root = self.right.take().unwrap();
        let t1 = self.left.take();
        let t2 = new_root.left.take();
        let t3 = new_root.right.take();
        swap_interval_data(self, &mut *new_root);

        new_root.left = t1;
        new_root.right = t2;
        new_root.update_height();
        new_root.update_max();

        self.right = t3;
        self.left = Some(new_root);
        self.update_height();
        self.update_max();
    }

    fn rotate_right(&mut self) {
        let mut new_root = self.left.take().unwrap();
        let t1 = new_root.left.take();
        let t2 = new_root.right.take();
        let t3 = self.right.take();
        swap_interval_data(self, &mut *new_root);

        new_root.left = t2;
        new_root.right = t3;
        new_root.update_height();
        new_root.update_max();

        self.left = t1;
        self.right = Some(new_root);
        self.update_height();
        self.update_max();
    }
}

fn swap_interval_data<N: Ord + Clone, D>(node_1: &mut Node<N, D>, node_2: &mut Node<N, D>) {
    mem::swap(&mut node_1.value, &mut node_2.value);
    mem::swap(&mut node_1.interval, &mut node_2.interval);
}

fn intersect<N: Ord + Clone>(range_1: &Interval<N>, range_2: &Interval<N>) -> bool {
    range_1.start < range_1.end
        && range_2.start < range_2.end
        && range_1.end > range_2.start
        && range_1.start < range_2.end
}

#[cfg(test)]
mod tests {
    use super::{Entry, IntervalTree, Node};
    use crate::utils::Interval;
    use std::cmp;
    use std::cmp::{max, min};
    use std::ops::Range;

    fn validate(node: &Node<i64, String>) {
        validate_height(node);
        validate_intervals(node);
        validate_string_metadata(node);
        if let Some(n) = node.left.as_ref() {
            validate(n)
        }
        if let Some(n) = node.right.as_ref() {
            validate(n)
        }
    }

    fn validate_height(node: &Node<i64, String>) {
        let left_height = node.left.as_ref().map_or(0, |n| n.height);
        let right_height = node.right.as_ref().map_or(0, |n| n.height);
        assert!((left_height - right_height).abs() <= 1);
        assert_eq!(node.height, cmp::max(left_height, right_height) + 1)
    }

    fn validate_intervals(node: &Node<i64, String>) {
        let mut reached_maximum: bool = false;
        if node.interval.end == node.max {
            reached_maximum = true;
        }
        if let Some(ref son) = node.left {
            assert!(
                son.max <= node.max,
                "left max invariant violated:\n{:?} -> {:?}",
                node,
                son
            );
            assert!(
                son.interval.start <= node.interval.start,
                "left ord invariant violated\n{:?} -> {:?}",
                node,
                son
            );
            if node.max == son.max {
                reached_maximum = true;
            }
        }
        if let Some(ref son) = node.right {
            assert!(
                son.max <= node.max,
                "right max invariant violated\n{:?} -> {:?}",
                node,
                son
            );
            assert!(
                son.interval.start >= node.interval.start,
                "right ord invariant violated\n{:?} -> {:?}",
                node,
                son
            );
            if node.max == son.max {
                reached_maximum = true;
            }
        }
        assert!(reached_maximum, "maximum invariant violated: {:?}", node);
    }

    fn validate_string_metadata(node: &Node<i64, String>) {
        let mut name: String = "".to_string();
        name.push_str(&node.interval.start.to_string());
        name.push(':');
        name.push_str(&node.interval.end.to_string());
        assert_eq!(name, node.value, "Invalid metadata for node {:?}", node);
    }

    fn insert_and_validate(tree: &mut IntervalTree<i64, String>, start: i64, end: i64) {
        let mut name: String = "".to_string();
        name.push_str(&start.to_string());
        name.push(':');
        name.push_str(&end.to_string());
        tree.insert(start..end, name);
        if let Some(ref n) = tree.root {
            validate(n);
        }
    }

    fn make_entry_tuples(intervals: Vec<Range<i64>>) -> Vec<(Interval<i64>, String)> {
        let mut entries = vec![];
        for interval in intervals {
            let mut data: String = "".to_string();
            data.push_str(&interval.start.to_string());
            data.push(':');
            data.push_str(&interval.end.to_string());
            entries.push((interval.into(), data));
        }
        entries.sort_by(|x1, x2| x1.1.cmp(&x2.1));
        entries
    }

    fn assert_intersections(
        tree: &IntervalTree<i64, String>,
        target: Range<i64>,
        expected_results: Vec<Range<i64>>,
    ) {
        let mut actual_entries: Vec<Entry<'_, i64, String>> = tree.find(&target).collect();
        println!("{:?}", actual_entries);
        actual_entries.sort_by(|x1, x2| x1.data.cmp(x2.data));
        let expected_entries = make_entry_tuples(expected_results);
        assert_eq!(actual_entries.len(), expected_entries.len());
        for (actual, expected) in actual_entries.iter().zip(expected_entries.iter()) {
            assert_eq!(*actual.interval, expected.0);
            assert_eq!(actual.data, &expected.1);
        }
    }

    fn assert_not_found(tree: &IntervalTree<i64, String>, target: Range<i64>) {
        assert_intersections(tree, target, vec![]);
    }

    // Clippy has a warning against `vec!` macros with a single range as argument.
    // Since we do actually want that here, we disable the warning.
    #[test]
    #[allow(clippy::single_range_in_vec_init)]
    fn test_insertion_and_intersection() {
        let mut tree: IntervalTree<i64, String> = IntervalTree::new();
        assert_eq!(tree.find(1..2).count(), 0);
        assert_eq!(tree.find_mut(1..2).count(), 0);
        tree.insert(50..51, "50:51".to_string());
        assert_not_found(&tree, 49..50);
        assert_intersections(&tree, 49..55, vec![50..51]);
        assert_not_found(&tree, 51..55);
        assert_not_found(&tree, 52..55);
        assert_not_found(&tree, 40..45);
        insert_and_validate(&mut tree, 80, 81);
        assert_intersections(&tree, 80..83, vec![80..81]);
        assert_intersections(&tree, 1..100, vec![50..51, 80..81]);
        assert_not_found(&tree, 82..83);
        insert_and_validate(&mut tree, 30, 35);
        assert_intersections(&tree, 25..33, vec![30..35]);
        assert_intersections(&tree, 1..100, vec![30..35, 50..51, 80..81]);
        assert_not_found(&tree, 42..43);
        assert_not_found(&tree, 35..36);
        assert_not_found(&tree, 22..29);
        insert_and_validate(&mut tree, 70, 77);
        assert_intersections(&tree, 75..79, vec![70..77]);
        assert_intersections(&tree, 1..100, vec![30..35, 50..51, 70..77, 80..81]);
        assert_not_found(&tree, 62..68);
        assert_intersections(&tree, 75..77, vec![70..77]);
        assert_not_found(&tree, 78..79);
        assert_intersections(&tree, 49..51, vec![50..51]);
        assert_intersections(&tree, 49..55, vec![50..51]);
        assert_not_found(&tree, 51..55);
        assert_not_found(&tree, 52..55);
        assert_not_found(&tree, 40..45);
        insert_and_validate(&mut tree, 101, 102);
        insert_and_validate(&mut tree, 103, 104);
        insert_and_validate(&mut tree, 105, 106);
        insert_and_validate(&mut tree, 107, 108);
        insert_and_validate(&mut tree, 111, 112);
        insert_and_validate(&mut tree, 113, 114);
        insert_and_validate(&mut tree, 115, 116);
        insert_and_validate(&mut tree, 117, 118);
        insert_and_validate(&mut tree, 119, 129);
        assert_not_found(&tree, 112..113);
        assert_not_found(&tree, 108..109);
        assert_intersections(&tree, 106..108, vec![107..108]);
        assert_intersections(&tree, 1..100, vec![30..35, 50..51, 70..77, 80..81]);
        assert_intersections(&tree, 1..101, vec![30..35, 50..51, 70..77, 80..81]);
        assert_intersections(
            &tree,
            1..102,
            vec![30..35, 50..51, 70..77, 80..81, 101..102],
        );
        assert_intersections(
            &tree,
            100..200,
            vec![
                101..102,
                103..104,
                105..106,
                107..108,
                111..112,
                113..114,
                115..116,
                117..118,
                119..129,
            ],
        );
    }

    #[test]
    fn test_insertion_and_intersection_2() {
        let mut tree: IntervalTree<i64, String> = IntervalTree::new();
        // interval size we'll insert into the tree
        let k = 10;
        for i in 100..200 {
            insert_and_validate(&mut tree, i, i + k);
        }
        for i in 90..210 {
            // "random" interval length we'll search against.
            // the exact formula doesn't matter; the point is
            // it will vary from 0.5 * k to 1.5 * k
            let l = k / 2 + i % k;
            let lower_bound = i;
            let upper_bound = i + l;
            let smallest_start = max(lower_bound - k + 1, 100);
            let largest_start = min(upper_bound, 200);
            let mut expected_intersections = vec![];
            for j in smallest_start..largest_start {
                expected_intersections.push(j..j + k);
            }
            assert_intersections(&tree, lower_bound..upper_bound, expected_intersections);
        }
    }

    #[test]
    fn zero_width_ranges() {
        let mut tree: IntervalTree<i64, String> = IntervalTree::new();
        tree.insert(10..10, "10:10".to_string());

        assert_not_found(&tree, 5..15);
        assert_not_found(&tree, 10..10);

        insert_and_validate(&mut tree, 50, 60);
        assert_not_found(&tree, 55..55);
    }

    #[test]
    fn from_iterator() {
        let tree: IntervalTree<i64, ()> = vec![(10..100, ()), (10..20, ()), (1..8, ())]
            .into_iter()
            .collect();
        assert_eq!(tree.find(&(0..1000)).count(), 3);
        let tree2: IntervalTree<_, _> = tree
            .find(&(11..30))
            .map(|e| (e.interval().clone(), *e.data()))
            .collect();
        assert_eq!(tree2.find(&(0..1000)).count(), 2);
    }

    #[test]
    fn iter_mut() {
        let mut tree: IntervalTree<i64, usize> = vec![(10..100, 0), (10..20, 0), (1..8, 0)]
            .into_iter()
            .collect();
        let q = Interval::new(11..30).unwrap();
        for mut e in tree.find_mut(q.clone()) {
            *e.data() += 1;
        }
        assert!(tree
            .find(0..100)
            .all(|e| if super::intersect(e.interval(), &q) {
                *e.data() == 1
            } else {
                *e.data() == 0
            }));
    }
}
