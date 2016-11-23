//! Interval tree, a data structure for efficiently storing and searching numeric intervals.
//!
//! This data structure uses the `std::ops:Range` type to define the intervals. The range bounds
//! may be specified by any type from the `num` crate that satisfies the `std::cmp::Ord` trait.
//! Upon inserting an interval may be associated with a data value. The interval data is stored in
//! an augmented AVL-tree which allows for efficient inserting and querying.
//!
//! # Example
//! ```
//! use bio::data_structures::interval_tree::IntervalTree;
//!
//! let mut tree = IntervalTree::new();
//! tree.insert((11..20), "Range_1").unwrap();
//! tree.insert((25..30), "Range_2").unwrap();
//! for r in tree.find(&(15..25)).unwrap() {
//!     assert_eq!(r.interval(), &(11..20));
//!     assert_eq!(r.data(), &"Range_1");
//! }
extern crate num;

use self::num::traits::Num;

use std::cmp;
use std::mem;
use std::fmt::Debug;
use std::ops::Range;

/// An interval tree for storing Ranges with data
#[derive(Debug, Clone)]
pub struct IntervalTree<N, D> {
    root: Option<Node<N, D>>,
}

/// An `Entry` is used by the IntervalTreeIterator to return references to the data in the tree
#[derive(PartialEq, Eq, Debug, Clone)]
pub struct Entry<'a, N: 'a, D: 'a> {
    data: &'a D,
    interval: &'a Range<N>,
}

/// A find query on the interval tree returns a `Iterator<
impl<'a, N: 'a, D: 'a> Entry<'a, N, D> {
    /// Get a reference to the data for this entry
    pub fn data(&self) -> &'a D {
        self.data
    }

    /// Get a reference to the interval for this entry
    pub fn interval(&self) -> &'a Range<N> {
        self.interval
    }
}

/// An IntervalTreeIterator is returned from by find and iterates over the entries
/// overlapping the query
pub struct IntervalTreeIterator<'a, N: 'a, D: 'a> {
    nodes: Vec<&'a Node<N, D>>,
    interval: Range<N>,
}

impl<'a, N: Debug + Num + Clone + Ord + 'a, D: Debug + 'a> Iterator
    for IntervalTreeIterator<'a, N, D> {
    type Item = Entry<'a, N, D>;

    fn next(&mut self) -> Option<Entry<'a, N, D>> {
        loop {
            let candidate = match self.nodes.pop() {
                None => return None,
                Some(node) => node,
            };

            // stop traversal if the query interval is beyond the current node and all children
            if self.interval.start < candidate.max {
                if let Some(ref left) = candidate.left {
                    self.nodes.push(left);
                }

                // don't traverse right if the query interval is completely before the current interval
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

impl<N: Debug + Num + Clone + Ord, D: Debug> IntervalTree<N, D> {
    /// Creates a new empty `IntervalTree`
    pub fn new() -> Self {
        IntervalTree { root: None }
    }

    /// Inserts a `Range` into the tree and associates it with `data`
    pub fn insert(&mut self, interval: Range<N>, data: D) -> Result<(), IntervalTreeError> {
        try!(validate(&interval));
        match self.root {
            Some(ref mut n) => n.insert(interval, data),
            None => self.root = Some(Node::new(interval, data)),
        };
        Ok(())
    }

    /// Uses the provided range to find overlapping intervals in the tree and returns an
    /// `IntervalTreeIterator`
    pub fn find(&self, interval: &Range<N>) -> Result<IntervalTreeIterator<N, D>, IntervalTreeError> {
        try!(validate(&interval));
        Ok(match self.root {
            Some(ref n) => n.find_iter(interval.clone()),
            None => {
                let empty_nodes = vec![];
                IntervalTreeIterator {
                    nodes: empty_nodes,
                    interval: interval.clone(),
                }
            }
        })
    }
}

fn validate<N: Ord + Debug>(interval: &Range<N>) -> Result<(), IntervalTreeError> {
    if interval.start > interval.end {
        Err(IntervalTreeError::InvalidRange)
    } else {
        Ok(())
    }
}

#[derive(Debug, Clone)]
struct Node<N, D> {
    // actual interval data
    interval: Range<N>,
    value: D,
    // tree metadata
    max: N,
    height: i64,
    left: Option<Box<Node<N, D>>>,
    right: Option<Box<Node<N, D>>>,
}

impl<N: Debug + Num + Clone + Ord, D: Debug> Node<N, D> {
    fn new(interval: Range<N>, data: D) -> Self {
        let max = interval.end.clone();
        Node {
            interval: interval,
            max: max,
            height: 1,
            value: data,
            left: None,
            right: None,
        }
    }

    fn insert(&mut self, interval: Range<N>, data: D) {
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

    fn find_iter(&self, interval: Range<N>) -> IntervalTreeIterator<N, D> {
        let nodes = vec![self];
        IntervalTreeIterator {
            nodes: nodes,
            interval: interval,
        }
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
                let mut right =
                    self.right.as_mut().expect("Invalid tree: leaf is taller than its sibling.");
                let right_left_h = right.left.as_ref().map_or(0, |n| n.height);
                let right_right_h = right.right.as_ref().map_or(0, |n| n.height);
                if right_left_h > right_right_h {
                    right.rotate_right();
                }
            }
            self.rotate_left();
        } else {
            {
                let mut left =
                    self.left.as_mut().expect("Invalid tree: leaf is taller than its sibling.");
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

fn swap_interval_data<N, D>(node_1: &mut Node<N, D>, node_2: &mut Node<N, D>) {
    mem::swap(&mut node_1.value, &mut node_2.value);
    mem::swap(&mut node_1.interval, &mut node_2.interval);
}

fn intersect<N: Ord>(range_1: &Range<N>, range_2: &Range<N>) -> bool {
    range_1.start < range_1.end && range_2.start < range_2.end &&
        range_1.end > range_2.start && range_1.start < range_2.end
}

quick_error! {
    #[derive(Debug)]
    pub enum IntervalTreeError {
        InvalidRange {
            description("Range for IntervalTree must have a positive width")
        }
    }
}



#[cfg(test)]
mod tests {
    use super::{Node, IntervalTree, Entry};
    use std::cmp;
    use std::cmp::{min, max};
    use std::ops::Range;

    fn validate(node: &Node<i64, String>) {
        validate_height(node);
        validate_intervals(node);
        validate_string_metadata(node);
        node.left.as_ref().map(|n| validate(n));
        node.right.as_ref().map(|n| validate(n));
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
            if !(son.max <= node.max) {
                panic!("left max invariant violated:\n{:?} -> {:?}", node, son);
            }
            if !(son.interval.start <= node.interval.start) {
                panic!("left ord invariant violated\n{:?} -> {:?}", node, son);
            }
            if node.max == son.max {
                reached_maximum = true;
            }
        }
        if let Some(ref son) = node.right {
            if !(son.max <= node.max) {
                panic!("right max invariant violated\n{:?} -> {:?}", node, son);
            }

            if !(son.interval.start >= node.interval.start) {
                panic!("right ord invariant violated\n{:?} -> {:?}", node, son);
            }
            if node.max == son.max {
                reached_maximum = true;
            }
        }
        if !reached_maximum {
            panic!("maximum invariant violated: {:?}", node);
        }
    }

    fn validate_string_metadata(node: &Node<i64, String>) {
        let mut name: String = "".to_string();
        name.push_str(&node.interval.start.to_string());
        name.push_str(":");
        name.push_str(&node.interval.end.to_string());
        if name != node.value {
            panic!("Invalid metadata for node {:?}", node);
        }
    }

    fn insert_and_validate(tree: &mut IntervalTree<i64, String>, start: i64, end: i64) {
        let mut name: String = "".to_string();
        name.push_str(&start.to_string());
        name.push_str(":");
        name.push_str(&end.to_string());
        tree.insert((start..end), name).unwrap();
        if let Some(ref n) = tree.root {
            validate(n);
        }
    }

    fn make_entry_tuples(intervals: Vec<Range<i64>>) -> Vec<(Range<i64>, String)> {
        let mut entries = vec![];
        for interval in intervals {
            let mut data: String = "".to_string();
            data.push_str(&interval.start.to_string());
            data.push_str(":");
            data.push_str(&interval.end.to_string());
            entries.push((interval, data));
        }
        entries.sort_by(|x1, x2| x1.1.cmp(&x2.1));
        entries
    }

    fn assert_intersections(tree: &IntervalTree<i64, String>,
                            target: Range<i64>,
                            expected_results: Vec<Range<i64>>) {
        let mut actual_entries: Vec<Entry<i64, String>> = tree.find(&target).unwrap().collect();
        println!("{:?}", actual_entries);
        actual_entries.sort_by(|x1, x2| x1.data.cmp(&x2.data));
        let expected_entries = make_entry_tuples(expected_results);
        assert_eq!(actual_entries.len(), expected_entries.len());
        for (actual, expected) in actual_entries.iter().zip(expected_entries.iter()) {
            assert_eq!(actual.interval, &expected.0);
            assert_eq!(actual.data, &expected.1);
        }
    }

    fn assert_not_found(tree: &IntervalTree<i64, String>, target: Range<i64>) {
        assert_intersections(tree, target, vec![]);
    }

    #[test]
    fn test_insertion_and_intersection() {
        let mut tree: IntervalTree<i64, String> = IntervalTree::new();
        tree.insert((50..51), "50:51".to_string()).unwrap();
        assert_not_found(&tree, (49..50));
        assert_intersections(&tree, (49..55), vec![(50..51)]);
        assert_not_found(&tree, (51..55));
        assert_not_found(&tree, (52..55));
        assert_not_found(&tree, (40..45));
        insert_and_validate(&mut tree, 80, 81);
        assert_intersections(&tree, (80..83), vec![(80..81)]);
        assert_intersections(&tree, (1..100), vec![(50..51), (80..81)]);
        assert_not_found(&tree, (82..83));
        insert_and_validate(&mut tree, 30, 35);
        assert_intersections(&tree, (25..33), vec![(30..35)]);
        assert_intersections(&tree, (1..100), vec![(30..35), (50..51), (80..81)]);
        assert_not_found(&tree, (42..43));
        assert_not_found(&tree, (35..36));
        assert_not_found(&tree, (22..29));
        insert_and_validate(&mut tree, 70, 77);
        assert_intersections(&tree, (75..79), vec![(70..77)]);
        assert_intersections(&tree,
                             (1..100),
                             vec![(30..35), (50..51), (70..77), (80..81)]);
        assert_not_found(&tree, (62..68));
        assert_intersections(&tree, (75..77), vec![(70..77)]);
        assert_not_found(&tree, (78..79));
        assert_intersections(&tree, (49..51), vec![(50..51)]);
        assert_intersections(&tree, (49..55), vec![(50..51)]);
        assert_not_found(&tree, (51..55));
        assert_not_found(&tree, (52..55));
        assert_not_found(&tree, (40..45));
        insert_and_validate(&mut tree, 101, 102);
        insert_and_validate(&mut tree, 103, 104);
        insert_and_validate(&mut tree, 105, 106);
        insert_and_validate(&mut tree, 107, 108);
        insert_and_validate(&mut tree, 111, 112);
        insert_and_validate(&mut tree, 113, 114);
        insert_and_validate(&mut tree, 115, 116);
        insert_and_validate(&mut tree, 117, 118);
        insert_and_validate(&mut tree, 119, 129);
        assert_not_found(&tree, (112..113));
        assert_not_found(&tree, (108..109));
        assert_intersections(&tree, (106..108), vec![(107..108)]);
        assert_intersections(&tree,
                             (1..100),
                             vec![(30..35), (50..51), (70..77), (80..81)]);
        assert_intersections(&tree,
                             (1..101),
                             vec![(30..35), (50..51), (70..77), (80..81)]);
        assert_intersections(&tree,
                             (1..102),
                             vec![(30..35), (50..51), (70..77), (80..81), (101..102)]);
        assert_intersections(&tree,
                             (100..200),
                             vec![(101..102), (103..104), (105..106), (107..108), (111..112),
                                  (113..114), (115..116), (117..118), (119..129)]);
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
            assert_intersections(&tree,
                                 (lower_bound..upper_bound),
                                 expected_intersections);
        }
    }

    #[test]
    fn zero_width_ranges() {
        let mut tree: IntervalTree<i64, String> = IntervalTree::new();
        tree.insert((10..10), "10:10".to_string()).unwrap();

        assert_not_found(&tree, (5..15));
        assert_not_found(&tree, (10..10));

        insert_and_validate(&mut tree, 50, 60);
        assert_not_found(&tree, (55..55));
    }
}
