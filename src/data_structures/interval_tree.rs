extern crate num;

use self::num::traits::Num;

use std::cmp;
use std::mem;
use std::fmt::Debug;
use std::ops::Range;

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct IntervalTree<N, D> {
    root: Option<Node<N, D>>,
}

impl<N: Debug + Num + Clone + Ord, D: Debug> IntervalTree<N, D> {
    pub fn new() -> Self {
        IntervalTree { root: None }
    }

    pub fn insert(&mut self, interval: Range<N>, data: D) {
        validate(&interval);
        match self.root {
            Some(ref mut n) => n.insert(interval, data),
            None => self.root = Some(Node::new(interval, data)),
        }
    }

    pub fn find_any(&self, interval: &Range<N>) -> Option<(&Range<N>, &D)> {
        validate(&interval);
        match self.root {
            Some(ref n) => {
                let mut result = vec![];
                n.find(interval, true, &mut result);
                result.pop()
            }
            None => None,
        }
    }

    pub fn find_all(&self, interval: &Range<N>) -> Vec<(&Range<N>, &D)> {
        validate(&interval);
        match self.root {
            Some(ref n) => {
                let mut result = vec![];
                n.find(interval, true, &mut result);
                result
            }
            None => vec![],
        }
    }
}

fn validate<N: Ord + Debug>(interval: &Range<N>) {
    if interval.start >= interval.end {
        panic!("Expected interval.end > interval.start, got: ({:?})",
               interval);
    }
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct Node<N, D> {
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
    pub fn new(interval: Range<N>, data: D) -> Self {
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

    pub fn insert(&mut self, interval: Range<N>, data: D) {
        if interval.start <= self.interval.start {
            if let Some(ref mut son) = self.left {
                son.insert(interval, data);
            } else {
                self.left = Some(Box::new(Node::new(interval, data)));
            }
        } else {
            if let Some(ref mut son) = self.right {
                son.insert(interval, data);
            } else {
                self.right = Some(Box::new(Node::new(interval, data)));
            }
        }
        self.repair();
    }

    pub fn find<'a>(&'a self,
                    interval: &Range<N>,
                    only_one: bool,
                    result: &mut Vec<(&'a Range<N>, &'a D)>) {
        if interval.end <= self.interval.start || self.interval.end <= interval.start {
            // no overlap with current node, recur into children
            if let Some(ref left) = self.left {
                if left.max > interval.start {
                    left.find(interval, only_one, result);
                    return;
                }
            }
            if let Some(ref right) = self.right {
                right.find(interval, only_one, result);
                return;
            }
        } else {
            result.push((&self.interval, &self.value));
            if only_one {
                return;
            } else {
                if let Some(ref left) = self.left {
                    left.find(interval, only_one, result);
                }
                if let Some(ref right) = self.right {
                    right.find(interval, only_one, result);
                }
            }
        }
    }

    fn update_height(&mut self) {
        let ref left_h = self.left.as_ref().map_or(0, |n| n.height);
        let ref right_h = self.right.as_ref().map_or(0, |n| n.height);
        self.height = 1 + cmp::max(*left_h, *right_h);
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
        let ref left_h = self.left.as_ref().map_or(0, |n| n.height);
        let ref right_h = self.right.as_ref().map_or(0, |n| n.height);
        // each case - update both height and max
        if (left_h - right_h).abs() <= 1 {
            self.update_height();
            self.update_max();
        } else if right_h > left_h {
            {
                let mut right =
                    self.right.as_mut().expect("Invalid tree: leaf is taller than its sibling.");
                let ref right_left_h = right.left.as_ref().map_or(0, |n| n.height);
                let ref right_right_h = right.right.as_ref().map_or(0, |n| n.height);
                if right_left_h > right_right_h {
                    right.rotate_right();
                }
            }
            self.rotate_left();
        } else {
            {
                let mut left =
                    self.left.as_mut().expect("Invalid tree: leaf is taller than its sibling.");
                let ref left_right_h = left.right.as_ref().map_or(0, |n| n.height);
                let ref left_left_h = left.left.as_ref().map_or(0, |n| n.height);
                if left_right_h > left_left_h {
                    println!("double rot");
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp;
    use std::ops::Range;

    fn validate(node: &Node<i64, String>) {
        validate_height(node);
        validate_intervals(node);
        validate_string_metadata(node);
        node.left.as_ref().map(|n| validate(n));
        node.right.as_ref().map(|n| validate(n));
    }

    fn validate_height(node: &Node<i64, String>) {
        let ref left_height = node.left.as_ref().map_or(0, |n| n.height);
        let ref right_height = node.right.as_ref().map_or(0, |n| n.height);
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
        tree.insert((start..end), name);
        if let Some(ref n) = tree.root {
            validate(n);
        }
    }

    fn assert_intersection(tree: &IntervalTree<i64, String>,
                           target: Range<i64>,
                           expected_bounds: Range<i64>) {
        let mut expected_data: String = "".to_string();
        expected_data.push_str(&expected_bounds.start.to_string());
        expected_data.push_str(":");
        expected_data.push_str(&expected_bounds.end.to_string());
        if let Some((actual_bounds, actual_data)) = tree.find_any(&target) {
            assert_eq!(&expected_data.to_string(), actual_data);
            assert_eq!(&expected_bounds, actual_bounds);
        } else {
            panic!("Expected to find {:?} in the tree", expected_bounds);
        }
    }


    #[test]
    fn test_insertion_and_intersection() {
        let mut tree: IntervalTree<i64, String> = IntervalTree::new();
        tree.insert((50..51), "50:51".to_string());
        assert_eq!(tree.find_any(&(49..50)), None);
        assert_intersection(&tree, (49..55), (50..51));
        assert_eq!(tree.find_any(&(51..55)), None);
        assert_eq!(tree.find_any(&(52..55)), None);
        assert_eq!(tree.find_any(&(40..45)), None);
        insert_and_validate(&mut tree, 80, 81);
        assert_intersection(&tree, (80..83), (80..81));
        assert_eq!(tree.find_any(&(82..83)), None);
        insert_and_validate(&mut tree, 30, 35);
        assert_intersection(&tree, (25..33), (30..35));
        assert_eq!(tree.find_any(&(42..43)), None);
        assert_eq!(tree.find_any(&(35..36)), None);
        assert_eq!(tree.find_any(&(22..29)), None);
        insert_and_validate(&mut tree, 70, 77);
        assert_intersection(&tree, (75..79), (70..77));
        assert_eq!(tree.find_any(&(62..68)), None);
        assert_intersection(&tree, (75..77), (70..77));
        assert_eq!(tree.find_any(&(78..79)), None);
        assert_intersection(&tree, (49..51), (50..51));
        assert_intersection(&tree, (49..55), (50..51));
        assert_eq!(tree.find_any(&(51..55)), None);
        assert_eq!(tree.find_any(&(52..55)), None);
        assert_eq!(tree.find_any(&(40..45)), None);
        insert_and_validate(&mut tree, 101, 102);
        insert_and_validate(&mut tree, 103, 104);
        insert_and_validate(&mut tree, 105, 106);
        insert_and_validate(&mut tree, 107, 108);
        insert_and_validate(&mut tree, 111, 112);
        insert_and_validate(&mut tree, 113, 114);
        insert_and_validate(&mut tree, 115, 116);
        insert_and_validate(&mut tree, 117, 118);
        insert_and_validate(&mut tree, 119, 129);
        assert_eq!(tree.find_any(&(112..113)), None);
        assert_eq!(tree.find_any(&(108..109)), None);
        assert_intersection(&tree, (106..108), (107..108));
    }
}
