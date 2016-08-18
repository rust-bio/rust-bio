extern crate num;

use self::num::traits::Num;

use std::cmp;
use std::mem;
use std::fmt::Debug;

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct IntervalTree<N, D> {
    root: Option<Node<N, D>>,
}

impl<N: Debug + Num + Clone + Ord, D: Debug> IntervalTree<N, D> {
    pub fn new() -> Self {
        IntervalTree { root: None }
    }

    pub fn insert(&mut self, start: N, end: N, data: D) {
        validate_interval(&start, &end);
        match self.root {
            Some(ref mut n) => n.insert(start, end, data),
            None => self.root = Some(Node::new(start, end, data)),
        }
    }

    pub fn find_intersection(&self, start: N, end: N) -> Option<&D> {
        validate_interval(&start, &end);
        match self.root {
            Some(ref n) => n.find_intersection(start, end),
            None => None,
        }
    }
}


fn validate_interval<N: Debug + Num + Ord>(start: &N, end: &N) {
    assert!(start <= end, "Invalid interval: [{:?}, {:?}]", start, end);
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct Node<N, D> {
    // actual interval data
    start: N,
    end: N,
    value: D,
    // tree metadata
    max: N,
    height: i64,
    left: Option<Box<Node<N, D>>>,
    right: Option<Box<Node<N, D>>>,
}

impl<N: Debug + Num + Clone + Ord, D: Debug> Node<N, D> {
    pub fn new(start: N, end: N, data: D) -> Self {
        Node {
            start: start,
            end: end.clone(),
            max: end,
            height: 1,
            value: data,
            left: None,
            right: None,
        }
    }

    pub fn insert(&mut self, start: N, end: N, data: D) {
        if start <= self.start {
            if let Some(ref mut son) = self.left {
                son.insert(start, end, data);
            } else {
                self.left = Some(Box::new(Node::new(start, end, data)));
            }
        } else {
            if let Some(ref mut son) = self.right {
                son.insert(start, end, data);
            } else {
                self.right = Some(Box::new(Node::new(start, end, data)));
            }
        }
        self.repair();
    }

    pub fn find_intersection(&self, start: N, end: N) -> Option<&D> {
        if end < self.start || self.end < start {
            // no overlap with current node, recur into children
            if let Some(ref left) = self.left {
                if left.max >= start {
                    return left.find_intersection(start, end);
                }
            }
            if let Some(ref right) = self.right {
                return right.find_intersection(start, end);
            }
        } else {
            // overlaps current node
            return Some(&self.value);
        }
        return None;
    }

    fn update_height(&mut self) {
        let ref left_h = self.left.as_ref().map_or(0, |n| n.height);
        let ref right_h = self.right.as_ref().map_or(0, |n| n.height);
        self.height = 1 + cmp::max(*left_h, *right_h);
    }

    fn update_max(&mut self) {
        self.max = self.end.clone();
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
            println!("will fix right tree");
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
            println!("will fix left tree");
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
    mem::swap(&mut node_1.start, &mut node_2.start);
    mem::swap(&mut node_1.end, &mut node_2.end);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp;

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
        if node.end == node.max {
            reached_maximum = true;
        }
        if let Some(ref son) = node.left {
            if !(son.max <= node.max) {
                panic!("left max invariant violated:\n{:?} -> {:?}", node, son);
            }
            if !(son.start <= node.start) {
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

            if !(son.start >= node.start) {
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
        name.push_str(&node.start.to_string());
        name.push_str(":");
        name.push_str(&node.end.to_string());
        if name != node.value {
            panic!("Invalid metadata for node {:?}", node);
        }
    }

    fn insert_and_validate(tree: &mut IntervalTree<i64, String>, start: i64, end: i64) {
        let mut name: String = "".to_string();
        name.push_str(&start.to_string());
        name.push_str(":");
        name.push_str(&end.to_string());
        tree.insert(start, end, name);
        if let Some(ref n) = tree.root {
            validate(n);
        }
    }

    #[test]
    fn test_insertion_and_intersection() {
        let mut tree: IntervalTree<i64, String> = IntervalTree::new();
        tree.insert(50, 51, "50:51".to_string());
        assert_eq!(tree.find_intersection(49, 50), Some(&"50:51".to_string()));
        assert_eq!(tree.find_intersection(49, 55), Some(&"50:51".to_string()));
        assert_eq!(tree.find_intersection(51, 55), Some(&"50:51".to_string()));
        assert_eq!(tree.find_intersection(52, 55), None);
        assert_eq!(tree.find_intersection(40, 45), None);
        insert_and_validate(&mut tree, 80, 81);
        assert_eq!(tree.find_intersection(81, 83), Some(&"80:81".to_string()));
        assert_eq!(tree.find_intersection(82, 83), None);
        insert_and_validate(&mut tree, 30, 35);
        assert_eq!(tree.find_intersection(25, 33), Some(&"30:35".to_string()));
        assert_eq!(tree.find_intersection(42, 43), None);
        assert_eq!(tree.find_intersection(35, 36), Some(&"30:35".to_string()));
        assert_eq!(tree.find_intersection(22, 29), None);
        insert_and_validate(&mut tree, 70, 77);
        assert_eq!(tree.find_intersection(75, 79), Some(&"70:77".to_string()));
        assert_eq!(tree.find_intersection(62, 68), None);
        assert_eq!(tree.find_intersection(75, 77), Some(&"70:77".to_string()));
        assert_eq!(tree.find_intersection(78, 79), None);
        assert_eq!(tree.find_intersection(49, 50), Some(&"50:51".to_string()));
        assert_eq!(tree.find_intersection(49, 55), Some(&"50:51".to_string()));
        assert_eq!(tree.find_intersection(51, 55), Some(&"50:51".to_string()));
        assert_eq!(tree.find_intersection(52, 55), None);
        assert_eq!(tree.find_intersection(40, 45), None);
        insert_and_validate(&mut tree, 101, 101);
        insert_and_validate(&mut tree, 103, 103);
        insert_and_validate(&mut tree, 105, 105);
        insert_and_validate(&mut tree, 107, 107);
        insert_and_validate(&mut tree, 109, 109);
        insert_and_validate(&mut tree, 111, 111);
        insert_and_validate(&mut tree, 113, 113);
        insert_and_validate(&mut tree, 115, 115);
        insert_and_validate(&mut tree, 117, 117);
        insert_and_validate(&mut tree, 119, 119);
        assert_eq!(tree.find_intersection(112, 112), None);
        assert_eq!(tree.find_intersection(108, 108), None);
        assert_eq!(tree.find_intersection(106, 108), Some(&"107:107".to_string()));
    }
}
