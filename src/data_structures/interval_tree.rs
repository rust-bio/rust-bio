// Copyright 2016 Brett Bowman.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::usize;
use std::collections::BTreeSet;

use data_structures::interval::Interval;

///
/// Definition of the IntervalTreeInterface
///
pub trait IntervalTreeInterface {
    /// Create an empty interval tree
    fn new() -> Self;

    /// Create a new interval tree from some existing intervals
    fn from_slice(ints: &[Interval]) -> Self;

    /// Empty out the tree by emptying out it's underlying storage
    fn clear(&mut self) -> ();

    /// Return whether a value is contained by an interval within the tree
    fn contains(&self, value: usize) -> bool;

    /// Return whether this tree contains any data
    fn is_empty(&self) -> bool;

    /// Return a new interval tree containing the gaps within this tree
    fn gaps(&self) -> Self;

    /// Add a new interval into the tree
    fn insert(&mut self, other: Interval) -> ();

    /// Return the number of discrete elements in the tree
    fn len(&self) -> usize;

    /// Return the first position covered by this tree
    fn start(&self) -> Option<usize>;

    /// Return the last position covered by this tree
    fn stop(&self) -> Option<usize>;
}

///
/// A simple interval-tree, built on the BTreeSet from the standard library
///
#[cfg_attr(feature = "serde_macros", derive(Serialize, Deserialize))]
pub struct IntervalTree {
    storage: BTreeSet<Interval>,
}

///
/// Implementation of the IntervalTreeInterface for IntervalTree
///
impl IntervalTreeInterface for IntervalTree {
    /// Create an empty interval tree
    fn new() -> Self {
        let storage: BTreeSet<Interval> = BTreeSet::new();
        
        IntervalTree {
            storage: storage, 
        }
    }

    /// Create a new interval tree from some existing intervals
    fn from_slice(ints: &[Interval]) -> Self {
        let storage: BTreeSet<Interval> = ints.iter().cloned().collect();
        
        IntervalTree {
            storage: storage, 
        }
    }

    /// Empty out the tree by emptying out it's underlying storage
    fn clear(&mut self) -> () {
        self.storage.clear();
    }

    /// Return whether a value is contained by an interval within the tree
    fn contains(&self, value: usize) -> bool {
        for int in self.storage.iter() {
            if int.start > value {
                return false;
            } else if int.contains(value) {
                return true;
            }
        }
        return false;
    }

    /// Return whether this tree contains any data
    fn is_empty(&self) -> bool {
        self.storage.is_empty()
    }

    /// Return a new interval tree containing the gaps within this tree
    fn gaps(&self) -> IntervalTree {
        let mut leaves: Vec<Interval> = Vec::new();

        // If this tree is empty, return a single node covering all values
        if self.is_empty() { 
            leaves.push(Interval::new(0, usize::MAX));
            return IntervalTree::from_slice(&leaves); 
        } 
        
        // If this tree is complete on the other hand, return an empty tree
        if self.storage.len() == 1 {
            let leaf = self.storage.iter().next().unwrap();
            if leaf.start == 0 && leaf.stop == usize::MAX {
                return IntervalTree::new();
            }
        }

        // Add gap intervals to preceed each non-gap interval
        let mut prev_stop = 0;
        for int in self.storage.iter() {
            if int.start > prev_stop {
                leaves.push(Interval::new(prev_stop, int.start));
            }
            prev_stop = int.stop;
        }
        
        // Add the final gap that succeeds the final non-gap interval
        if prev_stop < usize::MAX {
            leaves.push(Interval::new(prev_stop, usize::MAX));
        }

        IntervalTree::from_slice(&leaves)
    }

    /// Add a new interval into the tree
    fn insert(&mut self, other: Interval) -> () {
        let mut new_int = other;
        let mut to_remove: Vec<Interval> = Vec::new();

        // Identify and merge any overlapping intervals
        for int in self.storage.iter() {
            if new_int < *int && new_int.stop < int.start {
                break
            } else if other.overlaps(int) {
                new_int = new_int.union(int).unwrap();
                to_remove.push(*int);
            }
        }

        // Remove the overlapping intervals and insert the new union
        for int in to_remove.iter() {
            self.storage.remove(int);
        }
        self.storage.insert(new_int);
    }

    /// Return the number of discrete elements in the tree
    fn len(&self) -> usize {
        self.storage.len()
    }

    /// Return the first position covered by this tree
    fn start(&self) -> Option<usize> {
        if !self.is_empty() {
            Some(self.storage.iter().next().unwrap().start)
        } else {
            None
        }
    }

    /// Return the last position covered by this tree
    fn stop(&self) -> Option<usize> {
        if !self.is_empty() {
            Some(self.storage.iter().last().unwrap().stop)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::usize;
    use data_structures::interval::Interval;

    #[test]
    fn test_clear() {
        let mut tree = IntervalTree::new();
        assert_eq!(0, tree.len());
        tree.insert(Interval::new(2, 4));
        assert_eq!(1, tree.len());
        tree.clear();
        assert_eq!(0, tree.len());
    }
    
    #[test]
    fn test_contains() {
        let mut tree = IntervalTree::new();
        tree.insert(Interval::new(5, 7));
        assert_eq!(false, tree.contains(4));
        assert_eq!(true,  tree.contains(6));
        assert_eq!(false, tree.contains(8));
    }
    
    #[test]
    fn test_is_empty() {
        let mut tree = IntervalTree::new();
        assert_eq!(true, tree.is_empty());
        tree.insert(Interval::new(2, 4));
        assert_eq!(false, tree.is_empty());
        tree.clear();
        assert_eq!(true, tree.is_empty());
    }

    #[test]
    fn test_gaps1() {
        let tree = IntervalTree::new();
        let gaps = tree.gaps();
        assert_eq!(1, gaps.len());
        assert_eq!(0, gaps.start().unwrap());
        assert_eq!(usize::MAX, gaps.stop().unwrap());
    }

    #[test]
    fn test_gaps2() {
        let mut tree = IntervalTree::new();
        tree.insert(Interval::new(0, usize::MAX));
        let gaps = tree.gaps();
        assert_eq!(true, gaps.is_empty());
    }

    #[test]
    fn test_gaps3() {
        let mut tree = IntervalTree::new();
        tree.insert(Interval::new(5, 10));
        let gaps = tree.gaps();
        assert_eq!(2, gaps.len());
        assert_eq!(0, gaps.start().unwrap());
        assert_eq!(usize::MAX, gaps.stop().unwrap());
    }

    #[test]
    fn test_merging1() {
        let mut tree = IntervalTree::new();
        tree.insert(Interval::new(2, 4));
        tree.insert(Interval::new(4, 6));
        assert_eq!(1, tree.len());
        assert_eq!(2, tree.start().unwrap());
        assert_eq!(6, tree.stop().unwrap());
    }

    #[test]
    fn test_merging2() {
        let mut tree = IntervalTree::new();
        tree.insert(Interval::new(1, 3));
        tree.insert(Interval::new(5, 7));
        tree.insert(Interval::new(9, 11));
        assert_eq!(3, tree.len());
        tree.insert(Interval::new(3, 9));
        assert_eq!(1, tree.len());
    }
}    
