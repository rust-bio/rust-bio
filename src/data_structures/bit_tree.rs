// Copyright 2014-2015 Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! BIT-tree (Binary Indexed Trees, aka Fenwick Tree) maintains a prefix-sum or
//! prefix-max that can be efficiently queried and updated. From: Peter M. Fenwick (1994). "A new data structure for cumulative frequency tables". Software: Practice and Experience. 24 (3): 327–336.
//! Implementation outlined here: https://www.topcoder.com/community/data-science/data-science-tutorials/binary-indexed-trees/
//!
//! Time Complexity: O(log n) where `n = tree.len()`.
//! Memory Complexity: O(n) where `n = tree.len()`.
//! # Example for a max bit tree
//!
//! ```
//! use bio::data_structures::bit_tree::*;
//!
//! let mut bit = MaxBitTree::new(10);
//! bit.set(0, (1,0));
//! bit.set(1, (0,1));
//! bit.set(2, (2,2));
//! bit.set(3, (4,3));
//!
//! assert_eq!(bit.get(0), (1, 0));
//! assert_eq!(bit.get(1), (1, 0));
//! assert_eq!(bit.get(2), (2, 2));
//! assert_eq!(bit.get(3), (4, 3));
//! assert_eq!(bit.get(4), (4, 3));


use std::cmp::max;
use std::marker::PhantomData;
use std::ops::Add;

/// Fenwick tree prefix operator
pub trait PrefixOp<T> {
    fn operation(t1: T, t2: T) -> T;
}

/// In a max bit tree, get(i) will return the largest element e that has been added
/// to the bit tree with set(j, e), where j <= i. Initially all positions have
/// the value T::default(). Note that a set cannot be 'undone' by inserting
/// a smaller element at the same index.
/// Time Complexity: O(n) to build a new tree or O(log n) for get() and set() operations,
/// where `n = tree.len()`.

pub struct FenwickTree<T: Default + Ord, Op: PrefixOp<T>> {
    tree: Vec<T>,
    phantom: PhantomData<Op>,
}

impl<T: Ord + Default + Copy, Op: PrefixOp<T>> FenwickTree<T, Op> {
    /// Create a new bit tree with len elements
    pub fn new(len: usize) -> FenwickTree<T, Op> {
        let mut tree = Vec::new();

        // Pad length by one. The first element is unused.
        // Done this way to make the tree structure work correctly.
        for _ in 0..(len + 2) {
            tree.push(T::default());
        }

        FenwickTree {
            tree,
            phantom: PhantomData,
        }
    }

    /// Max bit tree: get(i) returns the largest element e that has been added
    /// to the bit tree with set(j, e), where j <= i.
    pub fn get(&self, idx: usize) -> T {
        let mut idx = idx + 1;
        let mut sum = T::default();
        while idx > 0 {
            sum = Op::operation(sum, self.tree[idx]);
            idx -= (idx as isize & -(idx as isize)) as usize;
        }

        sum
    }

    /// Set the value val at position idx. In max bit trees, val will
    /// be returned for any get(j) where j >= idx, if
    /// it is the maximum value inserted between 0 and j.
    /// Inserting a value val2 after inserting val1 where val1 > val2
    /// will have no effect.
    pub fn set(&mut self, idx: usize, val: T) {
        let mut idx = idx + 1;
        while idx < self.tree.len() {
            self.tree[idx] = Op::operation(self.tree[idx], val);
            idx += (idx as isize & -(idx as isize)) as usize;
        }
    }
}

pub struct MaxOp;
impl<T: Copy + Ord> PrefixOp<T> for MaxOp {
    fn operation(t1: T, t2: T) -> T {
        max(t1, t2)
    }
}

/// Fenwick tree specialized for prefix-max
pub type MaxBitTree<T> = FenwickTree<T, MaxOp>;

pub struct SumOp;
impl<T: Copy + Add> PrefixOp<T> for SumOp
where
    T: Add<Output = T>,
{
    fn operation(t1: T, t2: T) -> T {
        t1 + t2
    }
}

/// Fenwick tree specialized for prefix-sum
pub type SumBitTree<T> = FenwickTree<T, SumOp>;

#[cfg(test)]
mod test_bit_tree {
    use super::MaxBitTree;

    #[test]
    pub fn test_bit_tree() {
        let mut bit = MaxBitTree::new(10);

        bit.set(0, (1, 0));
        bit.set(1, (1, 1));
        bit.set(2, (2, 2));
        bit.set(3, (3, 3));
        bit.set(4, (2, 4));
        bit.set(5, (2, 5));
        bit.set(6, (4, 6));
        bit.set(7, (5, 7));

        assert_eq!(bit.get(0), (1, 0));
        assert_eq!(bit.get(1), (1, 1));
        assert_eq!(bit.get(2), (2, 2));
        assert_eq!(bit.get(3), (3, 3));
        assert_eq!(bit.get(4), (3, 3));
        assert_eq!(bit.get(5), (3, 3));
        assert_eq!(bit.get(6), (4, 6));
        assert_eq!(bit.get(7), (5, 7));
    }
}
