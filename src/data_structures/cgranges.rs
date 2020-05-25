use std::cmp::min;

use itertools::Itertools;

use crate::utils::Interval;

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct Entry<N: Ord + Clone + Copy, D> {
    data: D,
    interval: Interval<N>,
    max: N,
}

impl<N: Ord + Clone + Copy, D> Entry<N, D> {
    /// Get a reference to the data for this entry
    pub fn data(&self) -> &D {
        &self.data
    }

    /// Get a reference to the interval for this entry
    pub fn interval(&self) -> &Interval<N> {
        &self.interval
    }
}

impl<N: Ord + Clone + Copy, D> Default for IITree<N, D> {
    fn default() -> Self {
        IITree {
            entries: vec![],
            max_level: 0,
            indexed: false,
        }
    }
}

pub struct IITree<N: Ord + Clone + Copy, D> {
    entries: Vec<Entry<N, D>>,
    max_level: usize,
    indexed: bool,
}

impl<N: Ord + Clone + Copy, D: Clone> IITree<N, D> {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn from<I: Iterator<Item = (V, D)>, V>(iter: I) -> Self
    where
        V: Into<Interval<N>>,
    {
        let mut tree = Self::new();
        iter.for_each(|(interval, data)| tree.insert(interval, data));
        tree.index();
        tree
    }

    pub fn from_entries(entries: &[(Interval<N>, D)]) -> Self {
        let mut tree = Self {
            entries: entries
                .iter()
                .map(|(interval, data)| Entry {
                    interval: interval.clone(),
                    data: data.clone(),
                    max: interval.end,
                })
                .collect_vec(),
            max_level: 0,
            indexed: false,
        };
        tree.index();
        tree
    }

    pub fn insert<I: Into<Interval<N>>>(&mut self, interval: I, data: D) {
        let interval = interval.into();
        let max = interval.end;
        self.entries.push(Entry {
            interval,
            data,
            max: max,
        });
        self.indexed = false;
    }

    pub fn index(&mut self) {
        if !self.indexed {
            self.entries.sort_by_key(|e| e.interval.start);
            self.index_core();
            self.indexed = true;
        }
    }

    fn index_core(&mut self) {
        let a = &mut self.entries;
        if a.is_empty() {
            return;
        }

        let n = a.len();
        let mut last_i = 0;
        let mut last_value = a[0].max;
        (0..n).step_by(2).for_each(|i| {
            last_i = i;
            a[i].max = a[i].interval.end;
            last_value = a[i].max;
        });
        let mut k = 1;
        while (1 << k) < n {
            // process internal nodes in the bottom-up order
            let x = 1 << (k - 1);
            let i0 = (x << 1) - 1; // i0 is the first node
            let step = x << 2;
            for i in (i0..n).step_by(step) {
                // traverse all nodes at level k
                let end_left = a[i - x].max; // max value of the left child
                let end_right = if i + x < n { a[i + x].max } else { last_value }; // max value of the right child
                let end = max3(a[i].interval.end, end_left, end_right);
                a[i].max = end;
            }
            last_i = if (last_i >> k & 1) > 0 {
                last_i - x
            } else {
                last_i + x
            };
            if last_i < n && a[last_i].max > last_value {
                last_value = a[last_i].max
            }
            k += 1;
        }
        self.max_level = k - 1;
    }

    pub fn overlap(&self, start: N, end: N) -> Vec<&Entry<N, D>> {
        let mut buf = Vec::with_capacity(512);
        self.overlap_into(start, end, &mut buf);
        buf
    }

    pub fn overlap_into<'b, 'a: 'b>(&'a self, start: N, end: N, out: &mut Vec<&'b Entry<N, D>>) {
        if !self.indexed {
            panic!("This IITree has not been indexed yet. Call `index()` first.")
        }
        let n = self.entries.len() as usize;
        let a = &self.entries;
        out.clear();
        let mut stack = [StackCell::empty(); 64];
        // push the root; this is a top down traversal
        stack[0].k = self.max_level;
        stack[0].x = (1 << self.max_level) - 1;
        stack[0].w = false;
        let mut t = 1;
        while t > 0 {
            t -= 1;
            let StackCell { k, x, w } = stack[t];
            if k <= 3 {
                // we are in a small subtree; traverse every node in this subtree
                let i0 = x >> k << k;
                let i1 = min(i0 + (1 << (k + 1)) - 1, n);
                for i in i0..i1 {
                    if a[i].interval.start >= end {
                        break;
                    }
                    if start < a[i].interval.end {
                        // if overlap, append to out[]
                        out.push(&self.entries[i]);
                    }
                }
            } else if w == false {
                // if left child not processed
                let y = x - (1 << (k - 1)); // the left child of x; NB: y may be out of range (i.e. y>=n)
                stack[t].k = k;
                stack[t].x = x;
                stack[t].w = true; // re-add node x, but mark the left child having been processed
                t += 1;
                if y >= n || a[y].max > start {
                    // push the left child if y is out of range or may overlap with the query
                    stack[t].k = k - 1;
                    stack[t].x = y;
                    stack[t].w = false;
                    t += 1;
                }
            } else if x < n && a[x].interval.start < end {
                // need to push the right child
                if start < a[x].interval.end {
                    out.push(&self.entries[x]);
                }
                stack[t].k = k - 1;
                stack[t].x = x + (1 << (k - 1));
                stack[t].w = false;
                t += 1;
            }
        }
    }
}

fn max3<T: Ord>(a: T, b: T, c: T) -> T {
    a.max(b.max(c))
}

#[derive(Clone, Copy)]
struct StackCell {
    // node
    x: usize,
    // level
    k: usize,
    // false if left child hasn't been processed
    w: bool,
}

impl StackCell {
    fn empty() -> Self {
        Self {
            x: 0,
            k: 0,
            w: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_example() {
        let mut tree = IITree::new();
        tree.insert(12..34, 0);
        tree.insert(0..23, 1);
        tree.insert(34..56, 2);
        tree.index();
        let overlap = tree.overlap(22, 25);

        let e1 = Entry {
            interval: (0..23).into(),
            data: 1,
            max: 23,
        };
        let e2 = Entry {
            interval: (12..34).into(),
            data: 0,
            max: 56,
        };
        let expected = vec![&e1, &e2];
        assert_eq!(overlap, expected);
    }
}
