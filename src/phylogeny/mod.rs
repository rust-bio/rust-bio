//! This module contains algorithms related to phylogenetic trees.
//!
//! Currently, it contains two methods for phylogeny reconstruction:
//! - [UPGMA](https://en.wikipedia.org/wiki/UPGMA)
//! - [Neighbor Joining](https://en.wikipedia.org/wiki/Neighbor_joining)
//!
//! It also contains an implementation to compute the [Robinson-Foulds
//! metric](https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric).
//!
//! Each of these algorithms is also available as a command line tool in [`rust-bio-tools`](https://github.com/rust-bio/rust-bio-tools).

use std::{
    cmp::{max, min, Reverse},
    collections::{BTreeSet, BinaryHeap, HashSet},
    mem::swap,
};

use bio_types::{distancematrix::DistanceMatrix, phylogeny::Tree};
use itertools::Itertools;
use ordered_float::NotNan;
use petgraph::{
    graph::NodeIndex,
    visit::EdgeRef,
    EdgeDirection::{Incoming, Outgoing},
};

/// Helper function to join two nodes in a new node with the given distances.
fn join(t: &mut Tree, u: NodeIndex, u_dist: f32, v: NodeIndex, v_dist: f32) -> NodeIndex {
    let m = t.g.add_node("".into());
    t.g.add_edge(m, u, u_dist);
    t.g.add_edge(m, v, v_dist);
    m
}

/// Run the UPGMA algorithm.
///
/// Compute the rooted phylogeny given a `DistanceMatrix`.
pub fn upgma(distances: DistanceMatrix) -> Tree {
    let DistanceMatrix {
        names,
        mut distances,
        matrix_type: _,
    } = distances;
    let mut t = Tree::new();
    let n = distances.len();
    // Tuples of (subtree root, number of leafs in subtree, distance to leafs)
    let mut parts: Vec<Option<(NodeIndex, usize, f32)>> = names
        .into_iter()
        .map(|name| (t.g.add_node(name), 1, 0.).into())
        .collect();
    let mut heap = BinaryHeap::with_capacity(n * n);
    for (i, ds) in distances.iter().enumerate() {
        for (j, &d) in ds.iter().enumerate() {
            if i == j {
                continue;
            }
            heap.push(Reverse((NotNan::new(d).unwrap(), i, j)));
        }
    }

    while let Some(Reverse((d, i, j))) = heap.pop() {
        assert!(i != j);
        if parts[i].is_none() || parts[j].is_none() || distances[i][j] != *d {
            continue;
        }
        let (i, j) = (min(i, j), max(i, j));
        if let (Some((pi, si, di)), Some((pj, sj, dj))) = (parts[i].take(), parts[j].take()) {
            // Merge phylogeny i and j, adding theirs sizes. Store the result in the lower of the two indices.
            parts[i] = Some((
                join(&mut t, pi, *d / 2. - di, pj, *d / 2. - dj),
                si + sj,
                *d / 2.,
            ));
            // Update the distances to other nodes
            for k in 0..n {
                if k == i || k == j {
                    continue;
                }
                let dk =
                    (si as f32 * distances[i][k] + sj as f32 * distances[j][k]) / (si + sj) as f32;
                distances[i][k] = dk;
                distances[k][i] = dk;
                heap.push(Reverse((NotNan::new(dk).unwrap(), k, i)));
            }
        }
    }

    t
}

/// Run the Neighbor-Joining algorithm
///
/// Computes an unrooted phylogony given a `DistanceMatrix`.
///
/// The root vertex is arbitrarily chosen as the middle in between the last two remaining vertices.
pub fn neighbor_joining(distances: DistanceMatrix) -> Tree {
    let DistanceMatrix {
        names,
        mut distances,
        matrix_type: _,
    } = distances;
    let mut t = Tree::new();

    // Start with one node for each taxon.
    let mut parts: Vec<Option<NodeIndex>> = names
        .into_iter()
        .map(|name| Some(t.g.add_node(name)))
        .collect();

    // A list of 'active' indices, indicating roots of the remaining subtrees.
    let mut active: Vec<usize> = (0..distances.len()).collect();

    while active.len() > 2 {
        // The sum of all distances from index i to any other index.
        let sum_d = |i: usize| -> f32 { active.iter().map(|&k| distances[i][k]).sum::<f32>() };

        // Find i,j for which Q(i,j) is minimal.
        let q = |&(&i, &j): &(&usize, &usize)| -> NotNan<f32> {
            let r = NotNan::new((active.len() - 2) as f32 * distances[i][j] - sum_d(i) - sum_d(j))
                .unwrap();
            r
        };

        // Find minimal distance pair.
        let (&i, &j) = active
            .iter()
            .cartesian_product(active.iter())
            .filter(|&(&i, &j)| i != j)
            .min_by_key(q)
            .unwrap();

        let (i, j) = (min(i, j), max(i, j));

        // Compute distance from merged vertex to the nodes being merged.
        let di = distances[i][j] / 2. + (sum_d(i) - sum_d(j)) / (2. * (active.len() as f32 - 2.));
        let dj = distances[i][j] - di;

        // Remove j from positions considered in later iterations.
        active.remove(active.iter().position(|&x| x == j).unwrap());

        // Compute all other distances.
        active.iter().filter(|&&k| k != i).for_each(|&k| {
            let dk = (distances[i][k] + distances[j][k] - distances[i][j]) / 2.;
            distances[i][k] = dk;
            distances[k][i] = dk;
        });
        parts[i] = Some(join(
            &mut t,
            parts[i].take().unwrap(),
            di,
            parts[j].take().unwrap(),
            dj,
        ));
    }

    // Merge the two remaining vertices with the given distance.
    if let [i, j] = active[..] {
        let d = distances[i][j] / 2.;
        parts[i] = Some(join(
            &mut t,
            parts[i].take().unwrap(),
            d,
            parts[j].take().unwrap(),
            d,
        ));
    }

    t
}

/// Compute the Robingson-Foulds metric.
///
/// Also called the RF-distance, a distance metric between two phylogenetic trees.
///
/// Leaf nodes *must* be named. Names for internal nodes are ignored.
/// Panics when the two trees have different sets of leafs.
///
/// The root is never considered as a leaf node.
// The current implementation is relatively simple, and could probably be optimized further.
//
// It keeps a HashSet of all partitions, and each partition is stored as a sorted vector.
// The DFS keeps a BTreeSet to efficiently insert and merge sets of leafs.
pub fn robinson_foulds_distance(p: &Tree, q: &Tree) -> usize {
    type Partition = BTreeSet<String>;
    type Partitions = HashSet<Vec<String>>;

    let leafs_p = p.g.externals(Outgoing).map(|u| p.g[u].clone()).collect();
    let leafs_q = q.g.externals(Outgoing).map(|u| q.g[u].clone()).collect();
    assert_eq!(leafs_p, leafs_q);

    // Insert the set of names in a subtree to the hashset of all partitions seen so far.
    // Also inserts the complement of this subtree into the hashset.
    fn insert_part_and_complement(
        subtree: &BTreeSet<String>,
        universe: &BTreeSet<String>,
        parts: &mut Partitions,
    ) {
        parts.insert(subtree.iter().cloned().collect());
        parts.insert(universe.difference(subtree).cloned().collect());
    }

    // For each node, recurses into its children to collect all leafs in the
    // subtree. Then adds this set of leafs and the complement into `parts`.
    fn dfs(t: &Tree, node: NodeIndex, universe: &Partition, parts: &mut Partitions) -> Partition {
        let mut subtree = Partition::new();
        let mut iter = t.g.edges_directed(node, Outgoing).peekable();
        if iter.peek().is_none() {
            // For leafs, insert their name into the set of values seen so far.
            subtree.insert(t.g[node].to_string());
        } else {
            // For internal nodes, recurse and merge all leafs.
            for edge in iter {
                let mut childtree = dfs(t, edge.target(), universe, parts);
                // Merge the smaller of childtree and subtree into the larger
                if childtree.len() > subtree.len() {
                    swap(&mut subtree, &mut childtree);
                }
                for x in childtree.into_iter() {
                    subtree.insert(x);
                }
            }
        }
        insert_part_and_complement(&subtree, universe, parts);
        return subtree;
    }
    let (mut partitions_p, mut partitions_q) = (Partitions::new(), Partitions::new());
    // Find the partitions of p and q by running a DFS from all roots.
    p.g.externals(Incoming).for_each(|r| {
        dfs(p, r, &leafs_p, &mut partitions_p);
    });
    q.g.externals(Incoming).for_each(|r| {
        dfs(q, r, &leafs_q, &mut partitions_q);
    });
    if leafs_p.len() > 1 {
        assert_eq!(partitions_p.len(), 4 * leafs_p.len() - 4);
        assert_eq!(partitions_q.len(), 4 * leafs_p.len() - 4);
    }

    // Find the number of parts in one but not the other.
    let cnt = partitions_p.symmetric_difference(&partitions_q).count();
    assert!(cnt % 2 == 0);
    cnt / 2
}

#[cfg(test)]
mod tests {
    use crate::io::newick;

    use super::*;

    #[test]
    fn upgma_wiki() {
        // Sample from Wikipedia
        let d = DistanceMatrix::new(
            vec!["a", "b", "c", "d", "e"]
                .iter()
                .map(|&s| s.into())
                .collect(),
            vec![
                vec![0., 17., 21., 31., 23.],
                vec![17., 0., 30., 34., 21.],
                vec![32., 40., 0., 28., 39.],
                vec![31., 34., 28., 0., 43.],
                vec![23., 21., 39., 43., 0.],
            ],
        )
        .unwrap();
        let s = "(((a:8.5,b:8.5):2.5,e:11):5.5,(c:14,d:14):2.5);";
        assert_eq!(newick::to_string(&upgma(d)).unwrap(), s);
    }

    #[test]
    fn neighbor_joining_wiki() {
        // Sample from Wikipedia
        let d = DistanceMatrix::new(
            vec!["a", "b", "c", "d", "e"]
                .iter()
                .map(|&s| s.into())
                .collect(),
            vec![
                vec![0., 5., 9., 9., 8.],
                vec![5., 0., 10., 10., 9.],
                vec![9., 10., 0., 8., 7.],
                vec![9., 10., 8., 0., 3.],
                vec![8., 9., 7., 3., 0.],
            ],
        )
        .unwrap();
        let s = "((((a:2,b:3):3,c:4):2,d:2):0.5,e:0.5);";
        assert_eq!(newick::to_string(&neighbor_joining(d)).unwrap(), s);
    }

    #[test]
    fn robinson_foulds_distance_small() {
        let p1 = newick::from_string("a;").unwrap();
        let p2 = newick::from_string("(a);").unwrap();
        let p3 = newick::from_string("((a));").unwrap();
        assert_eq!(robinson_foulds_distance(&p1, &p2), 0);
        assert_eq!(robinson_foulds_distance(&p1, &p3), 0);
        assert_eq!(robinson_foulds_distance(&p2, &p3), 0);

        let p1 = newick::from_string("(a,b,c);").unwrap();
        let p2 = newick::from_string("((c,a),b);").unwrap();
        let p3 = newick::from_string("(c,(b,a));").unwrap();
        assert_eq!(robinson_foulds_distance(&p1, &p2), 0);
        assert_eq!(robinson_foulds_distance(&p1, &p3), 0);
        assert_eq!(robinson_foulds_distance(&p2, &p3), 0);

        let p1 = newick::from_string("((a,b),(c,d));").unwrap();
        let p2 = newick::from_string("(((a,b),c),d);").unwrap();
        let p3 = newick::from_string("((a,c),(b,d));").unwrap();
        assert_eq!(robinson_foulds_distance(&p1, &p2), 0);
        assert_eq!(robinson_foulds_distance(&p1, &p3), 2);
        assert_eq!(robinson_foulds_distance(&p2, &p3), 2);
    }

    #[test]
    #[should_panic]
    fn robinson_foulds_distance_different_labels() {
        let p1 = newick::from_string("((a,b),(c,d));").unwrap();
        let p2 = newick::from_string("(((a,B),c),d);").unwrap();
        robinson_foulds_distance(&p1, &p2);
    }
}
