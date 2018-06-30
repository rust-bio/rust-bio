// Copyright 2017 Brett Bowman
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of long
//! high-error reads (e.g. PacBio or Oxford Nanopore).  Represent the
//! multi-sequece alignment of many reads as a graph-structure, and use
//! Smith-Waterman alignments to the consensus-path to efficiently add
//! new sequences to the graph.  Because of this graph structure, both
//! traceback (i.e. consensus) is approximately O(N), and memory usage
//! is also O(N), where N is the number of nodes in the graph.
//!
//! For the original concept and theory, see:
//! * Lee, Christopher, Catherine Grasso, and Mark F. Sharlow. "Multiple sequence alignment using
//! partial order graphs." Bioinformatics 18.3 (2002): 452-464.
//! * Lee, Christopher. "Generating consensus sequences from partial order multiple sequence
//! alignment graphs." Bioinformatics 19.8 (2003): 999-1008.
//!
//! For a reference implementation that inspired this code, see poapy:
//! https://github.com/ljdursi/poapy

use std::str;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::cmp::max;
use std::collections::HashMap;

use alignment::{Alignment, AlignmentOperation, AlignmentMode};
use utils::TextSlice;

use petgraph::{Graph, Directed, Incoming};
use petgraph::graph::NodeIndex;
use petgraph::visit::{EdgeRef, IntoNeighborsDirected, Topo, Walker};
use petgraph::dot::Dot;

/// A Partial-Order Alignment Graph
pub struct POAGraph {
    // the POAGraph struct stores a DAG and labels for the sequences which
    // compose it
    graph: Graph<u8, i32, Directed>,
    cns: Vec<u8>,
    cns_path: Vec<NodeIndex>,
    node_idx: Vec<NodeIndex>,
    needs_sort: bool,
}

pub struct Aligner {
    traceback: Vec<Vec<i32>>,
    scoring: fn(u8, u8) -> i32,
}

impl Aligner {
    pub fn new() -> Self {
       let score_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
       Aligner {
            traceback: vec![vec![]],
            scoring: score_fn,
        }
    }

    pub fn custom(&mut self, g: &Graph<u8, i32, Directed, usize>, query: TextSlice) -> Alignment {
        // dimensions of the score/traceback matrices 
        let (m, n) = (g.node_count(), query.len());

        self.traceback = vec![vec![0; n + 1]; m + 1];
        let ops: Vec<AlignmentOperation> = vec![];

        for i in 0..(m + 1) {
            self.traceback[i][0] = -1 * i as i32;
        }
        for j in 0..(n + 1) {
            self.traceback[0][j] = -1 * j as i32;
        }

        // construct the score matrix (naive) 
        let mut topo = Topo::new(&g);
        while let Some(node) = topo.next(&g) {
            // reference base and index
            let r = g.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;

            // predecesors of this node
            let prevs: Vec<NodeIndex<usize>> = g.neighbors_directed(node, Incoming).collect();
            // query base and index
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and insertion scores for first reference base
                let (mat, del) = if prevs.len() == 0 {
                    (self.traceback[0][j - 1] + (self.scoring)(r, *q),
                     self.traceback[0][j] - 1i32)
                } else {
                    let mut mat_max = -999;
                    let mut del_max = -999;
                    for prev_n in 0..prevs.len() {
                        let i_p: usize = prevs[prev_n].index() + 1; // index of previous node

                        mat_max = max(mat_max, self.traceback[i_p][j - 1] + (self.scoring)(r, *q));
                        del_max = max(del_max, self.traceback[i_p][j] - 1i32);
                    }
                    (mat_max, del_max)
                };
                let score = max(mat, max(del, self.traceback[i][j - 1] - 1i32));
                self.traceback[i][j] = score;
            }
        }

        Alignment {
            score: self.traceback[m][n],
            xstart: 0, ystart: 0, xend: m, yend: n, ylen: m, xlen: n,
            operations: ops, mode: AlignmentMode::Custom,
        }
    }
}

impl POAGraph {
    /// Create new POAGraph instance initialized with a sequence
    ///
    /// # Arguments
    ///
    /// * `label` - sequence label for an initial sequence
    /// * `sequence` - TextSlice from which to initialize the POA
    ///
    pub fn new_from_sequence(label: &str, sequence: TextSlice) -> POAGraph {
        let mut poa = POAGraph::new();
        poa.add_unmatched_sequence(label, sequence);
        return poa;
    }

    /// Create an empty POAGraph 
    pub fn new() -> POAGraph {
        POAGraph {
            graph: Graph::new(),
            cns: Vec::new(),
            cns_path: Vec::new(),
            node_idx: Vec::new(),
            needs_sort: false,
        }
    }

    /// Add a new unaligned sequence to the underlying graph.
    /// Useful for both initializing the graph with its first sequence, as
    /// well as adding unaligned prefix or suffix sequence from partially
    /// aligned reads.
    ///
    /// # Arguments
    ///
    /// * `label` - the id of the sequence to be added
    /// * `sequence` - The sequence to be added
    ///
    pub fn add_unmatched_sequence(&mut self, _label: &str, sequence: TextSlice) {
        // this should return a POAGraph
        let mut graph = Graph::new();
        let mut prev: NodeIndex = graph.add_node(sequence[0]); 
        let mut node: NodeIndex;
        for i in 1..sequence.len() {
            node = graph.add_node(sequence[i]);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        self.graph = graph;
    }

    /// Calculate and return the current consensus sequence as an array reference
    pub fn consensus(&mut self) -> &[u8] {
        // If we've added no new nodes or edges since the last call, sort first
        if self.needs_sort {
            self.topological_sort();
        }

        self.cns = Vec::new();
        for node in self.consensus_path().to_vec() {
            self.cns.push(*self.graph.node_weight(node).unwrap() as u8);
        }

        &self.cns
    }

    /// Calculate and return the current consensus-path through the underlying
    /// graph structure and return it as an array reference of node indices
    pub fn consensus_path(&mut self) -> &[NodeIndex] {
        // If we've added no new nodes or edges since the last call, sort first
        if self.needs_sort {
            self.topological_sort();
        }

        // For each node find the best predecessor by edge-weight, breaking ties with path-weight
        let mut scores = HashMap::new();
        let mut next_in_path: HashMap<NodeIndex, Option<NodeIndex>> = HashMap::new();
        for node in self.node_idx.iter().rev() {
            let mut best_neighbor = None::<NodeIndex>;
            let mut best_weights = (0, 0); // (Edge-weight, Path-weight)
            for e_ref in self.graph.edges(*node) {
                let weight = *e_ref.weight();
                let target = e_ref.target();

                if (weight, *scores.entry(target).or_insert(0)) > best_weights {
                    best_neighbor = Some(target);
                    best_weights = (weight, *scores.entry(target).or_insert(0));
                }
            }

            scores.insert(*node, best_weights.0 + best_weights.1);
            next_in_path.insert(*node, best_neighbor);
        }

        // Find the start position by locating the highest scoring node
        let mut start_node = None::<NodeIndex>;
        let mut best_score = 0;
        for (node, score) in &scores {
            if score > &best_score {
                start_node = Some(*node);
                best_score = *score;
            }
        }

        // Traverse the graph from the start node, recording the path of nodes taken
        self.cns_path = Vec::new();
        let mut curr_node = start_node;
        loop {
            if curr_node.is_some() {
                let node = curr_node.unwrap();
                self.cns_path.push(node);
                curr_node = *next_in_path.get(&node).unwrap();
            } else {
                break;
            }
        }

        &self.cns_path
    }

    /// Align a sequence to the current graph and return the
    /// resulting alignment object
    ///
    /// # Arguments
    ///
    /// * `sequence` - The new sequence to align as a text slice
    ///
    pub fn align_sequence(&self, seq: &[u8]) -> Alignment {
        let mut aligner = Aligner::new();
        let mut g: Graph<u8, i32, Directed, usize> = Graph::with_capacity(25, 25);
        let a = g.add_node(b'G');
        let b = g.add_node(b'A');
        let c = g.add_node(b'T');
        let d = g.add_node(b'T');
        let e = g.add_node(b'A');
        let f = g.add_node(b'C');
        let h = g.add_node(b'A');

        g.add_edge(a, b, 0);
        g.add_edge(b, c, 0);
        g.add_edge(c, d, 0);
        g.add_edge(d, e, 0);
        g.add_edge(e, f, 0);
        g.add_edge(f, h, 0);

        aligner.custom(&g, seq)
   }

    /// Incorporate a new sequence into the graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment object of the new sequence to the graph 
    /// * `label` - The name of the new sequence being added to the graph
    /// * `seq` - The complete sequence of the read being incorporated
    ///
    pub fn incorporate_alignment(&mut self, aln: Alignment, label: &str, seq: TextSlice) {
        // this will replace self.graph
    }

    /// Sort the nodes in the underlying graph topologically,
    /// such that every node index preceeds all nodes it has connecting
    /// edges to, and succeeds all nodes that have edges connecting to
    /// it.  This guarantees stable, fast results from traversals
    /// and consensus operations
    pub fn topological_sort(&mut self) {
        // TODO: use petgraph's toposort
    }

    /// Write the current graph to a specified filepath in dot format for
    /// visualization, primarily for debugging / diagnostics
    ///
    /// # Arguments
    ///
    /// * `filename` - The filepath to write the dot file to, as a String
    ///
    pub fn write_dot(&self, filename: String) {
        let mut file = match File::create(&filename) {
            Err(why) => panic!("couldn't open file {}: {}", filename, why.description()),
            Ok(file) => file,
        };
        let g = self.graph.map(|ni, nw| *nw as char,
                               |ni, ew| ew);
        match file.write_all(Dot::new(&g).to_string().as_bytes()) {
            Err(why) => panic!("couldn't write to file {}: {}", filename, why.description()),
            _ => (),
        }
    }
}

#[cfg(test)]
mod tests {
    use alignment::poa::POAGraph;

//    #[test]
//    fn test_align_sequence() {
//        assert_eq!(true, true);
//    }

    #[test]
    fn test_init_graph() {
        // sanity check for String -> Graph
        let poa = POAGraph::new_from_sequence("seq1", b"123456789");
        assert!(poa.graph.is_directed());
        assert_eq!(poa.graph.node_count(), 9);
        assert_eq!(poa.graph.edge_count(), 8);
    }

    #[test]
    fn test_incorporate() {
        let _seq1 = b"AAAAAA";
        let _seq2 = b"ABBBBA";
        let _seq3 = b"ABCCBA";

    }

    #[test]
    fn test_alignment() {
        let _seq1 = b"PKMIVRPQKNETV";
        let _seq2 = b"THKMLVRNETIM";
        let mut poa = POAGraph::new_from_sequence("seq1", _seq1);
        let alignment = poa.align_sequence(b"GCATGCU");
//        poa.incorporate_alignment(alignment, "seq2", seq2);
        assert_eq!(alignment.score, 0);
    }

    #[test]
    fn test_incorporate_alignment() {
        let seq1 = b"CCGCTTTTCCGC";
        let seq2 = b"CCGCAAAACCGC";
        let mut poa = POAGraph::new_from_sequence("seq1", seq1);
        let alignment = poa.align_sequence(seq2);
        poa.incorporate_alignment(alignment, "seq2", seq2);
    }
}
