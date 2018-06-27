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
use std::collections::{HashMap, HashSet};

use alignment::{Alignment, AlignmentOperation, AlignmentMode};
use utils::TextSlice;

use petgraph::{Graph, Directed};
use petgraph::graph::{NodeIndex, EdgeIndex};
use petgraph::visit::EdgeRef;
use petgraph::dot::Dot;
use petgraph::algo::toposort;

/// A Partial-Order Alignment Graph
//#[derive(Default)]
pub struct POAGraph {
    graph: Graph<u8, u16, Directed>,
    cns: Vec<u8>,
    cns_path: Vec<NodeIndex>,
    node_idx: Vec<NodeIndex>,
    needs_sort: bool,
}

pub struct Aligner {
    traceback: Vec<Vec<u8>>,
    scoring: fn(u8, u8) -> i32,
}

impl Aligner {
    pub fn new() -> Self {
       let score_fn = |a: u8, b: u8| if a == b { 4i32 } else { -2i32 };
       Aligner {
            traceback: vec![vec![]],
            scoring: score_fn,
        }
    }

    pub fn custom(&mut self, x: &POAGraph, y: TextSlice) -> Alignment {
        // dimensions of the score/traceback matrices 
        let (m, n) = (x.graph.node_count(), y.len());
        
        self.traceback = vec![vec![]];
        // these iterators define the area of the traceback matrix
        //for x in 0..seq.len() {
            // traverse graph nodes in topological order
        //    for y in self.graph.nodes() {
                // instead of the predecesors from the familiar totally ordered
                // case, we instead consider the predecesor nodes from the
                // sequence graph. That's all.
        //        for pred in self.graph.edges() {
        //        score = min(score_fn(seq[i], self.graph[y]),
        //                    score_fn(seq[i-1], self.graph[y]) + gap,
        //                    preds); 

        //    }
        //}
        Alignment {
            score: 0, xstart: 0, ystart: 0, xend: 0, yend: 0, ylen: 0, xlen: 0,
            operations: vec![], mode: AlignmentMode::Custom, }
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
        let mut prev = graph.add_node(sequence[0]); 
        for i in 1..sequence.len() {
            let node = graph.add_node(sequence[i]);
            graph.add_edge(node, prev, 1);
            let prev = node;
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

        aligner.custom(self, seq)
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
mod poa {
    use alignment::poa::POAGraph;

//    #[test]
//    fn test_align_sequence() {
//        assert_eq!(true, true);
//    }

    #[test]
    fn test_init_graph() {
        let poa = POAGraph::new_from_sequence("seq1", b"ASDFGHJKL");
    }

    #[test]
    fn test_alignment() {
        let seq1 = b"PKMIVRPQKNETV";
        let seq2 = b"THKMLVRNETIM";
        let mut poa = POAGraph::new_from_sequence("seq1", seq1);
        let alignment = poa.align_sequence(seq2);
        poa.incorporate_alignment(alignment, "seq2", seq2);
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
