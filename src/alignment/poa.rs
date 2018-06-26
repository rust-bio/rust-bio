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
//! https://github.com/ljdursi/poapy.get

use std::str;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::collections::HashMap;
use std::collections::HashSet;

use alignment::Alignment;
use alignment::AlignmentOperation;
use alignment::pairwise::banded;
use utils::TextSlice;

use petgraph;
use petgraph::Graph;
use petgraph::graph::{NodeIndex, EdgeIndex};
use petgraph::visit::EdgeRef;
use petgraph::dot::Dot;

/// A Partial-Order Alignment Graph
#[derive(Default)]
pub struct POAGraph {
    graph: Graph<u8, u16, petgraph::Directed>,
    cns: Vec<u8>,
    cns_path: Vec<NodeIndex>,
    node_idx: Vec<NodeIndex>,
    needs_sort: bool,
}

impl POAGraph {
    /// Create new POAGraph instance, optionally initialized with a sequence
    ///
    /// # Arguments
    ///
    /// * `label` - optional sequence label for an initial sequence
    /// * `sequence` - opetional TextSlice from which to initialize the POA
    ///
    pub fn new_from_sequence(label: Option<&str>, sequence: Option<TextSlice>) -> POAGraph {
        let mut poa = POAGraph {
            graph: Graph::new(),
            cns: Vec::new(),
            cns_path: Vec::new(),
            node_idx: Vec::new(),
            needs_sort: false,
        };

        // Only initialize with a sequence if both a label and sequence were provided
        match (label, sequence) {
            (Some(lab), Some(seq)) => poa.add_unmatched_sequence(lab, seq),
            (_, _) => (None, None),
        };

        return poa;
    }

    /// Create new POAGraph instance, optionally initialized with a sequence
    pub fn new() -> POAGraph {
        POAGraph::new_from_sequence(None::<&str>, None::<TextSlice>)
    }

    /// Return the number of nodes in the underlying graph
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Return the number of edges (note: not edge-weights) in the underlying graph
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Add a new sequence node to the underlying graph
    ///
    /// # Arguments
    ///
    /// * `weight` - The character or nucleotide for this position
    ///
    pub fn add_node(&mut self, weight: u8) -> NodeIndex {
        let new_node = self.graph.add_node(weight);
        self.needs_sort = true;

        new_node
    }

    /// Add a new edge connecting two sequence nodes in the graph
    ///
    /// # Arguments
    ///
    /// * `in_node` - The 'source' of the new edge
    /// * `out_node` - The 'sink' or 'target' of the new edge
    /// * `weight` - The weight associated with the new edge
    ///
    pub fn add_edge(&mut self, in_node: NodeIndex, out_node: NodeIndex, weight: u16) -> EdgeIndex {
        let new_edge = self.graph.add_edge(in_node, out_node, weight);
        self.needs_sort = true;

        new_edge
    }

    /// Add a new sequence new unaligned sequence to the underlying graph.
    /// Useful for both initializing the graph with it's first sequence, as
    /// well as adding unaligned prefix or suffix sequence from partially
    /// aligned reads.
    ///
    /// # Arguments
    ///
    /// * `label` - the id of the sequence to be added
    /// * `sequence` - The sequence to be added
    ///
    pub fn add_unmatched_sequence(
        &mut self,
        _label: &str,
        sequence: TextSlice,
    ) -> (Option<NodeIndex>, Option<NodeIndex>) {
        // Add each character to the graph, and it's associated edges
        let mut first = None::<NodeIndex>;
        let mut last = None::<NodeIndex>;
        for b in sequence {
            let node = self.add_node(*b);
            if first.is_none() {
                first = Some(node);
            }
            if last.is_some() {
                self.add_edge(last.unwrap(), node, 1);
            }
            last = Some(node);
        }

        // Return the book-end positions of the new sequence
        return (first, last);
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

    /// Align a sequence to the current consensus and return the
    /// resulting alignment object
    ///
    /// # Arguments
    ///
    /// * `sequence` - The new sequence to align as a text slice
    ///
    pub fn align_sequence(&self, seq: &[u8]) -> Alignment {
        let score_fn = |a: u8, b: u8| if a == b { 4i32 } else { -2i32 };
        let mut aligner = banded::Aligner::new(-4, -2, &score_fn, 8, 10);

        aligner.global(seq, &self.cns)
    }

    /// Incorporate a new sequence into the graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment object of the new sequence to the current consensus
    /// * `label` - The name of the new sequence being added to the graph
    /// * `seq` - The complete sequence of the read being incorporated
    ///
    pub fn incorporate_alignment(&mut self, aln: Alignment, label: &str, seq: TextSlice) {
        // Ends of sequence may be unaligned, and incorporated separately
        let mut head_node = None::<NodeIndex>;
        let mut tail_node = None::<NodeIndex>;
        if aln.xstart > 0 {
            let (_start_node, end_node) = self.add_unmatched_sequence(&label, &seq[..aln.xstart]);
            head_node = end_node;
        }
        if aln.xend < seq.len() {
            let (start_node, _end_node) = self.add_unmatched_sequence(&label, &seq[aln.xend + 1..]);
            tail_node = start_node;
        }

        // Since each transition is an graph edge, we need to choose an initial source 'prev_node':
        //   1) Use the last node of the unmatched prefix as the initial "prev_node", start at 0
        //   2) Use the node at 0 as the initial "prev_node", start at 1
        let (astart, mut xpos, mut ypos, mut prev_node) = match head_node {
            Some(_n) => (0, aln.xstart, aln.ystart, head_node.unwrap()),
            None => (1, aln.xstart + 1, aln.ystart + 1, self.cns_path[aln.ystart]),
        };

        for apos in astart..aln.operations.len() {
            let op = aln.operations[apos];
            match op {
                AlignmentOperation::Match => {
                    let curr_node = self.cns_path[ypos];
                    match self.graph.find_edge(prev_node, curr_node) {
                        Some(e) => {
                            *self.graph.edge_weight_mut(e).unwrap() += 1;
                        }
                        None => {
                            self.add_edge(prev_node, curr_node, 1);
                        }
                    }
                    xpos += 1;
                    ypos += 1;
                    prev_node = self.cns_path[ypos - 1];
                }
                AlignmentOperation::Subst => {
                    let xchar = seq[xpos];
                    let out_edge = self.graph
                        .edges(prev_node)
                        .find(|e| *self.graph.node_weight(e.target()).unwrap() as u8 == xchar)
                        .map(|e| e.id());
                    match out_edge {
                        Some(e) => {
                            *self.graph.edge_weight_mut(e).unwrap() += 1;
                            prev_node = self.graph.edge_endpoints(e).unwrap().1;
                        }
                        None => {
                            let new_node = self.add_node(xchar);
                            self.add_edge(prev_node, new_node, 1);
                            prev_node = new_node;
                        }
                    }
                    xpos += 1;
                    ypos += 1;
                }
                AlignmentOperation::Del => {
                    ypos += 1;
                }
                AlignmentOperation::Ins => {
                    let xchar = seq[xpos];
                    let out_edge = self.graph
                        .edges(prev_node)
                        .filter(|e| e.target() != self.cns_path[ypos])
                        .find(|e| *self.graph.node_weight(e.target()).unwrap() as u8 == xchar)
                        .map(|e| e.id());
                    match out_edge {
                        Some(e) => {
                            *self.graph.edge_weight_mut(e).unwrap() += 1;
                            prev_node = self.graph.edge_endpoints(e).unwrap().1;
                        }
                        None => {
                            let new_node = self.add_node(xchar);
                            self.add_edge(prev_node, new_node, 1);
                            prev_node = new_node;
                        }
                    }
                    xpos += 1;
                }
                _ => ()  // Skip any clipped bases
            }

            // Finally, if the sequence had an unmatched suffix, connect it to the alignment's end
            if tail_node.is_some() {
                self.add_edge(prev_node, tail_node.unwrap(), 1);
            }
        }
    }

    /// Sort the nodes in the underlying graph topologically,
    /// such that every node index preceeds all nodes it has connecting
    /// edges to, and succeeds all nodes that have edges connecting to
    /// it.  This guarantees stable, fast results from traversals
    /// and consensus operations
    pub fn topological_sort(&mut self) {
        let mut sorted_indices: Vec<NodeIndex> = Vec::new();
        let mut completed: HashSet<NodeIndex> = HashSet::new();

        // Depth-first search for unsorted nodes from some defined start-point
        let dfs = |graph: &Graph<u8, u16, petgraph::Directed>,
                   start: NodeIndex,
                   completed: &mut HashSet<NodeIndex>,
                   sorted_indices: &mut Vec<NodeIndex>|
         -> () {
            let mut stack = vec![start];
            let mut started = HashSet::new();
            loop {
                let node = match stack.pop() {
                    Some(n) => n,
                    None => {
                        break;
                    }
                };

                if completed.contains(&node) {
                    continue;
                }

                if started.contains(&node) {
                    completed.insert(node);
                    sorted_indices.insert(0, node);
                    started.remove(&node);
                    continue;
                }

                let successors = graph
                    .edges(node)
                    .map(|e| e.target())
                    .filter(|n| !completed.contains(n));
                started.insert(node);
                stack.push(node);
                stack.extend(successors);
            }
        };

        // Loop over all uncompleted nodes, performing a depth-first search on each
        loop {
            if sorted_indices.len() == self.node_count() {
                break;
            }

            let mut found = None::<NodeIndex>;
            for node in self.graph.node_indices() {
                if !completed.contains(&node) {
                    found = Some(node);
                    break;
                }
            }
            assert!(found.is_some());
            dfs(
                &self.graph,
                found.unwrap(),
                &mut completed,
                &mut sorted_indices,
            );
        }

        // Finally, store the results of the sort
        self.node_idx = sorted_indices;
        self.needs_sort = false;
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
