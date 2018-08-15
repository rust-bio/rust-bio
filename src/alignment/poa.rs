// Copyright 2017 Brett Bowman
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of long homologous sequences.
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
use std::cmp::{max, Ordering};

use utils::TextSlice;

use petgraph::{Graph, Directed, Incoming};
use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;
use petgraph::dot::Dot;

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs

/// A Partial-Order Alignment Graph
pub struct POAGraph {
    // a POAGraph struct stores a reference DAG and labels the edges with the
    // count or names of the input sequences (TODO)
    graph: Graph<u8, i32, Directed, usize>,
}

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum Op {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}

pub struct Alignment {
    score: i32,
    operations: Vec<Op>,
}

pub struct Aligner {
    scoring: fn(u8, u8) -> i32,
}

#[derive(Debug, Clone)]
struct TracebackCell {
    score: i32,
    op: Op,
}

impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.score == other.score
    }
}

impl Eq for TracebackCell {}

impl Aligner {
    pub fn new() -> Self {
       let score_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
       Aligner {
            scoring: score_fn,
        }
    }

    /// Naive Needleman-Wunsch
    /// Populates the traceback matrix in $O(n^2)$ time using
    /// petgraph's constant time topological traversal
    pub fn global(&mut self,
                  g: &Graph<u8, i32, Directed, usize>,
                  query: TextSlice) -> Alignment {
        // dimensions of the traceback matrix
        let (m, n) = (g.node_count(), query.len());

        // initialize matrix
        let mut traceback: Vec<Vec<TracebackCell>> =
            vec![vec![TracebackCell { score: 0, op: Op::Match(None) }; n + 1]; m + 1];
        let mut ops: Vec<Op> = vec![];

        for i in 1..(m + 1) {
            // TODO: these should be -1 * distance from head node
            traceback[i][0] = TracebackCell { score: -1 * i as i32, op: Op::Del(None) };
        }
        for j in 1..(n + 1) {
            traceback[0][j] = TracebackCell { score: -1 * j as i32, op: Op::Ins(None) };
        }

        traceback[0][0] = TracebackCell { score: 0, op: Op::Match(None) };

        // construct the score matrix (naive) 
        let mut topo = Topo::new(&g);

        // store the last visited node in topological order so that
        // we can index into the end of the alignment when we backtrack
        let mut last: NodeIndex<usize> = NodeIndex::new(0);
        while let Some(node) = topo.next(&g) {
            // reference base and index
            let r = g.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> = g.neighbors_directed(node, Incoming).collect();
            // query base and its index in the DAG (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let (mat, del) = if prevs.len() == 0 {
                    (TracebackCell { score: traceback[0][j - 1].score + (self.scoring)(r, *q),
                            op: Op::Match(None) },
                     TracebackCell { score: traceback[0][j].score - 1i32, op: Op::Del(None) })
                } else {
                    let mut mat_max = TracebackCell { score: MIN_SCORE, op: Op::Match(None) };
                    let mut del_max = TracebackCell { score: MIN_SCORE, op: Op::Del(None) };
                    for prev_n in 0..prevs.len() {
                        let i_p: usize = prevs[prev_n].index() + 1; // index of previous node
                        mat_max = max(mat_max,
                            TracebackCell { score: traceback[i_p][j - 1].score + (self.scoring)(r, *q),
                                    op: Op::Match(Some((i_p - 1, i - 1)))});
                        del_max = max(del_max,
                            TracebackCell { score: traceback[i_p][j].score - 1i32,
                                   op: Op::Del(Some((i_p - 1, i)))});
                    }
                    (mat_max, del_max)
                };
                let score = max(mat, max(del, TracebackCell { score: traceback[i][j - 1].score - 1i32, 
                                                     op: Op::Ins(Some(i - 1)) }));
                traceback[i][j] = score;
            }
        }
        
        //dump_traceback(&traceback, g, query);

        // Now backtrack through the matrix to construct an optimal path
        let mut i = last.index() + 1;
        let mut j = n;

        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(traceback[i][j].op.clone());
            match traceback[i][j].op {
                Op::Match(Some((p, _))) => { i = p + 1; j = j - 1; },
                Op::Del(Some((p, _))) => { i = p + 1; },
                Op::Ins(Some(p)) => { i = p + 1; j = j - 1; },
                Op::Match(None) => { break; },
                Op::Del(None) => { j = j - 1; },
                Op::Ins(None) => { i = i - 1; },
            }
        }

        ops.reverse();

        Alignment {
            score: traceback[last.index() + 1][n].score,
            operations: ops,
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
    pub fn new(_label: &str, sequence: TextSlice) -> POAGraph {
        let graph = POAGraph::seq_to_graph(sequence);
        POAGraph {
            graph: graph,
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
    fn seq_to_graph(sequence: TextSlice) -> Graph<u8, i32, Directed, usize> {
        // this should return a POAGraph
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(sequence.len(), sequence.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(sequence[0]);
        let mut node: NodeIndex<usize>;
        for i in 1..sequence.len() {
            node = graph.add_node(sequence[i]);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        graph
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
        aligner.global(&self.graph, seq)
   }

    /// Incorporate a new sequence into the graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment object of the new sequence to the graph 
    /// * `label` - The name of the new sequence being added to the graph
    /// * `seq` - The complete sequence of the read being incorporated
    ///
    pub fn incorporate_alignment(&mut self, aln: Alignment, _label: &str, seq: TextSlice) {
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        for op in aln.operations {
            match op {
                Op::Match(None) => { i = i + 1; },
                Op::Match(Some((_, p))) => { 
                    let node = NodeIndex::new(p);
                    if seq[i] != self.graph.raw_nodes()[p].weight {
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                // where the previous node was newly added
                                self.graph.add_edge(prev, node, 1);
                            }
                        }
                        prev = NodeIndex::new(p);
                    }
                    i = i + 1;
                },
                Op::Ins(None) => { i = i + 1; },
                Op::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i = i + 1; },
                Op::Del(_) => { }, // we should only have to skip over deleted nodes
            }
        }
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
        let g = self.graph.map(|_, nw| *nw as char,
                               |_, ew| ew);
        match file.write_all(Dot::new(&g).to_string().as_bytes()) {
            Err(why) => panic!("couldn't write to file {}: {}", filename, why.description()),
            _ => (),
        }
    }
}

// print out a traceback matrix
#[allow(dead_code)]
fn dump_traceback(traceback: &Vec<Vec<TracebackCell>>,
         g: &Graph<u8, i32, Directed, usize>,
         query: TextSlice) {
    let (m, n) = (g.node_count(), query.len());
    print!(".\t");
    for i in 0..n {
        print!("{:?}\t", query[i]);
    }
    for i in 0..m {
        print!("\n{:?}\t", g.raw_nodes()[i].weight);
        for j in 0..n {
            print!("{}.\t", traceback[i+1][j+1].score);
        }
    }
    print!("\n");
}

#[cfg(test)]
mod tests {
    use alignment::poa::POAGraph;
    use petgraph::graph::NodeIndex;

    #[test]
    fn test_init_graph() {
        // sanity check for String -> Graph
        let poa = POAGraph::new("seq1", b"123456789");
        assert!(poa.graph.is_directed());
        assert_eq!(poa.graph.node_count(), 9);
        assert_eq!(poa.graph.edge_count(), 8);
    }

    #[test]
    fn test_alignment() {
        // examples from the POA paper
        //let _seq1 = b"PKMIVRPQKNETV";
        //let _seq2 = b"THKMLVRNETIM";
        let poa = POAGraph::new("seq1", b"GATTACA");
        let alignment = poa.align_sequence(b"GCATGCU");
        assert_eq!(alignment.score, 0);

        let alignment = poa.align_sequence(b"GCATGCUx");
        assert_eq!(alignment.score, -1);
 
        let alignment = poa.align_sequence(b"xCATGCU");
        assert_eq!(alignment.score, -2);
 
    }

    #[test]
    fn test_branched_alignment() {
        let seq1 = b"TTTTT";
        let seq2 = b"TTATT";
        let mut poa = POAGraph::new("seq1", seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.align_sequence(seq2);
        assert_eq!(alignment.score, 3);
    }

    #[test]
    fn test_alt_branched_alignment() {
        let seq1 = b"TTCCTTAA";
        let seq2 = b"TTTTGGAA";
        let mut poa = POAGraph::new("seq1", seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.align_sequence(seq2);
        poa.incorporate_alignment(alignment, "seq2", seq2);
        assert_eq!(poa.graph.edge_count(), 14);
        assert!(poa.graph.contains_edge(NodeIndex::new(5), NodeIndex::new(10)));
        assert!(poa.graph.contains_edge(NodeIndex::new(11), NodeIndex::new(6)));
    }

    #[test]
    fn test_insertion_on_branch() {
        let seq1 = b"TTCCGGTTTAA";
        let seq2 = b"TTGGTATGGGAA";
        let seq3 = b"TTGGTTTGCGAA";
        let mut poa = POAGraph::new("seq1", seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'C');
        let node2 = poa.graph.add_node(b'C');
        let node3 = poa.graph.add_node(b'C');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, node3, 1);
        poa.graph.add_edge(node3, tail, 1);
        let alignment = poa.align_sequence(seq2);
        assert_eq!(alignment.score, 2);
        poa.incorporate_alignment(alignment, "seq2", seq2);
        let alignment2 = poa.align_sequence(seq3);
        assert_eq!(alignment2.score, 10);
    }
}
