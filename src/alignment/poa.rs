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
use std::cmp::{max, Ordering};

use alignment::{AlignmentMode, AlignmentOperation};
use utils::TextSlice;

use petgraph::{Graph, Directed, Incoming};
use petgraph::graph::{Node, NodeIndex, EdgeIndex};
use petgraph::visit::Topo;
use petgraph::dot::Dot;

/// A Partial-Order Alignment Graph
pub struct POAGraph {
    // the POAGraph struct stores a DAG and labels for the sequences which
    // compose it
    graph: Graph<u8, i32, Directed, usize>,
    cns: Vec<u8>,
    cns_path: Vec<NodeIndex>,
    node_idx: Vec<NodeIndex>,
    head: NodeIndex<usize>,
    tail: NodeIndex<usize>,
}

#[derive(Debug, Clone)]
pub enum Op {
    Match(Option<usize>),
    Del(Option<usize>),
    Ins(Option<usize>),
}

pub struct Alignment {
    score: i32,
    operations: Vec<Op>,
}

pub struct Aligner {
    traceback: Vec<Vec<Cell>>,
    scoring: fn(u8, u8) -> i32,
}

#[derive(Debug, Clone)]
struct Cell {
    score: i32,
    op: Op,
}

impl Ord for Cell {
    fn cmp(&self, other: &Cell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Cell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Cell {
    fn eq(&self, other: &Cell) -> bool {
        self.score == other.score
    }
}

impl Eq for Cell {}

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

        // initialize matrix
        self.traceback = vec![vec![Cell { score: 0, op: Op::Del(None) } ; n + 1]; m + 1];
        let mut ops: Vec<Op> = vec![];

        for i in 0..(m + 1) {
            // TODO: these should be -1 * distance from head node
            self.traceback[i][0] = Cell { score: -1 * i as i32, op: Op::Del(Some(i)) };
        }
        for j in 0..(n + 1) {
            self.traceback[0][j] = Cell { score: -1 * j as i32, op: Op::Ins(Some(0)) };
        }

        self.traceback[0][0] = Cell { score: 0, op: Op::Match(Some(0)) };

        // construct the score matrix (naive) 
        let mut topo = Topo::new(&g);

        // we'll be storing the last visited node in topological order so that
        // we can index into the end of the alignment
        let mut last: NodeIndex<usize> = NodeIndex::new(0);
        while let Some(node) = topo.next(&g) {
            // reference base and index
            let r = g.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            last = node;
            // predecesors of this node
            let prevs: Vec<NodeIndex<usize>> = g.neighbors_directed(node, Incoming).collect();
            // query base and index
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for first reference base
                let (mat, del) = if prevs.len() == 0 {
                    (Cell { score: self.traceback[0][j - 1].score + (self.scoring)(r, *q),
                            op: Op::Match(Some(0)) },
                     Cell { score: self.traceback[0][j].score - 1i32, op: Op::Del(Some(0)) })
                } else {
                    let mut mat_max = Cell { score: -999, op: Op::Match(None) };
                    let mut del_max = Cell { score: -999, op: Op::Del(None) };
                    for prev_n in 0..prevs.len() {
                        let i_p: usize = prevs[prev_n].index() + 1; // index of previous node

                        mat_max = max(mat_max,
                            Cell { score: self.traceback[i_p][j - 1].score + (self.scoring)(r, *q),
                                    op: Op::Match(Some(i_p))});
                        del_max = max(del_max,
                            Cell { score: self.traceback[i_p][j].score - 1i32, op: Op::Del(Some(i_p)) });
                    }
                    (mat_max, del_max)
                };
                let score = max(mat, max(del, Cell { score: self.traceback[i][j - 1].score - 1i32, 
                                                     op: Op::Ins(Some(i)) }));
                self.traceback[i][j] = score;
            }
        }

        print!(".\t");
        for i in 0..n {
            print!("{:?}\t", query[i]);
        }
        for i in 0..m {
            print!("\n{:?}\t", g.raw_nodes()[i].weight);
            for j in 0..n {
                print!("{}\t", self.traceback[i+1][j+1].score);
            }
        }

        let mut i = last.index() + 1;
        let mut j = n;
        while (i > 0 && j > 0) {
            // push operation and edge corresponding to optimal route
            // nb. optimal route is ambiguous :o
//            let edges: Vec<EdgeIndex<usize>> = g.edges(node).collect();
//            let mut op;
            println!("\t{}, {}", i, j);
            ops.push(self.traceback[i][j].op.clone());
            match self.traceback[i][j].op {
                Op::Match(Some(p)) => { i = p; j = max(j - 1, 0); },
                Op::Del(Some(p)) => { i = p; },
                Op::Ins(Some(p)) => { i = p; j = max(j - 1, 0); },
                Op::Match(None) => { i = i - 1; j = max(j - 1, 0); },
                Op::Del(None) => { },
                Op::Ins(None) => { i = max(i - 1, 0) },
            }

            if (i == 0 && j == 0) {
//                ops.push(self.traceback[0][0].op.clone());
                break;
            }
        }

        ops.reverse();
        println!("{:?}", ops);

        Alignment {
            score: self.traceback[last.index() + 1][n].score,
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
            cns: Vec::new(),
            cns_path: Vec::new(),
            node_idx: Vec::new(),
            head: NodeIndex::new(0),
            tail: NodeIndex::new(0),
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
    pub fn seq_to_graph(sequence: TextSlice) -> Graph<u8, i32, Directed, usize> {
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
        aligner.custom(&self.graph, seq)
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
        let mut prev: Node<i16, usize>;
        let attach = |n: usize, i: usize| -> (Node<i16, usize>, Node<i16, usize>) {
            let new = self.graph.add_node(seq[i]);
            let node = self.graph.raw_nodes()[n];
            self.graph.add_edge(node, new, 1);
            (node, new)
        };
        for op in aln {
            match op {
                Op::Match(Some(n)) => { prev = attach(n, n).0; },
                Op::Ins(Some(n)) => { prev = attach(n, n).1; },
                Op::Del(Some(n)) => { prev = attach(n, n).1; },
            }
        }
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
        let g = self.graph.map(|_, nw| *nw as char,
                               |_, ew| ew);
        match file.write_all(Dot::new(&g).to_string().as_bytes()) {
            Err(why) => panic!("couldn't write to file {}: {}", filename, why.description()),
            _ => (),
        }
    }
}

#[cfg(test)]
mod tests {
    use alignment::poa::POAGraph;
    use petgraph::graph::NodeIndex;

//    #[test]
//    fn test_align_sequence() {
//        assert_eq!(true, true);
//    }

    #[test]
    fn test_init_graph() {
        // sanity check for String -> Graph
        let poa = POAGraph::new("seq1", b"123456789");
        assert!(poa.graph.is_directed());
        assert_eq!(poa.graph.node_count(), 9);
        assert_eq!(poa.graph.edge_count(), 8);
    }

//    #[test]
//    fn test_incorporate() {
//        let _seq1 = b"AAAAAA";
//        let _seq2 = b"ABBBBA";
//        let _seq3 = b"ABCCBA";
//    }

    #[test]
    fn test_alignment() {
        // examples from the POA paper
        //let _seq1 = b"PKMIVRPQKNETV";
        //let _seq2 = b"THKMLVRNETIM";
        let poa = POAGraph::new("seq1", b"GATTACA");
        let alignment = poa.align_sequence(b"GCATGCU");
        assert_eq!(alignment.score, 0);
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
//        let seq1 = b"CCGCTTTTCCGC";
//        let seq2 = b"CCGCAAAACCGC";
        let seq1 = b"TTTTT";
        let seq2 = b"TTCTT";
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
//        poa.incorporate_alignment(alignment, "seq2", seq2);
    }

    #[test]
    fn test_insertion_on_branch() {
        let seq1 = b"TTTT";
        let seq2 = b"TTCCCATT";
        let mut poa = POAGraph::new("seq1", seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'C');
        let node2 = poa.graph.add_node(b'C');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.align_sequence(seq2);
        assert_eq!(alignment.score, 4);
        poa.write_dot("/tmp/out.dot".to_string());
//        poa.incorporate_alignment(alignment, "seq2", seq2);
    }

}
