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

use utils::TextSlice;

use petgraph::{Graph, Directed, Incoming};
use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;
use petgraph::dot::Dot;

pub const MIN_SCORE: i32 = -858_993_459; // see alignment/pairwise/mod.rs

/// A Partial-Order Alignment Graph
pub struct POAGraph {
    // the POAGraph struct stores a DAG and labels for the sequences which
    // compose it
    graph: Graph<u8, i32, Directed, usize>,
}

// Unlike a total order alignment we have different possible successor nodes
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
            scoring: score_fn,
        }
    }

    // naive Needleman-Wunsch
    pub fn global(&mut self,
                  g: &Graph<u8, i32, Directed, usize>,
                  query: TextSlice) -> Alignment {
        // dimensions of the traceback matrix
        let (m, n) = (g.node_count(), query.len());

        // initialize matrix
        let mut traceback: Vec<Vec<Cell>> =
            vec![vec![Cell { score: 0, op: Op::Match(None) }; n + 1]; m + 1];
        let mut ops: Vec<Op> = vec![];

        for i in 1..(m + 1) {
            // TODO: these should be -1 * distance from head node
            traceback[i][0] = Cell { score: -1 * i as i32, op: Op::Del(None) };
        }
        for j in 1..(n + 1) {
            traceback[0][j] = Cell { score: -1 * j as i32, op: Op::Ins(None) };
        }

        traceback[0][0] = Cell { score: 0, op: Op::Match(None) };

        // construct the score matrix (naive) 
        let mut topo = Topo::new(&g);

        // store the last visited node in topological order so that
        // we can index into the end of the alignment
        let mut last: NodeIndex<usize> = NodeIndex::new(0);
        while let Some(node) = topo.next(&g) {
            // reference base and index
            let r = g.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            last = node;
            // predecesors of this node
            let prevs: Vec<NodeIndex<usize>> = g.neighbors_directed(node, Incoming).collect();
            // query base and index (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for first reference base
                let (mat, del) = if prevs.len() == 0 {
                    (Cell { score: traceback[0][j - 1].score + (self.scoring)(r, *q),
                            op: Op::Match(None) },
                     Cell { score: traceback[0][j].score - 1i32, op: Op::Del(None) })
                } else {
                    let mut mat_max = Cell { score: MIN_SCORE, op: Op::Match(None) };
                    let mut del_max = Cell { score: MIN_SCORE, op: Op::Del(None) };
                    for prev_n in 0..prevs.len() {
                        let i_p: usize = prevs[prev_n].index() + 1; // index of previous node
                        mat_max = max(mat_max,
                            Cell { score: traceback[i_p][j - 1].score + (self.scoring)(r, *q),
                                    op: Op::Match(Some((i_p - 1, i - 1)))});
                        del_max = max(del_max,
                            Cell { score: traceback[i_p][j].score - 1i32,
                                   op: Op::Del(Some((i_p - 1, i)))});
                    }
                    (mat_max, del_max)
                };
                let score = max(mat, max(del, Cell { score: traceback[i][j - 1].score - 1i32, 
                                                     op: Op::Ins(Some(i - 1)) }));
                traceback[i][j] = score;
            }
        }
        
//        debug(&traceback, g, query);

        // Now backtrack through the matrix to construct an optimal path
        let mut i = last.index() + 1;
        let mut j = n;
        //println!("last node: {:?} {:?}", i, j);
        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            //println!("\t{}, {} => {}", i, j, traceback[i][j].score);
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
//        let mut i: usize = seq.len();
        let mut i: usize = 0;
        for op in aln.operations {
            match op {
                Op::Match(None) => { i = i + 1; },
                Op::Match(Some((0, _))) => {
                    println!("(M) linking to zeroth node");
                    prev = NodeIndex::new(0);
                    i = i + 1;
                }
                Op::Match(Some((n, p))) => { 
                    let node = NodeIndex::new(n);
                    println!("(M) linking {:?} and {:?}", prev, node);
                    self.graph.add_edge(prev, node, 1);
                    prev = NodeIndex::new(p);
                    i = i + 1;
                },
                Op::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    println!("(I) linking {:?} and {:?}", prev, node);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i = i + 1; },
                Op::Del(Some((n, _))) => { 
                    println!("(D) skipping link {:?}", n);
                },
                Op::Ins(None) => {},
                Op::Del(None) => {},
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
fn debug(traceback: &Vec<Vec<Cell>>,
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
        // example from Wikipedia's Needleman-Wunsch distance article
        //let _seq1 = b"PKMIVRPQKNETV";
        //let _seq2 = b"THKMLVRNETIM";
        let mut poa = POAGraph::new("seq1", b"GATTACA");
        let alignment = poa.align_sequence(b"GCATGCU");
        assert_eq!(alignment.score, 0);
        poa.incorporate_alignment(alignment, "seq2", b"GCATCGU");
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
        assert_eq!(alignment.score, 3);
//        poa.incorporate_alignment(alignment, "seq2", seq2);
    }

    #[test]
    fn test_incorporation() {
//        example from POA paper
//        let seq1 = b"CCGCTTTTCCGC";
//        let seq2 = b"CCGCAAAACCGC";
        let seq1 = b"TTCCGGTTTAA";
        let seq2 = b"TTGGTTTGGGAA";
        let mut poa = POAGraph::new("seq1", seq1);
        let alignment = poa.align_sequence(seq2);
//        poa.write_dot("/tmp/inc1.dot".to_string());
        poa.incorporate_alignment(alignment, "seq2", seq2);
//        poa.write_dot("/tmp/inc2.dot".to_string());
//        assert!(false);
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
//        poa.write_dot("/tmp/qut1.dot".to_string());
        let alignment = poa.align_sequence(seq2);
        poa.incorporate_alignment(alignment, "seq2", seq2);
//        poa.write_dot("/tmp/qut2.dot".to_string());
//        assert_eq!(alignment.score, -4);
//        assert!(false);
    }

}
