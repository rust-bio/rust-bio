// Copyright 2017-2018 Brett Bowman, Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of multiple homologous sequences.
//!
//! For the original concept and theory, see:
//! * Lee, Christopher, Catherine Grasso, and Mark F. Sharlow. "Multiple sequence alignment using
//! partial order graphs." Bioinformatics 18.3 (2002): 452-464.
//! * Lee, Christopher. "Generating consensus sequences from partial order multiple sequence
//! alignment graphs." Bioinformatics 19.8 (2003): 999-1008.
//!
//! For a modern reference implementation, see poapy:
//! https://github.com/ljdursi/poapy
//!
//! # Example
//!
//! ```
//! use bio::alignment::poa::*;
//! let x = b"AAAAAAA";
//! let y = b"AABBBAA";
//! let z = b"AABCBAA";
//!
//! let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//! let mut aligner = Aligner::new(-1, &score);
//! aligner.add_sequence(x);
//! // z differs from x in 3 locations
//! assert_eq!(aligner.global(z).score, 1);
//! aligner.add_sequence(y);
//! // z differs from x and y's partial order alignment by 1 base
//! assert_eq!(aligner.global(z).score, 5);
//!
//! // an alignment can be recycled to add a sequence to the graph
//! let alignment = aligner.global(z);
//! aligner.add_alignment(alignment, z);
//! assert_eq!(aligner.global(z).score, 7); // a complete match
//! ```

use std::cmp::{max, Ordering};
use std::ops::Neg;

use utils::TextSlice;

use alignment::pairwise::{MatchFunc, Scoring};
use alignment::AlignmentMode;

use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}

pub struct Alignment {
    pub score: i32,
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
    mode: AlignmentMode,
}

#[derive(Debug, Clone)]
struct TracebackCell {
    score: i32,
    op: AlignmentOperation,
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

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

struct Traceback {
    rows: usize,
    cols: usize,
    matrix: Vec<Vec<TracebackCell>>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    ///
    fn with_capacity(m: usize, n: usize) -> Self {
        let mut matrix = vec![
            vec![
                TracebackCell {
                    score: 0,
                    op: AlignmentOperation::Match(None)
                };
                n + 1
            ];
            m + 1
        ];
        for (i, row) in matrix.iter_mut().enumerate().take(m + 1).skip(1) {
            // TODO: these should be -1 * distance from head node
            row[0] = TracebackCell {
                score: -(i as i32), // gap_open penalty
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..=n {
            matrix[0][j] = TracebackCell {
                score: -(j as i32),
                op: AlignmentOperation::Ins(None),
            };
        }

        Traceback {
            rows: m,
            cols: n,
            matrix: matrix,
        }
    }

    fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            matrix: Vec::new(),
        }
    }
}

/// A multiple sequence aligner builder
///
/// Uses builder pattern for constructing partial order alignments with method chainging
pub struct POA<F: MatchFunc> {
    aligner: Aligner<F>,
    refs: Vec<String>,
}

impl<F: MatchFunc> POA<F> {
    //    pub fn new() -> POA {
    //        POA { aligner: Aligner }
    //

    //    pub fn add_sequence() -> {
    //        POA { aligner: Aligner::new(self.scoring). }
    //    }
}

/// A global aligner on partially ordered graphs
///
/// Internally stores a directed acyclic graph datastructure that informs the topology of the
/// traceback matrix. A mutable Aligner may have additional sequences added to the internal
/// partially ordered graph.
///
/// The partial order is expressed with a directed acyclic graph and informs the traversal of the
/// traceback matrix.
///
pub struct Aligner<F: MatchFunc> {
    scoring: Scoring<F>,
    pub graph: Graph<u8, i32, Directed, usize>,
}

impl<F: MatchFunc> Aligner<F> {
    pub fn new(gap_open: i32, match_fn: F) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, 0, match_fn),
            //            traceback: Traceback::new(),
            graph: Graph::with_capacity(0, 0),
        }
    }
    /// Create a new aligner instance from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    /// * `gap_open` - the negative score assigned when branching from the reference graph
    /// * `match_fn` - the pairwise score for substitutions (see bio::scores)
    ///
    pub fn from_string(seq: TextSlice, gap_open: i32, match_fn: F) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }

        Aligner {
            scoring: Scoring::new(gap_open, 0, match_fn),
            graph: graph,
        }
    }

    /// Create a new aligner instance from the directed acyclic graph of another.
    ///
    /// # Arguments
    ///
    /// * `poa` - the partially ordered reference alignment
    /// * `gap_open` - the negative score assigned when branch from the reference graph
    /// * `match_fn` - the pairwise score for substitutions (see bio::scores)
    ///
    pub fn from_graph(poa: Graph<u8, i32, Directed, usize>, gap_open: i32, match_fn: F) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, 0, match_fn),
            graph: poa,
        }
    }

    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    ///
    pub fn global(&self, query: TextSlice) -> Alignment {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n);

        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        traceback.matrix[0][0] = TracebackCell {
            score: 0,
            op: AlignmentOperation::Match(None),
        };

        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);

        // store the last visited node in topological order so that
        // we can index into the end of the alignment when we backtrack
        let mut last: NodeIndex<usize> = NodeIndex::new(0);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            // query base and its index in the DAG (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    TracebackCell {
                        score: traceback.matrix[0][j - 1].score
                            + self.scoring.match_fn.score(r, *q),
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        score: MIN_SCORE,
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_n in 0..prevs.len() {
                        let i_p: usize = prevs[prev_n].index() + 1; // index of previous node
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: traceback.matrix[i_p][j - 1].score
                                        + self.scoring.match_fn.score(r, *q),
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.matrix[i_p][j].score + self.scoring.gap_open,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };

                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.matrix[i][j - 1].score + self.scoring.gap_open,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.matrix[i][j] = score;
            }
        }

        //dump_traceback(traceback, g, query);

        // Now backtrack through the matrix to construct an optimal path
        let mut i = last.index() + 1;
        let mut j = n;

        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(traceback.matrix[i][j].op.clone());
            match traceback.matrix[i][j].op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j = j - 1;
                }
                AlignmentOperation::Match(None) => {
                    break;
                }
                AlignmentOperation::Del(None) => {
                    j -= 1;
                }
                AlignmentOperation::Ins(None) => {
                    i -= 1;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: traceback.matrix[last.index() + 1][n].score,
            operations: ops,
            mode: AlignmentMode::Custom,
        }
    }

    /// Experimental: return sequence of traversed edges
    ///
    /// Only supported alignments for sequences that have already been added,
    /// so all operations must be Match.
    pub fn edges(&self, aln: Alignment) -> Vec<usize> {
        let mut path: Vec<usize> = vec![];
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        for op in aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(p);
                    let edge = self.graph.find_edge(prev, node).unwrap();
                    path.push(edge.index());
                    prev = NodeIndex::new(p);
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {}
                AlignmentOperation::Ins(Some(_)) => {}
                AlignmentOperation::Del(_) => {}
            }
        }
        path
    }

    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    ///
    pub fn add_alignment(&mut self, aln: Alignment, seq: TextSlice) {
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        for op in aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(p);
                    if (seq[i] != self.graph.raw_nodes()[p].weight) && (seq[i] != b'X') {
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
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }

    /// Align and incorporate a sequence to the partially ordered graph from a TextSlice
    ///
    /// # Arguments
    /// * `seq` - TextSlice to incorporate into the partial order aligner's graph
    ///
    pub fn add_sequence(&mut self, seq: TextSlice) {
        // if the Aligner has not been initialized with any sequence yet
        if self.graph.node_count() == 0 {
            self.graph = Graph::with_capacity(seq.len(), seq.len() - 1);
            let mut prev: NodeIndex<usize> = self.graph.add_node(seq[0]);
            let mut node: NodeIndex<usize>;
            for base in seq.iter().skip(1) {
                node = self.graph.add_node(*base);
                self.graph.add_edge(prev, node, 1);
                prev = node;
            }
        }
        let alignment = self.global(seq);
        self.add_alignment(alignment, seq);
    }
}

// print out a traceback matrix
#[allow(dead_code)]
fn dump_traceback(
    traceback: &[Vec<TracebackCell>],
    g: &Graph<u8, i32, Directed, usize>,
    query: TextSlice,
) {
    let (m, n) = (g.node_count(), query.len());
    print!(".\t");
    for base in query.iter().take(n) {
        print!("{:?}\t", *base);
    }
    for i in 0..m {
        print!("\n{:?}\t", g.raw_nodes()[i].weight);
        for j in 0..n {
            print!("{}.\t", traceback[i + 1][j + 1].score);
        }
    }
    print!("\n");
}

#[cfg(test)]
mod tests {
    use alignment::poa::Aligner;
    use petgraph::dot::Dot;
    use petgraph::graph::NodeIndex;
    use petgraph::{Directed, Graph};
    use std::error::Error;
    use std::fs::File;
    use std::io::Write;

    /// Write the current graph to a specified filepath in dot format for
    /// visualization, primarily for debugging / diagnostics
    ///
    /// # Arguments
    ///
    /// * `filename` - The filepath to write the dot file to, as a String
    ///
    #[allow(dead_code)]
    fn write_dot(graph: Graph<u8, i32, Directed, usize>, filename: String) {
        let mut file = match File::create(&filename) {
            Err(why) => panic!("couldn't open file {}: {}", filename, why.description()),
            Ok(file) => file,
        };
        let g = graph.map(|_, nw| *nw as char, |_, ew| ew);
        match file.write_all(Dot::new(&g).to_string().as_bytes()) {
            Err(why) => panic!("couldn't write to file {}: {}", filename, why.description()),
            _ => (),
        }
    }

    #[test]
    fn test_init_graph() {
        // sanity check for String -> Graph

        let match_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let poa = Aligner::from_string(b"123456789", -1, match_fn);
        assert!(poa.graph.is_directed());
        assert_eq!(poa.graph.node_count(), 9);
        assert_eq!(poa.graph.edge_count(), 8);
    }

    #[test]
    fn test_alignment() {
        let match_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // examples from the POA paper
        //let _seq1 = b"PKMIVRPQKNETV";
        //let _seq2 = b"THKMLVRNETIM";
        let mut poa = Aligner::from_string(b"GATTACA", -1, match_fn);
        let alignment = poa.global(b"GCATGCU");
        assert_eq!(alignment.score, 0);

        let alignment = poa.global(b"GCATGCUx");
        assert_eq!(alignment.score, -1);

        let alignment = poa.global(b"xCATGCU");
        assert_eq!(alignment.score, -2);
    }

    #[test]
    fn test_branched_alignment() {
        let match_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let seq1 = b"TTTTT";
        let seq2 = b"TTATT";
        let mut poa = Aligner::from_string(seq1, -1, match_fn);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.global(seq2);
        assert_eq!(alignment.score, 3);
    }

    #[test]
    fn test_alt_branched_alignment() {
        let match_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let seq1 = b"TTCCTTAA";
        let seq2 = b"TTTTGGAA";
        let mut poa = Aligner::from_string(seq1, -1, match_fn);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.global(seq2);
        poa.add_alignment(alignment, seq2);
        assert_eq!(poa.graph.edge_count(), 14);
        assert!(
            poa.graph
                .contains_edge(NodeIndex::new(5), NodeIndex::new(10))
        );
        assert!(
            poa.graph
                .contains_edge(NodeIndex::new(11), NodeIndex::new(6))
        );
    }

    #[test]
    fn test_insertion_on_branch() {
        let match_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let seq1 = b"TTCCGGTTTAA";
        let seq2 = b"TTGGTATGGGAA";
        let seq3 = b"TTGGTTTGCGAA";
        let mut poa = Aligner::from_string(seq1, -1, match_fn);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'C');
        let node2 = poa.graph.add_node(b'C');
        let node3 = poa.graph.add_node(b'C');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, node3, 1);
        poa.graph.add_edge(node3, tail, 1);
        let alignment = poa.global(seq2);
        assert_eq!(alignment.score, 2);
        poa.add_alignment(alignment, seq2);
        let alignment2 = poa.global(seq3);

        assert_eq!(alignment2.score, 10);
    }
}
