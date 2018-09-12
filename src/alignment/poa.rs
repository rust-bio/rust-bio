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

use std::cmp::{max, Ordering};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::str;

use utils::TextSlice;

use petgraph::dot::Dot;
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
    score: i32,
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
    mode: AlignmentMode,
}

impl MatchFunc for MatchParams {
    #[inline]
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

impl<F> MatchFunc for F
where
    F: Fn(u8, u8) -> i32,
{
    fn score(&self, a: u8, b: u8) -> i32 {
        (self)(a, b)
    }
}

struct Scoring<F: MatchFunc> {
    gap_open: i32, // indel penalty
    match_fn: F,
}

impl Scoring<MatchParams> {
    pub fn from_scores(
        gap_open: i32,
        gap_extend: i32,
        match_score: i32,
        mismatch_score: i32,
    ) -> Self {
        Scoring {
            gap_open,
            match_fn: MatchParams::new(match_score, mismatch_score),
            match_scores: Some((match_score, mismatch_score)),
        }
    }
}

pub struct Aligner<F: MatchFunc> {
    scoring: Scoring<F>,
    traceback: Traceback,
    graph: Graph<u8, i32, Directed, usize>,
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
    fn init(&mut self, m: usize, n: usize) -> Self {
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
        for i in 1..(m + 1) {
            // TODO: these should be -1 * distance from head node
            traceback[i][0] = TracebackCell {
                score: -1 * i as i32,
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..(n + 1) {
            traceback[0][j] = TracebackCell {
                score: -1 * j as i32,
                op: AlignmentOperation::Ins(None),
            };
        }

        Traceback {
            rows: m,
            cols: n,
            matrix: matrix,
        }
    }
}

impl<F: MatchFunc<N>> Aligner<F, N, E> {
    pub fn from_string(reference: TextSlice, gap_open: i32, match_fn: F<N>) -> Self {
        let mut poa: Graph<N, E, Directed, usize> =
            Graph::with_capacity(reference.len(), reference.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(reference[0]);
        let mut node: NodeIndex<usize>;
        for i in 1..reference.len() {
            node = poa.add_node(reference[i]);
            graph.add_edge(prev, node, 1);
            prev = node;
        }

        Aligner {
            scoring: Scoring::new(gap_open, match_fn),
            traceback: Traceback::new(reference.len() + 1, n),
            poa: poa,
        }
    }

    pub fn from_graph(poa: Graph<N, E, Directed, usize>, gap_open: i32, match_fn: F<N>) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, match_fn),
            traceback: Traceback::new(),
            poa: poa,
        }
    }

    /// Naive Needleman-Wunsch
    /// Populates the traceback matrix in $O(n^2)$ time using
    /// petgraph's constant time topological traversal
    pub fn global(&mut self, query: TextSlice) -> Alignment {
        // dimensions of the traceback matrix
        let (m, n) = (g.node_count(), query.len());
        self.traceback.init(m, n);

        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        self.traceback[0][0] = TracebackCell {
            score: 0,
            op: AlignmentOperation::Match(None),
        };

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
                let max_cell = if prevs.len() == 0 {
                    TracebackCell {
                        score: self.traceback[0][j - 1].score + (self.scoring)(r, *q),
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
                                    score: self.traceback[i_p][j - 1].score + (self.scoring)(r, *q),
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: self.traceback[i_p][j].score - 1i32,
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
                        score: self.traceback[i][j - 1].score - 1i32,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                self.traceback[i][j] = score;
            }
        }

        //dump_traceback(&self.traceback, g, query);

        // Now backtrack through the matrix to construct an optimal path
        let mut i = last.index() + 1;
        let mut j = n;

        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.traceback[i][j].op.clone());
            match self.traceback[i][j].op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j = j - 1;
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
                    j = j - 1;
                }
                AlignmentOperation::Ins(None) => {
                    i = i - 1;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: self.traceback[last.index() + 1][n].score,
            operations: ops,
        }
    }

    /// Incorporate a new sequence into a graph from an alignment
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
                AlignmentOperation::Match(None) => {
                    i = i + 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
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
                }
                AlignmentOperation::Ins(None) => {
                    i = i + 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i = i + 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
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
fn write_dot(graph: Graph, filename: String) {
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

// print out a traceback matrix
#[allow(dead_code)]
fn dump_traceback(
    traceback: &Vec<Vec<TracebackCell>>,
    g: &Graph<u8, i32, Directed, usize>,
    query: TextSlice,
) {
    let (m, n) = (g.node_count(), query.len());
    print!(".\t");
    for i in 0..n {
        print!("{:?}\t", query[i]);
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
