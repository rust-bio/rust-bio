// Copyright 2017-2025 Brett Bowman, Jeff Knaggs, Minindu Weerakoon
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of multiple homologous sequences.
//!
//! - time complexity: `O(N^2 * L^2)`, where `N` is the number of sequences and `L` is the length of each sequence.
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
//! use bio::alignment::pairwise::Scoring;
//! use bio::alignment::poa::*;
//!
//! let x = b"AAAAAAA";
//! let y = b"AABBBAA";
//! let z = b"AABCBAA";
//!
//! let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
//! let mut aligner = Aligner::new(scoring, x);
//! // z differs from x in 3 locations
//! assert_eq!(aligner.global(z).alignment().score, 1);
//! aligner.global(y).add_to_graph();
//! // z differs from x and y's partial order alignment by 1 base
//! assert_eq!(aligner.global(z).alignment().score, 5);
//! ```

use std::cmp::max;

use crate::utils::TextSlice;

use crate::alignment::pairwise::{MatchFunc, Scoring};

use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs
pub type POAGraph = Graph<u8, i32, Directed, usize>;

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
    Xclip(usize),
    Yclip(usize, usize), // to, from
}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Alignment {
    pub score: i32,
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
}

impl Alignment {
    /// Return the pretty formatted poa alignment as a String. The string
    /// contains sets of (number of queries + 1) lines of length (ncol). First line is for the
    /// consensus sequence, rest are for the queries. A '-' in the sequence
    /// indicates a blank (insertion/deletion).
    ///
    /// # Example
    ///
    /// If we align the strings,
    /// "ACCCCCTTTTTCCGG"
    /// "ACTTCCCTTTTTCCGG"
    /// "ACCGCCTTTTTCCGG"
    /// "ACCCCCTGTTTCAAGG"
    /// we will get the following output:
    ///
    /// ```c
    /// cons:   ACC---CCCTT-TTTCC--GG
    /// seq1:   ACC---CCCTT-TTTCC--GG
    /// seq2:   AC--TTCCCTT-TTTCC--GG
    /// seq3:   ACCG--CC-TT-TTTCC--GG
    /// seq4:   ACC---CCCT-GTTTC-AAGG
    /// ```
    pub fn pretty(
        &self,
        consensus: TextSlice,
        sequences: Vec<TextSlice>,
        graph: &POAGraph,
        ncol: usize,
    ) -> String {
        // initialize required variables
        let mut seq_indices: Vec<usize> = vec![0; sequences.len()];
        let mut seq_pretty: Vec<Vec<u8>> = vec![vec![]; sequences.len()];
        let mut con_index: usize = 0;
        let mut con_pretty: Vec<u8> = vec![];
        let mut topo = Topo::new(graph);
        // go through the nodes topologically
        while let Some(node) = topo.next(graph) {
            let topo_base = graph.raw_nodes()[node.index()].weight;
            let mut all_null = true;
            // go through the sequences checking if bases match with the current node base
            for current_seq in 0..sequences.len() {
                if seq_indices[current_seq] >= sequences[current_seq].len() {
                    seq_pretty[current_seq].push(b'-');
                } else if sequences[current_seq][seq_indices[current_seq]] == topo_base {
                    seq_pretty[current_seq].push(topo_base);
                    seq_indices[current_seq] += 1;
                    all_null = false;
                } else {
                    seq_pretty[current_seq].push(b'-');
                }
            }
            // do the same for the consensus
            if con_index >= consensus.len() {
                con_pretty.push(b'-');
            } else if consensus[con_index] == topo_base {
                con_pretty.push(topo_base);
                con_index += 1;
                all_null = false;
            } else {
                con_pretty.push(b'-');
            }
            if all_null {
                for current_seq in seq_pretty.iter_mut() {
                    current_seq.pop();
                }
                con_pretty.pop();
            }
        }
        let mut s = String::new();
        let mut idx = 0;
        use std::cmp::min;

        let ml = con_pretty.len();

        while idx < ml {
            let rng = idx..min(idx + ncol, ml);
            s.push_str(&format!(
                "cons:\t{}\n",
                &String::from_utf8_lossy(&con_pretty[rng.clone()])
            ));
            for (seq_index, seq) in seq_pretty.iter().enumerate() {
                s.push_str(&format!(
                    "seq{}:\t{}\n",
                    seq_index + 1,
                    &String::from_utf8_lossy(&seq[rng.clone()])
                ));
            }
            s.push('\n');
            idx += ncol;
        }
        s
    }
}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct Traceback {
    rows: usize,
    cols: usize,
    // store these for xclipping and yclipping
    best_in_last_row: usize,
    best_in_last_col: usize,
    best_overall: (usize, usize),
    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    // start and end of each row in matrix for banded poa
    start_end_vec: Vec<(usize, usize)>,
    // score matrix
    matrix: Vec<Vec<i32>>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    fn with_capacity(m: usize, n: usize) -> Self {
        // each row of matrix contain start end position and vec of traceback cells
        let matrix: Vec<Vec<i32>> = vec![vec![]; m + 1];
        let start_end_vec = vec![(0, n + 1); m + 1];
        Traceback {
            rows: m,
            cols: n,
            best_in_last_row: 0,
            best_in_last_col: 0,
            best_overall: (0, 0),
            last: NodeIndex::new(0),
            start_end_vec: start_end_vec,
            matrix,
        }
    }
    /// Populate the first row of the traceback matrix
    fn initialize_scores(&mut self, gap_open: i32, yclip_prefix: i32) {
        for j in 0..=self.cols {
            self.matrix[0].push(max((j as i32) * gap_open, yclip_prefix));
        }
        self.matrix[0][0] = 0;
    }

    fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            best_in_last_row: 0,
            best_in_last_col: 0,
            best_overall: (0, 0),
            last: NodeIndex::new(0),
            start_end_vec: Vec::new(),
            matrix: Vec::new(),
        }
    }

    // create a new row according to the parameters
    fn new_row(
        &mut self,
        row: usize,
        size: usize,
        gap_open: i32,
        xclip_prefix: i32,
        start: usize,
        end: usize,
    ) {
        self.start_end_vec[row] = (start, end);
        // when the row starts from the edge
        if start == 0 {
            self.matrix[row].push(max((row as i32) * gap_open, xclip_prefix));
        } else {
            self.matrix[row].push(MIN_SCORE);
        }
        for _ in 1..=size {
            self.matrix[row].push(MIN_SCORE);
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: i32) {
        // set the matrix cell if in band range
        if !(self.start_end_vec[i].0 > j || self.start_end_vec[i].1 < j) {
            let real_position = j - self.start_end_vec[i].0;
            self.matrix[i][real_position] = cell;
        }
    }

    fn get(&self, i: usize, j: usize) -> i32 {
        // get the matrix cell if in band range else return the appropriate values
        if !(self.start_end_vec[i].0 > j
            || self.start_end_vec[i].1 <= j
            || self.matrix[i].is_empty())
        {
            let real_position = j - self.start_end_vec[i].0;
            self.matrix[i][real_position]
        }
        // out of band
        else {
            MIN_SCORE
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
#[derive(Default, Clone, Debug)]
pub struct Aligner<F: MatchFunc> {
    traceback: Traceback,
    query: Vec<u8>,
    poa: Poa<F>,
}

impl<F: MatchFunc> Aligner<F> {
    /// Create new instance.
    pub fn new(scoring: Scoring<F>, reference: TextSlice) -> Self {
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(scoring, reference),
        }
    }

    /// Get the alignment of the last query to the graph and add to graph.
    pub fn add_to_graph(&mut self) -> &mut Self {
        let alignment = self.poa.recalculate_alignment(&self.traceback);
        self.poa.add_alignment(&alignment, &self.query);
        self
    }

    /// Return alignment of last added query against the graph.
    pub fn alignment(&self) -> Alignment {
        self.poa.recalculate_alignment(&self.traceback)
    }

    /// Add the alignment to the graph
    pub fn add_alignment(&mut self, alignment: &Alignment) -> &mut Self {
        self.poa.add_alignment(&alignment, &self.query);
        self
    }

    /// Globally align a given query against the graph.
    pub fn global(&mut self, query: TextSlice) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.scoring.xclip_prefix,
            self.poa.scoring.xclip_suffix,
            self.poa.scoring.yclip_prefix,
            self.poa.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.scoring.xclip_prefix = MIN_SCORE;
        self.poa.scoring.xclip_suffix = MIN_SCORE;
        self.poa.scoring.yclip_prefix = MIN_SCORE;
        self.poa.scoring.yclip_suffix = MIN_SCORE;

        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);

        // Set the clip penalties to the original values
        self.poa.scoring.xclip_prefix = clip_penalties[0];
        self.poa.scoring.xclip_suffix = clip_penalties[1];
        self.poa.scoring.yclip_prefix = clip_penalties[2];
        self.poa.scoring.yclip_suffix = clip_penalties[3];

        self
    }

    /// Semi-globally align a given query against the graph.
    pub fn semiglobal(&mut self, query: TextSlice) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.scoring.xclip_prefix,
            self.poa.scoring.xclip_suffix,
            self.poa.scoring.yclip_prefix,
            self.poa.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.scoring.xclip_prefix = MIN_SCORE;
        self.poa.scoring.xclip_suffix = MIN_SCORE;
        self.poa.scoring.yclip_prefix = 0;
        self.poa.scoring.yclip_suffix = 0;

        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);

        // Set the clip penalties to the original values
        self.poa.scoring.xclip_prefix = clip_penalties[0];
        self.poa.scoring.xclip_suffix = clip_penalties[1];
        self.poa.scoring.yclip_prefix = clip_penalties[2];
        self.poa.scoring.yclip_suffix = clip_penalties[3];

        self
    }

    /// Locally align a given query against the graph.
    pub fn local(&mut self, query: TextSlice) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.scoring.xclip_prefix,
            self.poa.scoring.xclip_suffix,
            self.poa.scoring.yclip_prefix,
            self.poa.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.scoring.xclip_prefix = 0;
        self.poa.scoring.xclip_suffix = 0;
        self.poa.scoring.yclip_prefix = 0;
        self.poa.scoring.yclip_suffix = 0;

        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);

        // Set the clip penalties to the original values
        self.poa.scoring.xclip_prefix = clip_penalties[0];
        self.poa.scoring.xclip_suffix = clip_penalties[1];
        self.poa.scoring.yclip_prefix = clip_penalties[2];
        self.poa.scoring.yclip_suffix = clip_penalties[3];

        self
    }

    /// Custom align a given query against the graph with custom xclip and yclip penalties.
    pub fn custom(&mut self, query: TextSlice) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);
        self
    }

    /// Globally align a given query against the graph with a band around the previous
    /// optimal score for speed.
    pub fn global_banded(&mut self, query: TextSlice, bandwidth: usize) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.global_banded(query, bandwidth);
        self
    }

    /// Return alignment graph.
    pub fn graph(&self) -> &POAGraph {
        &self.poa.graph
    }
    /// Return the consensus sequence generated from the POA graph.
    pub fn consensus(&self) -> Vec<u8> {
        let mut consensus: Vec<u8> = vec![];
        let max_index = self.poa.graph.node_count();
        let mut weight_score_next_vec: Vec<(i32, i32, usize)> = vec![(0, 0, 0); max_index + 1];
        let mut topo = Topo::new(&self.poa.graph);
        // go through the nodes topologically
        while let Some(node) = topo.next(&self.poa.graph) {
            let mut best_weight_score_next: (i32, i32, usize) = (0, 0, usize::MAX);
            let neighbour_nodes = self.poa.graph.neighbors_directed(node, Incoming);
            // go through the incoming neighbour nodes
            for neighbour_node in neighbour_nodes {
                let neighbour_index = neighbour_node.index();
                let neighbour_score = weight_score_next_vec[neighbour_index].1;
                let edges = self.poa.graph.edges_connecting(neighbour_node, node);
                let weight = edges.map(|edge| edge.weight()).sum();
                let current_node_score = weight + neighbour_score;
                // save the neighbour node with the highest weight and score as best
                if (weight, current_node_score, neighbour_index) > best_weight_score_next {
                    best_weight_score_next = (weight, current_node_score, neighbour_index);
                }
            }
            weight_score_next_vec[node.index()] = best_weight_score_next;
        }
        // get the index of the max scored node (end of consensus)
        let mut pos = weight_score_next_vec
            .iter()
            .enumerate()
            .max_by_key(|(_, &value)| value.1)
            .map(|(idx, _)| idx)
            .unwrap();
        // go through weight_score_next_vec appending to the consensus
        while pos != usize::MAX {
            consensus.push(self.poa.graph.raw_nodes()[pos].weight);
            pos = weight_score_next_vec[pos].2;
        }
        consensus.reverse();
        consensus
    }
}

/// A partially ordered alignment graph
///
/// A directed acyclic graph datastructure that represents the topology of a
/// traceback matrix.
#[derive(Default, Clone, Debug)]
pub struct Poa<F: MatchFunc> {
    scoring: Scoring<F>,
    pub graph: POAGraph,
}

impl<F: MatchFunc> Poa<F> {
    /// Create a new aligner instance from the directed acyclic graph of another.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `poa` - the partially ordered reference alignment
    pub fn new(scoring: Scoring<F>, graph: POAGraph) -> Self {
        Poa { scoring, graph }
    }

    /// Create a new POA graph from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    pub fn from_string(scoring: Scoring<F>, seq: TextSlice) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }

        Poa { scoring, graph }
    }
    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    pub fn custom(&self, query: TextSlice) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        // save score location of the max scoring node for the query for suffix clipping
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.scoring.gap_open, self.scoring.yclip_prefix);
        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        // required for semi global/ local alignment
        let mut max_score_last_column = i32::MIN; // score, row
        let mut max_score_overall = 0; // score, row, column
        while let Some(node) = topo.next(&self.graph) {
            let mut max_score_last_row = i32::MIN; // score, column
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1; // 0 index is for initialization so we start at 1
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            traceback.new_row(
                i,
                n + 1,
                self.scoring.gap_open,
                self.scoring.xclip_prefix,
                0,
                n + 1,
            );
            // get y clip min
            let y_clip_min = traceback.get(i, 0) + self.scoring.yclip_prefix;
            // query base and its index in the DAG (traceback matrix rows)
            for (query_index, query_base) in query.iter().enumerate() {
                let j = query_index + 1; // 0 index is initialized so we start at 1
                let max_cell = if prevs.is_empty() {
                    traceback.get(0, j - 1) + self.scoring.match_fn.score(r, *query_base)
                } else {
                    // get x clip min
                    let x_clip_min = traceback.get(0, j) + self.scoring.xclip_prefix;
                    let mut max_cell = max(MIN_SCORE, max(x_clip_min, y_clip_min));
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        max_cell = max(
                            max_cell,
                            max(
                                traceback.get(i_p, j - 1)
                                    + self.scoring.match_fn.score(r, *query_base),
                                traceback.get(i_p, j) + self.scoring.gap_open,
                            ),
                        );
                    }
                    max_cell
                };
                let score = max(max_cell, traceback.get(i, j - 1) + self.scoring.gap_open);
                if score > max_score_last_row {
                    max_score_last_row = score;
                    traceback.best_in_last_row = j;
                }
                if (score > max_score_last_column) && (query_index == query.len() - 1) {
                    max_score_last_column = score;
                    traceback.best_in_last_col = i;
                }
                if score > max_score_overall {
                    max_score_overall = score;
                    traceback.best_overall = (i, j);
                }
                traceback.set(i, j, score);
            }
        }
        traceback
    }
    /// A global Needleman-Wunsch aligner on partially ordered graphs with banding.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    /// * `bandwidth` - width of band, if too small, alignment may be suboptimal
    pub fn global_banded(&self, query: TextSlice, bandwidth: usize) -> Traceback {
        assert!(self.graph.node_count() != 0);

        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.scoring.gap_open, self.scoring.yclip_prefix);

        traceback.set(0, 0, 0);

        // construct the score matrix (O(n^2) space)
        // but this sucks, we want linear time!!!
        // at each row i we want to find the max scoring j
        // and band
        let mut max_scoring_j = 0;
        let mut max_score_for_row = MIN_SCORE;
        let mut topo = Topo::new(&self.graph);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1; // 0 index is for initialization so we start at 1
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            let start = if bandwidth > max_scoring_j {
                0
            } else {
                max_scoring_j - bandwidth
            };
            let end = max_scoring_j + bandwidth;
            traceback.new_row(
                i,
                (end - start) + 1,
                self.scoring.gap_open,
                self.scoring.xclip_prefix,
                start,
                end + 1,
            );
            for (query_index, query_base) in query.iter().enumerate().skip(start) {
                let j = query_index + 1; // 0 index is initialized so we start at 1
                if j > end {
                    break;
                }
                let max_cell = if prevs.is_empty() {
                    traceback.get(0, j - 1) + self.scoring.match_fn.score(r, *query_base)
                } else {
                    let mut max_cell = MIN_SCORE;
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        max_cell = max(
                            max_cell,
                            max(
                                traceback.get(i_p, j - 1)
                                    + self.scoring.match_fn.score(r, *query_base),
                                traceback.get(i_p, j) + self.scoring.gap_open,
                            ),
                        );
                    }
                    max_cell
                };

                let score = max(max_cell, traceback.get(i, j - 1) + self.scoring.gap_open);
                if score > max_score_for_row {
                    max_scoring_j = j;
                    max_score_for_row = score;
                }
                traceback.set(i, j, score);
            }
        }
        traceback
    }
    /// Computes the alignment using the provided traceback structure.
    /// This function requires the original graph used to generate the traceback score matrix.
    ///
    /// time complexity: `O(N * L)`, where `N` is the number of sequences and `L` is the length of each sequence.
    ///
    /// # Arguments
    /// * `traceback` - the traceback struct
    pub fn recalculate_alignment(&self, traceback: &Traceback) -> Alignment {
        // Get the alignment by backtracking and recalculating stuff
        let mut ops: Vec<AlignmentOperation> = vec![];
        // loop until we reach a node with no incoming
        let last_node = traceback.last.index() + 1;
        let last_query = traceback.cols;
        let final_score = traceback.get(last_node, last_query);

        let mut curr_node = traceback.last.index() + 1;
        let mut curr_query = traceback.cols;
        let xy_score = traceback.get(traceback.best_overall.0, traceback.best_overall.1)
            + self.scoring.xclip_suffix
            + self.scoring.yclip_suffix;
        let y_score =
            traceback.get(last_node, traceback.best_in_last_row) + self.scoring.yclip_suffix;
        let x_score =
            traceback.get(traceback.best_in_last_col, last_query) + self.scoring.xclip_suffix;
        // x and y suffix make sure not in last row and column
        if (xy_score >= final_score)
            && (xy_score >= x_score)
            && (xy_score >= y_score)
            && (traceback.best_overall.1 != last_query)
            && (traceback.best_overall.0 != last_node)
        {
            ops.push(AlignmentOperation::Xclip(traceback.best_overall.0));
            ops.push(AlignmentOperation::Yclip(
                traceback.best_overall.1,
                last_query,
            ));
            curr_node = traceback.best_overall.0;
            curr_query = traceback.best_overall.1;
        }
        // y suffix
        else if (y_score >= final_score)
            && (y_score >= x_score)
            && (traceback.best_in_last_row != last_query)
        {
            ops.push(AlignmentOperation::Yclip(
                traceback.best_in_last_row,
                last_query,
            ));
            curr_query = traceback.best_in_last_row;
        }
        // x suffix
        else if (x_score >= final_score) && (traceback.best_in_last_col != last_node) {
            ops.push(AlignmentOperation::Xclip(traceback.best_in_last_col));
            curr_node = traceback.best_in_last_col;
        }
        loop {
            let mut current_alignment_operation = AlignmentOperation::Match(None);
            let current_cell_score = traceback.get(curr_node, curr_query);
            let mut next_jump = curr_query;
            let mut next_node = 1;
            // check left if gap open difference with left
            let prevs: Vec<NodeIndex<usize>> = self
                .graph
                .neighbors_directed(NodeIndex::new(curr_node - 1), Incoming)
                .collect();
            let mut jump_up_score = MIN_SCORE;
            let mut jump_diagonal_score = MIN_SCORE;
            let jump_left_score = traceback.get(curr_node, curr_query - 1) + self.scoring.gap_open;
            if current_cell_score == jump_left_score {
                current_alignment_operation = AlignmentOperation::Ins(Some(curr_node - 1));
                next_node = curr_node;
                next_jump = curr_query - 1;
            } else {
                for prev in &prevs {
                    let prev_node = prev.index() + 1;
                    let diagonal_score = traceback.get(prev_node, curr_query - 1);
                    let top_score = traceback.get(prev_node, curr_query);
                    // Top
                    if current_cell_score == top_score + self.scoring.gap_open {
                        jump_up_score =
                            traceback.get(prev_node, curr_query) + self.scoring.gap_open;
                        current_alignment_operation = AlignmentOperation::Del(None);
                        next_jump = curr_query;
                        next_node = prev_node;
                    }
                    // Diagonal mismatch
                    else if current_cell_score
                        == diagonal_score + self.scoring.match_fn.score(0, 1)
                    {
                        jump_diagonal_score = traceback.get(prev_node, curr_query - 1)
                            + self.scoring.match_fn.score(0, 1);
                        current_alignment_operation =
                            AlignmentOperation::Match(Some((prev_node - 1, curr_node - 1)));
                        next_node = prev_node;
                        next_jump = curr_query - 1;
                    }
                    // Diagonal match
                    else if current_cell_score
                        == diagonal_score + self.scoring.match_fn.score(0, 0)
                    {
                        jump_diagonal_score = traceback.get(prev_node, curr_query - 1)
                            + self.scoring.match_fn.score(0, 0);
                        current_alignment_operation =
                            AlignmentOperation::Match(Some((prev_node - 1, curr_node - 1)));
                        next_node = prev_node;
                        next_jump = curr_query - 1;
                    }
                }
                // Last node to avoid Match(None) which should never occur
                if prevs.len() == 0 {
                    // match
                    if current_cell_score
                        == traceback.get(0, curr_query - 1) + self.scoring.match_fn.score(0, 0)
                    {
                        current_alignment_operation = AlignmentOperation::Match(None);
                        jump_diagonal_score =
                            traceback.get(0, curr_query - 1) + self.scoring.match_fn.score(0, 0);
                        next_node = 1;
                        next_jump = curr_query - 1;
                    }
                    // mismatch
                    if current_cell_score
                        == traceback.get(0, curr_query - 1) + self.scoring.match_fn.score(0, 1)
                    {
                        current_alignment_operation = AlignmentOperation::Match(None);
                        jump_diagonal_score =
                            traceback.get(0, curr_query - 1) + self.scoring.match_fn.score(0, 1);
                        next_node = 1;
                        next_jump = curr_query - 1;
                    }
                }
            }
            // max score of the three
            let max_score = max(jump_diagonal_score, max(jump_up_score, jump_left_score));
            // x clip
            if self.scoring.xclip_prefix >= max_score {
                next_node = 0;
                current_alignment_operation = AlignmentOperation::Xclip(0);
            }
            // y clip
            if self.scoring.yclip_prefix >= max(max_score, self.scoring.xclip_prefix) {
                next_jump = 0;
                current_alignment_operation = AlignmentOperation::Yclip(0, curr_query);
            }
            ops.push(current_alignment_operation);
            // iterate to next
            curr_query = next_jump;
            curr_node = next_node;
            // break point
            if prevs.len() == 0 || curr_query == 0 {
                //if at end but not at start of query add bunch of ins(None)
                if prevs.len() == 0 {
                    if curr_query > 0 {
                        for _ in 0..curr_query {
                            if self.scoring.yclip_prefix > MIN_SCORE {
                                ops.push(AlignmentOperation::Yclip(0, curr_query));
                                break;
                            } else {
                                ops.push(AlignmentOperation::Ins(None));
                            }
                        }
                    }
                } else {
                    //if at start of query but prev nodes available
                    if self.scoring.xclip_prefix > MIN_SCORE {
                        ops.push(AlignmentOperation::Xclip(0));
                    } else {
                        ops.push(AlignmentOperation::Del(None));
                    }
                }
                break;
            }
        }
        ops.reverse();
        Alignment {
            score: final_score as i32,
            operations: ops,
        }
    }
    /// Experimental: return sequence of traversed edges
    ///
    /// Only supports alignments for sequences that have already been added,
    /// so all operations must be Match.
    pub fn edges(&self, aln: Alignment) -> Vec<usize> {
        let mut path: Vec<usize> = vec![];
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut _i: usize = 0;
        for op in aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    _i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(p);
                    let edge = self.graph.find_edge(prev, node).unwrap();
                    path.push(edge.index());
                    prev = NodeIndex::new(p);
                    _i += 1;
                }
                AlignmentOperation::Ins(None) => {}
                AlignmentOperation::Ins(Some(_)) => {}
                AlignmentOperation::Del(_) => {}
                AlignmentOperation::Xclip(_) => {}
                AlignmentOperation::Yclip(_, _) => {}
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
    pub fn add_alignment(&mut self, aln: &Alignment, seq: TextSlice) {
        let head = Topo::new(&self.graph).next(&self.graph).unwrap();
        let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
        let mut i: usize = 0;
        let mut edge_not_connected: bool = false;
        for op in &aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    let node: NodeIndex<usize> = NodeIndex::new(head.index());
                    if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        if edge_not_connected {
                            self.graph.add_edge(prev, node, 1);
                        }
                        edge_not_connected = false;
                        prev = node;
                    }
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        edge_not_connected = false;
                    }
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
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
                                if prev.index() != head.index() && prev.index() != node.index() {
                                    self.graph.add_edge(prev, node, 1);
                                }
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    let node = self.graph.add_node(seq[i]);
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                    }
                    prev = node;
                    edge_not_connected = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes and xclip
                AlignmentOperation::Xclip(_) => {}
                AlignmentOperation::Yclip(_, r) => {
                    i = *r;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::pairwise::Scoring;
    use petgraph::dot::Dot;
    use petgraph::graph::NodeIndex;

    #[test]
    fn test_init_graph() {
        // sanity check for String -> Graph

        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let poa = Poa::from_string(scoring, b"123456789");
        assert!(poa.graph.is_directed());
        assert_eq!(poa.graph.node_count(), 9);
        assert_eq!(poa.graph.edge_count(), 8);
    }

    #[test]
    fn test_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        // examples from the POA paper
        //let _seq1 = b"PKMIVRPQKNETV";
        //let _seq2 = b"THKMLVRNETIM";
        let poa = Poa::from_string(scoring, b"GATTACA");
        let traceback = poa.custom(b"GCATGCU");
        let alignment = poa.recalculate_alignment(&traceback);
        assert_eq!(alignment.score, 0);

        let traceback = poa.custom(b"GCATGCUx");
        let alignment = poa.recalculate_alignment(&traceback);
        assert_eq!(alignment.score, -1);

        let traceback = poa.custom(b"xCATGCU");
        let alignment = poa.recalculate_alignment(&traceback);
        assert_eq!(alignment.score, -2);
    }

    #[test]
    fn test_branched_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let seq1 = b"TTTTT";
        let seq2 = b"TTATT";
        let mut poa = Poa::from_string(scoring, seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let traceback = poa.custom(seq2);
        let alignment = poa.recalculate_alignment(&traceback);
        assert_eq!(alignment.score, 3);
    }

    #[test]
    fn test_alt_branched_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

        let seq1 = b"TTCCTTAA";
        let seq2 = b"TTTTGGAA";
        let mut poa = Poa::from_string(scoring, seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let traceback = poa.custom(seq2);
        let alignment = poa.recalculate_alignment(&traceback);
        poa.add_alignment(&alignment, seq2);
        assert_eq!(poa.graph.edge_count(), 14);
        assert!(poa
            .graph
            .contains_edge(NodeIndex::new(5), NodeIndex::new(10)));
        assert!(poa
            .graph
            .contains_edge(NodeIndex::new(11), NodeIndex::new(6)));
    }

    #[test]
    fn test_insertion_on_branch() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

        let seq1 = b"TTCCGGTTTAA";
        let seq2 = b"TTGGTATGGGAA";
        let seq3 = b"TTGGTTTGCGAA";
        let mut poa = Poa::from_string(scoring, seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'C');
        let node2 = poa.graph.add_node(b'C');
        let node3 = poa.graph.add_node(b'C');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, node3, 1);
        poa.graph.add_edge(node3, tail, 1);
        let traceback = poa.custom(seq2);
        let alignment = poa.recalculate_alignment(&traceback);
        assert_eq!(alignment.score, 2);
        poa.add_alignment(&alignment, seq2);
        let traceback = poa.custom(seq3);
        let alignment2 = poa.recalculate_alignment(&traceback);

        assert_eq!(alignment2.score, 10);
    }

    #[test]
    fn test_poa_method_chaining() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"TTCCGGTTTAA");
        aligner
            .global(b"TTGGTATGGGAA")
            .add_to_graph()
            .global(b"TTGGTTTGCGAA");
        assert_eq!(aligner.alignment().score, 10);
    }

    #[test]
    fn test_global_banded() {
        // need strings long enough for the band to matter
        // create MSA with and without band and make sure they are the same
        let s1 = b"TGGCATGCTCAAGGACCGTTGAATACTATCTTAATGGACCGCAAGCTCCCTGAAGGTGGGCCACATTCGAGGGCC\
        CGGCCTCCACCTATTCCCAACGAAACTAGCATTAACATGGACAGGGGCGCATAAAACAGAGTTTCTCCTAATCCCCTTTCCCCTG\
        GAGTGCTAGTCAGAACCGCACATGTTGACGCTTTGGTCAGGTGTAGCCGATTCACTACCCGGGGTAGTACGAGTGGTAGCACCAT\
        GGTTAGCTTCTCCGGGATGTTCCGCGAAGAGAGCGGAGCGGGCGTGCACAAGCTCGGACAACCCTAGTGTGCATCAAATGCCATA\
        TGTTCTGCTTTGTCTGTGACTCACGCCCACGTTTGACATCACTCTTACTATCCAACGGGCCAAGCTTAGGAGGGGCGGACCTATT\
        GAACCATTAGAGGGGATCCTTCTGAAGTTAAGGCACAGCGTTGAGGGGCTATAGTCGATCCTCTTAGTAAATATAATGGACAGGT\
        CTTTACGACACAGTATGAATTAGTCCAATGGAGCCATTGTAATCGATGAAACTGTTATATCTGTTGGCCTAGTCGCAACGGTCTA\
        CATCGCTAGCGTAACGGTTAAGACCTCTTCCACGAGTGGGACACTCATAAAGCTCGCGGCCCTTACGATCTAGGGGAGCGCACTC\
        CGTAGTCAATCACGGCCAGCCGGTGTGCGCTAAGTTACGAAACAGTCACGAGCGATGAACCGTATGAAGAATGGACCCTTCTAAG\
        ATGTGAACACCTAGATGAGCAGCAAGACAATTGCTCTCGCCGACTCGTTCGAAAGTGTACCTCGAGAGCACAACACGCATTACCC\
        AGGTGACCGTGTATTGACTGCCCGTTACTCAGAAACCTTACAGTATTAATCGCCTAGTCTGTATAGTATTCATTCTGCCCGTGAC
        ATGCGGGAAGCCTGCTGAGATTGGCAGCGTCTTTGGAGGGTTACCAAGCGAGGACACGGGCAAATTGAGGTGT";
        let s2 = b"TGGCTACATGCTCAAGCATCGTTGAAGCTCATCTTAATGGACCGCAACGGCCGCCTGAAGGTGGGACACGTGACG\
        GGCGGGGGCCCGGCCTTAACCCATTCTCAAGCAACTAGCATACTGGACAGCGGCGCATATACAGAGAATCGCCTAAACCCACTTT\
        TGCCCTGAGTGCTAGTCAGTCCCCCACATCTGACACTTCGGTGGCGCACGTTTAGCAGCTTACACTACCCGGGGCAGTACGAGTG\
        CTAGCACGGTAGCCTCCGGAGGGCTGCGAGGAATAGAACGGAGAGGGCGTCCTCAAGCCGGACAACCCTAGTGTGCATCAAATGA\
        TGCCTGCTGATTTTCTGTGCATTTCACGCCCAATTCACAATCACTCCTACTATCCAACGGGCAAGCATAAGGGAGGGGGGGAGTA\
        CGTCTATTGCACCATTAGAGGGGTACTTCGAATTCGTTGAACTGAGATAGAGTCGATCCTCTTTGTATATAAACGCAGGTACTTT\
        GCTATAAGGTGAATTATTCAAATGGAGCCATTGTAATCGATGACAATGTTATACCTTTAGGCCTAGATCAACGGTCTCCATCGCA\
        AGCGTAACGATTATGACCCAACGAGTGGACACTCATAGAGCGGCCCTTACGAGCTAGGCGAGCGCAATCCGTGTGAATCACAGCC\
        AGACGGGATGTTGCGTTAAGCTACGAAACATCACGCGGTGAGCGTATGAATAATGGACCCGTCAAAATGTGGGCAGCGAGCAGCA\
        GGACAATTGCTCGGGTCCGGTAGCGACTCGTTCGAAACTGTAAACGTCGAGGCACAACACGCATTAGCCAGGTGAACATGTATTG\
        ACCGCCCCGTAGATACCTTACAGTATTAATCGCCTAGTCTGTATGATCTTCGTTCTGCCTGTGAACATGAGGGAAGCCTGCTTGA\
        GTTTGGCAGCGTCTTTGGGGTTTCCAAGCGAGCGACACGGGCAAATTGAGGTGT";
        let s3 = b"TGGCATGCTCAAGGAGTGCTGAAGCTCATTTTAATGGACCGCAACGGCCGCCTGAAGGTGGGGCACGTGACGGGC\
        GAGGGCCCGGCCTTAACCCATTCTCAAGAAACTCGTATACTGGACAGCGGCGCATATAGAGAGATTCTCCTAAACCCTCTTTTGC\
        CCTGACATGTGCTAGTCAGTCGCCCACATCTGAACACTTCGGCAGCGCACGTCTAGCAGCTTACACTACCGGGCGGGGCAGGTAC\
        GAGTGCTAGCACGGTAGCCTCTCCCGGAGTGCTGGGAATAGAAGGGAGAGGGCGTCCTCATGCCGGCGACCCTAGTGTGCATCAA\
        ATGAGATGCCTGCTGTGATTTTCACATTCACAATCACTCTTACCATCCCAACGGGACAAGCATAAGGGAGGGGGGGAGCTATTGA\
        ACCAAGAGGGGTCCTCCGGAATTCGTTGAGCTGCGATAGAGTCGATCCTCTTTGTATATAAACGCAGGTACTTTGCGATTAGGTG\
        AAGTATTCAAATGGAGCCATTGTAATCGATGACAATGTGATGCCTTTAGGCCTAGATCACGGTCTACATCGCGTAAGCGTAACGA\
        TTATGACCCAACGAGAGGCACACTCATAAAGCGCGGCCCTTACTAGCTAGGCGAGCTTAGCAATCCGTGCAATCACACCCAGACG\
        GGTTGAGCTAAGCTACGGAACACCACGCGATGAGCCGTATGAAGAATGGACCCGTCGAAAATGTGGACAGCGAGCATCAGGACAA\
        TTGCTCGGGTCCGCGACTCGTGCGGAACTGTAAACGTCGAGGCACAACACGATTAGCCAGGTGAACATGTAGACCGCCCCGTAGA\
        TATTTTACAGTATTAATCGCCTAGTCTGTATAGGATCTTCGTTCTGCCTGTGAACATGCGGGAAGCCTGCTTGAGATTGGCAGCG\
        TCTTTGGGCAAGCGAGGACACGGGCAAATCGAGGTGG";
        let scoring = Scoring::from_scores(-2, -2, 2, -4);
        let mut aligner_banded = Aligner::new(scoring, s1);
        aligner_banded.global_banded(s2, 25).add_to_graph();
        aligner_banded.global_banded(s3, 25);
        let scoring = Scoring::from_scores(-2, -2, 2, -4);
        let mut aligner_unbanded = Aligner::new(scoring, s1);
        aligner_unbanded.global(s2).add_to_graph();
        aligner_unbanded.global(s3);
        let alignment_banded = aligner_banded.alignment();
        let alignment_unbanded = aligner_unbanded.alignment();
        for (i, operation) in alignment_banded.operations.iter().enumerate() {
            assert_eq!(*operation, alignment_unbanded.operations[i]);
        }
    }

    #[test]
    fn test_edge_cases() {
        // case 1
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"BBA");
        aligner.global(b"AAA").add_to_graph();
        let g = aligner.graph().map(|_, n| (*n) as char, |_, e| *e);
        let dot = format!("{:?}", Dot::new(&g));
        let output = "digraph {\n    0 [ label = \"'B'\" ]\n    1 [ label = \"'B'\" ]\n    2 [ label = \"'A'\" ]\
        \n    3 [ label = \"'A'\" ]\n    4 [ label = \"'A'\" ]\n    0 -> 1 [ label = \"1\" ]\n    1 -> 2 [ label = \"1\" ]\
        \n    3 -> 4 [ label = \"1\" ]\n    4 -> 2 [ label = \"1\" ]\n}\n".to_string();
        assert_eq!(dot, output);
        // case 2
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"AAA");
        aligner.global(b"ABA").add_to_graph();
        let g = aligner.graph().map(|_, n| (*n) as char, |_, e| *e);
        let dot = format!("{:?}", Dot::new(&g));
        let output = "digraph {\n    0 [ label = \"'A'\" ]\n    1 [ label = \"'A'\" ]\n    2 [ label = \"'A'\" ]\
        \n    3 [ label = \"'B'\" ]\n    0 -> 1 [ label = \"1\" ]\n    1 -> 2 [ label = \"1\" ]\n    0 -> 3 [ label = \"1\" ]\
        \n    3 -> 2 [ label = \"1\" ]\n}\n".to_string();
        assert_eq!(dot, output);
        // case 3
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"BBBBBAAA");
        aligner.global(b"AAA").add_to_graph();
        let g = aligner.graph().map(|_, n| (*n) as char, |_, e| *e);
        let dot = format!("{:?}", Dot::new(&g));
        let output = "digraph {\n    0 [ label = \"'B'\" ]\n    1 [ label = \"'B'\" ]\n    2 [ label = \"'B'\" ]\
        \n    3 [ label = \"'B'\" ]\n    4 [ label = \"'B'\" ]\n    5 [ label = \"'A'\" ]\n    6 [ label = \"'A'\" ]\
        \n    7 [ label = \"'A'\" ]\n    0 -> 1 [ label = \"1\" ]\n    1 -> 2 [ label = \"1\" ]\n    2 -> 3 [ label = \"1\" ]\
        \n    3 -> 4 [ label = \"1\" ]\n    4 -> 5 [ label = \"1\" ]\n    5 -> 6 [ label = \"2\" ]\n    6 -> 7 [ label = \"2\" ]\n}\n".to_string();
        assert_eq!(dot, output);
        // case 4
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"AAA");
        aligner.global(b"BBBBBAAA").add_to_graph();
        let g = aligner.graph().map(|_, n| (*n) as char, |_, e| *e);
        let dot = format!("{:?}", Dot::new(&g));
        let output = "digraph {\n    0 [ label = \"'A'\" ]\n    1 [ label = \"'A'\" ]\n    2 [ label = \"'A'\" ]\
        \n    3 [ label = \"'B'\" ]\n    4 [ label = \"'B'\" ]\n    5 [ label = \"'B'\" ]\n    6 [ label = \"'B'\" ]\
        \n    7 [ label = \"'B'\" ]\n    0 -> 1 [ label = \"2\" ]\n    1 -> 2 [ label = \"2\" ]\n    3 -> 4 [ label = \"1\" ]\
        \n    4 -> 5 [ label = \"1\" ]\n    5 -> 6 [ label = \"1\" ]\n    6 -> 7 [ label = \"1\" ]\n    7 -> 0 [ label = \"1\" ]\n}\n".to_string();
        assert_eq!(dot, output);
    }

    #[test]
    fn test_consensus() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"GCATGCUx");
        aligner.global(b"GCATGCU").add_to_graph();
        aligner.global(b"xCATGCU").add_to_graph();
        assert_eq!(aligner.consensus(), b"GCATGCUx");
    }

    #[test]
    fn test_xclip_prefix_custom() {
        let x = b"GGGGGGATG";
        let y = b"ATG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(-5);

        let mut aligner = Aligner::new(scoring, x);
        let alignment = aligner.custom(y).alignment();

        assert_eq!(
            alignment.operations,
            [
                AlignmentOperation::Xclip(0),
                AlignmentOperation::Match(Some((5, 6))),
                AlignmentOperation::Match(Some((6, 7))),
                AlignmentOperation::Match(Some((7, 8)))
            ]
        );
    }

    #[test]
    fn test_yclip_prefix_custom() {
        let y = b"GGGGGGATG";
        let x = b"ATG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).yclip(-5);

        let mut aligner = Aligner::new(scoring, x);
        let alignment = aligner.custom(y).alignment();

        assert_eq!(
            alignment.operations,
            [
                AlignmentOperation::Yclip(0, 6),
                AlignmentOperation::Match(None),
                AlignmentOperation::Match(Some((0, 1))),
                AlignmentOperation::Match(Some((1, 2)))
            ]
        );
    }

    #[test]
    fn test_xclip_suffix_custom() {
        let x = b"GAAAA";
        let y = b"CG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(0).yclip(0);

        let mut aligner = Aligner::new(scoring, x);
        let alignment = aligner.custom(y).alignment();

        assert_eq!(
            alignment.operations,
            [
                AlignmentOperation::Yclip(0, 1),
                AlignmentOperation::Match(None),
                AlignmentOperation::Xclip(1)
            ]
        );
    }

    #[test]
    fn test_yclip_suffix_custom() {
        let y = b"GAAAA";
        let x = b"CG";

        let score = |a: u8, b: u8| if a == b { 3i32 } else { -3i32 };
        let scoring = Scoring::new(-5, -1, &score).yclip(-5).xclip(0);

        let mut aligner = Aligner::new(scoring, x);
        let alignment = aligner.custom(y).alignment();

        assert_eq!(
            alignment.operations,
            [
                AlignmentOperation::Xclip(0),
                AlignmentOperation::Match(Some((0, 1))),
                AlignmentOperation::Yclip(1, 5)
            ]
        );
    }
}
