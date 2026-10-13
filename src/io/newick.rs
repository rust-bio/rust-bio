// Copyright 2020 Franklin Delehelle.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A struct to read phylogenetic trees in the Newick format.
//!
//!  # Example
//!
//!  In this example, we parse a tree from a string and display all the taxons.
//!  See `petgraph` documentation for more details on how to handle the tree.
//!
//!  ```
//!  use bio::io::newick;
//!
//!  let tree = newick::from_string("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
//!  for taxon in tree.g.raw_nodes() {
//!      println!("{}", taxon.weight);
//!  }
//!  ```

use bio_types::phylogeny::{Tree, TreeGraph};
use pest::iterators::Pair;
use pest::Parser;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use thiserror::Error;

/// A `thiserror` error type gathering all the potential bad outcomes
#[derive(Debug, Error)]
pub enum Error {
    #[error("Error while opening {}: {}", filename.display(), source)]
    OpenFile {
        filename: PathBuf,
        source: std::io::Error,
    },

    #[error("Error while reading tree: {0}")]
    Read(#[from] std::io::Error),

    #[error("Tree contains invalid UTF-8: {0}")]
    InvalidContent(#[from] std::str::Utf8Error),

    #[error("Error while parsing tree: {0}")]
    ParsingError(#[from] pest::error::Error<crate::io::newick::Rule>),

    #[error("Tree has no root node")]
    EmptyTree,
}
type Result<T, E = Error> = std::result::Result<T, E>;

/// The parser is automagically derived from the `newick.pest` grammar
/// file
#[derive(Parser)]
#[grammar = "io/newick.pest"]
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct NewickParser;

/// A hidden, temporary datatype used to collect the parser result
/// before converting it to a `Tree`
enum TreeValue {
    Node {
        name: Option<String>,
        children: Option<Vec<TreeValue>>,
    },
    Link {
        weight: f32,
        node: Box<TreeValue>,
    },
}

/// Given a string representing a Newick tree, tries to parse it and
/// returns the `TreeValue` of the root
fn parse_newick_file(content: &str) -> Result<TreeValue> {
    fn parse_value(pair: Pair<Rule>) -> TreeValue {
        match pair.as_rule() {
            Rule::Leaf => {
                let name = pair.into_inner().next().unwrap().as_str();
                TreeValue::Node {
                    name: Some(name.into()),
                    children: None,
                }
            }
            Rule::Internal => {
                let mut inner_rules = pair.into_inner();
                let children = Some(
                    inner_rules
                        .next()
                        .unwrap()
                        .into_inner()
                        .map(parse_value)
                        .collect(),
                );
                let name = inner_rules.next().map(|clade| clade.as_str().into());
                TreeValue::Node { children, name }
            }

            Rule::Branch => {
                fn get_weight(mut inner: pest::iterators::Pairs<Rule>) -> f32 {
                    if let Some(weight) = inner.next() {
                        weight.as_str().parse::<f32>().unwrap()
                    } else {
                        f32::NAN
                    }
                }

                let mut inner = pair.into_inner();
                let (node, weight) = if let Some(next) = inner.next() {
                    match next.as_rule() {
                        Rule::SubTree => (parse_value(next), get_weight(inner)),
                        _ => (
                            TreeValue::Node {
                                name: None,
                                children: None,
                            },
                            next.as_str().parse::<f32>().unwrap(),
                        ),
                    }
                } else {
                    (
                        TreeValue::Node {
                            name: None,
                            children: None,
                        },
                        get_weight(inner),
                    )
                };

                TreeValue::Link {
                    weight,
                    node: Box::new(node),
                }
            }

            Rule::SubTree => parse_value(pair.into_inner().next().unwrap()),
            Rule::EOI
            | Rule::WHITESPACE
            | Rule::Tree
            | Rule::Length
            | Rule::BranchSet
            | Rule::float
            | Rule::safe
            | Rule::name => unreachable!(),
        }
    }

    let root = NewickParser::parse(Rule::Tree, content)
        .map_err(Error::ParsingError)?
        .next()
        .unwrap();

    Ok(parse_value(root))
}

/// Convert an intermediary `TreeValue` to the public `Tree` type
fn newick_to_graph(root: TreeValue) -> Result<Tree> {
    // The node index type belongs to `bio_types`' petgraph, which may be a
    // different (semver-incompatible) version from this crate's own petgraph;
    // naming it here would pick the wrong one, so it is left to inference.
    // `use<>` keeps the opaque type from capturing the borrow of `g` (the
    // default from edition 2024 on), which the recursion below needs.
    fn add_node(g: &mut TreeGraph, t: TreeValue) -> impl Copy + use<> {
        match t {
            TreeValue::Node { name, children } => {
                let node_id = g.add_node(name.unwrap_or_else(|| "N/A".into()));
                if let Some(children) = children {
                    for child in children {
                        match child {
                            TreeValue::Node { .. } => unimplemented!(),
                            TreeValue::Link { weight, node } => {
                                let child_id = add_node(g, *node);
                                g.add_edge(node_id, child_id, weight);
                            }
                        }
                    }
                };
                node_id
            }
            TreeValue::Link { .. } => unreachable!(),
        }
    }

    // The grammar accepts a bare `;` or `:0.5;` as a root branch without a
    // subtree: there is no node to build, and `add_node` would panic on it.
    if let TreeValue::Link { .. } = root {
        return Err(Error::EmptyTree);
    }

    let mut g = TreeGraph::new();
    add_node(&mut g, root);

    Ok(Tree { g })
}

/// Reads a tree from an `&str`-compatible type
pub fn from_string<S: AsRef<str>>(content: S) -> Result<Tree> {
    let raw_tree = parse_newick_file(content.as_ref())?;
    newick_to_graph(raw_tree)
}

/// Reads a tree from a file
pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Tree> {
    fs::File::open(&path)
        .map(read)
        .map_err(|e| Error::OpenFile {
            filename: path.as_ref().to_owned(),
            source: e,
        })?
}

/// Reads a tree from any type implementing `io::Read`
pub fn read<R: io::Read>(mut reader: R) -> Result<Tree> {
    let mut content_bytes = Vec::new();
    reader
        .read_to_end(&mut content_bytes)
        .map_err(Error::Read)?;
    let content_str = std::str::from_utf8(&content_bytes).map_err(Error::InvalidContent)?;
    from_string(content_str)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(tree: &Tree) -> Vec<String> {
        tree.g
            .raw_nodes()
            .iter()
            .map(|n| n.weight.clone())
            .collect()
    }

    #[test]
    fn test_parse_nested_tree_with_lengths() {
        let tree = from_string("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
        assert_eq!(names(&tree), vec!["F", "A", "B", "E", "C", "D"]);
        assert_eq!(tree.g.edge_count(), 5);
        let weights: Vec<f32> = tree.g.raw_edges().iter().map(|e| e.weight).collect();
        assert_eq!(weights, vec![0.1, 0.2, 0.3, 0.4, 0.5]);
    }

    #[test]
    fn test_missing_lengths_are_nan() {
        let tree = from_string("(A,B);").unwrap();
        assert_eq!(names(&tree), vec!["N/A", "A", "B"]);
        assert!(tree.g.raw_edges().iter().all(|e| e.weight.is_nan()));
    }

    #[test]
    fn test_single_node_trees_are_valid() {
        let tree = from_string("A;").unwrap();
        assert_eq!(names(&tree), vec!["A"]);
        assert_eq!(tree.g.edge_count(), 0);
        let tree = from_string("A:0.5;").unwrap();
        assert_eq!(names(&tree), vec!["A"]);
    }

    #[test]
    fn test_empty_root_is_an_error_not_a_panic() {
        for input in [";", ":0.5;", " ;"] {
            match from_string(input) {
                Err(Error::EmptyTree) => {}
                other => panic!(
                    "{:?}: expected EmptyTree, got {:?}",
                    input,
                    other.map(|t| names(&t))
                ),
            }
        }
    }

    #[test]
    fn test_syntax_error_is_reported() {
        assert!(matches!(from_string("(A,B"), Err(Error::ParsingError(_))));
        assert!(matches!(from_string("(A,B)"), Err(Error::ParsingError(_))));
    }

    #[test]
    fn test_read_from_reader() {
        let tree = read("(A:1,B:2)C;".as_bytes()).unwrap();
        assert_eq!(names(&tree), vec!["C", "A", "B"]);
    }
}
