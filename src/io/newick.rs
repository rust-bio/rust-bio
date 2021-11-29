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
//!  for taxon in tree.raw_nodes() {
//!      println!("{}", taxon.weight);
//!  }
//!  ```

use bio_types::phylogeny::{Tree, TreeGraph};
use pest::iterators::Pair;
use pest::Parser;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use petgraph::EdgeDirection::Incoming;
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

    #[error("Error while writing tree: {0}")]
    Write(std::io::Error),

    #[error("Tree contains invalid UTF-8: {0}")]
    InvalidContent(#[from] std::str::Utf8Error),

    #[error("Error while parsing tree: {0}")]
    ParsingError(#[from] pest::error::Error<crate::io::newick::Rule>),

    #[error("No unique root")]
    NoUniqueRoot,
}
type Result<T, E = Error> = std::result::Result<T, E>;

/// The parser is automagically derived from the `newick.pest` grammar
/// file
#[derive(Parser)]
#[grammar = "io/newick.pest"]
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
                let name = if let Some(clade) = inner_rules.next() {
                    Some(clade.as_str().into())
                } else {
                    None
                };
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

    let root = NewickParser::parse(Rule::Tree, &content)
        .map_err(Error::ParsingError)?
        .next()
        .unwrap();

    Ok(parse_value(root))
}

/// Convert an intermediary `TreeValue` to the public `Tree` type
fn newick_to_graph(root: TreeValue) -> Result<Tree> {
    fn add_node(g: &mut TreeGraph, t: TreeValue) -> NodeIndex {
        match t {
            TreeValue::Node { name, children } => {
                let node_id = g.add_node(name.unwrap_or("".into()).into());
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
pub fn read<R: io::Read>(reader: R) -> Result<Tree> {
    let content_bytes = reader
        .bytes()
        .collect::<Result<Vec<_>, _>>()
        .map_err(Error::Read)?;
    let content_str = std::str::from_utf8(&content_bytes).map_err(Error::InvalidContent)?;
    from_string(&content_str)
}

/// Convert the Tree to the Newick String representation.
pub fn to_string(t: &Tree) -> Result<String> {
    fn node_to_string(i: NodeIndex, g: &TreeGraph, s: &mut String) {
        let edges = g.edges_directed(i, petgraph::Outgoing).collect::<Vec<_>>();
        if !edges.is_empty() {
            *s += "(";
            let mut first = true;
            // Petgraph iterates over edges in reverse order of adding them.
            // Reversing here ensures that `to_string` is compatible with `from_string`.
            for edge in edges.iter().rev() {
                if first {
                    first = false;
                } else {
                    *s += ",";
                }
                node_to_string(edge.target(), g, s);
                if !edge.weight().is_nan() {
                    *s += ":";
                    *s += &edge.weight().to_string();
                }
            }
            *s += ")";
        }
        // `i` should never be an invalid index here.
        *s += g.node_weight(i).unwrap();
    }

    // Find the root of the tree.
    let roots = t.g.externals(Incoming);
    // Make sure there is exactly one vertex without incoming edges.
    let root = roots.next().ok_or(Error::NoUniqueRoot)?;
    if roots.next().is_some() {
        return Err(Error::NoUniqueRoot);
    }

    let mut s = String::new();
    node_to_string(root, &t.g, &mut s);
    Ok(s + ";")
}

/// Writes a tree to a file.
pub fn to_file<P: AsRef<Path>>(path: P, tree: &Tree) -> Result<()> {
    fs::File::open(&path)
        .map(|w| write(w, tree))
        .map_err(|e| Error::OpenFile {
            filename: path.as_ref().to_owned(),
            source: e,
        })?
}

/// Writes a tree to any type implementing `io::Write`.
pub fn write<W: io::Write>(mut writer: W, tree: &Tree) -> Result<()> {
    let s = to_string(&tree)?;
    writer.write_all(&s.into_bytes()).map_err(Error::Write)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test that to_string and from_string are inverses.
    #[test]
    fn newick_from_to_string() {
        // Test input from Wikipedia: https://en.wikipedia.org/wiki/Newick_format
        let strings = vec![
            // This is not supported currently:
            // The grammar allows it, but the value is ignored.
            //"(:0.1,:0.2,(:0.3,:0.4):0.5):1.0;", /* all have a distance to parent, including the root. */
            "A;",
            "(A,B);",
            "(,,(,));",                            //no nodes are named
            "(A,B,(C,D));",                        //leaf nodes are named
            "(A,B,(C,D)E)F;",                      //all nodes are named
            "(:0.1,:0.2,(:0.3,:0.4):0.5);",        //all but root node have a distance to parent
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",    //distances and leaf names (popular)
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",  //distances and all names
            "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", //a tree rooted on a leaf node (rare)
        ];
        for s in strings {
            let t = from_string(s).unwrap();
            assert_eq!(to_string(&t).unwrap(), s);
        }
    }
}
