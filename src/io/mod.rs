//! Readers and writers for common bioinformatics file formats.

pub mod bed;
pub mod bnx;
pub mod cmap;
pub mod fasta;
pub mod fastq;
pub mod gff;
#[cfg(feature = "phylogeny")]
pub mod newick;
pub mod om_utils;
pub mod xmap;
