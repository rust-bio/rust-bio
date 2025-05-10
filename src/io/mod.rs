//! Readers and writers for common bioinformatics file formats.

pub mod bed;
pub mod bedpe;
pub mod core;
pub mod fasta;
pub mod fastq;
pub mod fastx;
pub mod gff;
#[cfg(feature = "phylogeny")]
pub mod newick;
