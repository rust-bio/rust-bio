//! Readers and writers for common bioinformatics file formats.

#[path = "bed/bed.rs"]
pub mod bed;

#[path = "bed/common.rs"]
pub mod common;

pub mod bedpe;
pub mod fasta;
pub mod fastq;
pub mod fastx;
pub mod gff;
#[cfg(feature = "phylogeny")]
pub mod newick;
