//! Readers and writers for common bioinformatics file formats.

pub mod bed;
pub mod bedpe;
mod core_record;
pub mod fasta;
pub mod fastq;
pub mod fastx;
pub mod gff;
#[cfg(feature = "phylogeny")]
pub mod newick;
use core_record::Writer;
