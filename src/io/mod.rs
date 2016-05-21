//! Readers and writers for common bioinformatics file formats.


pub mod fastq;
pub mod fasta;
pub mod bed;
pub mod gff;

/// Strand information.
#[derive(Debug, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}
