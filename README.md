[![Build Status](https://travis-ci.org/rust-bio/rust-bio.svg?branch=master)](https://travis-ci.org/rust-bio/rust-bio)
[![](http://meritbadge.herokuapp.com/bio)](https://crates.io/crates/bio)

# Rust-Bio, a bioinformatics library for Rust.

This library provides implementations of many algorithms and data structures
that are useful for bioinformatics.
All provided implementations are rigorously tested via continuous
integration.

Currently, rust-bio provides

* most major pattern matching algorithms,
* a convenient alphabet implementation,
* pairwise alignment,
* suffix arrays,
* BWT and FM-Index,
* FMD-Index for finding supermaximal exact matches,
* a q-gram index,
* a rank/select data structure,
* FASTQ and FASTA and BED readers and writers,
* helper functions for combinatorics and dealing with log probabilities.

For reading and writing BAM and BCF files, have a look at https://github.com/christopher-schroeder/rust-htslib.

Please see the [homepage](https://rust-bio.github.io) for examples and documentation.

## Author

Johannes Köster (<koester@jimmy.harvard.edu>)

## Contributors

* Christopher Schröder (<christopher.schroeder@uni-due.de>)
* Peer Aramillo Irizar
* Fedor Gusev

The next name in this list could be you! If you are interested in joining the effort to build a general purpose Rust bioinformatics library, just introduce yourself [here](https://github.com/rust-bio/rust-bio/issues/3), or issue a pull request with your first contribution.

## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
