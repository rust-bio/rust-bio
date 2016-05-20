[![Crates.io](https://img.shields.io/crates/d/bio.svg?style=flat-square)](https://crates.io/crates/bio)
[![Crates.io](https://img.shields.io/crates/v/bio.svg?style=flat-square)](https://crates.io/crates/bio)
[![Crates.io](https://img.shields.io/crates/l/bio.svg?style=flat-square)](https://crates.io/crates/bio)
[![Travis](https://img.shields.io/travis/rust-bio/rust-bio/master.svg?style=flat-square)](https://travis-ci.org/rust-bio/rust-bio)

# Rust-Bio, a bioinformatics library for Rust.

This library provides implementations of many algorithms and data structures
that are useful for bioinformatics.
All provided implementations are rigorously tested via continuous
integration.

**Please see the [homepage](https://rust-bio.github.io) for examples and documentation.**

Currently, rust-bio provides

* most major pattern matching algorithms,
* a convenient alphabet implementation,
* pairwise alignment,
* suffix arrays,
* BWT and FM-Index,
* FMD-Index for finding supermaximal exact matches,
* a q-gram index,
* a rank/select data structure,
* [serde](https://github.com/serde-rs/serde) support for all data structures when built with `nightly` feature,
* FASTQ and FASTA and BED readers and writers,
* helper functions for combinatorics and dealing with log probabilities.

For reading and writing BAM and BCF files, have a look at https://github.com/christopher-schroeder/rust-htslib.

## Author

[Johannes Köster](https://github.com/johanneskoester)

## Contributors

* [Christopher Schröder](https://github.com/christopher-schroeder)
* [Peer Aramillo Irizar](https://github.com/parir)
* [Fedor Gusev](https://github.com/gusevfe)
* [Vadim Nazarov](https://github.com/vadimnazarov)
* [Brad Chapman](https://github.com/chapmanb)
* [Florian Gilcher](https://github.com/skade)
* [Erik Clarke](https://github.com/eclarke)
* [Rizky Luthfianto](https://github.com/rilut)
* [Adam Perry](https://github.com/dikaiosune)
* [Taylor Cramer](https://github.com/cramertj)
* [Andre Bogus](https://github.com/llogiq)
* Philipp Angerer

The next name in this list could be you! If you are interested in joining the effort to build a general purpose Rust bioinformatics library, just introduce yourself [here](https://github.com/rust-bio/rust-bio/issues/3), or issue a pull request with your first contribution.

## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
