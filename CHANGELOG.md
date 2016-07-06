# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.7.0] - 2016-07-06
- Reverse complement API has been refactored into plain functions.
- Reverse complement now supports the whole IUPAC alphabet.
- Various algorithms take now IntoTextIterator instead of only slices.
- Fasta reader and writer treat sequence names as strings.
- Refactoring of suffix array + fmindex API to provide more flexibility.

## [0.6.0] - 2016-05-09
- Type aliases for various text representations.
- Pattern matching algorithms take both iterators and slices where possible.
- logprobs::cumsum has been refactored to return an iterator.
- support for subtraction of logprobs.

## [0.5.0] - 2016-02-24
### Added
- Support for [serde](https://github.com/serde-rs/serde) serialization when used in combination with rust nightly (@dikaiosune).
