# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

# HEAD
- Added HMM implementation.

# [0.21.0] - 2018-06-19
- Added PSSM implementation (thanks to @hervold).

# [0.20.0] - 2018-06-01
- Refactored RankSelect API to consistently use u64.
- Use bv crate in suffix array implementation.

# [0.19.0] - 2018-05-25
- rank-0 and select-0 in RankSelect.
- use bv crate for RankSelect.

## [0.18.0] - 2018-05-04
- More flexible FASTA API.
- Fixed bug in KMP.

## [0.17.0] - 2018-02-22
- Bug fix in Ukkonen algorithm
- Convenience improvements to API

## [0.16.0] - 2018-01-05
- Pairwise alignment has been rewritten to support banded alignment and clips.
- Various minor API additions and improvements.
- Several small bug fixes. 

## [0.15.0] - 2017-11-20
- Add pair hidden markov model implementation to calculate the probability of two sequences being related.
- Various minor bug fixes and usability improvements.

## [0.14.2] - 2017-08-30
- Improved numerical stability of CDF construction.
- Speed improvements to occurrence array lookups in FM-index.
- Improved GFF/GTF variant format handling.
- Improved robustness of credible interval calculating in CDF.
- Bug fixes for log probability implementation.

## [0.14.1] - 2017-06-23
### Changed
- Replace nalgebra dependency with ndarray crate.

## [0.14.0] - 2017-06-15
### Changed
- GTF/GFF reader can now handle duplicate keys.
- Updated dependencies.
- RNA alphabet.
- Improved FASTQ reader.
- Fixes in alignment algorithm.


## [0.13.0] - 2017-05-09
### Changed
- fasta::IndexedReader now also provides an iterator.
- IntervalTree provides a mutable iterator.
- Various fixes to Fasta IO.
- Fixed calculation of expected FDR.

## [0.12.0] - 2017-04-03
### Changed
- Improved distance API.
- Moved Strand into utils.
- More robust gff/gtf parsing.


## [0.11.0] - 2017-02-16
### Changed
- Improved IntervalTree API.
- Updated dependencies.
- Speed improvements in alignment module.
- Improved test coverage.
- Speed improvements in fmindex module.

## [0.10.0] - 2016-11-02
### Added
- An interval tree implementation.
- Initial utilities for bayesian statistics.
### Changed
- Various small improvements to log-space probability API.

## [0.9.0] - 2016-08-18
### Added
- Implementation of discrete probability distributions via cumulative distribution functions.
### Changed
- Log-space probabilities have been refactored into newtypes.
- Performance improvements for FMIndex implementation.
- Improved documentation.

## [0.8.0] - 2016-07-20
### Changed
- Writers in the io module no longer take ownership of the given record.
- Various cosmetic changes.

## [0.7.0] - 2016-07-06
### Changed
- Reverse complement API has been refactored into plain functions.
- Reverse complement now supports the whole IUPAC alphabet.
- Various algorithms take now IntoTextIterator instead of only slices.
- Fasta reader and writer treat sequence names as strings.
- Refactoring of suffix array + fmindex API to provide more flexibility.

## [0.6.0] - 2016-05-09
### Changed
- Type aliases for various text representations.
- Pattern matching algorithms take both iterators and slices where possible.
- logprobs::cumsum has been refactored to return an iterator.
- support for subtraction of logprobs.

## [0.5.0] - 2016-02-24
### Added
- Support for [serde](https://github.com/serde-rs/serde) serialization when used in combination with rust nightly (@dikaiosune).
