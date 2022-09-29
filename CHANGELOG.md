# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).


## [1.0.0](https://www.github.com/rust-bio/rust-bio/compare/v0.42.0...v1.0.0) (2022-09-29)


### ⚠ BREAKING CHANGES

* update to latest ordered-float (#503)

### Features

* detect nested ORFs ([#501](https://www.github.com/rust-bio/rust-bio/issues/501)) ([19b6c36](https://www.github.com/rust-bio/rust-bio/commit/19b6c36bbea910eaaa0753a37e24b1389f31f527))


### Bug Fixes

* update to latest ordered-float ([#503](https://www.github.com/rust-bio/rust-bio/issues/503)) ([63fb752](https://www.github.com/rust-bio/rust-bio/commit/63fb752616bfb02791f10779cd37bc072487923e))

## [0.42.0](https://www.github.com/rust-bio/rust-bio/compare/v0.41.0...v0.42.0) (2022-08-30)


### ⚠ BREAKING CHANGES

* Update `strum` and `ordered-float` dependencies and change From<LogProb> into TryFrom<LogProb> for NotNan<f64>. (#491)

### Features

* Update `strum` and `ordered-float` dependencies and change From<LogProb> into TryFrom<LogProb> for NotNan<f64>. ([#491](https://www.github.com/rust-bio/rust-bio/issues/491)) ([57ccf8f](https://www.github.com/rust-bio/rust-bio/commit/57ccf8ff716416f7dbcba7f42a5e4369cea2fea0))


### Miscellaneous Chores

* widen range on statrs ([#499](https://www.github.com/rust-bio/rust-bio/issues/499)) ([0ff1e70](https://www.github.com/rust-bio/rust-bio/commit/0ff1e70cf5fee09f227c89316d62653c2414a8a0))

## [0.41.0](https://www.github.com/rust-bio/rust-bio/compare/v0.40.0...v0.41.0) (2022-03-30)


### Features

* adaptive integration of density functions using a binary search approach that tries to achieve good resolution around the maximum ([#486](https://www.github.com/rust-bio/rust-bio/issues/486)) ([207b76f](https://www.github.com/rust-bio/rust-bio/commit/207b76fa9bccce4236e3cda9c10e56be7a636a61))

## [0.40.0](https://www.github.com/rust-bio/rust-bio/compare/v0.39.2...v0.40.0) (2022-02-25)

### Features

* base specific hop parameters in homopoly-pair-hmm ([#480](https://www.github.com/rust-bio/rust-bio/issues/480)) ([cf75b6c](https://www.github.com/rust-bio/rust-bio/commit/cf75b6cb5280dde52b26f95f8ec9cd37706a642d))


### Miscellaneous Chores

* release 0.40.0 ([cf8ebc3](https://www.github.com/rust-bio/rust-bio/commit/cf8ebc36f0cef5cf3900a3aebe7bca7bd2f30e78))

### [0.39.2](https://www.github.com/rust-bio/rust-bio/compare/v0.39.1...v0.39.2) (2022-02-09)


### Bug Fixes

* Make QGramIndex use less memory; fix bug; improve tests&docs ([#471](https://www.github.com/rust-bio/rust-bio/issues/471)) ([48bac1c](https://www.github.com/rust-bio/rust-bio/commit/48bac1cc9efd236b8b904997b40af1f64cd5f255))
* overflow in `qgrams` for k=32 ([#478](https://www.github.com/rust-bio/rust-bio/issues/478)) ([8048eb8](https://www.github.com/rust-bio/rust-bio/commit/8048eb8ced2087659184f70bae9e0d91690aa212))

### [0.39.1](https://www.github.com/rust-bio/rust-bio/compare/v0.39.0...v0.39.1) (2022-01-12)


### Bug Fixes

* added code to ignore commented lines in a bed file  ([#474](https://www.github.com/rust-bio/rust-bio/issues/474)) ([d17f823](https://www.github.com/rust-bio/rust-bio/commit/d17f823de1466c0fa2f21ae1dbdd1298a36744e6))

## [0.39.0](https://www.github.com/rust-bio/rust-bio/compare/v0.38.0...v0.39.0) (2021-10-20)


### Features

* Backward search api ([#457](https://www.github.com/rust-bio/rust-bio/issues/457)) ([eb9d378](https://www.github.com/rust-bio/rust-bio/commit/eb9d378c70bc43d1dcf477e48cb167802799a546))

## [0.38.0](https://www.github.com/rust-bio/rust-bio/compare/v0.37.1...v0.38.0) (2021-10-04)


### Features

* adds iupac amino acid alpha ([#444](https://www.github.com/rust-bio/rust-bio/issues/444)) ([f4b7a7c](https://www.github.com/rust-bio/rust-bio/commit/f4b7a7c23aa7e165718c77e91f62f60fc431ba2c))
* Allow to compute bayesian model via the exploration of the marginal distribution. ([#453](https://www.github.com/rust-bio/rust-bio/issues/453)) ([d09033c](https://www.github.com/rust-bio/rust-bio/commit/d09033c1f42da89514ec87bf4aba5d1614e9810d))


### Bug Fixes

* backward search yielding potentially incorrect positions on FM-Index ([#454](https://www.github.com/rust-bio/rust-bio/issues/454)) ([#455](https://www.github.com/rust-bio/rust-bio/issues/455)) ([3489e6a](https://www.github.com/rust-bio/rust-bio/commit/3489e6ab493b5191b879b4405b868a4004bbb88e))

### [0.37.1](https://www.github.com/rust-bio/rust-bio/compare/v0.37.0...v0.37.1) (2021-08-23)


### Bug Fixes

* One Simple Trick to significantly improve the speed and memory usage of Occ ([#448](https://www.github.com/rust-bio/rust-bio/issues/448)) ([9aa79cb](https://www.github.com/rust-bio/rust-bio/commit/9aa79cbb76960af47178bf0c06d8507872fbe3b0))
* sampled suffix array ([#447](https://www.github.com/rust-bio/rust-bio/issues/447)) ([00f9846](https://www.github.com/rust-bio/rust-bio/commit/00f9846ba5cb717b3e5392301029f9c46ecea527))

## [0.37.0](https://www.github.com/rust-bio/rust-bio/compare/v0.36.0...v0.37.0) (2021-07-09)


### Features

* add method to retrieve all event posteriors from Bayesian model. ([e5feda8](https://www.github.com/rust-bio/rust-bio/commit/e5feda8ad4101cb992ca4ad82f4fc96e0e1574ba))

## [0.36.0](https://www.github.com/rust-bio/rust-bio/compare/v0.35.0...v0.36.0) (2021-07-06)


### Features

* Baum-Welch algorithm for Discrete HMM ([#432](https://www.github.com/rust-bio/rust-bio/issues/432)) ([eb8b8cb](https://www.github.com/rust-bio/rust-bio/commit/eb8b8cbad0016b0ab91861cb8d33f7fb624fb157))


## [0.35.0] - 2021-07-05
- Improved buffer control in Fasta and Fastq API (@natir).
- Fixed an indexing bug in ArrayBackedIntervalTree (@wabain).
- Fixed a corner case where the FASTX parser could have looped infinitely (@morsecodist).
- Fixed compiler warnings (@fxwiegand).
- Improved documentation for FASTA index (@mbhall88).


## [0.34.0] - 2021-05-04
- Bayesian model framework now relies on Hash instead of Ord for accessing events (@johanneskoester).
- Added wavelet matrix datastructure (@Identi, @tedil).


## [0.33.0] - 2021-03-09
- Fixed a floating point error in gcn_content (@tedil).
- Improved error messages in io module (@fxwiegand).
- Better memory usage of HomopolypairHMM (@tedil).
- Improved documentation (@TianShi2001, @m0ssc0de, @dcroote).
- Support for reading compact BED files (@manzt).
- API improvements for `Alphabet` and GFF reader (@tshauck).
- Switched to thiserror for error handling (@delehef).
- Removed unsafe code (@huonw).
- Fixed overflow in ShiftAnd algorithm (@dcroote).
- Various additional test cases (@dcroote).
- Extended API for SMEM computation on FMDIndex (@Identi).
- Added a parser for Newick phylogenetic trees (@delehef).


## [0.32.0] - 2020-07-28
This release mostly comprises of the documentation improvements made
by the recent Docathon. Lots of new doctests and a bunch of time and
memory complexity annotations. Big thanks to @tedil, @dcroote,
@TomKellyGenetics, @natir, @thomasmulvaney, @mbhall88, @jafors,
@HenningTimm, @luizirber, @dlaehnemann, @johanneskoester.
Further additions in this release are:
- homopolymer-error-aware pairHMM implementation (thanks to @tedil).
- SIMD-accelerated edit distance routines (thanks to @Daniel-Liu-c0deb0t).
- BitEnc derives more traits (thanks to @FelixMoelder).
- wider pinning of some dependencies' version numbers (thanks to @pmarks).


## [0.31.0] - 2020-06-02
- Bugfix for pHMM implementation (thanks to @tedil).
- Sorted array-backed interval trees (thanks to @tedil).


### [0.30.1] - 2020-05-13
- Improved occ counting speed for FM index (thanks to @thomasmulvaney)
- Various small bug fixes and code cleanups and linter fixes.

## [0.30.0] - 2019-11-14
- Bayesian models now allow to access internals.
- Various small bug fixes.

## [0.29.0] - 2019-09-27
- Migrate error handling to the snafu crate (this is an API breaking change).
- Fix edge cases in pairwise alignment.
- Fix error in backward search if symbol isn't found.

### [0.28.1] - 2019-06-28
- Fix select in RankSelect in cases where many superblocks have the same rank.

## [0.28.0] - 2019-06-19
- Myers bit-parallel pattern matching now supports arbitrarily long patterns via bit vectors (thanks to @markschl).
- Minor documentation updates (thanks to @anders-was-here).

## [0.27.0] - 2019-05-31
- Implement sequence-read-trait for FASTQ records.
- Cleanup dependencies.

### [0.26.1] - 2019-05-10
- Fix a bug in `select_1` and `select_0` that would lead to too large answers.

## [0.26.0] - 2019-05-09
- Added a trait system for computing Bayesian statistical models.
- Added an implementation of MSA via partial order alignment.
- Performance improvements to FASTQ reader.

## [0.25.0] - 2018-12-12
- Added `FQRead` and `FARead` traits to `FastaReader` and `FastqReader` to be more flexible with input types. This allows to use readers on gzipped and on plain text input interchangeably.
- Added an implementation of Bayes Factors and evidence scoring using the method of Kass and Raftery.

## [0.24.0] - 2018-11-26
- API overhaul to become more flexible when accepting text iterators. Now, anything that iterates over something can be borrowed as u8 is allowed.
- FMIndex and FMDIndex now also allow plain owned versions of BWT, Less and Occ. This should greatly simplify their usage.
- PairHMM and LogProb implementation has seen extensive performance improvements. Among that, (a) the usage of a fast approximation of exp() as presented by [Kopczynsi 2017](https://eldorado.tu-dortmund.de/bitstream/2003/36203/1/Dissertation_Kopczynski.pdf), and (b) banding of the pairHMM matrix with a given maximum edit distance.
- All IO records now support serde.

## [0.23.0] - 2018-11-06
- Generalized Myers pattern matching algorithm to arbitrary unsigned integer types (u64, u128) (thanks to @markschl).
- Implemented optional traceback and alignment output for Myers pattern matching algorithm (thanks to @markschl).
- Use Strand type from bio-types crate in BED module (thanks to @ingolia).
- Added an IntervalTree based data structure for looking up overlaps between annotation types (thanks to @ingolia).
- Various bug fixes.

## [0.22.0] - 2018-08-01
- Added HMM implementation (thanks to @holtgrewe).
- Moved Alignment types to `bio_types` crate (thanks to @pmarks).
- Ignore comment lines in GTF/GFF files (thanks to Yasunobu Okamura).
- API usability improvements.

## [0.21.0] - 2018-06-19
- Added PSSM implementation (thanks to @hervold).

## [0.20.0] - 2018-06-01
- Refactored RankSelect API to consistently use u64.
- Use bv crate in suffix array implementation.

## [0.19.0] - 2018-05-25
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
