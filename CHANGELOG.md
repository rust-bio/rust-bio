# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
- Type aliases for various text representations.
- Pattern matching algorithms take both iterators and slices where possible.
- logprobs::cumsum has been refactored to return an iterator.

## [0.5.0] - 2016-02-24
### Added
- Support for [serde](https://github.com/serde-rs/serde) serialization when used in combination with rust nightly (@dikaiosune).
