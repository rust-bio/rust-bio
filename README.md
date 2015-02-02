# Rust-Bio, a rusty bioinformatics library.

This library provides implementations of many algorithms and data structures
that are useful for bioinformatics.
All provided implementations are rigorously tested via continuous
integration (https://travis-ci.org/johanneskoester/rust-bio).

Currently, rust-bio provides

* pattern matching,
* pairwise alignment,
* suffix arrays,
* BWT and FM-Index,
* rank/select data structures.

Find out more at https://github.com/johanneskoester/rust-bio.

## Usage

To use rust-bio in your Rust project, add the following to your `Cargo.toml`

```toml
[dependencies]
bio = "*"
```

and import the crate from your source code:

```rust
extern crate bio;
// use e.g. a pattern matching algorithm
use bio::pattern_matching::bndm::BNDM;

let pattern = b"GAAAA";
let text = b"ACGGCTAGAAAAGGCTAGAAAA";
let bndm = bndm::BNDM::new(pattern);
let let matches = bndm.find_all(text);
```

For more information, please read the API documentation: https://johanneskoester.github.io/rust-bio/doc/bio/index.html

## Authors 

* Johannes Köster <johannes.koester@tu-dortmund.de>

## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
