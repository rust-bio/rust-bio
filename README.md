[![Build Status](https://travis-ci.org/johanneskoester/rust-bio.svg?branch=master)](https://travis-ci.org/johanneskoester/rust-bio)
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

## Resources

* Homepage: https://github.com/johanneskoester/rust-bio
* Documentation: https://johanneskoester.github.io/rust-bio
* Continuous integration tests: https://travis-ci.org/johanneskoester/rust-bio
* Roadmap: https://github.com/johanneskoester/rust-bio/issues/3

## Usage

We explain how to use Rust-Bio step-by-step. Users who already have experience with Rust can skip the first two steps.

### Step 1: Setting up Rust

Currently, Rust-Bio needs the development version (nightly) of Rust to compile properly, since it depends on features that are not in the stable branch yet. 
Rust nightly can be installed or updated by executing

```bash
curl -sSf https://static.rust-lang.org/rustup.sh | sh -s -- --channel=nightly
```
in your terminal.
For details or alternative ways of installation, have a look at the "Nightly" section of the Rust [download page](https://www.rust-lang.org/install.html).

### Step 2: Setting up a new Rust project

Since Rust-Bio is a library, you need to setup your own new Rust project to use Rust-Bio.
With Rust, projects and their dependencies are managed with the builtin package manager [Cargo](https://crates.io/).
To create a new Rust project, issue

```bash
cargo new hello_world --bin
cd hello_world
```
in your terminal. The flag `--bin` tells Cargo to create an executable project instead of a library.
In [this section](http://doc.rust-lang.org/nightly/book/hello-cargo.html#a-new-project) of the Rust docs, you find details about what Cargo just created for you.

Your new project can be compiled with
```bash
cargo build
```
and executed with
```bash
cargo run
```
If you are new to Rust, we suggest to proceed with [learning Rust](http://doc.rust-lang.org/nightly/book/learn-rust.html) via the Rust docs.

### Step 3: Use Rust-Bio from your project

To use Rust-Bio in your Rust project, add the following to your `Cargo.toml`

```toml
[dependencies]
bio = "*"
```

and import the crate from your source code:

```rust
extern crate bio;
```

### Example

An example usage of Rust-Bio is presented in the following:
```rust
// Import some modules
use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::bwt;
use bio::data_structures::fmindex::FMIndex;
use bio::io::fastq;

// Create an FM-Index for a given text.
let alphabet = alphabets::dna::iupac_alphabet();
let pos = suffix_array(text);
let bwt = bwt(text, &pos);
let fmindex = FMIndex::new(&bwt, 3, &alphabet);


// Iterate over a FASTQ file, use the alphabet to validate read
// sequences and search for exact matches in the FM-Index.
let reader = fastq::Reader::from_file("reads.fastq");
for record in reader.records() {
    let seq = record.seq();
    if alphabet.is_word(seq) {
        let interval = fmindex.backward_search(seq.iter());
        let positions = interval.occ(&pos);
    }
}
```
For more information and additional examples, please read the documentation: https://johanneskoester.github.io/rust-bio

## Benchmarks

Since Rust-Bio is based on a compiled language, similar performance to C/C++ based libraries can be expected. Indeed, we find the pattern matching algorithms of Rust-Bio to perform in the range of the C++ library Seqan:

| Algorithm | Rust-Bio | Seqan   |
| --------- | -------: | ------: |
| BNDM      | 77ms     | 82ms    |
| Horspool  | 122ms    | 127ms   |
| BOM       | 103ms    | 109ms   |
| Shift-And | 241ms    | 547ms   |

We measured 10000 iterations of searching pattern `GCGCGTACACACCGCCCG` in the sequence of the hg38 MT chromosome.
Initialization time of each algorithm for the given pattern was included in each iteration. Benchmarks were conducted with *Cargo bench* for Rust-Bio and *Python timeit* for Seqan on an Intel Core i5-3427U CPU.
Benchmarking Seqan from *Python timeit* entails an overhead of 1.46ms for calling a C++ binary.
Note that this benchmark only compares the two libraries to exemplify that Rust-Bio has comparable speed to C++ libraries: all used algorithms have their advantages for specific text and pattern structures and lengths (see [the Pattern Matching section in the documentation](http://johanneskoester.github.io/rust-bio/bio/pattern_matching/index.html)).

## Author

Johannes Köster (<koester@jimmy.harvard.edu>)

## Contributors

* Christopher Schröder (<christopher.schroeder@uni-due.de>)
* Peer Aramillo Irizar
* Fedor Gusev

The next name in this list could be you! If you are interested in joining the effort to build a general purpose Rust bioinformatics library, just send me an email, or issue a pull request with your first contribution.

## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
