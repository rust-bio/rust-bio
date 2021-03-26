# Contribute

Any contributions are welcome, from a simple bug report to full-blown new modules:

If you **find a bug** and don't have the time or in-depth knowledge to fix it, just [check if you can add info to an existing issue](https://github.com/rust-bio/rust-bio/issues) and otherwise [file a bug report](https://github.com/rust-bio/rust-bio/issues/new/choose) with as many infos as possible.

Pull requests are welcome if you want to contribute fixes, documentation, or new code.
Before making commits, it would be helpful to first install `pre-commit` to avoid failed continuous integration builds due to issues such as formatting:
1. Install `pre-commit` (see [pre-commit.com/#installation](https://pre-commit.com/#installation))
2. Run `pre-commit install` in the rust-bio base directory

Depending on your intended contribution frequency, you have two options for opening pull requests:
1. For one-time contributions, simply [fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the repository, apply your changes to a branch in your fork and then open a pull request.
2. If you plan on contributing more than once, become a contributor by saying hi [on the `rust-bio` Discord server](https://discord.gg/rssQABT), together with a short sentence saying who you are and mentioning what you want to contribute.
   We'll add you to the team.
   Then, you don't have to create a fork, but can simply push new branches into the main repository and open pull requests there.

If you want to contribute and don't know where to start, have a look at the [roadmap](https://github.com/rust-bio/rust-bio/issues/3), check if you could review a [pull request](https://github.com/rust-bio/rust-bio/pulls) or tackle an [existing issue](https://github.com/rust-bio/rust-bio/issues).


## Documentation guidelines

Every public function and module should have [documentation comments](https://doc.rust-lang.org/stable/rust-by-example/meta/doc.html).
Check out [which types of comments to use where](https://doc.rust-lang.org/stable/reference/comments.html#doc-comments).

In `rust-bio`, documentation comments should:
* [explain functionality](https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html)
* give at least one useful example of how to use it (best as [doctests](https://doc.rust-lang.org/rustdoc/documentation-tests.html), that run during testing, and using descriptive [`expect()`](https://doc.rust-lang.org/std/result/enum.Result.html#method.expect)
  statements for handling any `Err()`s that might occur)
* describe time and memory complexity listed (where applicable)
* cite and link sources and explanations for data structures, algorithms or code (where applicable)

For extra credit, feel free to familiarize yourself with:
* the Rust [documentation conventions](https://rust-lang.github.io/rfcs/1574-more-api-documentation-conventions.html#appendix-a-full-conventions-text)
* the Rust [API documentation guidelines](https://rust-lang.github.io/api-guidelines/documentation.html)

### variable names and mathematical expressions

* For mathematical expressions in doc strings (e.g. for "big-O" notation), you can use LaTeX syntax in between dollar signs (e.g. `$O(n)$`) in-line.
* For simple expressions that refer to a variable as an argument of a function or to a field of a struct, using backticks (`\`field_1 and var_x\``) is recommended.

In the following example, the public function `combinations` has two arguments `n` and `k`, which are referred to in the first line with backticks.
Then, a LaTeX math expressions is used to show the complexity.
The variable names *k* and *n* are not ambiguous (obviously referring to the two arguments), so we don't need to define again.

```rust
/// Calculate the number of combinations when choosing
/// `k` elements from `n` elements without replacement.
///
/// Time complexity: $O(\min(k, n - k))$
pub fn combinations(n: u64, k: u64) -> f64 {
    scaled_combinations(n, k, 1.0)
}
```

In the next example, none of the variables *n*, *k* and *w* used in the space complexity expression matches any field names of the struct `BitEnc`, so we need to define explicitly.
Simple expressions like `width = 7` and `4 * 7 = 28` are enclosed in backticks.

```rust
/// A sequence of bitencoded values.
///
/// Space complexity: $O(\dfrac{nw}{k})$; specifically $\Big\lceil\dfrac{nw}{k}\Big\rceil \times 32$
/// bits, where $w$ is `width`, $n$ is the length of the input
/// sequence, and $k = 32 - (32 \mod w)$ is the number of bits in each
/// 32-bit block that can be used to store values.
/// For values that are not a divider of 32, some bits will remain unused.
/// For example for `width = 7` only `4 * 7 = 28` bits are used.
#[derive(Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct BitEnc {
    storage: Vec<u32>,
    width: usize,
    mask: u32,
    len: usize,
    usable_bits_per_block: usize,
}
```

