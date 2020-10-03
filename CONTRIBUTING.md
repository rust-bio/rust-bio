# Contributing to rust-bio

- [General](#general)
- [Making a PR](#making-a-pr)
  - [Creating A Branch](#creating-a-branch)
  - [Install pre-commit hooks](#install-pre-commit-hooks)
  - [Applying Changes](#applying-changes)
  - [Merging](#merging)
- [Style Guidelines](#style-guidelines)
  - [Bypassing Formatting Warnings](#bypassing-formatting-warnings)

## General

Any contributions are welcome, from a simple bug report to full-blown new modules:

If you **find a bug** and don't have the time or in-depth knowledge to fix it, just [check if you can add info to an existing issue](https://github.com/rust-bio/rust-bio/issues) and otherwise [file a bug report](https://github.com/rust-bio/rust-bio/issues/new/choose) with as many infos as possible.
If you want to contribute fixes, documentation or new code, please [open a pull request](https://github.com/rust-bio/rust-bio/compare).


If you want to contribute and don't know where to start, have a look at the [roadmap](https://github.com/rust-bio/rust-bio/issues/3).



## Making a PR

If you are new to the GitHub Fork & Pull Request Workflow, please read [this](https://guides.github.com/introduction/flow/) and/or [this](https://gist.github.com/Chaser324/ce0505fbed06b947d962) guide first. Below we summarize details that are specific to **rust-bio**.

### Creating A Branch

You have two options to do this:
1. For one-time contributions, simply [fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the repository, apply your changes to a branch in your fork and then open a pull request.
2. If you plan on contributing more than once, become a contributor by saying hi in the [`Join the team!` issue](https://github.com/rust-bio/rust-bio/issues/27) or [on the `rust-bio` Discord server](https://discord.gg/rssQABT).
    We'll add you to the team.
    Then, you don't have to create a fork, but can simply push new branches into the main repository and open pull requests there.

### Install pre-commit hooks

**If you want to contribute more than once, we strongly suggest you install pre-commit hooks.** To do so:

1. [Follow the instructions](https://pre-commit.com) to install `pre-commit` if you haven't yet. (if you have `pip`, just `pip install pre-commit`)
2. Nagivate to `rust-bio`'s project directory, and run `pre-commit install`. This should update `.git/hooks/pre-commit`

### Applying Changes

Please follow our [style guidelines](#style-guidelines).

### Merging

Once you submit your pull request on GitHub, several automated tests should be run, and their results reported on the pull request.

`cargo test`, `cargo fmt` and `cargo clippy` must be passed before a PR can be merged. This is why we strongly suggest installing pre-commit hooks which warn you about any issues in advance

## Style Guidelines

Every public function and module should have [documentation comments](https://doc.rust-lang.org/stable/rust-by-example/meta/doc.html).
Check out [which types of comments to use where](https://doc.rust-lang.org/stable/reference/comments.html#doc-comments).
In `rust-bio`, documentation comments should:
* [explain functionality](https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html)
* give at least one useful example of how to use it (best as [doctests](https://doc.rust-lang.org/rustdoc/documentation-tests.html),
  that run during testing, and using descriptive [`expect()`](https://doc.rust-lang.org/std/result/enum.Result.html#method.expect)
  statements for handling any `Err()`s that might occur)
* describe time and memory complexity listed (where applicable)
* cite and link sources and explanations for data structures, algorithms or code (where applicable)

For extra credit, feel free to familiarize yourself with:
* the Rust [documentation conventions](https://rust-lang.github.io/rfcs/1574-more-api-documentation-conventions.html#appendix-a-full-conventions-text)
* the Rust [API documentation guidelines](https://rust-lang.github.io/api-guidelines/documentation.html)

### Bypassing Formatting Warnings

Formatting warnings emitted by `cargo check` and `cargo clippy` can be bypassed using the attribute `#[allow(<rule_name>)]`. Where you absolutely need to do so, please justify. For example:

```
#![allow(non_snake_case)]
// I, D, S are conventionally used for scoring matrices.
// Only one row of the matrices are stored, i.e. I, D, S are Vec<_> instead of Vec<Vec<_>> or similar.
// In addition, single-character variable names are abundant.
// Therefore, using uppercase letter helps to better signify their meaning.
#![allow(clippy::many_single_char_names)]
// Due to the strategy used by this set of algorithms, some variables have different meanings under different
// contexts, for example: `let mut s: i32; // S[i - 1][j - 1] (before updating S[i][j]) or S[i][j - 1] (after updating S[i][j])`
// Therefore, it is difficult to give them unchanged and descriptive names.
// The meanings of these variable names are explained in `Aligner.global`
```
