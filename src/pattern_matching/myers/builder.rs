use std::borrow::Borrow;
use std::collections::HashMap;

use super::long::Myers as MyersLong;
use super::{BitVec, Myers};

/// Builds a Myers instance, allowing to specify ambiguities.
///
/// # Example:
///
/// This example shows how recognition of IUPAC ambiguities in patterns can be implemented:
///
/// ```
/// # extern crate bio;
/// use bio::pattern_matching::myers::MyersBuilder;
///
/// # fn main() {
/// let ambigs = [
///     (b'M', &b"AC"[..]),
///     (b'R', &b"AG"[..]),
///     (b'W', &b"AT"[..]),
///     (b'S', &b"CG"[..]),
///     (b'Y', &b"CT"[..]),
///     (b'K', &b"GT"[..]),
///     (b'V', &b"ACGMRS"[..]),
///     (b'H', &b"ACTMWY"[..]),
///     (b'D', &b"AGTRWK"[..]),
///     (b'B', &b"CGTSYK"[..]),
///     (b'N', &b"ACGTMRWSYKVHDB"[..]),
/// ];
///
/// let mut builder = MyersBuilder::new();
///
/// for &(base, equivalents) in &ambigs {
///     builder.ambig(base, equivalents);
/// }
///
/// let text = b"GGATGNGCGCCATAG";
/// let pattern = b"TRANCGG";
/// //                *   * (mismatch)
///
/// let myers = builder.build_64(pattern);
/// assert_eq!(myers.distance(text), 2);
/// # }
/// ```
///
/// Note that only ambiguities in the pattern are recognized. The reverse is not true; ambiguities
/// in the search text are not matched by multiple symbols in the pattern. This would require
/// specifying additional ambiguities (`builder.ambig(b'A', b"MRWVHDN")`, etc...).
#[derive(Default, Clone, Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct MyersBuilder {
    ambigs: HashMap<u8, Vec<u8>>,
    wildcards: Vec<u8>,
}

impl MyersBuilder {
    pub fn new() -> MyersBuilder {
        Self::default()
    }

    /// Allows to specify ambiguous symbols in the pattern, which will match all
    /// given equivalent symbols in the text, as well as itself.
    ///
    /// The `ambig()` method can be called multiple times with the same symbol
    /// and (potentially) different equivalents, which are simply added to the
    /// existing list.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::MyersBuilder;
    ///
    /// # fn main() {
    /// let text = b"GGATGAGCGCCATAG";
    /// let pattern = b"TGAGCGN";
    ///
    /// let myers = MyersBuilder::new()
    ///     .ambig(b'N', b"ACGT")
    ///     .build_64(pattern);
    ///
    /// assert_eq!(myers.distance(text), 0);
    /// # }
    pub fn ambig<I, B>(&mut self, byte: u8, equivalents: I) -> &mut Self
    where
        I: IntoIterator<Item = B>,
        B: Borrow<u8>,
    {
        let eq = self.ambigs.entry(byte).or_default();
        eq.extend(equivalents.into_iter().map(|b| *b.borrow()));
        self
    }

    /// Allows to specify a wildcard symbol, that upon appearance in the search text
    /// shall be matched by any symbol of the pattern. Multiple wildcards are possible.
    /// For the inverse, that is, wildcards in the pattern matching any symbol in search
    /// text, use `ambig(byte, 0..255)`.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::MyersBuilder;
    ///
    /// # fn main() {
    /// let text = b"GGATGAGCG*CATAG";
    /// let pattern = b"TGAGCGT";
    ///
    /// let myers = MyersBuilder::new()
    ///     .text_wildcard(b'*')
    ///     .build_64(pattern);
    ///
    /// assert_eq!(myers.distance(text), 0);
    /// # }
    pub fn text_wildcard(&mut self, wildcard: u8) -> &mut Self {
        self.wildcards.push(wildcard);
        self
    }

    /// Creates a Myers instance given a pattern, using `u64` as bit vector type.
    /// Pattern length is restricted to at most 64 symbols.
    pub fn build_64<C, P>(&self, pattern: P) -> Myers<u64>
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build(pattern)
    }

    /// Creates a Myers instance given a pattern, using `u128` as bit vector type.
    /// Pattern length is restricted to at most 128 symbols.
    pub fn build_128<C, P>(&self, pattern: P) -> Myers<u128>
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build(pattern)
    }

    /// Creates a Myers instance given a pattern, using any desired type for bit vectors.
    /// Pattern length is restricted to the size of the bit vector `T`.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::{MyersBuilder, Myers};
    ///
    /// # fn main() {
    /// let myers: Myers<u32> = MyersBuilder::new()
    ///     .text_wildcard(b'*')
    ///     .build(b"TGAGCG*");
    /// // ...
    /// # }
    pub fn build<T, C, P>(&self, pattern: P) -> Myers<T>
    where
        T: BitVec,
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        Myers::new_ambig(pattern, Some(&self.ambigs), Some(&self.wildcards))
    }

    /// Creates a `long::Myers` instance given a pattern, using `u64` as bit vector type.
    /// Pattern length is not restricted regardless of the type of the bit vector.
    pub fn build_long_64<C, P>(&self, pattern: P) -> MyersLong<u64>
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build_long(pattern)
    }

    /// Creates a `long::Myers` instance given a pattern, using `u128` as bit vector type.
    /// Pattern length is not restricted regardless of the type of the bit vector.
    pub fn build_long_128<C, P>(&self, pattern: P) -> MyersLong<u128>
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build_long(pattern)
    }

    /// Creates a `long::Myers` instance given a pattern, using any desired type for bit vectors.
    /// Pattern length is not restricted regardless of the type of the bit vector.
    pub fn build_long<T, C, P>(&self, pattern: P) -> MyersLong<T>
    where
        T: BitVec,
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        MyersLong::new_ambig(pattern, Some(&self.ambigs), Some(&self.wildcards))
    }
}
