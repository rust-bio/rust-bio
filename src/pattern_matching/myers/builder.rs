use super::*;

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
///     (b'N', &b"ACGTMRWSYKVHDB"[..])
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
#[derive(Default, Clone, Eq, PartialEq)]
pub struct MyersBuilder {
    ambigs: HashMap<u8, Vec<u8>>,
    wildcards: Vec<u8>,
}

impl MyersBuilder {
    pub fn new() -> MyersBuilder {
        Self::default()
    }

    /// Allows to specify ambiguous symbols and their equivalents. Note that the ambiguous symbol
    /// will always be matched by itself. Explicitly including it in the equivalents is not
    /// necessary.
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
        let eq = equivalents
            .into_iter()
            .map(|b| *b.borrow())
            .chain(Some(byte))
            .collect();
        self.ambigs.insert(byte, eq);
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

    /// Creates a Myers instance given a pattern, using `u64` as bit vector type
    pub fn build_64<C, P>(&self, pattern: P) -> Myers<u64>
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build(pattern)
    }

    /// Creates a Myers instance given a pattern, using `u128` as bit vector type
    #[cfg(has_u128)]
    pub fn build_128<C, P>(&self, pattern: P) -> Myers<u128>
    where
        C: Borrow<u8>,
        P: IntoIterator<Item = C>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build(pattern)
    }

    /// Creates a Myers instance given a pattern, using any desired type for bit vectors
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
        let maxsize = T::DistType::from_usize(size_of::<T>() * 8).unwrap();
        let pattern = pattern.into_iter();
        let m = T::DistType::from_usize(pattern.len()).unwrap();
        assert!(m <= maxsize, "Pattern too long");
        assert!(m > T::DistType::zero(), "Pattern is empty");

        let mut peq = [T::zero(); 256];

        for (i, a) in pattern.enumerate() {
            let mask = T::one() << i;
            // equivalent
            peq[*a.borrow() as usize] |= mask;
            // ambiguities
            if let Some(equivalents) = self.ambigs.get(a.borrow()) {
                for &eq in equivalents {
                    peq[eq as usize] |= mask;
                }
            }
        }

        for &w in &self.wildcards {
            peq[w as usize] = T::max_value();
        }

        Myers {
            peq,
            bound: T::one() << (m.to_usize().unwrap() - 1),
            m,
            tb: Traceback::new(),
        }
    }
}
