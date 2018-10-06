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
///     (b'M', &b"ACM"[..]),
///     (b'R', &b"AGR"[..]),
///     (b'W', &b"ATW"[..]),
///     (b'S', &b"CGS"[..]),
///     (b'Y', &b"CTY"[..]),
///     (b'K', &b"GTK"[..]),
///     (b'V', &b"ACGMRSV"[..]),
///     (b'H', &b"ACTMWYH"[..]),
///     (b'D', &b"AGTRWKD"[..]),
///     (b'B', &b"CGTSYKB"[..]),
///     (b'N', &b"ACGTMRWSYKVHDBN"[..])
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
/// let myers = builder.build(pattern);
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

    /// Allows to specify ambiguous characters and their equivalents.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::MyersBuilder;
    ///
    /// # fn main() {
    /// let text = b"ACCGTGGATGAGCGCCATAG";
    /// let pattern =      b"TGAGCGN";
    ///
    /// let myers = MyersBuilder::new()
    ///     .ambig(b'N', b"ACGT")
    ///     .build(pattern);
    ///
    /// assert_eq!(myers.distance(text), 0);
    /// # }
    pub fn ambig<I, B>(&mut self, byte: u8, equivalents: I) -> &mut Self
    where
        I: IntoIterator<Item = B>,
        B: Borrow<u8>,
    {
        let eq = equivalents.into_iter().map(|b| *b.borrow()).collect();
        self.ambigs.insert(byte, eq);
        self
    }

    /// Allows to specify a wildcard character, that upon appearance in the search text
    /// shall be matched by any character of the pattern. Multiple wildcards are possible.
    /// For the inverse, that is, wildcards in the pattern matching any character in search
    /// text, use `ambig(byte, 0..255)`.
    ///
    /// # Example:
    ///
    /// ```
    /// # extern crate bio;
    /// use bio::pattern_matching::myers::MyersBuilder;
    ///
    /// # fn main() {
    /// let text = b"ACCGTGGATGAGCG*CATAG";
    /// let pattern =      b"TGAGCGT";
    ///
    /// let myers = MyersBuilder::new()
    ///     .text_wildcard(b'*')
    ///     .build(pattern);
    ///
    /// assert_eq!(myers.distance(text), 0);
    /// # }
    pub fn text_wildcard(&mut self, wildcard: u8) -> &mut Self {
        self.wildcards.push(wildcard);
        self
    }

    /// Creates a Myers instance given a pattern, using `u64` as bit vector type
    pub fn build<'a, P>(&self, pattern: P) -> Myers<u64>
    where
        P: IntoTextIterator<'a>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build_other(pattern)
    }

    /// Creates a Myers instance given a pattern, using `u128` as bit vector type
    #[cfg(has_u128)]
    pub fn build128<'a, P>(&self, pattern: P) -> Myers<u128>
    where
        P: IntoTextIterator<'a>,
        P::IntoIter: ExactSizeIterator,
    {
        self.build_other(pattern)
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
    ///     .build_other(b"TGAGCG*");
    /// // ...
    /// # }
    pub fn build_other<'a, T, P>(&self, pattern: P) -> Myers<T>
    where
        T: BitVec,
        P: IntoTextIterator<'a>,
        P::IntoIter: ExactSizeIterator,
    {
        let maxsize = T::DistType::from_usize(size_of::<T>() * 8).unwrap();
        let pattern = pattern.into_iter();
        let m = T::DistType::from_usize(pattern.len()).unwrap();
        assert!(m <= maxsize, "Pattern too long");
        assert!(m > T::DistType::zero(), "Pattern is empty");

        let mut peq = [T::zero(); 256];

        for (i, &a) in pattern.enumerate() {
            let mask = T::one() << i;
            // equivalent
            peq[a as usize] |= mask;
            // ambiguities
            if let Some(equivalents) = self.ambigs.get(&a) {
                for &eq in equivalents {
                    peq[eq as usize] |= mask;
                }
            }
        }

        for &w in &self.wildcards {
            peq[w as usize] = T::max_value();
        }

        Myers {
            peq: peq,
            bound: T::one() << (m.to_usize().unwrap() - 1),
            m: m,
            tb: Traceback::new(),
        }
    }
}
