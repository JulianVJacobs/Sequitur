use std::collections::HashMap;

/// Default minimum suffix length used by the Python prototype.
pub const DEFAULT_MIN_SUFFIX_LEN: usize = 3;

/// Distinguishes between suffix and (reverse) prefix entries when sorting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AffixKind {
    /// A suffix extracted from the tail of a read.
    Suffix,
    /// A prefix extracted from the head of a read (encoded with the `^` marker).
    Prefix,
}

impl AffixKind {
    #[inline]
    fn marker(self) -> char {
        match self {
            Self::Suffix => '$',
            Self::Prefix => '^',
        }
    }
}

/// A single entry in the affix array (suffix array plus mirrored prefixes).
///
/// The `text` field mirrors the Python implementation: the affix sequence is
/// followed by a marker (`$` for suffixes, `^` for prefixes) and the decimal
/// read index. This keeps lexicographic ordering identical across languages
/// while tracking the additional metadata in typed fields for downstream code.
#[derive(Debug, Clone)]
pub struct AffixEntry {
    text: String,
    marker_index: usize,
    kind: AffixKind,
    read_index: usize,
}

impl AffixEntry {
    fn new(sequence: &str, kind: AffixKind, read_index: usize, index_repr: &str) -> Self {
        let marker = kind.marker();
        let mut text = String::with_capacity(sequence.len() + 1 + index_repr.len());
        text.push_str(sequence);
        text.push(marker);
        text.push_str(index_repr);

        Self {
            text,
            marker_index: sequence.len(),
            kind,
            read_index,
        }
    }

    /// Full decorated key (`suffix$idx` or `prefix^idx`).
    pub fn key(&self) -> &str {
        &self.text
    }

    /// Returns the affix substring without the marker/index decoration.
    pub fn affix(&self) -> &str {
        &self.text[..self.marker_index]
    }

    /// Read identifier associated with this entry.
    pub fn read_index(&self) -> usize {
        self.read_index
    }

    /// Whether this entry originated from a suffix or prefix scan.
    pub fn kind(&self) -> AffixKind {
        self.kind
    }

    /// Length of the affix sequence prior to the marker.
    pub fn span(&self) -> usize {
        self.marker_index
    }
}

/// Sorted affix array enriched with constant-time lookups.
#[derive(Debug, Clone)]
pub struct AffixArray {
    entries: Vec<AffixEntry>,
    lookup: HashMap<String, usize>,
}

impl AffixArray {
    /// Construct the affix array mirroring the Python prototype semantics.
    pub fn build<I, S>(reads: I, min_suffix_len: usize) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let min_len = min_suffix_len.max(1);
        let mut entries: Vec<AffixEntry> = Vec::new();

        for (read_index, read) in reads.into_iter().enumerate() {
            let read = read.as_ref();
            let read_len = read.len();
            if read_len < min_len {
                continue;
            }

            let index_repr = read_index.to_string();
            let variants = read_len - min_len + 1;
            entries.reserve(variants * 2);

            for start in 0..=read_len - min_len {
                let suffix = &read[start..];
                entries.push(AffixEntry::new(
                    suffix,
                    AffixKind::Suffix,
                    read_index,
                    &index_repr,
                ));
            }

            for span in (min_len..=read_len).rev() {
                let prefix = &read[..span];
                entries.push(AffixEntry::new(
                    prefix,
                    AffixKind::Prefix,
                    read_index,
                    &index_repr,
                ));
            }
        }

        entries.sort_by(|lhs, rhs| lhs.key().cmp(rhs.key()));

        let mut lookup = HashMap::with_capacity(entries.len());
        for (idx, entry) in entries.iter().enumerate() {
            lookup.insert(entry.key().to_string(), idx);
        }

        Self { entries, lookup }
    }

    /// Immutable view of the sorted entries.
    pub fn entries(&self) -> &[AffixEntry] {
        &self.entries
    }

    /// Iterator over the entries in sorted order.
    pub fn iter(&self) -> impl Iterator<Item = &AffixEntry> {
        self.entries.iter()
    }

    /// Look up the position of an entry by its decorated key.
    pub fn index_of(&self, key: &str) -> Option<usize> {
        self.lookup.get(key).copied()
    }

    /// Access an entry by its position in the sorted array.
    pub fn get(&self, index: usize) -> Option<&AffixEntry> {
        self.entries.get(index)
    }

    /// Number of entries stored in the array.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when the array contains no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

impl Default for AffixArray {
    fn default() -> Self {
        Self {
            entries: Vec::new(),
            lookup: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_affix_array_for_single_read() {
        let reads = vec!["ACGT".to_string()];
        let array = AffixArray::build(reads, DEFAULT_MIN_SUFFIX_LEN);

        let keys: Vec<&str> = array.iter().map(|entry| entry.key()).collect();
        assert_eq!(keys, vec!["ACGT$0", "ACGT^0", "ACG^0", "CGT$0"]);

        let (first, second, third, fourth) = (
            array.get(0).unwrap(),
            array.get(1).unwrap(),
            array.get(2).unwrap(),
            array.get(3).unwrap(),
        );

        assert_eq!(first.affix(), "ACGT");
        assert_eq!(first.kind(), AffixKind::Suffix);
        assert_eq!(second.affix(), "ACGT");
        assert_eq!(second.kind(), AffixKind::Prefix);
        assert_eq!(third.affix(), "ACG");
        assert_eq!(third.kind(), AffixKind::Prefix);
        assert_eq!(fourth.affix(), "CGT");
        assert_eq!(fourth.kind(), AffixKind::Suffix);
    }

    #[test]
    fn skips_reads_shorter_than_threshold() {
        let reads = vec!["AC".to_string(), "ACGT".to_string()];
        let array = AffixArray::build(reads, DEFAULT_MIN_SUFFIX_LEN);
        assert_eq!(array.len(), 4);

        let keys: Vec<&str> = array.iter().map(|entry| entry.key()).collect();
        assert!(keys.iter().all(|key| key.ends_with('1')));
    }

    #[test]
    fn supports_variable_min_suffix_lengths() {
        let reads = vec!["AACGT".to_string()];
        let array = AffixArray::build(&reads, 2);
        assert_eq!(array.len(), 8);

        assert!(array.index_of("CGT$0").is_some());
        assert!(array.index_of("AACGT^0").is_some());
        assert!(array.index_of("GT$0").is_some());
        assert!(array.index_of("AA^0").is_some());
    }

    #[test]
    fn treats_zero_min_suffix_as_one() {
        let reads = vec!["AA".to_string()];
        let array = AffixArray::build(reads, 0);
        assert_eq!(array.len(), 4);
        assert!(array.index_of("AA$0").is_some());
        assert!(array.index_of("A^0").is_some());
    }
}
