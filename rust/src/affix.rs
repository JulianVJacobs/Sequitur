use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub struct AffixKey<'a> {
    pub affix: &'a str,
    pub kind: AffixKind,
}

#[derive(Debug, Clone)]
pub struct AffixMap<'a> {
    pub keys: Vec<AffixKey<'a>>,
    pub map: HashMap<AffixKey<'a>, Vec<usize>>,
}

impl<'a> AffixMap<'a> {
    pub fn build(reads: &'a [String], min_suffix_len: usize) -> Self {
        let min_len = min_suffix_len.max(1);
        let mut map: HashMap<AffixKey<'a>, Vec<usize>> = HashMap::new();

        for (read_index, read) in reads.iter().enumerate() {
            let read_len = read.len();
            if read_len < min_len {
                continue;
            }

            for start in 0..=read_len - min_len {
                let affix = &read[start..];
                let key = AffixKey {
                    affix,
                    kind: AffixKind::Suffix,
                };
                map.entry(key).or_insert_with(Vec::new).push(read_index);
            }

            for span in min_len..=read_len {
                let affix = &read[..span];
                let key = AffixKey {
                    affix,
                    kind: AffixKind::Prefix,
                };
                map.entry(key).or_insert_with(Vec::new).push(read_index);
            }
        }

        let mut keys: Vec<AffixKey<'a>> = map.keys().cloned().collect();
        keys.sort_unstable();

        AffixMap { keys, map }
    }

    pub fn iter(&self) -> impl Iterator<Item = &AffixKey<'a>> {
        self.keys.iter()
    }

    pub fn reads_for(&self, key: &AffixKey<'a>) -> Option<&Vec<usize>> {
        self.map.get(key)
    }
}
use std::cmp::Ordering;

const INVALID_ANCHOR: usize = usize::MAX;

/// Default minimum suffix length used by the Python prototype.
pub const DEFAULT_MIN_SUFFIX_LEN: usize = 3;

/// Distinguishes between suffix and (reverse) prefix entries when sorting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub enum AffixKind {
    /// A suffix extracted from the tail of a read.
    Suffix,
    /// A prefix extracted from the head of a read (encoded with the `^` marker).
    Prefix,
}

impl AffixKind {
    #[inline]
    #[allow(dead_code)]
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
#[derive(Debug, Clone, Copy)]
pub struct AffixEntry {
    read_index: usize,
    offset: usize,
    span: usize,
    kind: AffixKind,
}

impl AffixEntry {
    fn new_suffix(read_index: usize, offset: usize, span: usize) -> Self {
        Self {
            read_index,
            offset,
            span,
            kind: AffixKind::Suffix,
        }
    }

    fn new_prefix(read_index: usize, span: usize) -> Self {
        Self {
            read_index,
            offset: 0,
            span,
            kind: AffixKind::Prefix,
        }
    }

    pub fn read_index(&self) -> usize {
        self.read_index
    }

    pub fn kind(&self) -> AffixKind {
        self.kind
    }

    pub fn span(&self) -> usize {
        self.span
    }

    pub fn offset(&self) -> usize {
        self.offset
    }

    pub fn affix<'a>(&self, reads: &'a [String]) -> &'a str {
        let read = &reads[self.read_index];
        let end = self.offset + self.span;
        &read[self.offset..end]
    }
}

/// Sorted affix array enriched with constant-time lookups.
#[derive(Debug, Clone)]
pub struct AffixArray {
    entries: Vec<AffixEntry>,
    suffix_anchors: Vec<Vec<usize>>,
}

impl AffixArray {
    /// Construct the affix array mirroring the Python prototype semantics.
    pub fn build(reads: &[String], min_suffix_len: usize) -> Self {
        let min_len = min_suffix_len.max(1);
        let mut entries: Vec<AffixEntry> = Vec::new();
        let mut suffix_anchors: Vec<Vec<usize>> = vec![Vec::new(); reads.len()];

        for (read_index, read) in reads.iter().enumerate() {
            let read_len = read.len();
            if read_len < min_len {
                continue;
            }

            let suffix_slots = read_len - min_len + 1;
            suffix_anchors[read_index] = vec![INVALID_ANCHOR; suffix_slots];
            entries.reserve(suffix_slots * 2);

            for start in 0..=read_len - min_len {
                let span = read_len - start;
                entries.push(AffixEntry::new_suffix(read_index, start, span));
            }

            for span in (min_len..=read_len).rev() {
                entries.push(AffixEntry::new_prefix(read_index, span));
            }
        }

        entries.sort_by(|lhs, rhs| compare_entries(lhs, rhs, reads));

        for (idx, entry) in entries.iter().enumerate() {
            if entry.kind() == AffixKind::Suffix {
                if let Some(row) = suffix_anchors.get_mut(entry.read_index()) {
                    if entry.offset() < row.len() {
                        row[entry.offset()] = idx;
                    }
                }
            }
        }

        Self {
            entries,
            suffix_anchors,
        }
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
            suffix_anchors: Vec::new(),
        }
    }
}

fn compare_entries(lhs: &AffixEntry, rhs: &AffixEntry, reads: &[String]) -> Ordering {
    let lhs_seq = lhs.affix(reads).as_bytes();
    let rhs_seq = rhs.affix(reads).as_bytes();
    match lhs_seq.cmp(rhs_seq) {
        Ordering::Equal => {
            // Shorter sequence first
            match lhs_seq.len().cmp(&rhs_seq.len()) {
                Ordering::Equal => {
                    // Suffix before prefix
                    match Ord::cmp(&lhs.kind(), &rhs.kind()) {
                        Ordering::Equal => lhs.read_index().cmp(&rhs.read_index()),
                        ord => ord,
                    }
                }
                ord => ord,
            }
        }
        ord => ord,
    }
}

#[allow(dead_code)]
fn decorated_bytes<'a>(
    entry: &'a AffixEntry,
    reads: &'a [String],
) -> impl Iterator<Item = u8> + 'a {
    entry
        .affix(reads)
        .bytes()
        .chain(std::iter::once(entry.kind().marker() as u8))
        .chain(IndexDigits::new(entry.read_index()))
}

#[allow(dead_code)]
struct IndexDigits {
    buf: [u8; 20],
    pos: usize,
    end: usize,
}

impl IndexDigits {
    #[allow(dead_code)]
    fn new(mut value: usize) -> Self {
        let mut buf = [0u8; 20];
        let mut idx = buf.len();
        if value == 0 {
            idx -= 1;
            buf[idx] = b'0';
        } else {
            while value > 0 {
                idx -= 1;
                buf[idx] = b'0' + (value % 10) as u8;
                value /= 10;
            }
        }
        Self {
            buf,
            pos: idx,
            end: buf.len(),
        }
    }
}

impl Iterator for IndexDigits {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.end {
            None
        } else {
            let b = self.buf[self.pos];
            self.pos += 1;
            Some(b)
        }
    }
}

impl AffixArray {
    pub fn suffix_anchor(&self, read_index: usize, start: usize) -> Option<usize> {
        self.suffix_anchors
            .get(read_index)
            .and_then(|row| row.get(start))
            .and_then(|&idx| {
                if idx == INVALID_ANCHOR {
                    None
                } else {
                    Some(idx)
                }
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_affix_map_for_single_read() {
        let reads = vec!["ACGT".to_string()];
        let map = AffixMap::build(&reads, DEFAULT_MIN_SUFFIX_LEN);
        let keys: Vec<String> = map
            .iter()
            .map(|key| format!("{}{:?}", key.affix, key.kind))
            .collect();
        assert!(keys.contains(&"ACGTSuffix".to_string()));
        assert!(keys.contains(&"ACGTPrefix".to_string()));
        assert!(keys.contains(&"ACGPrefix".to_string()));
        assert!(keys.contains(&"CGTSuffix".to_string()));
        // Validate mapping
        for key in map.iter() {
            let indices = map.reads_for(key).unwrap();
            assert_eq!(indices, &vec![0]);
        }
    }

    #[test]
    fn skips_reads_shorter_than_threshold() {
        let reads = vec!["AC".to_string(), "ACGT".to_string()];
        let map = AffixMap::build(&reads, DEFAULT_MIN_SUFFIX_LEN);
        let keys: Vec<&AffixKey> = map.iter().collect();
        // Only keys for the second read (index 1) should exist
        for key in keys {
            let indices = map.reads_for(key).unwrap();
            assert!(indices.iter().all(|&idx| idx == 1));
        }
    }

    #[test]
    fn supports_variable_min_suffix_lengths() {
        let reads = vec!["AACGT".to_string()];
        let map = AffixMap::build(&reads, 2);
        // Should contain all expected affixes
        let expected_affixes = vec!["AACGT", "AACG", "ACGT", "CGT", "AAC", "GT", "AA"];
        for affix in expected_affixes {
            assert!(
                map.iter().any(|key| key.affix == affix
                    && (key.kind == AffixKind::Prefix || key.kind == AffixKind::Suffix)),
                "Missing affix: {}",
                affix
            );
        }
        // All keys should map to read index 0
        for key in map.iter() {
            let indices = map.reads_for(key).unwrap();
            assert_eq!(indices, &vec![0]);
        }
    }

    #[test]
    fn treats_zero_min_suffix_as_one() {
        let reads = vec!["AA".to_string()];
        let map = AffixMap::build(&reads, 0);
        // Should contain all expected affixes
        let expected_affixes = vec!["AA", "A"];
        for affix in expected_affixes {
            assert!(map.iter().any(|key| key.affix == affix));
        }
        // All keys should map to read index 0
        for key in map.iter() {
            let indices = map.reads_for(key).unwrap();
            assert_eq!(indices, &vec![0]);
        }
    }
}
