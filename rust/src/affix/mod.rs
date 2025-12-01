//! Unified affix structure interface supporting both array and trie implementations.

// Import from submodules
mod affix_array;
mod affix_trie;

// Re-export common types
pub use self::affix_array::{AffixKey, AffixKind, AffixMap, DEFAULT_MIN_SUFFIX_LEN};
pub use self::affix_trie::{AffixNode, AffixSlice, Direction, PrunedAffixTrie};

/// Unified affix structure that can be either an array or a trie.
#[derive(Debug, Clone)]
pub enum AffixStructure<'a> {
    /// HashMap-based affix array (original implementation)
    Array(AffixMap<'a>),
    /// Pruned and compressed trie (new optimized implementation)
    Trie(PrunedAffixTrie),
}

impl<'a> AffixStructure<'a> {
    /// Build an affix structure using the specified implementation.
    pub fn build(reads: &'a [String], min_suffix_len: usize, use_trie: bool) -> Self {
        if use_trie {
            Self::Trie(PrunedAffixTrie::build(reads, min_suffix_len))
        } else {
            Self::Array(AffixMap::build(reads, min_suffix_len))
        }
    }

    /// Build using the array implementation.
    pub fn build_array(reads: &'a [String], min_suffix_len: usize) -> Self {
        Self::Array(AffixMap::build(reads, min_suffix_len))
    }

    /// Build using the trie implementation.
    pub fn build_trie(reads: &'a [String], min_suffix_len: usize) -> Self {
        Self::Trie(PrunedAffixTrie::build(reads, min_suffix_len))
    }

    /// Check if this is a trie implementation.
    pub fn is_trie(&self) -> bool {
        matches!(self, Self::Trie(_))
    }

    /// Check if this is an array implementation.
    pub fn is_array(&self) -> bool {
        matches!(self, Self::Array(_))
    }

    /// Get the underlying AffixMap if this is an array, panics otherwise.
    pub fn as_array(&self) -> &AffixMap<'a> {
        match self {
            Self::Array(map) => map,
            Self::Trie(_) => panic!("Cannot get array from trie structure"),
        }
    }

    /// Get the underlying PrunedAffixTrie if this is a trie, panics otherwise.
    pub fn as_trie(&self) -> &PrunedAffixTrie {
        match self {
            Self::Trie(trie) => trie,
            Self::Array(_) => panic!("Cannot get trie from array structure"),
        }
    }
}
