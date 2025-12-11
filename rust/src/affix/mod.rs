//! Affix trie interface (array implementation removed).

// Import trie submodule
mod affix_trie;

// Re-export trie types
pub use self::affix_trie::{AffixNode, AffixSlice, Direction, PrunedAffixTrie};

/// Default minimum suffix length (migrated from the removed array implementation).
pub const DEFAULT_MIN_SUFFIX_LEN: usize = 3;

// Note: The previous `AffixStructure` enum and array-based helpers were removed
// in favor of the trie-only `PrunedAffixTrie`. Callers should use
// `PrunedAffixTrie::build(reads, min_k, max_diff)` to construct the trie.
