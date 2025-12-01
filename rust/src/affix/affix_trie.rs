//! Compressed pruned affix trie for efficient overlap candidate generation.
//!
//! This implementation uses a map-like structure where affixes are represented as
//! (read_idx, start, end) slices rather than storing full strings. Nodes are
//! progressively promoted from Leaf → Chain → Delta as branching is detected.

use std::collections::{HashMap, HashSet};

/// A key representing an affix as a slice reference into a read.
pub type AffixSlice = (usize, usize, usize); // (read_idx, start, end)

/// Direction of affix extension (growing toward full read).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Prefix direction: growing rightward (adding chars to end)
    Prefix,
    /// Suffix direction: growing leftward (adding chars to start)
    Suffix,
}

/// Tracks which reads terminate at a specific affix (for containment/duplicate detection).
/// Optimized for the common case of no roots (empty 99%+ of the time).
/// Memory: Empty=8B, Single=16B, Multiple=32B+data (vs always 24B for Vec<usize>)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RootsList {
    /// No reads terminate here (common case: 99%+)
    Empty,
    /// Exactly one read terminates here (rare: single duplicate or contained read)
    Single(usize),
    /// Multiple reads terminate here (extremely rare: multiple exact duplicates)
    Multiple(Vec<usize>),
}

impl RootsList {
    fn new() -> Self {
        Self::Empty
    }

    fn add(&mut self, read_idx: usize) {
        match self {
            Self::Empty => *self = Self::Single(read_idx),
            Self::Single(existing) if *existing != read_idx => {
                *self = Self::Multiple(vec![*existing, read_idx]);
            }
            Self::Single(_) => {} // Already contains this read
            Self::Multiple(vec) => {
                if !vec.contains(&read_idx) {
                    vec.push(read_idx);
                }
            }
        }
    }

    fn contains(&self, read_idx: usize) -> bool {
        match self {
            Self::Empty => false,
            Self::Single(id) => *id == read_idx,
            Self::Multiple(vec) => vec.contains(&read_idx),
        }
    }

    fn is_empty(&self) -> bool {
        matches!(self, Self::Empty)
    }

    fn iter(&self) -> RootsListIter<'_> {
        RootsListIter {
            roots: self,
            index: 0,
        }
    }
}

pub struct RootsListIter<'a> {
    roots: &'a RootsList,
    index: usize,
}

impl<'a> Iterator for RootsListIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match self.roots {
            RootsList::Empty => None,
            RootsList::Single(id) if self.index == 0 => {
                self.index += 1;
                Some(*id)
            }
            RootsList::Single(_) => None,
            RootsList::Multiple(vec) => {
                if self.index < vec.len() {
                    let result = vec[self.index];
                    self.index += 1;
                    Some(result)
                } else {
                    None
                }
            }
        }
    }
}

/// A node in the affix trie, progressively promoted based on branching.
#[derive(Debug, Clone)]
pub enum AffixNode {
    /// Single extension in one or both directions - a collapsible chain link.
    Chain {
        prefix_link: Option<AffixSlice>,
        suffix_link: Option<AffixSlice>,
    },

    /// Multiple extensions or contains full read(s) - non-collapsible branching point.
    Delta {
        prefix_extensions: Vec<AffixSlice>,
        suffix_extensions: Vec<AffixSlice>,
        roots: RootsList, // Full reads that end at this exact affix (containment edge case)
    },

    /// Terminal node with no extensions.
    Leaf,
}

/// Pruned and compressed affix trie for overlap detection.
#[derive(Debug, Clone)]
pub struct PrunedAffixTrie {
    /// All nodes in the trie (key = affix slice reference)
    pub nodes: HashMap<AffixSlice, AffixNode>,

    /// Root k-mers (affixes of length min_k)
    pub roots: Vec<AffixSlice>,

    /// Minimum k-mer length
    pub min_k: usize,

    /// Global canonicalization: sequence → first key encountered
    sequence_to_key: HashMap<String, AffixSlice>,

    /// Track which reads contribute to each affix sequence
    affix_contributors: HashMap<String, HashSet<usize>>,
}

impl PrunedAffixTrie {
    /// Build a pruned and compressed affix trie from reads.
    pub fn build(reads: &[String], min_k: usize) -> Self {
        let min_len = min_k.max(1);
        let mut trie = Self {
            nodes: HashMap::new(),
            roots: Vec::new(),
            min_k: min_len,
            sequence_to_key: HashMap::new(),
            affix_contributors: HashMap::new(),
        };

        // Phase 1a: Insert all affixes with canonicalization
        for (read_idx, read) in reads.iter().enumerate() {
            let read_len = read.len();
            if read_len < min_len {
                continue;
            }

            // Generate PREFIXES: start=0, end varies from min_len to read_len
            for end in min_len..=read_len {
                let sequence = &read[0..end];

                // Track that this read has this affix
                trie.affix_contributors
                    .entry(sequence.to_string())
                    .or_default()
                    .insert(read_idx);

                // Canonicalize: reuse existing key if sequence seen before
                let key = if let Some(&existing_key) = trie.sequence_to_key.get(sequence) {
                    existing_key
                } else {
                    let new_key = (read_idx, 0, end);
                    trie.nodes.insert(new_key, AffixNode::Leaf);
                    trie.sequence_to_key.insert(sequence.to_string(), new_key);

                    // Track k-mers (minimum length affixes)
                    if end == min_len {
                        trie.roots.push(new_key);
                    }

                    new_key
                };

                // Phase 1b: Link to PREFIX parent (longer prefix = one char added right/end)
                if end < read_len {
                    let parent_seq = &read[0..end + 1];
                    if let Some(&parent_key) = trie.sequence_to_key.get(parent_seq) {
                        Self::insert_extension(&mut trie.nodes, parent_key, key, Direction::Prefix);
                    }
                }
            }

            // Generate SUFFIXES: start varies from 0 to read_len-min_len, end=read_len
            // Skip suffix that equals the full read (already added as prefix)
            for start in 1..=read_len.saturating_sub(min_len) {
                let sequence = &read[start..read_len];

                // Track that this read has this affix
                trie.affix_contributors
                    .entry(sequence.to_string())
                    .or_default()
                    .insert(read_idx);

                // Canonicalize: reuse existing key if sequence seen before
                let key = if let Some(&existing_key) = trie.sequence_to_key.get(sequence) {
                    existing_key
                } else {
                    let new_key = (read_idx, start, read_len);
                    trie.nodes.insert(new_key, AffixNode::Leaf);
                    trie.sequence_to_key.insert(sequence.to_string(), new_key);

                    // Track k-mers (minimum length affixes)
                    if read_len - start == min_len {
                        trie.roots.push(new_key);
                    }

                    new_key
                };

                // Phase 1b: Link to SUFFIX parent (longer suffix = one char added left/start)
                if start > 1 {
                    let parent_seq = &read[start - 1..read_len];
                    if let Some(&parent_key) = trie.sequence_to_key.get(parent_seq) {
                        Self::insert_extension(&mut trie.nodes, parent_key, key, Direction::Suffix);
                    }
                }
            }

            // Phase 1c: Mark full read
            let full_key = (read_idx, 0, read_len);
            Self::mark_full_read(&mut trie.nodes, full_key, read_idx);
        }

        // Phase 2: Prune nodes without cross-read connections
        // Retain:
        // - Delta nodes (branching points and full reads)
        // - Leaf nodes that appear in multiple reads (cross-read overlaps)
        // - Chains with cross-read links

        trie.nodes.retain(|&key, node| match node {
            AffixNode::Delta { .. } => true, // Always keep branching points and full reads (containment)
            AffixNode::Leaf => {
                // Keep if it appears in multiple reads (cross-read overlap)
                let affix_str = &reads[key.0][key.1..key.2];
                trie.affix_contributors
                    .get(affix_str)
                    .map_or(false, |set| set.len() > 1)
            }
            AffixNode::Chain { .. } => Self::has_cross_read_link(key, node),
        });

        trie
    }

    /// Insert an extension, progressively promoting nodes as branching is detected.
    fn insert_extension(
        nodes: &mut HashMap<AffixSlice, AffixNode>,
        parent_key: AffixSlice,
        child_key: AffixSlice,
        direction: Direction,
    ) {
        let parent_node = nodes.get_mut(&parent_key).expect("Parent node must exist");

        match parent_node {
            // Leaf becomes Chain
            AffixNode::Leaf => {
                *parent_node = AffixNode::Chain {
                    prefix_link: if direction == Direction::Prefix {
                        Some(child_key)
                    } else {
                        None
                    },
                    suffix_link: if direction == Direction::Suffix {
                        Some(child_key)
                    } else {
                        None
                    },
                };
            }

            // Chain stays Chain or promotes to Delta
            AffixNode::Chain {
                prefix_link,
                suffix_link,
            } => {
                match direction {
                    Direction::Prefix => {
                        if prefix_link.is_none() {
                            *prefix_link = Some(child_key);
                        } else if *prefix_link != Some(child_key) {
                            // BRANCHING! Promote to Delta
                            let old_prefix = prefix_link.unwrap();
                            let old_suffix = *suffix_link;
                            *parent_node = AffixNode::Delta {
                                prefix_extensions: vec![old_prefix, child_key],
                                suffix_extensions: old_suffix.into_iter().collect(),
                                roots: RootsList::new(),
                            };
                        }
                    }
                    Direction::Suffix => {
                        if suffix_link.is_none() {
                            *suffix_link = Some(child_key);
                        } else if *suffix_link != Some(child_key) {
                            // BRANCHING! Promote to Delta
                            let old_suffix = suffix_link.unwrap();
                            let old_prefix = *prefix_link;
                            *parent_node = AffixNode::Delta {
                                prefix_extensions: old_prefix.into_iter().collect(),
                                suffix_extensions: vec![old_suffix, child_key],
                                roots: RootsList::new(),
                            };
                        }
                    }
                }
            }

            // Delta just adds unique extension
            AffixNode::Delta {
                prefix_extensions,
                suffix_extensions,
                ..
            } => match direction {
                Direction::Prefix => {
                    if !prefix_extensions.contains(&child_key) {
                        prefix_extensions.push(child_key);
                    }
                }
                Direction::Suffix => {
                    if !suffix_extensions.contains(&child_key) {
                        suffix_extensions.push(child_key);
                    }
                }
            },
        }
    }

    /// Mark a complete read at this affix key.
    fn mark_full_read(
        nodes: &mut HashMap<AffixSlice, AffixNode>,
        full_key: AffixSlice,
        read_idx: usize,
    ) {
        let node = nodes.get_mut(&full_key).expect("Full read key must exist");

        match node {
            // Isolated full read - promote to Delta to track containment
            AffixNode::Leaf => {
                let mut roots = RootsList::new();
                roots.add(read_idx);
                *node = AffixNode::Delta {
                    prefix_extensions: vec![],
                    suffix_extensions: vec![],
                    roots,
                };
            }

            // Full read that's also a chain - promote to Delta with roots
            AffixNode::Chain {
                prefix_link,
                suffix_link,
            } => {
                let old_prefix = *prefix_link;
                let old_suffix = *suffix_link;
                let mut roots = RootsList::new();
                roots.add(read_idx);
                *node = AffixNode::Delta {
                    prefix_extensions: old_prefix.into_iter().collect(),
                    suffix_extensions: old_suffix.into_iter().collect(),
                    roots,
                };
            }

            // Already a Delta - just add to roots
            AffixNode::Delta { roots, .. } => {
                roots.add(read_idx);
            }
        }
    }

    /// Check if a chain node links to a different read.
    fn has_cross_read_link(key: AffixSlice, node: &AffixNode) -> bool {
        let own_read = key.0;

        if let AffixNode::Chain {
            prefix_link,
            suffix_link,
        } = node
        {
            if let Some(prefix) = prefix_link {
                if prefix.0 != own_read {
                    return true;
                }
            }
            if let Some(suffix) = suffix_link {
                if suffix.0 != own_read {
                    return true;
                }
            }
        }

        false
    }

    /// Generate all overlap candidates from the trie.
    /// Returns a vector of (suffix_idx, prefix_idx, affix_key) tuples.
    ///
    /// Overlaps are detected at Delta nodes where different reads meet,
    /// and at FullRead nodes that are referenced by multiple reads.
    pub fn overlap_candidates(&self, reads: &[String]) -> Vec<(usize, usize, String)> {
        let mut candidates = Vec::new();

        for (&key, node) in &self.nodes {
            match node {
                // Deltas with roots = full reads overlapping with extensions
                AffixNode::Delta {
                    prefix_extensions,
                    suffix_extensions,
                    roots,
                } => {
                    // Collect all reads at this delta
                    let mut reads_here = HashSet::new();
                    reads_here.insert(key.0);

                    for &ext in prefix_extensions {
                        reads_here.insert(ext.0);
                    }
                    for &ext in suffix_extensions {
                        reads_here.insert(ext.0);
                    }
                    for root in roots.iter() {
                        reads_here.insert(root);
                    }

                    // All pairs of different reads = overlap at this affix
                    let reads_vec: Vec<usize> = reads_here.into_iter().collect();
                    for &read1 in &reads_vec {
                        for &read2 in &reads_vec {
                            if read1 != read2 {
                                // Resolve the affix string from the key
                                let affix_str = reads[key.0][key.1..key.2].to_string();
                                candidates.push((read1, read2, affix_str));
                            }
                        }
                    }
                }

                // Leaf nodes: check for canonicalized affixes present in multiple reads
                AffixNode::Leaf => {
                    let affix_str = &reads[key.0][key.1..key.2];
                    if let Some(contributors) = self.affix_contributors.get(affix_str) {
                        if contributors.len() > 1 {
                            let ids: Vec<usize> = contributors.iter().copied().collect();
                            for &read1 in &ids {
                                for &read2 in &ids {
                                    if read1 != read2 {
                                        candidates.push((read1, read2, affix_str.to_string()));
                                    }
                                }
                            }
                        }
                    }
                }

                // Chains don't generate overlaps (no branching)
                _ => {}
            }
        }

        candidates
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_trie_for_simple_reads() {
        let reads = vec!["ACGT".to_string(), "GTAA".to_string()];
        let trie = PrunedAffixTrie::build(&reads, 2);

        // Should have nodes for overlapping affixes
        assert!(!trie.nodes.is_empty());

        // Should find the GT overlap
        let candidates = trie.overlap_candidates(&reads);
        let has_overlap = candidates
            .iter()
            .any(|(src, dst, affix)| *src != *dst && affix.contains("GT"));
        assert!(has_overlap, "Should find overlap between ACGT and GTAA");
    }

    #[test]
    fn prunes_nodes_without_overlaps() {
        let reads = vec!["AAAA".to_string(), "TTTT".to_string()];
        let trie = PrunedAffixTrie::build(&reads, 2);

        // Should have few or no nodes since no overlaps exist
        let candidates = trie.overlap_candidates(&reads);
        assert!(
            candidates.is_empty(),
            "Should have no overlap candidates for non-overlapping reads"
        );
    }

    #[test]
    fn handles_empty_reads() {
        let reads: Vec<String> = vec![];
        let trie = PrunedAffixTrie::build(&reads, 3);

        assert!(trie.nodes.is_empty());
        assert!(trie.roots.is_empty());
        assert!(trie.overlap_candidates(&reads).is_empty());
    }

    #[test]
    fn handles_contained_reads() {
        // read1 is a prefix of read2 (containment via prefix)
        let reads = vec!["ACGTACGT".to_string(), "ACGTACGTAA".to_string()];
        let trie = PrunedAffixTrie::build(&reads, 3);

        // Should find overlaps
        let candidates = trie.overlap_candidates(&reads);
        assert!(
            !candidates.is_empty(),
            "Should find overlaps for contained reads"
        );

        // Should detect that read1 is a full read (contained in read2 as a prefix)
        let full_read_key = (0, 0, 8); // read1[0:8]
        assert!(
            matches!(trie.nodes.get(&full_read_key), Some(AffixNode::Delta { roots, .. }) if roots.contains(0)),
            "read1 should be marked as a full read in a Delta node"
        );
    }

    #[test]
    fn canonicalizes_duplicate_sequences() {
        // Two reads where a suffix of read1 equals a prefix of read2
        // read1 suffix "ACG" matches read2 prefix "ACG"
        let reads = vec!["TTACG".to_string(), "ACGAA".to_string()];
        let trie = PrunedAffixTrie::build(&reads, 3);

        // The shared "ACG" affix should be canonicalized
        let acg_key = trie.sequence_to_key.get("ACG");
        assert!(
            acg_key.is_some(),
            "Shared affix should be in canonicalization map"
        );

        // Should detect that ACG appears in both reads
        assert!(
            trie.affix_contributors
                .get("ACG")
                .map_or(false, |c| c.len() > 1),
            "ACG should have multiple contributors"
        );

        // Should create overlaps at the shared affix
        let candidates = trie.overlap_candidates(&reads);
        let has_acg_overlap = candidates
            .iter()
            .any(|(r1, r2, affix)| *r1 != *r2 && affix == "ACG");
        assert!(has_acg_overlap, "Should find overlap at shared affix ACG");
    }
}
