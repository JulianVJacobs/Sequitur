//! Overlap graph construction primitives translated from the Python reference.

use std::collections::HashMap;

use sprs::CsMat;
use strsim::damerau_levenshtein;

use crate::affix::{AffixKind, AffixMap, AffixStructure};

/// Mapping from source read index to weighted successor indices.
pub type Adjacency = Vec<HashMap<usize, usize>>;

/// Companion mapping that records the overlap span associated with each edge.
pub type OverlapLengths = Vec<HashMap<usize, usize>>;

/// Configuration options that govern overlap graph construction.
#[derive(Debug, Clone, Copy)]
pub struct OverlapConfig {
    /// Maximum normalised Damerau–Levenshtein distance permitted for an edge.
    pub max_diff: f32,
    /// Minimum suffix length used when generating candidate overlaps.
    pub min_suffix_len: usize,
    /// Enable a threaded implementation (currently downgraded to sequential).
    pub use_threads: bool,
    /// Desired worker count when threading is enabled.
    pub max_workers: usize,
    /// Use pruned trie implementation instead of affix array (default: true).
    pub use_trie: bool,
    /// Enable parabolic early-out optimization for overlap detection.
    pub use_parabolic_cutoff: bool,
    /// Patience for parabolic early-out (consecutive increases before cutoff).
    pub parabolic_patience: usize,
}

impl Default for OverlapConfig {
    fn default() -> Self {
        Self {
            max_diff: 0.25,
            min_suffix_len: crate::affix::DEFAULT_MIN_SUFFIX_LEN,
            use_threads: false,
            max_workers: 1,
            use_trie: true, // Trie is the default (optimized implementation)
            use_parabolic_cutoff: true,
            parabolic_patience: 3,
        }
    }
}

/// Compute the normalised Damerau–Levenshtein distance and overlap weight.
pub fn normalised_damerau_levenshtein_distance(
    suffix: &str,
    prefix: &str,
) -> Option<(f32, f32, usize)> {
    let overlap_len = suffix.len().min(prefix.len());
    if overlap_len == 0 {
        return None;
    }

    let suffix_window = &suffix[suffix.len() - overlap_len..];
    let prefix_window = &prefix[..overlap_len];
    let distance = damerau_levenshtein(suffix_window, prefix_window);

    if distance > overlap_len {
        return None;
    }

    let diff = distance as f32 / overlap_len as f32;
    let score = 1.0 - diff;
    Some((diff, score, overlap_len))
}

/// Build the weighted overlap graph for the supplied reads using the unified affix interface.
/// Returns sparse matrices directly to avoid intermediate HashMap allocation.
/// Supports both array and trie implementations based on config.use_trie.
pub fn create_overlap_graph_unified<'a>(
    reads: &'a [String],
    config: OverlapConfig,
) -> (CsMat<usize>, CsMat<usize>) {
    let min_suffix_len = config.min_suffix_len.max(1);
    let affix_structure = AffixStructure::build(reads, min_suffix_len, config.use_trie);

    match affix_structure {
        AffixStructure::Array(affix_map) => {
            let (_map, adj, ovl) = create_overlap_graph_from_array(reads, affix_map, config);
            (adj, ovl)
        }
        AffixStructure::Trie(affix_trie) => {
            create_overlap_graph_from_trie(reads, &affix_trie, config)
        }
    }
}

/// Build overlap graph from trie implementation (optimized).
fn create_overlap_graph_from_trie(
    reads: &[String],
    affix_trie: &crate::affix::PrunedAffixTrie,
    config: OverlapConfig,
) -> (CsMat<usize>, CsMat<usize>) {
    if reads.is_empty() {
        let empty1 = CsMat::<usize>::zero((0, 0));
        let empty2 = CsMat::<usize>::zero((0, 0));
        return (empty1, empty2);
    }

    let n = reads.len();
    let candidates = affix_trie.overlap_candidates(reads);

    // Build best scores and spans for each read pair
    let mut best_scores: HashMap<(usize, usize), f32> = HashMap::new();
    let mut best_spans: HashMap<(usize, usize), usize> = HashMap::new();

    for (suffix_idx, prefix_idx, shared_affix) in &candidates {
        let suffix_read = &reads[*suffix_idx];
        let prefix_read = &reads[*prefix_idx];

        // Verify the overlap with edit distance
        if let Some((score, overlap_len)) =
            verify_overlap_with_parabolic(suffix_read, prefix_read, shared_affix, config)
        {
            let key = (*suffix_idx, *prefix_idx);
            let prev_score = best_scores.get(&key).copied().unwrap_or(f32::MIN);
            let prev_span = best_spans.get(&key).copied().unwrap_or(0);

            if score > prev_score || (score == prev_score && overlap_len > prev_span) {
                best_scores.insert(key, score);
                best_spans.insert(key, overlap_len);
            }
        }
    }

    // Convert to CSR sparse matrices
    let mut a_indptr: Vec<usize> = Vec::with_capacity(n + 1);
    let mut a_indices: Vec<usize> = Vec::new();
    let mut a_data: Vec<usize> = Vec::new();
    let mut o_indptr: Vec<usize> = Vec::with_capacity(n + 1);
    let mut o_indices: Vec<usize> = Vec::new();
    let mut o_data: Vec<usize> = Vec::new();

    a_indptr.push(0);
    o_indptr.push(0);

    for src_idx in 0..n {
        let mut edges: Vec<(usize, f32)> = best_scores
            .iter()
            .filter_map(|((src, dst), &score)| {
                if *src == src_idx {
                    Some((*dst, score))
                } else {
                    None
                }
            })
            .collect();
        edges.sort_by_key(|(dst, _)| *dst);

        for (dst, score) in edges {
            let overlap_len = best_spans.get(&(src_idx, dst)).copied().unwrap_or(0);
            a_indices.push(dst);
            a_data.push(score as usize);
            o_indices.push(dst);
            o_data.push(overlap_len);
        }

        a_indptr.push(a_indices.len());
        o_indptr.push(o_indices.len());
    }

    let adjacency = CsMat::new((n, n), a_indptr, a_indices, a_data);
    let overlaps = CsMat::new((n, n), o_indptr, o_indices, o_data);

    (adjacency, overlaps)
}

/// Verify overlap between suffix and prefix with parabolic early-out optimization.
fn verify_overlap_with_parabolic(
    suffix_read: &str,
    prefix_read: &str,
    _shared_affix: &str,
    config: OverlapConfig,
) -> Option<(f32, usize)> {
    let min_len = config.min_suffix_len.max(1);
    let max_span = suffix_read.len().min(prefix_read.len());

    let mut best_score = f32::MIN;
    let mut best_overlap = 0;

    let mut min_diff_seen = f32::MAX;
    let mut increasing_count = 0;

    for span in min_len..=max_span {
        if span > suffix_read.len() || span > prefix_read.len() {
            break;
        }

        let suffix_window = &suffix_read[suffix_read.len() - span..];
        let prefix_window = &prefix_read[..span];

        let distance = damerau_levenshtein(suffix_window, prefix_window);
        let float_diff = distance as f32 / span as f32;

        // Parabolic early-out logic
        if config.use_parabolic_cutoff && span > min_len + 5 {
            if float_diff < min_diff_seen {
                min_diff_seen = float_diff;
                increasing_count = 0;
            } else if float_diff > config.max_diff {
                increasing_count += 1;
                if increasing_count >= config.parabolic_patience {
                    break; // Passed the valley, no point continuing
                }
            }
        }

        if float_diff < config.max_diff {
            let score = (span as i32 - distance as i32) as f32;
            if score > best_score || (score == best_score && span > best_overlap) {
                best_score = score;
                best_overlap = span;
            }
        }
    }

    if best_score > f32::MIN {
        Some((best_score, best_overlap))
    } else {
        None
    }
}

/// Build the weighted overlap graph for the supplied reads (legacy array-based).
/// Returns sparse matrices directly to avoid intermediate HashMap allocation.
pub fn create_overlap_graph<'a>(
    reads: &'a [String],
    affix_map: Option<AffixMap<'a>>,
    config: OverlapConfig,
) -> (AffixMap<'a>, CsMat<usize>, CsMat<usize>) {
    let (affix_map, adj, ovl) = create_overlap_graph_from_array(
        reads,
        affix_map.unwrap_or_else(|| AffixMap::build(reads, config.min_suffix_len.max(1))),
        config,
    );
    (affix_map, adj, ovl)
}

/// Build overlap graph from array implementation (original algorithm).
fn create_overlap_graph_from_array<'a>(
    reads: &'a [String],
    affix_map: AffixMap<'a>,
    config: OverlapConfig,
) -> (AffixMap<'a>, CsMat<usize>, CsMat<usize>) {
    #[cfg(not(feature = "parallel"))]
    if config.use_threads {
        use log::warn;
        warn!(
            "Threaded overlap construction requested, but the 'parallel' feature is not enabled; falling back to sequential mode"
        );
    }

    let min_suffix_len = config.min_suffix_len.max(1);

    if reads.is_empty() {
        let empty1 = CsMat::<usize>::zero((0, 0));
        let empty2 = CsMat::<usize>::zero((0, 0));
        return (affix_map, empty1, empty2);
    }

    let n = reads.len();
    let mut a_indptr: Vec<usize> = Vec::with_capacity(n + 1);
    let mut a_indices: Vec<usize> = Vec::new();
    let mut a_data: Vec<usize> = Vec::new();
    let mut o_indptr: Vec<usize> = Vec::with_capacity(n + 1);
    let mut o_indices: Vec<usize> = Vec::new();
    let mut o_data: Vec<usize> = Vec::new();

    a_indptr.push(0);
    o_indptr.push(0);

    #[cfg(feature = "parallel")]
    use rayon::prelude::*;

    let process_read = |(read_idx, read): (usize, &String)| {
        if read.len() < min_suffix_len {
            return (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        }
        let read_len = read.len();
        let mut best_scores: HashMap<usize, f32> = HashMap::new();
        let mut best_spans: HashMap<usize, usize> = HashMap::new();

        for span in min_suffix_len..read_len {
            let suffix_seq = &read[read_len - span..];
            for key in affix_map.keys.iter() {
                if key.kind != AffixKind::Prefix {
                    continue;
                }
                static EMPTY: [usize; 0] = [];
                let prefix_read_indices = affix_map
                    .reads_for(key)
                    .map(|v| v.as_slice())
                    .unwrap_or(&EMPTY);
                if prefix_read_indices.contains(&read_idx) {
                    continue;
                }
                let prefix_seq = key.affix;
                let overlap_len = suffix_seq.len().min(prefix_seq.len());
                if overlap_len == 0 {
                    continue;
                }
                let suffix_overlap = &suffix_seq[suffix_seq.len() - overlap_len..];
                let prefix_overlap = &prefix_seq[..overlap_len];
                let int_edit_distance = damerau_levenshtein(suffix_overlap, prefix_overlap);
                let float_diff = int_edit_distance as f32 / overlap_len as f32;
                if float_diff >= config.max_diff {
                    continue;
                }
                let score = overlap_len as i32 - int_edit_distance as i32;
                for &dst in prefix_read_indices {
                    let prev_score = best_scores.get(&dst).copied().unwrap_or(f32::MIN);
                    let prev_span = best_spans.get(&dst).copied().unwrap_or(0);
                    if score as f32 > prev_score
                        || (score as f32 == prev_score && overlap_len > prev_span)
                    {
                        best_scores.insert(dst, score as f32);
                        best_spans.insert(dst, overlap_len);
                    }
                }
            }
        }

        let mut a_idx = Vec::new();
        let mut a_dat = Vec::new();
        let mut o_idx = Vec::new();
        let mut o_dat = Vec::new();
        let mut pairs: Vec<(usize, f32)> = best_scores.into_iter().collect();
        pairs.sort_by_key(|(dst, _)| *dst);
        for (dst, score) in pairs.into_iter() {
            let overlap_len = best_spans.get(&dst).copied().unwrap_or(0);
            a_idx.push(dst);
            a_dat.push(score as usize);
            o_idx.push(dst);
            o_dat.push(overlap_len);
        }
        (a_idx, a_dat, o_idx, o_dat)
    };

    let results: Vec<(Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>)> = if config.use_threads {
        #[cfg(feature = "parallel")]
        {
            use rayon::ThreadPoolBuilder;
            let pool = ThreadPoolBuilder::new()
                .num_threads(config.max_workers)
                .build()
                .expect("Failed to build rayon thread pool");
            pool.install(|| reads.par_iter().enumerate().map(process_read).collect())
        }
        #[cfg(not(feature = "parallel"))]
        {
            reads.iter().enumerate().map(process_read).collect()
        }
    } else {
        reads.iter().enumerate().map(process_read).collect()
    };

    for (a_idx, a_dat, o_idx, o_dat) in results {
        for idx in a_idx {
            a_indices.push(idx);
        }
        for dat in a_dat {
            a_data.push(dat);
        }
        for idx in o_idx {
            o_indices.push(idx);
        }
        for dat in o_dat {
            o_data.push(dat);
        }
        a_indptr.push(a_indices.len());
        o_indptr.push(o_indices.len());
    }

    let adjacency = CsMat::new((n, n), a_indptr, a_indices, a_data);
    let overlaps = CsMat::new((n, n), o_indptr, o_indices, o_data);

    (affix_map, adjacency, overlaps)
}
/// Compute normalised per-edge confidences (softmax over outgoing weights).
pub fn compute_edge_confidences(adjacency: &Adjacency) -> Vec<Vec<(usize, f64)>> {
    adjacency
        .iter()
        .map(|edges| {
            if edges.is_empty() {
                return Vec::new();
            }

            let max_weight = edges.values().copied().max().unwrap_or(0) as f64;
            let mut ordered: Vec<(usize, usize)> =
                edges.iter().map(|(&dst, &weight)| (dst, weight)).collect();
            ordered.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

            let exps: Vec<f64> = ordered
                .iter()
                .map(|(_, weight)| ((*weight as f64) - max_weight).exp())
                .collect();
            let sum: f64 = exps.iter().sum();

            ordered
                .into_iter()
                .zip(exps.into_iter())
                .map(|((dst, _weight), value)| (dst, if sum == 0.0 { 0.0 } else { value / sum }))
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn computes_normalised_distance() {
        let (diff, score, overlap) =
            normalised_damerau_levenshtein_distance("ACGT", "ACGA").expect("diff");
        assert!((diff - 0.25).abs() < f32::EPSILON);
        assert_eq!(score, 1.0 - diff);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn builds_overlap_graph_for_simple_pair() {
        let reads = vec!["ACGT".to_string(), "GTAA".to_string()];
        let config = OverlapConfig {
            max_diff: 0.5,
            min_suffix_len: 2,
            ..Default::default()
        };

        let (_affixes, adjacency, overlaps) = create_overlap_graph(&reads, None, config);

        let a_csc = adjacency.to_csc();
        let o_csc = overlaps.to_csc();

        assert_eq!(a_csc.rows(), 2);
        assert_eq!(a_csc.cols(), 2);
        assert_eq!(o_csc.rows(), 2);
        assert_eq!(o_csc.cols(), 2);

        let edge = a_csc.get(0, 1).copied();
        let span = o_csc.get(0, 1).copied();

        assert_eq!(edge, Some(2));
        assert_eq!(span, Some(2));
    }

    #[test]
    fn computes_edge_confidences_per_source() {
        let mut adjacency: Adjacency = vec![HashMap::new(); 3];
        adjacency[0].insert(1, 3);
        adjacency[0].insert(2, 1);
        adjacency[1].insert(2, 4);

        let confidences = compute_edge_confidences(&adjacency);
        assert_eq!(confidences.len(), 3);

        let row0_sum: f64 = confidences[0].iter().map(|(_, p)| *p).sum();
        assert!((row0_sum - 1.0).abs() < 1e-6);
        assert!(confidences[0][0].1 > confidences[0][1].1);

        let row1_sum: f64 = confidences[1].iter().map(|(_, p)| *p).sum();
        assert!((row1_sum - 1.0).abs() < 1e-6);
        assert_eq!(confidences[1][0].0, 2);
        assert!((confidences[1][0].1 - 1.0).abs() < 1e-12);

        assert!(confidences[2].is_empty());
    }
}
