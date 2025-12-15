//! Overlap graph construction primitives translated from the Python reference.

use std::collections::HashMap;

use sprs::CsMat;
use strsim::damerau_levenshtein;

use crate::affix::PrunedAffixTrie;

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
    /// Exponent for error penalty in quality-adjusted scoring (1.0=linear, 2.0=quadratic).
    pub error_penalty_exponent: f32,
    /// Detect and remove low-quality overlaps using knee-point detection on score distribution.
    pub detect_score_cliff: bool,
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
            error_penalty_exponent: 2.0, // Quadratic penalty by default
            detect_score_cliff: true,    // Knee detection enabled by default
        }
    }
}

/// Detect the knee (elbow) point in a sorted score distribution.
/// Returns the index where scores drop dramatically, using the distance-from-line method.
/// Scores should be sorted in descending order.
fn detect_knee_point(scores: &[f32]) -> usize {
    if scores.is_empty() {
        return 0;
    }

    if scores.len() == 1 {
        return 0; // Only one score, no knee
    }

    let n = scores.len();
    let start_score = scores[0];
    let end_score = scores[n - 1];

    // Handle edge case where all scores are the same
    if (start_score - end_score).abs() < 1e-6 {
        return n - 1; // Linear, keep all points
    }

    // For just 2 points, check if drop is dramatic (>50% reduction)
    if n == 2 {
        let ratio = scores[1] / scores[0];
        if ratio < 0.5 {
            return 0; // Keep only the first (highest) score
        } else {
            return 1; // Keep both
        }
    }

    // Find point with maximum distance from line connecting first and last point
    let mut max_distance = 0.0f32;
    let mut knee_idx = n - 1;

    for i in 1..n - 1 {
        let score = scores[i];
        // Point on line at rank i (linear interpolation)
        let line_score = start_score + (end_score - start_score) * (i as f32 / (n - 1) as f32);
        // Distance from actual score to line
        let distance = (score - line_score).abs();
        if distance > max_distance {
            max_distance = distance;
            knee_idx = i;
        }
    }

    knee_idx
}

/// Compute the normalised Damerau–Levenshtein distance and overlap weight.
pub fn normalised_damerau_levenshtein_distance(
    suffix: &str,
    prefix: &str,
) -> Option<(f32, f32, usize)> {
    let overlap_len = suffix.len().min(prefix.len());
    let longer_len = suffix.len().max(prefix.len());
    if overlap_len == 0 {
        return None;
    }

    let suffix_window = &suffix[suffix.len() - overlap_len..];
    let prefix_window = &prefix[..overlap_len];
    let distance = damerau_levenshtein(suffix_window, prefix_window);

    if distance > overlap_len {
        return None;
    }

    // Normalize by longer read length to penalize insertions/deletions
    let diff = distance as f32 / longer_len as f32;
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
    let trie = PrunedAffixTrie::build(reads, min_suffix_len, config.max_diff);
    create_overlap_graph_from_trie(reads, &trie, config)
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
    // Trie now includes fuzzy matches in extension vectors automatically
    log::info!("Collecting & verifying candidates from trie...");
    // Use the trie-side verified candidate emitter (simple verifier)
    let verified = affix_trie.overlap_verified_candidates_simple(
        reads,
        config.max_diff,
        config.min_suffix_len,
    );
    log::info!("Found {} verified edges", verified.len());

    // Debug: log all verified edges
    for (suffix_idx, prefix_idx, score, span) in &verified {
        log::debug!(
            "[VERIFIED] {} -> {}: score={:.2}, span={} ",
            suffix_idx,
            prefix_idx,
            score,
            span
        );
    }

    // Build from verified tuples

    // Build per-source rows from verified tuples
    #[cfg(feature = "parallel")]
    let mut per_source: Vec<HashMap<usize, (f32, usize)>> =
        (0..n).map(|_| HashMap::new()).collect();

    #[cfg(feature = "parallel")]
    {
        for (s, d, score, span) in &verified {
            let row = &mut per_source[*s];
            let prev = row.get(d).copied().unwrap_or((f32::MIN, 0));
            if *score > prev.0 || (*score == prev.0 && *span > prev.1) {
                row.insert(*d, (*score, *span));
            }
        }
    }

    #[cfg(not(feature = "parallel"))]
    {
        // Streaming emission from verified tuples grouped by source
        log::info!(
            "Streaming CSR emission from {} verified edges...",
            verified.len()
        );

        // Prepare CSR buffers and build directly from grouped verified tuples
        let mut a_indptr: Vec<usize> = Vec::with_capacity(n + 1);
        let mut a_indices: Vec<usize> = Vec::new();
        let mut a_data: Vec<usize> = Vec::new();
        let mut o_indptr: Vec<usize> = Vec::with_capacity(n + 1);
        let mut o_indices: Vec<usize> = Vec::new();
        let mut o_data: Vec<usize> = Vec::new();

        a_indptr.push(0);
        o_indptr.push(0);

        // Sort verified tuples by source
        let mut verified_sorted = verified.clone();
        verified_sorted.sort_by_key(|(s, _d, _score, _span)| *s);

        let mut i = 0usize;
        while i < verified_sorted.len() {
            let src = verified_sorted[i].0;
            // collect best per-dst for this src
            let mut row_map: HashMap<usize, (f32, usize)> = HashMap::new();
            while i < verified_sorted.len() && verified_sorted[i].0 == src {
                let (s, d, score, span) = verified_sorted[i];
                let prev = row_map.get(&d).copied().unwrap_or((f32::MIN, 0));
                if score > prev.0 || (score == prev.0 && span > prev.1) {
                    row_map.insert(d, (score, span));
                }
                i += 1;
            }

            // Convert to edge list and apply knee-detection
            let mut edges: Vec<(usize, f32)> =
                row_map.iter().map(|(&dst, &(sco, _))| (dst, sco)).collect();
            let edges_before = edges.len();
            if config.detect_score_cliff && !edges.is_empty() {
                edges.sort_by(|(_, s1), (_, s2)| s2.partial_cmp(s1).unwrap());
                let scores: Vec<f32> = edges.iter().map(|(_, s)| *s).collect();
                let knee_idx = detect_knee_point(&scores);
                edges.truncate(knee_idx + 1);
            }
            edges.sort_by_key(|(dst, _)| *dst);

            for (dst, score) in edges {
                let span = row_map.get(&dst).copied().map(|(_s, sp)| sp).unwrap_or(0);
                a_indices.push(dst);
                a_data.push(score as usize);
                o_indices.push(dst);
                o_data.push(span);
            }

            a_indptr.push(a_indices.len());
            o_indptr.push(o_indices.len());
        }

        let adjacency = CsMat::new((n, n), a_indptr, a_indices, a_data);
        let overlaps = CsMat::new((n, n), o_indptr, o_indices, o_data);

        return (adjacency, overlaps);
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
        // Collect outgoing edges from the per-row map
        let mut edges: Vec<(usize, f32)> = per_source[src_idx]
            .iter()
            .map(|(&dst, &(score, _span))| (dst, score))
            .collect();

        let edges_before = edges.len();

        // Apply knee-point detection if configured
        if config.detect_score_cliff && !edges.is_empty() {
            // Sort by score descending to detect the cliff
            edges.sort_by(|(_, s1), (_, s2)| s2.partial_cmp(s1).unwrap());
            let scores: Vec<f32> = edges.iter().map(|(_, s)| *s).collect();
            let knee_idx = detect_knee_point(&scores);

            log::debug!(
                "[SCORE_CLIFF] Row {}: {} edges before, knee at index {}, keeping {} edges",
                src_idx,
                edges_before,
                knee_idx,
                knee_idx + 1
            );
            log::debug!("[SCORE_CLIFF] Row {} scores: {:?}", src_idx, scores);

            // Truncate to keep only overlaps up to and including the knee point
            edges.truncate(knee_idx + 1);

            if edges.len() < edges_before {
                log::debug!(
                    "[SCORE_CLIFF] Row {}: removed {} edges due to score cliff",
                    src_idx,
                    edges_before - edges.len()
                );
            }
        }

        edges.sort_by_key(|(dst, _)| *dst);

        for (dst, score) in edges {
            let overlap_len = per_source[src_idx]
                .get(&dst)
                .copied()
                .map(|(_s, span)| span)
                .unwrap_or(0);
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

/// Verify overlap from trie anchor by extending bidirectionally from exact match.
/// This leverages the shared_affix hint to avoid O(m²) full scan.
/// Public for benchmarking purposes.
pub fn verify_overlap_from_anchor(
    suffix_read: &str,
    prefix_read: &str,
    shared_affix: &str,
    config: OverlapConfig,
) -> Option<(f32, usize)> {
    let anchor_len = shared_affix.len();
    if anchor_len == 0 {
        // No anchor, do simple full scan
        return verify_overlap_simple(suffix_read, prefix_read, config);
    }

    let min_len = config.min_suffix_len.max(1);
    let max_span = suffix_read.len().min(prefix_read.len());

    let mut best_score = f32::MIN;
    let mut best_overlap = 0;

    // Start search from anchor and extend outward
    let start_span = anchor_len.max(min_len);
    let end_span = max_span.min(start_span + 20); // Only check ±20bp around anchor

    for span in start_span..=end_span {
        if span > suffix_read.len() || span > prefix_read.len() {
            break;
        }

        let suffix_window = &suffix_read[suffix_read.len() - span..];
        let prefix_window = &prefix_read[..span];

        let distance = damerau_levenshtein(suffix_window, prefix_window);
        let float_diff = distance as f32 / span as f32;

        if float_diff <= config.max_diff {
            // Quality-adjusted scoring: penalize errors based on error rate
            let error_rate = distance as f32 / span as f32;
            let quality_factor = (1.0 - error_rate).powf(config.error_penalty_exponent);
            let score = span as f32 * quality_factor;
            if score > best_score || (score == best_score && span > best_overlap) {
                best_score = score;
                best_overlap = span;
            }
        } else if span > anchor_len + 5 {
            // Error rate increasing beyond anchor, stop extending
            break;
        }
    }

    if best_score > f32::MIN {
        Some((best_score, best_overlap))
    } else {
        None
    }
}

/// Simple full-scan verification without optimizations.
/// Used as fallback when no anchor hint exists.
fn verify_overlap_simple(
    suffix_read: &str,
    prefix_read: &str,
    config: OverlapConfig,
) -> Option<(f32, usize)> {
    let min_len = config.min_suffix_len.max(1);
    let max_span = suffix_read.len().min(prefix_read.len());

    let mut best_score = f32::MIN;
    let mut best_overlap = 0;

    for span in min_len..=max_span {
        if span > suffix_read.len() || span > prefix_read.len() {
            break;
        }

        let suffix_window = &suffix_read[suffix_read.len() - span..];
        let prefix_window = &prefix_read[..span];

        let distance = damerau_levenshtein(suffix_window, prefix_window);
        let float_diff = distance as f32 / span as f32;

        if float_diff <= config.max_diff {
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
    _affix_map: Option<()>,
    config: OverlapConfig,
) -> (PrunedAffixTrie, CsMat<usize>, CsMat<usize>) {
    let trie = PrunedAffixTrie::build(reads, config.min_suffix_len.max(1), config.max_diff);
    let (adj, ovl) = create_overlap_graph_from_trie(reads, &trie, config);
    (trie, adj, ovl)
}

/// Build overlap graph from array implementation (original algorithm).
// The original array-based implementation was removed. Use the trie-based
// construction path implemented above.
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
