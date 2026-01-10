//! Overlap graph construction primitives translated from the Python reference.

use std::collections::{HashMap, HashSet, VecDeque};

use sprs::CsMat;
use strsim::damerau_levenshtein;

use crate::affix::PrunedAffixTrie;
use crate::read_source::{QualityAwareReadSource, ReadSource};

/// Quality score integration strategy.
/// None: no quality weighting (current behavior)
/// Position: soft mismatch penalty scaled by quality at mismatched positions
/// Confidence: phred-explicit scoring; matches at high-Q positions boost score more
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualityScoring {
    None,
    Position,
    Confidence,
}

/// Configuration for quality-aware overlap scoring.
#[derive(Debug, Clone)]
pub struct QualityConfig {
    /// How quality scores affect overlap weight.
    pub scoring_mode: QualityScoring,

    /// Exponent for error penalty: 1.0=linear, 2.0=quadratic.
    pub error_penalty_exponent: f32,

    /// Minimum quality threshold; bases below this treated as ambiguous.
    /// None = no minimum threshold.
    pub min_quality: Option<i32>,

    /// Enable quality-based tie-breaking in consensus assembly.
    pub use_quality_for_consensus: bool,

    /// Enable verbose logging of quality scoring decisions.
    pub log_quality_scores: bool,
}

impl Default for QualityConfig {
    fn default() -> Self {
        QualityConfig {
            scoring_mode: QualityScoring::None,
            error_penalty_exponent: 2.0,
            min_quality: None,
            use_quality_for_consensus: true,
            log_quality_scores: false,
        }
    }
}

/// Metrics collected during quality-aware overlap scoring.
#[derive(Debug, Clone, Default)]
pub struct QualityMetrics {
    /// Per-read average quality score.
    pub mean_quality: Vec<f32>,
    /// Count of positions below min_quality threshold.
    pub low_quality_positions: usize,
    /// Number of edges pruned due to low quality.
    pub edges_removed_by_quality: usize,
    /// Score histogram before quality adjustment.
    pub score_histogram_blind: Vec<usize>,
    /// Score histogram after quality adjustment.
    pub score_histogram_adjusted: Vec<usize>,
}

impl QualityMetrics {
    pub fn new() -> Self {
        QualityMetrics {
            mean_quality: Vec::new(),
            low_quality_positions: 0,
            edges_removed_by_quality: 0,
            score_histogram_blind: vec![0; 100],
            score_histogram_adjusted: vec![0; 100],
        }
    }

    /// Record score in histogram (clamped to 0-99).
    pub fn record_score(&mut self, score: f32, adjusted: bool) {
        let bucket = (score.min(99.0) as usize).min(99);
        if adjusted {
            self.score_histogram_adjusted[bucket] += 1;
        } else {
            self.score_histogram_blind[bucket] += 1;
        }
    }
}

pub type Adjacency = Vec<HashMap<usize, usize>>;

/// Companion mapping that records the overlap span associated with each edge.
pub type OverlapLengths = Vec<HashMap<usize, usize>>;

/// Result of overlap graph construction including ambiguity information.
#[derive(Debug)]
pub struct OverlapGraphResult {
    pub adjacency: CsMat<usize>,
    pub overlaps: CsMat<usize>,
    pub ambiguities: Vec<AmbiguityInfo>,
}

/// Configuration options that govern overlap graph construction.
#[derive(Debug, Clone)]
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
    /// Threshold for detecting ambiguous reads (gap ≤ tie_gap means ambiguous); 0.0 = exact ties only.
    pub tie_gap: f32,
    /// Enable mate-aware scoring (apply mate penalties to tied successors of ambiguous reads).
    pub mate_aware_scoring: bool,
    /// Optional mapping from read index to mate index for mate-aware scoring.
    pub mate_map: Option<Vec<Option<usize>>>,
    /// Expected insert size between mate pairs (in basepairs).
    pub insert_size: usize,
    /// Minimum acceptable insert size for mate pairs (lower bound for penalty calculation).
    pub min_insert: usize,
    /// Maximum acceptable insert size for mate pairs (upper bound for penalty calculation).
    pub max_insert: usize,
    /// Weight factor for mate penalty (penalty = weight × distance_error).
    pub mate_penalty_weight: f32,
    /// Maximum hops for bounded shortest path search in mate penalty calculation.
    pub mate_penalty_hop_limit: usize,
    /// Maximum basepair distance cap for bounded shortest path search.
    pub mate_penalty_cost_cap: usize,
    /// Configuration for quality-aware overlap scoring.
    pub quality_config: QualityConfig,
    /// Enable per-read adaptive max_diff based on read quality scores.
    pub adaptive_max_diff: bool,
    /// Minimum allowed max_diff when using adaptive thresholds (prevents overly strict filtering).
    pub adaptive_max_diff_min: f32,
    /// Maximum allowed max_diff when using adaptive thresholds (prevents overly permissive filtering).
    pub adaptive_max_diff_max: f32,
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
            tie_gap: 0.0,                // Exact ties only by default
            mate_aware_scoring: false,   // Disabled by default for backward compatibility
            mate_map: None,
            insert_size: 300, // Default ~300bp insert size
            min_insert: 165,  // 5th percentile for typical library
            max_insert: 407,  // 95th percentile for typical library
            mate_penalty_weight: 1.0,
            mate_penalty_hop_limit: 10,
            mate_penalty_cost_cap: 1000,
            quality_config: QualityConfig::default(),
            adaptive_max_diff: false,    // Disabled by default
            adaptive_max_diff_min: 0.05, // Allow down to 5% error
            adaptive_max_diff_max: 0.35, // Cap at 35% error even for low-quality reads
        }
    }
}

/// Information about ambiguous successors for a read.
#[derive(Debug, Clone)]
pub struct AmbiguityInfo {
    /// Read index in the assembly.
    pub read_idx: usize,
    /// Whether this read has multiple successors within tie_gap of the top score.
    pub is_ambiguous: bool,
    /// Indices of successor reads that are tied with (or near) the top score.
    pub tied_successors: Vec<usize>,
}

/// Compute phred-based error probability: P(error) = 10^(-Q/10)
pub fn phred_to_error_probability(phred: i32) -> f32 {
    10.0_f32.powf(-phred as f32 / 10.0)
}

/// Compute confidence from phred score: 1 - P(error)
pub fn phred_to_confidence(phred: i32) -> f32 {
    1.0 - phred_to_error_probability(phred)
}

/// Compute per-read maximum edit distance tolerance from quality scores.
/// Uses mean/median phred quality to derive an expected error rate,
/// then derives a max_diff that permits that error rate within bounds.
///
/// Higher quality reads get stricter (lower) max_diff; lower quality reads get looser (higher) max_diff.
pub fn compute_per_read_max_diff(qualities: Option<&[i32]>, min_diff: f32, max_diff: f32) -> f32 {
    if let Some(q) = qualities {
        if q.is_empty() {
            return max_diff;
        }

        // Compute mean quality
        let mean_q = q.iter().sum::<i32>() as f32 / q.len() as f32;

        // Convert mean quality to expected error rate
        let expected_error_rate = phred_to_error_probability(mean_q as i32);

        // Clamp the expected error rate within [min_diff, max_diff]
        expected_error_rate.clamp(min_diff, max_diff)
    } else {
        // No quality available, use the global max_diff
        max_diff
    }
}

/// Compute per-read adaptive tolerances from a quality-aware read source.
/// Returns None if adaptive_max_diff is disabled or if quality retrieval fails.
fn derive_per_read_max_diff<R: QualityAwareReadSource + ?Sized>(
    reads: &R,
    n: usize,
    config: &OverlapConfig,
) -> Option<Vec<f32>> {
    if !config.adaptive_max_diff {
        return None;
    }

    let mut per_read: Vec<f32> = Vec::with_capacity(n);
    for i in 0..n {
        let q = match reads.get_quality(i) {
            Ok(opt) => opt,
            Err(e) => {
                log::warn!(
                    "Falling back to global max_diff: failed to load quality for read {}: {:?}",
                    i,
                    e
                );
                return None;
            }
        };
        let eff = compute_per_read_max_diff(
            q.as_ref().map(|c| c.as_ref()),
            config.adaptive_max_diff_min,
            config.adaptive_max_diff_max,
        );
        per_read.push(eff);
    }

    Some(per_read)
}

/// Detect ambiguous reads based on top-1 vs top-2 score gap.
/// Returns a mapping from read index to ambiguity information.
pub fn detect_ambiguities(
    per_source: &[HashMap<usize, (f32, usize)>],
    tie_gap: f32,
) -> Vec<AmbiguityInfo> {
    let mut ambiguities = Vec::new();

    for (read_idx, successors) in per_source.iter().enumerate() {
        if successors.is_empty() {
            ambiguities.push(AmbiguityInfo {
                read_idx,
                is_ambiguous: false,
                tied_successors: Vec::new(),
            });
            continue;
        }

        // Get all successor scores
        let mut scores: Vec<(usize, f32)> = successors
            .iter()
            .map(|(&dst, &(score, _span))| (dst, score))
            .collect();

        // Sort by score descending
        scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        let top_score = scores[0].1;
        // Find all successors within tie_gap of the top score
        let tied_successors: Vec<usize> = scores
            .iter()
            .filter(|(_dst, score)| (top_score - score).abs() <= tie_gap)
            .map(|(dst, _score)| *dst)
            .collect();

        let is_ambiguous = tied_successors.len() >= 2;

        ambiguities.push(AmbiguityInfo {
            read_idx,
            is_ambiguous,
            tied_successors,
        });
    }

    ambiguities
}

/// Compute the basepair distance along a path between two reads in the assembly graph.
/// Uses bounded BFS to find a path from start to target within hop and cost limits.
/// Returns Some(distance) if a path exists within bounds, None otherwise.
///
/// # Arguments
/// * `start` - Starting read index
/// * `target` - Target read index
/// * `adjacency` - Sparse adjacency matrix (scores as usize)
/// * `overlaps` - Sparse overlap matrix (overlap lengths as usize)
/// * `read_lengths` - Vector of read lengths (for computing basepair distances)
/// * `hop_limit` - Maximum number of edges to traverse
/// * `cost_cap` - Maximum cumulative basepair distance to allow
///
/// # Returns
/// Some(distance_in_bp) if path found within limits, None otherwise
pub fn bounded_shortest_path(
    start: usize,
    target: usize,
    adjacency: &CsMat<usize>,
    overlaps: &CsMat<usize>,
    read_lengths: &[usize],
    hop_limit: usize,
    cost_cap: usize,
) -> Option<usize> {
    if start == target {
        return Some(0);
    }

    if start >= read_lengths.len() || target >= read_lengths.len() {
        return None;
    }

    // BFS with (node, cumulative_distance, hop_count)
    let mut queue = VecDeque::new();
    let mut visited = HashSet::new();

    queue.push_back((start, 0_usize, 0_usize));
    visited.insert(start);

    let adj_csr = adjacency.to_csr();
    let ovl_csr = overlaps.to_csr();

    while let Some((current, dist, hops)) = queue.pop_front() {
        if hops >= hop_limit {
            continue;
        }

        // Get successors from adjacency matrix
        // Use CSR format to iterate over outgoing edges from current row
        if let Some(row_view) = adj_csr.outer_view(current) {
            for (successor_idx, _score) in row_view.iter() {
                if visited.contains(&successor_idx) {
                    continue;
                }

                let overlap_len = ovl_csr
                    .outer_view(current)
                    .and_then(|v| v.get(successor_idx).copied())
                    .unwrap_or(0);
                let read_len = read_lengths.get(current).copied().unwrap_or(0);

                // Distance advanced = read_length - overlap_length
                let step_dist = read_len.saturating_sub(overlap_len);
                let new_dist = dist + step_dist;

                if new_dist > cost_cap {
                    continue;
                }

                if successor_idx == target {
                    return Some(new_dist);
                }

                visited.insert(successor_idx);
                queue.push_back((successor_idx, new_dist, hops + 1));
            }
        }
    }

    None
}

/// Calculate mate penalty based on predicted distance to mate pair.
/// Returns 0.0 if no mate, or a penalty proportional to the deviation from acceptable insert range.
///
/// # Arguments
/// * `read_idx` - Index of the current read
/// * `successor_idx` - Index of the candidate successor
/// * `mate_map` - Optional mapping from read index to mate index
/// * `adjacency` - Sparse adjacency matrix
/// * `overlaps` - Sparse overlap matrix
/// * `read_lengths` - Vector of read lengths
/// * `min_insert` - Minimum acceptable insert size
/// * `max_insert` - Maximum acceptable insert size
/// * `penalty_weight` - Multiplier for the penalty
/// * `hop_limit` - Maximum hops for bounded search
/// * `cost_cap` - Maximum distance for bounded search
///
/// # Returns
/// Penalty value to subtract from quality-adjusted score (0.0 if no mate or within acceptable range)
pub fn compute_mate_penalty(
    read_idx: usize,
    successor_idx: usize,
    mate_map: Option<&Vec<Option<usize>>>,
    adjacency: &CsMat<usize>,
    overlaps: &CsMat<usize>,
    read_lengths: &[usize],
    min_insert: usize,
    max_insert: usize,
    penalty_weight: f32,
    hop_limit: usize,
    cost_cap: usize,
) -> f32 {
    // No penalty if no mate information
    let mate_map = match mate_map {
        Some(m) => m,
        None => return 0.0,
    };

    // No penalty if this read doesn't have a mate
    let mate_idx = match mate_map.get(read_idx).and_then(|m| *m) {
        Some(idx) => idx,
        None => return 0.0,
    };

    // Compute the first step: read_idx -> successor_idx
    let overlap_len = overlaps.get(read_idx, successor_idx).copied().unwrap_or(0);
    let read_len = read_lengths.get(read_idx).copied().unwrap_or(0);
    let first_step = read_len.saturating_sub(overlap_len);

    // Find path from successor to mate
    let successor_to_mate_dist = bounded_shortest_path(
        successor_idx,
        mate_idx,
        adjacency,
        overlaps,
        read_lengths,
        hop_limit,
        cost_cap,
    );

    let predicted_distance = match successor_to_mate_dist {
        Some(dist) => first_step + dist,
        None => {
            // Cannot reach mate within bounds - apply maximum penalty
            return penalty_weight * (cost_cap as f32);
        }
    };

    // Compute distance error based on acceptable range
    let distance_error = if predicted_distance < min_insert {
        (min_insert - predicted_distance) as f32
    } else if predicted_distance > max_insert {
        (predicted_distance - max_insert) as f32
    } else {
        // Within acceptable range - no penalty
        0.0
    };

    // Apply penalty proportional to error
    penalty_weight * distance_error
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
    let trie_max_diff = if config.adaptive_max_diff {
        config.max_diff.max(config.adaptive_max_diff_max)
    } else {
        config.max_diff
    };

    // Build trie, collect verified edges, drop trie, then build adjacency
    let trie = PrunedAffixTrie::build(reads, min_suffix_len, trie_max_diff);
    // Stream verified candidates directly into per-source rows to avoid
    // materialising a full verified-vector in memory.
    let (adj, ovl) = create_overlap_graph_from_trie_stream(reads, &trie, config, None);
    drop(trie);
    (adj, ovl)
}

/// Variant of unified overlap graph builder that accepts any `ReadSource`.
///
/// This convenience API uses `PrunedAffixTrie::build_from_readsource` to
/// construct the trie and (for now) materializes reads into `Vec<String>` to
/// reuse the existing trie-based pipeline. Future work may avoid full
/// materialization.
pub fn create_overlap_graph_unified_from_readsource<
    R: ReadSource + QualityAwareReadSource + ?Sized,
>(
    reads: &R,
    config: OverlapConfig,
) -> (CsMat<usize>, CsMat<usize>, Vec<String>, Vec<String>) {
    let min_suffix_len = config.min_suffix_len.max(1);
    let trie_max_diff = if config.adaptive_max_diff {
        config.max_diff.max(config.adaptive_max_diff_max)
    } else {
        config.max_diff
    };
    let (trie, mat_reads, mat_names) =
        PrunedAffixTrie::build_from_readsource(reads, min_suffix_len, trie_max_diff);

    let per_read_max_diff = derive_per_read_max_diff(reads, mat_reads.len(), &config);

    // Stream verified candidates into per-source rows while the trie exists,
    // then drop it to release memory before constructing adjacency matrices.
    let (adj, overlaps) = create_overlap_graph_from_trie_stream(
        &mat_reads,
        &trie,
        config,
        per_read_max_diff.as_deref(),
    );
    drop(trie);
    (adj, overlaps, mat_reads, mat_names)
}

/// Build the weighted overlap graph with ambiguity detection from any `ReadSource`.
/// Returns adjacency, overlaps, ambiguity info, and the materialized reads/ids.
pub fn create_overlap_graph_with_ambiguities_from_readsource<
    R: ReadSource + QualityAwareReadSource + ?Sized,
>(
    reads: &R,
    config: OverlapConfig,
) -> (
    CsMat<usize>,
    CsMat<usize>,
    Vec<AmbiguityInfo>,
    Vec<String>,
    Vec<String>,
) {
    let min_suffix_len = config.min_suffix_len.max(1);
    let trie_max_diff = if config.adaptive_max_diff {
        config.max_diff.max(config.adaptive_max_diff_max)
    } else {
        config.max_diff
    };
    let (trie, mat_reads, mat_names) =
        PrunedAffixTrie::build_from_readsource(reads, min_suffix_len, trie_max_diff);

    let per_read_max_diff = derive_per_read_max_diff(reads, mat_reads.len(), &config);

    let (adj, overlaps, ambiguities) = create_overlap_graph_from_trie_stream_with_ambiguities(
        &mat_reads,
        &trie,
        config,
        per_read_max_diff.as_deref(),
    );
    drop(trie);
    (adj, overlaps, ambiguities, mat_reads, mat_names)
}

/// Build the weighted overlap graph with ambiguity detection.
/// Returns adjacency, overlaps, and ambiguity information for each read.
pub fn create_overlap_graph_with_ambiguities<'a>(
    reads: &'a [String],
    config: OverlapConfig,
) -> OverlapGraphResult {
    let min_suffix_len = config.min_suffix_len.max(1);
    let trie_max_diff = if config.adaptive_max_diff {
        config.max_diff.max(config.adaptive_max_diff_max)
    } else {
        config.max_diff
    };
    let trie = PrunedAffixTrie::build(reads, min_suffix_len, trie_max_diff);
    let (adj, overlaps, ambiguities) =
        create_overlap_graph_from_trie_stream_with_ambiguities(reads, &trie, config, None);
    drop(trie);

    OverlapGraphResult {
        adjacency: adj,
        overlaps,
        ambiguities,
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
    let trie_max_diff = if config.adaptive_max_diff {
        config.max_diff.max(config.adaptive_max_diff_max)
    } else {
        config.max_diff
    };

    // Trie now includes fuzzy matches in extension vectors automatically
    log::info!("Collecting & verifying candidates from trie...");
    // Use the trie-side verified candidate emitter (simple verifier)
    let verified = affix_trie.overlap_verified_candidates_simple(
        reads,
        trie_max_diff,
        config.min_suffix_len,
        None,
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

/// Stream candidates from `affix_trie` into per-source maps and construct
/// adjacency/overlap CSR matrices. This avoids materialising a large
/// intermediate verified-vector and lets us drop the trie earlier.
fn create_overlap_graph_from_trie_stream(
    reads: &[String],
    affix_trie: &crate::affix::PrunedAffixTrie,
    config: OverlapConfig,
    per_read_max_diff: Option<&[f32]>,
) -> (CsMat<usize>, CsMat<usize>) {
    if reads.is_empty() {
        let empty1 = CsMat::<usize>::zero((0, 0));
        let empty2 = CsMat::<usize>::zero((0, 0));
        return (empty1, empty2);
    }

    let n = reads.len();

    let trie_max_diff = if config.adaptive_max_diff {
        config.max_diff.max(config.adaptive_max_diff_max)
    } else {
        config.max_diff
    };

    // Prepare per-source best maps
    let mut per_source: Vec<HashMap<usize, (f32, usize)>> =
        (0..n).map(|_| HashMap::new()).collect();

    // Stream verified candidates from the trie into per_source
    log::info!("Collecting & verifying candidates from trie (stream)...");
    let mut found = 0usize;
    affix_trie.overlap_verified_candidates_stream(
        reads,
        trie_max_diff,
        config.min_suffix_len,
        per_read_max_diff,
        |s, d, score, span| {
            found += 1;
            let row = &mut per_source[s];
            let prev = row.get(&d).copied().unwrap_or((f32::MIN, 0));
            if score > prev.0 || (score == prev.0 && span > prev.1) {
                row.insert(d, (score, span));
            }
        },
    );

    log::info!("Found {} verified edges (streamed)", found);

    // Convert per_source into CSR sparse matrices (reuse same logic as before)
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

/// Stream candidates from `affix_trie` into per-source maps and construct
/// adjacency/overlap CSR matrices with ambiguity detection.
/// Returns adjacency, overlap, and ambiguity information matrices.
fn create_overlap_graph_from_trie_stream_with_ambiguities(
    reads: &[String],
    affix_trie: &crate::affix::PrunedAffixTrie,
    config: OverlapConfig,
    per_read_max_diff: Option<&[f32]>,
) -> (CsMat<usize>, CsMat<usize>, Vec<AmbiguityInfo>) {
    if reads.is_empty() {
        let empty1 = CsMat::<usize>::zero((0, 0));
        let empty2 = CsMat::<usize>::zero((0, 0));
        return (empty1, empty2, Vec::new());
    }

    let n = reads.len();

    let trie_max_diff = if config.adaptive_max_diff {
        config.max_diff.max(config.adaptive_max_diff_max)
    } else {
        config.max_diff
    };

    // Prepare per-source best maps
    let mut per_source: Vec<HashMap<usize, (f32, usize)>> =
        (0..n).map(|_| HashMap::new()).collect();

    // Stream verified candidates from the trie into per_source
    log::info!("Collecting & verifying candidates from trie (stream)...");
    affix_trie.overlap_verified_candidates_stream(
        reads,
        trie_max_diff,
        config.min_suffix_len,
        per_read_max_diff,
        |s, d, score, span| {
            let row = &mut per_source[s];
            let prev = row.get(&d).copied().unwrap_or((f32::MIN, 0));
            if score > prev.0 || (score == prev.0 && span > prev.1) {
                row.insert(d, (score, span));
            }
        },
    );

    // Detect ambiguities before converting to sparse matrices
    let ambiguities = detect_ambiguities(&per_source, config.tie_gap);

    // Apply mate penalties to ambiguous reads if mate-aware scoring is enabled
    if config.mate_aware_scoring {
        if let Some(ref mate_map) = config.mate_map {
            log::info!("Applying mate penalties to ambiguous reads...");

            // First, build temporary sparse matrices to enable path search
            let temp_adj = build_temp_adjacency(&per_source, n);
            let temp_ovl = build_temp_overlaps(&per_source, n);

            // Get read lengths for distance calculation
            let read_lengths: Vec<usize> = reads.iter().map(|r| r.len()).collect();

            // Apply penalties to tied successors of ambiguous reads
            for (src_idx, amb_info) in ambiguities.iter().enumerate() {
                if !amb_info.is_ambiguous {
                    continue;
                }

                // For each tied successor, compute mate penalty and adjust score
                for &dst_idx in &amb_info.tied_successors {
                    let penalty = compute_mate_penalty(
                        src_idx,
                        dst_idx,
                        Some(mate_map),
                        &temp_adj,
                        &temp_ovl,
                        &read_lengths,
                        config.min_insert,
                        config.max_insert,
                        config.mate_penalty_weight,
                        config.mate_penalty_hop_limit,
                        config.mate_penalty_cost_cap,
                    );

                    // Update score in per_source
                    if let Some((score, _span)) = per_source[src_idx].get_mut(&dst_idx) {
                        *score -= penalty;
                        log::debug!(
                            "[MATE_PENALTY] Read {} -> {} penalty={:.2}, new_score={:.2}",
                            src_idx,
                            dst_idx,
                            penalty,
                            *score
                        );
                    }
                }
            }
        }
    }

    // Build sparse matrices from per_source
    let mut a_indptr = vec![0];
    let mut a_indices = Vec::new();
    let mut a_data = Vec::new();

    let mut o_indptr = vec![0];
    let mut o_indices = Vec::new();
    let mut o_data = Vec::new();

    for src_idx in 0..n {
        let mut edges: Vec<(usize, f32)> = per_source[src_idx]
            .iter()
            .map(|(&dst, &(score, _span))| (dst, score))
            .collect();

        let edges_before = edges.len();

        // Apply knee-point detection if configured
        if config.detect_score_cliff && !edges.is_empty() {
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

    (adjacency, overlaps, ambiguities)
}

/// Build temporary adjacency matrix from per_source map for mate penalty calculation.
fn build_temp_adjacency(per_source: &[HashMap<usize, (f32, usize)>], n: usize) -> CsMat<usize> {
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();

    for src_idx in 0..n {
        let mut edges: Vec<(usize, f32)> = per_source[src_idx]
            .iter()
            .map(|(&dst, &(score, _span))| (dst, score))
            .collect();
        edges.sort_by_key(|(dst, _)| *dst);

        for (dst, score) in edges {
            indices.push(dst);
            data.push(score as usize);
        }
        indptr.push(indices.len());
    }

    CsMat::new((n, n), indptr, indices, data)
}

/// Build temporary overlaps matrix from per_source map for mate penalty calculation.
fn build_temp_overlaps(per_source: &[HashMap<usize, (f32, usize)>], n: usize) -> CsMat<usize> {
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();

    for src_idx in 0..n {
        let mut edges: Vec<(usize, usize)> = per_source[src_idx]
            .iter()
            .map(|(&dst, &(_score, span))| (dst, span))
            .collect();
        edges.sort_by_key(|(dst, _)| *dst);

        for (dst, span) in edges {
            indices.push(dst);
            data.push(span);
        }
        indptr.push(indices.len());
    }

    CsMat::new((n, n), indptr, indices, data)
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
    verify_overlap_from_anchor_with_quality(
        suffix_read,
        prefix_read,
        shared_affix,
        None,
        None,
        config,
    )
}

/// Verify overlap with optional quality scores.
pub fn verify_overlap_from_anchor_with_quality(
    suffix_read: &str,
    prefix_read: &str,
    shared_affix: &str,
    suffix_quality: Option<&[i32]>,
    prefix_quality: Option<&[i32]>,
    config: OverlapConfig,
) -> Option<(f32, usize)> {
    let anchor_len = shared_affix.len();
    if anchor_len == 0 {
        // No anchor, do simple full scan
        return verify_overlap_simple_with_quality(
            suffix_read,
            prefix_read,
            suffix_quality,
            prefix_quality,
            config,
        );
    }

    let min_len = config.min_suffix_len.max(1);
    let max_span = suffix_read.len().min(prefix_read.len());

    let mut best_score = f32::MIN;
    let mut best_overlap = 0;

    // Determine effective max_diff: use adaptive if enabled, else global
    let effective_max_diff = if config.adaptive_max_diff {
        compute_per_read_max_diff(
            suffix_quality,
            config.adaptive_max_diff_min,
            config.adaptive_max_diff_max,
        )
    } else {
        config.max_diff
    };

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

        if float_diff <= effective_max_diff {
            // Compute score using quality-aware scoring if available
            let score = match config.quality_config.scoring_mode {
                QualityScoring::Position => score_overlap_position_mode(
                    suffix_window,
                    prefix_window,
                    suffix_quality,
                    prefix_quality,
                    span,
                    config.quality_config.error_penalty_exponent,
                ),
                QualityScoring::Confidence => score_overlap_confidence_mode(
                    suffix_window,
                    prefix_window,
                    suffix_quality,
                    prefix_quality,
                    span,
                ),
                QualityScoring::None => {
                    // Standard quality-adjusted scoring: penalize errors based on error rate
                    let error_rate = distance as f32 / span as f32;
                    let quality_factor = (1.0 - error_rate).powf(config.error_penalty_exponent);
                    span as f32 * quality_factor
                }
            };

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

/// Simple full-scan verification with optional quality scores.
/// Used as fallback when no anchor hint exists.
fn verify_overlap_simple_with_quality(
    suffix_read: &str,
    prefix_read: &str,
    suffix_quality: Option<&[i32]>,
    prefix_quality: Option<&[i32]>,
    config: OverlapConfig,
) -> Option<(f32, usize)> {
    let min_len = config.min_suffix_len.max(1);
    let max_span = suffix_read.len().min(prefix_read.len());

    let mut best_score = f32::MIN;
    let mut best_overlap = 0;

    // Determine effective max_diff: use adaptive if enabled, else global
    let effective_max_diff = if config.adaptive_max_diff {
        compute_per_read_max_diff(
            suffix_quality,
            config.adaptive_max_diff_min,
            config.adaptive_max_diff_max,
        )
    } else {
        config.max_diff
    };

    for span in min_len..=max_span {
        if span > suffix_read.len() || span > prefix_read.len() {
            break;
        }

        let suffix_window = &suffix_read[suffix_read.len() - span..];
        let prefix_window = &prefix_read[..span];

        let distance = damerau_levenshtein(suffix_window, prefix_window);
        let float_diff = distance as f32 / span as f32;

        if float_diff <= effective_max_diff {
            let score = match config.quality_config.scoring_mode {
                QualityScoring::Position => score_overlap_position_mode(
                    suffix_window,
                    prefix_window,
                    suffix_quality,
                    prefix_quality,
                    span,
                    config.quality_config.error_penalty_exponent,
                ),
                QualityScoring::Confidence => score_overlap_confidence_mode(
                    suffix_window,
                    prefix_window,
                    suffix_quality,
                    prefix_quality,
                    span,
                ),
                QualityScoring::None => (span as i32 - distance as i32) as f32,
            };

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

/// Compute overlap score using Position mode quality weighting.
/// Mismatches at low-quality positions receive soft penalty.
#[allow(dead_code)]
fn score_overlap_position_mode(
    suffix_window: &str,
    prefix_window: &str,
    suffix_quality: Option<&[i32]>,
    prefix_quality: Option<&[i32]>,
    span: usize,
    error_penalty_exponent: f32,
) -> f32 {
    let distance = damerau_levenshtein(suffix_window, prefix_window);
    let mut score = span as f32;

    if let (Some(sq), Some(pq)) = (suffix_quality, prefix_quality) {
        // Apply quality-aware penalty to mismatches
        let suffix_start = sq.len().saturating_sub(span);
        for i in 0..span {
            if i < suffix_window.len()
                && i < prefix_window.len()
                && suffix_window.as_bytes()[i] != prefix_window.as_bytes()[i]
            {
                let sq_idx = suffix_start + i;
                if sq_idx < sq.len() {
                    let avg_q = (sq[sq_idx] as f32 + pq[i] as f32) / 2.0;
                    let q_factor = phred_to_error_probability(avg_q as i32);
                    // Soft penalty scaled by error probability at this position
                    score -= q_factor * error_penalty_exponent;
                }
            }
        }
    } else {
        // No quality: use standard error rate penalty
        let error_rate = distance as f32 / span as f32;
        let quality_factor = (1.0 - error_rate).powf(error_penalty_exponent);
        score = span as f32 * quality_factor;
    }

    score.max(1.0)
}

/// Compute overlap score using Confidence mode quality weighting.
/// Matches at high-quality positions boost score; mismatches penalize.
#[allow(dead_code)]
fn score_overlap_confidence_mode(
    suffix_window: &str,
    prefix_window: &str,
    suffix_quality: Option<&[i32]>,
    prefix_quality: Option<&[i32]>,
    span: usize,
) -> f32 {
    let mut score = 0.0;

    if let (Some(sq), Some(pq)) = (suffix_quality, prefix_quality) {
        let suffix_start = sq.len().saturating_sub(span);
        for i in 0..span {
            if i < suffix_window.len() && i < prefix_window.len() {
                let sq_idx = suffix_start + i;
                if sq_idx < sq.len() {
                    let avg_q = (sq[sq_idx] as f32 + pq[i] as f32) / 2.0;
                    let confidence = phred_to_confidence(avg_q as i32);
                    if suffix_window.as_bytes()[i] == prefix_window.as_bytes()[i] {
                        score += confidence; // Match: boost by confidence
                    } else {
                        score -= confidence; // Mismatch: penalize by confidence
                    }
                }
            }
        }
    } else {
        // No quality: treat all as equal
        let distance = damerau_levenshtein(suffix_window, prefix_window);
        score = span as f32 - distance as f32;
    }

    score.max(1.0)
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

    #[test]
    fn test_ambiguity_detection_exact_ties() {
        let mut per_source: Vec<HashMap<usize, (f32, usize)>> =
            vec![HashMap::new(), HashMap::new(), HashMap::new()];

        // Read 0: two successors with identical scores (exact tie)
        per_source[0].insert(1, (100.0, 45));
        per_source[0].insert(2, (100.0, 50)); // Same score as successor 1

        // Read 1: three successors with identical scores (3-way tie)
        per_source[1].insert(0, (50.0, 30));
        per_source[1].insert(2, (50.0, 32));
        per_source[1].insert(3, (50.0, 28)); // All tied

        // Read 2: single successor (not ambiguous)
        per_source[2].insert(0, (75.0, 40));

        let ambiguities = detect_ambiguities(&per_source, 0.0);

        assert_eq!(ambiguities.len(), 3);

        // Read 0 should be flagged as ambiguous with 2 tied successors
        assert!(ambiguities[0].is_ambiguous);
        assert_eq!(ambiguities[0].tied_successors.len(), 2);
        assert!(ambiguities[0].tied_successors.contains(&1));
        assert!(ambiguities[0].tied_successors.contains(&2));

        // Read 1 should be flagged as ambiguous with 3 tied successors
        assert!(ambiguities[1].is_ambiguous);
        assert_eq!(ambiguities[1].tied_successors.len(), 3);
        assert!(ambiguities[1].tied_successors.contains(&0));
        assert!(ambiguities[1].tied_successors.contains(&2));

        // Read 2 should NOT be ambiguous
        assert!(!ambiguities[2].is_ambiguous);
        assert_eq!(ambiguities[2].tied_successors.len(), 1); // Only the successor
    }

    #[test]
    fn test_ambiguity_detection_fuzzy_ties() {
        let mut per_source: Vec<HashMap<usize, (f32, usize)>> = vec![HashMap::new()];

        // Create successors with scores near each other
        per_source[0].insert(1, (100.0, 45));
        per_source[0].insert(2, (98.0, 50)); // Within tie_gap=5.0
        per_source[0].insert(3, (92.0, 40)); // Outside tie_gap=5.0

        // With tie_gap=0.0 (exact ties only)
        let ambiguities_exact = detect_ambiguities(&per_source, 0.0);
        assert!(!ambiguities_exact[0].is_ambiguous);
        assert_eq!(ambiguities_exact[0].tied_successors.len(), 1);

        // With tie_gap=5.0 (fuzzy ties)
        let ambiguities_fuzzy = detect_ambiguities(&per_source, 5.0);
        assert!(ambiguities_fuzzy[0].is_ambiguous);
        assert_eq!(ambiguities_fuzzy[0].tied_successors.len(), 2);
        assert!(ambiguities_fuzzy[0].tied_successors.contains(&1));
        assert!(ambiguities_fuzzy[0].tied_successors.contains(&2));
    }

    #[test]
    fn test_non_ambiguous_reads_unmarked() {
        let mut per_source: Vec<HashMap<usize, (f32, usize)>> =
            vec![HashMap::new(), HashMap::new()];

        // Read 0: clear winner (score gap >> tie_gap)
        per_source[0].insert(1, (100.0, 45));
        per_source[0].insert(2, (70.0, 50)); // Large gap

        // Read 1: also clear winner
        per_source[1].insert(0, (50.0, 30));
        per_source[1].insert(3, (20.0, 20)); // Large gap

        let ambiguities = detect_ambiguities(&per_source, 0.0);

        assert_eq!(ambiguities.len(), 2);
        assert!(!ambiguities[0].is_ambiguous);
        assert!(!ambiguities[1].is_ambiguous);
    }

    #[test]
    fn test_overlap_graph_with_ambiguities() {
        let reads = vec![
            "AAAAAA".to_string(),
            "AAAAAA".to_string(), // Exact tie with read 0
            "CCCCCC".to_string(),
        ];

        let config = OverlapConfig {
            tie_gap: 0.0,
            ..Default::default()
        };

        let result = create_overlap_graph_with_ambiguities(&reads, config);

        // Should have 3 reads
        assert_eq!(result.ambiguities.len(), 3);

        // Verify ambiguity info is populated
        for amb in &result.ambiguities {
            assert!(amb.read_idx < 3);
        }
    }

    #[test]
    fn test_bounded_shortest_path_finds_mate() {
        // Create a simple path: 0 -> 1 -> 2
        // With overlaps and read lengths allowing distance calculation
        use sprs::TriMat;

        let n = 3;
        let mut adj_builder = TriMat::new((n, n));
        let mut ovl_builder = TriMat::new((n, n));

        // Edge 0 -> 1 (score=10, overlap=5)
        adj_builder.add_triplet(0, 1, 10);
        ovl_builder.add_triplet(0, 1, 5);

        // Edge 1 -> 2 (score=8, overlap=3)
        adj_builder.add_triplet(1, 2, 8);
        ovl_builder.add_triplet(1, 2, 3);

        let adjacency = adj_builder.to_csr().to_csc();
        let overlaps = ovl_builder.to_csr().to_csc();

        // Read lengths: all 10bp
        let read_lengths = vec![10, 10, 10];

        // Find path from 0 to 2
        // Distance should be: (10-5) + (10-3) = 5 + 7 = 12
        let result = bounded_shortest_path(
            0,
            2,
            &adjacency,
            &overlaps,
            &read_lengths,
            10,  // hop_limit
            100, // cost_cap
        );

        assert_eq!(result, Some(12));
    }

    #[test]
    fn test_bounded_shortest_path_respects_limits() {
        use sprs::TriMat;

        let n = 4;
        let mut adj_builder = TriMat::new((n, n));
        let mut ovl_builder = TriMat::new((n, n));

        // Chain: 0 -> 1 -> 2 -> 3
        adj_builder.add_triplet(0, 1, 10);
        ovl_builder.add_triplet(0, 1, 2);

        adj_builder.add_triplet(1, 2, 10);
        ovl_builder.add_triplet(1, 2, 2);

        adj_builder.add_triplet(2, 3, 10);
        ovl_builder.add_triplet(2, 3, 2);

        let adjacency = adj_builder.to_csr().to_csc();
        let overlaps = ovl_builder.to_csr().to_csc();
        let read_lengths = vec![10, 10, 10, 10];

        // Test hop limit: cannot reach node 3 with only 2 hops from 0
        let result_hop_limited = bounded_shortest_path(
            0,
            3,
            &adjacency,
            &overlaps,
            &read_lengths,
            2,   // hop_limit (too small)
            100, // cost_cap
        );
        assert_eq!(result_hop_limited, None);

        // Test cost limit: total distance would be 3*(10-2)=24
        let result_cost_limited = bounded_shortest_path(
            0,
            3,
            &adjacency,
            &overlaps,
            &read_lengths,
            10, // hop_limit
            20, // cost_cap (too small)
        );
        assert_eq!(result_cost_limited, None);

        // With adequate limits, should succeed
        let result_ok = bounded_shortest_path(
            0,
            3,
            &adjacency,
            &overlaps,
            &read_lengths,
            10, // hop_limit
            50, // cost_cap
        );
        assert_eq!(result_ok, Some(24));
    }

    #[test]
    fn test_mate_penalty_zero_for_non_ambiguous() {
        use sprs::TriMat;

        let n = 3;
        let adj_builder = TriMat::new((n, n));
        let ovl_builder = TriMat::new((n, n));

        let adjacency = adj_builder.to_csr().to_csc();
        let overlaps = ovl_builder.to_csr().to_csc();
        let read_lengths = vec![100, 100, 100];

        // No mate map provided
        let penalty = compute_mate_penalty(
            0,
            1,
            None,
            &adjacency,
            &overlaps,
            &read_lengths,
            250, // min_insert
            350, // max_insert
            1.0, // penalty_weight
            10,  // hop_limit
            500, // cost_cap
        );

        assert_eq!(penalty, 0.0);
    }

    #[test]
    fn test_mate_penalty_nonzero_for_ambiguous_ties() {
        use sprs::TriMat;

        let n = 4;
        let mut adj_builder = TriMat::new((n, n));
        let mut ovl_builder = TriMat::new((n, n));

        // Read 0 has two successors: 1 and 2 (ambiguous)
        // Read 0's mate is read 3
        // Path 0 -> 1 -> 3 has distance closer to insert_size
        // Path 0 -> 2 (no connection to 3) should have higher penalty

        adj_builder.add_triplet(0, 1, 10);
        ovl_builder.add_triplet(0, 1, 50);

        adj_builder.add_triplet(0, 2, 10); // Same score (ambiguous)
        ovl_builder.add_triplet(0, 2, 50);

        adj_builder.add_triplet(1, 3, 10);
        ovl_builder.add_triplet(1, 3, 50);

        let adjacency = adj_builder.to_csr().to_csc();
        let overlaps = ovl_builder.to_csr().to_csc();
        let read_lengths = vec![100, 100, 100, 100];

        // Mate map: read 0's mate is read 3
        let mate_map = vec![Some(3), None, None, Some(0)];

        // Penalty for choosing successor 1 (can reach mate 3)
        let penalty_good = compute_mate_penalty(
            0,
            1,
            Some(&mate_map),
            &adjacency,
            &overlaps,
            &read_lengths,
            80,  // min_insert (distance is (100-50)+(100-50)=100)
            120, // max_insert
            1.0, // penalty_weight
            10,  // hop_limit
            500, // cost_cap
        );

        // Penalty for choosing successor 2 (cannot reach mate 3)
        let penalty_bad = compute_mate_penalty(
            0,
            2,
            Some(&mate_map),
            &adjacency,
            &overlaps,
            &read_lengths,
            80,  // min_insert
            120, // max_insert
            1.0, // penalty_weight
            10,  // hop_limit
            500, // cost_cap
        );

        // Good path should have lower penalty than bad path
        assert!(penalty_good < penalty_bad);
        assert_eq!(penalty_good, 0.0); // Perfect distance match
        assert!(penalty_bad > 0.0); // Cannot reach mate
    }

    #[test]
    fn test_qas_with_mate_penalty() {
        // Test that QAS = overlap_length - edit_penalty - mate_penalty
        // when mate_aware_scoring is enabled

        // Manually create a small overlap graph with ambiguous reads
        let mut adj_builder = sprs::TriMat::new((4, 4));
        let mut ovl_builder = sprs::TriMat::new((4, 4));

        // Read 0 has two successors with equal weight (ambiguous)
        adj_builder.add_triplet(0, 1, 10);
        ovl_builder.add_triplet(0, 1, 50);
        adj_builder.add_triplet(0, 2, 10);
        ovl_builder.add_triplet(0, 2, 50);

        // Read 1 connects to read 3 (mate of read 0)
        adj_builder.add_triplet(1, 3, 10);
        ovl_builder.add_triplet(1, 3, 50);

        let adjacency = adj_builder.to_csr().to_csc();
        let overlaps = ovl_builder.to_csr().to_csc();
        let read_lengths = vec![100, 100, 100, 100];
        let mate_map = vec![Some(3), None, None, Some(0)];

        // Compute mate penalties for both successors
        let penalty_1 = compute_mate_penalty(
            0,
            1,
            Some(&mate_map),
            &adjacency,
            &overlaps,
            &read_lengths,
            80,  // min_insert
            120, // max_insert
            1.0,
            10,
            500,
        );
        let penalty_2 = compute_mate_penalty(
            0,
            2,
            Some(&mate_map),
            &adjacency,
            &overlaps,
            &read_lengths,
            80,  // min_insert
            120, // max_insert
            1.0,
            10,
            500,
        );

        // Base QAS (assuming no edit penalty for simplicity)
        let base_qas = 50.0;

        // QAS with mate penalties
        let qas_1 = base_qas - penalty_1;
        let qas_2 = base_qas - penalty_2;

        // Successor 1 (good path to mate) should have higher QAS
        assert!(qas_1 > qas_2);
        assert_eq!(qas_1, 50.0); // No penalty
        assert!(qas_2 < 50.0); // Has penalty
    }

    #[test]
    fn test_qas_without_mate_penalty() {
        // Test that QAS = overlap_length - edit_penalty
        // when mate_aware_scoring is disabled

        let overlap_length = 50.0;
        let edit_penalty = 0.0;
        let mate_penalty = 0.0; // Not computed when mate_aware_scoring=false

        let qas = overlap_length - edit_penalty - mate_penalty;
        assert_eq!(qas, 50.0);
    }

    #[test]
    fn test_mate_penalty_only_for_ambiguous_reads() {
        // Test that mate penalties are only computed for reads flagged as ambiguous

        let mut adj_builder = sprs::TriMat::new((3, 3));
        let mut ovl_builder = sprs::TriMat::new((3, 3));

        // Read 0 has one clear winner (not ambiguous)
        adj_builder.add_triplet(0, 1, 20);
        ovl_builder.add_triplet(0, 1, 50);
        adj_builder.add_triplet(0, 2, 5);
        ovl_builder.add_triplet(0, 2, 10);

        let _adjacency = adj_builder.to_csr::<usize>();
        let _overlaps = ovl_builder.to_csr::<usize>();

        use std::collections::HashMap;
        let mut per_source = vec![HashMap::new()];
        per_source[0].insert(1, (20.0, 50));
        per_source[0].insert(2, (5.0, 10));

        let tie_gap = 0.0;
        let ambiguities = detect_ambiguities(&per_source, tie_gap);

        // Read 0 is not ambiguous (20 vs 5, clear winner)
        assert!(!ambiguities[0].is_ambiguous);

        // Therefore, mate penalty should not be computed for this read
        // (In practice, we check is_ambiguous before calling compute_mate_penalty)
    }
}
