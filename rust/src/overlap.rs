//! Overlap graph construction primitives translated from the Python reference.

use std::collections::HashMap;

use log::warn;
use sprs::CsMat;
use strsim::damerau_levenshtein;

use crate::suffix::{AffixArray, AffixKind};

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
}

impl Default for OverlapConfig {
    fn default() -> Self {
        Self {
            max_diff: 0.25,
            min_suffix_len: crate::suffix::DEFAULT_MIN_SUFFIX_LEN,
            use_threads: false,
            max_workers: 1,
        }
    }
}

/// Compute the normalised Damerau–Levenshtein distance and overlap weight.
pub fn normalised_damerau_levenshtein_distance(
    suffix: &str,
    prefix: &str,
) -> Option<(f32, usize, usize)> {
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
    let score = overlap_len - distance;
    Some((diff, score, overlap_len))
}

/// Build the weighted overlap graph for the supplied reads.
/// Returns sparse matrices directly to avoid intermediate HashMap allocation.
pub fn create_overlap_graph(
    reads: &[String],
    mut affix_array: Option<AffixArray>,
    config: OverlapConfig,
) -> (AffixArray, CsMat<usize>, CsMat<usize>) {
    if config.use_threads {
        warn!(
            "Threaded overlap construction is not yet implemented; falling back to sequential mode"
        );
    }

    let min_suffix_len = config.min_suffix_len.max(1);

    let affix_array = affix_array
        .take()
        .unwrap_or_else(|| AffixArray::build(reads, min_suffix_len));

    if reads.is_empty() {
        let empty1 = CsMat::<usize>::zero((0, 0));
        let empty2 = CsMat::<usize>::zero((0, 0));
        return (affix_array, empty1, empty2);
    }

    // CSR builders
    let n = reads.len();
    let mut a_indptr: Vec<usize> = Vec::with_capacity(n + 1);
    let mut a_indices: Vec<usize> = Vec::new();
    let mut a_data: Vec<usize> = Vec::new();
    let mut o_indptr: Vec<usize> = Vec::with_capacity(n + 1);
    let mut o_indices: Vec<usize> = Vec::new();
    let mut o_data: Vec<usize> = Vec::new();

    a_indptr.push(0);
    o_indptr.push(0);

    for (read_idx, read) in reads.iter().enumerate() {
        if read.len() < min_suffix_len {
            continue;
        }

        let mut up_cut = false;
        let mut down_cut = false;
        let read_len = read.len();

        // Track best edge per target for this source row to avoid duplicate triplets
        let mut best_scores: HashMap<usize, usize> = HashMap::new();
        let mut best_spans: HashMap<usize, usize> = HashMap::new();

        for span in min_suffix_len..read_len {
            let suffix_seq = &read[read_len - span..];
            let start = read_len - span;
            let Some(anchor) = affix_array.suffix_anchor(read_idx, start) else {
                continue;
            };

            let mut forward = anchor + 1;
            while forward < affix_array.len() {
                let entry = affix_array.get(forward).expect("in-bounds forward entry");

                if entry.kind() == AffixKind::Suffix || entry.read_index() == read_idx {
                    forward += 1;
                    continue;
                }

                if entry.kind() != AffixKind::Prefix {
                    break;
                }

                if let Some((diff, score, overlap_len)) =
                    normalised_damerau_levenshtein_distance(suffix_seq, entry.affix(reads))
                {
                    if diff > config.max_diff {
                        down_cut = true;
                        break;
                    }

                    let dst = entry.read_index();
                    let prev = best_scores.get(&dst).copied().unwrap_or(0);
                    if score > prev {
                        best_scores.insert(dst, score);
                        best_spans.insert(dst, overlap_len);
                    }
                }

                forward += 1;
            }

            let mut backward = anchor as isize - 1;
            while backward >= 0 {
                let idx = backward as usize;
                let entry = affix_array.get(idx).expect("in-bounds backward entry");

                if entry.kind() == AffixKind::Suffix || entry.read_index() == read_idx {
                    backward -= 1;
                    continue;
                }

                if entry.kind() != AffixKind::Prefix {
                    break;
                }

                if let Some((diff, score, overlap_len)) =
                    normalised_damerau_levenshtein_distance(suffix_seq, entry.affix(reads))
                {
                    if diff > config.max_diff {
                        up_cut = true;
                        break;
                    }

                    let dst = entry.read_index();
                    let prev = best_scores.get(&dst).copied().unwrap_or(0);
                    if score > prev {
                        best_scores.insert(dst, score);
                        best_spans.insert(dst, overlap_len);
                    }
                }

                backward -= 1;
            }

            if up_cut && down_cut {
                break;
            }
        }

        // Commit best edges for this source row to the sparse matrices (sorted by dst)
        let mut pairs: Vec<(usize, usize)> = best_scores.into_iter().collect();
        pairs.sort_by_key(|(dst, _)| *dst);
        for (dst, score) in pairs.into_iter() {
            let overlap_len = best_spans.get(&dst).copied().unwrap_or(0);
            a_indices.push(dst);
            a_data.push(score);
            o_indices.push(dst);
            o_data.push(overlap_len);
        }
        a_indptr.push(a_indices.len());
        o_indptr.push(o_indices.len());
    }

    let adjacency = CsMat::new((n, n), a_indptr, a_indices, a_data);
    let overlaps = CsMat::new((n, n), o_indptr, o_indices, o_data);

    (affix_array, adjacency, overlaps)
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
        assert_eq!(score, 3);
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
