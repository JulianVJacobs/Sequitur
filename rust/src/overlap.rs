//! Overlap graph construction primitives translated from the Python reference.

use std::collections::HashMap;

use strsim::damerau_levenshtein;
use log::warn;

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

fn update_edge(
	adjacency: &mut HashMap<usize, usize>,
	overlaps: &mut HashMap<usize, usize>,
	target: usize,
	score: usize,
	overlap_len: usize,
) {
	match adjacency.get(&target) {
		Some(&current) if current > score => return,
		Some(&current) if current == score => {
			if let Some(&existing_overlap) = overlaps.get(&target) {
				if existing_overlap >= overlap_len {
					return;
				}
			}
		}
		_ => {}
	}

	adjacency.insert(target, score);
	overlaps.insert(target, overlap_len);
}

/// Build the weighted overlap graph for the supplied reads.
pub fn create_overlap_graph(
	reads: &[String],
	mut affix_array: Option<AffixArray>,
	config: OverlapConfig,
) -> (AffixArray, Adjacency, OverlapLengths) {
	if config.use_threads {
		warn!("Threaded overlap construction is not yet implemented; falling back to sequential mode");
	}

	let min_suffix_len = config.min_suffix_len.max(1);

	let affix_array = affix_array
		.take()
		.unwrap_or_else(|| AffixArray::build(reads.iter().map(|s| s.as_str()), min_suffix_len));

	if reads.is_empty() {
		return (affix_array, Vec::new(), Vec::new());
	}

	let mut adjacency: Adjacency = vec![HashMap::new(); reads.len()];
	let mut overlap_lengths: OverlapLengths = vec![HashMap::new(); reads.len()];

	for (read_idx, read) in reads.iter().enumerate() {
		if read.len() < min_suffix_len {
			continue;
		}

		let mut up_cut = false;
		let mut down_cut = false;
		let read_len = read.len();

		for span in min_suffix_len..read_len {
			let suffix_seq = &read[read_len - span..];
			let key = format!("{}${}", suffix_seq, read_idx);
			let Some(anchor) = affix_array.index_of(&key) else {
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
					normalised_damerau_levenshtein_distance(suffix_seq, entry.affix())
				{
					if diff > config.max_diff {
						down_cut = true;
						break;
					}

					update_edge(
						&mut adjacency[read_idx],
						&mut overlap_lengths[read_idx],
						entry.read_index(),
						score,
						overlap_len,
					);
				}

				forward += 1;
			}

			let mut backward = anchor as isize - 1;
			while backward >= 0 {
				let idx = backward as usize;
				let entry = affix_array
					.get(idx)
					.expect("in-bounds backward entry");

				if entry.kind() == AffixKind::Suffix || entry.read_index() == read_idx {
					backward -= 1;
					continue;
				}

				if entry.kind() != AffixKind::Prefix {
					break;
				}

				if let Some((diff, score, overlap_len)) =
					normalised_damerau_levenshtein_distance(suffix_seq, entry.affix())
				{
					if diff > config.max_diff {
						up_cut = true;
						break;
					}

					update_edge(
						&mut adjacency[read_idx],
						&mut overlap_lengths[read_idx],
						entry.read_index(),
						score,
						overlap_len,
					);
				}

				backward -= 1;
			}

			if up_cut && down_cut {
				break;
			}
		}
	}

	(affix_array, adjacency, overlap_lengths)
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
			let mut ordered: Vec<(usize, usize)> = edges.iter().map(|(&dst, &weight)| (dst, weight)).collect();
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

		assert_eq!(adjacency.len(), 2);
		assert_eq!(overlaps.len(), 2);

		let edge = adjacency[0].get(&1).copied();
		let span = overlaps[0].get(&1).copied();

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
