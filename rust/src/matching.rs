//! Sequence reconstruction helpers.

use std::collections::{HashMap, HashSet};

use sprs::{indexing::SpIndex, CsMat, TriMat};

use crate::overlap::{Adjacency, OverlapLengths};

fn column_sums(adjacency: &Adjacency) -> Vec<usize> {
	let size = adjacency.len();
	let mut sums = vec![0; size];
	for edges in adjacency {
		for (&col, &weight) in edges {
			sums[col] += weight;
		}
	}
	sums
}

fn row_sums(adjacency: &Adjacency) -> Vec<usize> {
	adjacency
		.iter()
		.map(|edges| edges.values().copied().sum())
		.collect()
}

fn argmin_index(values: &[usize]) -> Option<usize> {
	values
		.iter()
		.enumerate()
		.min_by_key(|&(_, value)| value)
		.map(|(idx, _)| idx)
}

fn select_next_read(
	current: usize,
	visited: &HashSet<usize>,
	adjacency: &Adjacency,
	overlaps: &OverlapLengths,
) -> Option<usize> {
	let mut candidates: Vec<(usize, usize, usize)> = adjacency[current]
		.iter()
		.filter(|(&target, _)| !visited.contains(&target))
		.map(|(&target, &weight)| {
			let overlap = overlaps[current].get(&target).copied().unwrap_or_default();
			(target, weight, overlap)
		})
		.collect();
	if candidates.is_empty() {
		return None;
	}
	candidates.sort_by(|a, b| {
		b.1.cmp(&a.1)
			.then_with(|| b.2.cmp(&a.2))
			.then_with(|| a.0.cmp(&b.0))
	});
	Some(candidates[0].0)
}

/// Construct a COO sparse matrix mirroring the Python implementation.
pub fn adjacency_to_sparse(
	adjacency: &Adjacency,
	size: Option<usize>,
) -> TriMat<usize> {
	let mut max_dim = adjacency.len();
	for (row, cols) in adjacency.iter().enumerate() {
		for (&col, _) in cols {
			max_dim = max_dim.max(row + 1).max(col + 1);
		}
	}

	let dim = size.unwrap_or(max_dim);
	let mut matrix = TriMat::new((dim, dim));

	for (row, cols) in adjacency.iter().enumerate() {
		for (&col, &value) in cols {
			matrix.add_triplet(row, col, value);
		}
	}

	matrix
}

/// Convert the COO matrix into a compressed sparse column representation.
pub fn adjacency_to_csc(
	adjacency: &Adjacency,
	size: Option<usize>,
) -> CsMat<usize> {
	adjacency_to_sparse(adjacency, size).to_csc()
}

/// Relabel column indices in-place according to the provided mapping.
pub fn relabel_columns(matrix: &mut TriMat<usize>, cols_map: &HashMap<usize, usize>) {
	let triplets: Vec<(usize, usize, usize)> = matrix
		.triplet_iter()
		.map(|(value, (row, col))| (row.index(), col.index(), *value))
		.collect();

	let mut updated = TriMat::with_capacity(matrix.shape(), triplets.len());
	for (row, col, value) in triplets {
		let mapped_col = cols_map.get(&col).copied().unwrap_or(col);
		updated.add_triplet(row, mapped_col, value);
	}

	*matrix = updated;
}

/// Relabel row indices in-place according to the provided mapping.
pub fn relabel_rows(matrix: &mut TriMat<usize>, rows_map: &HashMap<usize, usize>) {
	let triplets: Vec<(usize, usize, usize)> = matrix
		.triplet_iter()
		.map(|(value, (row, col))| (row.index(), col.index(), *value))
		.collect();

	let mut updated = TriMat::with_capacity(matrix.shape(), triplets.len());
	for (row, col, value) in triplets {
		let mapped_row = rows_map.get(&row).copied().unwrap_or(row);
		updated.add_triplet(mapped_row, col, value);
	}

	*matrix = updated;
}

/// Placeholder for the lower diagonal traversal; not yet ported from Python.
pub fn find_lower_diagonal_path(
	matrix: &CsMat<usize>,
	overlaps: &OverlapLengths,
	reads: &[String],
	qualities: Option<&[Vec<i32>]>,
) -> String {
	if reads.is_empty() {
		return String::new();
	}

	let csr = matrix.to_csr();
	let (row_count, _) = csr.shape();
	let mut adjacency: Adjacency = vec![HashMap::new(); row_count];
	for (row_idx, row_vec) in csr.outer_iterator().enumerate() {
		for (col_idx, weight) in row_vec.iter() {
			if *weight > 0 {
				adjacency[row_idx].insert(col_idx, *weight);
			}
		}
	}

	let col_sums = column_sums(&adjacency);
	let row_sums = row_sums(&adjacency);
	let mut start = argmin_index(&col_sums).unwrap_or(0);
	if col_sums.get(start).copied().unwrap_or(0) != 0 {
		if let Some(row_idx) = argmin_index(&row_sums) {
			start = row_idx;
		}
	}
	start = start.min(reads.len().saturating_sub(1));

	let mut path = Vec::with_capacity(reads.len());
	let mut visited = HashSet::with_capacity(reads.len());
	path.push(start);
	visited.insert(start);

	let mut current = start;
	while visited.len() < reads.len() {
		if let Some(next) = select_next_read(current, &visited, &adjacency, overlaps) {
			path.push(next);
			visited.insert(next);
			current = next;
		} else {
			break;
		}
	}

	if visited.len() < reads.len() {
		let mut remaining: Vec<usize> = (0..reads.len())
			.filter(|idx| !visited.contains(idx))
			.collect();
		remaining.sort_by_key(|idx| col_sums.get(*idx).copied().unwrap_or(usize::MAX));
		for idx in remaining {
			path.push(idx);
		}
	}

	if path.is_empty() {
		return String::new();
	}

	let mut sequence: Vec<u8> = Vec::new();
	sequence.extend_from_slice(reads[path[0]].as_bytes());

	for window in path.windows(2) {
		let prev = window[0];
		let next = window[1];
		let overlap_len = overlaps[prev].get(&next).copied().unwrap_or(0);
		let overlap_len = overlap_len.min(reads[prev].len()).min(reads[next].len());

		if let Some(all_qualities) = qualities {
			let prev_q = &all_qualities[prev];
			let next_q = &all_qualities[next];
			let prev_start = prev_q.len().saturating_sub(overlap_len);
			for idx in 0..overlap_len {
				let seq_idx = sequence.len().saturating_sub(overlap_len) + idx;
				let prev_quality = prev_q.get(prev_start + idx).copied().unwrap_or(0);
				let next_quality = next_q.get(idx).copied().unwrap_or(0);
				if next_quality > prev_quality {
					if let Some(&base) = reads[next].as_bytes().get(idx) {
						if seq_idx < sequence.len() {
							sequence[seq_idx] = base;
						}
					}
				}
			}
		}

		if overlap_len < reads[next].len() {
			sequence.extend_from_slice(&reads[next].as_bytes()[overlap_len..]);
		}
	}

	String::from_utf8(sequence).expect("assembled sequence should be valid UTF-8")
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn builds_sparse_matrix_from_adjacency() {
		let mut adjacency: Adjacency = vec![HashMap::new(); 3];
		adjacency[0].insert(1, 2);
		adjacency[1].insert(2, 3);

		let matrix = adjacency_to_sparse(&adjacency, None);
		assert_eq!(matrix.nnz(), 2);

		let csc = matrix.to_csc::<usize>();
		assert_eq!(csc.rows(), 3);
		assert_eq!(csc.cols(), 3);
	}

	#[test]
	fn relabels_columns_and_rows() {
		let mut adjacency: Adjacency = vec![HashMap::new(); 2];
		adjacency[0].insert(1, 1);

		let mut matrix = adjacency_to_sparse(&adjacency, None);

		let mut cols_map = HashMap::new();
		cols_map.insert(1, 0);
		relabel_columns(&mut matrix, &cols_map);

		let mut rows_map = HashMap::new();
		rows_map.insert(0, 1);
		relabel_rows(&mut matrix, &rows_map);

		let rows = matrix.row_inds();
		let cols = matrix.col_inds();
		assert_eq!(rows, &[1]);
		assert_eq!(cols, &[0]);
	}

	#[test]
	fn reconstructs_simple_path_without_qualities() {
		let mut adjacency: Adjacency = vec![HashMap::new(); 2];
		adjacency[0].insert(1, 2);
		let mut overlaps: OverlapLengths = vec![HashMap::new(); 2];
		overlaps[0].insert(1, 2);

		let reads = vec!["ACGT".to_string(), "GTAA".to_string()];
		let matrix = adjacency_to_csc(&adjacency, Some(2));

		let assembled = find_lower_diagonal_path(&matrix, &overlaps, &reads, None);
		assert_eq!(assembled, "ACGTAA");
	}

	#[test]
	fn reconstructs_using_qualities() {
		let mut adjacency: Adjacency = vec![HashMap::new(); 2];
		adjacency[0].insert(1, 2);
		let mut overlaps: OverlapLengths = vec![HashMap::new(); 2];
		overlaps[0].insert(1, 2);

		let reads = vec!["ACGT".to_string(), "GTAA".to_string()];
		let qualities = vec![vec![10, 10, 5, 5], vec![5, 5, 30, 30]];
		let matrix = adjacency_to_csc(&adjacency, Some(2));

		let assembled = find_lower_diagonal_path(&matrix, &overlaps, &reads, Some(&qualities));
		assert_eq!(assembled, "ACGTAA");
	}
}
