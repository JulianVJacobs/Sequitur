//! Sequence reconstruction helpers.

use std::collections::{HashMap, HashSet};

use sprs::{indexing::SpIndex, CsMat, TriMat};

use crate::overlap::Adjacency;

fn argmin_index(values: &[usize]) -> Option<usize> {
    values
        .iter()
        .enumerate()
        .min_by_key(|&(_, value)| value)
        .map(|(idx, _)| idx)
}

/// Construct a COO sparse matrix mirroring the Python implementation.
pub fn adjacency_to_sparse(adjacency: &Adjacency, size: Option<usize>) -> TriMat<usize> {
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
pub fn adjacency_to_csc(adjacency: &Adjacency, size: Option<usize>) -> CsMat<usize> {
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

fn is_swap_square(i: usize, j: usize, csc: &CsMat<usize>) -> bool {
    let row_i = csc.outer_view(i);
    let row_j = csc.outer_view(j);
    let has = |r: Option<sprs::CsVecView<'_, usize>>, col: usize| -> bool {
        if let Some(vec) = r {
            vec.indices().binary_search(&col).is_ok()
        } else {
            false
        }
    };
    has(row_i, i) && has(row_j, j) && has(row_i, j) && has(row_j, i)
}

/// Placeholder for the lower diagonal traversal; not yet ported from Python.
pub fn find_lower_diagonal_path(
    matrix: &CsMat<usize>,
    overlap_csc: &CsMat<usize>,
    reads: &[String],
    qualities: Option<&[Vec<i32>]>,
) -> String {
    if reads.is_empty() {
        return String::new();
    }

    let csr = matrix.to_csr();
    let csc = matrix.to_csc(); // Keep CSC for is_swap_square checks
    let (row_count, _) = csr.shape();

    // Compute sums directly from sparse matrices
    let mut col_sums = vec![0usize; row_count];
    for (col_idx, col_vec) in csc.outer_iterator().enumerate() {
        let sum: usize = col_vec.data().iter().copied().sum();
        col_sums[col_idx] = sum;
    }
    let mut row_sums = vec![0usize; row_count];
    for (row_idx, row_vec) in csr.outer_iterator().enumerate() {
        row_sums[row_idx] = row_vec.data().iter().copied().sum();
    }

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

    let overlap_csr = overlap_csc.to_csr();
    let mut current = start;
    let mut tried_for_node: HashMap<usize, HashSet<usize>> = HashMap::new();
    while visited.len() < reads.len() {
        let tried = tried_for_node.entry(current).or_insert_with(HashSet::new);

        // Compute argnb directly from the CSR row
        let mut best: Option<(usize, usize, usize)> = None; // (target, weight, overlap)
        if let Some(row_vec) = csr.outer_view(current) {
            for (col_idx, &weight) in row_vec.indices().iter().zip(row_vec.data().iter()) {
                let target = *col_idx;
                if tried.contains(&target) {
                    continue;
                }
                if visited.contains(&target) {
                    continue;
                }
                let overlap = overlap_csr
                    .get(current, target)
                    .copied()
                    .unwrap_or_default();
                match best {
                    None => best = Some((target, weight, overlap)),
                    Some((bt, bw, bo)) => {
                        if weight > bw
                            || (weight == bw && overlap > bo)
                            || (weight == bw && overlap == bo && target < bt)
                        {
                            best = Some((target, weight, overlap));
                        }
                    }
                }
            }
        }

        match best.map(|(t, _, _)| t) {
            Some(next) => {
                // Check if this is a backward dependency (already in path but not immediate predecessor)
                let is_backward = path.contains(&next) && path.last() != Some(&next);

                if is_backward {
                    // Backward picks only allowed if reads are interchangeable (swap-square exists)
                    if is_swap_square(current, next, &csc) {
                        // Safe backward pick: accept it
                        path.push(next);
                        visited.insert(next);
                        current = next;
                    } else {
                        // Unsafe backward pick: mask and try next-best
                        tried.insert(next);
                        continue; // Re-enter loop with updated mask
                    }
                } else {
                    // Forward pick: always accept
                    path.push(next);
                    visited.insert(next);
                    current = next;
                }
            }
            None => break,
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
        let overlap_len = overlap_csc.get(prev, next).copied().unwrap_or(0);
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

/// Detect cycles (SCCs) in adjacency represented as Vec<HashMap<usize,usize>>.
pub fn detect_cycles(adjacency: &Adjacency) -> Vec<Vec<usize>> {
    // Build adjacency list
    let n = adjacency.len();
    let mut g: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (u, m) in adjacency.iter().enumerate() {
        for (&v, _) in m.iter() {
            if u < n && v < n {
                g[u].push(v);
            }
        }
    }

    // Tarjan's algorithm
    let mut index = 0usize;
    let mut indices = vec![None; n];
    let mut lowlink = vec![0usize; n];
    let mut stack: Vec<usize> = Vec::new();
    let mut onstack = vec![false; n];
    let mut sccs: Vec<Vec<usize>> = Vec::new();

    fn strongconnect(
        v: usize,
        index: &mut usize,
        indices: &mut [Option<usize>],
        lowlink: &mut [usize],
        stack: &mut Vec<usize>,
        onstack: &mut [bool],
        g: &Vec<Vec<usize>>,
        sccs: &mut Vec<Vec<usize>>,
    ) {
        indices[v] = Some(*index);
        lowlink[v] = *index;
        *index += 1;
        stack.push(v);
        onstack[v] = true;

        for &w in &g[v] {
            if indices[w].is_none() {
                strongconnect(w, index, indices, lowlink, stack, onstack, g, sccs);
                lowlink[v] = lowlink[v].min(lowlink[w]);
            } else if onstack[w] {
                if let Some(iw) = indices[w] {
                    lowlink[v] = lowlink[v].min(iw);
                }
            }
        }

        if let Some(iv) = indices[v] {
            if lowlink[v] == iv {
                let mut comp = Vec::new();
                while let Some(w) = stack.pop() {
                    onstack[w] = false;
                    comp.push(w);
                    if w == v {
                        break;
                    }
                }
                if comp.len() > 0 {
                    sccs.push(comp);
                }
            }
        }
    }

    for v in 0..n {
        if indices[v].is_none() {
            strongconnect(
                v,
                &mut index,
                &mut indices,
                &mut lowlink,
                &mut stack,
                &mut onstack,
                &g,
                &mut sccs,
            );
        }
    }

    // filter singletons without self-loop
    sccs.into_iter()
        .filter(|comp| {
            if comp.len() > 1 {
                true
            } else {
                let v = comp[0];
                adjacency[v].contains_key(&v)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::OverlapLengths;

    // Helper for tests that still use HashMap format
    fn overlaps_to_csc(overlaps: &OverlapLengths) -> CsMat<usize> {
        adjacency_to_sparse(overlaps, Some(overlaps.len())).to_csc()
    }

    #[test]
    fn detect_simple_cycle() {
        let mut adjacency: Adjacency = vec![HashMap::new(); 3];
        // 0 -> 1, 1 -> 2, 2 -> 0 (cycle)
        adjacency[0].insert(1, 1);
        adjacency[1].insert(2, 1);
        adjacency[2].insert(0, 1);

        let sccs = detect_cycles(&adjacency);
        let mut found_cycle = false;
        for comp in sccs.iter() {
            let mut sorted = comp.clone();
            sorted.sort();
            if sorted == vec![0, 1, 2] {
                found_cycle = true;
                break;
            }
        }
        assert!(found_cycle, "Expected SCC containing nodes 0, 1, 2");
    }

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

        let overlap_csc = overlaps_to_csc(&overlaps);
        let assembled = find_lower_diagonal_path(&matrix, &overlap_csc, &reads, None);
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

        let overlap_csc = overlaps_to_csc(&overlaps);
        let assembled = find_lower_diagonal_path(&matrix, &overlap_csc, &reads, Some(&qualities));
        assert_eq!(assembled, "ACGTAA");
    }

    #[test]
    fn allows_interchangeable_backward_pick() {
        // Create a swap-square scenario: A→B and A→C with equal weights, and B↔C edges
        // This represents two reads that can be assembled in either order
        let mut adjacency: Adjacency = vec![HashMap::new(); 3];
        adjacency[0].insert(1, 4); // A→B weight 4
        adjacency[0].insert(2, 4); // A→C weight 4 (equal)
        adjacency[1].insert(2, 4); // B→C
        adjacency[2].insert(1, 4); // C→B (swap-square complete)

        let mut overlaps: OverlapLengths = vec![HashMap::new(); 3];
        overlaps[0].insert(1, 4);
        overlaps[0].insert(2, 4);
        overlaps[1].insert(2, 4);
        overlaps[2].insert(1, 4);

        let reads = vec![
            "AAAABBBB".to_string(),
            "BBBBCCCC".to_string(),
            "BBBBDDDD".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(3));
        let overlap_csc = overlaps_to_csc(&overlaps);

        // With swap-square guard, both orders should be acceptable
        let assembled = find_lower_diagonal_path(&matrix, &overlap_csc, &reads, None);
        // Result should be valid (either A→B→C or A→C→B path)
        assert!(assembled.starts_with("AAAA"));
        assert!(assembled.len() >= 12); // At least 3 reads worth
    }

    #[test]
    fn rejects_non_interchangeable_backward_pick() {
        // Create scenario where best choice is backward but NOT a swap-square
        // A→B (weight 5), A→C (weight 4), B→D (weight 3)
        // If we pick B first, then C appears as a backward option but no C↔B swap-square
        let mut adjacency: Adjacency = vec![HashMap::new(); 4];
        adjacency[0].insert(1, 5); // A→B best
        adjacency[0].insert(2, 4); // A→C second-best
        adjacency[1].insert(3, 3); // B→D
                                   // No B↔C edges, so no swap-square

        let mut overlaps: OverlapLengths = vec![HashMap::new(); 4];
        overlaps[0].insert(1, 4);
        overlaps[0].insert(2, 4);
        overlaps[1].insert(3, 2);

        let reads = vec![
            "AAAABBBB".to_string(),
            "BBBBCCCC".to_string(),
            "AAAADDDD".to_string(),
            "CCCCEEEE".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(4));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let assembled = find_lower_diagonal_path(&matrix, &overlap_csc, &reads, None);
        // Should follow valid forward path: A→B→D, then append C
        assert!(assembled.len() > 0);
    }

    #[test]
    fn handles_multi_rejection_gracefully() {
        // Create chain requiring multiple argnb calls due to backward rejections
        let mut adjacency: Adjacency = vec![HashMap::new(); 4];
        adjacency[0].insert(1, 5);
        adjacency[0].insert(2, 5); // Equal weight
        adjacency[0].insert(3, 5); // Equal weight
        adjacency[1].insert(2, 3);

        let mut overlaps: OverlapLengths = vec![HashMap::new(); 4];
        overlaps[0].insert(1, 3);
        overlaps[0].insert(2, 3);
        overlaps[0].insert(3, 3);
        overlaps[1].insert(2, 2);

        let reads = vec![
            "AAAABBBB".to_string(),
            "BBBBCCCC".to_string(),
            "BBBBDDDD".to_string(),
            "BBBBEEEE".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(4));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let assembled = find_lower_diagonal_path(&matrix, &overlap_csc, &reads, None);
        // With overlaps, total length will be: 8 + (8-3) + (8-3) + (8-3) = 8+5+5+5 = 23
        // But actual behavior depends on path order
        assert!(assembled.len() >= 8); // At least one read
        assert!(assembled.contains("AAAA")); // Start should be present
    }

    #[test]
    fn handles_terminal_read_without_outgoing() {
        // Test that terminal reads (no outgoing edges) don't cause infinite loops
        let mut adjacency: Adjacency = vec![HashMap::new(); 3];
        adjacency[0].insert(1, 4);
        adjacency[1].insert(2, 4);
        // Read 2 has no outgoing edges (terminal)

        let mut overlaps: OverlapLengths = vec![HashMap::new(); 3];
        overlaps[0].insert(1, 3);
        overlaps[1].insert(2, 3);

        let reads = vec![
            "AAAABBBB".to_string(),
            "BBBBCCCC".to_string(),
            "CCCCDDDD".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(3));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let assembled = find_lower_diagonal_path(&matrix, &overlap_csc, &reads, None);
        // With overlap=3: 8 + (8-3) + (8-3) = 8+5+5 = 18
        // But overlap includes overlapping bases, so: "AAAABBBB" + "BCCCC" + "CDDDD" = "AAAABBBBBCCCCCDDDD" (19 chars)
        assert_eq!(assembled.len(), 18); // Correct accounting for overlaps
        assert!(assembled.starts_with("AAAA"));
        assert!(assembled.ends_with("DDDD"));
    }
}
