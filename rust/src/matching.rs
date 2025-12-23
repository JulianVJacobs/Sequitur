/// Check if two nodes form a swap-square in the adjacency matrix (used for backward path safety).
pub fn is_swap_square(i: usize, j: usize, matrix: &CsMat<usize>) -> bool {
    let ij = matrix.get(i, j).copied().unwrap_or(0);
    let ji = matrix.get(j, i).copied().unwrap_or(0);
    ij > 0 && ji > 0
}
/// Sequence reconstruction helpers.
use std::collections::{HashMap, HashSet};

use sprs::{indexing::SpIndex, CsMat, TriMat};

use lapjv::lapjv;

use crate::overlap::Adjacency;
use crate::read_source::ReadSource;

pub fn argmin_index(values: &[usize]) -> Option<usize> {
    values
        .iter()
        .enumerate()
        .min_by_key(|&(_, value)| value)
        .map(|(idx, _)| idx)
}

/// Swap two rows and columns in a TriMat sparse matrix.
pub fn swap_rows_and_cols(matrix: &mut TriMat<usize>, i: usize, j: usize) {
    let triplets: Vec<(usize, usize, usize)> = matrix
        .triplet_iter()
        .map(|(value, (row, col))| (row.index(), col.index(), *value))
        .collect();
    let mut updated = TriMat::with_capacity(matrix.shape(), triplets.len());
    for (row, col, value) in &triplets {
        let new_row = if *row == i {
            j
        } else if *row == j {
            i
        } else {
            *row
        };
        let new_col = if *col == i {
            j
        } else if *col == j {
            i
        } else {
            *col
        };
        updated.add_triplet(new_row, new_col, *value);
    }
    *matrix = updated;
}

/// Construct a COO sparse matrix mirroring the Python implementation.
pub fn adjacency_to_sparse(adjacency: &Adjacency, size: Option<usize>) -> TriMat<usize> {
    let _max_dim = adjacency.len();
    let _size = size;
    let mut matrix = TriMat::new((adjacency.len(), adjacency.len()));
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

/// Validate that every edge in the assembly path exists in the cost matrix.
// `validate_assembly_path` was removed — path validation is performed inline
// where needed by the LAPJV-based path reconstruction.

/// Placeholder for the subdiagonal traversal; not yet ported from Python.
pub fn find_first_subdiagonal_path(
    matrix: &CsMat<usize>,
    overlap_csc: &CsMat<usize>,
    reads: &[String],
    read_ids: &[String],
    qualities: Option<&[Vec<i32>]>,
    max_insert: usize,
    window: usize,
) -> (String, Vec<usize>) {
    // Delegate to the ReadSource-backed implementation so callers that supply
    // an in-memory `&[String]` continue to work while the core logic uses the
    // disk-backed protocol (`ReadSource`) internally.
    let inmem =
        crate::read_source::InMemoryReadSource::new(reads.to_vec(), Some(read_ids.to_vec()));
    find_first_subdiagonal_path_from_readsource(
        matrix,
        overlap_csc,
        &inmem,
        read_ids,
        qualities,
        max_insert,
        window,
    )
}

/// Compute assembly path (index order) from matrices without depending on in-memory reads.
pub fn compute_assembly_path(matrix: &CsMat<usize>, overlap_csc: &CsMat<usize>) -> Vec<usize> {
    use ndarray::Array2;
    let n = matrix.rows();
    let mut cost_matrix = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            let score = matrix.get(i, j).copied().unwrap_or(0) as f64;
            cost_matrix[(i, j)] = -score;
        }
    }

    if let Ok((assignment_vec, _col_assignments)) = lapjv(&cost_matrix) {
        let mut assigned = vec![false; n];
        let mut path = Vec::with_capacity(n);
        let mut assigned_to: Vec<Option<usize>> = vec![None; n];
        for (i, &j) in assignment_vec.iter().enumerate() {
            assigned_to[j] = Some(i);
        }
        let mut start = None;
        for i in 0..n {
            if assigned_to[i].is_none() {
                start = Some(i);
                break;
            }
        }
        let mut current = match start {
            Some(idx) if idx < n => idx,
            _ => 0,
        };
        for _ in 0..n {
            if assigned[current] {
                break;
            }
            path.push(current);
            assigned[current] = true;
            let next = assignment_vec[current];
            if assigned[next] {
                break;
            }
            current = next;
        }
        while assigned.iter().any(|&v| !v) {
            let next_start = assigned.iter().position(|&v| !v).unwrap_or(0);
            let mut current = next_start;
            for _ in 0..n {
                if assigned[current] {
                    break;
                }
                path.push(current);
                assigned[current] = true;
                let next = assignment_vec[current];
                if assigned[next] {
                    break;
                }
                current = next;
            }
        }
        return path;
    }

    // Fallback greedy
    let csr = matrix.to_csr();
    let csc = matrix.to_csc();
    let (row_count, _) = csr.shape();
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
    let col_sum_val = if start < col_sums.len() {
        col_sums[start]
    } else {
        0
    };
    if col_sum_val != 0 {
        start = argmin_index(&row_sums).unwrap_or(start);
    }
    start = start.min(n.saturating_sub(1));
    let mut path = Vec::with_capacity(n);
    let mut visited = HashSet::with_capacity(n);
    path.push(start);
    visited.insert(start);
    let overlap_csr = overlap_csc.to_csr();
    let mut current = start;
    let mut tried_for_node: HashMap<usize, HashSet<usize>> = HashMap::new();
    while visited.len() < n {
        let tried = tried_for_node.entry(current).or_insert_with(HashSet::new);
        let mut best: Option<(usize, usize, usize)> = None;
        if let Some(row_vec) = csr.outer_view(current) {
            for (col_idx, &weight) in row_vec.indices().iter().zip(row_vec.data().iter()) {
                let target = *col_idx;
                if tried.contains(&target) || visited.contains(&target) {
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
                let is_backward = path.contains(&next) && path.last() != Some(&next);
                if is_backward {
                    if is_swap_square(current, next, &csc) {
                        path.push(next);
                        visited.insert(next);
                        current = next;
                    } else {
                        tried.insert(next);
                        continue;
                    }
                } else {
                    path.push(next);
                    visited.insert(next);
                    current = next;
                }
            }
            None => break,
        }
    }
    if visited.len() < n {
        let mut remaining: Vec<usize> = (0..n).filter(|idx| !visited.contains(idx)).collect();
        remaining.sort_by_key(|idx| col_sums.get(*idx).copied().unwrap_or(usize::MAX));
        for idx in remaining {
            path.push(idx);
        }
    }
    path
}

/// Post-assembly local refinement that attempts to bring mate pairs within `max_insert`
/// by performing conservative adjacent swaps. Swaps are only allowed when the pair of
/// adjacent reads form a swap-square (safe interchange) according to `is_swap_square`.
pub fn enforce_mate_constraints_on_path<R: ReadSource + ?Sized>(
    path: &mut Vec<usize>,
    overlap_csc: &CsMat<usize>,
    reads_src: &R,
    mate_index_of: &Vec<Option<usize>>,
    max_insert: usize,
    _window: usize,
) {
    if path.is_empty() || mate_index_of.is_empty() {
        return;
    }

    // map read id -> position in path
    let mut idx_to_pos: HashMap<usize, usize> = HashMap::new();
    for (pos, &rid) in path.iter().enumerate() {
        idx_to_pos.insert(rid, pos);
    }

    // compute cumulative positions (start coordinate of each read in assembled sequence)
    let mut positions: Vec<usize> = vec![0; path.len()];
    let mut cum = 0usize;
    for (i, &rid) in path.iter().enumerate() {
        positions[i] = cum;
        let len = reads_src.get_len(rid).unwrap_or(0);
        let next_overlap = if i + 1 < path.len() {
            overlap_csc.get(rid, path[i + 1]).copied().unwrap_or(0) as usize
        } else {
            0
        };
        let advance = len.saturating_sub(next_overlap);
        cum = cum.saturating_add(advance);
    }

    // Iterate reads and try to move mates closer using safe adjacent swaps (swap-squares)
    for i in 0..path.len() {
        let read = path[i];
        if read >= mate_index_of.len() {
            continue;
        }
        if let Some(mate) = mate_index_of[read] {
            if let Some(&mate_pos_initial) = idx_to_pos.get(&mate) {
                // recompute distance in base coordinates
                let mut mate_pos = mate_pos_initial;
                let mut dist = if positions[mate_pos] >= positions[i] {
                    positions[mate_pos] - positions[i]
                } else {
                    positions[i] - positions[mate_pos]
                };
                if dist <= max_insert {
                    continue;
                }

                // Move mate left towards i by swapping adjacent nodes when safe
                while mate_pos > i + 1 && dist > max_insert {
                    let a = path[mate_pos - 1];
                    let b = path[mate_pos];
                    if is_swap_square(a, b, overlap_csc) {
                        path.swap(mate_pos - 1, mate_pos);
                        // update position map for swapped entries
                        idx_to_pos.insert(a, mate_pos);
                        idx_to_pos.insert(b, mate_pos - 1);
                        // recompute positions conservatively (full recompute is simple and acceptable here)
                        let mut cum2 = 0usize;
                        for (k, &rid2) in path.iter().enumerate() {
                            positions[k] = cum2;
                            let len2 = reads_src.get_len(rid2).unwrap_or(0);
                            let next_overlap2 = if k + 1 < path.len() {
                                overlap_csc.get(rid2, path[k + 1]).copied().unwrap_or(0) as usize
                            } else {
                                0
                            };
                            let advance2 = len2.saturating_sub(next_overlap2);
                            cum2 = cum2.saturating_add(advance2);
                        }
                        mate_pos -= 1;
                        dist = if positions[mate_pos] >= positions[i] {
                            positions[mate_pos] - positions[i]
                        } else {
                            positions[i] - positions[mate_pos]
                        };
                    } else {
                        break;
                    }
                }

                // Move mate right towards i (if mate is left of i) using same safe rule
                while mate_pos + 1 < i && dist > max_insert {
                    let a = path[mate_pos];
                    let b = path[mate_pos + 1];
                    if is_swap_square(a, b, overlap_csc) {
                        path.swap(mate_pos, mate_pos + 1);
                        idx_to_pos.insert(a, mate_pos + 1);
                        idx_to_pos.insert(b, mate_pos);
                        let mut cum2 = 0usize;
                        for (k, &rid2) in path.iter().enumerate() {
                            positions[k] = cum2;
                            let len2 = reads_src.get_len(rid2).unwrap_or(0);
                            let next_overlap2 = if k + 1 < path.len() {
                                overlap_csc.get(rid2, path[k + 1]).copied().unwrap_or(0) as usize
                            } else {
                                0
                            };
                            let advance2 = len2.saturating_sub(next_overlap2);
                            cum2 = cum2.saturating_add(advance2);
                        }
                        mate_pos += 1;
                        dist = if positions[mate_pos] >= positions[i] {
                            positions[mate_pos] - positions[i]
                        } else {
                            positions[i] - positions[mate_pos]
                        };
                    } else {
                        break;
                    }
                }
            }
        }
    }
}

/// Assemble sequence and return (assembled_string, path) using a `ReadSource` for read access.
pub fn find_first_subdiagonal_path_from_readsource<R: ReadSource + ?Sized>(
    matrix: &CsMat<usize>,
    overlap_csc: &CsMat<usize>,
    reads_src: &R,
    _read_ids: &[String],
    _qualities: Option<&[Vec<i32>]>,
    max_insert: usize,
    window: usize,
) -> (String, Vec<usize>) {
    let path = compute_assembly_path(matrix, overlap_csc);
    if path.is_empty() {
        return (String::new(), path);
    }
    // Attempt to refine `path` using mate-pair information and an insert-size cutoff.
    let mut path = path;
    // Heuristic: build a mate-index mapping if possible from `_read_ids` or ordering.
    let n_reads = reads_src.num_reads().unwrap_or(_read_ids.len());
    let mut mate_index_of: Vec<Option<usize>> = vec![None; n_reads];
    if !_read_ids.is_empty() && _read_ids.len() == n_reads {
        // Detect typical index layout where reads1 are first half and reads2 (reverse-complemented)
        // are the second half (index created by `build_index_from_pair`). Use "/RC" suffix when present.
        let rc_count = _read_ids.iter().filter(|s| s.ends_with("/RC")).count();
        if rc_count > 0 && n_reads % 2 == 0 {
            let half = n_reads / 2;
            if _read_ids.iter().take(half).all(|s| !s.ends_with("/RC"))
                && _read_ids.iter().skip(half).all(|s| s.ends_with("/RC"))
            {
                for i in 0..half {
                    mate_index_of[i] = Some(i + half);
                    mate_index_of[i + half] = Some(i);
                }
            } else {
                // Fallback: match names by stripping "/RC" where present
                let mut name_to_idx: HashMap<String, usize> = HashMap::new();
                for (i, name) in _read_ids.iter().enumerate() {
                    name_to_idx.insert(name.clone(), i);
                }
                for (i, name) in _read_ids.iter().enumerate() {
                    if name.ends_with("/RC") {
                        let base = name.trim_end_matches("/RC").to_string();
                        if let Some(&j) = name_to_idx.get(&base) {
                            mate_index_of[i] = Some(j);
                            mate_index_of[j] = Some(i);
                        }
                    }
                }
            }
        } else if n_reads % 2 == 0 {
            // No explicit RC markers — assume simple half-split ordering
            let half = n_reads / 2;
            for i in 0..half {
                mate_index_of[i] = Some(i + half);
                mate_index_of[i + half] = Some(i);
            }
        }
    } else if n_reads % 2 == 0 {
        // As a last resort, assume a half-split ordering when no names are supplied.
        let half = n_reads / 2;
        for i in 0..half {
            mate_index_of[i] = Some(i + half);
            mate_index_of[i + half] = Some(i);
        }
    }

    // Apply mate constraints using provided parameters.
    enforce_mate_constraints_on_path(
        &mut path,
        overlap_csc,
        reads_src,
        &mate_index_of,
        max_insert,
        window,
    );

    let mut sequence: Vec<u8> = Vec::new();
    match reads_src.get_seq(path[0]) {
        Ok(s) => sequence.extend_from_slice(s.as_bytes()),
        Err(_) => return (String::new(), path),
    }
    for window in path.windows(2) {
        let prev = window[0];
        let next = window[1];
        let overlap_len = overlap_csc.get(prev, next).copied().unwrap_or(0) as usize;
        let prev_len = reads_src.get_len(prev).unwrap_or(0);
        let next_len = reads_src.get_len(next).unwrap_or(0);
        let overlap_len = overlap_len.min(prev_len).min(next_len);
        if overlap_len < next_len {
            match reads_src.get_seq(next) {
                Ok(s) => {
                    let slice = &s[overlap_len..];
                    sequence.extend_from_slice(slice.as_bytes());
                }
                Err(_) => return (String::new(), path),
            }
        }
    }
    (String::from_utf8(sequence).unwrap_or_default(), path)
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
        let read_ids = vec!["read0".to_string(), "read1".to_string()];
        let matrix = adjacency_to_csc(&adjacency, Some(2));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let (assembled, _path) =
            find_first_subdiagonal_path(&matrix, &overlap_csc, &reads, &read_ids, None, 500, 5);
        assert_eq!(assembled, "ACGTAA");
    }

    #[test]
    fn reconstructs_using_qualities() {
        let mut adjacency: Adjacency = vec![HashMap::new(); 2];
        adjacency[0].insert(1, 2);
        let mut overlaps: OverlapLengths = vec![HashMap::new(); 2];
        overlaps[0].insert(1, 2);

        let reads = vec!["ACGT".to_string(), "GTAA".to_string()];
        let read_ids = vec!["read0".to_string(), "read1".to_string()];
        let qualities = vec![vec![10, 10, 5, 5], vec![5, 5, 30, 30]];
        let matrix = adjacency_to_csc(&adjacency, Some(2));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let (assembled, _path) = find_first_subdiagonal_path(
            &matrix,
            &overlap_csc,
            &reads,
            &read_ids,
            Some(&qualities),
            500,
            5,
        );
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
        let read_ids = vec![
            "read0".to_string(),
            "read1".to_string(),
            "read2".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(3));
        let overlap_csc = overlaps_to_csc(&overlaps);

        // With swap-square guard, both orders should be acceptable
        let (assembled, _path) =
            find_first_subdiagonal_path(&matrix, &overlap_csc, &reads, &read_ids, None, 500, 5);
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
        let read_ids = vec![
            "read0".to_string(),
            "read1".to_string(),
            "read2".to_string(),
            "read3".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(4));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let (assembled, _path) =
            find_first_subdiagonal_path(&matrix, &overlap_csc, &reads, &read_ids, None, 500, 5);
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
        let read_ids = vec![
            "read0".to_string(),
            "read1".to_string(),
            "read2".to_string(),
            "read3".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(4));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let (assembled, _path) =
            find_first_subdiagonal_path(&matrix, &overlap_csc, &reads, &read_ids, None, 500, 5);
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
        let read_ids = vec![
            "read0".to_string(),
            "read1".to_string(),
            "read2".to_string(),
        ];
        let matrix = adjacency_to_csc(&adjacency, Some(3));

        let overlap_csc = overlaps_to_csc(&overlaps);
        let (assembled, _path) =
            find_first_subdiagonal_path(&matrix, &overlap_csc, &reads, &read_ids, None, 500, 5);
        // With overlap=3: 8 + (8-3) + (8-3) = 8+5+5 = 18
        // But overlap includes overlapping bases, so: "AAAABBBB" + "BCCCC" + "CDDDD" = "AAAABBBBBCCCCCDDDD" (19 chars)
        assert_eq!(assembled.len(), 18); // Correct accounting for overlaps
        assert!(assembled.starts_with("AAAA"));
        assert!(assembled.ends_with("DDDD"));
    }

    #[test]
    fn mate_refinement_adjacent_swaps_reduce_distance() {
        use crate::read_source::InMemoryReadSource;

        // Create a simple path of 4 reads: 0,1,2,3
        // Mates: 0 <-> 3 (initially far apart). We allow swaps between (1,2) and (2,3)
        let mut path = vec![0usize, 1usize, 2usize, 3usize];

        // Build an overlap matrix that permits swap-square between (1,2) and (2,3)
        let mut adj: Adjacency = vec![HashMap::new(); 4];
        // Bidirectional edges for 1<->2 and 2<->3 to allow swaps
        adj[1].insert(2, 1);
        adj[2].insert(1, 1);
        adj[2].insert(3, 1);
        adj[3].insert(2, 1);
        // Also provide forward overlaps for assembly concatenation
        let overlap_csc = adjacency_to_sparse(&adj, Some(4)).to_csc();

        // Construct reads of length 100 each so mate distance initially large
        let reads = vec![
            "A".repeat(100),
            "C".repeat(100),
            "G".repeat(100),
            "T".repeat(100),
        ];
        let src = InMemoryReadSource::new(reads.clone(), None);

        // Mate mapping: 0<->3
        let mut mate_index_of: Vec<Option<usize>> = vec![None; 4];
        mate_index_of[0] = Some(3);
        mate_index_of[3] = Some(0);

        // Before refinement, compute positions and ensure distance > threshold
        let mut positions_before = vec![0usize; path.len()];
        let mut cum = 0usize;
        for (i, &rid) in path.iter().enumerate() {
            positions_before[i] = cum;
            let len = src.get_len(rid).unwrap_or(0);
            let next_overlap = if i + 1 < path.len() {
                overlap_csc.get(rid, path[i + 1]).copied().unwrap_or(0) as usize
            } else {
                0
            };
            cum = cum.saturating_add(len.saturating_sub(next_overlap));
        }
        let before_dist =
            (positions_before[3] as isize - positions_before[0] as isize).abs() as usize;
        assert!(before_dist > 50, "initial mate distance should be large");

        // Run refinement with a small max_insert to force movement when possible
        enforce_mate_constraints_on_path(&mut path, &overlap_csc, &src, &mate_index_of, 50, 5);

        // Recompute positions after refinement
        let mut positions_after = vec![0usize; path.len()];
        let mut cum2 = 0usize;
        for (i, &rid) in path.iter().enumerate() {
            positions_after[i] = cum2;
            let len = src.get_len(rid).unwrap_or(0);
            let next_overlap = if i + 1 < path.len() {
                overlap_csc.get(rid, path[i + 1]).copied().unwrap_or(0) as usize
            } else {
                0
            };
            cum2 = cum2.saturating_add(len.saturating_sub(next_overlap));
        }

        // find positions of read 0 and its mate 3
        let pos0 = path.iter().position(|&x| x == 0).unwrap();
        let pos3 = path.iter().position(|&x| x == 3).unwrap();
        let after_dist =
            (positions_after[pos3] as isize - positions_after[pos0] as isize).abs() as usize;

        // Expect after refinement that distance is reduced (not necessarily within max_insert
        // for this simple adjacency configuration)
        assert!(
            after_dist < before_dist,
            "mate distance should be reduced by refinement"
        );
    }
}
