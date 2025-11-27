/// Given a square list and adjacency matrix, select swap configurations to maximise first subdiagonal sum.
/// Greedily applies swaps for non-overlapping squares; for overlapping, a global optimisation is needed.
pub fn optimise_first_subdiagonal_sum(matrix: &mut CsMat<usize>, squares: &[SwapSquare]) {
    use crate::matching::swap_rows_and_cols;
    use std::collections::HashSet;
    let mut swapped: HashSet<(usize, usize)> = HashSet::new();
    let shape = matrix.shape();
    let mut tri = sprs::TriMat::new(shape);
    for (row_idx, row_vec) in matrix.outer_iterator().enumerate() {
        for (col_idx, &value) in row_vec.iter() {
            tri.add_triplet(row_idx, col_idx, value);
        }
    }
    for sq in squares {
        let i = sq.i;
        let j = sq.j;
        if swapped.contains(&(i, j)) || swapped.contains(&(j, i)) {
            continue;
        }
        let val_ii = matrix.get(i, i).copied().unwrap_or(0);
        let val_jj = matrix.get(j, j).copied().unwrap_or(0);
        let val_ij = matrix.get(i, j).copied().unwrap_or(0);
        let val_ji = matrix.get(j, i).copied().unwrap_or(0);
        if (val_ij + val_ji) > (val_ii + val_jj) {
            swap_rows_and_cols(&mut tri, i, j);
            swapped.insert((i, j));
        }
    }
    *matrix = tri.to_csc();
    // For overlapping squares, a global optimisation (e.g. max-weight matching) should be used
}
/// Alternative path and cycle detection via swap-square analysis.
use sprs::CsMat;
use std::collections::{HashMap, HashSet, VecDeque};

#[derive(Debug, Clone)]
pub struct SwapSquare {
    pub i: usize,
    pub j: usize,
    pub delta: f64,
}

#[derive(Debug, Clone)]
pub struct AlternativesAnalysis {
    pub squares: Vec<SwapSquare>,
    pub components: Vec<Vec<usize>>,
    pub cycles: Vec<Vec<usize>>,
    pub chains: Vec<Vec<usize>>,
    pub ambiguity_count: usize,
}

/// Detect ambiguous square submatrices (swap squares) in the overlap graph.
///
/// For each pair of reads (i, j), check if the following predecessor overlaps are non-zero:
///     matrix[i, i-1], matrix[j, j-1], matrix[i, j-1], matrix[j, i-1]
/// This forms a square region in the first subdiagonal, indicating possible alternative paths.
///
/// Returns a list of SwapSquare where delta is the score change if swapped:
///     delta = (matrix[i, j-1] + matrix[j, i-1]) - (matrix[i, i-1] + matrix[j, j-1])
///
/// If score_gap is provided, only return squares with |delta| <= score_gap.
pub fn detect_swap_squares(matrix: &CsMat<usize>, score_gap: Option<f64>) -> Vec<SwapSquare> {
    let (n, _) = matrix.shape();

    // Build fast lookup for non-zero entries
    let mut lookup: HashMap<(usize, usize), usize> = HashMap::new();
    for (row_idx, row_vec) in matrix.outer_iterator().enumerate() {
        for (col_idx, &value) in row_vec.iter() {
            if value > 0 {
                lookup.insert((row_idx, col_idx), value);
            }
        }
    }

    let mut squares = Vec::new();

    // Check all pairs of reads for predecessor-based swap squares
    for i in 1..n {
        for j in (i + 1)..n {
            // Check predecessor overlaps
            let val_ii1 = lookup.get(&(i, i - 1)).copied().unwrap_or(0);
            let val_jj1 = lookup.get(&(j, j - 1)).copied().unwrap_or(0);
            let val_ij1 = lookup.get(&(i, j - 1)).copied().unwrap_or(0);
            let val_ji1 = lookup.get(&(j, i - 1)).copied().unwrap_or(0);

            if val_ii1 > 0 && val_jj1 > 0 && val_ij1 > 0 && val_ji1 > 0 {
                let delta = (val_ij1 + val_ji1) as f64 - (val_ii1 + val_jj1) as f64;

                if score_gap.is_none() || delta.abs() <= score_gap.unwrap() {
                    squares.push(SwapSquare { i, j, delta });
                }
            }
        }
    }

    squares
}

/// Build adjacency graph from swap squares.
///
/// Nodes are row/col indices, edges connect (i,j) if they form a swap square.
pub fn build_swap_graph(squares: &[SwapSquare]) -> HashMap<usize, HashSet<usize>> {
    let mut graph: HashMap<usize, HashSet<usize>> = HashMap::new();

    for square in squares {
        graph
            .entry(square.i)
            .or_insert_with(HashSet::new)
            .insert(square.j);
        graph
            .entry(square.j)
            .or_insert_with(HashSet::new)
            .insert(square.i);
    }

    graph
}

/// Find connected components in the swap graph using BFS.
pub fn find_connected_components(graph: &HashMap<usize, HashSet<usize>>) -> Vec<Vec<usize>> {
    let mut visited: HashSet<usize> = HashSet::new();
    let mut components = Vec::new();

    let all_nodes: Vec<usize> = graph.keys().copied().collect();

    for &start in &all_nodes {
        if visited.contains(&start) {
            continue;
        }

        // BFS to find component
        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited.insert(start);

        while let Some(node) = queue.pop_front() {
            component.push(node);

            if let Some(neighbors) = graph.get(&node) {
                for &neighbor in neighbors {
                    if !visited.contains(&neighbor) {
                        visited.insert(neighbor);
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        component.sort_unstable();
        components.push(component);
    }

    components
}

/// Check if a connected component forms a cycle using DFS.
///
/// A component is a cycle if:
/// - Size >= 3, AND
/// - There exists a back edge (path from any node back to itself)
pub fn is_cycle(component: &[usize], graph: &HashMap<usize, HashSet<usize>>) -> bool {
    if component.len() < 3 {
        return false;
    }

    let component_set: HashSet<usize> = component.iter().copied().collect();
    let mut visited = HashSet::new();
    let mut rec_stack = HashSet::new();

    fn has_cycle_dfs(
        node: usize,
        parent: Option<usize>,
        visited: &mut HashSet<usize>,
        rec_stack: &mut HashSet<usize>,
        graph: &HashMap<usize, HashSet<usize>>,
        component_set: &HashSet<usize>,
    ) -> bool {
        visited.insert(node);
        rec_stack.insert(node);

        if let Some(neighbors) = graph.get(&node) {
            for &neighbor in neighbors {
                if !component_set.contains(&neighbor) {
                    continue;
                }

                if !visited.contains(&neighbor) {
                    if has_cycle_dfs(
                        neighbor,
                        Some(node),
                        visited,
                        rec_stack,
                        graph,
                        component_set,
                    ) {
                        return true;
                    }
                } else if Some(neighbor) != parent && rec_stack.contains(&neighbor) {
                    // Back edge found (not to immediate parent)
                    return true;
                }
            }
        }

        rec_stack.remove(&node);
        false
    }

    has_cycle_dfs(
        component[0],
        None,
        &mut visited,
        &mut rec_stack,
        graph,
        &component_set,
    )
}

/// Analyse alternative assembly paths via swap square detection.
///
/// Returns AlternativesAnalysis containing:
/// - squares: List of detected swap squares
/// - components: Connected components in the swap graph
/// - cycles: Components that form cycles
/// - chains: Components that are linear chains
/// - ambiguity_count: Total number of swappable positions
pub fn analyse_alternatives(matrix: &CsMat<usize>, score_gap: Option<f64>) -> AlternativesAnalysis {
    let squares = detect_swap_squares(matrix, score_gap);

    if squares.is_empty() {
        return AlternativesAnalysis {
            squares: vec![],
            components: vec![],
            cycles: vec![],
            chains: vec![],
            ambiguity_count: 0,
        };
    }

    let graph = build_swap_graph(&squares);
    let components = find_connected_components(&graph);

    let mut cycles = Vec::new();
    let mut chains = Vec::new();

    for component in &components {
        if is_cycle(component, &graph) {
            cycles.push(component.clone());
        } else {
            chains.push(component.clone());
        }
    }

    let ambiguity_count: usize = components.iter().map(|c| c.len()).sum();

    AlternativesAnalysis {
        squares,
        components,
        cycles,
        chains,
        ambiguity_count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    #[test]
    fn detects_swap_square() {
        // Create a matrix with a swap square at (1,2) on the first subdiagonal
        let mut matrix = TriMat::new((4, 4));
        // Predecessor overlaps for i=1, j=2
        matrix.add_triplet(1, 0, 8); // (i, i-1)
        matrix.add_triplet(2, 1, 7); // (j, j-1)
        matrix.add_triplet(1, 1, 6); // (i, j-1)
        matrix.add_triplet(2, 0, 9); // (j, i-1)

        let csr = matrix.to_csr();
        let squares = detect_swap_squares(&csr, None);

        assert_eq!(squares.len(), 1);
        assert_eq!(squares[0].i, 1);
        assert_eq!(squares[0].j, 2);
        // delta = (6 + 9) - (8 + 7) = 15 - 15 = 0
        assert!((squares[0].delta - 0.0).abs() < 0.001);
    }

    #[test]
    fn detects_cycle_in_component() {
        // Create a 3-node cycle: 0 <-> 1 <-> 2 <-> 0
        let squares = vec![
            SwapSquare {
                i: 0,
                j: 1,
                delta: 0.0,
            },
            SwapSquare {
                i: 1,
                j: 2,
                delta: 0.0,
            },
            SwapSquare {
                i: 2,
                j: 0,
                delta: 0.0,
            },
        ];

        let graph = build_swap_graph(&squares);
        let components = find_connected_components(&graph);

        assert_eq!(components.len(), 1);
        assert_eq!(components[0].len(), 3);
        assert!(is_cycle(&components[0], &graph));
    }

    #[test]
    fn detects_chain_not_cycle() {
        // Create a linear chain: 0 <-> 1 <-> 2
        let squares = vec![
            SwapSquare {
                i: 0,
                j: 1,
                delta: 0.0,
            },
            SwapSquare {
                i: 1,
                j: 2,
                delta: 0.0,
            },
        ];

        let graph = build_swap_graph(&squares);
        let components = find_connected_components(&graph);

        assert_eq!(components.len(), 1);
        assert_eq!(components[0].len(), 3);
        assert!(!is_cycle(&components[0], &graph));
    }
}
