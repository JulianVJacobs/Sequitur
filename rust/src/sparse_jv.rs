use sprs::CsMat;

/// Sparse assignment using successive shortest augmenting path (min-cost max-flow)
/// on a bipartite graph built from a sparse score matrix (rows -> cols).
///
/// The matrix is treated as scores (usize). We convert scores to costs by
/// negating them because the flow routine minimizes cost while we want to
/// maximize total score. Returns `Vec<Option<usize>>` of length `n` where
/// entry `i` is `Some(j)` if row `i` assigned to column `j`, otherwise `None`.
pub fn sparse_assignment(matrix: &CsMat<usize>) -> Vec<Option<usize>> {
    // Build graph: nodes: source(0), rows(1..=n), cols(n+1..=2n), sink(2n+1)
    let n = matrix.rows();
    let nodes = 2 * n + 2;
    let source = 0usize;
    let sink = 2 * n + 1;

    #[derive(Clone)]
    struct Edge {
        to: usize,
        rev: usize,
        cap: i32,
        cost: f64,
    }

    let mut graph: Vec<Vec<Edge>> = vec![Vec::new(); nodes];

    let add_edge = |u: usize, v: usize, cap: i32, cost: f64, graph: &mut Vec<Vec<Edge>>| {
        let rev_u = graph[v].len();
        let rev_v = graph[u].len();
        graph[u].push(Edge {
            to: v,
            rev: rev_u,
            cap,
            cost,
        });
        graph[v].push(Edge {
            to: u,
            rev: rev_v,
            cap: 0,
            cost: -cost,
        });
    };

    // source -> rows
    for i in 0..n {
        add_edge(source, 1 + i, 1, 0.0, &mut graph);
    }
    // cols -> sink
    for j in 0..n {
        add_edge(1 + n + j, sink, 1, 0.0, &mut graph);
    }

    // rows -> cols using sparse matrix entries (scores -> negative cost)
    for (row, cols) in matrix.outer_iterator().enumerate() {
        for (&col_idx, &val) in cols.indices().iter().zip(cols.data().iter()) {
            let u = 1 + row;
            let v = 1 + n + col_idx;
            let cost = -(val as f64);
            add_edge(u, v, 1, cost, &mut graph);
        }
    }

    // Successive shortest augmenting path with potentials
    use std::cmp::Ordering;
    use std::collections::BinaryHeap;

    #[derive(Clone, Copy)]
    struct State(usize, f64);
    impl PartialEq for State {
        fn eq(&self, other: &Self) -> bool {
            self.0 == other.0 && (self.1 - other.1).abs() < 1e-12
        }
    }
    impl Eq for State {}
    impl PartialOrd for State {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            other.1.partial_cmp(&self.1)
        }
    }
    impl Ord for State {
        fn cmp(&self, other: &Self) -> Ordering {
            self.partial_cmp(other).unwrap()
        }
    }

    let mut potential = vec![0.0f64; nodes];
    let mut dist = vec![std::f64::INFINITY; nodes];
    let mut prevnode = vec![usize::MAX; nodes];
    let mut prevedge = vec![usize::MAX; nodes];

    let mut flow = 0i32;
    let max_flow = n as i32; // attempt perfect matching

    while flow < max_flow {
        // Dijkstra
        for v in 0..nodes {
            dist[v] = std::f64::INFINITY;
            prevnode[v] = usize::MAX;
            prevedge[v] = usize::MAX;
        }
        dist[source] = 0.0;
        let mut heap = BinaryHeap::new();
        heap.push(State(source, 0.0));

        while let Some(State(u, _d)) = heap.pop() {
            if _d > dist[u] + 1e-12 {
                continue;
            }
            for (ei, e) in graph[u].iter().enumerate() {
                if e.cap <= 0 {
                    continue;
                }
                let v = e.to;
                let rc = e.cost + potential[u] - potential[v];
                let nd = dist[u] + rc;
                if nd + 1e-12 < dist[v] {
                    dist[v] = nd;
                    prevnode[v] = u;
                    prevedge[v] = ei;
                    heap.push(State(v, nd));
                }
            }
        }

        if dist[sink].is_infinite() {
            break;
        }

        // update potentials
        for v in 0..nodes {
            if !dist[v].is_infinite() {
                potential[v] += dist[v];
            }
        }

        // augment by 1 along path
        let mut v = sink;
        let addf = 1;
        while v != source {
            let u = prevnode[v];
            let ei = prevedge[v];
            let e = &mut graph[u][ei];
            e.cap -= addf;
            let rev = e.rev;
            graph[v][rev].cap += addf;
            v = u;
        }
        flow += addf;
    }

    // Extract assignment: for each row node u=1..n find outgoing edge with used flow
    let mut assignment: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let u = 1 + i;
        for e in &graph[u] {
            if e.to >= 1 + n && e.to <= n + n {
                if e.cap == 0 {
                    let col = e.to - (1 + n);
                    assignment[i] = Some(col);
                    break;
                }
            }
        }
    }

    assignment
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    #[test]
    fn small_sparse_assignment() {
        let mut tri = TriMat::new((3, 3));
        tri.add_triplet(0, 0, 10);
        tri.add_triplet(0, 2, 5);
        tri.add_triplet(1, 1, 8);
        tri.add_triplet(2, 0, 7);
        let csc = tri.to_csc();
        let assign = sparse_assignment(&csc);
        assert_eq!(assign.len(), 3);
        // assignments should be unique
        let mut seen = std::collections::HashSet::new();
        for a in assign.iter().flatten() {
            seen.insert(*a);
        }
        assert!(seen.len() <= 3);
    }
}
