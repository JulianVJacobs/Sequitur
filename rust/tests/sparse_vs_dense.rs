use lapjv::lapjv;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use sprs::TriMat;

use sequitur::sparse_assignment;

#[test]
fn compare_sparse_and_dense_small_random() {
    let n = 8usize;
    let mut tri = TriMat::new((n, n));
    let mut rng = StdRng::seed_from_u64(42);
    for i in 0..n {
        for j in 0..n {
            if rng.gen_bool(0.6) {
                let val = rng.gen_range(1..20);
                tri.add_triplet(i, j, val);
            }
        }
    }
    let csc = tri.to_csc();

    // dense lapjv expects cost matrix (f64) with cost = -score
    let mut cost = ndarray::Array2::<f64>::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            let score = csc.get(i, j).copied().unwrap_or(0) as f64;
            cost[(i, j)] = -score;
        }
    }

    let (lap_assign, _) = lapjv(&cost).expect("lapjv should succeed");
    let lap_score: usize = lap_assign
        .iter()
        .enumerate()
        .map(|(i, &j)| csc.get(i, j).copied().unwrap_or(0))
        .sum();

    let sparse_assign = sparse_assignment(&csc);
    let sparse_score: usize = sparse_assign
        .iter()
        .enumerate()
        .map(|(i, o)| o.map(|j| csc.get(i, j).copied().unwrap_or(0)).unwrap_or(0))
        .sum();

    // For this dense-ish random test, scores should match (both find max assignment)
    assert_eq!(lap_score, sparse_score);
}
