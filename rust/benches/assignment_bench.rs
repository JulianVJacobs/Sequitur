use criterion::{criterion_group, criterion_main, Criterion};
use lapjv::lapjv;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use sequitur::sparse_assignment;
use sprs::TriMat;

fn make_random_sparse(n: usize, degree: usize, seed: u64) -> sprs::CsMat<usize> {
    let mut tri = TriMat::new((n, n));
    let mut rng = StdRng::seed_from_u64(seed);
    for i in 0..n {
        for _ in 0..degree {
            let j = rng.gen_range(0..n);
            let val = rng.gen_range(1..50);
            tri.add_triplet(i, j, val);
        }
    }
    tri.to_csc()
}

fn bench_assignment(c: &mut Criterion) {
    let n = 200usize;
    let deg = 10usize;
    let mat = make_random_sparse(n, deg, 123);

    c.bench_function("sparse_assignment", |b| {
        b.iter(|| {
            let _ = sparse_assignment(&mat);
        })
    });

    // build dense cost for lapjv
    let mut cost = ndarray::Array2::<f64>::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            let score = mat.get(i, j).copied().unwrap_or(0) as f64;
            cost[(i, j)] = -score;
        }
    }
    c.bench_function("lapjv_dense", |b| {
        b.iter(|| {
            let _ = lapjv(&cost);
        })
    });
}

criterion_group!(benches, bench_assignment);
criterion_main!(benches);
