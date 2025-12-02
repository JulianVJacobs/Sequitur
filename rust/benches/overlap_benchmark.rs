use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use sequitur_rs::{create_overlap_graph, create_overlap_graph_unified, AffixMap, OverlapConfig};
use std::time::Duration;

/// Benchmarks comparing trie vs array-based affix structures.
/// For comprehensive end-to-end benchmarks, see sequitur_benchmark.rs

fn generate_synthetic_reads(n: usize, read_len: usize, overlap_len: usize) -> Vec<String> {
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    let mut rng = StdRng::seed_from_u64(42);
    let bases = b"ACGT";
    let mut reads = Vec::with_capacity(n);

    // Generate first read
    let mut current: Vec<u8> = (0..read_len).map(|_| bases[rng.gen_range(0..4)]).collect();
    reads.push(String::from_utf8(current.clone()).unwrap());

    // Generate subsequent reads with controlled overlap
    for _ in 1..n {
        let mut next = Vec::with_capacity(read_len);
        // Copy overlap from previous read
        next.extend_from_slice(&current[read_len - overlap_len..]);
        // Generate new bases
        for _ in overlap_len..read_len {
            next.push(bases[rng.gen_range(0..4)]);
        }
        reads.push(String::from_utf8(next.clone()).unwrap());
        current = next;
    }

    reads
}

fn bench_trie_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("trie_construction");
    group.measurement_time(Duration::from_secs(10));

    for n in [100, 500, 1000].iter() {
        let reads = generate_synthetic_reads(*n, 150, 20);

        group.bench_with_input(BenchmarkId::new("build_trie", n), &reads, |b, reads| {
            b.iter(|| {
                use sequitur_rs::affix::PrunedAffixTrie;
                let _trie = PrunedAffixTrie::build(black_box(reads), 3, 0.25);
            });
        });
    }

    group.finish();
}

fn bench_array_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("array_construction");
    group.measurement_time(Duration::from_secs(10));

    for n in [100, 500, 1000].iter() {
        let reads = generate_synthetic_reads(*n, 150, 20);

        group.bench_with_input(BenchmarkId::new("build_array", n), &reads, |b, reads| {
            b.iter(|| {
                let _array = AffixMap::build(black_box(reads), 3);
            });
        });
    }

    group.finish();
}

fn bench_overlap_graph_trie(c: &mut Criterion) {
    let mut group = c.benchmark_group("overlap_graph_trie");
    group.measurement_time(Duration::from_secs(15));
    group.sample_size(20);

    for n in [100, 500, 1000].iter() {
        let reads = generate_synthetic_reads(*n, 150, 20);

        group.bench_with_input(BenchmarkId::new("full_pipeline", n), &reads, |b, reads| {
            b.iter(|| {
                let config = OverlapConfig {
                    use_trie: true,
                    max_diff: 0.25,
                    min_suffix_len: 3,
                    ..Default::default()
                };
                let (_adj, _ovl) = create_overlap_graph_unified(black_box(reads), config);
            });
        });
    }

    group.finish();
}

fn bench_overlap_graph_array(c: &mut Criterion) {
    let mut group = c.benchmark_group("overlap_graph_array");
    group.measurement_time(Duration::from_secs(15));
    group.sample_size(20);

    for n in [100, 500].iter() {
        let reads = generate_synthetic_reads(*n, 150, 20);

        group.bench_with_input(BenchmarkId::new("full_pipeline", n), &reads, |b, reads| {
            b.iter(|| {
                let config = OverlapConfig {
                    use_trie: false,
                    max_diff: 0.25,
                    min_suffix_len: 3,
                    ..Default::default()
                };
                let affix_map = AffixMap::build(reads, 3);
                let (_map, _adj, _ovl) =
                    sequitur_rs::create_overlap_graph(black_box(reads), Some(affix_map), config);
            });
        });
    }

    group.finish();
}

// Removed: bench_candidate_generation, bench_verification_strategies,
// bench_different_read_lengths, bench_error_rates
// These are now in sequitur_benchmark.rs for comprehensive analysis

criterion_group!(
    benches,
    bench_trie_construction,
    bench_array_construction,
    bench_overlap_graph_trie,
    bench_overlap_graph_array,
);

criterion_main!(benches);
