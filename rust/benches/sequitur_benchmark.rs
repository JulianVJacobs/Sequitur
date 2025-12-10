use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use sequitur_rs::affix::PrunedAffixTrie;
use sequitur_rs::find_first_subdiagonal_path;
use sequitur_rs::overlap::{
    create_overlap_graph_unified, verify_overlap_from_anchor, OverlapConfig,
};
use std::time::Duration;

/// Generate synthetic reads with controlled overlaps
fn generate_synthetic_reads(n: usize, read_len: usize, overlap_len: usize) -> Vec<String> {
    let mut rng = StdRng::seed_from_u64(42);
    let bases = ['A', 'C', 'G', 'T'];
    let mut reads = Vec::with_capacity(n);

    for _ in 0..n {
        let read: String = (0..read_len).map(|_| bases[rng.gen_range(0..4)]).collect();
        reads.push(read);
    }

    // Ensure some reads have controlled overlaps
    for i in 0..n.min(n / 2) {
        if i + 1 < n {
            let suffix = &reads[i][reads[i].len().saturating_sub(overlap_len)..];
            let prefix_rest = &reads[i + 1][overlap_len..];
            reads[i + 1] = format!("{}{}", suffix, prefix_rest);
        }
    }

    reads
}

fn bench_trie_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("trie_construction");

    for n in [100, 500, 1000] {
        let reads = generate_synthetic_reads(n, 150, 20);
        group.bench_with_input(BenchmarkId::new("build_trie", n), &reads, |b, reads| {
            b.iter(|| PrunedAffixTrie::build(black_box(reads), 3, 0.25));
        });
    }

    group.finish();
}

fn bench_overlap_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("overlap_construction");
    group.measurement_time(Duration::from_secs(10));

    for n in [100, 500, 1000] {
        let reads = generate_synthetic_reads(n, 150, 20);
        let trie = PrunedAffixTrie::build(&reads, 3, 0.25);
        let config = OverlapConfig::default();

        group.bench_with_input(
            BenchmarkId::new("build_overlap_graph", n),
            &reads,
            |b, reads| {
                let trie = &trie;
                b.iter(|| {
                    let candidates = trie.overlap_candidates(black_box(reads));
                    // Build overlap graph from candidates (without trie construction)
                    create_overlap_graph_from_candidates(black_box(reads), candidates, config)
                });
            },
        );
    }

    group.finish();
}

// Helper function to build overlap graph from pre-computed candidates
fn create_overlap_graph_from_candidates(
    reads: &[String],
    candidates: Vec<(usize, usize, String)>,
    config: OverlapConfig,
) -> (sprs::CsMat<usize>, sprs::CsMat<usize>) {
    use std::collections::HashMap;

    if reads.is_empty() {
        let empty1 = sprs::CsMat::<usize>::zero((0, 0));
        let empty2 = sprs::CsMat::<usize>::zero((0, 0));
        return (empty1, empty2);
    }

    let n = reads.len();
    let mut best_scores: HashMap<(usize, usize), f32> = HashMap::new();
    let mut best_spans: HashMap<(usize, usize), usize> = HashMap::new();

    for (suffix_idx, prefix_idx, shared_affix) in &candidates {
        let suffix_read = &reads[*suffix_idx];
        let prefix_read = &reads[*prefix_idx];

        if let Some((score, overlap_len)) =
            verify_overlap_from_anchor(suffix_read, prefix_read, shared_affix, config)
        {
            let key = (*suffix_idx, *prefix_idx);
            let current_best = best_scores.get(&key).copied().unwrap_or(0.0);
            if score > current_best {
                best_scores.insert(key, score);
                best_spans.insert(key, overlap_len);
            }
        }
    }

    // Build sparse matrices
    let mut adj_triplets_mat = sprs::TriMat::new((n, n));
    let mut ovl_triplets_mat = sprs::TriMat::new((n, n));

    for ((i, j), &score) in &best_scores {
        let weight = (score * 100.0).round() as usize;
        adj_triplets_mat.add_triplet(*i, *j, weight);
        if let Some(&span) = best_spans.get(&(*i, *j)) {
            ovl_triplets_mat.add_triplet(*i, *j, span);
        }
    }

    let adj_matrix = adj_triplets_mat.to_csr();
    let ovl_matrix = ovl_triplets_mat.to_csr();

    (adj_matrix, ovl_matrix)
}

fn bench_candidate_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("candidate_generation");
    group.measurement_time(Duration::from_secs(10));

    for n in [100, 500, 1000] {
        let reads = generate_synthetic_reads(n, 150, 20);
        let trie = PrunedAffixTrie::build(&reads, 3, 0.25);

        group.bench_with_input(BenchmarkId::new("trie_candidates", n), &trie, |b, trie| {
            b.iter(|| trie.overlap_candidates(black_box(&reads)));
        });
    }

    group.finish();
}

fn bench_verification(c: &mut Criterion) {
    let mut group = c.benchmark_group("verification");
    group.measurement_time(Duration::from_secs(10));

    let reads = generate_synthetic_reads(500, 150, 20);
    let trie = PrunedAffixTrie::build(&reads, 3, 0.25);
    let candidates = trie.overlap_candidates(&reads);
    let config = OverlapConfig::default();

    // Take first 100 candidates for verification benchmark
    let test_cases: Vec<_> = candidates
        .iter()
        .take(100)
        .map(|(suffix_idx, prefix_idx, shared_affix)| {
            (
                reads[*suffix_idx].clone(),
                reads[*prefix_idx].clone(),
                shared_affix.clone(),
            )
        })
        .collect();

    group.bench_function("anchor_verification", |b| {
        b.iter(|| {
            for (suffix_read, prefix_read, shared_affix) in test_cases.iter() {
                let _ = verify_overlap_from_anchor(
                    black_box(suffix_read),
                    black_box(prefix_read),
                    black_box(shared_affix),
                    config,
                );
            }
        });
    });

    group.finish();
}

fn bench_full_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_pipeline");
    group.measurement_time(Duration::from_secs(15));
    group.sample_size(20);

    for n in [100, 500, 1000] {
        let reads = generate_synthetic_reads(n, 150, 20);
        let cfg = OverlapConfig::default();
        group.bench_with_input(
            BenchmarkId::new("trie_fuzzy_kmer", n),
            &reads,
            |b, reads| {
                b.iter(|| create_overlap_graph_unified(black_box(reads), cfg));
            },
        );
    }

    group.finish();
}

fn bench_read_length_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("read_length_scaling");
    group.measurement_time(Duration::from_secs(10));

    let n = 200;
    for read_len in [100, 200, 500, 1000] {
        let reads = generate_synthetic_reads(n, read_len, 20);
        group.bench_with_input(BenchmarkId::new("length", read_len), &reads, |b, reads| {
            b.iter(|| create_overlap_graph_unified(black_box(reads), OverlapConfig::default()));
        });
    }

    group.finish();
}

fn bench_matching_assembly(c: &mut Criterion) {
    let mut group = c.benchmark_group("matching_assembly");
    group.measurement_time(Duration::from_secs(10));

    for n in [100, 500, 1000] {
        let reads = generate_synthetic_reads(n, 150, 20);
        let config = OverlapConfig::default();

        // Pre-compute overlap graph
        let (adj_matrix, ovl_matrix) = create_overlap_graph_unified(&reads, config);
        let adj_csc = adj_matrix.to_csc();
        let ovl_csc = ovl_matrix.to_csc();
        let read_ids: Vec<String> = (0..reads.len()).map(|i| format!("read{}", i)).collect();

        group.bench_with_input(BenchmarkId::new("path_finding", n), &reads, |b, reads| {
            b.iter(|| {
                find_first_subdiagonal_path(
                    black_box(&adj_csc),
                    black_box(&ovl_csc),
                    black_box(reads),
                    black_box(&read_ids),
                    None,
                )
            });
        });
    }

    group.finish();
}

fn bench_end_to_end(c: &mut Criterion) {
    let mut group = c.benchmark_group("end_to_end");
    group.measurement_time(Duration::from_secs(15));
    group.sample_size(20);

    for n in [100, 500, 1000] {
        let reads = generate_synthetic_reads(n, 150, 20);
        let config = OverlapConfig::default();

        group.bench_with_input(BenchmarkId::new("full_assembly", n), &reads, |b, reads| {
            b.iter(|| {
                // Complete pipeline: trie -> overlap graph -> matching -> assembly
                let (adj_matrix, ovl_matrix) =
                    create_overlap_graph_unified(black_box(reads), config);
                let adj_csc = adj_matrix.to_csc();
                let ovl_csc = ovl_matrix.to_csc();
                let read_ids: Vec<String> = (0..reads.len()).map(|i| format!("read{}", i)).collect();
                find_first_subdiagonal_path(&adj_csc, &ovl_csc, reads, &read_ids, None)
            });
        });
    }

    group.finish();
}

fn bench_error_rate_impact(c: &mut Criterion) {
    let mut group = c.benchmark_group("error_rate_impact");
    group.measurement_time(Duration::from_secs(10));

    let n = 200;
    let reads = generate_synthetic_reads(n, 150, 20);

    for max_diff in [0.05, 0.15, 0.25, 0.35] {
        let config = OverlapConfig {
            max_diff,
            ..Default::default()
        };

        group.bench_with_input(
            BenchmarkId::new("max_diff", (max_diff * 100.0) as u32),
            &config,
            |b, config| {
                b.iter(|| create_overlap_graph_unified(black_box(&reads), *config));
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_trie_construction,
    bench_overlap_construction,
    bench_candidate_generation,
    bench_verification,
    bench_matching_assembly,
    bench_end_to_end,
    bench_full_pipeline,
    bench_read_length_scaling,
    bench_error_rate_impact
);
criterion_main!(benches);
