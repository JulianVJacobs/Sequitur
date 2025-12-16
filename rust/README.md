Sequitur Rust port
===================

This folder contains the initial scaffold for porting Sequitur's core pipeline to Rust.

Goals
-----
- Provide a low-level, memory-efficient implementation of the Sequitur core algorithms:
  - suffix-array construction (streaming-friendly)
  - normalized Damerau–Levenshtein overlap scoring
  - sparse adjacency representation (coordinate lists / CSR)
  - matching/reconstruction pipeline that enforces at-most-one-predecessor/one-successor

Why Rust?
---------
Rust gives precise control over memory layout and ownership without a GC. That lets us:
- preallocate buffers and reuse them across passes;
- use compact numeric types and tightly-packed sparse structures;
- parallelise safely with Rayon while avoiding data races;
- isolate unsafe code for aggressive performance wins.

How to build (host must have Rust toolchain installed)
-----------------------------------------------------

# build
cargo build --release

# run the placeholder CLI (currently prints a message)
cargo run -- --help

Structure
---------
- Cargo.toml  - dependencies and project metadata
- src/main.rs - CLI entrypoint (placeholder)
- src/lib.rs  - library surface and TODO module skeletons

Next steps
----------
1. Implement suffix array and a lightweight disk-backed index for very large read sets.
2. Implement Damerau–Levenshtein scoring with early exit and memory pooling.
3. Build coordinate-list adjacency representation and conversion to CSR for matching.
4. Implement a memory-efficient matching algorithm (greedy, then refine) with an optional exact solver.

Disk-backed Read Index
----------------------
Sequitur's Rust prototype includes a lightweight on-disk read index to avoid holding all reads in memory. Key points:

- Index files: `<base>.seqs` (concatenated sequence bytes) and `<base>.sidx.json` (JSON array of `{name,offset,len}` entries).
- Indexer: a small binary `index_reads` is available; build it with `cargo build --release` and run `./target/release/index_reads --input <combined.fastq> --output /tmp/index_base`.
- CLI: pass the base path to the index with `--read-index /tmp/index_base` when running `sequitur`.
- Orientation: `sequitur` applies reverse-complement to `reads2` by default. Ensure the index was produced from the same read orientation (recommended: concatenate `reads_1` then `reads_2` without pre-applying RC).

Example
```
cd rust
cargo build --release
# create index from concatenated reads (reads1 followed by reads2)
cat ../tests/synthetic/swap_test/data/reads_1.fastq ../tests/synthetic/swap_test/data/reads_2.fastq > /tmp/swap_combined.fastq
./target/release/index_reads --input /tmp/swap_combined.fastq --output /tmp/swap_index

# run sequitur using the on-disk index
./target/release/sequitur ../tests/synthetic/swap_test/data/reads_1.fastq ../tests/synthetic/swap_test/data/reads_2.fastq \
  --read-index /tmp/swap_index --output-fasta ../tests/synthetic/swap_test/results/assembled_index.fasta
```

Limitations / notes:
- Prototype implementation may materialise reads during trie construction; future work will stream trie construction to avoid materialisation.
- Current index format is intentionally simple (JSON + concatenated sequences); compression-aware formats and a compact binary index are planned.


Notes
-----
This is an initial scaffold to get the Rust port started. It intentionally keeps the first check-in small and compiles as-is to allow incremental development and CI integration.
