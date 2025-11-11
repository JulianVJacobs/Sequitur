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

Notes
-----
This is an initial scaffold to get the Rust port started. It intentionally keeps the first check-in small and compiles as-is to allow incremental development and CI integration.
