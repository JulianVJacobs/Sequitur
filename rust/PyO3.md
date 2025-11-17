# Python bindings for sequitur_rs (PyO3 prototype)

This document explains the minimal steps to build and use the Rust-based Sequitur functions from Python using PyO3/maturin.

Prerequisites (dev container / local):
- Rust toolchain (rustc, cargo)
- Python 3.8+ and pip
- maturin (recommended) — install with `pip install maturin` or `cargo install maturin`

Quick dev workflow

1) Build and install into the active Python environment (development mode):

```bash
# from repo root
cd rust
maturin develop --release
```

`maturin develop` will build a wheel and `pip install` it into the active Python environment.

2) Use from Python / Jupyter

```python
import sequitur_rs as sequitur

reads = ["ACGT", "GTAA"]
result = sequitur.assemble_from_reads(reads)
print(result.best)

analysis = sequitur.analyse_reads(reads)
print("Adjacency rows:", analysis["adjacency"])
```

Notes
- The current prototype exposes `assemble_from_reads(reads: List[str]) -> AssemblyResult` for the best assembly and `analyse_reads(reads: List[str]) -> dict` for overlap inspection (weights, confidences, cycles, and overlap spans).
- The prototype is intentionally minimal; the next steps are:
  - expose the AffixArray and OverlapGraph types as Python classes (opaque Rust-backed objects),
  - return richer `AssemblyResult` with `alternatives` and `confidence`, and
  - use numpy / scipy-compatible buffers for heavy data transfer if needed.

Troubleshooting
- If maturin is not installed, either install it (`pip install maturin`) or build with `cargo build --release` and then create a Python wrapper manually (more work).
- On Linux manylinux wheels you may need to use `maturin build` in CI to produce portable wheels.

CI
- CI should run `cargo test` and then `maturin build` / `maturin develop` to ensure the Python extension builds cleanly before running Python tests.
