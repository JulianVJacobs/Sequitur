**Notebook Context**: `python/sequitur.ipynb` and `python/results.ipynb` document algorithmic intent; they may import deprecated Python code—update to use `sequitur` when possible.

**Notebook Rules**:
 - For full instructions and best practices, refer to `python/README.md`.
 
 **Documentation Rules:**
 Refer to the root `README.md` for documentation rules regarding markdown file naming, placement, and indexing. The root README.md serves as the canonical index for all further instructions and rules.
**Sequitur Agent Guide**
**Terminal Command Directory Rule**: Always run shell commands using the pattern `(cd /desired/path && command ...)` to ensure the correct working directory is used. This avoids confusion and errors from manual `cd` changes, and always falls back to the workspace root if not specified. Never rely on the current terminal directory state.
- **Architecture**: Rust implementation in `rust/` is the primary codebase; Python in `python/` is deprecated for core algorithms but maintained for prototyping and notebooks—prefer Rust with PyO3 bindings for production; see `docs/RUST_ARCHITECTURE.md` for rationale and migration guide.
- **Alternative Paths**: `alternative_paths.rs` implements swap-square detection; pure Python version in `sequitur_core/alternative_paths.py` is deprecated—direct new code to use `sequitur.analyse_alternative_paths()`.
- **Documentation**: All feature docs live in `docs/` with `docs/README.md` as index; implementation-specific docs stay in language subdirs (`python/README.md`, `rust/README.md`); update root `README.md` only for quick-start or roadmap changes.
- **Pipeline Flow**: Reads move through `build_suffix_array` → `create_bipartite_adjacency_matrix` → `adjacency_to_sparse` → `find_lower_diagonal_path`; keep these interfaces stable because CLI, experiments, and notebooks rely on them.
- **Suffix Array Contract**: `build_suffix_array` decorates affixes as `suffix${idx}` / `prefix^{idx}`—Rust's `AffixArray` expects the same format; any change must update both languages plus lookup logic.
- **Overlap Construction**: `create_bipartite_adjacency_matrix` filters overlaps by `max_diff` and keeps the highest weight per edge; tightening thresholds or score semantics will ripple into matching correctness tests.
- **Threading Toggle**: Python overlap building optionally uses threads (`--threads`), but deterministic single-threaded runs are the baseline for debugging; avoid introducing non-determinism in scoring.
- **Matching Semantics**: Rust `matching.rs` is canonical. Use `argnb` (next-best masked selection) with swap-square safety (see `alternative_paths.rs`) to avoid swap loops. `find_lower_diagonal_path` assumes square COO/CSC matrices with identical row/col orderings and uses fallback greedy assembly on failure—validate both code paths when editing traversal heuristics.
- **Quality Scores**: When qualities are supplied, overlapping bases are resolved by per-position quality comparisons; maintain array alignment (`quality_map[idx][-overlap_len:]`) when altering overlap math.
- **CLI Entry**: Both `rust/src/main.rs` and `python/sequitur.py` ingest paired FASTQ files, reverse-complement reads2, and optionally write metrics; Rust CLI is primary, Python CLI wraps Rust when available.
- **Python Setup**: Work inside `python/`; create `.venv`, install `numpy scipy biopython`; pure Python core (`sequitur_core/`) is deprecated—including matching—use the Rust library (`sequitur`) via PyO3 for production (`maturin develop`).
- **Experiments**: `python/experiments/natural_language.py` reproduces toy datasets; it depends on the same core functions and is used for regression checks alongside notebooks.
- **Rust Modules**: `suffix.rs`, `overlap.rs`, `matching.rs`, `alternative_paths.rs` are the canonical implementations; keep public APIs stable and exported via `lib.rs` for PyO3 bindings.
- **Rust CLI**: `cargo run -- <reads1> <reads2> --output-fasta out.fasta` ingests `.fastq[.gz]` or `.fasta[.gz]` pairs, reverse-complements `reads2`, supports `--analyse-alternatives` and `--reference` for validation.
- **Rust Tests**: `cargo test` exercises all modules including alternative path detection; add Rust tests for new features—they're the authoritative test suite.
	For all test directory structure, dataset generation, and test running instructions, refer to `tests/README.md`.
- **Python ↔ Rust Bridge**: `rust/src/python_bindings.rs` exposes `assemble_from_reads`, `analyse_reads`, and `analyse_alternative_paths` via PyO3; maintain API compatibility when updating Rust signatures.
- **Integration Workflow**: `scripts/integration/run_integration.sh` compares Rust and Python assemblies over `tests/fixtures/*.fastq`; ensure the Rust binary exists before running.
- **Results Storage**: Both pipelines write FASTA outputs under `tests/results/`; avoid hard-coding paths so integration tests discover artifacts automatically.
- **Logging & Metrics**: Rust uses `env_logger`; Python prints timing stats; both support `--metrics-csv` for benchmarking—keep CSV schema consistent across implementations.
- **External Dependencies**: Rust uses `strsim`, `sprs`, `bio`, `serde_json`; Python (deprecated core) uses Biopython + SciPy—account for these in build scripts.
- **Notebook Context**: `python/sequitur.ipynb` and `python/results.ipynb` document algorithmic intent; they may import deprecated Python code—update to use `sequitur` when possible.
- **Data Fixtures**: FASTQ samples live in `tests/fixtures/`; use them for deterministic regression tests instead of crafting new inline data.
- **Reverse Complements**: Both CLIs reverse-complement reads from the second FASTQ; maintain this behaviour in loaders and dataset handlers.
- **Error Handling**: Rust returns `anyhow::Result`; Python (deprecated) raises exceptions; PyO3 bindings convert Rust errors to Python exceptions—ensure error messages are actionable.
- **Performance Notes**: Rust is 30x faster than pure Python for alternative path analysis; always recommend Rust bindings over pure Python in new code.
- **Extending Features**: New algorithms go in Rust first (`rust/src/`), then expose via PyO3 if Python access needed; avoid duplicating logic in pure Python.
- **Release Checklist**: Run `cargo test`, `cargo build --release`, test PyO3 bindings with `maturin develop`, verify integration tests pass; pure Python tests are optional legacy validation.
