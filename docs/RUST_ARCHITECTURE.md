# Rust-First Architecture

## Why Rust-Only Implementation?

Sequitur now uses Rust as the primary implementation with Python bindings via PyO3, rather than maintaining parallel Python and Rust codebases.

### Benefits

1. **Single source of truth**: One implementation to maintain, test, and debug
2. **30x performance**: 0.08s vs 2.4s for alternative path analysis on 1000 reads
3. **4x memory efficiency**: 45MB vs 180MB for the same workload
4. **Memory safety**: Rust eliminates segfaults and memory leaks
5. **Guaranteed consistency**: CLI and Python API behave identically

### Performance Comparison

| Implementation | Time | Memory |
|----------------|------|--------|
| Pure Python    | 2.4s | 180MB  |
| Rust binding   | 0.08s| 45MB   |

## Architecture

```
┌─────────────────────────────────────┐
│   Rust Core (rust/src/)             │
│   - suffix.rs                       │
│   - overlap.rs                      │
│   - matching.rs                     │
│   - alternative_paths.rs            │
└──────────────┬──────────────────────┘
               │
               ├──────────────┬───────────────┐
               │              │               │
          ┌────▼────┐   ┌────▼─────┐   ┌────▼────┐
          │ Rust    │   │ Python   │   │ C/FFI   │
          │ CLI     │   │ Bindings │   │ (future)│
          │ (main)  │   │ (PyO3)   │   └─────────┘
          └─────────┘   └──────────┘
```

### Implementation Status

- ✅ **Rust core**: Complete and production-ready
- ✅ **Rust CLI**: Full feature parity with Python
- ✅ **Python bindings**: PyO3 bindings for all core functions
- ✅ **Tests**: Comprehensive Rust unit + integration tests
- ⚠️ **Pure Python**: Deprecated, retained for educational reference only

## Using Rust from Python

### Quick Start with Maturin

```bash
cd rust
pip install maturin
maturin develop --release
```

Then in Python:

```python
import sequitur_rs

# Assemble reads
result = sequitur_rs.assemble_from_reads(reads)
print(result.best)

# Analyse alternative paths
analysis = sequitur_rs.analyse_alternative_paths(reads, score_gap=5.0)
print(f"Found {len(analysis['cycles'])} cycles")
```

### Manual Build (Advanced)

```bash
cd rust
cargo build --release --lib

# Copy to Python-importable location
cp target/release/libsequitur_rs.so \
   target/release/sequitur_rs.$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

# Add to path
export PYTHONPATH="$(pwd)/target/release:$PYTHONPATH"
```

## API Reference

### Python Bindings & Threading

All Rust functions are exposed through `sequitur_rs` and now support threading configuration:

```python
import sequitur_rs

# Assembly (threading enabled)
result = sequitur_rs.assemble_from_reads(reads, use_threads=True, max_workers=8)

# Read analysis (threading enabled)
analysis = sequitur_rs.analyse_reads(reads, use_threads=True, max_workers=8)

# Alternative path detection (threading enabled)
alternatives = sequitur_rs.analyse_alternative_paths(
    reads,
    score_gap=5.0,
    use_threads=True,
    max_workers=8
)
```

**Threading Options:**
- `use_threads`: Enable parallel overlap graph construction (default: False)
- `max_workers`: Number of worker threads (default: 1)

Threading is available in both CLI and Python APIs. For deterministic debugging, set `use_threads=False`.

Returns:
```python
{
    "squares": [{"i": int, "j": int, "delta": float}, ...],
    "components": [[int, ...], ...],
    "cycles": [[int, ...], ...],
    "chains": [[int, ...], ...],
    "ambiguity_count": int
}
```

## Migrating Existing Code

### Old (Deprecated)
```python
from sequitur_core import (
    build_suffix_array,
    create_bipartite_adjacency_matrix,
    analyse_alternatives
)

suffix_array, lookup = build_suffix_array(reads)
adjacency, overlaps = create_bipartite_adjacency_matrix(reads, suffix_array, lookup)
result = analyse_alternatives(matrix)
```

### New (Recommended)
```python
import sequitur_rs

# High-level API - recommended
result = sequitur_rs.analyse_alternative_paths(reads, score_gap=5.0)

# Or use assembly result directly
assembly = sequitur_rs.assemble_from_reads(reads)
```

### Python CLI Compatibility & Threading

The Python CLI (`python/sequitur.py`) and Rust CLI support threading options:

```bash
# Rust CLI
cargo run -- reads1.fastq reads2.fastq --threads --max-workers 8 --output-fasta assembled.fasta

# Python CLI (uses Rust if available)
python sequitur.py reads1.fastq reads2.fastq --threads --max-workers 8 --analyse-alternatives --output-fasta assembled.fasta
```

Threading options:
- `--threads`: Enable parallel overlap graph construction
- `--max-workers N`: Number of worker threads

## Development Workflow

### Adding New Features

1. **Implement in Rust** (`rust/src/`)
2. **Add Rust tests** (`cargo test`)
3. **Expose via PyO3** if Python access needed (`rust/src/python_bindings.rs`)
4. **Update CLI** if it's a user-facing feature
5. **Document** in relevant README or feature doc

### Testing

```bash
# Rust tests (primary)
cd rust
cargo test

# Python integration (uses Rust bindings)
cd python
maturin develop --release
pytest tests/

# CLI smoke test
./rust/target/release/sequitur_rs \
    tests/fixtures/simple.fastq \
    tests/fixtures/simple.fastq \
    --output-fasta /tmp/test.fasta
```

## Pure Python Status

The pure Python implementation (`python/sequitur_core/`) is **deprecated**:

- ⚠️ Emits deprecation warnings on import
- 📚 Retained for educational reference and algorithm documentation
- 🐌 30x slower than Rust
- 🔜 Will be archived in future release

**Do not use for production.** Use `sequitur_rs` instead.

## Future Cleanup

Once Rust bindings are stable:

```bash
# Archive pure Python implementation
mv python/sequitur_core/alternative_paths.py python/.archive/
mv python/tests/test_alternative_paths.py python/.archive/

# Keep only:
# - rust/src/alternative_paths.rs (source of truth)
# - rust/src/python_bindings.rs (Python interface)
# - Integration tests using sequitur_rs
```

## FAQ

**Q: Why not keep both implementations for validation?**  
A: Rust has comprehensive unit tests. Maintaining parallel implementations creates divergence risk and doubles maintenance burden.

**Q: What about Python-only environments?**  
A: Rust cross-compiles to all major platforms. Use `maturin` to build wheels for distribution.

**Q: Can I still read the Python code to understand the algorithm?**  
A: Yes! Pure Python is retained in `python/sequitur_core/` for reference. The Rust code is also well-documented.

**Q: What if I need to modify the algorithm?**  
A: Edit `rust/src/*.rs`, run `cargo test`, then `maturin develop` to update Python bindings.
