# Documentation Index

Feature-specific and architecture deep-dives for Sequitur.  
For quick-start and project overview, see the [root README](../README.md).

## Feature Documentation

* **[ALTERNATIVE_PATHS.md](ALTERNATIVE_PATHS.md)** — Alternative assembly path detection: swap squares, cycles, and chain analysis
* **Threading & Parallelization** — See [RUST_ARCHITECTURE.md](RUST_ARCHITECTURE.md) for threading config, CLI and Python API usage, and migration guide

## Architecture

- **[RUST_ARCHITECTURE.md](RUST_ARCHITECTURE.md)** — Why Rust-only implementation, API reference, migration guide

## Implementation Guides

- **[python/README.md](../python/README.md)** — Python setup, CLI, notebooks
- **[rust/README.md](../rust/README.md)** — Rust setup, CLI, Cargo tasks
- **[rust/README_PyO3.md](../rust/PyO3.md)** — Building Python bindings with Maturin

# Or view directly in VS Code / GitHub
```

## Maintenance

When adding new features:
1. Update relevant implementation README (`python/` or `rust/`)
2. Add feature documentation to `docs/` if it's a major feature
3. Update this index to reference the new documentation
4. Update root README if it affects quick start or roadmap
