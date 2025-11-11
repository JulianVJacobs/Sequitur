# Sequitur

Sequitur is an experimental genome assembler that blends De Bruijn-graph efficiency
with overlap-layout-consensus accuracy. It now lives in this repository as
two parallel tracks:

- a feature-complete Python reference implementation (used for validation and rapid
	iteration);
- a memory-conscious Rust port aimed at large-scale datasets.

![Sequitur overview diagram](./images/sequitur_process_correct.png "Sequitur process overview")

## Repository layout

- `python/` – reference Python implementation, notebooks, and experiment drivers.
- `rust/` – in-progress Rust port with Cargo project scaffolding.
- `images/` – figures and poster assets referenced in documentation.
- `julian jacobs 1605267 masters dissertation.pdf` – the original thesis describing
	the Sequitur methodology in detail.

Archived Python scripts are preserved under `python/.archive/` for historical context.

## Quick start

### Python reference pipeline

```bash
cd python
python -m venv .venv
source .venv/bin/activate
pip install numpy scipy biopython fastDamerauLevenshtein networkx pylcs
PYTHONPATH=. python sequitur.py \
	data/input/example.1.fastq \
	data/input/example.2.fastq \
	--reference data/input/example.fasta \
	--output-fasta data/output/example.sequitur.fasta \
	--metrics-csv data/output/example.metrics.csv
```

See `python/README.md` for environment notes, experiment runners, and notebook usage.

### Rust prototype

```bash
cd rust
cargo build --release
cargo run -- --help
```

Implementation of the Rust modules is ongoing; the CLI currently acts as a placeholder.

## Roadmap

- Port suffix-array construction, overlap scoring, and matching to Rust (tracked in
	`rust/` and the issue/todo list).
- Reconcile results between Python and Rust on toy datasets and realistic read sets.
- Package reproducible benchmarks and automate CI checks once the Rust path is feature
	complete.

## Background

If you are new to the project, start with the dissertation PDF and the notebooks in
`python/` for a detailed walkthrough of the algorithm, experiments, and lessons learned.
The natural-language toy datasets provide an intuitive introduction before tackling
real FASTQ data.

Contributions are welcome—open an issue or start a discussion if you plan to add new
features or fix bugs.
