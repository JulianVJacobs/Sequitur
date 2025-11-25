Sequitur Python Reference
=========================

This directory holds the original Python prototype of the Sequitur assembler along with
supporting notebooks and benchmarking scripts. It stays feature-complete with the
work described in the dissertation and serves as the ground truth while the Rust port
is developed.

Directory layout
----------------

- `sequitur_core/` – reusable library code (suffix array, overlap graph, matching).
- `sequitur.py` – command line entry point that assembles paired FASTQ reads.
- `experiments/` – reproducible drivers for the natural-language toy datasets and
  other exploratory runs migrated from the notebook.
- `scripts/` – legacy scripts kept for reference (De Bruijn baseline, earlier
  prototypes, cluster helpers).
- `sequitur.ipynb`, `results.ipynb` – notebooks with exploratory analysis and
  figures used throughout the research.
- `results.csv` – cached outputs from earlier experiments.

Python environment
------------------

The code targets Python 3.10+. Create a virtual environment and install the required
packages before running the CLI or notebooks:

```bash
cd python
python -m venv .venv
source .venv/bin/activate
pip install numpy scipy biopython fastDamerauLevenshtein networkx pylcs
deactivate  # optional when you are done
```

The experiments also rely on `hungarian-algorithm`, `memory-profiler`, and `scipy`'s
sparse matrix support. Install them as needed, or capture the full environment with
`pip freeze > requirements.txt` once you are satisfied with the version set.

Running the assembler
---------------------

From the `python/` directory:

```bash
source .venv/bin/activate
PYTHONPATH=. python sequitur.py \
  data/input/reads.1.fastq \
  data/input/reads.2.fastq \
  --reference data/input/reference.fasta \
  --output-fasta data/output/assembly.fasta \
  --metrics-csv data/output/metrics.csv
```

The script expects paired FASTQ reads and will write the assembled contig, along with
optional metrics comparing the reconstruction to a supplied reference.

Running experiments
-------------------

The natural-language toy experiments can be reproduced with:

```bash
source .venv/bin/activate
PYTHONPATH=. python experiments/natural_language.py --compare-debruijn
```

This generates CSV reports under `data/output/` (mirroring the notebook setup).

Working with notebooks
----------------------

## Notebook Import and Dependency Instructions

When working with Jupyter notebooks in this directory:

- **Place all imports in the first code cell** of the notebook for clarity and reproducibility.
- **Install any non-standard dependencies** (not in the Python standard library) using an inline pip command in the first cell, e.g.:
  ```python
  %pip install matplotlib networkx plotly
  ```
- Do **not repeat imports** in subsequent cells; reference only the first cell for all imports.

Refer to this section for best practices when creating or editing notebooks in the Sequitur Python environment.

Launch Jupyter with the virtual environment activated:

```bash
source .venv/bin/activate
jupyter notebook
```

Open `sequitur.ipynb` or `results.ipynb` to explore the detailed methodology,
plots, and intermediate validations.

Relationship to the Rust port
-----------------------------

The Rust implementation (in `../rust/`) will eventually replace the performance-critical
parts of this Python code. Until feature parity is achieved, keep this reference version
up-to-date and use it for validation against the new port.
