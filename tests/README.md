# Sequitur Test Datasets

This directory contains all test datasets for the Sequitur assembler.

## Structure

```
tests/
├── datasets.yaml           # Master registry of all datasets
├── synthetic/              # Synthetic test datasets
│   ├── {dataset_id}/
│   │   ├── data/           # Reference and reads
│   │   └── results/        # Assembly outputs (gitignored)
├── real/                   # Real sequencing datasets
│   └── {dataset_id}/
│       ├── data/           # Downloaded data (gitignored for large files)
│       └── results/        # Assembly outputs (gitignored)
└── fixtures/               # DEPRECATED - old test location (to be removed)
```

## Dataset Organization

Each dataset follows a consistent structure:
- **data/** contains reference genome and paired-end reads
- **results/** contains assembly outputs (assembly.fasta, alternatives.json, etc.)
- Results directories are gitignored to avoid committing large generated files

## Synthetic Datasets

### Small Tests (in repository)
- **simple** - Minimal 16bp test case
- **swap_test** - Swap-square detection test
- **ambiguous** - Ambiguous path resolution
- **repeats** - Repeat handling
- **reverse_complement** - RC handling
- **with_errors** - Sequencing error tolerance

### Large Tests (generated on-demand)
- **bee_movie** - Natural language text encoded as DNA (86KB text → 344KB DNA)
  - Regenerate with: See datasets.yaml for command
  - Provides realistic assembly test with controlled ground truth

## Real Datasets

Downloaded from public repositories (SRA, NCBI):
- **SRR390728** - RNA-Seq DLBCL cell line
- **SRR11117158** - Raphanus sativus transcriptome

## Usage

### Run Assembly on a Dataset

```bash
# Synthetic dataset
./rust/target/release/sequitur_rs \
  tests/synthetic/{dataset}/data/reads_1.fastq \
  tests/synthetic/{dataset}/data/reads_2.fastq \
  --output-fasta tests/synthetic/{dataset}/results/assembly.fasta \
  --reference tests/synthetic/{dataset}/data/reference.fasta \
  --analyse-alternatives \
  --alternatives-json tests/synthetic/{dataset}/results/alternatives.json

# Validate synthetic assembly
python3 scripts/synthetic/validate_assembly.py \
  tests/synthetic/{dataset}/data/reference.fasta \
  tests/synthetic/{dataset}/results/assembly.fasta
```

### Regenerate Large Synthetic Datasets

```bash
# See datasets.yaml for exact parameters
python3 scripts/synthetic/text_to_fastq.py \
  tests/synthetic/bee_movie/source.txt \
  tests/synthetic/bee_movie/data \
  --read-length 250 --coverage 50 --insert-size 350 --seed 42 --ensure-boundaries
```

### Download Real Datasets

```bash
# Use provided download script (if available)
bash scripts/download_dataset.sh {dataset_id}

# Or manually download using URLs in datasets.yaml
```

## Migration from Old Structure

The old `tests/fixtures/` structure has been reorganized:
- Small synthetic tests → `tests/synthetic/{dataset}/data/`
- Real data → `tests/real/{dataset}/data/`
- All results → `tests/{synthetic|real}/{dataset}/results/`

Old fixture files can be removed after verification.
