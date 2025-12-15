#!/usr/bin/env python3
"""Generate a deterministic, "perfect" read set from the original bee_movie reference.

Outputs to `tests/synthetic/bee_movie_high_overlap_simple/data/`:
- copies the original reference.fasta
- writes `reads_1.fastq` and `reads_2.fastq` with deterministic tiled pairs
- writes `manifest.json` with parameters and coverage summary

Defaults: `read_length=100`, `insert_size=200`, `step=50` (deterministic overlap).
"""

from pathlib import Path
import importlib.util
import json


def load_datasets_cli():
    path = Path(__file__).resolve().parent / 'datasets_cli.py'
    spec = importlib.util.spec_from_file_location('datasets_cli', str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    mod = load_datasets_cli()
    # Source reference (original bee_movie_high_overlap)
    src_ref = Path('tests') / 'synthetic' / 'bee_movie_high_overlap' / 'data' / 'reference.fasta'
    if not src_ref.exists():
        raise SystemExit(f"Source reference not found: {src_ref}")

    outdir = Path('tests') / 'synthetic' / 'bee_movie_high_overlap_simple' / 'data'
    outdir.mkdir(parents=True, exist_ok=True)

    # Read source reference
    seqs = mod.read_fasta(src_ref)
    ref_seq = ''.join(s for _, s in seqs)

    # Copy reference into dataset dir (preserve header)
    mod.write_fasta(ref_seq, outdir / 'reference.fasta', header='bee_movie_high_overlap|original')

    # Parameters (deterministic)
    read_length = 100
    insert_size = 200
    step = 50

    fwd_reads = []
    rev_reads = []
    ref_len = len(ref_seq)

    # Tile reads deterministically
    for start in range(0, max(1, ref_len), step):
        # forward
        f = ref_seq[start:start + read_length]
        if len(f) < read_length:
            f = f + 'N' * (read_length - len(f))
        # mate start
        r_start = start + insert_size
        if r_start < 0 or r_start >= ref_len:
            r = 'N' * read_length
        else:
            r = ref_seq[r_start:r_start + read_length]
            if len(r) < read_length:
                r = r + 'N' * (read_length - len(r))
        # reverse-complement mate (standard paired-end orientation)
        r_rc = mod.reverse_complement(r)
        fwd_reads.append(f)
        rev_reads.append(r_rc)

    # Write FASTQ files
    mod.write_fastq(fwd_reads, outdir / 'reads_1.fastq', '1')
    mod.write_fastq(rev_reads, outdir / 'reads_2.fastq', '2')

    # Compute naive coverage from record positions
    covered = [0] * ref_len
    for i, start in enumerate(range(0, max(1, ref_len), step)):
        # mark forward
        for p in range(start, min(ref_len, start + read_length)):
            covered[p] += 1
        # mark mate
        r_start = start + insert_size
        for p in range(r_start, min(ref_len, r_start + read_length)):
            if 0 <= p < ref_len:
                covered[p] += 1

    covered_positions = sum(1 for c in covered if c > 0)
    coverage_percent = covered_positions / ref_len * 100 if ref_len else 0.0

    manifest = {
        'source_reference': str(src_ref),
        'out_reference': str(outdir / 'reference.fasta'),
        'read_length': read_length,
        'insert_size': insert_size,
        'step': step,
        'ref_length': ref_len,
        'reads_per_file': len(fwd_reads),
        'covered_positions': covered_positions,
        'coverage_percent': coverage_percent
    }

    with open(outdir / 'manifest.json', 'w') as mf:
        json.dump(manifest, mf, indent=2)

    print(f"Wrote perfect reads to {outdir}")
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
