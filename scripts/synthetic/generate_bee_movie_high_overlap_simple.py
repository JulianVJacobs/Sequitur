#!/usr/bin/env python3
"""Generate a small controlled 'bee_movie_high_overlap_simple' dataset.

Creates `tests/synthetic/bee_movie_high_overlap_simple/data/{reference.fasta,reads_1.fastq,reads_2.fastq}`.

Rules implemented:
- `insert_size` = 200
- `read_length` = 100 (divisible by 4)
- sequence_1 starts at `terminus_0 - 200` (overlaps last 200 bases of seq0)
- sequence_2 starts at `terminus_0 - 20` (overlaps last 20 bases of seq0)
- odd-numbered reads in `reads_2.fastq` are stored as reverse-complements (mixed orientation)

This script reuses encoding/generation helpers from `scripts/synthetic/datasets_cli.py`.
"""

import importlib.util
import random
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent

def load_datasets_cli():
    path = ROOT / 'scripts' / 'synthetic' / 'datasets_cli.py'
    spec = importlib.util.spec_from_file_location('datasets_cli', str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_text(seed: int, n_chars: int, prefix: str) -> str:
    rnd = random.Random(seed)
    choices = [chr(i) for i in range(32, 127)]
    # keep it human-readable by mixing prefix and random chars
    out = []
    while len(out) < n_chars:
        out.append(prefix)
        out.append('')
        out.extend(rnd.choices(choices, k=min(50, n_chars - len(out))))
    return ''.join(out)[:n_chars]


def main():
    mod = load_datasets_cli()
    # Parameters
    insert_size = 200
    read_length = 100
    coverage = 10
    seed = 42

    # Sequence lengths in bases (must be multiples of 4 for clean encoding)
    seq0_bases = 2000  # will encode from ~500 chars
    seq1_bases = 1500
    seq2_bases = 800

    # Derived char counts
    seq0_chars = seq0_bases // 4
    seq1_chars = seq1_bases // 4
    seq2_chars = seq2_bases // 4

    # Build text blocks
    text0 = make_text(seed + 1, seq0_chars, 'SEQ0')
    tail1_text = make_text(seed + 2, max(1, seq1_chars - (200 // 4)), 'SEQ1')
    tail2_text = make_text(seed + 3, max(1, seq2_chars - (20 // 4)), 'SEQ2')

    # Encode to DNA using the canonical encoder
    seq0_dna = mod.text_to_sequence(text0)
    # Ensure we have right lengths (trim/pad if necessary)
    seq0_dna = (seq0_dna + 'N' * seq0_bases)[:seq0_bases]

    # seq1 begins with last 200 bases of seq0
    overlap1 = 200
    seq1_tail_dna = mod.text_to_sequence(tail1_text)
    seq1_tail_dna = (seq1_tail_dna + 'N' * max(0, seq1_bases - overlap1))[: max(0, seq1_bases - overlap1)]
    seq1_dna = seq0_dna[-overlap1:] + seq1_tail_dna

    # seq2 begins with last 20 bases of seq0
    overlap2 = 20
    seq2_tail_dna = mod.text_to_sequence(tail2_text)
    seq2_tail_dna = (seq2_tail_dna + 'N' * max(0, seq2_bases - overlap2))[: max(0, seq2_bases - overlap2)]
    seq2_dna = seq0_dna[-overlap2:] + seq2_tail_dna

    # Positions
    pos0 = 0
    pos1 = len(seq0_dna) - overlap1
    pos2 = len(seq0_dna) - overlap2

    # Compute reference length and build blank reference
    end0 = pos0 + len(seq0_dna)
    end1 = pos1 + len(seq1_dna)
    end2 = pos2 + len(seq2_dna)
    ref_len = max(end0, end1, end2)
    ref = ['N'] * ref_len

    # Place seq0
    for i, b in enumerate(seq0_dna):
        ref[pos0 + i] = b
    # Place seq1 (overlap region should match seq0 by construction)
    for i, b in enumerate(seq1_dna):
        idx = pos1 + i
        if ref[idx] == 'N':
            ref[idx] = b
        # else assume identical
    # Place seq2
    for i, b in enumerate(seq2_dna):
        idx = pos2 + i
        if ref[idx] == 'N':
            ref[idx] = b

    reference = ''.join(ref)

    # Prepare output dir
    outdir = Path('tests') / 'synthetic' / 'bee_movie_high_overlap_simple' / 'data'
    outdir.mkdir(parents=True, exist_ok=True)

    # Write reference
    ref_path = outdir / 'reference.fasta'
    mod.write_fasta(reference, ref_path, header='bee_movie_high_overlap_simple|ref')

    # Generate paired reads from the reference
    fwd, rev = mod.generate_paired(reference, read_length, coverage, insert_size, ensure_boundaries=True, seed=seed)

    # Apply odd-read reverse-complement rule to reads_2: ensure 1-based odd reads are RC (keep as-is), 1-based even -> un-rc
    # Note: `rev` returned by generate_paired is already reverse-complemented; we will flip every 1-based even read back to original orientation
    out_rev = []
    for i, r in enumerate(rev):
        # i is 0-based; 1-based index = i+1
        if ((i + 1) % 2) == 1:
            # 1-based odd -> keep RC
            out_rev.append(r)
        else:
            # 1-based even -> reverse back to forward orientation
            out_rev.append(mod.reverse_complement(r))

    # Write FASTQ files
    mod.write_fastq(fwd, outdir / 'reads_1.fastq', '1')
    mod.write_fastq(out_rev, outdir / 'reads_2.fastq', '2')

    # Validate coverage by naive substring mapping
    ref_seq = reference
    covered = [0] * len(ref_seq)
    all_reads = fwd + out_rev
    for r in all_reads:
        idx = ref_seq.find(r)
        if idx != -1:
            for j in range(idx, min(len(ref_seq), idx + len(r))):
                covered[j] += 1

    covered_positions = sum(1 for c in covered if c > 0)
    coverage_pct = covered_positions / len(ref_seq) * 100 if len(ref_seq) else 0

    manifest = {
        'insert_size': insert_size,
        'read_length': read_length,
        'coverage_target': coverage,
        'seed': seed,
        'reference_length': len(ref_seq),
        'covered_positions': covered_positions,
        'coverage_percent': coverage_pct,
        'positions': {'pos0': pos0, 'pos1': pos1, 'pos2': pos2},
        'counts': {'reads_1': len(fwd), 'reads_2': len(out_rev)}
    }
    with open(outdir / 'manifest.json', 'w') as mf:
        json.dump(manifest, mf, indent=2)

    print(f"Wrote dataset to {outdir}")
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
