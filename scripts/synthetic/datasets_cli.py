#!/usr/bin/env python3
"""Unified synthetic dataset utility CLI for Sequitur.

Subcommands:
  encode      Convert text to DNA FASTA and paired FASTQ reads
  decode      Decode DNA FASTA back to text
  generate    Generate dataset by id from tests/datasets.yaml (synthetic only)
  validate    Validate assembly against reference using edlib
  diagnose    Diagnose overlaps in paired FASTQ reads
  list        List synthetic datasets from tests/datasets.yaml

Rationale:
  Consolidates previously separate scripts:
    - text_to_fastq.py (encode)
    - decode_dna_to_text.py (decode)
    - validate_assembly.py (validate)
    - diagnose_overlaps.py (diagnose)
    - generate_synthetic_datasets.sh (generate all)

South African English spelling for user-facing text (e.g. "optimise").
"""

import argparse
import sys
import random
import json
from pathlib import Path
from typing import List, Tuple, Dict, Any

try:
    import edlib  # only needed for validate
except ImportError:
    edlib = None

# ----------------------------- Encoding / Decoding -----------------------------

BASES = ['A', 'C', 'G', 'T']

def char_to_dna(char: str) -> str:
    code = ord(char)
    # Allow space (32) and line break (10) to be encoded
    if code == 10:  # line break '\n'
        idx = 127 - 32  # use next available codon after printable ASCII
        d0 = idx % 4
        d1 = (idx // 4) % 4
        d2 = (idx // 16) % 4
        d3 = (idx // 64) % 4
        return BASES[d3] + BASES[d2] + BASES[d1] + BASES[d0]
    elif code == 32:  # space
        idx = code - 32
        d0 = idx % 4
        d1 = (idx // 4) % 4
        d2 = (idx // 16) % 4
        d3 = (idx // 64) % 4
        return BASES[d3] + BASES[d2] + BASES[d1] + BASES[d0]
    elif code < 32 or code > 126:
        code = ord('?')
        idx = code - 32
        d0 = idx % 4
        d1 = (idx // 4) % 4
        d2 = (idx // 16) % 4
        d3 = (idx // 64) % 4
        return BASES[d3] + BASES[d2] + BASES[d1] + BASES[d0]
    else:
        idx = code - 32
        d0 = idx % 4
        d1 = (idx // 4) % 4
        d2 = (idx // 16) % 4
        d3 = (idx // 64) % 4
        return BASES[d3] + BASES[d2] + BASES[d1] + BASES[d0]

def text_to_sequence(text: str) -> str:
    return ''.join(char_to_dna(c) for c in text)

def build_reverse_map() -> Dict[str, str]:
    reverse = {}
    # Add space
    c = chr(32)
    dna = char_to_dna(c)
    reverse.setdefault(dna, c)
    # Add line break
    c = chr(10)
    dna = char_to_dna(c)
    reverse.setdefault(dna, c)
    for code in range(32, 127):
        c = chr(code)
        dna = char_to_dna(c)
        reverse.setdefault(dna, c)
    # Add line break codon (using next available codon)
    dna = char_to_dna(chr(10))
    reverse.setdefault(dna, chr(10))
    return reverse

def decode_sequence(dna_seq: str, reverse_map: Dict[str, str]) -> str:
    out = []
    i = 0
    while i < len(dna_seq):
        if i + 4 <= len(dna_seq):
            codon = dna_seq[i:i+4]
            out.append(reverse_map.get(codon, '?'))
            i += 4
        else:
            out.append('?')
            i += 1
    return ''.join(out)

# ----------------------------- FASTA / FASTQ IO -----------------------------

def write_fasta(sequence: str, path: Path, header: str) -> None:
    with open(path, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")

def write_fastq(reads: List[str], path: Path, tag: str) -> None:
    with open(path, 'w') as f:
        for i, read in enumerate(reads):
            f.write(f"@text_read_{i}/{tag}\n{read}\n+\n" + 'I' * len(read) + "\n")

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    seqs = []
    header = None
    buf = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    seqs.append((header, ''.join(buf)))
                header = line[1:]
                buf = []
            else:
                buf.append(line)
    if header is not None:
        seqs.append((header, ''.join(buf)))
    return seqs

def read_fastq(path: Path, max_reads: int = None) -> List[str]:
    out = []
    with open(path, 'r') as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()
            out.append(seq)
            if max_reads and len(out) >= max_reads:
                break
    return out

# ----------------------------- Generation logic -----------------------------

def reverse_complement(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def generate_paired(sequence: str, read_length: int, coverage: int, insert_size: int, ensure_boundaries: bool, seed: int) -> Tuple[List[str], List[str]]:
    random.seed(seed)
    seq_len = len(sequence)
    if seq_len < insert_size:
        raise ValueError(f"Sequence length {seq_len} shorter than insert size {insert_size}")
    total_bases = seq_len * coverage
    num_pairs = int(total_bases / (2 * read_length))
    forward = []
    reverse = []
    max_start = seq_len - insert_size - read_length
    for _ in range(num_pairs):
        if max_start < 0:
            start = random.randint(0, max(0, seq_len - read_length))
        else:
            start = random.randint(0, max_start)
        f_end = min(start + read_length, seq_len)
        fwd = sequence[start:f_end]
        r_start = min(start + insert_size, seq_len - read_length)
        r_end = min(r_start + read_length, seq_len)
        rev_seq = sequence[r_start:r_end]
        rev_rc = reverse_complement(rev_seq)
        if len(fwd) < read_length:
            fwd += 'N' * (read_length - len(fwd))
        if len(rev_rc) < read_length:
            rev_rc += 'N' * (read_length - len(rev_rc))
        forward.append(fwd)
        reverse.append(rev_rc)
    if ensure_boundaries:
        def make_pair(start: int) -> Tuple[str, str]:
            f_end = min(start + read_length, seq_len)
            fwd = sequence[start:f_end]
            r_start = min(start + insert_size, max(0, seq_len - read_length))
            r_end = min(r_start + read_length, seq_len)
            rev_seq = sequence[r_start:r_end]
            rev_rc = reverse_complement(rev_seq)
            if len(fwd) < read_length:
                fwd += 'N' * (read_length - len(fwd))
            if len(rev_rc) < read_length:
                rev_rc += 'N' * (read_length - len(rev_rc))
            return fwd, rev_rc
        f0, r0 = make_pair(0)
        end_start = max(0, seq_len - insert_size - read_length)
        fe, re = make_pair(end_start)
        forward.insert(0, f0)
        reverse.insert(0, r0)
        forward.append(fe)
        reverse.append(re)
    return forward, reverse

# ----------------------------- Overlap diagnosis -----------------------------

def find_overlap(seq1: str, seq2: str, min_overlap: int) -> int:
    max_o = min(len(seq1), len(seq2))
    for l in range(max_o, min_overlap - 1, -1):
        if seq1[-l:] == seq2[:l]:
            return l
    return 0

def diagnose(reads1: List[str], reads2: List[str], sample: int, no_revcomp: bool) -> Dict[str, Any]:
    if not no_revcomp:
        reads2 = [reverse_complement(r) for r in reads2]
    sample1 = reads1[:min(sample, len(reads1))]
    sample2 = reads2[:min(sample, len(reads2))]
    mate = []
    for i in range(min(len(sample1), len(sample2))):
        o = find_overlap(sample1[i], sample2[i], 10)
        if o > 0:
            mate.append(o)
    cross = []
    for i, r1 in enumerate(sample1):
        for j, r2 in enumerate(sample2):
            o = find_overlap(r1, r2, 20)
            if o >= 50:
                cross.append((i, j, o))
    return {
        'sample_forward': len(sample1),
        'sample_reverse': len(sample2),
        'mate_pairs_checked': min(len(sample1), len(sample2)),
        'mate_pairs_with_overlap': len(mate),
        'mate_avg_overlap': (sum(mate) / len(mate)) if mate else 0.0,
        'cross_significant': len(cross),
        'cross_top5': sorted(cross, key=lambda x: -x[2])[:5]
    }

# ----------------------------- Validation -----------------------------

def validate(ref: Path, asm: Path, detailed: bool) -> int:
    if edlib is None:
        print("Error: edlib not installed. Run: pip install edlib", file=sys.stderr)
        return 1
    ref_seq = ''.join(s for _, s in read_fasta(ref))
    asm_seq = ''.join(s for _, s in read_fasta(asm))
    result = edlib.align(ref_seq, asm_seq, mode="NW", task="path")
    dist = result['editDistance']
    cigar = result['cigar']
    max_len = max(len(ref_seq), len(asm_seq))
    accuracy = 1 - (dist / max_len) if max_len else 0
    len_ratio = (len(asm_seq) / len(ref_seq)) if ref_seq else 0
    print(f"Edit distance: {dist}\nAccuracy: {accuracy*100:.2f}%\nLength ratio (asm/ref): {len_ratio:.3f}")
    if detailed and cigar:
        print(f"CIGAR: {cigar}")
    if accuracy >= 0.95:
        print("✓ PASS: accuracy ≥ 95%")
        return 0
    elif accuracy >= 0.90:
        print("⚠ WARN: accuracy 90–95%")
        return 0
    else:
        print("✗ FAIL: accuracy < 90%")
        return 1

# ----------------------------- Dataset registry -----------------------------

def load_registry(path: Path) -> Dict[str, Any]:
    import yaml  # lazy import
    with open(path, 'r') as f:
        data = yaml.safe_load(f)
    synthetic = {}
    for entry in data.get('synthetic', []):
        synthetic[entry['id']] = entry
    return {'synthetic': synthetic, 'raw': data}

def list_datasets(registry: Dict[str, Any]) -> None:
    print("Synthetic datasets:")
    for did, meta in registry['synthetic'].items():
        keep = 'cached' if meta.get('keep_in_repo') else 'on-demand'
        desc = meta.get('description', '')
        print(f"  - {did:25s} [{keep}] {desc}")

def generate_dataset(dataset_id: str, registry: Dict[str, Any]) -> None:
    meta = registry['synthetic'].get(dataset_id)
    if not meta:
        print(f"Error: dataset id '{dataset_id}' not found", file=sys.stderr)
        sys.exit(1)
    gen = meta.get('generator')
    if not gen:
        print(f"Error: dataset '{dataset_id}' has no generator spec", file=sys.stderr)
        sys.exit(1)
    params = gen.get('params', {})
    source_path = Path(meta['paths']['source']) if 'source' in meta['paths'] else None
    if source_path and not source_path.exists():
        print(f"Error: source text file missing: {source_path}", file=sys.stderr)
        sys.exit(1)
    # Read source text (or synthesise placeholder)
    if source_path:
        text = source_path.read_text()
    else:
        text = f"SYNTHETIC DATASET {dataset_id}"
    dna = text_to_sequence(text)
    data_dir = Path(meta['paths']['data_dir'])
    data_dir.mkdir(parents=True, exist_ok=True)
    ref_path = data_dir / 'reference.fasta'
    write_fasta(dna, ref_path, header=f"{dataset_id}|source")
    read_len = int(params.get('read_length', 150))
    coverage = int(params.get('coverage', 30))
    insert_size = int(params.get('insert_size', 500))
    seed = int(params.get('seed', 42))
    ensure_boundaries = bool(params.get('ensure_boundaries', False))
    fwd, rev = generate_paired(dna, read_len, coverage, insert_size, ensure_boundaries, seed)
    write_fastq(fwd, data_dir / 'reads_1.fastq', '1')
    write_fastq(rev, data_dir / 'reads_2.fastq', '2')
    print(f"Generated dataset '{dataset_id}' in {data_dir}")

# ----------------------------- CLI parsing -----------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description='Unified synthetic dataset utilities')
    sub = p.add_subparsers(dest='cmd', required=True)
    # encode
    pe = sub.add_parser('encode', help='Convert text to DNA + paired reads')
    pe.add_argument('input', type=Path)
    pe.add_argument('output_dir', type=Path)
    pe.add_argument('--read-length', type=int, default=150)
    pe.add_argument('--coverage', type=int, default=30)
    pe.add_argument('--insert-size', type=int, default=500)
    pe.add_argument('--seed', type=int, default=42)
    pe.add_argument('--ensure-boundaries', action='store_true')
    # decode
    pd = sub.add_parser('decode', help='Decode DNA FASTA to text')
    pd.add_argument('input', type=Path)
    pd.add_argument('output', type=Path)
    # validate
    pv = sub.add_parser('validate', help='Validate assembly vs reference')
    pv.add_argument('reference', type=Path)
    pv.add_argument('assembly', type=Path)
    pv.add_argument('--detailed', action='store_true')
    # diagnose
    pg = sub.add_parser('diagnose', help='Diagnose overlaps in paired FASTQ reads')
    pg.add_argument('reads1', type=Path)
    pg.add_argument('reads2', type=Path)
    pg.add_argument('--sample', type=int, default=100)
    pg.add_argument('--no-revcomp', action='store_true')
    # list
    pl = sub.add_parser('list', help='List synthetic datasets from registry')
    pl.add_argument('--registry', type=Path, default=Path('tests/datasets.yaml'))
    # generate
    pgx = sub.add_parser('generate', help='Generate dataset by id from registry')
    pgx.add_argument('dataset_id')
    pgx.add_argument('--registry', type=Path, default=Path('tests/datasets.yaml'))
    return p

def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    cmd = args.cmd
    if cmd == 'encode':
        text = args.input.read_text()
        dna = text_to_sequence(text)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_fasta(dna, args.output_dir / 'reference.fasta', header=f"text_ref|{args.input.name}")
        fwd, rev = generate_paired(dna, args.read_length, args.coverage, args.insert_size, args.ensure_boundaries, args.seed)
        write_fastq(fwd, args.output_dir / 'reads_1.fastq', '1')
        write_fastq(rev, args.output_dir / 'reads_2.fastq', '2')
        print(f"Encode complete. Reads + reference in {args.output_dir}")
        return 0
    elif cmd == 'decode':
        seqs = read_fasta(args.input)
        reverse_map = build_reverse_map()
        with open(args.output, 'w') as out:
            for i, (h, s) in enumerate(seqs):
                txt = decode_sequence(s, reverse_map)
                if i:
                    out.write("\n" + '='*60 + f"\nSEQUENCE {i+1}: {h}\n" + '='*60 + "\n")
                out.write(txt)
        print(f"Decoded {len(seqs)} sequences to {args.output}")
        return 0
    elif cmd == 'validate':
        return validate(args.reference, args.assembly, args.detailed)
    elif cmd == 'diagnose':
        r1 = read_fastq(args.reads1, max_reads=args.sample * 2)
        r2 = read_fastq(args.reads2, max_reads=args.sample * 2)
        stats = diagnose(r1, r2, args.sample, args.no_revcomp)
        print(json.dumps(stats, indent=2))
        return 0
    elif cmd == 'list':
        reg = load_registry(args.registry)
        list_datasets(reg)
        return 0
    elif cmd == 'generate':
        reg = load_registry(args.registry)
        generate_dataset(args.dataset_id, reg)
        return 0
    else:
        parser.print_help()
        return 1

if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
