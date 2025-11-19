#!/usr/bin/env python3
"""
Convert plain text files to synthetic FASTQ read pairs for assembly testing.

Text is mapped to DNA alphabet to make it compatible with genomic tools while
remaining human-readable for debugging. Fragments text into overlapping reads
with configurable coverage.

Usage:
    python text_to_fastq.py input.txt output_prefix [--read-length 150] [--coverage 30]

Output:
    output_prefix_1.fastq  (forward reads)
    output_prefix_2.fastq  (reverse reads)
    output_prefix_ref.fasta (reference text as single sequence)
"""

import argparse
import sys
import random
from pathlib import Path


def char_to_dna(char):
    """Map printable ASCII to DNA bases with bijective (1-to-1) encoding."""
    # Use fixed 4-base codons (4^4 = 256) to represent printable ASCII (95 chars)
    # ASCII 32-126 = 95 printable chars

    code = ord(char)
    if code < 32 or code > 126:
        code = ord('?')

    idx = code - 32  # 0-94
    bases = ['A', 'C', 'G', 'T']

    # Convert idx into 4 base-4 digits: d3 d2 d1 d0 (most-significant first)
    d0 = idx % 4
    d1 = (idx // 4) % 4
    d2 = (idx // 16) % 4
    d3 = (idx // 64) % 4

    return bases[d3] + bases[d2] + bases[d1] + bases[d0]


def text_to_sequence(text):
    """Convert text to DNA sequence."""
    return ''.join(char_to_dna(c) for c in text)


def sequence_to_text(seq):
    """Best-effort reverse: DNA sequence back to text (for debugging)."""
    reverse_map = {}
    for code in range(32, 127):  # printable ASCII
        c = chr(code)
        dna = char_to_dna(c)
        reverse_map[dna] = c

    result = []
    i = 0
    # Decode fixed 4-base codons
    while i < len(seq):
        if i + 4 <= len(seq):
            codon = seq[i:i+4]
            result.append(reverse_map.get(codon, '?'))
            i += 4
        else:
            # leftover bases not forming a full codon
            result.append('?')
            i += 1
    
    return ''.join(result)


def generate_paired_reads(sequence, read_length, coverage, insert_size=500):
    """
    Fragment sequence into overlapping paired-end reads.
    
    Args:
        sequence: DNA sequence (text converted to DNA)
        read_length: Length of each read
        coverage: Target coverage depth
        insert_size: Distance between forward and reverse read starts
    
    Returns:
        (forward_reads, reverse_reads) as lists of strings
    """
    seq_len = len(sequence)
    if seq_len < insert_size:
        raise ValueError(f"Sequence too short ({seq_len}) for insert size {insert_size}")
    
    # Calculate number of read pairs needed
    total_bases_needed = seq_len * coverage
    num_pairs = int(total_bases_needed / (2 * read_length))
    
    forward_reads = []
    reverse_reads = []
    
    # Generate reads with random starts (uniform coverage)
    for i in range(num_pairs):
        # Random start position ensuring both reads fit
        max_start = seq_len - insert_size - read_length
        if max_start < 0:
            # Sequence too short, use overlapping positions
            start = random.randint(0, max(0, seq_len - read_length))
        else:
            start = random.randint(0, max_start)
        
        # Forward read
        fwd_end = min(start + read_length, seq_len)
        fwd = sequence[start:fwd_end]
        
        # Reverse read (from other end of insert, reverse complemented)
        rev_start = min(start + insert_size, seq_len - read_length)
        rev_end = min(rev_start + read_length, seq_len)
        rev_seq = sequence[rev_start:rev_end]
        rev_rc = reverse_complement(rev_seq)
        
        # Pad if needed (shouldn't happen with proper sizing)
        if len(fwd) < read_length:
            fwd += 'N' * (read_length - len(fwd))
        if len(rev_rc) < read_length:
            rev_rc += 'N' * (read_length - len(rev_rc))
        
        forward_reads.append(fwd)
        reverse_reads.append(rev_rc)
    
    return forward_reads, reverse_reads


def reverse_complement(seq):
    """Reverse complement DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def write_fastq(reads, output_path, read_type='1'):
    """Write reads to FASTQ file with constant high quality scores."""
    with open(output_path, 'w') as f:
        for i, read in enumerate(reads):
            qual = 'I' * len(read)  # Phred+33 quality score 40 (99.99% accuracy)
            f.write(f"@text_read_{i}/{read_type}\n")
            f.write(f"{read}\n")
            f.write("+\n")
            f.write(f"{qual}\n")


def write_fasta(sequence, output_path, header="text_reference"):
    """Write sequence to FASTA file."""
    with open(output_path, 'w') as f:
        f.write(f">{header}\n")
        # Wrap at 80 characters
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Convert text to synthetic FASTQ reads for assembly testing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('input', type=Path, help='Input text file')
    parser.add_argument('output_prefix', type=str, help='Output file prefix')
    parser.add_argument('--read-length', type=int, default=150,
                       help='Read length in bases (default: 150)')
    parser.add_argument('--coverage', type=int, default=30,
                       help='Target coverage depth (default: 30)')
    parser.add_argument('--insert-size', type=int, default=500,
                       help='Insert size for paired reads (default: 500)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--ensure-boundaries', action='store_true',
                       help='Ensure at least one read pair starts at position 0 and one at the end')
    
    args = parser.parse_args()
    
    # Set random seed
    random.seed(args.seed)
    
    # Read input text
    print(f"Reading text from {args.input}")
    with open(args.input, 'r') as f:
        text = f.read()
    
    print(f"Text length: {len(text)} characters")
    
    # Convert to DNA
    print("Converting text to DNA sequence...")
    sequence = text_to_sequence(text)
    print(f"DNA sequence length: {len(sequence)} bases ({len(sequence)//4} codons)")
    
    # Generate reads
    print(f"Generating reads (length={args.read_length}, coverage={args.coverage}x)...")
    forward_reads, reverse_reads = generate_paired_reads(
        sequence, args.read_length, args.coverage, args.insert_size
    )
    # Optionally ensure boundary coverage: add deterministic pairs at start and near the end
    if args.ensure_boundaries:
        seq_len = len(sequence)
        # helper to create a pair at a given start
        def make_pair(start):
            # forward
            fwd_end = min(start + args.read_length, seq_len)
            fwd = sequence[start:fwd_end]
            # reverse mate
            rev_start = min(start + args.insert_size, max(0, seq_len - args.read_length))
            rev_end = min(rev_start + args.read_length, seq_len)
            rev_seq = sequence[rev_start:rev_end]
            rev_rc = reverse_complement(rev_seq)
            if len(fwd) < args.read_length:
                fwd += 'N' * (args.read_length - len(fwd))
            if len(rev_rc) < args.read_length:
                rev_rc += 'N' * (args.read_length - len(rev_rc))
            return fwd, rev_rc

        # pair at the very start
        f0, r0 = make_pair(0)
        # pair near the end: choose start so the forward read and mate fit near the tail
        end_start = max(0, seq_len - args.insert_size - args.read_length)
        fend, rend = make_pair(end_start)
        # prepend start pair and append end pair for deterministic coverage
        forward_reads.insert(0, f0)
        reverse_reads.insert(0, r0)
        forward_reads.append(fend)
        reverse_reads.append(rend)
    print(f"Generated {len(forward_reads)} read pairs")
    
    # Write outputs
    fwd_path = f"{args.output_prefix}_1.fastq"
    rev_path = f"{args.output_prefix}_2.fastq"
    ref_path = f"{args.output_prefix}_ref.fasta"
    
    print(f"Writing forward reads to {fwd_path}")
    write_fastq(forward_reads, fwd_path, read_type='1')
    
    print(f"Writing reverse reads to {rev_path}")
    write_fastq(reverse_reads, rev_path, read_type='2')
    
    print(f"Writing reference sequence to {ref_path}")
    write_fasta(sequence, ref_path, header=f"text_ref|{args.input.name}")
    
    print("\nDone! Test assembly with:")
    print(f"  sequitur_rs {fwd_path} {rev_path} --output-fasta assembled.fasta")
    print(f"\nValidate with:")
    print(f"  python scripts/synthetic/validate_assembly.py {ref_path} assembled.fasta")


if __name__ == '__main__':
    main()
