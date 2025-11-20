#!/usr/bin/env python3
"""
Diagnose overlap graph issues by sampling reads and checking what overlaps are detected.

Usage:
    python diagnose_overlaps.py reads1.fastq reads2.fastq [--sample 100]
"""

import argparse
import sys
from pathlib import Path


def read_fastq(path, max_reads=None):
    """Read FASTQ, return list of sequences."""
    sequences = []
    with open(path, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()
            sequences.append(seq.upper())
            if max_reads and len(sequences) >= max_reads:
                break
    return sequences


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def find_overlap(seq1, seq2, min_overlap=10):
    """Find overlap length between suffix of seq1 and prefix of seq2."""
    max_overlap = min(len(seq1), len(seq2))
    
    for overlap_len in range(max_overlap, min_overlap - 1, -1):
        if seq1[-overlap_len:] == seq2[:overlap_len]:
            return overlap_len
    return 0


def analyze_read_set(reads, label, sample_size=100):
    """Analyze a set of reads for self-overlaps and duplicates."""
    print(f"\n{'='*60}")
    print(f"Analyzing {label}")
    print(f"{'='*60}")
    print(f"Total reads: {len(reads)}")
    
    # Check for exact duplicates
    unique_reads = set(reads)
    print(f"Unique reads: {len(unique_reads)}")
    if len(unique_reads) < len(reads):
        print(f"⚠ WARNING: {len(reads) - len(unique_reads)} duplicate reads detected!")
    
    # Sample reads for overlap analysis
    sample = reads[:min(sample_size, len(reads))]
    
    # Find overlaps within sample
    overlaps = []
    for i, read1 in enumerate(sample):
        for j, read2 in enumerate(sample):
            if i != j:
                overlap = find_overlap(read1, read2, min_overlap=20)
                if overlap >= 50:  # Report significant overlaps
                    overlaps.append((i, j, overlap))
    
    print(f"\nSample of {len(sample)} reads:")
    print(f"  Significant overlaps (≥50bp): {len(overlaps)}")
    
    if overlaps:
        print(f"  Top 5 overlaps:")
        for i, j, overlap in sorted(overlaps, key=lambda x: -x[2])[:5]:
            print(f"    Read {i} → Read {j}: {overlap} bp")
    
    # Check if reads overlap with themselves (should be rare)
    self_overlaps = [(i, j, o) for i, j, o in overlaps if i == j]
    if self_overlaps:
        print(f"  ⚠ Self-overlaps detected: {len(self_overlaps)}")


def analyze_cross_overlaps(reads1, reads2, sample_size=100):
    """Check overlaps between read1 and read2 sets."""
    print(f"\n{'='*60}")
    print(f"Cross-set overlap analysis (Read1 → Read2)")
    print(f"{'='*60}")
    
    sample1 = reads1[:min(sample_size, len(reads1))]
    sample2 = reads2[:min(sample_size, len(reads2))]
    
    overlaps = []
    for i, read1 in enumerate(sample1):
        for j, read2 in enumerate(sample2):
            overlap = find_overlap(read1, read2, min_overlap=20)
            if overlap >= 50:
                overlaps.append((i, j, overlap))
    
    print(f"Sample: {len(sample1)} reads1 × {len(sample2)} reads2")
    print(f"Significant overlaps (≥50bp): {len(overlaps)}")
    
    if overlaps:
        print(f"Top 5 cross-set overlaps:")
        for i, j, overlap in sorted(overlaps, key=lambda x: -x[2])[:5]:
            print(f"  Read1[{i}] → Read2[{j}]: {overlap} bp")


def main():
    parser = argparse.ArgumentParser(
        description='Diagnose overlap graph structure',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('reads1', type=Path, help='Forward reads FASTQ')
    parser.add_argument('reads2', type=Path, help='Reverse reads FASTQ')
    parser.add_argument('--sample', type=int, default=100,
                       help='Number of reads to sample for analysis')
    parser.add_argument('--no-revcomp', action='store_true',
                       help='Skip reverse complement of reads2')
    
    args = parser.parse_args()
    
    if not args.reads1.exists():
        print(f"Error: {args.reads1} not found", file=sys.stderr)
        return 1
    
    if not args.reads2.exists():
        print(f"Error: {args.reads2} not found", file=sys.stderr)
        return 1
    
    print(f"Reading {args.reads1}...")
    reads1 = read_fastq(args.reads1, max_reads=args.sample * 2)
    
    print(f"Reading {args.reads2}...")
    reads2_raw = read_fastq(args.reads2, max_reads=args.sample * 2)
    
    if args.no_revcomp:
        reads2 = reads2_raw
        print("Skipping reverse complement (--no-revcomp)")
    else:
        print("Applying reverse complement to reads2...")
        reads2 = [reverse_complement(seq) for seq in reads2_raw]
    
    # Analyze each set
    analyze_read_set(reads1, "Forward reads (reads1)", args.sample)
    analyze_read_set(reads2, "Reverse reads (reads2, rev-comp'd)", args.sample)
    
    # Check mate-pair overlaps (reads1[i] should overlap reads2[i])
    print(f"\n{'='*60}")
    print(f"Mate-pair overlap check")
    print(f"{'='*60}")
    mate_overlaps = []
    for i in range(min(args.sample, len(reads1), len(reads2))):
        overlap = find_overlap(reads1[i], reads2[i], min_overlap=10)
        if overlap > 0:
            mate_overlaps.append((i, overlap))
    
    print(f"Checked {min(args.sample, len(reads1), len(reads2))} mate pairs")
    print(f"Pairs with overlap: {len(mate_overlaps)}")
    if mate_overlaps:
        avg_overlap = sum(o for _, o in mate_overlaps) / len(mate_overlaps)
        print(f"Average mate overlap: {avg_overlap:.1f} bp")
        print(f"Top 5 mate overlaps:")
        for i, overlap in sorted(mate_overlaps, key=lambda x: -x[1])[:5]:
            print(f"  Pair {i}: {overlap} bp")
    
    # Cross-set overlaps
    analyze_cross_overlaps(reads1, reads2, args.sample)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
