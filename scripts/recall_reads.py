#!/usr/bin/env python3
"""
Recall specific reads by matrix index from the original FASTQ files.

Usage:
    python recall_reads.py reads.json 0 10 25
    python recall_reads.py reads.json --all
    python recall_reads.py reads.json --range 10 20
"""

import json
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def load_read_map(json_path):
    """Load the read map JSON file."""
    with open(json_path) as f:
        return json.load(f)


def load_fastq_sequences(fastq_path):
    """Load all sequences from a FASTQ file into a list."""
    sequences = []
    for record in SeqIO.parse(fastq_path, "fastq"):
        sequences.append(str(record.seq).upper())
    return sequences


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def recall_reads(read_map_json, matrix_indices):
    """
    Recall specific reads by their matrix indices.
    
    Args:
        read_map_json: Path to the read map JSON file
        matrix_indices: List of matrix indices to recall
    
    Returns:
        List of dicts with matrix_index, source, source_index, and sequence
    """
    metadata = load_read_map(read_map_json)
    
    reads1_path = metadata["reads1_path"]
    reads2_path = metadata["reads2_path"]
    reads1_count = metadata["reads1_count"]
    reverse_complemented = metadata["reverse_complemented"]
    read_map = metadata["read_map"]
    
    # Load sequences from both files
    print(f"Loading sequences from {reads1_path}...", file=sys.stderr)
    reads1 = load_fastq_sequences(reads1_path)
    
    print(f"Loading sequences from {reads2_path}...", file=sys.stderr)
    reads2 = load_fastq_sequences(reads2_path)
    
    if reverse_complemented:
        print("Applying reverse complement to reads2...", file=sys.stderr)
        reads2 = [reverse_complement(seq) for seq in reads2]
    
    results = []
    for matrix_idx in matrix_indices:
        if matrix_idx >= len(read_map):
            print(f"Warning: Matrix index {matrix_idx} out of range (max: {len(read_map)-1})", file=sys.stderr)
            continue
        
        entry = read_map[matrix_idx]
        source = entry["source"]
        source_idx = entry["source_index"]
        
        if source == "reads1":
            sequence = reads1[source_idx]
        else:
            sequence = reads2[source_idx]
        
        results.append({
            "matrix_index": matrix_idx,
            "source": source,
            "source_index": source_idx,
            "sequence": sequence,
            "length": len(sequence)
        })
    
    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: python recall_reads.py reads.json <indices...>", file=sys.stderr)
        print("       python recall_reads.py reads.json --all", file=sys.stderr)
        print("       python recall_reads.py reads.json --range START END", file=sys.stderr)
        sys.exit(1)
    
    read_map_json = sys.argv[1]
    
    if not Path(read_map_json).exists():
        print(f"Error: {read_map_json} not found", file=sys.stderr)
        sys.exit(1)
    
    # Determine which indices to recall
    if len(sys.argv) == 3 and sys.argv[2] == "--all":
        metadata = load_read_map(read_map_json)
        matrix_indices = list(range(len(metadata["read_map"])))
    elif len(sys.argv) == 5 and sys.argv[2] == "--range":
        start = int(sys.argv[3])
        end = int(sys.argv[4])
        matrix_indices = list(range(start, end))
    else:
        matrix_indices = [int(idx) for idx in sys.argv[2:]]
    
    results = recall_reads(read_map_json, matrix_indices)
    
    # Print results in a readable format
    for result in results:
        print(f"Matrix Index: {result['matrix_index']}")
        print(f"Source: {result['source']} (index {result['source_index']})")
        print(f"Length: {result['length']}")
        print(f"Sequence: {result['sequence']}")
        print()


if __name__ == "__main__":
    main()
