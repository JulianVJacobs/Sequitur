#!/usr/bin/env python3
"""Example: Detect alternative assembly paths in a dataset.

This script demonstrates how to use the alternative path detection API
directly without going through the CLI.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from Bio import SeqIO
from sequitur_core import (
    adjacency_to_sparse,
    analyse_alternatives,
    build_suffix_array,
    create_bipartite_adjacency_matrix,
)


def main():
    # Example reads with potential ambiguity
    reads = [
        "AAAABBBBCCCC",
        "BBBBCCCCDDDD",
        "CCCCDDDDEEEE",
        "AAAADDDDEEEE",
    ]
    
    print(f"Analysing {len(reads)} reads for alternative paths...")
    print()
    
    # Build overlap graph
    suffix_array, suffix_lookup = build_suffix_array(reads, min_suffix_len=3)
    adjacency, overlaps = create_bipartite_adjacency_matrix(
        reads,
        suffix_array,
        suffix_lookup,
        max_diff=0.25,
        min_suffix_len=3,
    )
    
    # Convert to sparse matrix
    matrix = adjacency_to_sparse(adjacency, size=len(reads)).tocoo()
    
    print(f"Overlap matrix: {matrix.shape[0]}×{matrix.shape[1]} with {matrix.nnz} non-zeros")
    print()
    
    # Analyse alternatives
    result = analyse_alternatives(matrix, score_gap=None)
    
    print("Alternative Path Analysis:")
    print(f"  Swap squares detected: {len(result['squares'])}")
    print(f"  Ambiguous components:  {len(result['components'])}")
    print(f"  Cycles detected:       {len(result['cycles'])}")
    print(f"  Linear chains:         {len(result['chains'])}")
    print(f"  Total ambiguous pos:   {result['ambiguity_count']}")
    print()
    
    if result['squares']:
        print("Swap squares (i, j, delta):")
        for i, j, delta in result['squares'][:10]:  # Show first 10
            print(f"  ({i}, {j}): delta = {delta:.2f}")
    
    if result['cycles']:
        print("\nCycles detected:")
        for idx, cycle in enumerate(result['cycles'], 1):
            print(f"  Cycle {idx}: positions {cycle}")
    
    if result['chains']:
        print("\nLinear chains:")
        for idx, chain in enumerate(result['chains'], 1):
            print(f"  Chain {idx}: positions {chain}")


if __name__ == "__main__":
    main()
