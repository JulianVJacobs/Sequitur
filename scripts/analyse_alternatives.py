#!/usr/bin/env python3
"""
Analyse per-read alternatives exported from Sequitur assembly.

Usage:
    python analyse_alternatives.py alternatives.json
    python analyse_alternatives.py alternatives.json --read 42
    python analyse_alternatives.py alternatives.json --reads 10 20 30
    python analyse_alternatives.py alternatives.json --top 10
"""

import json
import sys
from pathlib import Path


def load_alternatives(json_path):
    """Load the per-read alternatives JSON file."""
    with open(json_path) as f:
        return json.load(f)


def analyse_read(alternatives, read_index):
    """Analyse alternatives for a specific read."""
    read_data = None
    for entry in alternatives:
        if entry["read_index"] == read_index:
            read_data = entry
            break
    
    if read_data is None:
        print(f"Read {read_index} has no successors (dead end or not in graph)")
        return
    
    successors = read_data["successors"]
    print(f"Read {read_index} has {len(successors)} possible successors:")
    print()
    
    for i, succ in enumerate(successors, 1):
        print(f"  {i}. Target: {succ['target']}")
        print(f"     Score: {succ['score']}")
        print(f"     Overlap Length: {succ['overlap_length']}")
        print()
    
    if len(successors) > 1:
        print("AMBIGUITY DETECTED:")
        best_score = successors[0]["score"]
        ties = [s for s in successors if s["score"] == best_score]
        if len(ties) > 1:
            print(f"  {len(ties)} successors tied with score {best_score}")
            for s in ties:
                print(f"    -> {s['target']} (overlap: {s['overlap_length']})")
        else:
            print(f"  Best: {successors[0]['target']} (score: {best_score})")
            print(f"  Next: {successors[1]['target']} (score: {successors[1]['score']})")
            print(f"  Score gap: {best_score - successors[1]['score']}")


def summarize_alternatives(alternatives):
    """Provide summary statistics about alternatives."""
    total_reads_with_successors = len(alternatives)
    total_successors = sum(len(entry["successors"]) for entry in alternatives)
    
    ambiguous_reads = [entry for entry in alternatives if len(entry["successors"]) > 1]
    highly_ambiguous = [entry for entry in alternatives if len(entry["successors"]) > 2]
    
    print("=== Assembly Graph Summary ===")
    print(f"Reads with successors: {total_reads_with_successors}")
    print(f"Total successor edges: {total_successors}")
    print(f"Average successors per read: {total_successors / total_reads_with_successors:.2f}")
    print()
    print(f"Ambiguous reads (>1 successor): {len(ambiguous_reads)}")
    print(f"Highly ambiguous (>2 successors): {len(highly_ambiguous)}")
    print()
    
    if ambiguous_reads:
        print("=== Top Ambiguous Reads ===")
        ambiguous_reads.sort(key=lambda x: len(x["successors"]), reverse=True)
        for entry in ambiguous_reads[:10]:
            read_idx = entry["read_index"]
            num_succ = len(entry["successors"])
            scores = [s["score"] for s in entry["successors"]]
            print(f"  Read {read_idx}: {num_succ} successors, scores: {scores[:5]}{'...' if len(scores) > 5 else ''}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyse_alternatives.py alternatives.json [--read INDEX | --reads INDEX... | --top N]", file=sys.stderr)
        sys.exit(1)
    
    alternatives_path = sys.argv[1]
    
    if not Path(alternatives_path).exists():
        print(f"Error: {alternatives_path} not found", file=sys.stderr)
        sys.exit(1)
    
    alternatives = load_alternatives(alternatives_path)
    
    if len(sys.argv) == 2:
        # Just summarize
        summarize_alternatives(alternatives)
    elif sys.argv[2] == "--read":
        # Analyse single read
        read_index = int(sys.argv[3])
        analyse_read(alternatives, read_index)
    elif sys.argv[2] == "--reads":
        # Analyse multiple reads
        read_indices = [int(idx) for idx in sys.argv[3:]]
        for read_index in read_indices:
            analyse_read(alternatives, read_index)
            print("=" * 60)
            print()
    elif sys.argv[2] == "--top":
        # Show top N most ambiguous reads
        n = int(sys.argv[3])
        ambiguous = [e for e in alternatives if len(e["successors"]) > 1]
        ambiguous.sort(key=lambda x: len(x["successors"]), reverse=True)
        
        print(f"=== Top {n} Most Ambiguous Reads ===")
        for entry in ambiguous[:n]:
            print()
            analyse_read(alternatives, entry["read_index"])
            print("=" * 60)


if __name__ == "__main__":
    main()
