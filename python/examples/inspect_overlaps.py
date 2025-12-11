#!/usr/bin/env python3
"""Quick helper to inspect overlap confidences and cycles using the Rust backend."""

from __future__ import annotations

import json
from pathlib import Path

import sequitur as sequitur


def main() -> None:
    reads = ["ACGT", "GTAA", "AACG", "CGTT"]
    analysis = sequitur.analyse_reads(reads)

    print("Adjacency (src -> dst, weight, overlap):")
    for src_idx, edges in enumerate(analysis["adjacency"]):
        if not edges:
            print(f"  read {src_idx}: no outgoing edges")
            continue
        for dst, weight, overlap in edges:
            print(f"  read {src_idx} -> {dst}: weight={weight}, overlap={overlap}")

    print("\nConfidence per outgoing edge:")
    for src_idx, edges in enumerate(analysis["confidences"]):
        if not edges:
            print(f"  read {src_idx}: no outgoing edges")
            continue
        for dst, probability in edges:
            print(f"  read {src_idx} -> {dst}: confidence={probability:.3f}")

    if analysis["cycles"]:
        print("\nDetected cycles (strongly connected components):")
        for comp in analysis["cycles"]:
            print("  nodes:", comp)
    else:
        print("\nNo cycles detected in the overlap graph.")

    # Persist the analysis in case it is helpful for notebook work.
    output_path = Path("analysis.json")
    output_path.write_text(json.dumps(analysis, indent=2))
    print(f"\nAnalysis written to {output_path.resolve()}")


if __name__ == "__main__":
    main()
