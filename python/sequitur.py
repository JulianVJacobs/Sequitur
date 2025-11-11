"""Command line entrypoint for the Sequitur assembly prototype."""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fastDamerauLevenshtein import damerauLevenshtein as damerau_levenshtein_distance

from sequitur_core import (
    adjacency_to_sparse,
    build_suffix_array,
    create_bipartite_adjacency_matrix,
    find_lower_diagonal_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Sequitur assembler")
    parser.add_argument("reads1", type=Path, help="Path to the first FASTQ file")
    parser.add_argument("reads2", type=Path, help="Path to the second FASTQ file")
    parser.add_argument(
        "--reference",
        type=Path,
        help="Reference FASTA used to score the reconstruction",
    )
    parser.add_argument(
        "--output-fasta",
        type=Path,
        help="Optional output FASTA path for the assembled sequence",
    )
    parser.add_argument(
        "--metrics-csv",
        type=Path,
        help="Optional CSV file to append assembly metrics",
    )
    parser.add_argument(
        "--max-diff",
        type=float,
        default=0.25,
        help="Maximum normalised Damerau-Levenshtein difference",
    )
    parser.add_argument(
        "--min-suffix",
        type=int,
        default=3,
        help="Minimum suffix length to consider for overlaps",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Worker threads for overlap construction (1 disables threading)",
    )
    return parser.parse_args()


def load_reads(read1_path: Path, read2_path: Path) -> Tuple[List[str], List[List[int]]]:
    sequences: List[str] = []
    qualities: List[List[int]] = []

    for record in SeqIO.parse(read1_path, "fastq"):
        sequences.append(str(record.seq.upper()))
        qualities.append(list(record.letter_annotations["phred_quality"]))

    for record in SeqIO.parse(read2_path, "fastq"):
        seq = record.seq.upper().reverse_complement()
        sequences.append(str(seq))
        qualities.append(list(record.letter_annotations["phred_quality"][::-1]))

    return sequences, qualities


def load_reference(reference_path: Path | None) -> str | None:
    if reference_path is None:
        return None
    record = next(SeqIO.parse(reference_path, "fasta"), None)
    return str(record.seq) if record is not None else None


def ensure_header(csv_path: Path, header: Sequence[str]) -> None:
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open("w", encoding="utf-8") as handle:
            handle.write(",".join(header) + "\n")


def append_metrics(csv_path: Path, row: Sequence[str], header: Sequence[str]) -> None:
    ensure_header(csv_path, header)
    with csv_path.open("a", encoding="utf-8") as handle:
        handle.write(",".join(row) + "\n")


def assemble(args: argparse.Namespace) -> None:
    reads, qualities = load_reads(args.reads1, args.reads2)
    if not reads:
        raise RuntimeError("No reads were parsed from the supplied FASTQ files")

    reference = load_reference(args.reference)

    suffix_start = time.time()
    suffix_array, suffix_lookup = build_suffix_array(
        reads, min_suffix_len=args.min_suffix
    )
    suffix_time = time.time() - suffix_start

    overlap_start = time.time()
    adjacency, overlaps = create_bipartite_adjacency_matrix(
        reads,
        suffix_array,
        suffix_lookup,
        max_diff=args.max_diff,
        min_suffix_len=args.min_suffix,
        use_threads=args.threads > 1,
        max_workers=max(args.threads, 1),
    )
    overlap_time = time.time() - overlap_start

    matrix = adjacency_to_sparse(adjacency, size=len(reads)).tocoo()
    rows = list(range(len(reads)))
    cols = list(range(len(reads)))

    reads_map: Dict[int, str] = dict(enumerate(reads))
    quality_map: Dict[int, List[int]] = dict(enumerate(qualities))

    if len(reads) < 2 or matrix.nnz == 0:
        assembled = reads[0]
        assembly_time = 0.0
    else:
        assembled, assembly_time = find_lower_diagonal_path(
            matrix,
            overlaps,
            reads_map,
            cols,
            rows,
            quality_map=quality_map,
            do_time=True,
            file=str(args.output_fasta) if args.output_fasta else None,
        )

    if args.output_fasta and (len(reads) < 2 or matrix.nnz == 0):
        SeqIO.write(
            [SeqRecord(Seq(assembled), id="sequitur", description="")],
            args.output_fasta,
            "fasta",
        )

    edit_distance = None
    if reference is not None:
        edit_distance = damerau_levenshtein_distance(assembled, reference, similarity=False)

    print("Sequitur assembly complete")
    print(f"Reads processed       : {len(reads)}")
    print(f"Suffix array time     : {suffix_time:.3f}s")
    print(f"Overlap matrix time   : {overlap_time:.3f}s")
    print(f"Assembly time         : {assembly_time:.3f}s")
    if reference is not None and edit_distance is not None:
        print(f"Edit distance to ref  : {edit_distance}")

    if args.metrics_csv is not None:
        header = (
            "edit_distance",
            "target_length",
            "assembly_length",
            "suffix_array_time",
            "overlap_time",
            "assembly_time",
        )
        row = (
            str(edit_distance) if edit_distance is not None else "",
            str(len(reference)) if reference is not None else "",
            str(len(assembled)),
            f"{suffix_time:.6f}",
            f"{overlap_time:.6f}",
            f"{assembly_time:.6f}",
        )
        append_metrics(args.metrics_csv, row, header)


def main() -> None:
    args = parse_args()
    assemble(args)


if __name__ == "__main__":
    main()
