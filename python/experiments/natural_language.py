"""Natural language toy experiments for the Sequitur assembler."""

from __future__ import annotations

import argparse
import random
import time
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from fastDamerauLevenshtein import damerauLevenshtein as damerau_levenshtein_distance

from sequitur_core import (
    adjacency_to_sparse,
    build_suffix_array,
    create_bipartite_adjacency_matrix,
    create_de_bruijn_graph,
    eulerian_path,
    find_longest_overlap,
    find_lower_diagonal_path,
)

ToyDataset = Tuple[str, List[str]]
TOY_DATASETS: Tuple[ToyDataset, ...] = (
    (
        "betty_bought_butter_the_butter_was_bitter_betty_bought_better_butter_to_make_the_bitter_butter_better",
        [
            "betty_bought_butter_th",
            "tter_the_butter_was_",
            "he_butter_was_bitter_",
            "as_bitter_betty_bought",
            "tty_bought_better_butter_t",
            "ught_better_butter_to",
            "r_butter_to_make_the_",
            "ke_the_bitter_butter_better",
        ],
    ),
    (
        "you say hello world, i bellow go to hell",
        [
            "you say hel",
            " say hello wo",
            "lo world, i be",
            "ld, i bellow go t",
            "ow go to hell",
        ],
    ),
    (
        "she_sells_sea_shells_on_the_sea_shore",
        [
            "she_sells_s",
            "lls_sea_shel",
            "ea_shells_o",
            "shells_on_the_s",
            "he_sea_s",
            "ea_shore",
        ],
    ),
)


def run_sequitur(
    reads: Sequence[str],
    *,
    max_diff: float,
    min_suffix_len: int,
) -> Tuple[str, float, float, float]:
    """Execute Sequitur on the provided reads and return timings."""

    suffix_start = time.time()
    suffix_array, suffix_lookup = build_suffix_array(reads, min_suffix_len=min_suffix_len)
    suffix_time = time.time() - suffix_start

    overlap_start = time.time()
    adjacency, overlaps = create_bipartite_adjacency_matrix(
        reads,
        suffix_array,
        suffix_lookup,
        max_diff=max_diff,
        min_suffix_len=min_suffix_len,
    )
    overlap_time = time.time() - overlap_start

    if not adjacency:
        return "", suffix_time, overlap_time, 0.0

    matrix = adjacency_to_sparse(adjacency, size=len(reads)).tocoo()
    rows = list(range(len(reads)))
    cols = list(range(len(reads)))
    assembled, assembly_time = find_lower_diagonal_path(
        matrix,
        overlaps,
        dict(enumerate(reads)),
        cols,
        rows,
        do_time=True,
    )
    return assembled, suffix_time, overlap_time, assembly_time


def run_debruijn(
    reads: Sequence[str],
    target: str,
) -> Tuple[str, float, float, int]:
    """Assemble using the reference De Bruijn implementation."""

    k_min, k_max = find_longest_overlap(reads)
    best_edit = float("inf")
    best_seq = ""
    dbg_time_total = 0.0
    euler_time_total = 0.0

    for k in range(k_min, k_max + 1):
        graph, dbg_time = create_de_bruijn_graph(k, reads, do_time=True)
        seq, euler_time = eulerian_path(graph, do_time=True)
        dbg_time_total += dbg_time
        euler_time_total += euler_time
        if not seq:
            continue
        edit = damerau_levenshtein_distance(seq, target, similarity=False)
        if edit < best_edit:
            best_edit = edit
            best_seq = seq
    return best_seq, dbg_time_total, euler_time_total, int(best_edit if best_edit != float("inf") else -1)


def experiment(
    dataset: ToyDataset,
    *,
    seeds: Iterable[int],
    repeats: int,
    max_diff: float,
    min_suffix_len: int,
) -> List[Tuple[int, int, float, float, float, int]]:
    """Run Sequitur on shuffled versions of the toy dataset."""

    target, reads = dataset
    results: List[Tuple[int, int, float, float, float, int]] = []

    for seed in seeds:
        for _ in range(repeats):
            shuffled = list(reads)
            random.Random(seed).shuffle(shuffled)
            assembled, t_suffix, t_overlap, t_assembly = run_sequitur(
                shuffled,
                max_diff=max_diff,
                min_suffix_len=min_suffix_len,
            )
            edit = damerau_levenshtein_distance(assembled, target, similarity=False)
            results.append((seed, len(shuffled), t_suffix, t_overlap, t_assembly, edit))
    return results


def write_results(
    path: Path,
    header: Sequence[str],
    rows: Iterable[Sequence[object]],
) -> None:
    if not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(",".join(header) + "\n")
        for row in rows:
            handle.write(",".join(str(field) for field in row) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sequitur toy natural-language experiment")
    parser.add_argument(
        "--dataset",
        type=int,
        default=-1,
        help="Dataset index (default all)",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="*",
        default=list(range(5)),
        help="Random seeds for shuffling",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=5,
        help="Repeats per seed",
    )
    parser.add_argument(
        "--max-diff",
        type=float,
        default=0.25,
        help="Max normalised edit distance",
    )
    parser.add_argument(
        "--min-suffix",
        type=int,
        default=3,
        help="Minimum suffix length",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/natural_language_sequences.sequitur.csv"),
        help="CSV output path",
    )
    parser.add_argument(
        "--compare-debruijn",
        action="store_true",
        help="Also run the De Bruijn baseline and write a comparison CSV",
    )
    parser.add_argument(
        "--dbg-output",
        type=Path,
        default=Path("data/output/natural_language_sequences.debruijn.csv"),
        help="CSV output path for De Bruijn baseline",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    indices = range(len(TOY_DATASETS)) if args.dataset < 0 else [args.dataset]
    header = (
        "dataset",
        "seed",
        "read_count",
        "suffix_array_time",
        "overlap_time",
        "assembly_time",
        "edit_distance",
    )
    all_rows: List[Tuple[int, int, int, float, float, float, int]] = []
    dbg_header = (
        "dataset",
        "seed",
        "read_count",
        "dbg_time",
        "euler_time",
        "edit_distance",
    )
    dbg_rows: List[Tuple[int, int, int, float, float, int]] = []
    for idx in indices:
        target, reads = TOY_DATASETS[idx]
        print(f"Running Sequitur on dataset {idx} ({len(reads)} reads)")
        rows = experiment(
            (target, reads),
            seeds=args.seeds,
            repeats=args.repeats,
            max_diff=args.max_diff,
            min_suffix_len=args.min_suffix,
        )
        for seed, read_count, t_suffix, t_overlap, t_assembly, edit in rows:
            all_rows.append((idx, seed, read_count, t_suffix, t_overlap, t_assembly, edit))
        if args.compare_debruijn:
            print(f"Running De Bruijn baseline on dataset {idx}")
            for seed in args.seeds:
                shuffled = list(reads)
                random.Random(seed).shuffle(shuffled)
                seq, t_dbg, t_euler, edit = run_debruijn(shuffled, target)
                dbg_rows.append((idx, seed, len(shuffled), t_dbg, t_euler, edit))
    write_results(args.output, header, all_rows)
    print(f"Wrote {len(all_rows)} Sequitur rows to {args.output}")
    if args.compare_debruijn:
        write_results(args.dbg_output, dbg_header, dbg_rows)
        print(f"Wrote {len(dbg_rows)} De Bruijn rows to {args.dbg_output}")


if __name__ == "__main__":
    main()
