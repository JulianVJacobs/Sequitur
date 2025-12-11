"""
Overlap detection utilities for the Sequitur assembler (DEPRECATED).

WARNING: This pure Python implementation is deprecated and much slower than the Rust version.
For production, use the Rust library via PyO3 bindings (`sequitur`).

Threading logic here is legacy and not recommended for new code.
See docs/RUST_ARCHITECTURE.md for migration guide and usage.
"""

from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Iterable, Mapping, MutableMapping, Sequence, Tuple

from fastDamerauLevenshtein import (
    damerauLevenshtein as damerau_levenshtein_distance,
)

SuffixArray = Sequence[str]
SuffixLookup = Mapping[str, int]
Adjacency = Dict[int, Dict[int, int]]
OverlapLengths = Dict[int, Dict[int, int]]


def normalised_damerau_levenshtein_distance(
    suffix: str,
    prefix: str,
) -> Tuple[float, int]:
    """Return the normalised Damerau–Levenshtein distance and match count.

    The metric is normalised by the shorter string so it can compare a suffix of
    one read against a prefix of another. The companion integer is the number of
    matching positions (``m - d``) which we treat as the overlap weight.
    """

    overlap_len = min(len(suffix), len(prefix))
    if overlap_len == 0:
        return 0.0, 0
    distance = damerau_levenshtein_distance(
        suffix[-overlap_len:], prefix[:overlap_len], similarity=False
    )
    return distance / overlap_len, overlap_len - distance


def build_suffix_array(
    reads: Iterable[str],
    min_suffix_len: int = 3,
) -> Tuple[SuffixArray, SuffixLookup]:
    """Construct a suffix array enriched with read identifiers.

    Each entry is either ``suffix${read_idx}`` or ``prefix^{read_idx}``, mirroring
    the notebook implementation so downstream code can locate neighbouring
    prefixes/suffixes efficiently.
    """

    reads_list = list(reads)
    suffixes: list[str] = []

    for index, read in enumerate(reads_list):
        if len(read) < min_suffix_len:
            continue
        decorated = f"{read}${index}"
        base_len = len(read)
        for start in range(base_len - min_suffix_len + 1):
            suffixes.append(decorated[start:])

        prefix_view = decorated.replace("$", "^")
        pivot = prefix_view.index("^")
        for offset in range(base_len - min_suffix_len + 1):
            prefix = prefix_view[: pivot - offset] + prefix_view[pivot:]
            suffixes.append(prefix)

    suffixes.sort()
    lookup = {entry: idx for idx, entry in enumerate(suffixes)}
    return suffixes, lookup


def create_bipartite_adjacency_matrix(
    reads: Iterable[str],
    suffix_array: SuffixArray | None = None,
    suffix_lookup: SuffixLookup | None = None,
    *,
    max_diff: float = 0.25,
    min_suffix_len: int = 3,
    use_threads: bool = False,
    max_workers: int = 16,
    # DEPRECATED: Threading in Python is legacy. Use Rust for performance.
) -> Tuple[Adjacency, OverlapLengths]:
    """Build the weighted successor/predecessor mapping for Sequitur.

    Parameters
    ----------
    reads
        Iterable of input reads. Iterated exactly once and cached into a list.
    suffix_array / suffix_lookup
        Optional pre-computed suffix array and lookup produced by
        :func:`build_suffix_array`.
    max_diff
        Normalised Damerau–Levenshtein ratio threshold. Smaller is a better
        match; values above this limit are ignored.
    min_suffix_len
        Minimum overlap length to consider when generating suffix candidates.
    use_threads / max_workers
        Enable a thread pool version of the suffix scan. Useful for large read
        sets when Python's GIL is not a bottleneck (the edit distance routine is
        C-backed), but the sequential path is deterministic and easier to
        reason about.
    """

    reads_list = list(reads)
    if not reads_list:
        return {}, {}

    if suffix_array is None or suffix_lookup is None:
        suffix_array, suffix_lookup = build_suffix_array(
            reads_list, min_suffix_len=min_suffix_len
        )
    if not suffix_array:
        return {i: {} for i in range(len(reads_list))}, {
            i: {} for i in range(len(reads_list))
        }

    adjacency: Adjacency = {i: {} for i in range(len(reads_list))}
    overlap_lengths: OverlapLengths = {i: {} for i in range(len(reads_list))}

    def update_edge(source: int, target: int, score: int, overlap_len: int) -> None:
        current = adjacency[source].get(target)
        if current is None or score > current or (
            score == current and overlap_len > overlap_lengths[source][target]
        ):
            adjacency[source][target] = score
            overlap_lengths[source][target] = overlap_len

    def process_suffix(suffix_idx: int) -> None:
        suffix = reads_list[suffix_idx]
        if len(suffix) < min_suffix_len:
            return
        _up_cut = False
        _down_cut = False

        for span in range(min_suffix_len, len(suffix)):
            marker = f"{suffix[-span:]}${suffix_idx}"
            anchor = suffix_lookup.get(marker)
            if anchor is None:
                continue

            forward = anchor + 1
            while forward < len(suffix_array) and (
                "$" in suffix_array[forward]
                or suffix_array[forward].endswith(str(suffix_idx))
            ):
                forward += 1

            while forward < len(suffix_array):
                candidate = suffix_array[forward]
                if "^" not in candidate:
                    break
                prefix, prefix_idx_str = candidate.split("^")
                prefix_idx = int(prefix_idx_str)
                if prefix_idx == suffix_idx:
                    forward += 1
                    continue
                diff, score = normalised_damerau_levenshtein_distance(
                    suffix[-span:], prefix
                )
                if diff > max_diff:
                    _down_cut = True
                    break
                overlap = min(len(suffix[-span:]), len(prefix))
                update_edge(suffix_idx, prefix_idx, score, overlap)

                forward += 1
                while forward < len(suffix_array) and (
                    "$" in suffix_array[forward]
                    or suffix_array[forward].endswith(str(suffix_idx))
                ):
                    forward += 1

            backward = anchor - 1
            while backward >= 0 and (
                "$" in suffix_array[backward]
                or suffix_array[backward].endswith(str(suffix_idx))
            ):
                backward -= 1

            while backward >= 0:
                candidate = suffix_array[backward]
                if "^" not in candidate:
                    break
                prefix, prefix_idx_str = candidate.split("^")
                prefix_idx = int(prefix_idx_str)
                if prefix_idx == suffix_idx:
                    backward -= 1
                    continue
                diff, score = normalised_damerau_levenshtein_distance(
                    suffix[-span:], prefix
                )
                if diff > max_diff:
                    _up_cut = True
                    break
                overlap = min(len(suffix[-span:]), len(prefix))
                update_edge(suffix_idx, prefix_idx, score, overlap)

                backward -= 1
                while backward >= 0 and (
                    "$" in suffix_array[backward]
                    or suffix_array[backward].endswith(str(suffix_idx))
                ):
                    backward -= 1

            if _up_cut and _down_cut:
                break

    import warnings
    if use_threads:
        warnings.warn(
            "Python threading is deprecated. Use Rust implementation for performance.",
            DeprecationWarning,
        )
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            list(executor.map(process_suffix, range(len(reads_list))))
    else:
        for idx in range(len(reads_list)):
            process_suffix(idx)

    return adjacency, overlap_lengths
