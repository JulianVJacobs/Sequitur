"""Sparse matching utilities for Sequitur reconstruction."""

from __future__ import annotations

import time
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.sparse import coo_matrix, vstack, hstack, csr_matrix

OverlapMatrix = Mapping[int, Mapping[int, int]]


def adjacency_to_sparse(
    adjacency: Mapping[int, Mapping[int, int]],
    *,
    size: int | None = None,
) -> coo_matrix:
    """Convert nested adjacency dictionaries into a COO matrix (transposed)."""

    rows: List[int] = []
    cols: List[int] = []
    data: List[int] = []

    for row, col_dict in adjacency.items():
        for col, value in col_dict.items():
            rows.append(row)
            cols.append(col)
            data.append(value)

    if not rows:
        shape = (size or 0, size or 0)
        return coo_matrix(shape, dtype=int)

    max_dim = max(max(rows), max(cols)) + 1
    dim = size or max_dim
    matrix = coo_matrix((np.array(data), (np.array(rows), np.array(cols))), shape=(dim, dim))
    return matrix.T


def move_col(matrix: coo_matrix, cols_map: Mapping[int, int]) -> None:
    """Relabel COO column indices in-place."""

    for idx in range(len(matrix.col)):
        matrix.col[idx] = cols_map[matrix.col[idx]]


def move_row(matrix: coo_matrix, rows_map: Mapping[int, int]) -> None:
    """Relabel COO row indices in-place."""

    for idx in range(len(matrix.row)):
        matrix.row[idx] = rows_map[matrix.row[idx]]


def find_lower_diagonal_path(
    matrix: coo_matrix,
    overlaps: OverlapMatrix,
    reads_map: Mapping[int, str] | Iterable[Tuple[int, str]],
    cols: Sequence[int],
    rows: Sequence[int],
    *,
    quality_map: Mapping[int, Sequence[int]] | Iterable[Tuple[int, Sequence[int]]] | None = None,
    file: str | None = None,
    do_time: bool = False,
) -> Tuple[str, float] | str:
    """Reconstruct a sequence by traversing the lower diagonal of the matching."""

    if do_time:
        start = time.time()

    argpen = lambda arr: np.argpartition(arr, -2)[-2]

    new_cols = list(cols)
    matrix = matrix.tocoo()
    if matrix.sum(axis=0).min() == 0:
        new_cols = [c for c in new_cols if c not in [new_cols[matrix.sum(axis=0).argmin()]]] + [
            new_cols[matrix.sum(axis=0).argmin()]
        ]
    if matrix.sum(axis=1).min() == 0:
        if matrix.sum(axis=1).argmin() == matrix.sum(axis=0).argmin():
            new_cols = [new_cols[-1]] + [
                c
                for c in new_cols[:-1]
                if c not in [cols[matrix.getrow(rows.index(new_cols[-1])).argmax()]]
            ] + [cols[matrix.getrow(rows.index(new_cols[-1])).argmax()]]
        else:
            new_cols = [rows[matrix.sum(axis=1).argmin()]] + [
                c for c in new_cols if c not in [rows[matrix.sum(axis=1).argmin()]]
            ]

    cols_map = {cols.index(c): new_cols.index(c) for c in cols}
    move_col(matrix, cols_map)
    cols = new_cols

    new_rows = list(cols)
    rows_map = {rows.index(r): new_rows.index(r) for r in rows}
    move_row(matrix, rows_map)
    rows = new_rows

    i = len(rows)
    j = len(cols) - 1
    k = matrix.sum(axis=1).argmin() if matrix.sum(axis=1).min() == 0 else None

    while j > (k if matrix.sum(axis=1).min() == 0 else 0):
        if k is not None and matrix.getrow(rows.index(cols[j])).argmax() == k:
            cols_, c_ = [], 0
            while j + c_ + 1 < len(rows):
                c_ += 1
                if len(matrix.getrow(j + c_).nonzero()[1]) > 1:
                    cols_ = np.argpartition(matrix.getrow(j + c_).toarray().flatten(), -2)[::-1][:2]
                    if cols[cols_[1]] in cols[:j] and matrix.getcol(cols_[1]).argmax() == j + c_:
                        break
            if j + c_ + 1 == len(cols):
                new_cols = cols[: k + 1] + cols[j:] + cols[k + 1 : j]
            else:
                new_cols = (
                    cols[: k + 1]
                    + cols[j : j + c_]
                    + [c for c in cols[k + 1 : j] if c not in [cols[min(cols_)]]]
                    + [cols[min(cols_)]]
                    + cols[j + c_ :]
                )
            cols_map = {cols.index(c): new_cols.index(c) for c in cols}
            move_col(matrix, cols_map)
            cols = new_cols

            if j + c_ + 1 == len(rows):
                new_rows = list(cols)
            else:
                new_rows = (
                    cols[: k + c_ + 1]
                    + [r for r in rows[k : j + c_] if r not in cols[: k + c_ + 1] + cols[j + c_ :]]
                    + cols[j + c_ :]
                )
            rows_map = {rows.index(r): new_rows.index(r) for r in rows}
            move_row(matrix, rows_map)
            rows = new_rows

            i, j, k = j + c_ + 1, j + c_, k + c_
        else:
            # safe accessors to avoid IndexError / empty-data accesses
            def row_nonzero_count(idx: int) -> int:
                if idx < 0 or idx >= matrix.shape[0]:
                    return 0
                return len(matrix.getrow(idx).nonzero()[1])

            def row_argmax(idx: int) -> int:
                if idx < 0 or idx >= matrix.shape[0]:
                    return -1
                return matrix.getrow(idx).argmax()

            def row_getcol_data(idx: int, col_idx: int) -> int:
                if idx < 0 or idx >= matrix.shape[0] or col_idx < 0 or col_idx >= matrix.shape[1]:
                    return 0
                data = matrix.getrow(idx).getcol(col_idx).data
                return int(data[0]) if len(data) > 0 else 0

            cmax = row_argmax(rows.index(cols[j]))
            if row_nonzero_count(rows.index(cols[j])) > 1:
                cpen = argpen(matrix.getrow(rows.index(cols[j])).toarray().flatten())
                if cmax > j:
                    # check cmax+1 is in-bounds and has a meaningful second-best
                    if (
                        row_nonzero_count(cmax + 1) > 1
                        and row_getcol_data(cmax + 1, argpen(matrix.getrow(cmax + 1).toarray().flatten()))
                        >= row_getcol_data(rows.index(cols[j]), cpen)
                    ):
                        crange = [argpen(matrix.getrow(cmax).toarray().flatten()), cmax]
                    else:
                        crange = [cpen]
                else:
                    crange = [cmax]
            else:
                crange = [cmax]
            while crange[0] > j:
                if len(matrix.getrow(crange[0]).nonzero()[1]) > 1:
                    crange = [argpen(matrix.getrow(crange[0]).toarray().flatten())] + crange
                else:
                    crange = [matrix.getrow(crange[0]).argmax()] + crange
                if crange[0] == j:
                    crange = [matrix.getrow(crange[1]).argmax()] + crange[1:]

            new_cols = (
                [c for c in cols[:j] if c not in [cols[cr] for cr in crange]]
                + [cols[cr] for cr in crange]
                + [c for c in cols[j:] if c not in [cols[cr] for cr in crange]]
            )
            cols_map = {cols.index(c): new_cols.index(c) for c in cols}
            move_col(matrix, cols_map)
            cols = new_cols

            new_rows = [r for r in rows[:i] if r not in cols[j:]] + list(cols[j:])
            rows_map = {rows.index(r): new_rows.index(r) for r in rows}
            move_row(matrix, rows_map)
            rows = new_rows
        j -= 1
        i -= 1

    seq = ""
    reads_map = dict(reads_map)
    overlap_lookup: Dict[int, Dict[int, int]] = {
        src: dict(dst_map) for src, dst_map in overlaps.items()
    }

    def get_overlap(c_idx: int, r_idx: int) -> int:
        return overlap_lookup.get(c_idx, {}).get(r_idx, 0)

    def fallback_assemble() -> str:
        # Simple greedy fallback: start at read 0 (if present) and greedily follow highest overlap
        seq_parts: List[str] = []
        if not reads_map:
            return ""
        start = 0 if 0 in reads_map else next(iter(reads_map.keys()))
        seq_parts.append(reads_map[start])
        used = {start}
        current = start
        while True:
            nxt = overlaps.get(current, {})
            # choose unused target with maximum overlap length
            candidates = [(olen, tgt) for tgt, olen in nxt.items() if tgt not in used]
            if not candidates:
                break
            olen, tgt = max(candidates)
            if olen and olen < len(reads_map[tgt]):
                seq_parts.append(reads_map[tgt][olen:])
            else:
                seq_parts.append(reads_map[tgt])
            used.add(tgt)
            current = tgt
        # append any remaining reads
        for idx in sorted(set(reads_map.keys()) - used):
            seq_parts.append(reads_map[idx])
        return "".join(seq_parts)

    try:
        if quality_map is None:
            for r_idx, c_idx in zip(rows[1:], cols[:-1]):
                overlap_len = get_overlap(c_idx, r_idx)
                if overlap_len > 0:
                    seq += reads_map[c_idx][:-overlap_len]
                else:
                    seq += reads_map[c_idx]
            seq += reads_map.get(cols[-1], "")
        else:
            quality_lookup = {idx: list(q) for idx, q in dict(quality_map).items()}
            first = True
            for c_idx, r_idx in zip(cols, rows):
                overlap_len = get_overlap(c_idx, r_idx)
                if overlap_len <= 0:
                    # no overlap info: append full read
                    if first:
                        seq += reads_map[c_idx]
                        first = False
                    else:
                        seq += reads_map[c_idx]
                    continue
                if first:
                    seq += reads_map[c_idx][:-overlap_len]
                    first = False
                prefix = reads_map[c_idx][-overlap_len:]
                prefix_q = quality_lookup.get(c_idx, [])[-overlap_len:]
                suffix = reads_map[r_idx][:-overlap_len]
                suffix_q = quality_lookup.get(r_idx, [])[-overlap_len:]
                # pad quality arrays if needed
                while len(prefix_q) < overlap_len:
                    prefix_q.insert(0, 0)
                while len(suffix_q) < overlap_len:
                    suffix_q.insert(0, 0)
                for base_idx in range(overlap_len):
                    seq += prefix[base_idx] if prefix_q[base_idx] >= suffix_q[base_idx] else suffix[base_idx]
            last_overlap = get_overlap(cols[-1], rows[-1])
            if last_overlap > 0:
                seq += reads_map.get(rows[-1], "")[-last_overlap:]
            else:
                seq += reads_map.get(rows[-1], "")
    except Exception as e:
        # On any unexpected indexing errors, fall back to a simpler greedy assembler
        import sys

        print(f"Warning: traversal failed with {e}; using fallback assembler", file=sys.stderr)
        seq = fallback_assemble()

    if do_time:
        elapsed = time.time() - start
        if file is not None:
            SeqIO.write(
                [SeqRecord(Seq(seq), id="sequitur", description="")],
                file,
                "fasta",
            )
        return seq, elapsed
    return seq
