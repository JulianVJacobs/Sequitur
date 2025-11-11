"""Reference De Bruijn assembler used for Sequitur comparisons."""

from __future__ import annotations

import time
from typing import Iterable, Sequence, Tuple, Union

import networkx as nx
import pylcs


def find_longest_overlap(reads: Sequence[str]) -> Tuple[int, int]:
    """Return the minimum and maximum pairwise overlap across the reads."""

    if not reads:
        return 0, 0
    overlaps = []
    unique_reads = list(reads)
    for read in unique_reads:
        others = list(set(unique_reads).difference({read}))
        if not others:
            continue
        overlaps.extend(pylcs.lcs2_of_list(read, others))
    if not overlaps:
        return 0, 0
    return min(overlaps), max(overlaps)


def create_de_bruijn_graph(
    k: int,
    sequences: Iterable[str],
    *,
    do_time: bool = False,
) -> Union[nx.DiGraph, Tuple[nx.DiGraph, float]]:
    """Build a simple k-mer De Bruijn graph from the supplied sequences."""

    if do_time:
        start = time.time()
    graph = nx.DiGraph()
    for sequence in sequences:
        if len(sequence) < k:
            continue
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i : i + k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            if graph.has_edge(prefix, suffix):
                graph[prefix][suffix]["weight"] += 1
            else:
                graph.add_edge(prefix, suffix, weight=1)
    if do_time:
        return graph, time.time() - start
    return graph


def eulerian_path(
    graph: nx.DiGraph,
    *,
    do_time: bool = False,
) -> Union[str, Tuple[str, float]]:
    """Reconstruct a contig using an Eulerian path if one exists."""

    if do_time:
        start = time.time()
    if nx.has_eulerian_path(graph):
        path_nodes = list(nx.eulerian_path(graph))
        if not path_nodes:
            result = ""
        else:
            contig = [path_nodes[0][0]]
            contig.extend(edge[1][-1] for edge in path_nodes)
            result = "".join(contig)
        if do_time:
            return result, time.time() - start
        return result

    # Fallback: greedily follow available edges
    graph_copy = graph.copy()
    result_nodes = []
    start_node = next(
        (node for node in graph_copy.nodes if graph_copy.out_degree(node) > 0),
        None,
    )
    if start_node is None:
        result = ""
    else:
        current = start_node
        while True:
            result_nodes.append(current)
            successors = list(graph_copy.successors(current))
            if not successors:
                break
            nxt = successors[0]
            graph_copy.remove_edge(current, nxt)
            current = nxt
        result = "".join(result_nodes)

    if do_time:
        return result, time.time() - start
    return result
