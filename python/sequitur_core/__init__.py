"""Core utilities for the Sequitur assembly prototype."""

from .overlap import (
    normalised_damerau_levenshtein_distance,
    build_suffix_array,
    create_bipartite_adjacency_matrix,
)
from .debruijn import (
    find_longest_overlap,
    create_de_bruijn_graph,
    eulerian_path,
)
from .matching import (
    adjacency_to_sparse,
    move_col,
    move_row,
    find_lower_diagonal_path,
)
from .alternative_paths import (
    detect_swap_squares,
    build_swap_graph,
    find_connected_components,
    is_cycle,
    analyse_alternatives,
)

__all__ = [
    "normalised_damerau_levenshtein_distance",
    "build_suffix_array",
    "create_bipartite_adjacency_matrix",
    "find_longest_overlap",
    "create_de_bruijn_graph",
    "eulerian_path",
    "adjacency_to_sparse",
    "move_col",
    "move_row",
    "find_lower_diagonal_path",
    "detect_swap_squares",
    "build_swap_graph",
    "find_connected_components",
    "is_cycle",
    "analyse_alternatives",
]
