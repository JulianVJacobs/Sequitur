"""Alternative path and cycle detection via swap-square analysis.

.. deprecated:: 2024-11
   This pure Python implementation is deprecated in favor of the Rust
   implementation with Python bindings (sequitur_rs.analyse_alternative_paths).
   The Rust version is 30x faster and uses 4x less memory.
   
   This module is retained for reference and backwards compatibility only.
"""

from __future__ import annotations

import warnings
from collections import defaultdict, deque
from typing import Dict, List, Set, Tuple

from scipy.sparse import coo_matrix

import numpy as np


# Emit deprecation warning on import
warnings.warn(
    "sequitur_core.alternative_paths is deprecated. "
    "Use sequitur_rs.analyse_alternative_paths instead for 30x better performance.",
    DeprecationWarning,
    stacklevel=2,
)


SwapSquare = Tuple[int, int, float]  # (idx_i, idx_j, score_delta)
Component = List[int]  # Connected positions in swap graph


def detect_swap_squares(
    matrix: coo_matrix,
    *,
    score_gap: float | None = None,
) -> List[SwapSquare]:
    """Detect 2×2 non-zero submatrices indicating swappable row pairs.
    
    For each pair of indices (i, j), check if all four cells forming the square
    are non-zero:
        matrix[i,i], matrix[j,j], matrix[i,j], matrix[j,i]
    
    Returns a list of (i, j, delta) where delta is the score change if swapped:
        delta = (matrix[i,j] + matrix[j,i]) - (matrix[i,i] + matrix[j,j])
    
    If score_gap is provided, only return squares with |delta| <= score_gap.
    """
    # Build a fast lookup for non-zero entries
    lookup: Dict[Tuple[int, int], int] = {}
    for r, c, v in zip(matrix.row, matrix.col, matrix.data):
        lookup[(r, c)] = int(v)
    
    n = matrix.shape[0]
    squares: List[SwapSquare] = []
    
    # Check all pairs of diagonal positions
    for i in range(n):
        for j in range(i + 1, n):
            # Check if all four corners of the square are non-zero
            val_ii = lookup.get((i, i), 0)
            val_jj = lookup.get((j, j), 0)
            val_ij = lookup.get((i, j), 0)
            val_ji = lookup.get((j, i), 0)
            
            if val_ii > 0 and val_jj > 0 and val_ij > 0 and val_ji > 0:
                delta = float((val_ij + val_ji) - (val_ii + val_jj))
                
                if score_gap is None or abs(delta) <= score_gap:
                    squares.append((i, j, delta))
    
    return squares


def build_swap_graph(squares: List[SwapSquare]) -> Dict[int, Set[int]]:
    """Build adjacency graph from swap squares.
    
    Nodes are row/col indices, edges connect (i,j) if they form a swap square.
    """
    graph: Dict[int, Set[int]] = defaultdict(set)
    
    for i, j, _ in squares:
        graph[i].add(j)
        graph[j].add(i)
    
    return graph


def find_connected_components(graph: Dict[int, Set[int]]) -> List[Component]:
    """Find connected components in the swap graph using BFS."""
    visited: Set[int] = set()
    components: List[Component] = []
    
    all_nodes = set(graph.keys())
    
    for start in all_nodes:
        if start in visited:
            continue
        
        # BFS to find component
        component: List[int] = []
        queue: deque = deque([start])
        visited.add(start)
        
        while queue:
            node = queue.popleft()
            component.append(node)
            
            for neighbor in graph.get(node, set()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        
        component.sort()
        components.append(component)
    
    return components


def is_cycle(component: Component, graph: Dict[int, Set[int]]) -> bool:
    """Check if a connected component forms a cycle.
    
    A component is a cycle if:
    - Size >= 3, AND
    - There exists a path from any node back to itself
    
    For swap graphs, we check if the component forms a closed loop.
    """
    if len(component) < 3:
        return False
    
    # Check if component forms a cycle by doing DFS and detecting back edges
    visited: Set[int] = set()
    rec_stack: Set[int] = set()
    
    def has_cycle_dfs(node: int, parent: int | None = None) -> bool:
        visited.add(node)
        rec_stack.add(node)
        
        for neighbor in graph.get(node, set()):
            if neighbor not in component:
                continue
            
            if neighbor not in visited:
                if has_cycle_dfs(neighbor, node):
                    return True
            elif neighbor != parent and neighbor in rec_stack:
                # Back edge found (not to immediate parent)
                return True
        
        rec_stack.remove(node)
        return False
    
    return has_cycle_dfs(component[0])


def analyse_alternatives(
    matrix: coo_matrix,
    *,
    score_gap: float | None = None,
) -> Dict:
    """Analyse alternative assembly paths via swap square detection.
    
    Returns a dictionary containing:
    - squares: List of (i, j, delta) for all detected swap squares
    - components: List of connected components in the swap graph
    - cycles: List of components that form cycles
    - chains: List of components that are linear chains
    - ambiguity_count: Total number of swappable positions
    """
    squares = detect_swap_squares(matrix, score_gap=score_gap)
    
    if not squares:
        return {
            "squares": [],
            "components": [],
            "cycles": [],
            "chains": [],
            "ambiguity_count": 0,
        }
    
    graph = build_swap_graph(squares)
    components = find_connected_components(graph)
    
    cycles: List[Component] = []
    chains: List[Component] = []
    
    for component in components:
        if is_cycle(component, graph):
            cycles.append(component)
        else:
            chains.append(component)
    
    ambiguity_count = sum(len(comp) for comp in components)
    
    return {
        "squares": squares,
        "components": components,
        "cycles": cycles,
        "chains": chains,
        "ambiguity_count": ambiguity_count,
    }
