"""Tests for alternative path detection."""

import numpy as np
from scipy.sparse import coo_matrix

from sequitur_core.alternative_paths import (
    analyse_alternatives,
    build_swap_graph,
    detect_swap_squares,
    find_connected_components,
    is_cycle,
)


def test_detects_simple_swap_square():
    """Test detection of a single swap square."""
    # Create a matrix with a swap square at (0,1)
    row = [0, 1, 0, 1, 2, 3]
    col = [0, 1, 1, 0, 2, 3]
    data = [10, 8, 7, 9, 5, 6]
    matrix = coo_matrix((data, (row, col)), shape=(4, 4))
    
    squares = detect_swap_squares(matrix)
    
    assert len(squares) == 1
    assert squares[0][0] == 0
    assert squares[0][1] == 1
    # delta = (7 + 9) - (10 + 8) = -2
    assert abs(squares[0][2] - (-2.0)) < 0.001


def test_detects_cycle_in_component():
    """Test detection of a 3-node cycle."""
    # Create a graph with a cycle: 0 <-> 1 <-> 2 <-> 0
    # This requires all 6 swap squares for a full 3-cycle
    row = [0, 1, 2, 0, 1, 2, 0, 1, 2, 1, 2, 0]
    col = [0, 1, 2, 1, 2, 0, 2, 0, 1, 1, 2, 0]
    data = [10, 10, 10, 8, 8, 8, 9, 9, 9, 10, 10, 10]
    matrix = coo_matrix((data, (row, col)), shape=(3, 3))
    
    squares = detect_swap_squares(matrix)
    graph = build_swap_graph(squares)
    components = find_connected_components(graph)
    
    assert len(components) == 1
    assert len(components[0]) == 3
    assert is_cycle(components[0], graph)


def test_detects_chain_not_cycle():
    """Test detection of a linear chain without a cycle."""
    # Create a linear chain: 0 <-> 1 <-> 2
    row = [0, 1, 2, 0, 1, 1, 2]
    col = [0, 1, 2, 1, 0, 2, 1]
    data = [10, 10, 10, 8, 8, 8, 8]
    matrix = coo_matrix((data, (row, col)), shape=(3, 3))
    
    squares = detect_swap_squares(matrix)
    graph = build_swap_graph(squares)
    components = find_connected_components(graph)
    
    assert len(components) == 1
    assert len(components[0]) == 3
    assert not is_cycle(components[0], graph)


def test_analyse_alternatives_full():
    """Test the full analysis pipeline."""
    # Create a matrix with both a chain and isolated squares
    row = [0, 1, 2, 3, 0, 1, 1, 2]
    col = [0, 1, 2, 3, 1, 0, 2, 1]
    data = [10, 10, 10, 5, 8, 8, 8, 8]
    matrix = coo_matrix((data, (row, col)), shape=(4, 4))
    
    result = analyse_alternatives(matrix)
    
    assert len(result["squares"]) > 0
    assert len(result["components"]) > 0
    assert result["ambiguity_count"] > 0


def test_score_gap_filtering():
    """Test that score_gap filters squares correctly."""
    row = [0, 1, 0, 1, 2, 3]
    col = [0, 1, 1, 0, 2, 3]
    data = [10, 8, 20, 1, 5, 6]  # Delta = (20 + 1) - (10 + 8) = 3
    matrix = coo_matrix((data, (row, col)), shape=(4, 4))
    
    # Without filtering
    squares_all = detect_swap_squares(matrix, score_gap=None)
    assert len(squares_all) == 1
    
    # With strict filtering (delta is 3)
    squares_strict = detect_swap_squares(matrix, score_gap=2.0)
    assert len(squares_strict) == 0
    
    # With permissive filtering
    squares_permissive = detect_swap_squares(matrix, score_gap=5.0)
    assert len(squares_permissive) == 1


if __name__ == "__main__":
    test_detects_simple_swap_square()
    test_detects_cycle_in_component()
    test_detects_chain_not_cycle()
    test_analyse_alternatives_full()
    test_score_gap_filtering()
    print("All tests passed!")
