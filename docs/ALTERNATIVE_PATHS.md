# Alternative Path Detection

This feature detects and reports alternative assembly paths and cycles in the Sequitur overlap graph using swap-square analysis of the lower-diagonal matrix.

## Theory

After solving the maximum bipartite matching to create a lower-diagonal matrix ordering, alternative assembly paths can exist when:

1. **Swap Squares**: A 2×2 submatrix with all four corners non-zero indicates two rows/columns can be swapped while maintaining a valid lower-diagonal matrix:
   ```
   matrix[i,i] ≠ 0, matrix[j,j] ≠ 0
   matrix[i,j] ≠ 0, matrix[j,i] ≠ 0
   ```

2. **Cycles**: Connected chains of swap squares that form a closed loop, where positions can be freely rotated.

3. **Chains**: Linear sequences of swappable positions without cycles.

## Usage

### Python

```bash
python sequitur.py reads1.fastq reads2.fastq \
  --analyse-alternatives \
  --score-gap 5.0 \
  --alternatives-json results/alternatives.json \
  --output-fasta results/assembled.fasta
```

Options:
- `--analyse-alternatives`: Enable alternative path detection
- `--score-gap FLOAT`: Only report swap squares with score delta ≤ gap (optional)
- `--alternatives-json PATH`: Write analysis results to JSON file

### Rust

```bash
./sequitur_rs reads1.fastq reads2.fastq \
  --analyse-alternatives \
  --score-gap 5.0 \
  --alternatives-json results/alternatives.json \
  --output-fasta results/assembled.fasta \
  --verbose
```

## Output Format

JSON output structure:
```json
{
  "squares": [
    {"i": 0, "j": 1, "delta": -2.5}
  ],
  "components": [[0, 1, 2]],
  "cycles": [[0, 1, 2]],
  "chains": [[3, 4]],
  "ambiguity_count": 5
}
```

- **squares**: List of swap squares with score delta `(M[i,j] + M[j,i]) - (M[i,i] + M[j,j])`
- **components**: Connected components in the swap graph
- **cycles**: Components that form closed loops (size ≥ 3)
- **chains**: Linear swappable segments
- **ambiguity_count**: Total positions involved in ambiguous regions

## Implementation

### Python
- `sequitur_core/alternative_paths.py`: Core detection algorithms
- `sequitur.py`: CLI integration
- `tests/test_alternative_paths.py`: Unit tests

### Rust
- `src/alternative_paths.rs`: Core detection algorithms  
- `src/main.rs`: CLI integration
- Built-in module tests

## Complexity

For a sparse matrix with *n* positions and *m* non-zero entries:
- **Swap square detection**: O(n² · d) where d is average diagonal density
- **Component detection**: O(n + e) where e is number of swap edges
- **Cycle detection**: O(n + e) via DFS

Sparse matrix representation keeps this efficient even for large assemblies.

## Swap-Square Guard in Path Traversal

The Rust assembler uses swap-square detection to enable implicit backtracking during path construction.

### Motivation

When building a path through the overlap graph, traditional greedy selection (argmax) can fail when:
1. Multiple overlaps have equal or near-equal scores
2. The best choice leads to a dead-end while a slightly worse choice continues successfully
3. The assembly needs to "reconsider" a previous decision without explicit recursion

### Algorithm

The path-building loop uses a **masked next-best selector (`argnb`)** with swap-square safety:

1. **Forward Selection**: For each current node, select the highest-scoring unvisited target using `argnb`
2. **Backward Detection**: Check if the selected target is already in the path (but not the immediate predecessor)
3. **Interchangeability Guard**: 
   - If backward: check `is_swap_square(current, target, matrix)` 
   - Swap-square means the two reads can be assembled in either order → safe to allow backward pick
   - No swap-square means reads are ordered-dependent → reject and mask target
4. **Masking & Retry**: If rejected, add target to `tried` set and call `argnb` again to get next-best candidate
5. **Termination**: When no more candidates available, append remaining reads by column sum

### Benefits

- **Memory Efficient**: Masked re-scan (O(k·deg)) instead of allocating sorted candidate lists
- **Deterministic**: Explicit tie-breakers (weight desc → overlap desc → index asc)
- **No Recursion**: Handles ambiguity without backtracking stack or global re-analysis
- **Conservative**: Only allows backward picks when proven safe (swap-square exists)

### Example

Given overlap graph where read A has edges to B and C with equal scores, and B↔C form a swap-square:
```
A → B (weight=10, overlap=5)
A → C (weight=10, overlap=5)
B ↔ C (swap-square)
```

1. Start at A, pick B (arbitrary tie-break by index)
2. From B, consider C (backward dependency)
3. Check `is_swap_square(B, C)` → true (swap-square exists)
4. Accept C, continue path: A → B → C

Without swap-square guard, the algorithm would reject C and potentially terminate early or produce suboptimal assembly.

### Implementation

- **Function**: `find_lower_diagonal_path` in `rust/src/matching.rs`
- **Helper**: `is_swap_square(i, j, csc)` checks all four corners non-zero
- **Selector**: `argnb(row_edges, mask, overlaps, visited)` returns best unmasked candidate
- **Tests**: `matching::tests::allows_interchangeable_backward_pick`, `rejects_non_interchangeable_backward_pick`

