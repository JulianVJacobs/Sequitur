# Swap Square Test Design

This test is designed to ensure the swap square optimisation logic is exercised and validated.

## Actual Read Sequences
```
0 read0: TTAAAGTTAAAAAAAAAAAAAAGGGAAACCC
1 read3: CCCTTTTTAAAGTTTTTTTTTGGGGAAACCC
2 read1: GAAACCCCCCCCCCCCCCCCTTTTTAAAGTT
3 read2: GGGAAACCCGGGGGGGGGGGGGGGGGGGGGG
```

## Expected Adjacency Matrix (Overlap Scores)

|   |  0  |  1  |  2  |  3  |
|---|-----|-----|-----|-----|
| 0 |  31 |   0 |   8 |   0 |
| 1 |   3 |  31 |  14 |   0 |
| 2 |   7 |   7 |  31 |   0 |
| 3 |   9 |   9 |   0 |  31 |

Where:
- Diagonal (31) = read length (self-overlap, forbidden in cost matrix)
- 0 = no overlap exists between these reads
- Non-zero off-diagonal = overlap score (quality-adjusted)

Key overlaps:
- 0→2: score 8 (overlap length 7)
- 1→0: score 3
- 1→2: score 14 (overlap length 14)
- 2→0: score 7, 2→1: score 7 (swap square: rows 1 and 2 both overlap)
- 3→0: score 9, 3→1: score 9

## FASTQ Files
### reads_1.fastq
```
@read0
TTAAAGTTAAAAAAAAAAAAAAGGGAAACCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3
CCCTTTTTAAAGTTTTTTTTTGGGGAAACCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

### reads_2.fastq
```
@read1
GAAACCCCCCCCCCCCCCCCTTTTTAAAGTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
GGGAAACCCGGGGGGGGGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## Expected Cost Matrix (for LAPJV)

After negating scores and marking missing edges as INF:

|   |  0  |  1  |  2  |  3  |
|---|-----|-----|-----|-----|
| 0 |  31 | inf |  -8 | inf |
| 1 |  -3 |  31 | -14 | inf |
| 2 |  -7 |  -7 |  31 | inf |
| 3 |  -9 |  -9 | inf |  31 |

- Diagonal = 31 (or INF to forbid self-assignments)
- Missing edges = INF (forbidden)
- Real overlaps = negative score (higher overlap → lower cost)

## Expected Assembly Path

LAPJV should find the minimum-cost perfect matching. With the above costs, the optimal path is likely `[0, 2, 1, 3]` assembling:
- 0 → 2 (cost -8)
- 2 → 1 (cost -7)  
- 1 self or 3 terminus

## Expected Swap Square
- Rows 1 and 2 both have overlaps to columns 0 and 1
- This creates potential ambiguity in assembly order

## Usage
Run the CLI with these files and --optimise-diagonal. The optimiser should detect and swap the ambiguous region for maximal diagonal sum.

---
This file documents the rationale and expected behaviour for the swap square test.