# Swap Square Test Design

This test is designed to ensure the swap square optimisation logic is exercised and validated.

## Matrix Structure
We want a 4x4 overlap matrix with a swap square between reads 1 and 2:

|   | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| 0 | X | A | 0 | 0 |
| 1 | B | X | 0 | 0 |
| 2 | 0 | 0 | X | C |
| 3 | 0 | 0 | D | X |

Where:
- X = self overlap (irrelevant)
- A, B = overlaps between reads 0 and 1 (swap square)
- C, D = overlaps between reads 2 and 3 (no swap square)

## Read Construction
To create a swap square, reads 0 and 1 should overlap each other in both directions with similar scores. Reads 2 and 3 should overlap only in one direction.

Example:
- read0: AAAACCCC
- read1: CCCCAAAA
- read2: GGGGTTTT
- read3: TTTTGGGG

This produces:
- read0 overlaps read1 (suffix/prefix: CCCC)
- read1 overlaps read0 (suffix/prefix: AAAA)
- read2 overlaps read3 (suffix/prefix: TTTT)
- read3 overlaps read2 (suffix/prefix: GGGG)

## FASTQ Files
### reads_1.fastq
@read0
AAAACCCC
+
IIIIIIII
@read2
GGGGTTTT
+
IIIIIIII

### reads_2.fastq
@read1
CCCCAAAA
+
IIIIIIII
@read3
TTTTGGGG
+
IIIIIIII

## Expected Swap Square
- Swap square between read0/read1
- No swap square between read2/read3

## Usage
Run the CLI with these files and --optimise-diagonal. The optimiser should detect and swap the ambiguous region for maximal diagonal sum.

---
This file documents the rationale and expected behaviour for the swap square test.