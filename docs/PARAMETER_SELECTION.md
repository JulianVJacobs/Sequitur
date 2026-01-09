# Deterministic Parameter Selection

## Problem

Assembly quality depends on several CLI parameters that historically required manual tuning:
- `--quality-mode`, `--quality-exponent`, `--min-quality` (quality scoring)
- `--tie-gap` (ambiguity tolerance)
- `--max-diff` (overlap error threshold)
- `--insert-size`, `--min-insert`, `--max-insert` (mate pair constraints)
- `--mate-penalty-weight` (mate-aware tie-breaking)

Manual selection involves guessing and iterative refinement, which is inefficient and non-reproducible.

## Solution: Dataset Profiling

Use the `profile_dataset` tool to analyze your input data and receive evidence-based parameter recommendations.

### Usage

```bash
# Profile your dataset
cargo run --release --bin profile_dataset -- \
  reads1.fastq.gz reads2.fastq.gz \
  --sample-size 10000

# Or use the compiled binary
./target/release/profile_dataset reads1.fastq.gz reads2.fastq.gz
```

### What It Analyzes

**Quality Score Distribution:**
- Mean, median, quartiles (Q25/Q75)
- Fraction of bases below Q20 and Q30 thresholds
- Inter-quartile range (IQR) as measure of quality variance

**Read Length Statistics:**
- Mean, min, max for both read sets
- Used for insert size estimation

**Insert Size Estimation:**
- Heuristic: typical paired-end overlapping reads have insert ≈ 1.7× read length
- Calculates 5th-95th percentile bounds assuming 15% coefficient of variation

### Decision Logic

#### Quality Mode Selection
```
if avg_quality > 35.0:
    mode = "position"
    exponent = 1.5
    min_quality = Q25 (first quartile)
    # High-quality Illumina (HiSeq, MiSeq, NovaSeq S4)
    # Quality scores are reliable; use strict per-position weighting

elif avg_quality > 25.0:
    mode = "position"  
    exponent = 2.0
    min_quality = 20
    # Medium-quality data (NovaSeq S1/S2, older platforms)
    # Quality provides useful signal with some noise

else:
    mode = "none"
    exponent = 2.0
    min_quality = null
    # Low or unreliable quality
    # Fall back to edit distance only
```

#### Tie Gap (Ambiguity Tolerance)
```
quality_IQR = Q75 - Q25

if quality_IQR < 5.0:
    tie_gap = 0.0
    # Uniform quality → confident scoring → exact ties only

elif quality_IQR < 10.0:
    tie_gap = 5.0
    # Moderate variation → allow small score gaps

else:
    tie_gap = 10.0
    # High variation → permissive tie-breaking needed
```

#### Max Diff (Error Threshold)
```
if avg_quality > 35.0:
    max_diff = 0.15  # Expect <15% error rate

elif avg_quality > 25.0:
    max_diff = 0.25  # Allow 25% error

else:
    max_diff = 0.30  # More permissive for noisy data
```

#### Insert Size Bounds
```
estimated_mean = avg_read_length × 1.7
estimated_std = estimated_mean × 0.15

insert_size = estimated_mean
min_insert = estimated_mean - 2×std  (5th percentile)
max_insert = estimated_mean + 2×std  (95th percentile)
```

### Example Output

```
=== Dataset Profile ===

Reads1:
  Sampled: 10000
  Length: 150.0 bp (range: 150-150)
  Quality: mean=37.2, median=38, Q25=34, Q75=40
  Low quality: 1243 bases <Q20 (0.8%), 12405 bases <Q30 (8.3%)

Reads2:
  Sampled: 10000
  Length: 150.0 bp (range: 150-150)
  Quality: mean=36.8, median=37, Q25=33, Q75=40
  Low quality: 1456 bases <Q20 (1.0%), 14201 bases <Q30 (9.5%)

Insert Size (estimated):
  Mean: 255 bp
  Median: 255 bp
  Std: 38.2 bp
  5th-95th percentile: 179-331 bp

=== Recommended Parameters ===

--quality-mode position
--quality-exponent 1.5
--min-quality 33
--tie-gap 5.0
--max-diff 0.15
--insert-size 255
--min-insert 179
--max-insert 331

=== Rationale ===
High average quality (Q37.0) detected → using Position mode with strict error penalty
Quality scores will modulate mismatch penalties at each position
Moderate quality variation (IQR=6.5) → tie_gap=5.0 to handle ambiguity
```

## Advanced: Insert Size from Alignment

For more accurate insert size estimation, align a sample to a reference:

```bash
# Extract sample
seqtk sample -s100 reads1.fastq.gz 10000 > sample_r1.fq
seqtk sample -s100 reads2.fastq.gz 10000 > sample_r2.fq

# Align with minimap2
minimap2 -ax sr reference.fa sample_r1.fq sample_r2.fq \
  | samtools view -bS - \
  | samtools sort -o sample.bam -
samtools index sample.bam

# Extract insert sizes
samtools view -f 0x2 sample.bam \
  | awk '{if ($9>0) print $9}' \
  | sort -n \
  | awk 'BEGIN{c=0;sum=0;}{a[c++]=$1;sum+=$1;}END{
      mean=sum/c;
      p5=a[int(c*0.05)];
      p50=a[int(c*0.5)];
      p95=a[int(c*0.95)];
      print "mean="mean" median="p50" min="p5" max="p95;
    }'
```

Pass these values to sequitur:
```bash
--insert-size <mean> --min-insert <p5> --max-insert <p95>
```

## When to Override Recommendations

**Known library characteristics:**
If you know your library's exact insert size distribution from library prep protocols, use those values directly.

**Reference-based validation:**
If you have a reference genome:
1. Run assembly with recommended parameters
2. Align assembly to reference (minimap2)
3. Check identity/coverage in PAF output
4. Adjust `max_diff` or `tie_gap` if needed

**Iterative refinement:**
For critical applications, run a parameter sweep:
```bash
for qmode in none position confidence; do
  for tie in 0.0 5.0 10.0; do
    sequitur reads1.fq reads2.fq \
      --quality-mode $qmode \
      --tie-gap $tie \
      --output-fasta asm_${qmode}_${tie}.fa
  done
done
# Compare assemblies against reference or via N50/contiguity metrics
```

## Integration with Slurm Scripts

Update your job script to use profiling:

```bash
# Profile dataset first
PROFILE_OUTPUT=$(./profile_dataset "$READ1" "$READ2")

# Extract recommended parameters
QUALITY_MODE=$(echo "$PROFILE_OUTPUT" | grep "^--quality-mode" | awk '{print $2}')
TIE_GAP=$(echo "$PROFILE_OUTPUT" | grep "^--tie-gap" | awk '{print $2}')
# ... etc

# Run assembly with profiled parameters
sequitur "$READ1" "$READ2" \
  --quality-mode "$QUALITY_MODE" \
  --tie-gap "$TIE_GAP" \
  # ... other recommended params
```

Or save profile to JSON and parse:
```bash
./profile_dataset reads1.fq reads2.fq --json > profile.json
QUALITY_MODE=$(jq -r '.suggested_quality_mode' profile.json)
```

## Summary

**Deterministic approach:**
1. Run `profile_dataset` on your inputs
2. Use recommended parameters directly
3. (Optional) Validate against reference if available
4. Document profiling output for reproducibility

**Benefits:**
- Evidence-based parameter selection
- Reproducible across datasets
- Adapts to data quality automatically
- No manual guessing or iteration required
