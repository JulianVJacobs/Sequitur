# Quality Score Integration in Sequitur

## Overview

Sequitur now supports phred quality scores from FASTQ files to improve overlap detection accuracy. This document describes the quality integration system, its design, and how to use it.

## What are Phred Quality Scores?

Phred scores (Q) quantify the confidence of a base call using:

$$Q = -10 \log_{10}(P)$$

where $P$ is the probability of an incorrect base call.

**Common interpretations:**
- Q=10: 1 in 10 error probability (90% accuracy)
- Q=20: 1 in 100 error probability (99% accuracy)
- Q=30: 1 in 1000 error probability (99.9% accuracy)
- Q=40: 1 in 10000 error probability (99.99% accuracy)

## What Phred Scores Tell Us

**What they capture:**
- Sequencer confidence in a single base call
- Position-level confidence from signal strength
- Per-read base-by-base quality profile

**What they do NOT capture:**
- Biological variation vs. sequencing error (they only measure sequencer confidence, not truth)
- Systematic bias (homopolymer runs, GC extremes, repetitive sequences)
- Correlation between adjacent positions
- Platform-specific bias (sequencers can be over/under-calibrated)

## Integration Strategy

### Three Scoring Modes

1. **None (default)**: No quality weighting
   - Uses current behavior: score = overlap_length × (1 - error_rate)^exponent
   - Backward compatible; works with FASTA input

2. **Position**: Soft mismatch penalty scaled by quality
   - Mismatches at low-quality positions receive reduced penalty
   - Mismatches at high-quality positions receive full penalty
   - Formula: score -= error_probability × exponent per mismatch

3. **Confidence**: Phred-explicit scoring
   - Matches at high-quality positions boost score more
   - Mismatches at high-quality positions penalize more
   - Formula: score += confidence per match, score -= confidence per mismatch

### How Quality Affects Assembly

**Example: Resolving conflicting overlaps**

Consider two 50bp overlaps with 2 mismatches each:

- **Overlap A**: Mismatches at Q=5 (67% error) and Q=30 (0.1% error)
  - Quality-blind score: 48 (penalizes both equally)
  - Quality-aware score: ~49 (Q=5 mismatch is likely noise; Q=30 is likely real SNP)

- **Overlap B**: Both mismatches at Q=30 (0.1% error each)
  - Quality-blind score: 48
  - Quality-aware score: ~47 (both high-confidence, likely true divergence)

With quality weighting, Overlap A scores higher, suggesting it's more likely correct.

## Configuration

### CLI Flags

```bash
# Quality-disabled (default)
cargo run -- reads1.fastq reads2.fastq

# Quality-aware with Position mode
cargo run -- reads1.fastq reads2.fastq --quality-mode position

# Quality-aware with Confidence mode
cargo run -- reads1.fastq reads2.fastq --quality-mode confidence

# Set quality penalty exponent (higher = stricter)
cargo run -- reads1.fastq reads2.fastq --quality-mode confidence --quality-exponent 1.5

# Set minimum quality threshold (reject bases below Q15)
cargo run -- reads1.fastq reads2.fastq --quality-mode position --min-quality 15

# Enable quality logging
cargo run -- reads1.fastq reads2.fastq --quality-mode confidence --log-quality-scores
```

### Programmatic API

```rust
use sequitur::overlap::{OverlapConfig, QualityConfig, QualityScoring};
use sequitur::read_source::load_fastq_pair_with_quality;

let reads = load_fastq_pair_with_quality("reads1.fastq", "reads2.fastq")?;

let mut config = OverlapConfig::default();
config.quality_config = QualityConfig {
    scoring_mode: QualityScoring::Confidence,
    error_penalty_exponent: 1.5,
    min_quality: Some(10),
    use_quality_for_consensus: true,
    log_quality_scores: true,
};

let result = create_overlap_graph_with_quality(&reads, config)?;
```

## Known Limitations & Assumptions

### Sequencer Calibration

Phred scores assume the sequencer's base-caller has been trained on known-truth sequences. In practice:

- **Modern Illumina**: Generally well-calibrated (within ±2 points)
- **PacBio early versions**: Often under-calibrated (reported Q too low)
- **Oxford Nanopore**: Variable calibration across runs

**Recommendation**: If mixing reads from multiple platforms, consider a recalibration step (not yet implemented in Sequitur).

### Independence Assumption

Phred scores assume base-call errors are independent. In reality:

- Homopolymer runs show correlated errors
- Systematic biases (e.g., AT-rich regions) increase error rates across stretches
- Adapter sequences have different error profiles

**Implication**: Quality-aware scoring improves on average, but may miss systematic patterns.

### Biological Divergence vs. Sequencing Error

High-quality mismatches suggest:
1. Real biological variation (SNPs, indels)
2. Paralog amplification
3. Contamination
4. Polyploidy

Sequitur's quality scoring cannot distinguish these cases. Users should validate high-divergence assemblies.

## Quality Metrics Output

When `--log-quality-scores` is enabled, Sequitur outputs:

```
Quality scoring (mode=Confidence): X edges removed, Y low-Q positions detected
Mean quality: 35.2 (range: 3-40)
Score histogram (blind): [count at Q=0-1, count at Q=1-2, ...]
Score histogram (adjusted): [...]
```

This allows you to inspect the quality distribution and verify that quality had expected impact.

## Backward Compatibility

- **FASTA input**: Quality defaults to None mode (no quality weighting)
- **FASTQ without quality weighting**: Identical to v0.1.0 behavior
- **Mixed FASTA+FASTQ**: FASTA reads have no quality; FASTQ reads use supplied quality

No configuration changes required for existing workflows.

## Future Enhancements

### Phase 2: Quality Recalibration
- Detect over/under-calibration by comparing observed error rate to reported phred
- Optional recalibration to normalize scores across platforms

### Phase 3: Systematic Bias Detection
- Fit position-dependent error rates
- Flag high-error regions (homopolymers, GC extremes)
- Optional trimming of unreliable regions

### Phase 4: Contig Polishing
- Use quality scores to resolve ambiguous bases in consensus
- Bayesian-style posterior probability estimation

## Testing & Validation

For regression tests, Sequitur includes:

1. **Synthetic FASTQ**: Datasets with known quality gradients
2. **Golden edges**: Pre-computed overlap matrices for reference
3. **Parallel test**: Running same input with/without quality to verify topology match
4. **Benchmark**: Tracking quality-scoring overhead (target: <5% vs. blind mode)

Run tests with:
```bash
cargo test quality
cargo test --release -- --nocapture  # Show metrics output
```

## References

- Ewing, B., & Green, P. (1998). Base-calling of automated sequencer traces using Phred. Genome Research, 8(3), 175-185.
- Cock, P. J., et al. (2009). The Sanger FASTQ file format for sequences with quality scores. Nucleic Acids Research, 38(6), e181.
- Phrap documentation: http://phrap.org/

## Example Workflow

```bash
# 1. Assemble with quality awareness
cargo run --release -- \
    reads1.fastq reads2.fastq \
    --quality-mode confidence \
    --quality-exponent 1.5 \
    --output-fasta assembly.fasta \
    --metrics-csv metrics.json

# 2. Inspect quality impact
cat metrics.json | jq .quality

# 3. Compare to quality-blind run
cargo run --release -- \
    reads1.fastq reads2.fastq \
    --output-fasta assembly_blind.fasta

# 4. Validate that quality improved specificity
# (fewer false positives, same true positives)
minimap2 -c assembly.fasta assembly_blind.fasta
```
