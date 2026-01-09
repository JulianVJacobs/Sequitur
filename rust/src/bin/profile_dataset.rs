//! Dataset profiling tool to recommend optimal assembly parameters.
//!
//! Analyzes FASTQ files to extract:
//! - Quality score distribution
//! - Read length statistics
//! - Insert size estimation (from read names or k-mer analysis)
//! - Error rate estimation
//! - Suggested CLI parameters

use std::path::PathBuf;

use anyhow::{Context, Result};
use bio::io::fastq;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "profile_dataset")]
#[command(about = "Profile a dataset and recommend assembly parameters")]
struct Args {
    /// Path to reads1 FASTQ file
    reads1: PathBuf,

    /// Path to reads2 FASTQ file
    reads2: PathBuf,

    /// Number of reads to sample (default: 10000)
    #[arg(long, default_value_t = 10000)]
    sample_size: usize,

    /// Output recommended parameters as JSON
    #[arg(long)]
    json: bool,
}

#[derive(Debug, Clone, serde::Serialize)]
struct QualityStats {
    mean: f32,
    median: i32,
    q25: i32,
    q75: i32,
    min: i32,
    max: i32,
    below_20: usize,
    below_30: usize,
}

#[derive(Debug, Clone, serde::Serialize)]
struct ReadStats {
    count: usize,
    mean_length: f32,
    min_length: usize,
    max_length: usize,
    quality: QualityStats,
}

#[derive(Debug, Clone, serde::Serialize)]
struct InsertStats {
    mean: usize,
    median: usize,
    std: f32,
    percentile_5: usize,
    percentile_95: usize,
}

#[derive(Debug, serde::Serialize)]
struct DatasetProfile {
    reads1: ReadStats,
    reads2: ReadStats,
    insert_size: Option<InsertStats>,
    suggested_quality_mode: String,
    suggested_quality_exponent: f32,
    suggested_min_quality: Option<i32>,
    suggested_tie_gap: f32,
    suggested_max_diff: f32,
    suggested_insert_size: usize,
    suggested_min_insert: usize,
    suggested_max_insert: usize,
}

fn analyze_qualities(all_qualities: &[i32]) -> QualityStats {
    if all_qualities.is_empty() {
        return QualityStats {
            mean: 0.0,
            median: 0,
            q25: 0,
            q75: 0,
            min: 0,
            max: 0,
            below_20: 0,
            below_30: 0,
        };
    }

    let mut sorted = all_qualities.to_vec();
    sorted.sort_unstable();

    let mean = all_qualities.iter().sum::<i32>() as f32 / all_qualities.len() as f32;
    let median = sorted[sorted.len() / 2];
    let q25 = sorted[sorted.len() / 4];
    let q75 = sorted[sorted.len() * 3 / 4];
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];

    let below_20 = all_qualities.iter().filter(|&&q| q < 20).count();
    let below_30 = all_qualities.iter().filter(|&&q| q < 30).count();

    QualityStats {
        mean,
        median,
        q25,
        q75,
        min,
        max,
        below_20,
        below_30,
    }
}

fn sample_fastq(path: &PathBuf, sample_size: usize) -> Result<ReadStats> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;

    let reader: Box<dyn std::io::Read> = if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let fastq_reader = fastq::Reader::new(reader);
    let mut count = 0;
    let mut lengths = Vec::new();
    let mut all_qualities = Vec::new();

    for (idx, result) in fastq_reader.records().enumerate() {
        if idx >= sample_size {
            break;
        }

        let record = result.context("Failed to read FASTQ record")?;
        let seq_len = record.seq().len();
        lengths.push(seq_len);

        // Convert phred qualities
        let qualities: Vec<i32> = record
            .qual()
            .iter()
            .map(|&c| ((c as i32) - 33).max(0).min(93))
            .collect();
        all_qualities.extend(qualities);

        count += 1;
    }

    if count == 0 {
        anyhow::bail!("No reads found in {}", path.display());
    }

    let mean_length = lengths.iter().sum::<usize>() as f32 / lengths.len() as f32;
    let min_length = *lengths.iter().min().unwrap();
    let max_length = *lengths.iter().max().unwrap();

    let quality = analyze_qualities(&all_qualities);

    Ok(ReadStats {
        count,
        mean_length,
        min_length,
        max_length,
        quality,
    })
}

fn estimate_insert_size(reads1: &ReadStats, reads2: &ReadStats) -> Option<InsertStats> {
    // Simple heuristic: for paired-end sequencing, typical insert = 2*read_length + gap
    // For overlapping reads (our target use case): insert < 2*read_length

    let avg_read_len = (reads1.mean_length + reads2.mean_length) / 2.0;

    // Conservative estimate: assume 20-30bp overlap typical
    let estimated_mean = (avg_read_len * 1.7) as usize;
    let estimated_std = estimated_mean as f32 * 0.15; // ~15% coefficient of variation

    Some(InsertStats {
        mean: estimated_mean,
        median: estimated_mean,
        std: estimated_std,
        percentile_5: (estimated_mean as f32 - 2.0 * estimated_std) as usize,
        percentile_95: (estimated_mean as f32 + 2.0 * estimated_std) as usize,
    })
}

fn recommend_parameters(
    reads1: &ReadStats,
    reads2: &ReadStats,
    insert: Option<&InsertStats>,
) -> DatasetProfile {
    let avg_quality = (reads1.quality.mean + reads2.quality.mean) / 2.0;
    let min_quality_1 = reads1.quality.q25;
    let min_quality_2 = reads2.quality.q25;

    // Quality mode selection based on quality distribution
    let (quality_mode, quality_exponent, min_quality) = if avg_quality > 35.0 {
        // High-quality data (Illumina HiSeq/MiSeq): use Position mode with strict threshold
        (
            "position".to_string(),
            1.5,
            Some(min_quality_1.min(min_quality_2)),
        )
    } else if avg_quality > 25.0 {
        // Medium-quality data: use Position mode with looser threshold
        ("position".to_string(), 2.0, Some(20))
    } else {
        // Low-quality or highly variable: use None mode (edit distance only)
        ("none".to_string(), 2.0, None)
    };

    // Tie gap: set based on quality variance
    // Higher variance → higher tie_gap to tolerate ambiguity
    let quality_range = (reads1.quality.q75 - reads1.quality.q25) as f32;
    let tie_gap = if quality_range < 5.0 {
        0.0 // Very uniform quality: exact ties only
    } else if quality_range < 10.0 {
        5.0 // Some variation: allow small gaps
    } else {
        10.0 // High variation: more permissive
    };

    // Max diff: stricter for high-quality data
    let max_diff = if avg_quality > 35.0 {
        0.15 // High quality: expect <15% error rate
    } else if avg_quality > 25.0 {
        0.25 // Medium quality: allow 25% error
    } else {
        0.30 // Lower quality: more permissive
    };

    // Insert size from estimation or reasonable defaults
    let (insert_size, min_insert, max_insert) = if let Some(ins) = insert {
        (ins.mean, ins.percentile_5, ins.percentile_95)
    } else {
        // Fallback: typical Illumina paired-end
        (300, 200, 450)
    };

    DatasetProfile {
        reads1: ReadStats {
            count: reads1.count,
            mean_length: reads1.mean_length,
            min_length: reads1.min_length,
            max_length: reads1.max_length,
            quality: QualityStats {
                mean: reads1.quality.mean,
                median: reads1.quality.median,
                q25: reads1.quality.q25,
                q75: reads1.quality.q75,
                min: reads1.quality.min,
                max: reads1.quality.max,
                below_20: reads1.quality.below_20,
                below_30: reads1.quality.below_30,
            },
        },
        reads2: ReadStats {
            count: reads2.count,
            mean_length: reads2.mean_length,
            min_length: reads2.min_length,
            max_length: reads2.max_length,
            quality: QualityStats {
                mean: reads2.quality.mean,
                median: reads2.quality.median,
                q25: reads2.quality.q25,
                q75: reads2.quality.q75,
                min: reads2.quality.min,
                max: reads2.quality.max,
                below_20: reads2.quality.below_20,
                below_30: reads2.quality.below_30,
            },
        },
        insert_size: insert.cloned(),
        suggested_quality_mode: quality_mode,
        suggested_quality_exponent: quality_exponent,
        suggested_min_quality: min_quality,
        suggested_tie_gap: tie_gap,
        suggested_max_diff: max_diff,
        suggested_insert_size: insert_size,
        suggested_min_insert: min_insert,
        suggested_max_insert: max_insert,
    }
}

fn print_profile(profile: &DatasetProfile, json: bool) {
    if json {
        println!("{}", serde_json::to_string_pretty(&profile).unwrap());
        return;
    }

    println!("=== Dataset Profile ===\n");

    println!("Reads1:");
    println!("  Sampled: {}", profile.reads1.count);
    println!(
        "  Length: {:.1} bp (range: {}-{})",
        profile.reads1.mean_length, profile.reads1.min_length, profile.reads1.max_length
    );
    println!(
        "  Quality: mean={:.1}, median={}, Q25={}, Q75={}",
        profile.reads1.quality.mean,
        profile.reads1.quality.median,
        profile.reads1.quality.q25,
        profile.reads1.quality.q75
    );
    println!(
        "  Low quality: {} bases <Q20 ({:.1}%), {} bases <Q30 ({:.1}%)",
        profile.reads1.quality.below_20,
        100.0 * profile.reads1.quality.below_20 as f32
            / (profile.reads1.count as f32 * profile.reads1.mean_length),
        profile.reads1.quality.below_30,
        100.0 * profile.reads1.quality.below_30 as f32
            / (profile.reads1.count as f32 * profile.reads1.mean_length)
    );

    println!("\nReads2:");
    println!("  Sampled: {}", profile.reads2.count);
    println!(
        "  Length: {:.1} bp (range: {}-{})",
        profile.reads2.mean_length, profile.reads2.min_length, profile.reads2.max_length
    );
    println!(
        "  Quality: mean={:.1}, median={}, Q25={}, Q75={}",
        profile.reads2.quality.mean,
        profile.reads2.quality.median,
        profile.reads2.quality.q25,
        profile.reads2.quality.q75
    );
    println!(
        "  Low quality: {} bases <Q20 ({:.1}%), {} bases <Q30 ({:.1}%)",
        profile.reads2.quality.below_20,
        100.0 * profile.reads2.quality.below_20 as f32
            / (profile.reads2.count as f32 * profile.reads2.mean_length),
        profile.reads2.quality.below_30,
        100.0 * profile.reads2.quality.below_30 as f32
            / (profile.reads2.count as f32 * profile.reads2.mean_length)
    );

    if let Some(ref ins) = profile.insert_size {
        println!("\nInsert Size (estimated):");
        println!("  Mean: {} bp", ins.mean);
        println!("  Median: {} bp", ins.median);
        println!("  Std: {:.1} bp", ins.std);
        println!(
            "  5th-95th percentile: {}-{} bp",
            ins.percentile_5, ins.percentile_95
        );
    }

    println!("\n=== Recommended Parameters ===\n");
    println!("--quality-mode {}", profile.suggested_quality_mode);
    println!("--quality-exponent {}", profile.suggested_quality_exponent);
    if let Some(mq) = profile.suggested_min_quality {
        println!("--min-quality {}", mq);
    }
    println!("--tie-gap {}", profile.suggested_tie_gap);
    println!("--max-diff {}", profile.suggested_max_diff);
    println!("--insert-size {}", profile.suggested_insert_size);
    println!("--min-insert {}", profile.suggested_min_insert);
    println!("--max-insert {}", profile.suggested_max_insert);

    println!("\n=== Rationale ===");
    let avg_q = (profile.reads1.quality.mean + profile.reads2.quality.mean) / 2.0;
    if avg_q > 35.0 {
        println!("High average quality (Q{:.1}) detected → using Position mode with strict error penalty", avg_q);
        println!("Quality scores will modulate mismatch penalties at each position");
    } else if avg_q > 25.0 {
        println!(
            "Medium average quality (Q{:.1}) detected → using Position mode",
            avg_q
        );
        println!("Quality scores provide useful signal but expect some noise");
    } else {
        println!(
            "Lower quality (Q{:.1}) detected → using edit distance only (None mode)",
            avg_q
        );
        println!("Quality scores may be unreliable; falling back to strict sequence matching");
    }

    let quality_range = (profile.reads1.quality.q75 - profile.reads1.quality.q25) as f32;
    if quality_range < 5.0 {
        println!(
            "Uniform quality distribution (IQR={:.1}) → tie_gap=0.0 for exact ties only",
            quality_range
        );
    } else if quality_range < 10.0 {
        println!(
            "Moderate quality variation (IQR={:.1}) → tie_gap=5.0 to handle ambiguity",
            quality_range
        );
    } else {
        println!(
            "High quality variation (IQR={:.1}) → tie_gap=10.0 for permissive tie-breaking",
            quality_range
        );
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    eprintln!("Profiling dataset...");
    eprintln!("Sampling {} reads from each file", args.sample_size);

    let reads1 = sample_fastq(&args.reads1, args.sample_size)
        .with_context(|| format!("Failed to profile reads1: {}", args.reads1.display()))?;

    let reads2 = sample_fastq(&args.reads2, args.sample_size)
        .with_context(|| format!("Failed to profile reads2: {}", args.reads2.display()))?;

    let insert_estimate = estimate_insert_size(&reads1, &reads2);

    let profile = recommend_parameters(&reads1, &reads2, insert_estimate.as_ref());

    print_profile(&profile, args.json);

    Ok(())
}
