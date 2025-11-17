use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use bio::alignment::distance::levenshtein;
use bio::alphabets::dna;
use bio::io::{fasta, fastq};
use clap::Parser;
use flate2::read::MultiGzDecoder;

use sequitur_rs::{
    adjacency_to_csc, analyse_alternatives, create_overlap_graph, find_lower_diagonal_path,
    AffixArray, OverlapConfig,
};

/// Sequitur Rust CLI (prototype)
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// FASTQ/FASTA file for read set 1
    reads1: String,

    /// FASTQ/FASTA file for read set 2 (reverse-complemented during ingest)
    reads2: String,

    /// Optional reference FASTA used for post-assembly confirmation
    #[arg(long)]
    reference: Option<String>,
    /// Optional output FASTA path for the assembled sequence
    #[arg(long)]
    output_fasta: Option<String>,

    /// Optional CSV file to append assembly metrics (reserved for future use)
    #[arg(long)]
    metrics_csv: Option<String>,

    /// Detect and report alternative assembly paths and cycles
    #[arg(long)]
    analyse_alternatives: bool,

    /// Maximum score gap for alternative paths (default: no filter)
    #[arg(long)]
    score_gap: Option<f64>,

    /// Output JSON file for alternative path analysis
    #[arg(long)]
    alternatives_json: Option<String>,

    /// Verbose output (default: quiet)
    #[arg(long, short = 'v')]
    verbose: bool,

    /// Skip reverse complement of reads2 (for non-genomic or unpaired data)
    #[arg(long)]
    no_revcomp: bool,
}

fn main() {
    env_logger::init();
    let args = Args::parse();

    if args.verbose {
        println!("Sequitur Rust prototype");
        println!("reads1: {}", args.reads1);
        println!("reads2: {}", args.reads2);
        if let Some(refp) = &args.reference {
            println!("reference: {}", refp);
        }
    }

    if let Err(error) = run_pipeline(
        &args.reads1,
        &args.reads2,
        args.output_fasta.as_deref(),
        args.reference.as_deref(),
        args.analyse_alternatives,
        args.score_gap,
        args.alternatives_json.as_deref(),
        args.verbose,
        args.no_revcomp,
    ) {
        eprintln!("Assembly failed: {error:?}");
        std::process::exit(1);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SequenceFormat {
    Fastq,
    Fasta,
    Lines,
}

fn is_gzip(path: &Path) -> bool {
    path.extension()
        .map(|ext| ext.eq_ignore_ascii_case("gz") || ext.eq_ignore_ascii_case("bgz"))
        .unwrap_or(false)
}

fn infer_format(path: &Path) -> SequenceFormat {
    let mut ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    if ext == "gz" || ext == "bgz" {
        if let Some(stem) = path.file_stem() {
            let stem_path = Path::new(stem);
            ext = stem_path
                .extension()
                .and_then(|e| e.to_str())
                .map(|s| s.to_ascii_lowercase())
                .unwrap_or_default();
        }
    }

    match ext.as_str() {
        "fastq" | "fq" => SequenceFormat::Fastq,
        "fasta" | "fa" | "fna" => SequenceFormat::Fasta,
        _ => SequenceFormat::Lines,
    }
}

fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    if is_gzip(path) {
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn uppercase_sequence(bytes: &[u8]) -> Result<String> {
    let upper = bytes
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect::<Vec<u8>>();
    String::from_utf8(upper).map_err(|_| anyhow!("Encountered non-UTF-8 symbols in sequence data"))
}

fn read_sequences(path: &Path) -> Result<Vec<String>> {
    let format = infer_format(path);
    let reader = open_reader(path)?;

    match format {
        SequenceFormat::Fastq => {
            let fastq_reader = fastq::Reader::new(reader);
            let mut sequences = Vec::new();
            for record in fastq_reader.records() {
                let record = record.with_context(|| {
                    format!("Error reading FASTQ record from {}", path.display())
                })?;
                let seq = uppercase_sequence(record.seq())?;
                sequences.push(seq);
            }
            Ok(sequences)
        }
        SequenceFormat::Fasta => {
            let fasta_reader = fasta::Reader::new(reader);
            let mut sequences = Vec::new();
            for record in fasta_reader.records() {
                let record = record.with_context(|| {
                    format!("Error reading FASTA record from {}", path.display())
                })?;
                let seq = uppercase_sequence(record.seq())?;
                sequences.push(seq);
            }
            Ok(sequences)
        }
        SequenceFormat::Lines => {
            let mut sequences = Vec::new();
            let mut buf_reader = reader;
            loop {
                let mut line = String::new();
                let bytes = buf_reader.read_line(&mut line)?;
                if bytes == 0 {
                    break;
                }
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                sequences.push(trimmed.to_ascii_uppercase());
            }
            Ok(sequences)
        }
    }
}

fn reverse_complement_all(reads: Vec<String>) -> Result<Vec<String>> {
    reads
        .into_iter()
        .map(|seq| {
            let rc = dna::revcomp(seq.as_bytes());
            String::from_utf8(rc).map_err(|_| anyhow!("Reverse complement produced invalid UTF-8"))
        })
        .collect()
}

fn load_reference(path: &Path) -> Result<Option<String>> {
    if !path.exists() {
        bail!("Reference path {} does not exist", path.display());
    }
    let sequences = read_sequences(path)?;
    if sequences.is_empty() {
        return Ok(None);
    }
    let mut combined = String::new();
    for seq in sequences {
        combined.push_str(&seq);
    }
    Ok(Some(combined))
}

fn run_pipeline(
    reads1_path: &str,
    reads2_path: &str,
    output_fasta: Option<&str>,
    reference_path: Option<&str>,
    analyse_alts: bool,
    score_gap: Option<f64>,
    alternatives_json: Option<&str>,
    verbose: bool,
    no_revcomp: bool,
) -> Result<String> {
    let reads1 = read_sequences(Path::new(reads1_path))
        .with_context(|| format!("Failed to parse reads from {}", reads1_path))?;
    let reads2_raw = read_sequences(Path::new(reads2_path))
        .with_context(|| format!("Failed to parse reads from {}", reads2_path))?;

    let reads2 = if no_revcomp {
        if verbose {
            println!("Skipping reverse complement (--no-revcomp flag)");
        }
        reads2_raw
    } else {
        if verbose {
            println!("Applying reverse complement to reads2");
        }
        reverse_complement_all(reads2_raw)?
    };

    let reads1_count = reads1.len();
    let reads2_count = reads2.len();

    let mut reads = Vec::with_capacity(reads1_count + reads2_count);
    reads.extend(reads1);
    reads.extend(reads2);

    if reads.is_empty() {
        bail!(
            "No reads were parsed from {} or {}",
            reads1_path,
            reads2_path
        );
    }

    if verbose {
        eprintln!(
            "Loaded {} reads ({} from reads1, {} reverse-complemented from reads2)",
            reads.len(),
            reads1_count,
            reads2_count
        );
    }

    let reference = if let Some(path) = reference_path {
        load_reference(Path::new(path))?
    } else {
        None
    };

    let affix = AffixArray::build(&reads, 3);
    if verbose {
        eprintln!("Built affix array with {} entries", affix.len());
    }

    let config = OverlapConfig::default();
    let (_affix_array, adjacency_matrix, overlap_matrix) =
        create_overlap_graph(&reads, Some(affix), config);

    let adjacency_csc = adjacency_matrix.to_csc();
    let overlap_csc = overlap_matrix.to_csc();

    if verbose {
        eprintln!("Created adjacency CSC with {} edges", adjacency_csc.nnz());
    }

    let assembled = find_lower_diagonal_path(&adjacency_csc, &overlap_csc, &reads, None);

    // Analyse alternative paths if requested
    if analyse_alts {
        let analysis = analyse_alternatives(&adjacency_csc, score_gap);

        if verbose || alternatives_json.is_none() {
            eprintln!("\nAlternative Path Analysis:");
            eprintln!("Swap squares detected : {}", analysis.squares.len());
            eprintln!("Ambiguous components  : {}", analysis.components.len());
            eprintln!("Cycles detected       : {}", analysis.cycles.len());
            eprintln!("Linear chains         : {}", analysis.chains.len());
            eprintln!("Total ambiguous pos   : {}", analysis.ambiguity_count);

            if !analysis.cycles.is_empty() {
                eprintln!("\nCycle positions:");
                for (idx, cycle) in analysis.cycles.iter().enumerate() {
                    eprintln!("  Cycle {}: {:?}", idx + 1, cycle);
                }
            }

            if !analysis.chains.is_empty() {
                eprintln!("\nChain positions:");
                for (idx, chain) in analysis.chains.iter().enumerate() {
                    eprintln!("  Chain {}: {:?}", idx + 1, chain);
                }
            }
        }

        if let Some(json_path) = alternatives_json {
            use serde_json::json;

            let squares_json: Vec<_> = analysis
                .squares
                .iter()
                .map(|sq| json!({"i": sq.i, "j": sq.j, "delta": sq.delta}))
                .collect();

            let output = json!({
                "squares": squares_json,
                "components": analysis.components,
                "cycles": analysis.cycles,
                "chains": analysis.chains,
                "ambiguity_count": analysis.ambiguity_count,
            });

            if let Some(parent) = Path::new(json_path).parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent)?;
                }
            }

            let mut file = File::create(json_path)?;
            writeln!(file, "{}", serde_json::to_string_pretty(&output)?)?;

            if verbose {
                eprintln!("Alternative analysis written to {}", json_path);
            }
        }
    }

    if let Some(path) = output_fasta {
        if let Some(parent) = Path::new(path).parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let mut fh = File::create(path)?;
        let header = format!(
            "assembled_from_{}_{}",
            Path::new(reads1_path)
                .file_name()
                .unwrap_or_else(|| "reads1".as_ref())
                .to_string_lossy(),
            Path::new(reads2_path)
                .file_name()
                .unwrap_or_else(|| "reads2".as_ref())
                .to_string_lossy()
        );
        writeln!(fh, ">{header}")?;
        writeln!(fh, "{assembled}")?;
    } else if verbose {
        println!(
            ">assembled_from_{}_{}\n{}",
            Path::new(reads1_path)
                .file_name()
                .unwrap_or_else(|| "reads1".as_ref())
                .to_string_lossy(),
            Path::new(reads2_path)
                .file_name()
                .unwrap_or_else(|| "reads2".as_ref())
                .to_string_lossy(),
            assembled
        );
    }

    if let Some(reference_seq) = reference {
        const MAX_DISTANCE_LEN: usize = 20_000;
        if assembled == reference_seq {
            println!(
                "Assembled contig matches the reference sequence exactly ({} bp).",
                assembled.len()
            );
        } else if assembled.len() <= MAX_DISTANCE_LEN && reference_seq.len() <= MAX_DISTANCE_LEN {
            let distance = levenshtein(assembled.as_bytes(), reference_seq.as_bytes());
            println!(
                "Edit distance to reference (len {} vs {}): {}",
                assembled.len(),
                reference_seq.len(),
                distance
            );
        } else {
            println!(
                "Reference check skipped: assembled length {} or reference length {} exceeds {} bp threshold.",
                assembled.len(),
                reference_seq.len(),
                MAX_DISTANCE_LEN
            );
        }
    }

    Ok(assembled)
}

#[cfg(test)]
mod smoke {
    use super::*;

    #[test]
    fn smoke_run() {
        use std::io::Write;

        let tmp1 = tempfile::NamedTempFile::new().expect("tmpfile");
        let tmp2 = tempfile::NamedTempFile::new().expect("tmpfile");
        writeln!(tmp1.as_file(), "ACGT").unwrap();
        writeln!(tmp2.as_file(), "GTAA").unwrap();

        let res = run_pipeline(
            tmp1.path().to_str().unwrap(),
            tmp2.path().to_str().unwrap(),
            None,
            None,
            false,
            None,
            None,
            false,
            false,
        );
        assert!(res.is_ok());
    }
}
