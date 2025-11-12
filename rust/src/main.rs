use clap::Parser;

/// Sequitur Rust CLI (prototype)
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// FASTQ file for read set 1
    reads1: String,

    /// FASTQ file for read set 2
    reads2: String,

    /// Optional reference FASTA
    #[arg(long)]
    reference: Option<String>,
    /// Optional output FASTA path for the assembled sequence
    #[arg(long)]
    output_fasta: Option<String>,

    /// Optional CSV file to append assembly metrics
    #[arg(long)]
    metrics_csv: Option<String>,
    /// Verbose output (default: quiet)
    #[arg(long, short = 'v')]
    verbose: bool,
}

fn main() {
    env_logger::init();
    let args = Args::parse();

    if args.verbose {
        println!("Sequitur Rust prototype");
        println!("reads1: {}", args.reads1);
        println!("reads2: {}", args.reads2);
        if let Some(refp) = args.reference.clone() {
            println!("reference: {}", refp);
        }
    }

    // Run the pipeline for reads1 and reads2 sequentially (smoke-test behaviour).
    if let Err(e) = run_pipeline(&args.reads1, args.output_fasta.as_deref(), args.verbose) {
        eprintln!("Error running pipeline on {}: {:#}", args.reads1, e);
        std::process::exit(1);
    }

    if let Err(e) = run_pipeline(&args.reads2, args.output_fasta.as_deref(), args.verbose) {
        eprintln!("Error running pipeline on {}: {:#}", args.reads2, e);
        std::process::exit(1);
    }
}

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use sequitur_rs::{
    AffixArray,
    create_overlap_graph,
    adjacency_to_csc,
    find_lower_diagonal_path,
    OverlapConfig,
};

/// Robust FASTA/FASTQ reader that handles multiple records.
fn read_records<P: AsRef<Path>>(path: P) -> io::Result<Vec<String>> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    // Peek first non-empty byte to decide format
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    if first_line.trim().is_empty() {
        // try to find a non-empty
        loop {
            first_line.clear();
            if reader.read_line(&mut first_line)? == 0 {
                return Ok(Vec::new());
            }
            if !first_line.trim().is_empty() {
                break;
            }
        }
    }

    let mut reads: Vec<String> = Vec::new();
    if first_line.starts_with('>') {
        // FASTA: accumulate multi-line sequences per header
        let mut seq = String::new();
        for line in reader.lines() {
            let l = line?;
            if l.starts_with('>') {
                if !seq.is_empty() {
                    reads.push(seq.clone());
                    seq.clear();
                }
                continue;
            }
            let s = l.trim();
            if s.is_empty() {
                continue;
            }
            seq.push_str(s);
        }
        if !seq.is_empty() {
            reads.push(seq);
        }
    } else if first_line.starts_with('@') {
        // FASTQ: read groups of 4 lines
        // we already consumed the first header line; read sequence, plus, quality
        let mut seq = String::new();
        if reader.read_line(&mut seq)? > 0 {
            reads.push(seq.trim().to_string());
        }
        // consume remaining records
        loop {
            let mut header = String::new();
            if reader.read_line(&mut header)? == 0 {
                break;
            }
            if header.trim().is_empty() {
                continue;
            }
            if !header.starts_with('@') {
                // not a header; stop
                break;
            }
            let mut sequence = String::new();
            let mut plus = String::new();
            let mut qual = String::new();
            if reader.read_line(&mut sequence)? == 0 {
                break;
            }
            if reader.read_line(&mut plus)? == 0 {
                break;
            }
            if reader.read_line(&mut qual)? == 0 {
                break;
            }
            reads.push(sequence.trim().to_string());
        }
    } else {
        // one-sequence-per-line style
        if !first_line.trim().is_empty() {
            reads.push(first_line.trim().to_string());
        }
        for line in reader.lines() {
            let l = line?;
            if l.trim().is_empty() {
                continue;
            }
            reads.push(l.trim().to_string());
        }
    }

    Ok(reads)
}

fn run_pipeline(reads_path: &str, output_fasta: Option<&str>, verbose: bool) -> anyhow::Result<String> {
    let reads = read_records(reads_path)?;
    if verbose {
        eprintln!("Read {} sequences from {}", reads.len(), reads_path);
    }

    let affix = AffixArray::build(reads.iter().map(|s| s.as_str()), 3);
    if verbose {
        eprintln!("Built affix array with {} entries", affix.len());
    }

    let config = OverlapConfig::default();
    let (_affix_array, adjacency, overlaps) = create_overlap_graph(&reads, Some(affix), config);
    if verbose {
        eprintln!("Created adjacency for {} reads", adjacency.len());
    }

    let csc = adjacency_to_csc(&adjacency, None);
    let assembled = find_lower_diagonal_path(&csc, &overlaps, &reads, None);

    if let Some(path) = output_fasta {
        std::fs::create_dir_all(Path::new(path).parent().unwrap_or(Path::new(".")))?;
        let mut fh = File::create(path)?;
        writeln!(fh, ">assembled_from_{}", Path::new(reads_path).file_name().unwrap_or_default().to_string_lossy())?;
        writeln!(fh, "{}", assembled)?;
    } else if verbose {
        println!(">assembled_from_{}\n{}", Path::new(reads_path).file_name().unwrap_or_default().to_string_lossy(), assembled);
    }

    Ok(assembled)
}

#[cfg(test)]
mod smoke {
    use super::*;

    #[test]
    fn smoke_run() {
        // small in-memory smoke test: create two short reads and run pipeline via temporary file
        use std::io::Write;
        let tmp = tempfile::NamedTempFile::new().expect("tmpfile");
        writeln!(tmp.as_file(), "ACGT").unwrap();
        writeln!(tmp.as_file(), "GTAA").unwrap();
        let path = tmp.path().to_string_lossy().to_string();
        let res = run_pipeline(&path, None, false);
        assert!(res.is_ok());
    }
}
