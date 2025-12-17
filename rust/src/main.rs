use env_logger;
use log::{debug, info};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use bio::alignment::distance::levenshtein;
use bio::io::{fasta, fastq};
use clap::{Parser, Subcommand};
use flate2::read::MultiGzDecoder;

use sequitur::{
    create_overlap_graph_unified, create_overlap_graph_unified_from_readsource,
    find_first_subdiagonal_path, OverlapConfig,
};

/// Sequitur Rust CLI (prototype)
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Optional subcommand. If omitted, the CLI behaves like the previous flat interface (assemble).
    #[command(subcommand)]
    command: Option<Commands>,

    #[command(flatten)]
    assemble: AssembleArgs,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Build a read index from two read files and write the index base to `--output-index`
    Index(IndexArgs),
}

#[derive(Parser, Debug)]
struct IndexArgs {
    /// FASTQ/FASTA file for read set 1
    reads1: String,

    /// FASTQ/FASTA file for read set 2 (reverse-complemented during ingest)
    reads2: String,

    /// Output base path for the index (will create `<base>.seqs` and `<base>.sidx.json`)
    #[arg(long, required = true)]
    output_index: String,

    /// Skip reverse complement of reads2 when building the index
    #[arg(long)]
    no_revcomp: bool,
}

#[derive(Parser, Debug)]
struct AssembleArgs {
    /// Optional NDJSON file (.jsonl) to export read index-to-sequence mapping
    #[arg(long)]
    export_read_map: Option<String>,

    /// Enable threaded overlap graph construction (default: off)
    #[arg(long, default_value_t = false)]
    threads: bool,

    /// Number of worker threads for overlap graph construction (default: max available - 1)
    #[arg(long, default_value_t = num_cpus::get() - 1)]
    max_workers: usize,

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

    /// Optional output file for assembly graph (NDJSON .jsonl format with edge list and node attributes)
    #[arg(long)]
    export_graph_jsonl: Option<String>,

    /// Optional CSV file to append assembly metrics (reserved for future use)
    #[arg(long)]
    metrics_csv: Option<String>,

    /// Detect and remove low-quality overlaps by finding the cliff (knee) in score distribution (default: true)
    #[arg(long, default_value_t = true)]
    detect_score_cliff: bool,

    /// Output NDJSON file (.jsonl) for alternative path analysis
    #[arg(long)]
    alternatives_jsonl: Option<String>,

    /// Verbose/info output (default: quiet)
    #[arg(long, short = 'v', alias = "info")]
    verbose: bool,

    /// Debug output
    #[arg(long)]
    debug: bool,

    /// Trace output
    #[arg(long)]
    trace: bool,

    /// Wrap assembled FASTA lines to this width (0 = no-wrap)
    #[arg(long, default_value_t = 60)]
    fasta_line_width: usize,

    /// Skip reverse complement of reads2 (for non-genomic or unpaired data)
    #[arg(long)]
    no_revcomp: bool,

    /// Optimise swap squares for maximal diagonal sum (deterministic)
    #[arg(long)]
    optimise_diagonal: bool,

    /// Use the original array-based affix structure instead of the trie (trie is default)
    #[arg(long)]
    use_array: bool,

    /// Maximum normalised edit fraction for overlaps (e.g., 0.25). <0.1 disables fuzzy k-mer augmentation.
    #[arg(long, default_value_t = 0.25)]
    max_diff: f32,

    /// Minimum suffix/overlap length to consider (filters tiny overlaps)
    #[arg(long, default_value_t = sequitur::DEFAULT_MIN_SUFFIX_LEN)]
    min_suffix_len: usize,

    /// Optional path to a prebuilt read index (.seqs/.sidx.json) to use for read access.
    /// If the flag is provided with no value, a temporary index will be created and
    /// deleted after the run unless `--output-index` is also supplied (in which case
    /// the index will be persisted to that path).
    #[arg(long, num_args = 0..=1)]
    read_index: Option<Option<String>>,

    /// Exponent for error penalty in quality-adjusted scoring (1.0=linear, 2.0=quadratic)
    #[arg(long, default_value_t = 2.0)]
    error_penalty_exponent: f32,
    /// Optional output base path to persist an index when auto-building (`--output-index base`)
    #[arg(long)]
    output_index: Option<String>,
}

fn main() {
    let cli = Cli::parse();
    // If an explicit subcommand was supplied, handle it first.
    if let Some(Commands::Index(idx)) = &cli.command {
        // Call the library index builder and exit.
        let idx_base = std::path::Path::new(&idx.output_index);
        if let Err(e) = sequitur::read_source::build_index_from_pair(
            std::path::Path::new(&idx.reads1),
            std::path::Path::new(&idx.reads2),
            idx_base,
            idx.no_revcomp,
        ) {
            eprintln!("Failed to build read index: {:?}", e);
            std::process::exit(1);
        }
        println!("Wrote index base: {}", idx.output_index);
        std::process::exit(0);
    }

    // Default to assemble behavior when no subcommand is provided
    let args = cli.assemble;
    // Set log level based on CLI flags
    let log_level = if args.trace {
        "trace"
    } else if args.debug {
        "debug"
    } else if args.verbose {
        "info"
    } else {
        "error"
    };
    unsafe {
        std::env::set_var("RUST_LOG", log_level);
    }
    env_logger::init();

    info!("Sequitur Rust prototype");
    info!("reads1: {}", args.reads1);
    info!("reads2: {}", args.reads2);
    if let Some(refp) = &args.reference {
        info!("reference: {}", refp);
    }

    // Normalize the `--read-index` flag into `ReadIndexArg` (handled below by caller).
    let read_index_mode = match &args.read_index {
        None => ReadIndexArg::Absent,
        Some(None) => ReadIndexArg::Temp,
        Some(Some(s)) => ReadIndexArg::Base(s.clone()),
    };

    if let Err(error) = run_pipeline(
        &args.reads1,
        &args.reads2,
        args.output_fasta.as_deref(),
        args.reference.as_deref(),
        args.detect_score_cliff,
        args.alternatives_jsonl.as_deref(),
        args.verbose,
        args.fasta_line_width,
        args.no_revcomp,
        args.export_graph_jsonl.as_deref(),
        args.optimise_diagonal,
        args.threads,
        args.max_workers,
        args.use_array,
        args.max_diff,
        args.min_suffix_len,
        read_index_mode,
        args.output_index.as_deref(),
        args.export_read_map.as_deref(),
        args.error_penalty_exponent,
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

fn read_sequences(path: &Path) -> Result<(Vec<String>, Vec<String>)> {
    let format = infer_format(path);
    let reader = open_reader(path)?;

    match format {
        SequenceFormat::Fastq => {
            let fastq_reader = fastq::Reader::new(reader);
            let mut sequences = Vec::new();
            let mut read_ids = Vec::new();
            for record in fastq_reader.records() {
                let record = record.with_context(|| {
                    format!("Error reading FASTQ record from {}", path.display())
                })?;
                let seq = uppercase_sequence(record.seq())?;
                let id = record.id().to_string();
                sequences.push(seq);
                read_ids.push(id);
            }
            Ok((sequences, read_ids))
        }
        SequenceFormat::Fasta => {
            let fasta_reader = fasta::Reader::new(reader);
            let mut sequences = Vec::new();
            let mut read_ids = Vec::new();
            for record in fasta_reader.records() {
                let record = record.with_context(|| {
                    format!("Error reading FASTA record from {}", path.display())
                })?;
                let seq = uppercase_sequence(record.seq())?;
                let id = record.id().to_string();
                sequences.push(seq);
                read_ids.push(id);
            }
            Ok((sequences, read_ids))
        }
        SequenceFormat::Lines => {
            let mut sequences = Vec::new();
            let mut read_ids = Vec::new();
            let mut buf_reader = reader;
            let mut idx = 0;
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
                read_ids.push(format!("line_{}", idx));
                idx += 1;
            }
            Ok((sequences, read_ids))
        }
    }
}

fn load_reference(path: &Path) -> Result<Option<String>> {
    if !path.exists() {
        bail!("Reference path {} does not exist", path.display());
    }
    let (sequences, _read_ids) = read_sequences(path)?;
    if sequences.is_empty() {
        return Ok(None);
    }
    let mut combined = String::new();
    for seq in sequences {
        combined.push_str(&seq);
    }
    Ok(Some(combined))
}

#[derive(Debug, Clone)]
enum ReadIndexArg {
    Absent,
    Temp,
    Base(String),
}

fn run_pipeline(
    reads1_path: &str,
    reads2_path: &str,
    output_fasta: Option<&str>,
    reference_path: Option<&str>,
    detect_score_cliff: bool,
    alternatives_jsonl: Option<&str>,
    verbose: bool,
    fasta_line_width: usize,
    no_revcomp: bool,
    export_graph_jsonl: Option<&str>,
    optimise_diagonal: bool,
    use_threads: bool,
    max_workers: usize,
    use_array: bool,
    max_diff_cli: f32,
    min_suffix_len_cli: usize,
    read_index: ReadIndexArg,
    output_index: Option<&str>,
    export_read_map: Option<&str>,
    error_penalty_exponent: f32,
) -> Result<String> {
    // If an index is requested (explicit base, temp request, or output_index supplied)
    // avoid reading and reverse-complementing the reads into memory here to prevent
    // double-orientation effects; the index-backed branch will materialize reads
    // when needed. Otherwise load reads into memory as before.
    let index_requested = !matches!(read_index, ReadIndexArg::Absent) || output_index.is_some();

    let mut reads: Vec<String> = Vec::new();
    let mut read_ids: Vec<String> = Vec::new();
    let mut reads1_count: usize = 0;
    let mut reads2_count: usize = 0;

    if !index_requested {
        let (reads1, reads1_ids) = read_sequences(Path::new(reads1_path))
            .with_context(|| format!("Failed to parse reads from {}", reads1_path))?;
        let (reads2_raw, reads2_ids_raw) = read_sequences(Path::new(reads2_path))
            .with_context(|| format!("Failed to parse reads from {}", reads2_path))?;

        // In-memory path: do NOT reverse-complement reads2 here. The index builder
        // is responsible for reverse-complementing reads2 when requested, and the
        // pipeline will consume indexed materialization as-is. Keep in-memory reads
        // in their original orientation to avoid double-RC confusion.
        let reads2 = reads2_raw;

        // Keep original read IDs (no "/RC" annotation) for in-memory mode.
        let reads2_ids: Vec<String> = reads2_ids_raw;

        reads1_count = reads1.len();
        reads2_count = reads2.len();
        reads = Vec::with_capacity(reads1_count + reads2_count);
        read_ids = Vec::with_capacity(reads1_count + reads2_count);
        reads.extend(reads1.iter().cloned());
        reads.extend(reads2.iter().cloned());
        read_ids.extend(reads1_ids.iter().cloned());
        read_ids.extend(reads2_ids.iter().cloned());
    }

    // Export read map as NDJSON if requested
    if let Some(json_path) = export_read_map {
        use serde_json::json;
        use std::io::Write;
        if let Some(parent) = Path::new(json_path).parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let mut file = std::fs::File::create(json_path)?;
        // Write metadata as first line
        let metadata = json!({
            "type": "metadata",
            "reads1_path": reads1_path,
            "reads2_path": reads2_path,
            "reads1_count": reads1_count,
            "reads2_count": reads2_count,
            "reverse_complemented": !no_revcomp
        });
        writeln!(file, "{}", serde_json::to_string(&metadata)?)?;
        // Write each read as a separate line
        for idx in 0..reads.len() {
            let read_entry = if idx < reads1_count {
                json!({
                    "type": "read",
                    "matrix_index": idx,
                    "source": "reads1",
                    "source_index": idx
                })
            } else {
                json!({
                    "type": "read",
                    "matrix_index": idx,
                    "source": "reads2",
                    "source_index": idx - reads1_count
                })
            };
            writeln!(file, "{}", serde_json::to_string(&read_entry)?)?;
        }
        info!("Read map exported to {} (NDJSON format)", json_path);
    }

    let config = OverlapConfig {
        use_threads,
        max_workers,
        use_trie: !use_array,
        max_diff: max_diff_cli,
        min_suffix_len: min_suffix_len_cli,
        error_penalty_exponent,
        detect_score_cliff,
        ..OverlapConfig::default()
    };
    if config.use_trie {
        info!("Using pruned affix trie (unified path)");
    } else {
        info!("Using affix array (legacy path)");
    }
    if config.detect_score_cliff {
        info!("Score-cliff detection enabled: will remove low-quality overlaps using knee-point analysis");
    }
    info!("Creating overlap graph (unified)...");
    // Determine whether to use an on-disk index (either from `--read-index` or `--output-index`).
    let (mut adjacency_matrix, overlap_matrix) = {
        // `_temp_index_dir` will hold a TempDir when we auto-create a temporary index
        // so that the directory (and files) live for the duration of this scope.
        // Leading underscore suppresses unused-variable warnings while keeping it alive.
        let mut _temp_index_dir: Option<tempfile::TempDir> = None;
        // Choose a candidate base: prefer explicit `--read-index <base>`, otherwise
        // `--output-index`. If `--read-index` was provided with no value, a temporary
        // index will be created (unless `--output-index` is supplied, which causes
        // the index to be persisted to that path per user preference).
        let mut chosen_base: Option<std::path::PathBuf> = None;
        match &read_index {
            ReadIndexArg::Base(s) => chosen_base = Some(std::path::PathBuf::from(s.clone())),
            ReadIndexArg::Temp => {
                if let Some(o) = output_index {
                    chosen_base = Some(std::path::PathBuf::from(o));
                } else {
                    let td = tempfile::tempdir()
                        .map_err(|e| anyhow!("Failed to create tempdir: {:?}", e))?;
                    let tmpbase = td.path().join("sequitur_index");
                    info!("Building temporary index at {:?}", tmpbase);
                    _temp_index_dir = Some(td);
                    chosen_base = Some(tmpbase);
                }
            }
            ReadIndexArg::Absent => {
                if let Some(o) = output_index {
                    chosen_base = Some(std::path::PathBuf::from(o));
                }
            }
        }

        if let Some(base) = chosen_base {
            let seqs_path = base.with_extension("seqs");
            let sidx_path = seqs_path.with_extension("sidx.json");
            if !seqs_path.exists() || !sidx_path.exists() {
                // Build index at `base`. If the user requested persistence via
                // `--output-index` then `base` will reflect that path; if the user
                // requested a temporary index (via `--read-index` with no value),
                // `base` was already set above to a path under a tempdir.
                info!("Index files not found for {:?}, building index", base);
                sequitur::read_source::build_index_from_pair(
                    std::path::Path::new(reads1_path),
                    std::path::Path::new(reads2_path),
                    &base,
                    no_revcomp,
                )
                .map_err(|e| anyhow!("Failed to build read index: {:?}", e))?;
            }

            let base_str = base.to_string_lossy().to_string();
            match sequitur::read_source::BinaryIndexReadSource::open(&base_str) {
                Ok(idx_src) => {
                    let (adj, overlaps, mat_reads, mat_names) =
                        create_overlap_graph_unified_from_readsource(&idx_src, config);
                    reads = mat_reads;
                    read_ids = mat_names;
                    (adj, overlaps)
                }
                Err(e) => {
                    log::warn!(
                        "Failed to open read index {:?}, falling back to in-memory: {:?}",
                        base_str,
                        e
                    );
                    create_overlap_graph_unified(&reads, config)
                }
            }
        } else {
            create_overlap_graph_unified(&reads, config)
        }
    };

    info!("Overlap graph created.");

    let reference = if let Some(path) = reference_path {
        load_reference(Path::new(path))?
    } else {
        None
    };

    if optimise_diagonal {
        use sequitur::alternative_paths::{detect_swap_squares, optimise_first_subdiagonal_sum};
        let squares = detect_swap_squares(&adjacency_matrix, None);
        optimise_first_subdiagonal_sum(&mut adjacency_matrix, &squares);
        info!(
            "Applied maximal diagonal optimisation over {} swap squares.",
            squares.len()
        );
    }

    if let Some(graph_path) = export_graph_jsonl {
        use serde_json::json;
        use std::io::Write;
        if let Some(parent) = std::path::Path::new(graph_path).parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let mut nodes = Vec::new();
        for (idx, read) in reads.iter().enumerate() {
            nodes.push(json!({"id": idx, "sequence": read}));
        }
        let mut file = std::fs::File::create(graph_path)?;
        // Write metadata as first line
        let metadata = json!({"type": "metadata", "node_count": nodes.len()});
        writeln!(file, "{}", serde_json::to_string(&metadata)?)?;
        // Write each node as a separate line
        for node in nodes {
            writeln!(
                file,
                "{}",
                serde_json::to_string(&json!({"type": "node", "data": node}))?
            )?;
        }
        // Write each edge as a separate line
        for (row_idx, row_vec) in adjacency_matrix.outer_iterator().enumerate() {
            for (col_idx, &weight) in row_vec.indices().iter().zip(row_vec.data().iter()) {
                let edge =
                    json!({"type": "edge", "source": row_idx, "target": col_idx, "weight": weight});
                writeln!(file, "{}", serde_json::to_string(&edge)?)?;
            }
        }
        info!("Assembly graph written to {} (NDJSON format)", graph_path);
    }

    let adjacency_csc = adjacency_matrix.to_csc();
    let overlap_csc = overlap_matrix.to_csc();

    debug!(
        "[DEBUG] Created adjacency CSC with {} edges",
        adjacency_csc.nnz()
    );

    let (assembled, assembly_path) =
        find_first_subdiagonal_path(&adjacency_csc, &overlap_csc, &reads, &read_ids, None);
    debug!("[DEBUG] Assembly sequence length: {}", assembled.len());
    debug!("[DEBUG] Assembly path order: {:?}", assembly_path);

    // Export combined alternatives JSONL if requested (per-read successors with assembly order)
    if let Some(json_path) = alternatives_jsonl {
        use sequitur::extract_read_alternatives;
        use serde_json::json;

        let per_read = extract_read_alternatives(
            &adjacency_matrix,
            &overlap_matrix,
            &read_ids,
            Some(&assembly_path),
        );
        if let Some(parent) = Path::new(json_path).parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let mut file = File::create(json_path)?;

        // Write each per_read entry as a separate line (now includes assembly_order field)
        for entry in per_read {
            writeln!(
                file,
                "{}",
                serde_json::to_string(&json!({"type": "per_read", "data": entry}))?
            )?;
        }
        info!(
            "[ALT] Alternatives with assembly order written to {} (NDJSON format)",
            json_path
        );
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
        if fasta_line_width == 0 {
            writeln!(fh, "{assembled}")?;
        } else {
            let mut i = 0;
            let seq_len = assembled.len();
            while i < seq_len {
                let end = std::cmp::min(i + fasta_line_width, seq_len);
                writeln!(fh, "{}", &assembled[i..end])?;
                i = end;
            }
        }
    } else if verbose {
        let header = format!(
            ">assembled_from_{}_{}",
            Path::new(reads1_path)
                .file_name()
                .unwrap_or_else(|| "reads1".as_ref())
                .to_string_lossy(),
            Path::new(reads2_path)
                .file_name()
                .unwrap_or_else(|| "reads2".as_ref())
                .to_string_lossy()
        );
        info!("{header}");
        if fasta_line_width == 0 {
            debug!("{assembled}");
        } else {
            let mut i = 0;
            let seq_len = assembled.len();
            while i < seq_len {
                let end = std::cmp::min(i + fasta_line_width, seq_len);
                debug!("{}", &assembled[i..end]);
                i = end;
            }
        }
    }

    if let Some(reference_seq) = reference {
        const MAX_DISTANCE_LEN: usize = 20_000;
        if assembled == reference_seq {
            info!(
                "Assembled contig matches the reference sequence exactly ({} bp).",
                assembled.len()
            );
        } else if assembled.len() <= MAX_DISTANCE_LEN && reference_seq.len() <= MAX_DISTANCE_LEN {
            let distance = levenshtein(assembled.as_bytes(), reference_seq.as_bytes());
            info!(
                "Edit distance to reference (len {} vs {}): {}",
                assembled.len(),
                reference_seq.len(),
                distance
            );
        } else {
            info!(
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
            true, // detect_score_cliff
            None,
            false,
            60,
            false,
            None,
            false,
            false, // use_threads
            1,     // max_workers
            false, // use_array
            0.25,
            sequitur::DEFAULT_MIN_SUFFIX_LEN,
            ReadIndexArg::Absent,
            None, // output_index
            None, // export_read_map
            2.0,  // error_penalty_exponent
        );
        assert!(res.is_ok());
    }
}
