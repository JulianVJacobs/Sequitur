use std::fs::OpenOptions;
use std::io::Write;
use std::path::PathBuf;

use bio::io::fastq;
use clap::Parser;
use serde::Serialize;

#[derive(Parser)]
struct Args {
    /// Input FASTQ/FASTA path
    #[arg(long, short)]
    input: PathBuf,

    /// Output base path (default: input with .seqs)
    #[arg(long)]
    output: Option<PathBuf>,
}

#[derive(Serialize)]
struct IndexRecord {
    name: String,
    offset: u64,
    len: u64,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let input = args.input;
    let out_base = args.output.unwrap_or_else(|| input.with_extension("seqs"));
    let seqs_tmp = out_base.with_extension("seqs.tmp");
    let sidx_tmp = out_base.with_extension("sidx.json.tmp");

    let mut seqs_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&seqs_tmp)?;
    let rdr = fastq::Reader::from_file(&input)?;

    let mut records = Vec::new();
    let mut offset: u64 = 0;
    for result in rdr.records() {
        let rec = result?;
        let seq = rec.seq();
        let name = rec.id().to_string();
        seqs_file.write_all(seq)?;
        records.push(IndexRecord {
            name,
            offset,
            len: seq.len() as u64,
        });
        offset += seq.len() as u64;
    }
    seqs_file.flush()?;
    seqs_file.sync_all()?;

    // write index JSON
    let mut sidx_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&sidx_tmp)?;
    serde_json::to_writer(&mut sidx_file, &records)?;
    sidx_file.flush()?;
    sidx_file.sync_all()?;

    // atomic rename
    let seqs_final = out_base.with_extension("seqs");
    let sidx_final = out_base.with_extension("sidx.json");
    std::fs::rename(&seqs_tmp, &seqs_final)?;
    std::fs::rename(&sidx_tmp, &sidx_final)?;

    println!(
        "Wrote {} records to {} and {}",
        records.len(),
        seqs_final.display(),
        sidx_final.display()
    );
    Ok(())
}
