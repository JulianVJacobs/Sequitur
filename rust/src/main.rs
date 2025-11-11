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
}

fn main() {
    env_logger::init();
    let args = Args::parse();

    println!("Sequitur Rust prototype");
    println!("reads1: {}", args.reads1);
    println!("reads2: {}", args.reads2);
    if let Some(refp) = args.reference {
        println!("reference: {}", refp);
    }

    // TODO: Wire up to library functions in src/lib.rs
}
