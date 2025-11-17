#!/usr/bin/env bash
set -euo pipefail

# Benchmark Sequitur assembly on synthetic datasets.
# Runs Rust CLI against tiered datasets, captures timing and memory metrics.
#
# Usage:
#   scripts/synthetic/benchmark_synthetic.sh [--rust-binary PATH] [--tiers small,medium] [--output-dir DIR]
#
# Outputs:
#   - Assembled FASTA files in output directory
#   - Benchmark CSV with timing and memory stats
#   - Optional alternatives JSON if --analyse-alternatives passed
#
# Requirements:
#   - Rust sequitur_rs binary built (rust/target/release/sequitur_rs)
#   - Synthetic datasets generated (tests/fixtures/synthetic/)
#   - /usr/bin/time (GNU time for memory stats) or fallback to shell time

RUST_BIN="rust/target/release/sequitur_rs"
OUT_DIR="tests/results/synthetic"
TIERS="small,medium"
ANALYSE_ALT=0
METRICS_CSV="$OUT_DIR/benchmark_metrics.csv"

log(){ printf "[benchmark] %s\n" "$*"; }
err(){ printf "[benchmark][error] %s\n" "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case $1 in
    --rust-binary) RUST_BIN=$2; shift 2;;
    --tiers) TIERS=$2; shift 2;;
    --output-dir) OUT_DIR=$2; shift 2;;
    --analyse-alternatives) ANALYSE_ALT=1; shift;;
    *) err "Unknown argument: $1";;
  esac
done

[[ -f "$RUST_BIN" ]] || err "Rust binary not found: $RUST_BIN (run: cd rust && cargo build --release)"
[[ -f "tests/fixtures/synthetic/synthetic_datasets.yaml" ]] || err "Synthetic datasets not generated. Run: scripts/synthetic/generate_synthetic_datasets.sh"

mkdir -p "$OUT_DIR"

# Check if GNU time available (for memory tracking)
if command -v /usr/bin/time >/dev/null 2>&1; then
  TIME_CMD="/usr/bin/time -f 'TIME:%E MEM:%M'"
  HAS_GNU_TIME=1
else
  TIME_CMD="time"
  HAS_GNU_TIME=0
  log "GNU time not found, memory stats unavailable (install: apk add time)"
fi

# CSV header
echo "tier,read_count,assembly_time_sec,peak_memory_kb,assembled_length,output_file" > "$METRICS_CSV"

IFS=',' read -r -a ACTIVE_TIERS <<< "$TIERS"

for tier in "${ACTIVE_TIERS[@]}"; do
  r1="tests/fixtures/synthetic/${tier}_reads_1.fq"
  r2="tests/fixtures/synthetic/${tier}_reads_2.fq"
  [[ -f "$r1" && -f "$r2" ]] || { log "Skipping missing tier '$tier'"; continue; }
  
  # Count reads (lines / 4)
  read_count=$(($(wc -l < "$r1") / 4))
  
  log "Benchmarking tier '$tier' with $read_count read pairs"
  
  out_fasta="$OUT_DIR/${tier}_assembled.fasta"
  alt_json="$OUT_DIR/${tier}_alternatives.json"
  
  # Build command
  cmd=("$RUST_BIN" "$r1" "$r2" --output-fasta "$out_fasta")
  [[ $ANALYSE_ALT -eq 1 ]] && cmd+=(--analyse-alternatives --alternatives-json "$alt_json")
  
  # Run with timing
  start=$(date +%s.%N)
  if [[ $HAS_GNU_TIME -eq 1 ]]; then
    timing_out=$(mktemp)
    /usr/bin/time -f 'TIME:%E MEM:%M' -o "$timing_out" "${cmd[@]}" 2>&1 | tee "$OUT_DIR/${tier}_log.txt" || true
    timing=$(cat "$timing_out")
    rm -f "$timing_out"
    # Parse TIME:MM:SS.ms MEM:kb (BusyBox-compatible)
    time_str=$(echo "$timing" | sed -n 's/.*TIME:\([0-9:.]*\).*/\1/p')
    mem_kb=$(echo "$timing" | sed -n 's/.*MEM:\([0-9]*\).*/\1/p')
    [[ -z "$time_str" ]] && time_str="0:00.00"
    [[ -z "$mem_kb" ]] && mem_kb="0"
    # Convert MM:SS.ms to seconds
    if [[ "$time_str" =~ ([0-9]+):([0-9.]+) ]]; then
      time_sec=$(awk "BEGIN {print ${BASH_REMATCH[1]} * 60 + ${BASH_REMATCH[2]}}")
    else
      time_sec="$time_str"
    fi
  else
    # Fallback: shell time (no memory)
    { time "${cmd[@]}" 2>&1; } 2>&1 | tee "$OUT_DIR/${tier}_log.txt"
    end=$(date +%s.%N)
    time_sec=$(awk "BEGIN {print $end - $start}")
    mem_kb="N/A"
  fi
  
  # Get assembled length
  if [[ -f "$out_fasta" ]]; then
    asm_len=$(grep -v '^>' "$out_fasta" | tr -d '\n' | wc -c)
  else
    asm_len="FAILED"
  fi
  
  log "  Time: ${time_sec}s, Memory: ${mem_kb} KB, Length: $asm_len bp"
  echo "$tier,$read_count,$time_sec,$mem_kb,$asm_len,$out_fasta" >> "$METRICS_CSV"
done

log "Benchmark complete. Results: $METRICS_CSV"
log "Assembled FASTAs in: $OUT_DIR"
cat "$METRICS_CSV"
