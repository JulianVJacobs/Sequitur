#!/usr/bin/env bash
set -euo pipefail

# Synthetic dataset generation using mason2 if available, else Python fallback.
# Generates tiered datasets from a subset of the Streptococcus pneumoniae reference.
# Output FASTQ read pairs and a metadata YAML manifest.
#
# Default tiers:
#   small  : 10k fragments / ~10k read pairs
#   medium : 50k fragments / ~50k read pairs
#   large  : 100k fragments / ~100k read pairs (ensure memory headroom)
#
# Usage:
#   scripts/synthetic/generate_synthetic_datasets.sh \
#       [--ref-length 200000] [--tiers small,medium] [--read-length 150] [--force]
#
# Reference used (subset): tests/fixtures/real/GCF_000329985.1_ASM32998v1_genomic.fna.gz
# Mason2 required binaries: mason_frag_sequencing, mason_simulator
# Fallback: random Python generator if mason2 is not found on PATH.
# Metadata: tests/fixtures/synthetic/synthetic_datasets.yaml

REF_GZ="tests/fixtures/real/GCF_000329985.1_ASM32998v1_genomic.fna.gz"
OUT_DIR="tests/fixtures/synthetic"
META_FILE="$OUT_DIR/synthetic_datasets.yaml"
REF_LEN=200000
TIERS="small,medium"
READ_LEN=150
FORCE=0

log(){ printf "[generate] %s\n" "$*"; }
err(){ printf "[generate][error] %s\n" "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case $1 in
    --ref-length) REF_LEN=$2; shift 2;;
    --tiers) TIERS=$2; shift 2;;
    --read-length) READ_LEN=$2; shift 2;;
    --force) FORCE=1; shift;;
    *) err "Unknown argument: $1";;
  esac
done

[[ -f "$REF_GZ" ]] || err "Reference gzip not found at $REF_GZ"
mkdir -p "$OUT_DIR"

# Build truncated reference FASTA (preserve header + wrap at 80 cols)
TRUNC_REF="$OUT_DIR/ref_subset_${REF_LEN}.fasta"
if [[ ! -f "$TRUNC_REF" || $FORCE -eq 1 ]]; then
  log "Creating truncated reference ($REF_LEN bp) -> $TRUNC_REF"
  gzip -dc "$REF_GZ" | awk -v len="$REF_LEN" 'NR==1{header=$0;next} /^[>]/{next} {seq=seq $0} END{print header; subseq=substr(seq,1,len); for(i=1;i<=length(subseq);i+=80) print substr(subseq,i,80)}' > "$TRUNC_REF"
fi

# Tier fragment counts
declare -A FRAG_COUNTS=( [small]=10000 [medium]=50000 [large]=100000 )

IFS=',' read -r -a ACTIVE_TIERS <<< "$TIERS"

MASON_SIM=$(command -v mason_simulator || true)
USE_FALLBACK=0
if [[ -z "$MASON_SIM" ]]; then
  log "Mason2 not found, using Python fallback generator."
  USE_FALLBACK=1
fi

# Initialize metadata file
echo "datasets:" > "$META_FILE"
echo "  reference_subset: $TRUNC_REF" >> "$META_FILE"
echo "  reference_length: $REF_LEN" >> "$META_FILE"
echo "  read_length: $READ_LEN" >> "$META_FILE"
echo "  generated_at: $(date -Iseconds)" >> "$META_FILE"
echo "  tiers:" >> "$META_FILE"

for tier in "${ACTIVE_TIERS[@]}"; do
  [[ -n "${FRAG_COUNTS[$tier]:-}" ]] || err "Unknown tier '$tier'"
  count="${FRAG_COUNTS[$tier]}"
  log "Generating tier '$tier' with $count read pairs"
  prefix="$OUT_DIR/${tier}" 
  if [[ $USE_FALLBACK -eq 0 ]]; then
    # Mason2 direct simulation from reference
    r1="${prefix}_reads_1.fq"; r2="${prefix}_reads_2.fq"
    if [[ $FORCE -eq 1 || ! -f "$r1" ]]; then
      log "  mason_simulator -> ${r1} / ${r2}"
      "$MASON_SIM" -ir "$TRUNC_REF" -n "$count" --seq-technology illumina \
        --fragment-mean-size 400 --fragment-size-std-dev 50 \
        --illumina-read-length "$READ_LEN" \
        -o "$r1" -or "$r2" -q >/dev/null 2>&1
    fi
    [[ -f "$r1" && -f "$r2" ]] || err "Expected mason output $r1 / $r2 missing"
  else
    # Fallback Python random generator
    r1="${prefix}_reads_1.fastq"; r2="${prefix}_reads_2.fastq"
    if [[ $FORCE -eq 1 || ! -f "$r1" ]]; then
      log "  Python fallback generating $r1 / $r2"
      python3 - "$TRUNC_REF" "$count" "$READ_LEN" "$r1" "$r2" <<'PYFALLBACK'
import sys, random
ref_path, count, read_len, r1_path, r2_path = sys.argv[1:6]
count=int(count); read_len=int(read_len)
with open(ref_path) as fh:
    fh.readline()  # header
    seq=''.join(line.strip() for line in fh if not line.startswith('>'))
ql='I'*read_len
def rand_read():
    start=random.randint(0, max(0,len(seq)-read_len))
    return seq[start:start+read_len]
def revcomp(s):
    comp={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    return ''.join(comp.get(b,'N') for b in reversed(s))
with open(r1_path,'w') as o1, open(r2_path,'w') as o2:
    for i in range(count):
        r1=rand_read(); r2=revcomp(rand_read())
        o1.write(f"@synthetic_{i}/1\n{r1}\n+\n{ql}\n")
        o2.write(f"@synthetic_{i}/2\n{r2}\n+\n{ql}\n")
PYFALLBACK
    fi
  fi
  echo "    - id: $tier" >> "$META_FILE"
  echo "      fragments: $count" >> "$META_FILE"
  echo "      read_pair_files:" >> "$META_FILE"
  echo "        - $r1" >> "$META_FILE"
  echo "        - $r2" >> "$META_FILE"
done

log "Synthetic dataset generation complete. Metadata: $META_FILE"
log "To benchmark: scripts/synthetic/benchmark_synthetic.sh (after creation)"