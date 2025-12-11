#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FIXTURES=("$ROOT/tests/fixtures"/*.fastq)
VERBOSE_FLAG=""
if [ "${VERBOSE:-}" = "1" ]; then
  VERBOSE_FLAG="--verbose"
fi
RUST_BIN="$ROOT/rust/target/debug/sequitur"
# Prefer project venv if present
if [ -x "$ROOT/python/.venv/bin/python" ]; then
  PYTHON_BIN="$ROOT/python/.venv/bin/python"
else
  PYTHON_BIN="$(command -v python3 || true)"
fi

for FIXTURE in "${FIXTURES[@]}"; do
  echo "Running integration test on fixture: $FIXTURE"

  # Build Rust binary if not present
  if [ ! -x "$RUST_BIN" ]; then
    echo "Building Rust binary..."
    (cd "$ROOT/rust" && cargo build)
  fi

  BASE=$(basename "$FIXTURE" .fastq)
  RUST_OUT="$ROOT/tests/results/${BASE}_rust.fasta"
  PY_OUT="$ROOT/tests/results/${BASE}_py.fasta"
  mkdir -p "$(dirname "$RUST_OUT")"
  "$RUST_BIN" "$FIXTURE" "$FIXTURE" --output-fasta "$RUST_OUT" $VERBOSE_FLAG || true

  if [ -n "$PYTHON_BIN" ] && "$PYTHON_BIN" -c "import sys; import Bio" 2>/dev/null; then
    echo "Running Python reference..."
    "$PYTHON_BIN" "$ROOT/python/sequitur.py" "$FIXTURE" "$FIXTURE" --output-fasta "$PY_OUT" || true
  else
    echo "Python environment with Biopython not available; skipping Python run."
  fi

  RUST_SEQ=$(grep -v '^>' "$RUST_OUT" | tr -d '\n' || true)
  PY_SEQ=$(grep -v '^>' "$PY_OUT" | tr -d '\n' || true)

  echo "Rust assembled: $RUST_SEQ"
  echo "Python assembled: $PY_SEQ"

  if [ -n "$PY_SEQ" ] && [ "$RUST_SEQ" != "$PY_SEQ" ]; then
    echo "Mismatch between Rust and Python assembled sequences for $FIXTURE"
    exit 2
  fi
done

echo "Integration completed" 
