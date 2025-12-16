```markdown
# Read Index (.seqs + .sidx.json)

This document describes the lightweight on-disk read index used by the Rust `sequitur` CLI.

Purpose
-------
- Allow large read datasets to remain on-disk and be referenced via random access (mmap/pread) instead of loading all reads into memory.
- Provide a compact, mmap-friendly layout for fast access during overlap/assembly passes.

Index files
-----------
- `<base>.seqs` — concatenated raw sequence bytes for all reads (no FASTQ headers or qualities). Designed to be mmapped for fast random access.
- `<base>.sidx.json` — JSON array of index entries with fields `{ "name": string, "offset": integer, "len": integer }` describing the read name and the byte range inside the `.seqs` file.

Indexer
-------
An indexer binary is provided under the Rust workspace as `index_reads`.

Example (build + index):

```bash
cd rust
cargo build --release
# concatenate reads1 then reads2 (non-RC) and make an index
cat ../tests/synthetic/swap_test/data/reads_1.fastq ../tests/synthetic/swap_test/data/reads_2.fastq > /tmp/swap_combined.fastq
./target/release/index_reads --input /tmp/swap_combined.fastq --output /tmp/swap_index
```

CLI usage
---------
Pass the index base (path without extension) to `sequitur` with `--read-index`:

```bash
./target/release/sequitur reads_1.fastq reads_2.fastq --read-index /tmp/swap_index \
    --output-fasta results/assembled.fasta --alternatives-jsonl results/alternatives.jsonl
```

Orientation note
----------------
By default `sequitur` applies a reverse-complement to the second read file (`reads2`) during loading. To ensure parity between an indexed run and the in-memory run, one of the following must hold:

- Build the index from the same input orientation that `sequitur` will use (recommended: concatenate `reads_1.fastq` then `reads_2.fastq` without pre-applying RC), or
- If your index was created from pre-RC (or pre-RC reads2), reindex using the non-RC files, or be aware the code will apply RC to the reads2 partition at runtime if needed.

Current limitations
-------------------
- The current prototype materialises read sequences during trie construction to reuse existing code paths; future work will allow streaming trie construction that avoids full materialisation.
- The `BinaryIndexReadSource` implementation currently expects an mmapped `.seqs` file; a non-mmap `pread` fallback and compressed/BGZF-aware formats are planned.
- The `.sidx.json` is intentionally simple (JSON) for prototyping. A compact binary index format may be added later.

Implementation notes
--------------------
- The ReadSource trait abstracts read access (num_reads, get_name, get_seq, get_subseq). The Rust implementation includes `InMemoryReadSource` and `BinaryIndexReadSource` (mmap + LRU cache).
- When using `--read-index`, `sequitur` will attempt to open `<base>.seqs` and `<base>.sidx.json` and will fall back to the in-memory path on failure.

If you add or change the index format, update this document and the indexer (`rust/src/bin/index_reads.rs`).

```
