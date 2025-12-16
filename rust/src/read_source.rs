use std::borrow::Cow;
use std::fs::File;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use bio::alphabets::dna;
use bio::io::{fasta, fastq};
use flate2::read::MultiGzDecoder;
use memmap2::Mmap;
use std::io::BufReader;

/// Errors returned by ReadSource implementations.
#[derive(thiserror::Error, Debug)]
pub enum ReadSourceError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Not found: {0}")]
    NotFound(String),
    #[error("Invalid range")]
    InvalidRange,
    #[error("Corrupt index or unsupported index format")]
    CorruptIndex,
    #[error("Unsupported compressed input for random access")]
    UnsupportedCompressed,
}

/// Trait that abstracts random-access read retrieval.
/// Implementations should be `Send + Sync` so they can be shared across threads.
pub trait ReadSource: Send + Sync {
    fn num_reads(&self) -> Result<usize, ReadSourceError>;
    fn get_name(&self, id: usize) -> Result<Option<String>, ReadSourceError>;
    fn resolve_name(&self, name: &str) -> Result<Option<usize>, ReadSourceError>;
    fn get_len(&self, id: usize) -> Result<usize, ReadSourceError>;
    fn get_seq(&self, id: usize) -> Result<Cow<'_, str>, ReadSourceError>;
    fn get_subseq(
        &self,
        id: usize,
        start: usize,
        len: usize,
    ) -> Result<Cow<'_, str>, ReadSourceError>;
}

/// A simple in-memory adapter over Vec<String> for compatibility.
pub struct InMemoryReadSource {
    reads: Arc<Vec<String>>,
    names: Option<Arc<Vec<String>>>,
}

impl InMemoryReadSource {
    pub fn new(reads: Vec<String>, names: Option<Vec<String>>) -> Self {
        InMemoryReadSource {
            reads: Arc::new(reads),
            names: names.map(Arc::new),
        }
    }
}

impl ReadSource for InMemoryReadSource {
    fn num_reads(&self) -> Result<usize, ReadSourceError> {
        Ok(self.reads.len())
    }

    fn get_name(&self, id: usize) -> Result<Option<String>, ReadSourceError> {
        if let Some(names) = &self.names {
            Ok(names.get(id).cloned())
        } else {
            Ok(None)
        }
    }

    fn resolve_name(&self, _name: &str) -> Result<Option<usize>, ReadSourceError> {
        // Linear scan fallback; callers should avoid using this on large datasets.
        if let Some(names) = &self.names {
            for (i, n) in names.iter().enumerate() {
                if n == _name {
                    return Ok(Some(i));
                }
            }
        }
        Ok(None)
    }

    fn get_len(&self, id: usize) -> Result<usize, ReadSourceError> {
        self.reads
            .get(id)
            .map(|s| s.len())
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))
    }

    fn get_seq(&self, id: usize) -> Result<Cow<'_, str>, ReadSourceError> {
        self.reads
            .get(id)
            .map(|s| Cow::Borrowed(s.as_str()))
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))
    }

    fn get_subseq(
        &self,
        id: usize,
        start: usize,
        len: usize,
    ) -> Result<Cow<'_, str>, ReadSourceError> {
        let s = self
            .reads
            .get(id)
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))?;
        if start.checked_add(len).map_or(false, |end| end <= s.len()) {
            Ok(Cow::Owned(s[start..start + len].to_string()))
        } else {
            Err(ReadSourceError::InvalidRange)
        }
    }
}

/// Index entry for binary index. Kept minimal for now.
#[derive(Debug, Clone)]
pub struct IndexEntry {
    pub data_offset: u64,
    pub data_len: u32,
    pub name_idx: u32,
    pub flags: u32,
}

/// A minimal BinaryIndexReadSource skeleton: attempts to mmap a `.seqs` file and read a simple JSON index.
#[allow(dead_code)]
pub struct BinaryIndexReadSource {
    seqs_path: PathBuf,
    sidx_path: PathBuf,
    mmap: Option<Mmap>,
    entries: Vec<IndexEntry>,
    names: Vec<String>,
    // Simple LRU cache of recently accessed read sequences (owned Strings inside Arc)
    cache: Mutex<lru::LruCache<usize, Arc<String>>>,
}

impl BinaryIndexReadSource {
    /// Open a pair of files: data file (".seqs") and index file (".sidx.json").
    /// For this initial implementation the index file is expected to be JSON written by the indexer.
    pub fn open<P: AsRef<Path>>(base: P) -> Result<Self, ReadSourceError> {
        let base = base.as_ref();
        let seqs_path = if base.extension().is_some() && base.extension().unwrap() == "seqs" {
            base.to_path_buf()
        } else {
            base.with_extension("seqs")
        };
        let sidx_path = seqs_path.with_extension("sidx.json");

        if !seqs_path.exists() || !sidx_path.exists() {
            return Err(ReadSourceError::NotFound(format!(
                "index files for {:?}",
                base
            )));
        }

        let file = File::open(&seqs_path)?;
        let mmap = unsafe { Mmap::map(&file).ok() };

        // Read JSON index (simple array of {name,offset,len}) for initial implementation.
        let idx_file = File::open(&sidx_path)?;
        let sidx: serde_json::Value =
            serde_json::from_reader(idx_file).map_err(|_| ReadSourceError::CorruptIndex)?;
        let mut entries = Vec::new();
        let mut names = Vec::new();
        if let Some(arr) = sidx.as_array() {
            for v in arr.iter() {
                let name = v
                    .get("name")
                    .and_then(|x| x.as_str())
                    .unwrap_or("")
                    .to_string();
                let offset = v.get("offset").and_then(|x| x.as_u64()).unwrap_or(0);
                let len = v.get("len").and_then(|x| x.as_u64()).unwrap_or(0) as u32;
                let entry = IndexEntry {
                    data_offset: offset,
                    data_len: len,
                    name_idx: names.len() as u32,
                    flags: 0,
                };
                names.push(name);
                entries.push(entry);
            }
        } else {
            return Err(ReadSourceError::CorruptIndex);
        }

        // Default cache capacity: 1024 entries
        let cache = Mutex::new(lru::LruCache::new(NonZeroUsize::new(1024).unwrap()));

        Ok(BinaryIndexReadSource {
            seqs_path,
            sidx_path,
            mmap,
            entries,
            names,
            cache,
        })
    }
}

/// Build an index from two sequence files (reads1 followed by reads2).
/// If `no_revcomp` is false, reads2 sequences will be reverse-complemented
/// before being written into the `.seqs` file and index.
pub fn build_index_from_pair<P: AsRef<Path>>(
    reads1: P,
    reads2: P,
    out_base: P,
    no_revcomp: bool,
) -> Result<(), ReadSourceError> {
    use std::fs::OpenOptions;
    use std::io::Write;

    let reads1 = reads1.as_ref();
    let reads2 = reads2.as_ref();
    let out_base = out_base.as_ref();
    let seqs_tmp = out_base.with_extension("seqs.tmp");
    let sidx_tmp = out_base.with_extension("sidx.json.tmp");

    let mut seqs_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&seqs_tmp)
        .map_err(ReadSourceError::Io)?;

    let mut records: Vec<serde_json::Value> = Vec::new();
    let mut offset: u64 = 0;

    // helper to open potentially gzipped files
    let open_buf = |p: &Path| -> Result<BufReader<Box<dyn std::io::Read>>, ReadSourceError> {
        let f = std::fs::File::open(p).map_err(ReadSourceError::Io)?;
        if p.extension()
            .and_then(|e| e.to_str())
            .map(|s| s.eq_ignore_ascii_case("gz"))
            .unwrap_or(false)
        {
            let dec = MultiGzDecoder::new(f);
            Ok(BufReader::new(Box::new(dec)))
        } else {
            Ok(BufReader::new(Box::new(f)))
        }
    };

    let mut write_seq_record = |name: String, seq_bytes: &[u8]| -> Result<(), ReadSourceError> {
        seqs_file
            .write_all(seq_bytes)
            .map_err(ReadSourceError::Io)?;
        records.push(serde_json::json!({"name": name, "offset": offset, "len": seq_bytes.len() }));
        offset += seq_bytes.len() as u64;
        Ok(())
    };

    // Process a file that may be FASTQ or FASTA
    let process_path =
        |p: &Path,
         is_reads2: bool,
         write_seq_record: &mut dyn FnMut(String, &[u8]) -> Result<(), ReadSourceError>|
         -> Result<(), ReadSourceError> {
            // determine format by extension (simple heuristic)
            let ext = p
                .extension()
                .and_then(|e| e.to_str())
                .map(|s| s.to_ascii_lowercase())
                .unwrap_or_default();
            let reader = open_buf(p)?;
            if ext == "fastq"
                || ext == "fq"
                || (ext == "gz"
                    && p.file_stem()
                        .and_then(|s| Path::new(s).extension())
                        .and_then(|e| e.to_str())
                        .map(|s| s.eq_ignore_ascii_case("fastq") || s.eq_ignore_ascii_case("fq"))
                        .unwrap_or(false))
            {
                let fq = fastq::Reader::new(reader);
                for res in fq.records() {
                    let rec = res.map_err(|_| {
                        ReadSourceError::Io(std::io::Error::new(
                            std::io::ErrorKind::Other,
                            "fastq read error",
                        ))
                    })?;
                    let mut seq = rec.seq().to_vec();
                    if is_reads2 && !no_revcomp {
                        seq = dna::revcomp(&seq);
                    }
                    let name = rec.id().to_string();
                    write_seq_record(name, &seq)?;
                }
                Ok(())
            } else {
                // treat as FASTA or plain lines using fasta reader (will handle fasta/gz)
                let fa = fasta::Reader::new(reader);
                for res in fa.records() {
                    let rec = res.map_err(|_| {
                        ReadSourceError::Io(std::io::Error::new(
                            std::io::ErrorKind::Other,
                            "fasta read error",
                        ))
                    })?;
                    let mut seq = rec.seq().to_vec();
                    if is_reads2 && !no_revcomp {
                        seq = dna::revcomp(&seq);
                    }
                    let name = rec.id().to_string();
                    write_seq_record(name, &seq)?;
                }
                Ok(())
            }
        };

    // Process reads1 then reads2
    process_path(reads1, false, &mut write_seq_record)?;
    process_path(reads2, true, &mut write_seq_record)?;

    seqs_file.flush().map_err(ReadSourceError::Io)?;
    seqs_file.sync_all().map_err(ReadSourceError::Io)?;

    // write JSON index
    let mut sidx_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&sidx_tmp)
        .map_err(ReadSourceError::Io)?;
    serde_json::to_writer(&mut sidx_file, &records).map_err(|_| {
        ReadSourceError::Io(std::io::Error::new(
            std::io::ErrorKind::Other,
            "failed to write index json",
        ))
    })?;
    sidx_file.flush().map_err(ReadSourceError::Io)?;
    sidx_file.sync_all().map_err(ReadSourceError::Io)?;

    // atomic rename
    let seqs_final = out_base.with_extension("seqs");
    let sidx_final = out_base.with_extension("sidx.json");
    std::fs::rename(&seqs_tmp, &seqs_final).map_err(ReadSourceError::Io)?;
    std::fs::rename(&sidx_tmp, &sidx_final).map_err(ReadSourceError::Io)?;

    Ok(())
}

impl BinaryIndexReadSource {
    /// Prefetch a batch of read ids into the LRU cache. Best-effort; errors are returned
    /// if underlying IO fails. This is intended to amortize random IO during verification.
    pub fn prefetch(&self, ids: &[usize]) -> Result<(), ReadSourceError> {
        // If no mmap, nothing to prefetch in this minimal implementation
        let m = match &self.mmap {
            Some(m) => m,
            None => return Ok(()),
        };

        let mut cache = self.cache.lock().map_err(|_| {
            ReadSourceError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                "cache lock poisoned",
            ))
        })?;

        for &id in ids {
            if cache.get(&id).is_some() {
                continue;
            }
            let e = self
                .entries
                .get(id)
                .ok_or(ReadSourceError::NotFound(format!("id {}", id)))?;
            let start = e.data_offset as usize;
            let end = start + e.data_len as usize;
            if end <= m.len() {
                let slice = &m[start..end];
                let s = std::str::from_utf8(slice).map_err(|_| ReadSourceError::CorruptIndex)?;
                cache.put(id, Arc::new(s.to_string()));
            } else {
                return Err(ReadSourceError::CorruptIndex);
            }
        }

        Ok(())
    }
}

impl ReadSource for BinaryIndexReadSource {
    fn num_reads(&self) -> Result<usize, ReadSourceError> {
        Ok(self.entries.len())
    }

    fn get_name(&self, id: usize) -> Result<Option<String>, ReadSourceError> {
        Ok(self.names.get(id).cloned())
    }

    fn resolve_name(&self, name: &str) -> Result<Option<usize>, ReadSourceError> {
        for (i, n) in self.names.iter().enumerate() {
            if n == name {
                return Ok(Some(i));
            }
        }
        Ok(None)
    }

    fn get_len(&self, id: usize) -> Result<usize, ReadSourceError> {
        self.entries
            .get(id)
            .map(|e| e.data_len as usize)
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))
    }

    fn get_seq(&self, id: usize) -> Result<Cow<'_, str>, ReadSourceError> {
        // First check cache
        if let Ok(mut cache) = self.cache.lock() {
            if let Some(arc) = cache.get(&id) {
                return Ok(Cow::Owned(arc.as_ref().clone()));
            }
        }

        let e = self
            .entries
            .get(id)
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))?;
        if let Some(m) = &self.mmap {
            let start = e.data_offset as usize;
            let end = start + e.data_len as usize;
            if end <= m.len() {
                // Safety: assume sequences are valid UTF-8 ASCII ACGT.
                let slice = &m[start..end];
                let s = std::str::from_utf8(slice).map_err(|_| ReadSourceError::CorruptIndex)?;
                // Optionally populate cache with owned copy
                if let Ok(mut cache) = self.cache.lock() {
                    cache.put(id, Arc::new(s.to_string()));
                }
                return Ok(Cow::Owned(s.to_string()));
            } else {
                return Err(ReadSourceError::CorruptIndex);
            }
        }
        // Fallback: not implemented for non-mmap in this minimal implementation.
        Err(ReadSourceError::Io(std::io::Error::new(
            std::io::ErrorKind::Other,
            "non-mmap fallback not implemented",
        )))
    }

    fn get_subseq(
        &self,
        id: usize,
        start: usize,
        len: usize,
    ) -> Result<Cow<'_, str>, ReadSourceError> {
        // Try cache first
        if let Ok(mut cache) = self.cache.lock() {
            if let Some(arc) = cache.get(&id) {
                let s = arc.as_str();
                if start.checked_add(len).map_or(false, |end| end <= s.len()) {
                    return Ok(Cow::Owned(s[start..start + len].to_string()));
                } else {
                    return Err(ReadSourceError::InvalidRange);
                }
            }
        }

        // Fall back to reading full sequence and slicing
        let full = self.get_seq(id)?;
        if start
            .checked_add(len)
            .map_or(false, |end| end <= full.len())
        {
            Ok(Cow::Owned(full[start..start + len].to_string()))
        } else {
            Err(ReadSourceError::InvalidRange)
        }
    }
}

/// Small helper to construct an in-memory adapter from Vec<String> easily.
impl From<Vec<String>> for InMemoryReadSource {
    fn from(v: Vec<String>) -> Self {
        InMemoryReadSource::new(v, None)
    }
}

// Implement ReadSource for string slices so existing code that passes
// `&[String]` can be used directly with ReadSource-based APIs without
// requiring construction of an adapter.
impl ReadSource for [String] {
    fn num_reads(&self) -> Result<usize, ReadSourceError> {
        Ok(self.len())
    }

    fn get_name(&self, id: usize) -> Result<Option<String>, ReadSourceError> {
        Ok(self.get(id).map(|s| s.clone()))
    }

    fn resolve_name(&self, name: &str) -> Result<Option<usize>, ReadSourceError> {
        for (i, n) in self.iter().enumerate() {
            if n == name {
                return Ok(Some(i));
            }
        }
        Ok(None)
    }

    fn get_len(&self, id: usize) -> Result<usize, ReadSourceError> {
        self.get(id)
            .map(|s| s.len())
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))
    }

    fn get_seq(&self, id: usize) -> Result<Cow<'_, str>, ReadSourceError> {
        self.get(id)
            .map(|s| Cow::Borrowed(s.as_str()))
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))
    }

    fn get_subseq(
        &self,
        id: usize,
        start: usize,
        len: usize,
    ) -> Result<Cow<'_, str>, ReadSourceError> {
        let s = self
            .get(id)
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))?;
        if start.checked_add(len).map_or(false, |end| end <= s.len()) {
            Ok(Cow::Owned(s[start..start + len].to_string()))
        } else {
            Err(ReadSourceError::InvalidRange)
        }
    }
}
