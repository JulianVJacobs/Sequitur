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
    #[error("Quality validation error: {0}")]
    QualityValidation(String),
}

/// Quality score status: Present, Missing, or Invalid
#[derive(Debug, Clone)]
pub enum QualityStatus {
    Present { mean_q: f32, min_q: i32, max_q: i32 },
    Missing,
    Invalid { reason: String },
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

/// Optional extension trait for quality-aware read sources.
/// Provides access to phred quality scores (0-93 scale).
/// Default implementation returns None, indicating no quality available.
pub trait QualityAwareReadSource: ReadSource {
    /// Returns quality scores for a read (phred scale 0-93).
    /// Returns None if quality not available for this read.
    /// Precondition: returned slice must have same length as sequence.
    fn get_quality(&self, id: usize) -> Result<Option<Cow<'_, [i32]>>, ReadSourceError> {
        let _ = id;
        Ok(None)
    }

    /// Get quality scores for a windowed region (memory-efficient for large reads).
    /// Returns quality scores for bases [start..start+len).
    /// Returns None if quality not available.
    fn get_quality_window(
        &self,
        id: usize,
        start: usize,
        len: usize,
    ) -> Result<Option<Vec<i32>>, ReadSourceError> {
        // Default: fall back to full quality and slice
        if let Some(full_q) = self.get_quality(id)? {
            if start
                .checked_add(len)
                .map_or(false, |end| end <= full_q.len())
            {
                Ok(Some(full_q[start..start + len].to_vec()))
            } else {
                Err(ReadSourceError::InvalidRange)
            }
        } else {
            Ok(None)
        }
    }

    /// Validate all quality arrays match their sequence lengths.
    /// Returns QualityStatus describing validation outcome.
    fn validate_quality_alignment(&self) -> Result<QualityStatus, ReadSourceError> {
        let num_reads = self.num_reads()?;
        let mut quality_present = false;
        let mut all_qualities: Vec<i32> = Vec::new();

        for i in 0..num_reads {
            let seq_len = self.get_len(i)?;
            if let Some(q) = self.get_quality(i)? {
                quality_present = true;
                if q.len() != seq_len {
                    return Ok(QualityStatus::Invalid {
                        reason: format!(
                            "Read {}: quality len {} ≠ seq len {}",
                            i,
                            q.len(),
                            seq_len
                        ),
                    });
                }
                // Check phred range
                for &phred in q.iter() {
                    if phred < 0 || phred > 93 {
                        return Ok(QualityStatus::Invalid {
                            reason: format!("Read {}: quality {} out of range [0, 93]", i, phred),
                        });
                    }
                }
                all_qualities.extend(q.iter());
            }
        }

        if !quality_present {
            Ok(QualityStatus::Missing)
        } else {
            let mean_q = if all_qualities.is_empty() {
                0.0
            } else {
                all_qualities.iter().sum::<i32>() as f32 / all_qualities.len() as f32
            };
            let min_q = *all_qualities.iter().min().unwrap_or(&0);
            let max_q = *all_qualities.iter().max().unwrap_or(&0);

            Ok(QualityStatus::Present {
                mean_q,
                min_q,
                max_q,
            })
        }
    }
}

/// A simple in-memory adapter over Vec<String> for compatibility.
pub struct InMemoryReadSource {
    reads: Arc<Vec<String>>,
    names: Option<Arc<Vec<String>>>,
    qualities: Option<Arc<Vec<Vec<i32>>>>,
}

impl InMemoryReadSource {
    pub fn new(reads: Vec<String>, names: Option<Vec<String>>) -> Self {
        InMemoryReadSource {
            reads: Arc::new(reads),
            names: names.map(Arc::new),
            qualities: None,
        }
    }

    pub fn with_qualities(
        reads: Vec<String>,
        names: Option<Vec<String>>,
        qualities: Vec<Vec<i32>>,
    ) -> Self {
        InMemoryReadSource {
            reads: Arc::new(reads),
            names: names.map(Arc::new),
            qualities: Some(Arc::new(qualities)),
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

/// Implement QualityAwareReadSource for InMemoryReadSource.
impl QualityAwareReadSource for InMemoryReadSource {
    fn get_quality(&self, id: usize) -> Result<Option<Cow<'_, [i32]>>, ReadSourceError> {
        if let Some(qualities) = &self.qualities {
            Ok(qualities.get(id).map(|q| Cow::Borrowed(q.as_slice())))
        } else {
            Ok(None)
        }
    }

    fn get_quality_window(
        &self,
        id: usize,
        start: usize,
        len: usize,
    ) -> Result<Option<Vec<i32>>, ReadSourceError> {
        if let Some(qualities) = &self.qualities {
            if let Some(q) = qualities.get(id) {
                if start.checked_add(len).map_or(false, |end| end <= q.len()) {
                    Ok(Some(q[start..start + len].to_vec()))
                } else {
                    Err(ReadSourceError::InvalidRange)
                }
            } else {
                Err(ReadSourceError::NotFound(format!("id {}", id)))
            }
        } else {
            Ok(None)
        }
    }
}

/// Index entry for binary index, with optional quality metadata.
#[derive(Debug, Clone)]
pub struct IndexEntry {
    pub data_offset: u64,
    pub data_len: u32,
    pub name_idx: u32,
    pub flags: u32,
    pub qual_offset: Option<u64>,
    pub qual_len: Option<u32>,
}

/// A minimal BinaryIndexReadSource skeleton: attempts to mmap a `.seqs` file and read a simple JSON index.
#[allow(dead_code)]
pub struct BinaryIndexReadSource {
    seqs_path: PathBuf,
    sidx_path: PathBuf,
    mmap: Option<Mmap>,
    mmap_quals: Option<Mmap>,
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
        let quals_path = seqs_path.with_extension("quals");

        if !seqs_path.exists() || !sidx_path.exists() {
            return Err(ReadSourceError::NotFound(format!(
                "index files for {:?}",
                base
            )));
        }

        let file = File::open(&seqs_path)?;
        let mmap = unsafe { Mmap::map(&file).ok() };
        let mmap_quals = if quals_path.exists() {
            let qf = File::open(&quals_path)?;
            unsafe { Mmap::map(&qf).ok() }
        } else {
            None
        };

        // Read JSON index (simple array of {name,offset,len,qual_offset?,qual_len?}) for initial implementation.
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
                let qual_offset = v.get("qual_offset").and_then(|x| x.as_u64());
                let qual_len = v.get("qual_len").and_then(|x| x.as_u64()).map(|n| n as u32);
                let entry = IndexEntry {
                    data_offset: offset,
                    data_len: len,
                    name_idx: names.len() as u32,
                    flags: 0,
                    qual_offset,
                    qual_len,
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
            mmap_quals,
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

    // Prepare a temporary qualities file to store phred values
    let quals_tmp = out_base.with_extension("quals.tmp");
    let mut quals_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&quals_tmp)
        .map_err(ReadSourceError::Io)?;

    let mut records: Vec<serde_json::Value> = Vec::new();
    let mut offset: u64 = 0;
    let mut qual_offset: u64 = 0;

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

    // Process a file that may be FASTQ or FASTA
    let mut process_path = |p: &Path, is_reads2: bool| -> Result<(), ReadSourceError> {
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
                // Sequence (apply RC for reads2 if requested)
                let mut seq = rec.seq().to_vec();
                if is_reads2 && !no_revcomp {
                    seq = dna::revcomp(&seq);
                }
                // Phred qualities as 0-93 u8; reverse if RC was applied
                let mut qual_bytes: Vec<u8> = rec
                    .qual()
                    .iter()
                    .map(|&c| {
                        let phred = (c as i32) - 33;
                        phred.max(0).min(93) as u8
                    })
                    .collect();
                if is_reads2 && !no_revcomp {
                    qual_bytes.reverse();
                }
                let name = if is_reads2 && !no_revcomp {
                    format!("{}/RC", rec.id())
                } else {
                    rec.id().to_string()
                };
                // Write sequence
                seqs_file.write_all(&seq).map_err(ReadSourceError::Io)?;
                let seq_off = offset;
                offset += seq.len() as u64;
                // Write quality
                quals_file
                    .write_all(&qual_bytes)
                    .map_err(ReadSourceError::Io)?;
                let qoff = qual_offset;
                qual_offset += qual_bytes.len() as u64;
                // Push record with both offsets
                records.push(serde_json::json!({
                    "name": name,
                    "offset": seq_off,
                    "len": seq.len(),
                    "qual_offset": qoff,
                    "qual_len": qual_bytes.len()
                }));
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
                let name = if is_reads2 && !no_revcomp {
                    format!("{}/RC", rec.id())
                } else {
                    rec.id().to_string()
                };
                // Write sequence only (no quality for FASTA)
                seqs_file.write_all(&seq).map_err(ReadSourceError::Io)?;
                let seq_off = offset;
                offset += seq.len() as u64;
                records.push(serde_json::json!({
                    "name": name,
                    "offset": seq_off,
                    "len": seq.len(),
                    "qual_offset": null,
                    "qual_len": null
                }));
            }
            Ok(())
        }
    };

    // Process reads1 then reads2
    process_path(reads1, false)?;
    process_path(reads2, true)?;

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
    // Note: if FASTA inputs were used, quals_tmp may be empty; create file unconditionally
    let quals_tmp = out_base.with_extension("quals.tmp");
    if std::fs::metadata(&quals_tmp).is_ok() {
        let quals_final = out_base.with_extension("quals");
        std::fs::rename(&quals_tmp, &quals_final).map_err(ReadSourceError::Io)?;
    }
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

impl QualityAwareReadSource for BinaryIndexReadSource {
    fn get_quality(&self, id: usize) -> Result<Option<Cow<'_, [i32]>>, ReadSourceError> {
        let m = match &self.mmap_quals {
            Some(m) => m,
            None => return Ok(None),
        };
        let e = self
            .entries
            .get(id)
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))?;
        let (qoff, qlen) = match (e.qual_offset, e.qual_len) {
            (Some(qoff), Some(qlen)) => (qoff as usize, qlen as usize),
            _ => return Ok(None),
        };
        if qoff.checked_add(qlen).map_or(false, |end| end <= m.len()) {
            let slice = &m[qoff..qoff + qlen];
            let v: Vec<i32> = slice.iter().map(|&b| b as i32).collect();
            Ok(Some(Cow::Owned(v)))
        } else {
            Err(ReadSourceError::CorruptIndex)
        }
    }

    fn get_quality_window(
        &self,
        id: usize,
        start: usize,
        len: usize,
    ) -> Result<Option<Vec<i32>>, ReadSourceError> {
        let m = match &self.mmap_quals {
            Some(m) => m,
            None => return Ok(None),
        };
        let e = self
            .entries
            .get(id)
            .ok_or(ReadSourceError::NotFound(format!("id {}", id)))?;
        let (qoff, qlen) = match (e.qual_offset, e.qual_len) {
            (Some(qoff), Some(qlen)) => (qoff as usize, qlen as usize),
            _ => return Ok(None),
        };
        if start.checked_add(len).map_or(false, |end| end <= qlen) {
            let s = qoff + start;
            let slice = &m[s..s + len];
            Ok(Some(slice.iter().map(|&b| b as i32).collect()))
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

/// Load FASTQ sequences and quality scores, with optional reverse-complement for reads2.
/// Returns a tuple of (sequences, quality_scores, read_names).
pub fn load_fastq_with_quality<P: AsRef<Path>>(
    path: P,
    reverse_complement: bool,
) -> Result<(Vec<String>, Vec<Vec<i32>>, Vec<String>), ReadSourceError> {
    let path = path.as_ref();
    let file = File::open(path)?;

    let reader: Box<dyn std::io::Read> = if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let buf_reader = BufReader::new(reader);
    let fastq_reader = fastq::Reader::new(buf_reader);

    let mut sequences = Vec::new();
    let mut qualities = Vec::new();
    let mut names = Vec::new();

    for result in fastq_reader.records() {
        let record = result.map_err(|_| {
            ReadSourceError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                "FASTQ read error",
            ))
        })?;

        let mut seq = record.seq().to_vec();
        if reverse_complement {
            seq = dna::revcomp(&seq);
        }

        // Convert phred quality scores from ASCII (33-126) to integers (0-93)
        let qual_scores: Vec<i32> = record
            .qual()
            .iter()
            .map(|&c| {
                let phred = (c as i32) - 33;
                phred.max(0).min(93)
            })
            .collect();

        let name = if reverse_complement {
            format!("{}/RC", record.id())
        } else {
            record.id().to_string()
        };

        sequences.push(String::from_utf8_lossy(&seq).into_owned());
        qualities.push(qual_scores);
        names.push(name);
    }

    Ok((sequences, qualities, names))
}

/// Load a pair of FASTQ files as an InMemoryReadSource with quality support.
/// Reads1 are kept as-is; reads2 are reverse-complemented.
pub fn load_fastq_pair_with_quality<P: AsRef<Path>>(
    reads1_path: P,
    reads2_path: P,
) -> Result<Arc<dyn QualityAwareReadSource>, ReadSourceError> {
    let (seqs1, quals1, names1) = load_fastq_with_quality(&reads1_path, false)?;
    let (seqs2, quals2, names2) = load_fastq_with_quality(&reads2_path, true)?;

    let mut sequences = seqs1;
    sequences.extend(seqs2);

    let mut qualities = quals1;
    qualities.extend(quals2);

    let mut names = names1;
    names.extend(names2);

    Ok(Arc::new(InMemoryReadSource::with_qualities(
        sequences,
        Some(names),
        qualities,
    )))
}
