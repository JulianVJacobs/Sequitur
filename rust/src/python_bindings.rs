use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::affix::AffixMap;
use crate::alternative_paths::analyse_alternatives;
use crate::matching::detect_cycles;
use crate::matching::find_first_subdiagonal_path;
use crate::overlap::compute_edge_confidences;
use crate::overlap::{create_overlap_graph, Adjacency, OverlapConfig};
use pyo3::types::{PyDict, PyList};
use sprs::CsMat;
use std::collections::HashMap;

/// Local helper for Python bindings that still need HashMap format
fn sparse_to_adjacency_local(matrix: &CsMat<usize>) -> Adjacency {
    let (rows, _) = matrix.shape();
    let mut adjacency: Adjacency = vec![HashMap::new(); rows];

    let csr = matrix.to_csr();
    for (row_idx, row_vec) in csr.outer_iterator().enumerate() {
        for (col_idx, &val) in row_vec.indices().iter().zip(row_vec.data().iter()) {
            adjacency[row_idx].insert(*col_idx, val);
        }
    }

    adjacency
}

#[pyclass]
pub struct AssemblyResult {
    #[pyo3(get)]
    pub best: String,
}

#[pymethods]
impl AssemblyResult {
    #[new]
    fn new(best: String) -> Self {
        AssemblyResult { best }
    }
}

#[pyfunction]
fn assemble_from_reads(
    reads: Vec<String>,
    use_threads: Option<bool>,
    max_workers: Option<usize>,
    use_trie: Option<bool>,
) -> PyResult<AssemblyResult> {
    if reads.is_empty() {
        return Err(PyValueError::new_err("reads list is empty"));
    }

    // Build affix map
    let affix_map = AffixMap::build(&reads, 3);

    // Create overlap graph with threading config
    let config = OverlapConfig {
        use_threads: use_threads.unwrap_or(false),
        max_workers: max_workers.unwrap_or(1),
        use_trie: use_trie.unwrap_or(true),
        ..OverlapConfig::default()
    };
    let (_affix_out, adjacency_matrix, overlap_matrix) =
        create_overlap_graph(&reads, Some(affix_map), config);

    // Convert to CSC
    let adjacency_csc = adjacency_matrix.to_csc();
    let overlap_csc = overlap_matrix.to_csc();

    // Call assembler (qualities not provided)
    let seq = find_first_subdiagonal_path(&adjacency_csc, &overlap_csc, &reads, None);

    Ok(AssemblyResult::new(seq))
}

#[pyfunction]
fn analyse_reads(
    py: Python,
    reads: Vec<String>,
    use_threads: Option<bool>,
    max_workers: Option<usize>,
    use_trie: Option<bool>,
) -> PyResult<PyObject> {
    if reads.is_empty() {
        return Err(PyValueError::new_err("reads list is empty"));
    }
    let affix_map = AffixMap::build(&reads, 3);
    let config = OverlapConfig {
        use_threads: use_threads.unwrap_or(false),
        max_workers: max_workers.unwrap_or(1),
        use_trie: use_trie.unwrap_or(true),
        ..OverlapConfig::default()
    };
    let (_affix_out, adjacency_matrix, overlap_matrix) =
        create_overlap_graph(&reads, Some(affix_map), config);

    let adjacency = sparse_to_adjacency_local(&adjacency_matrix);
    let overlap_lengths = sparse_to_adjacency_local(&overlap_matrix);
    let confidences = compute_edge_confidences(&adjacency);
    let cycles = detect_cycles(&adjacency);

    let py_adj = PyList::empty(py);
    for (src_idx, src_edges) in adjacency.iter().enumerate() {
        let py_row = PyList::empty(py);
        for (&dst, &w) in src_edges.iter() {
            // include overlap length if present in overlap_lengths
            let overlap_len = overlap_lengths
                .get(src_idx)
                .and_then(|m| m.get(&dst))
                .copied()
                .unwrap_or(0usize);
            py_row.append((dst, w, overlap_len))?;
        }
        py_adj.append(py_row)?;
    }

    let py_conf = PyList::empty(py);
    for src in confidences.iter() {
        let py_row = PyList::empty(py);
        for (dst, p) in src.iter() {
            py_row.append((dst, *p))?;
        }
        py_conf.append(py_row)?;
    }

    let py_cycles = PyList::empty(py);
    for comp in cycles.iter() {
        let py_comp = PyList::empty(py);
        for v in comp.iter() {
            py_comp.append(*v)?;
        }
        py_cycles.append(py_comp)?;
    }

    let dict = PyDict::new(py);
    dict.set_item("adjacency", py_adj)?;
    dict.set_item("confidences", py_conf)?;
    dict.set_item("cycles", py_cycles)?;
    Ok(dict.into())
}

#[pyfunction]
fn analyse_alternative_paths(
    py: Python,
    reads: Vec<String>,
    score_gap: Option<f64>,
    use_threads: Option<bool>,
    max_workers: Option<usize>,
    use_trie: Option<bool>,
) -> PyResult<PyObject> {
    if reads.is_empty() {
        return Err(PyValueError::new_err("reads list is empty"));
    }

    // Build overlap graph
    let affix_map = AffixMap::build(&reads, 3);
    let config = OverlapConfig {
        use_threads: use_threads.unwrap_or(false),
        max_workers: max_workers.unwrap_or(1),
        use_trie: use_trie.unwrap_or(true),
        ..OverlapConfig::default()
    };
    let (_affix_out, adjacency_matrix, _overlap_matrix) =
        create_overlap_graph(&reads, Some(affix_map), config);

    // Convert to CSC
    let adjacency_csc = adjacency_matrix.to_csc();

    // Analyse alternatives
    let analysis = analyse_alternatives(&adjacency_csc, score_gap);

    // Convert to Python objects
    let py_squares = PyList::empty(py);
    for square in analysis.squares.iter() {
        let py_dict = PyDict::new(py);
        py_dict.set_item("i", square.i)?;
        py_dict.set_item("j", square.j)?;
        py_dict.set_item("delta", square.delta)?;
        py_squares.append(py_dict)?;
    }

    let py_components = PyList::empty(py);
    for component in analysis.components.iter() {
        let py_comp = PyList::empty(py);
        for &pos in component.iter() {
            py_comp.append(pos)?;
        }
        py_components.append(py_comp)?;
    }

    let py_cycles = PyList::empty(py);
    for cycle in analysis.cycles.iter() {
        let py_cycle = PyList::empty(py);
        for &pos in cycle.iter() {
            py_cycle.append(pos)?;
        }
        py_cycles.append(py_cycle)?;
    }

    let py_chains = PyList::empty(py);
    for chain in analysis.chains.iter() {
        let py_chain = PyList::empty(py);
        for &pos in chain.iter() {
            py_chain.append(pos)?;
        }
        py_chains.append(py_chain)?;
    }

    let dict = PyDict::new(py);
    dict.set_item("squares", py_squares)?;
    dict.set_item("components", py_components)?;
    dict.set_item("cycles", py_cycles)?;
    dict.set_item("chains", py_chains)?;
    dict.set_item("ambiguity_count", analysis.ambiguity_count)?;

    Ok(dict.into())
}

#[pymodule]
fn sequitur_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<AssemblyResult>()?;
    m.add_function(wrap_pyfunction!(assemble_from_reads, m)?)?;
    m.add_function(wrap_pyfunction!(analyse_reads, m)?)?;
    m.add_function(wrap_pyfunction!(analyse_alternative_paths, m)?)?;
    Ok(())
}
