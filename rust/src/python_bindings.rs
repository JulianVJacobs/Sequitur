use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use crate::suffix::AffixArray;
use crate::overlap::{create_overlap_graph, OverlapConfig};
use crate::matching::{adjacency_to_csc, find_lower_diagonal_path};
use crate::matching::{detect_cycles};
use crate::overlap::{compute_edge_confidences};
use pyo3::types::{PyDict, PyList};

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
fn assemble_from_reads(reads: Vec<String>) -> PyResult<AssemblyResult> {
    if reads.is_empty() {
        return Err(PyValueError::new_err("reads list is empty"));
    }

    // Build affix array
    let affix = AffixArray::build(reads.iter().map(|s| s.as_str()), 3);

    // Create overlap graph with default config
    let config = OverlapConfig::default();
    let (_affix_out, adjacency, overlap_lengths) = create_overlap_graph(&reads, Some(affix), config);

    // Convert adjacency to sparse CSC matrix
    let mat = adjacency_to_csc(&adjacency, None);

    // Call assembler (qualities not provided)
    let seq = find_lower_diagonal_path(&mat, &overlap_lengths, &reads, None);

    Ok(AssemblyResult::new(seq))
}

#[pyfunction]
fn analyse_reads(py: Python, reads: Vec<String>) -> PyResult<PyObject> {
    if reads.is_empty() {
        return Err(PyValueError::new_err("reads list is empty"));
    }
    let affix = AffixArray::build(reads.iter().map(|s| s.as_str()), 3);
    let config = OverlapConfig::default();
    let (_affix_out, adjacency, overlap_lengths) = create_overlap_graph(&reads, Some(affix), config);

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

#[pymodule]
fn sequitur_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<AssemblyResult>()?;
    m.add_function(wrap_pyfunction!(assemble_from_reads, m)?)?;
    m.add_function(wrap_pyfunction!(analyse_reads, m)?)?;
    Ok(())
}
