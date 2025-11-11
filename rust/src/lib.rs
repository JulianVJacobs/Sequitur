//! sequitur_rs library skeleton
//!
//! This crate will expose the minimal core primitives needed for the Sequitur
//! algorithm: suffix-array construction, overlap scoring and adjacency
//! representation, and the reconstruction algorithm.

pub mod suffix;
pub mod overlap;
pub mod matching;

pub use suffix::{
    AffixArray, AffixEntry, AffixKind, DEFAULT_MIN_SUFFIX_LEN,
};
pub use overlap::{
    create_overlap_graph,
    normalised_damerau_levenshtein_distance,
    Adjacency,
    OverlapConfig,
    OverlapLengths,
};
pub use matching::{
    adjacency_to_csc,
    adjacency_to_sparse,
    relabel_columns,
    relabel_rows,
    find_lower_diagonal_path,
};

// TODO: export a compact public API mirroring the Python core
// e.g. pub fn build_suffix_array(...), pub fn create_adjacency(...)

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
