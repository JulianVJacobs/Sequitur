//! sequitur_rs library skeleton
//!
//! This crate will expose the minimal core primitives needed for the Sequitur
//! algorithm: suffix-array construction, overlap scoring and adjacency
//! representation, and the reconstruction algorithm.

pub mod affix;
pub mod alternative_paths;
pub mod matching;
pub mod overlap;
pub mod read_source;

pub use affix::{AffixNode, AffixSlice, Direction, PrunedAffixTrie, DEFAULT_MIN_SUFFIX_LEN};
pub mod python_bindings;
pub use alternative_paths::{
    analyse_alternatives, build_swap_graph, detect_swap_squares, extract_read_alternatives,
    find_connected_components, is_cycle, AlternativesAnalysis, ReadAlternative, ReadAlternatives,
    SwapSquare,
};
pub use matching::{
    adjacency_to_csc, adjacency_to_sparse, find_first_subdiagonal_path,
    find_first_subdiagonal_path_with_options, relabel_columns, relabel_rows, AssemblyOptions,
    AssemblyResult,
};
pub use overlap::{
    compute_edge_confidences, create_overlap_graph, create_overlap_graph_unified,
    create_overlap_graph_unified_from_readsource, create_overlap_graph_with_ambiguities,
    create_overlap_graph_with_ambiguities_from_readsource, normalised_damerau_levenshtein_distance,
    Adjacency, AmbiguityInfo, OverlapConfig, OverlapGraphResult, OverlapLengths,
};
pub use read_source::{
    build_index_from_pair, BinaryIndexReadSource, InMemoryReadSource, QualityAwareReadSource,
    QualityStatus, ReadSource, ReadSourceError,
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
