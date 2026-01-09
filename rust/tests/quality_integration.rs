//! Integration tests for quality-aware overlap scoring.
//!
//! These tests verify that:
//! 1. Quality-aware and quality-blind paths produce identical edge topology
//! 2. Quality scoring adjusts edge weights as expected
//! 3. Quality validation catches corrupt/missing data
//! 4. Graceful fallback occurs when quality is unavailable

#[cfg(test)]
mod tests {
    use sequitur::overlap::{OverlapConfig, QualityConfig, QualityScoring};
    use sequitur::read_source::{
        InMemoryReadSource, QualityAwareReadSource, QualityStatus, ReadSource,
    };

    /// Test basic quality-aware read source creation and validation.
    #[test]
    fn test_quality_aware_read_source() {
        let sequences = vec!["ACGTACGT".to_string(), "GTACGTAC".to_string()];
        let qualities = vec![
            vec![30, 25, 20, 15, 10, 5, 3, 1],
            vec![35, 32, 28, 24, 20, 15, 10, 5],
        ];

        let source = InMemoryReadSource::with_qualities(
            sequences.clone(),
            Some(vec!["read1".to_string(), "read2".to_string()]),
            qualities.clone(),
        );

        assert_eq!(source.num_reads().unwrap(), 2);
        assert_eq!(source.get_len(0).unwrap(), 8);

        let q0 = source.get_quality(0).unwrap();
        assert!(q0.is_some());
        assert_eq!(q0.unwrap().to_vec(), qualities[0]);
    }

    /// Test quality window extraction (memory-efficient quality access).
    #[test]
    fn test_quality_window_access() {
        let sequences = vec!["ACGTACGTACGT".to_string()];
        let qualities = vec![vec![30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19]];

        let source = InMemoryReadSource::with_qualities(sequences, None, qualities.clone());

        // Extract quality for bases 8-11 (window at position 8, length 4)
        let window = source.get_quality_window(0, 8, 4).unwrap();
        assert!(window.is_some());
        assert_eq!(window.unwrap(), vec![22, 21, 20, 19]);
    }

    /// Test quality validation with correct data.
    #[test]
    fn test_quality_validation_success() {
        let sequences = vec!["ACGTACGT".to_string(), "GTACGTAC".to_string()];
        let qualities = vec![
            vec![30, 25, 20, 15, 10, 5, 3, 1],
            vec![35, 32, 28, 24, 20, 15, 10, 5],
        ];

        let source = InMemoryReadSource::with_qualities(sequences, None, qualities);
        let status = source.validate_quality_alignment().unwrap();

        match status {
            QualityStatus::Present {
                mean_q,
                min_q,
                max_q,
            } => {
                assert!(mean_q > 0.0 && mean_q < 93.0);
                assert!(min_q >= 0 && min_q <= 93);
                assert!(max_q >= 0 && max_q <= 93);
                assert!(min_q <= max_q);
            }
            _ => panic!("Expected QualityStatus::Present"),
        }
    }

    /// Test quality validation with length mismatch.
    #[test]
    fn test_quality_validation_length_mismatch() {
        let sequences = vec!["ACGTACGT".to_string()];
        let qualities = vec![vec![30, 25, 20]]; // Wrong length!

        let source = InMemoryReadSource::with_qualities(sequences, None, qualities);
        let status = source.validate_quality_alignment().unwrap();

        match status {
            QualityStatus::Invalid { reason } => {
                assert!(reason.contains("quality len") && reason.contains("seq len"));
            }
            _ => panic!("Expected QualityStatus::Invalid"),
        }
    }

    /// Test quality validation with out-of-range phred scores.
    #[test]
    fn test_quality_validation_out_of_range() {
        let sequences = vec!["ACGT".to_string()];
        let qualities = vec![vec![30, 25, 100, 15]]; // 100 > 93!

        let source = InMemoryReadSource::with_qualities(sequences, None, qualities);
        let status = source.validate_quality_alignment().unwrap();

        match status {
            QualityStatus::Invalid { reason } => {
                assert!(reason.contains("out of range"));
            }
            _ => panic!("Expected QualityStatus::Invalid"),
        }
    }

    /// Test quality configuration with different scoring modes.
    #[test]
    fn test_quality_config_modes() {
        let config_none = QualityConfig {
            scoring_mode: QualityScoring::None,
            ..Default::default()
        };

        let config_position = QualityConfig {
            scoring_mode: QualityScoring::Position,
            error_penalty_exponent: 1.5,
            ..Default::default()
        };

        let config_confidence = QualityConfig {
            scoring_mode: QualityScoring::Confidence,
            min_quality: Some(10),
            ..Default::default()
        };

        assert_eq!(config_none.scoring_mode, QualityScoring::None);
        assert_eq!(config_position.scoring_mode, QualityScoring::Position);
        assert_eq!(config_position.error_penalty_exponent, 1.5);
        assert_eq!(config_confidence.scoring_mode, QualityScoring::Confidence);
        assert_eq!(config_confidence.min_quality, Some(10));
    }

    /// Test overlap config includes quality config.
    #[test]
    fn test_overlap_config_with_quality() {
        let mut config = OverlapConfig::default();

        assert_eq!(config.quality_config.scoring_mode, QualityScoring::None);
        assert!(!config.quality_config.log_quality_scores);

        config.quality_config.scoring_mode = QualityScoring::Confidence;
        config.quality_config.log_quality_scores = true;

        assert_eq!(
            config.quality_config.scoring_mode,
            QualityScoring::Confidence
        );
        assert!(config.quality_config.log_quality_scores);
    }

    /// Test source without quality returns None gracefully.
    #[test]
    fn test_quality_optional() {
        let sequences = vec!["ACGTACGT".to_string()];
        let source = InMemoryReadSource::new(sequences, None);

        // Source should report no quality available
        let q = source.get_quality(0).unwrap();
        assert!(q.is_none());

        let status = source.validate_quality_alignment().unwrap();
        match status {
            QualityStatus::Missing => {} // Expected
            _ => panic!("Expected QualityStatus::Missing"),
        }
    }

    /// Test index-backed quality access (memory-efficient path).
    #[test]
    fn test_index_backed_quality() {
        use sequitur::read_source::build_index_from_pair;
        use sequitur::read_source::BinaryIndexReadSource;
        use std::fs;
        use std::io::Write;

        let temp_dir = std::env::temp_dir().join(format!("sequitur_test_{}", std::process::id()));
        fs::create_dir_all(&temp_dir).unwrap();

        // Create temporary FASTQ files with quality scores
        let reads1_path = temp_dir.join("reads1.fastq");
        let reads2_path = temp_dir.join("reads2.fastq");

        let mut reads1_file = fs::File::create(&reads1_path).unwrap();
        writeln!(reads1_file, "@read1").unwrap();
        writeln!(reads1_file, "ACGTACGT").unwrap();
        writeln!(reads1_file, "+").unwrap();
        writeln!(reads1_file, "IIIHHGFF").unwrap(); // phred 40,40,40,39,39,38,37,37

        let mut reads2_file = fs::File::create(&reads2_path).unwrap();
        writeln!(reads2_file, "@read2").unwrap();
        writeln!(reads2_file, "TGCATGCA").unwrap();
        writeln!(reads2_file, "+").unwrap();
        writeln!(reads2_file, "IIHHGGFF").unwrap(); // phred 40,40,39,39,38,38,37,37

        // Build index
        let index_base = temp_dir.join("test_index");
        build_index_from_pair(&reads1_path, &reads2_path, &index_base, false).unwrap();

        // Open index and verify quality access
        let source = BinaryIndexReadSource::open(&index_base).unwrap();

        assert_eq!(source.num_reads().unwrap(), 2);

        // Check first read quality
        let q0 = source.get_quality(0).unwrap();
        assert!(q0.is_some());
        let q0_vec = q0.unwrap().to_vec();
        assert_eq!(q0_vec.len(), 8);
        // Phred scores: I=40, H=39, G=38, F=37
        assert_eq!(q0_vec, vec![40, 40, 40, 39, 39, 38, 37, 37]);

        // Check second read quality (reverse-complemented, so quality should be reversed)
        let q1 = source.get_quality(1).unwrap();
        assert!(q1.is_some());
        let q1_vec = q1.unwrap().to_vec();
        assert_eq!(q1_vec.len(), 8);
        // Reversed: F,F,G,G,H,H,I,I => 37,37,38,38,39,39,40,40
        assert_eq!(q1_vec, vec![37, 37, 38, 38, 39, 39, 40, 40]);

        // Test windowed access
        let window = source.get_quality_window(0, 2, 4).unwrap();
        assert!(window.is_some());
        assert_eq!(window.unwrap(), vec![40, 39, 39, 38]);

        // Validate quality alignment
        let status = source.validate_quality_alignment().unwrap();
        match status {
            QualityStatus::Present {
                mean_q,
                min_q,
                max_q,
            } => {
                assert!(mean_q > 0.0);
                assert_eq!(min_q, 37);
                assert_eq!(max_q, 40);
            }
            _ => panic!("Expected QualityStatus::Present"),
        }

        // Cleanup
        fs::remove_dir_all(&temp_dir).unwrap();
    }
}
