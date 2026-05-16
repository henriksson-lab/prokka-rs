//! Plain-text annotation statistics (`.txt`) writer.
//!
//! Reproduces the small `prefix.txt` summary file Perl Prokka writes at
//! lines 1394-1413: a single `organism:` line, contig count, total bases,
//! and one `<feature_type>: <count>` line for each feature type (sorted
//! alphabetically by feature type).

use std::io::Write;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::AnnotationResult;

/// Write the annotation statistics text file (`.txt`).
///
/// Feature counts are sorted alphabetically by feature-type name to match
/// Perl Prokka's `sort keys %count` (line 1408). Replicates Perl Prokka
/// lines 1394-1413.
pub fn write_txt(
    writer: &mut impl Write,
    result: &AnnotationResult,
    config: &ProkkaConfig,
) -> Result<(), ProkkaError> {
    let plasmid = config.plasmid.as_deref().unwrap_or("");

    writeln!(
        writer,
        "organism: {} {} {} {}",
        config.genus, config.species, config.strain, plasmid
    )?;
    writeln!(writer, "contigs: {}", result.stats.num_contigs)?;
    writeln!(writer, "bases: {}", result.stats.total_bp)?;

    // Feature counts, sorted by name
    let mut counts: Vec<_> = result.stats.feature_counts.iter().collect();
    counts.sort_by_key(|(k, _)| (*k).clone());
    for (ftype, count) in counts {
        writeln!(writer, "{}: {}", ftype, count)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::IndexMap;

    #[test]
    fn test_txt_format() {
        let mut feature_counts = IndexMap::new();
        feature_counts.insert("CDS".to_string(), 100usize);
        feature_counts.insert("tRNA".to_string(), 50usize);

        let result = crate::model::AnnotationResult {
            contigs: vec![],
            stats: crate::model::AnnotationStats {
                num_contigs: 3,
                total_bp: 5000000,
                feature_counts,
            },
            log: Vec::new(),
        };
        let config = ProkkaConfig {
            genus: "Escherichia".into(),
            species: "coli".into(),
            strain: "K12".into(),
            plasmid: Some("pUC19".into()),
            ..Default::default()
        };

        let mut buf = Vec::new();
        write_txt(&mut buf, &result, &config).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("organism: Escherichia coli K12 pUC19"));
        assert!(output.contains("contigs: 3"));
        assert!(output.contains("bases: 5000000"));
        assert!(output.contains("CDS: 100"));
        assert!(output.contains("tRNA: 50"));
    }
}
