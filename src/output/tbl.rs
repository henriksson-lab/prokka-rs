//! NCBI five-column feature table (`.tbl`) writer.
//!
//! This is the format consumed by `tbl2asn` to produce the GenBank/Sequin
//! outputs. Each contig is introduced by a `>Feature <id>` line, followed
//! by one block per feature: a `left<TAB>right<TAB>type` coordinate line
//! (with `left > right` indicating the minus strand) and any number of
//! triple-tab-indented qualifier lines (`<TAB><TAB><TAB>key<TAB>value`).

use std::io::Write;

use crate::error::ProkkaError;
use crate::model::{AnnotationResult, Strand};

/// Write the NCBI five-column feature table (`.tbl`).
///
/// Features inside a contig are sorted by `(start asc, end desc, has Parent
/// asc)` to match the Perl Prokka comparator at line 1329. Coordinates are
/// written as `start<TAB>end` for the forward strand and `end<TAB>start`
/// for the reverse strand, as expected by `tbl2asn`. Tags whose name starts
/// with an uppercase letter (GFF-specific things like `ID`, `Parent`,
/// `Name`) are skipped, except for `EC_number` — matching Perl Prokka line
/// 1345.
///
/// Replicates Perl Prokka lines 1328-1349.
pub fn write_tbl(
    writer: &mut impl Write,
    result: &AnnotationResult,
) -> Result<(), ProkkaError> {
    for contig in &result.contigs {
        writeln!(writer, ">Feature {}", contig.id)?;

        let mut sorted_features: Vec<_> = contig.features.iter().collect();
        sorted_features.sort_by(|a, b| {
            a.start.cmp(&b.start)
                .then(b.end.cmp(&a.end))
                .then(a.has_tag("Parent").cmp(&b.has_tag("Parent")))
        });

        for feature in &sorted_features {
            // Coordinates: L R type (swapped for reverse strand)
            let (left, right) = match feature.strand {
                Strand::Forward => (feature.start, feature.end),
                Strand::Reverse => (feature.end, feature.start),
            };
            writeln!(writer, "{}\t{}\t{}", left, right, feature.feature_type.as_str())?;

            // Tags: skip GFF-specific tags (start with uppercase, except EC_number)
            for (key, values) in &feature.tags {
                // Skip tags starting with uppercase unless it's EC_number
                let first_char = key.chars().next().unwrap_or('a');
                if first_char.is_ascii_uppercase() && key != "EC_number" {
                    continue;
                }
                for value in values {
                    writeln!(writer, "\t\t\t{}\t{}", key, value)?;
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::*;

    fn make_result(features: Vec<SeqFeature>) -> AnnotationResult {
        AnnotationResult {
            contigs: vec![Contig {
                id: "ctg1".into(),
                sequence: b"ACGTACGT".to_vec(),
                features,
            }],
            stats: AnnotationStats::default(),
            log: Vec::new(),
        }
    }

    #[test]
    fn test_tbl_forward_strand() {
        let mut f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "test".into(), 100, 400, Strand::Forward,
        );
        f.add_tag("product", "test protein");
        f.add_tag("locus_tag", "T_00001");

        let result = make_result(vec![f]);
        let mut buf = Vec::new();
        write_tbl(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains(">Feature ctg1\n"));
        assert!(output.contains("100\t400\tCDS\n")); // forward: start < end
        assert!(output.contains("\t\t\tproduct\ttest protein\n"));
        assert!(output.contains("\t\t\tlocus_tag\tT_00001\n"));
    }

    #[test]
    fn test_tbl_reverse_strand() {
        let f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "test".into(), 100, 400, Strand::Reverse,
        );
        let result = make_result(vec![f]);
        let mut buf = Vec::new();
        write_tbl(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // Reverse strand: end comes first
        assert!(output.contains("400\t100\tCDS\n"));
    }

    #[test]
    fn test_tbl_skips_uppercase_tags_except_ec() {
        let mut f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "test".into(), 1, 100, Strand::Forward,
        );
        f.add_tag("ID", "should_skip");
        f.add_tag("Parent", "should_skip");
        f.add_tag("Name", "should_skip");
        f.add_tag("EC_number", "1.2.3.4");
        f.add_tag("product", "test");

        let result = make_result(vec![f]);
        let mut buf = Vec::new();
        write_tbl(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(!output.contains("ID\t"));
        assert!(!output.contains("Parent\t"));
        assert!(!output.contains("Name\t"));
        assert!(output.contains("EC_number\t1.2.3.4"));
        assert!(output.contains("product\ttest"));
    }
}
