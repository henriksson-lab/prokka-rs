use std::io::Write;

use crate::error::ProkkaError;
use crate::model::AnnotationResult;

/// Write TSV summary output.
///
/// Replicates Perl Prokka lines 1311-1379.
/// Header: locus_tag ftype length_bp gene EC_number COG product
pub fn write_tsv(
    writer: &mut impl Write,
    result: &AnnotationResult,
) -> Result<(), ProkkaError> {
    writeln!(writer, "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct")?;

    for contig in &result.contigs {
        let mut sorted_features: Vec<_> = contig.features.iter().collect();
        sorted_features.sort_by(|a, b| {
            a.start.cmp(&b.start)
                .then(b.end.cmp(&a.end))
                .then(a.has_tag("Parent").cmp(&b.has_tag("Parent")))
        });

        for feature in &sorted_features {
            let locus_tag = feature.get_tag("locus_tag").unwrap_or("");
            let ftype = feature.feature_type.as_str();
            let length_bp = feature.length();
            let gene = feature.get_tag("gene").unwrap_or("");
            let ec = feature.get_tag("EC_number").unwrap_or("");

            // Extract COG from db_xref tag
            // Find COG among possibly multiple db_xref values
            let cog = feature.tags.get("db_xref")
                .and_then(|values| {
                    values.iter()
                        .find_map(|v| v.strip_prefix("COG:"))
                })
                .unwrap_or("");

            let product = feature.get_tag("product").unwrap_or("");

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                locus_tag, ftype, length_bp, gene, ec, cog, product
            )?;
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
                sequence: vec![],
                features,
            }],
            stats: AnnotationStats::default(),
            log: Vec::new(),
        }
    }

    #[test]
    fn test_tsv_header() {
        let result = make_result(vec![]);
        let mut buf = Vec::new();
        write_tsv(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert_eq!(
            output.lines().next().unwrap(),
            "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct"
        );
    }

    #[test]
    fn test_tsv_cds_fields() {
        let mut f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "t".into(), 1, 300, Strand::Forward,
        );
        f.add_tag("locus_tag", "T_00001");
        f.add_tag("product", "kinase");
        f.add_tag("gene", "abc");
        f.add_tag("EC_number", "1.2.3.4");
        f.add_tag("db_xref", "COG:COG0001");

        let result = make_result(vec![f]);
        let mut buf = Vec::new();
        write_tsv(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let data_line = output.lines().nth(1).unwrap();
        let fields: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(fields[0], "T_00001");
        assert_eq!(fields[1], "CDS");
        assert_eq!(fields[2], "300");
        assert_eq!(fields[3], "abc");
        assert_eq!(fields[4], "1.2.3.4");
        assert_eq!(fields[5], "COG0001");
        assert_eq!(fields[6], "kinase");
    }

    #[test]
    fn test_tsv_cog_from_multiple_db_xref() {
        let mut f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "t".into(), 1, 100, Strand::Forward,
        );
        f.add_tag("locus_tag", "T_00001");
        f.add_tag("product", "test");
        f.add_tag("db_xref", "GI:12345");     // first db_xref, not COG
        f.add_tag("db_xref", "COG:COG0042");  // second db_xref, is COG

        let result = make_result(vec![f]);
        let mut buf = Vec::new();
        write_tsv(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let data_line = output.lines().nth(1).unwrap();
        let fields: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(fields[5], "COG0042"); // should find COG in second db_xref
    }
}
