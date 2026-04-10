use std::collections::HashMap;

use crate::model::Contig;

/// De-duplicate colliding /gene names across all contigs.
///
/// Groups CDS features by gene name (stripping existing `_N` suffix),
/// and renames duplicates as `gene_1`, `gene_2`, etc.
///
/// Replicates Perl Prokka lines 1201-1226.
pub fn deduplicate_genes(contigs: &mut [Contig]) {
    // Collect all gene names and their (contig_idx, feature_idx) locations
    let mut gene_locs: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    for (ci, contig) in contigs.iter().enumerate() {
        for (fi, feature) in contig.features.iter().enumerate() {
            if feature.feature_type != crate::model::FeatureType::CDS {
                continue;
            }
            if let Some(gene) = feature.get_tag("gene") {
                // Strip existing _N suffix for grouping
                let base = gene
                    .rfind('_')
                    .and_then(|pos| {
                        if gene[pos + 1..].chars().all(|c| c.is_ascii_digit()) {
                            Some(&gene[..pos])
                        } else {
                            None
                        }
                    })
                    .unwrap_or(gene)
                    .to_string();

                gene_locs
                    .entry(base)
                    .or_default()
                    .push((ci, fi));
            }
        }
    }

    // Rename duplicates
    for (gene, locs) in &gene_locs {
        if locs.len() <= 1 {
            continue;
        }
        for (n, &(ci, fi)) in locs.iter().enumerate() {
            let feature = &mut contigs[ci].features[fi];
            feature.remove_tag("gene");
            feature.add_tag("gene", &format!("{}_{}", gene, n + 1));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{FeatureType, SeqFeature, Strand};

    fn cds_with_gene(start: usize, gene: &str) -> SeqFeature {
        let mut f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "test".into(), start, start + 300, Strand::Forward,
        );
        f.add_tag("gene", gene);
        f
    }

    #[test]
    fn test_no_duplicates_unchanged() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        contigs[0].features = vec![
            cds_with_gene(100, "abc"),
            cds_with_gene(500, "def"),
            cds_with_gene(900, "ghi"),
        ];
        deduplicate_genes(&mut contigs);
        assert_eq!(contigs[0].features[0].get_tag("gene"), Some("abc"));
        assert_eq!(contigs[0].features[1].get_tag("gene"), Some("def"));
        assert_eq!(contigs[0].features[2].get_tag("gene"), Some("ghi"));
    }

    #[test]
    fn test_duplicates_renamed() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        contigs[0].features = vec![
            cds_with_gene(100, "abc"),
            cds_with_gene(500, "abc"),
            cds_with_gene(900, "abc"),
        ];
        deduplicate_genes(&mut contigs);
        assert_eq!(contigs[0].features[0].get_tag("gene"), Some("abc_1"));
        assert_eq!(contigs[0].features[1].get_tag("gene"), Some("abc_2"));
        assert_eq!(contigs[0].features[2].get_tag("gene"), Some("abc_3"));
    }

    #[test]
    fn test_existing_suffix_stripped_for_grouping() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        contigs[0].features = vec![
            cds_with_gene(100, "abc_1"),
            cds_with_gene(500, "abc_2"),
        ];
        deduplicate_genes(&mut contigs);
        // Both have base "abc", so they get renumbered
        assert_eq!(contigs[0].features[0].get_tag("gene"), Some("abc_1"));
        assert_eq!(contigs[0].features[1].get_tag("gene"), Some("abc_2"));
    }

    #[test]
    fn test_no_gene_tag_ignored() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        let f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "test".into(), 100, 400, Strand::Forward,
        );
        // No gene tag
        contigs[0].features = vec![f];
        deduplicate_genes(&mut contigs);
        assert!(!contigs[0].features[0].has_tag("gene"));
    }

    #[test]
    fn test_cross_contig_duplicates() {
        let mut contigs = vec![
            Contig::new("ctg1".into(), b"ACGT".to_vec()),
            Contig::new("ctg2".into(), b"ACGT".to_vec()),
        ];
        contigs[0].features = vec![cds_with_gene(100, "abc")];
        contigs[1].features = vec![cds_with_gene(100, "abc")];
        deduplicate_genes(&mut contigs);
        // Same gene across contigs should be deduplicated
        assert_eq!(contigs[0].features[0].get_tag("gene"), Some("abc_1"));
        assert_eq!(contigs[1].features[0].get_tag("gene"), Some("abc_2"));
    }
}
