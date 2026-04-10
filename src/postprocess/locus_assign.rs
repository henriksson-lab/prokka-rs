use crate::config::ProkkaConfig;
use crate::model::{Contig, FeatureType, SeqFeature};

/// Assign locus_tags and protein_ids to all CDS and RNA features.
///
/// Optionally creates sister `gene` and `mRNA` features.
///
/// Replicates Perl Prokka lines 1231-1293.
pub fn assign_locus_tags(contigs: &mut [Contig], config: &ProkkaConfig) {
    let locustag = config.locustag.as_deref().unwrap_or("PROKKA");
    let increment = config.increment;
    let mut counter: u32 = 0;

    // We collect new features to add (gene/mRNA sisters) separately to avoid borrow issues
    let mut new_features_per_contig: Vec<Vec<SeqFeature>> = vec![Vec::new(); contigs.len()];

    for (ci, contig) in contigs.iter_mut().enumerate() {
        // Sort features by start position
        contig.features.sort_by(|a, b| a.start.cmp(&b.start));

        for feature in &mut contig.features {
            // Only assign to CDS and RNA features
            let is_cds_or_rna = matches!(
                feature.feature_type,
                FeatureType::CDS
                    | FeatureType::TRNA
                    | FeatureType::TMRNA
                    | FeatureType::RRNA
                    | FeatureType::MiscRNA
            );
            if !is_cds_or_rna {
                continue;
            }

            counter += 1;
            let id = format!("{}_{:05}", locustag, counter * increment);

            feature.add_tag("ID", &id);
            feature.add_tag("locus_tag", &id);

            // CDS features must have a valid /protein_id for tbl2asn
            if let Some(ref centre) = config.centre {
                if feature.feature_type == FeatureType::CDS {
                    feature.add_tag("protein_id", &format!("gnl|{}|{}", centre, id));
                }
            }

            // Create sister gene feature if --addgenes
            if config.addgenes {
                let gene_id = format!("{}_gene", id);
                let mut gene = SeqFeature::new(
                    FeatureType::Gene,
                    feature.seq_id.clone(),
                    "prokka-rs".to_string(),
                    feature.start,
                    feature.end,
                    feature.strand,
                );
                gene.add_tag("locus_tag", &id);
                gene.add_tag("ID", &gene_id);

                // Copy /gene from the CDS
                if let Some(g) = feature.get_tag("gene") {
                    gene.add_tag("gene", g);
                }

                // Add Parent tag to CDS pointing to gene
                feature.add_tag("Parent", &gene_id);

                new_features_per_contig[ci].push(gene);
            }

            // Create sister mRNA feature if --addmrna
            if config.addmrna {
                let mrna_id = format!("{}_mRNA", id);
                let mut mrna = SeqFeature::new(
                    FeatureType::MRNA,
                    feature.seq_id.clone(),
                    "prokka-rs".to_string(),
                    feature.start,
                    feature.end,
                    feature.strand,
                );
                mrna.add_tag("locus_tag", &id);
                mrna.add_tag("ID", &mrna_id);

                if let Some(g) = feature.get_tag("gene") {
                    mrna.add_tag("gene", g);
                }

                feature.add_tag("Parent", &mrna_id);

                new_features_per_contig[ci].push(mrna);
            }
        }
    }

    // Add new features to contigs
    for (ci, new_features) in new_features_per_contig.into_iter().enumerate() {
        contigs[ci].features.extend(new_features);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Strand;

    fn make_cds(start: usize, end: usize) -> SeqFeature {
        SeqFeature::new(FeatureType::CDS, "ctg1".into(), "test".into(), start, end, Strand::Forward)
    }

    fn make_trna(start: usize, end: usize) -> SeqFeature {
        SeqFeature::new(FeatureType::TRNA, "ctg1".into(), "test".into(), start, end, Strand::Forward)
    }

    #[test]
    fn test_basic_locus_tags() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        contigs[0].features = vec![
            make_cds(100, 400),
            make_cds(500, 800),
            make_trna(900, 970),
        ];
        let config = ProkkaConfig {
            locustag: Some("TEST".into()),
            increment: 1,
            ..Default::default()
        };
        assign_locus_tags(&mut contigs, &config);

        assert_eq!(contigs[0].features[0].get_tag("locus_tag"), Some("TEST_00001"));
        assert_eq!(contigs[0].features[1].get_tag("locus_tag"), Some("TEST_00002"));
        assert_eq!(contigs[0].features[2].get_tag("locus_tag"), Some("TEST_00003"));
    }

    #[test]
    fn test_increment() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        contigs[0].features = vec![make_cds(100, 400), make_cds(500, 800)];
        let config = ProkkaConfig {
            locustag: Some("X".into()),
            increment: 10,
            ..Default::default()
        };
        assign_locus_tags(&mut contigs, &config);

        assert_eq!(contigs[0].features[0].get_tag("locus_tag"), Some("X_00010"));
        assert_eq!(contigs[0].features[1].get_tag("locus_tag"), Some("X_00020"));
    }

    #[test]
    fn test_addgenes_creates_sister_features() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        let mut cds = make_cds(100, 400);
        cds.add_tag("gene", "abc");
        contigs[0].features = vec![cds];

        let config = ProkkaConfig {
            locustag: Some("T".into()),
            addgenes: true,
            ..Default::default()
        };
        assign_locus_tags(&mut contigs, &config);

        // Should have CDS + sister gene = 2 features
        assert_eq!(contigs[0].features.len(), 2);
        assert_eq!(contigs[0].features[1].feature_type, FeatureType::Gene);
        assert_eq!(contigs[0].features[1].get_tag("locus_tag"), Some("T_00001"));
        assert_eq!(contigs[0].features[1].get_tag("gene"), Some("abc"));
        // CDS should have Parent tag
        assert!(contigs[0].features[0].has_tag("Parent"));
    }

    #[test]
    fn test_centre_adds_protein_id() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        contigs[0].features = vec![make_cds(100, 400)];
        let config = ProkkaConfig {
            locustag: Some("T".into()),
            centre: Some("UoN".into()),
            ..Default::default()
        };
        assign_locus_tags(&mut contigs, &config);

        assert_eq!(
            contigs[0].features[0].get_tag("protein_id"),
            Some("gnl|UoN|T_00001")
        );
    }

    #[test]
    fn test_non_cds_rna_skipped() {
        let mut contigs = vec![Contig::new("ctg1".into(), b"ACGT".to_vec())];
        // Gene and RepeatRegion should NOT get locus tags
        contigs[0].features = vec![
            SeqFeature::new(FeatureType::Gene, "ctg1".into(), "t".into(), 100, 400, Strand::Forward),
            SeqFeature::new(FeatureType::RepeatRegion, "ctg1".into(), "t".into(), 500, 800, Strand::Forward),
        ];
        let config = ProkkaConfig {
            locustag: Some("T".into()),
            ..Default::default()
        };
        assign_locus_tags(&mut contigs, &config);

        // Neither should get locus_tag
        assert!(!contigs[0].features[0].has_tag("locus_tag"));
        assert!(!contigs[0].features[1].has_tag("locus_tag"));
    }
}
