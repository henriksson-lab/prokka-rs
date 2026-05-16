//! Post-prediction filter that removes CDS overlapping a previously called
//! RNA or CRISPR feature.
//!
//! Mirrors the inline check in Perl Prokka (lines ~759-774): Prodigal will
//! occasionally call ORFs inside tRNA/rRNA/ncRNA/CRISPR regions, and those
//! CDS are excluded unless the user passes `--cdsrnaolap` (mitochondrial
//! genomes are packed, so overlaps are allowed there).

use crate::model::SeqFeature;

/// Drop CDS features that overlap any RNA or CRISPR feature.
///
/// Operates in place: any CDS whose interval intersects an entry in
/// `rna_features` (same contig) is removed. When `allow_overlap` is true
/// (corresponding to Perl Prokka's `--cdsrnaolap` flag), the input is
/// returned unchanged.
pub fn filter_overlapping_cds(
    cds_features: &mut Vec<SeqFeature>,
    rna_features: &[SeqFeature],
    allow_overlap: bool,
) {
    if allow_overlap {
        return;
    }
    cds_features.retain(|cds| {
        !rna_features.iter().any(|rna| cds.overlaps(rna))
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{FeatureType, Strand};

    fn make_feature(ftype: FeatureType, start: usize, end: usize) -> SeqFeature {
        SeqFeature::new(ftype, "ctg1".into(), "test".into(), start, end, Strand::Forward)
    }

    #[test]
    fn test_overlapping_cds_removed() {
        let mut cds = vec![
            make_feature(FeatureType::CDS, 100, 500),
            make_feature(FeatureType::CDS, 600, 900),
            make_feature(FeatureType::CDS, 1000, 1200),
        ];
        let rna = vec![
            make_feature(FeatureType::TRNA, 450, 550), // overlaps first CDS
        ];
        filter_overlapping_cds(&mut cds, &rna, false);
        assert_eq!(cds.len(), 2);
        assert_eq!(cds[0].start, 600);
        assert_eq!(cds[1].start, 1000);
    }

    #[test]
    fn test_no_overlap_keeps_all() {
        let mut cds = vec![
            make_feature(FeatureType::CDS, 100, 200),
            make_feature(FeatureType::CDS, 300, 400),
        ];
        let rna = vec![
            make_feature(FeatureType::RRNA, 500, 600),
        ];
        filter_overlapping_cds(&mut cds, &rna, false);
        assert_eq!(cds.len(), 2);
    }

    #[test]
    fn test_allow_overlap_keeps_all() {
        let mut cds = vec![
            make_feature(FeatureType::CDS, 100, 500),
        ];
        let rna = vec![
            make_feature(FeatureType::TRNA, 100, 500), // fully overlapping
        ];
        filter_overlapping_cds(&mut cds, &rna, true);
        assert_eq!(cds.len(), 1); // kept because allow_overlap=true
    }

    #[test]
    fn test_crispr_excludes_cds() {
        let mut cds = vec![
            make_feature(FeatureType::CDS, 11000, 12000),
        ];
        let rna = vec![
            make_feature(FeatureType::RepeatRegion, 11491, 12592),
        ];
        filter_overlapping_cds(&mut cds, &rna, false);
        assert_eq!(cds.len(), 0); // CDS overlaps CRISPR
    }

    #[test]
    fn test_different_contigs_no_overlap() {
        let mut cds = vec![
            SeqFeature::new(FeatureType::CDS, "ctg1".into(), "t".into(), 100, 500, Strand::Forward),
        ];
        let rna = vec![
            SeqFeature::new(FeatureType::TRNA, "ctg2".into(), "t".into(), 100, 500, Strand::Forward),
        ];
        filter_overlapping_cds(&mut cds, &rna, false);
        assert_eq!(cds.len(), 1); // different contigs, no overlap
    }
}
