//! CDS (coding sequence) prediction via the native `prodigal-rs` crate.
//!
//! Implements pipeline step 6: replaces the Perl call to the `prodigal`
//! binary (Perl Prokka lines ~721-786) with an in-process port of
//! Prodigal v2.6.3. Decides between `single` (trained) and `meta` mode
//! the same way the Perl pipeline does.

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{Contig, FeatureType, SeqFeature, Strand};

/// Predict CDS features using prodigal-rs.
///
/// Replicates Perl Prokka lines 721-786:
/// - Mode selection: `single` if total_bp >= 100000 and not metagenome, else `meta`
/// - In single mode: train on all contigs (concatenated), then predict per contig
/// - In meta mode: predict each contig independently
/// - Flags: closed_ends=true, mask_n_runs=true
pub fn predict_cds(
    contigs: &[Contig],
    config: &ProkkaConfig,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    let gcode = config.effective_gcode();
    let total_bp: usize = contigs.iter().map(|c| c.len()).sum();

    let prodigal_config = prodigal_rs::ProdigalConfig {
        translation_table: gcode,
        closed_ends: true,
        mask_n_runs: true,
        force_non_sd: false,
    };

    // Mode selection: single if >= 100000bp and not metagenome, else meta
    // If a training file is provided, always use single mode with that file
    let use_single = config.prodigaltf.is_some()
        || (total_bp >= 100_000 && !config.metagenome);

    let mut all_features = Vec::new();

    if use_single {
        // Single mode: train on all contigs (or load training file), then predict per contig
        let training = if let Some(ref tf_path) = config.prodigaltf {
            // Load user-supplied Prodigal training file
            prodigal_rs::TrainingData::load(tf_path)
                .map_err(|e| ProkkaError::Prodigal(format!("Cannot load training file: {}", e)))?
        } else {
            // Train on concatenated contigs with stop-codon spacers
            let spacer = b"TTAATTAATTAA";
            let mut concat = Vec::with_capacity(total_bp + contigs.len() * spacer.len());
            for (i, contig) in contigs.iter().enumerate() {
                if i > 0 {
                    concat.extend_from_slice(spacer);
                }
                concat.extend_from_slice(&contig.sequence);
            }
            prodigal_rs::train_with(&concat, &prodigal_config)
                .map_err(|e| ProkkaError::Prodigal(e.to_string()))?
        };

        for contig in contigs {
            if contig.is_empty() {
                continue;
            }
            match prodigal_rs::predict_with(&contig.sequence, &training, &prodigal_config) {
                Ok(genes) => {
                    for gene in genes {
                        all_features.push(prodigal_gene_to_feature(gene, &contig.id, gcode));
                    }
                }
                Err(e) => {
                    // Skip contigs that are too short for prediction
                    if !config.quiet {
                        eprintln!("Warning: Prodigal skipped contig {}: {}", contig.id, e);
                    }
                }
            }
        }
    } else {
        // Meta mode: predict each contig independently
        for contig in contigs {
            if contig.is_empty() {
                continue;
            }
            match prodigal_rs::predict_meta_with(&contig.sequence, &prodigal_config) {
                Ok(genes) => {
                    for gene in genes {
                        all_features.push(prodigal_gene_to_feature(gene, &contig.id, gcode));
                    }
                }
                Err(e) => {
                    if !config.quiet {
                        eprintln!("Warning: Prodigal skipped contig {}: {}", contig.id, e);
                    }
                }
            }
        }
    }

    Ok(all_features)
}

/// Convert a `prodigal_rs::PredictedGene` into a Prokka [`SeqFeature`].
///
/// Sets `feature.source` to `Prodigal:prodigal-rs-<crate-version>` and
/// adds the standard `inference=ab initio prediction:<source>` tag, mirroring
/// the source-tag format used by Perl Prokka.
fn prodigal_gene_to_feature(
    gene: prodigal_rs::PredictedGene,
    seq_id: &str,
    _gcode: u8,
) -> SeqFeature {
    let strand = match gene.strand {
        prodigal_rs::Strand::Forward => Strand::Forward,
        prodigal_rs::Strand::Reverse => Strand::Reverse,
    };

    let version = env!("CARGO_PKG_VERSION");
    let tool = format!("Prodigal:prodigal-rs-{}", version);

    let mut feature = SeqFeature::new(
        FeatureType::CDS,
        seq_id.to_string(),
        tool.clone(),
        gene.begin,
        gene.end,
        strand,
    );

    feature.add_tag("inference", &format!("ab initio prediction:{}", tool));

    feature
}

/// Decide whether Prodigal should run in `single` (trained) or `meta` mode.
///
/// Matches Perl Prokka line 727:
/// `single` if a user training file is supplied, or if total bp >= 100 000
/// and `--metagenome` was not set; otherwise `meta`. Returns `true` for
/// single mode and `false` for meta mode.
pub fn select_mode(config: &ProkkaConfig, total_bp: usize) -> bool {
    config.prodigaltf.is_some() || (total_bp >= 100_000 && !config.metagenome)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mode_single_for_large_genome() {
        let config = ProkkaConfig::default(); // metagenome=false, prodigaltf=None
        assert!(select_mode(&config, 200_000));  // large -> single
        assert!(!select_mode(&config, 50_000));   // small -> meta
        assert!(!select_mode(&config, 99_999));   // just under threshold -> meta
        assert!(select_mode(&config, 100_000));   // exactly at threshold -> single
    }

    #[test]
    fn test_mode_meta_for_metagenome() {
        let config = ProkkaConfig {
            metagenome: true,
            ..Default::default()
        };
        assert!(!select_mode(&config, 10_000_000)); // metagenome overrides size
    }

    #[test]
    fn test_mode_single_with_training_file() {
        let config = ProkkaConfig {
            prodigaltf: Some(std::path::PathBuf::from("/tmp/fake.trn")),
            metagenome: true, // even with metagenome, training file forces single
            ..Default::default()
        };
        assert!(select_mode(&config, 1000)); // training file -> always single
    }

    #[test]
    fn test_predict_cds_on_small_sequence() {
        // A very small sequence should still work (meta mode)
        let contigs = vec![Contig::new(
            "test".into(),
            b"ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG".to_vec(),
        )];
        let config = ProkkaConfig {
            quiet: true,
            ..Default::default()
        };
        // This may return 0 genes for such a tiny sequence, but should not panic
        let result = predict_cds(&contigs, &config);
        assert!(result.is_ok());
    }
}
