use std::io::BufReader;
use std::path::Path;

use crate::annotate::database::{parse_annotation_header, AnnotationDb};
use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::Contig;

/// Annotate unannotated CDS features via HMMER profile search.
///
/// Loads HMM profiles and searches each against the set of unannotated
/// CDS proteins. Maps to hmmscan semantics (best HMM per sequence).
///
/// Replicates Perl Prokka's HMMER annotation pipeline.
pub fn annotate_hmmer(
    contigs: &mut [Contig],
    db: &AnnotationDb,
    config: &ProkkaConfig,
) -> Result<usize, ProkkaError> {
    let db_path = Path::new(&db.path);
    if !db_path.exists() {
        return Ok(0);
    }

    let evalue_threshold = db.evalue.unwrap_or(config.evalue);
    let gcode = config.effective_gcode();

    // Load HMMs — try pressed format first, then ASCII
    let h3m_path = format!("{}.h3m", db.path);
    let hmms = if Path::new(&h3m_path).exists() {
        hmmer::io::binary::read_all_pressed(&h3m_path)
            .map_err(|e| ProkkaError::Hmmer(format!("Failed to read pressed HMMs: {}", e)))?
    } else {
        let file = std::fs::File::open(db_path)
            .map_err(|e| ProkkaError::Hmmer(format!("Cannot open HMM file: {}", e)))?;
        let mut reader = BufReader::new(file);
        hmmer::io::hmm_file::read_all_hmms(&mut reader)
            .map_err(|e| ProkkaError::Hmmer(format!("Failed to read HMMs: {}", e)))?
    };

    if hmms.is_empty() {
        return Ok(0);
    }

    // Collect unannotated CDS proteins
    let abc = hmmer::alphabet::Alphabet::amino();

    // Build digital sequences for all unannotated CDS
    struct CdsRef {
        contig_idx: usize,
        feature_idx: usize,
    }
    let mut digital_seqs = Vec::new();
    let mut cds_refs = Vec::new();

    for (ci, contig) in contigs.iter().enumerate() {
        for (fi, feature) in contig.features.iter().enumerate() {
            if feature.feature_type != crate::model::FeatureType::CDS {
                continue;
            }
            if feature.has_tag("product") {
                continue;
            }

            let dna = contig.extract(feature.start, feature.end, feature.strand);
            let protein = crate::annotate::blast::translate_dna(&dna, gcode);
            if protein.is_empty() {
                continue;
            }

            let seq_name = format!("cds_{}_{}", ci, fi);
            match hmmer::alphabet::DigitalSequence::new(&seq_name, "", &protein, &abc) {
                Ok(dseq) => {
                    digital_seqs.push(dseq);
                    cds_refs.push(CdsRef {
                        contig_idx: ci,
                        feature_idx: fi,
                    });
                }
                Err(_) => continue,
            }
        }
    }

    if digital_seqs.is_empty() {
        return Ok(0);
    }

    // Initialize logsum tables (required once for hmmer-pure-rs)
    hmmer::core::logsum::flogsum_init();

    let pipeline_config = hmmer::pipeline::PipelineConfig {
        e_value: evalue_threshold,
        noali: true,
        ..hmmer::pipeline::PipelineConfig::default()
    };

    // For hmmscan semantics: for each HMM, search against all sequences.
    // Track best HMM hit per sequence.
    struct BestHit {
        hmm_name: String,
        hmm_desc: String,
        evalue: f64,
    }

    let mut best_hits: Vec<Option<BestHit>> = (0..digital_seqs.len()).map(|_| None).collect();

    for hmm in &hmms {
        let (hits, _stats) = hmmer::pipeline::hmmsearch(hmm, &digital_seqs, &pipeline_config);

        for hit in hits {
            // Find which sequence this hit corresponds to
            if let Some(seq_idx) = digital_seqs.iter().position(|s| s.name == hit.name) {
                let is_better = match &best_hits[seq_idx] {
                    None => true,
                    Some(prev) => hit.evalue < prev.evalue,
                };
                if is_better {
                    best_hits[seq_idx] = Some(BestHit {
                        hmm_name: hmm.name.clone(),
                        hmm_desc: hmm.desc.clone().unwrap_or_default(),
                        evalue: hit.evalue,
                    });
                }
            }
        }
    }

    // Apply annotations
    let mut num_annotated = 0;

    for (seq_idx, best) in best_hits.into_iter().enumerate() {
        if let Some(hit) = best {
            let cref = &cds_refs[seq_idx];
            let feature = &mut contigs[cref.contig_idx].features[cref.feature_idx];

            let ann = parse_annotation_header(&hit.hmm_desc);

            let product = if ann.product.is_empty() {
                hit.hmm_desc.clone()
            } else {
                ann.product
            };

            feature.add_tag("product", &product);
            if !ann.ec_number.is_empty() {
                feature.add_tag("EC_number", &ann.ec_number);
            }
            if !ann.gene.is_empty() {
                feature.add_tag("gene", &ann.gene);
            }
            if !ann.cog.is_empty() {
                feature.add_tag("db_xref", &format!("COG:{}", ann.cog));
            }
            if !db.source_prefix.is_empty() {
                feature.add_tag(
                    "inference",
                    &format!("{}{}", db.source_prefix, hit.hmm_name),
                );
            }

            num_annotated += 1;
        }
    }

    Ok(num_annotated)
}
