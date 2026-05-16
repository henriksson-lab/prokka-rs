//! HMMER profile annotation of CDS features.
//!
//! Loads HMM profiles from an annotation database and searches them against
//! the translations of any CDS still missing a `/product` tag, recording the
//! best (lowest e-value) hit per sequence. Hit descriptions are parsed as
//! `~~~`-delimited Prokka deflines.
//!
//! Replaces Perl Prokka's external `hmmscan` invocation (`$HMMER3CMD` at
//! line 73 / annotation loop near line 1075).

use std::path::Path;

use crate::annotate::database::{parse_annotation_header, AnnotationDb};
use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::Contig;
use hmmer_pure_rs::alphabet::Alphabet;
use hmmer_pure_rs::bg::Bg;
use hmmer_pure_rs::hmmfile;
use hmmer_pure_rs::logsum;
use hmmer_pure_rs::pipeline::Pipeline;
use hmmer_pure_rs::profile::{self, Profile, P7_LOCAL};
use hmmer_pure_rs::sequence::Sequence;
use hmmer_pure_rs::simd::oprofile::OProfile;
use hmmer_pure_rs::tophits::{TopHits, P7_IS_REPORTED};

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

    // newhmmer's public search API reads the source HMM file and builds the
    // required profiles in-process. Prokka only schedules HMM DBs that have
    // pressed sidecars, but the ASCII .hmm path remains the source of truth.
    let hmms = hmmfile::read_hmm_file(db_path)
        .map_err(|e| ProkkaError::Hmmer(format!("Failed to read HMMs: {}", e)))?;

    if hmms.is_empty() {
        return Ok(0);
    }

    // Collect unannotated CDS proteins
    let abc = Alphabet::amino();
    let bg = Bg::new(&abc);

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

            let dsq = abc.digitize(&protein);
            let n = dsq.len().saturating_sub(2);
            if n == 0 {
                continue;
            }

            let seq_name = format!("cds_{}_{}", ci, fi);
            digital_seqs.push(Sequence {
                name: seq_name,
                acc: String::new(),
                desc: String::new(),
                dsq,
                n,
                l: n,
            });
            cds_refs.push(CdsRef {
                contig_idx: ci,
                feature_idx: fi,
            });
        }
    }

    if digital_seqs.is_empty() {
        return Ok(0);
    }

    logsum::p7_flogsuminit();
    hmmer_pure_rs::util::simd_env::init();

    // For hmmscan semantics: for each HMM, search against all sequences.
    // Track best HMM hit per sequence.
    struct BestHit {
        hmm_name: String,
        hmm_desc: String,
        evalue: f64,
    }

    let mut best_hits: Vec<Option<BestHit>> = (0..digital_seqs.len()).map(|_| None).collect();

    let z = hmms.len() as f64;
    for hmm in &hmms {
        let mut local_bg = bg.clone();
        local_bg.set_filter(hmm.m, &hmm.compo);

        for (seq_idx, seq) in digital_seqs.iter().enumerate() {
            local_bg.set_length(seq.n);

            let mut gm = Profile::new(hmm.m, &abc);
            profile::profile_config(hmm, &local_bg, &mut gm, seq.n as i32, P7_LOCAL);
            let mut om = OProfile::convert(&gm);

            let mut pli = Pipeline::new();
            pli.e_value_threshold = evalue_threshold;
            pli.inc_e = evalue_threshold;
            pli.do_alignment = false;
            pli.do_alignment_display = false;
            pli.new_model(&gm);

            let mut th = TopHits::new();
            if !pli.run(&mut gm, &mut om, &local_bg, hmm, seq, &mut th) {
                continue;
            }

            th.threshold(&pli, z, z);
            let Some(hit) = th
                .hits
                .iter()
                .filter(|hit| hit.flags & P7_IS_REPORTED != 0)
                .min_by(|a, b| {
                    let a_eval = z * a.lnp.exp();
                    let b_eval = z * b.lnp.exp();
                    a_eval.total_cmp(&b_eval)
                })
            else {
                continue;
            };

            let evalue = z * hit.lnp.exp();
            if evalue > evalue_threshold {
                continue;
            }

            let is_better = match &best_hits[seq_idx] {
                None => true,
                Some(prev) => evalue < prev.evalue,
            };
            if is_better {
                best_hits[seq_idx] = Some(BestHit {
                    hmm_name: hmm.name.clone(),
                    hmm_desc: hmm.desc.clone().unwrap_or_default(),
                    evalue,
                });
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
                if hit.hmm_desc.contains("~~~") {
                    crate::postprocess::product::HYPO.to_string()
                } else {
                    hit.hmm_desc.clone()
                }
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
