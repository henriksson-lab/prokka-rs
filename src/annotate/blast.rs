use std::path::Path;

use crate::annotate::database::{parse_annotation_header, AnnotationDb};
use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::Contig;

/// Annotate unannotated CDS features via BLAST protein search.
///
/// Matches Perl Prokka's approach: all queries searched against a pre-indexed
/// database with `-num_descriptions 1 -num_alignments 1 -seg no`.
/// Parallelizes across queries via a shared rayon thread pool.
pub fn annotate_blast(
    contigs: &mut [Contig],
    db: &AnnotationDb,
    config: &ProkkaConfig,
    pool: Option<&rayon::ThreadPool>,
) -> Result<usize, ProkkaError> {
    let db_path = Path::new(&db.path);
    if !db_path.exists() {
        return Ok(0);
    }

    let evalue_threshold = db.evalue.unwrap_or(config.evalue);
    let coverage_threshold = db.min_coverage.unwrap_or(config.coverage);
    let gcode = config.effective_gcode();

    // Load FASTA and build indexed BlastDb
    let db_file = std::fs::File::open(db_path)
        .map_err(|e| ProkkaError::Blast(format!("Cannot open database {}: {}", db.path, e)))?;
    let records = blast_rs::input::parse_fasta(db_file);
    if records.is_empty() {
        return Ok(0);
    }

    let tmp_dir = tempfile::tempdir()
        .map_err(|e| ProkkaError::Blast(format!("Cannot create temp dir: {}", e)))?;
    let db_base = tmp_dir.path().join("blastdb");

    let mut builder = blast_rs::BlastDbBuilder::new(
        blast_rs::db::DbType::Protein,
        "prokka_annotation_db",
    );
    for rec in &records {
        builder.add(blast_rs::SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&db_base)
        .map_err(|e| ProkkaError::Blast(format!("Cannot build BlastDb: {}", e)))?;

    let blast_db = blast_rs::BlastDb::open(&db_base)
        .map_err(|e| ProkkaError::Blast(format!("Cannot open BlastDb: {}", e)))?;

    // Collect all unannotated CDS queries
    struct CdsRef {
        contig_idx: usize,
        feature_idx: usize,
    }
    let mut refs: Vec<CdsRef> = Vec::new();
    let mut proteins: Vec<Vec<u8>> = Vec::new();

    for (ci, contig) in contigs.iter().enumerate() {
        for (fi, feature) in contig.features.iter().enumerate() {
            if feature.feature_type != crate::model::FeatureType::CDS {
                continue;
            }
            if feature.has_tag("product") {
                continue;
            }
            let dna = contig.extract(feature.start, feature.end, feature.strand);
            let protein = translate_dna(&dna, gcode);
            if protein.is_empty() {
                continue;
            }
            refs.push(CdsRef { contig_idx: ci, feature_idx: fi });
            proteins.push(protein);
        }
    }

    if proteins.is_empty() {
        return Ok(0);
    }

    // Per-query blastp with parallel dispatch across queries.
    // Each blastp uses two-hit algorithm + num_threads=1.
    let params = blast_rs::SearchParams::blastp()
        .evalue(evalue_threshold)
        .max_target_seqs(1)
        .filter_low_complexity(false)
        .num_threads(1);

    let search_one = |protein: &[u8]| -> Option<(String, String)> {
        let results = blast_rs::api::blastp(&blast_db, protein, &params);
        results.iter().find(|r| {
            r.hsps.iter().any(|hsp| {
                let query_cov = (hsp.alignment_length as f64 / protein.len() as f64) * 100.0;
                hsp.evalue <= evalue_threshold && query_cov >= coverage_threshold
            })
        }).map(|hit| (hit.subject_accession.clone(), hit.subject_title.clone()))
    };

    let search_results: Vec<Option<(usize, String, String)>> = match pool {
        Some(pool) => {
            use rayon::prelude::*;
            pool.install(|| {
                proteins.par_iter().enumerate()
                    .map(|(qi, p)| search_one(p).map(|(acc, title)| (qi, acc, title)))
                    .collect()
            })
        }
        None => {
            proteins.iter().enumerate()
                .map(|(qi, p)| search_one(p).map(|(acc, title)| (qi, acc, title)))
                .collect()
        }
    };

    // Apply annotations
    let mut num_annotated = 0;
    for result in search_results.into_iter().flatten() {
        let (qi, hit_accession, hit_title) = result;

        let clean_title = clean_blast_title(&hit_title);
        let ann = parse_annotation_header(&clean_title);

        let product = if ann.product.is_empty() {
            if hit_title.contains("~~~") {
                crate::postprocess::product::HYPO.to_string()
            } else {
                hit_title.clone()
            }
        } else {
            ann.product
        };

        let r = &refs[qi];
        let feature = &mut contigs[r.contig_idx].features[r.feature_idx];
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
            let accession = hit_accession.trim_start_matches(|c: char| !c.is_ascii_alphanumeric());
            feature.add_tag(
                "inference",
                &format!("{}{}", db.source_prefix, accession),
            );
        }
        num_annotated += 1;
    }

    Ok(num_annotated)
}

/// Clean blast-rs subject_title to extract the Prokka defline.
fn clean_blast_title(raw: &str) -> String {
    if let Some(tilde_pos) = raw.find("~~~") {
        let before_tilde = &raw[..tilde_pos];
        if let Some(space_pos) = before_tilde.rfind(' ') {
            return raw[space_pos + 1..].to_string();
        }
        return raw.to_string();
    }
    let start = raw.find(|c: char| c.is_ascii_graphic()).unwrap_or(0);
    raw[start..].to_string()
}

/// Translate a DNA sequence to protein using the given genetic code.
pub fn translate_dna(dna: &[u8], gcode: u8) -> Vec<u8> {
    crate::codon_table::translate_dna(dna, gcode)
}
