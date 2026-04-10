use std::path::Path;

use crate::annotate::database::{parse_annotation_header, AnnotationDb};
use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::Contig;

use blast_rs::matrix::AA_SIZE;

/// Create a simplified protein scoring matrix compatible with NCBIstdaa encoding.
/// +5 on diagonal for standard amino acids (indices 1-20), -2 elsewhere.
/// This matches the approach used in blast-rs main.rs.
fn make_simple_protein_matrix() -> [[i32; AA_SIZE]; AA_SIZE] {
    let mut m = [[-2i32; AA_SIZE]; AA_SIZE];
    for (i, row) in m.iter_mut().enumerate().take(21).skip(1) {
        row[i] = 5;
    }
    m
}

/// Annotate unannotated CDS features via BLAST protein search.
///
/// For each unannotated CDS, translates to protein and searches against
/// the BLAST database using blast-rs. Only CDS without a 'product' tag
/// are processed.
///
/// Replicates Perl Prokka lines 1075-1161.
pub fn annotate_blast(
    contigs: &mut [Contig],
    db: &AnnotationDb,
    config: &ProkkaConfig,
) -> Result<usize, ProkkaError> {
    let db_path = Path::new(&db.path);
    if !db_path.exists() {
        return Ok(0);
    }

    let evalue_threshold = db.evalue.unwrap_or(config.evalue);
    let coverage_threshold = db.min_coverage.unwrap_or(config.coverage);
    let gcode = config.effective_gcode();

    // Load database subjects from FASTA
    let db_file = std::fs::File::open(db_path)
        .map_err(|e| ProkkaError::Blast(format!("Cannot open database {}: {}", db.path, e)))?;
    let subjects = blast_rs::input::parse_fasta(db_file);

    if subjects.is_empty() {
        return Ok(0);
    }

    // Pre-encode all subjects
    let encoded_subjects: Vec<(String, String, Vec<u8>)> = subjects
        .iter()
        .map(|s| {
            let encoded: Vec<u8> = s
                .sequence
                .iter()
                .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
                .collect();
            (s.id.clone(), s.defline.clone(), encoded)
        })
        .collect();

    // Calculate total database length for E-value computation
    let db_total_len: i64 = encoded_subjects.iter().map(|(_, _, s)| s.len() as i64).sum();
    let num_db_seqs = encoded_subjects.len() as i32;

    // Karlin-Altschul statistics for BLOSUM62 with gap_open=11, gap_extend=1
    let kbp = match blast_rs::stat::lookup_protein_params(11, 1) {
        Some(p) => blast_rs::stat::KarlinBlk {
            lambda: p.lambda,
            k: p.k,
            log_k: p.k.ln(),
            h: p.h,
        },
        None => blast_rs::stat::protein_ungapped_kbp(),
    };

    // NOTE: blast_rs::matrix::BLOSUM62 has incorrect indexing for NCBIstdaa encoding.
    // Use simplified identity-like matrix matching blast-rs main.rs behavior:
    // +5 on diagonal for standard amino acids (indices 1-20), -2 elsewhere.
    let matrix = make_simple_protein_matrix();

    let mut num_annotated = 0;

    // First pass: collect unannotated CDS info (to avoid borrow conflicts)
    for contig in contigs.iter_mut() {
        // Collect indices and protein translations of unannotated CDS
        let mut cds_data: Vec<(usize, Vec<u8>)> = Vec::new();

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

            cds_data.push((fi, protein));
        }

        // Second pass: search and annotate
        for (fi, protein) in cds_data {
            let query_encoded: Vec<u8> = protein
                .iter()
                .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
                .collect();

            let len_adj = blast_rs::stat::compute_length_adjustment(
                query_encoded.len() as i32,
                db_total_len,
                num_db_seqs,
                &kbp,
            );
            let search_space = blast_rs::stat::compute_search_space(
                query_encoded.len() as i64,
                db_total_len,
                num_db_seqs,
                len_adj,
            );

            let mut best_hit: Option<(f64, String, String)> = None;

            for (subj_id, subj_desc, subj_encoded) in &encoded_subjects {
                let hits = blast_rs::protein_lookup::protein_gapped_scan(
                    &query_encoded,
                    subj_encoded,
                    &matrix,
                    3,    // word_size
                    11.0, // threshold
                    30,   // ungap_x_dropoff
                    11,   // gap_open
                    1,    // gap_extend
                    50,   // gap_x_dropoff
                    0,    // ungap_cutoff
                );

                for hit in hits {
                    let query_cov =
                        (hit.align_length as f64 / query_encoded.len() as f64) * 100.0;
                    if query_cov < coverage_threshold {
                        continue;
                    }

                    let eval = kbp.raw_to_evalue(hit.score, search_space);
                    if eval > evalue_threshold {
                        continue;
                    }

                    if best_hit.is_none() || eval < best_hit.as_ref().unwrap().0 {
                        best_hit = Some((eval, subj_id.clone(), subj_desc.clone()));
                    }
                }
            }

            // Apply best hit annotation
            if let Some((_eval, hit_id, hit_desc)) = best_hit {
                let ann = parse_annotation_header(&hit_desc);

                let product = if ann.product.is_empty() {
                    hit_desc.clone()
                } else {
                    ann.product
                };

                let feature = &mut contig.features[fi];
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
                        &format!("{}{}", db.source_prefix, hit_id),
                    );
                }

                num_annotated += 1;
            }
        }
    }

    Ok(num_annotated)
}

/// Translate a DNA sequence to protein using the given genetic code.
///
/// Delegates to `crate::codon_table::translate_dna()` which implements
/// all NCBI genetic code tables (1-25).
pub fn translate_dna(dna: &[u8], gcode: u8) -> Vec<u8> {
    crate::codon_table::translate_dna(dna, gcode)
}

// Translation tests are in crate::codon_table::tests
