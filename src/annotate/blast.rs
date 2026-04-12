use std::path::Path;

use crate::annotate::database::{parse_annotation_header, AnnotationDb};
use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::Contig;

/// Annotate unannotated CDS features via BLAST protein search.
///
/// Builds an indexed BlastDb from the FASTA database, then runs each
/// unannotated CDS query through `blastp()` for fast indexed search.
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

    // Load FASTA and build indexed BlastDb
    let db_file = std::fs::File::open(db_path)
        .map_err(|e| ProkkaError::Blast(format!("Cannot open database {}: {}", db.path, e)))?;
    let records = blast_rs::input::parse_fasta(db_file);

    if records.is_empty() {
        return Ok(0);
    }

    // Build indexed protein database in a temp directory
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

    // Open the indexed database
    let blast_db = blast_rs::BlastDb::open(&db_base)
        .map_err(|e| ProkkaError::Blast(format!("Cannot open BlastDb: {}", e)))?;

    // Configure search parameters
    let mut params = blast_rs::SearchParams::blastp();
    params.evalue_threshold = evalue_threshold;

    let mut num_annotated = 0;

    for contig in contigs.iter_mut() {
        // Collect unannotated CDS proteins (to avoid borrow conflicts)
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

        // Search each protein against indexed database
        for (fi, protein) in cds_data {
            let results = blast_rs::api::blastp(&blast_db, &protein, &params);

            // Find best hit passing coverage threshold
            let best = results.iter().find(|r| {
                r.hsps.iter().any(|hsp| {
                    let query_cov = (hsp.alignment_length as f64 / protein.len() as f64) * 100.0;
                    hsp.evalue <= evalue_threshold && query_cov >= coverage_threshold
                })
            });

            if let Some(hit) = best {
                // blast-rs subject_title may have binary header prefix bytes.
                // Strip non-ASCII-printable prefix and parse the defline.
                let clean_title = clean_blast_title(&hit.subject_title);
                let ann = parse_annotation_header(&clean_title);

                let product = if ann.product.is_empty() {
                    if hit.subject_title.contains("~~~") {
                        crate::postprocess::product::HYPO.to_string()
                    } else {
                        hit.subject_title.clone()
                    }
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
                    let accession = hit.subject_accession.trim_start_matches(|c: char| !c.is_ascii_alphanumeric());
                    feature.add_tag(
                        "inference",
                        &format!("{}{}", db.source_prefix, accession),
                    );
                }

                num_annotated += 1;
            }
        }
    }

    Ok(num_annotated)
}

/// Clean blast-rs subject_title which may have binary header prefix bytes.
///
/// The .phr format stores headers with ASN.1 encoding. blast-rs returns
/// raw bytes as a string, so we strip any leading non-printable or non-ASCII
/// characters to get the clean defline.
/// Clean blast-rs subject_title to extract the Prokka defline.
///
/// blast-rs returns raw .phr header bytes as a string, which may include
/// binary ASN.1 prefix and a duplicate accession. The original FASTA
/// defline format is: `>ACCESSION EC~~~gene~~~product~~~COG`
///
/// We need to extract just: `EC~~~gene~~~product~~~COG`
fn clean_blast_title(raw: &str) -> String {
    // Find the ~~~-delimited part. The raw title may be:
    //   "<binary>ACCESSION ACCESSION EC~~~gene~~~product~~~COG"
    // or for entries with empty EC:
    //   "<binary>ACCESSION ACCESSION ~~~gene~~~product~~~COG"
    if let Some(tilde_pos) = raw.find("~~~") {
        // Walk backwards from first ~~~ to find the space before the EC/~~~ field.
        // The EC field (or empty) is between the last space before ~~~ and ~~~.
        let before_tilde = &raw[..tilde_pos];
        if let Some(space_pos) = before_tilde.rfind(' ') {
            // Everything from after that space: "EC~~~gene~~~product~~~COG"
            return raw[space_pos + 1..].to_string();
        }
        // No space found — the whole thing starts with EC~~~
        return raw.to_string();
    }
    // No ~~~ — return as-is after stripping non-printable prefix
    let start = raw.find(|c: char| c.is_ascii_graphic()).unwrap_or(0);
    raw[start..].to_string()
}

/// Translate a DNA sequence to protein using the given genetic code.
///
/// Delegates to `crate::codon_table::translate_dna()` which implements
/// all NCBI genetic code tables (1-25).
pub fn translate_dna(dna: &[u8], gcode: u8) -> Vec<u8> {
    crate::codon_table::translate_dna(dna, gcode)
}

// Translation tests are in crate::codon_table::tests
