use std::io::Write;
use std::path::Path;

/// Get current time as HH:MM:SS string (matching Perl Prokka's log format).
fn current_time_str() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    // Convert to local time approximation (UTC for now)
    let day_secs = (secs % 86400) as u32;
    let h = day_secs / 3600;
    let m = (day_secs % 3600) / 60;
    let s = day_secs % 60;
    format!("{:02}:{:02}:{:02}", h, m, s)
}

/// Convert days since Unix epoch to (year, month, day).
pub fn days_to_date(days: i64) -> (i64, u32, u32) {
    // Algorithm from http://howardhinnant.github.io/date_algorithms.html
    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u32;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d)
}

use crate::annotate::database::{AnnotationDb, DbFormat};
use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{AnnotationResult, AnnotationStats, Contig};
use crate::predict::overlap::filter_overlapping_cds;

/// Build the annotation database hierarchy in priority order.
///
/// Replicates Perl Prokka lines 940-1063.
fn build_annotation_databases(config: &ProkkaConfig) -> Vec<AnnotationDb> {
    let mut databases = Vec::new();
    let kingdom = config.kingdom.as_str();
    let dbdir = &config.dbdir;

    // User-supplied proteins (highest priority)
    if let Some(ref proteins) = config.proteins {
        databases.push(AnnotationDb {
            path: proteins.display().to_string(),
            source_prefix: if config.compliant {
                String::new()
            } else {
                let src = proteins.file_name()
                    .map(|f| f.to_string_lossy().to_string())
                    .unwrap_or_default();
                format!("similar to AA sequence:{}:", src)
            },
            format: DbFormat::Blast,
            evalue: None,
            min_coverage: None,
        });
    }

    // User-supplied HMMs
    if let Some(ref hmms) = config.hmms {
        let src = hmms.file_name()
            .map(|f| f.to_string_lossy().replace(".hmm", ""))
            .unwrap_or_default();
        databases.push(AnnotationDb {
            path: hmms.display().to_string(),
            source_prefix: if config.compliant {
                String::new()
            } else {
                format!("protein motif:{}:", src)
            },
            format: DbFormat::Hmmer3,
            evalue: None,
            min_coverage: None,
        });
    }

    // Genus-specific database
    if config.usegenus {
        let genus_db = dbdir.join("genus").join(&config.genus);
        if genus_db.exists() {
            databases.push(AnnotationDb {
                path: genus_db.display().to_string(),
                source_prefix: "similar to AA sequence:RefSeq:".to_string(),
                format: DbFormat::Blast,
                evalue: None,
                min_coverage: None,
            });
        }
    }

    // IS database (ISfinder)
    let is_db = dbdir.join("kingdom").join(kingdom).join("IS");
    if is_db.exists() {
        databases.push(AnnotationDb {
            path: is_db.display().to_string(),
            source_prefix: "similar to AA sequence:ISfinder:".to_string(),
            format: DbFormat::Blast,
            evalue: Some(1e-30),
            min_coverage: Some(90.0),
        });
    }

    // AMR database
    let amr_db = dbdir.join("kingdom").join(kingdom).join("AMR");
    if amr_db.exists() {
        databases.push(AnnotationDb {
            path: amr_db.display().to_string(),
            source_prefix: "similar to AA sequence:BARRGD:".to_string(),
            format: DbFormat::Blast,
            evalue: Some(1e-300),
            min_coverage: Some(90.0),
        });
    }

    // UniProtKB/Swiss-Prot
    let sprot_db = dbdir.join("kingdom").join(kingdom).join("sprot");
    if sprot_db.exists() {
        databases.push(AnnotationDb {
            path: sprot_db.display().to_string(),
            source_prefix: "similar to AA sequence:UniProtKB:".to_string(),
            format: DbFormat::Blast,
            evalue: None,
            min_coverage: None,
        });
    }

    // HMM databases (HAMAP, etc.)
    if kingdom != "Viruses" {
        let hmm_dir = dbdir.join("hmm");
        if hmm_dir.exists() {
            if let Ok(entries) = std::fs::read_dir(&hmm_dir) {
                for entry in entries.flatten() {
                    let path = entry.path();
                    if let Some(ext) = path.extension() {
                        if ext == "hmm" {
                            let name = path.file_stem()
                                .map(|s| s.to_string_lossy().to_string())
                                .unwrap_or_default();
                            databases.push(AnnotationDb {
                                path: path.display().to_string(),
                                source_prefix: format!("protein motif:{}:", name),
                                format: DbFormat::Hmmer3,
                                evalue: None,
                                min_coverage: None,
                            });
                        }
                    }
                }
            }
        }
    }

    databases
}

/// Run the full Prokka annotation pipeline.
///
/// This is the main library entry point. It takes contigs and returns
/// an AnnotationResult without requiring file I/O.
pub fn annotate(
    mut contigs: Vec<Contig>,
    config: &ProkkaConfig,
    fna_path: Option<&Path>,
) -> Result<AnnotationResult, ProkkaError> {
    let mut log = Vec::new();

    log.push(format!(
        "Annotating as >>> {} <<<",
        config.kingdom.as_str()
    ));
    log.push(format!("Loaded {} contigs", contigs.len()));

    let total_bp: usize = contigs.iter().map(|c| c.len()).sum();
    log.push(format!("Total {} bp", total_bp));

    let contig_lengths: Vec<(String, usize)> = contigs
        .iter()
        .map(|c| (c.id.clone(), c.len()))
        .collect();

    // ---- Phase 3: Feature prediction ----

    // 1. tRNA/tmRNA prediction (Aragorn)
    let mut trna_features = Vec::new();
    if !config.notrna {
        if let Some(fna) = fna_path {
            match crate::predict::trna::predict_trna(fna, config, &contig_lengths) {
                Ok(features) => {
                    log.push(format!("Found {} tRNAs", features.len()));
                    trna_features = features;
                }
                Err(ProkkaError::ToolNotFound { .. }) => {
                    log.push("Aragorn not found, skipping tRNA prediction.".into());
                }
                Err(e) => return Err(e),
            }
        }
    } else {
        log.push("Skipping tRNA search at user request.".into());
    }

    // 2. rRNA prediction (Barrnap)
    let mut rrna_features = Vec::new();
    if !config.norrna && config.kingdom.barrnap_mode().is_some() {
        if let Some(fna) = fna_path {
            match crate::predict::rrna::predict_rrna(fna, config) {
                Ok(features) => {
                    log.push(format!("Found {} rRNAs", features.len()));
                    rrna_features = features;
                }
                Err(ProkkaError::ToolNotFound { .. }) => {
                    log.push("Barrnap not found, skipping rRNA prediction.".into());
                }
                Err(e) => return Err(e),
            }
        }
    }

    // 3. ncRNA prediction (Infernal/cmscan, if --rfam)
    let mut ncrna_features = Vec::new();
    if config.rfam {
        if let Some(fna) = fna_path {
            match crate::predict::ncrna::predict_ncrna(fna, config, total_bp) {
                Ok(features) => {
                    log.push(format!("Found {} ncRNAs", features.len()));
                    ncrna_features = features;
                }
                Err(ProkkaError::ToolNotFound { .. }) => {
                    log.push("cmscan not found, skipping ncRNA prediction.".into());
                }
                Err(e) => return Err(e),
            }
        }
    }

    // 4. CRISPR detection (minced)
    let mut crispr_features = Vec::new();
    if let Some(fna) = fna_path {
        match crate::predict::crispr::predict_crispr(fna) {
            Ok(features) => {
                log.push(format!("Found {} CRISPRs", features.len()));
                crispr_features = features;
            }
            Err(ProkkaError::ToolNotFound { .. }) => {
                log.push("minced not found, skipping CRISPR detection.".into());
            }
            Err(e) => return Err(e),
        }
    }

    // Collect all RNA features for overlap exclusion
    let all_rna: Vec<_> = trna_features.iter()
        .chain(rrna_features.iter())
        .chain(crispr_features.iter())
        .cloned()
        .collect();

    log.push(format!(
        "Total of {} tRNA + rRNA features",
        trna_features.len() + rrna_features.len()
    ));

    // 5. CDS prediction (Prodigal)
    log.push("Predicting coding sequences".into());
    let mut cds_features = crate::predict::cds::predict_cds(&contigs, config)?;
    log.push(format!("Found {} CDS", cds_features.len()));

    // Exclude CDS overlapping RNA/CRISPR features
    let allow_overlap = config.cdsrnaolap || config.kingdom.allow_cds_rna_overlap();
    let pre_filter = cds_features.len();
    filter_overlapping_cds(&mut cds_features, &all_rna, allow_overlap);
    if pre_filter != cds_features.len() {
        log.push(format!(
            "Excluded {} CDS overlapping RNA features",
            pre_filter - cds_features.len()
        ));
    }

    // ---- Assign features to contigs (needed before SignalP) ----
    let mut all_features = Vec::new();
    all_features.extend(trna_features);
    all_features.extend(rrna_features);
    all_features.extend(ncrna_features);
    all_features.extend(crispr_features);
    all_features.extend(cds_features);

    for feature in all_features {
        if let Some(contig) = contigs.iter_mut().find(|c| c.id == feature.seq_id) {
            contig.features.push(feature);
        }
    }

    // Sort features by start position within each contig
    for contig in &mut contigs {
        contig.features.sort_by(|a, b| {
            a.start.cmp(&b.start)
                .then(b.end.cmp(&a.end))
        });
    }

    // 6. Signal peptide detection (SignalP, if --gram)
    if config.gram.is_some() {
        if let Some(fna) = fna_path {
            let outdir = fna.parent().unwrap_or(Path::new("."));
            match crate::predict::signalp::predict_signalp(outdir, &contigs, config) {
                Ok(features) => {
                    if !features.is_empty() {
                        log.push(format!("Found {} signal peptides", features.len()));
                        for f in features {
                            if let Some(contig) = contigs.iter_mut().find(|c| c.id == f.seq_id) {
                                contig.features.push(f);
                            }
                        }
                    }
                }
                Err(ProkkaError::ToolNotFound { .. }) => {
                    log.push("SignalP not found, skipping signal peptide detection.".into());
                }
                Err(e) => {
                    log.push(format!("SignalP error: {}", e));
                }
            }
        }
    }

    // ---- Phase 4: CDS Annotation ----
    if !config.noanno {
        log.push("Annotating CDS, please be patient.".into());

        // Build the annotation database hierarchy
        let databases = build_annotation_databases(config);
        let num_cds: usize = contigs.iter()
            .flat_map(|c| c.features.iter())
            .filter(|f| f.feature_type == crate::model::FeatureType::CDS)
            .count();

        for db in &databases {
            // Skip HMMs in --fast mode
            if config.fast && db.format == crate::annotate::database::DbFormat::Hmmer3 {
                log.push(format!("In --fast mode, skipping HMM search against {}", db.path));
                continue;
            }

            // Count remaining unannotated CDS
            let unannotated: usize = contigs.iter()
                .flat_map(|c| c.features.iter())
                .filter(|f| f.feature_type == crate::model::FeatureType::CDS && !f.has_tag("product"))
                .count();

            if unannotated == 0 {
                break;
            }

            log.push(format!(
                "There are still {} unannotated CDS left (started with {})",
                unannotated, num_cds
            ));

            let annotated = match db.format {
                crate::annotate::database::DbFormat::Blast => {
                    log.push(format!("Will use BLAST to search against {}", db.path));
                    crate::annotate::blast::annotate_blast(&mut contigs, db, config)?
                }
                crate::annotate::database::DbFormat::Hmmer3 => {
                    log.push(format!("Will use HMMER to search against {}", db.path));
                    crate::annotate::hmmer::annotate_hmmer(&mut contigs, db, config)?
                }
            };

            log.push(format!("Annotated {} CDS from {}", annotated, db.path));
        }

    }

    // ---- Phase 5: Post-processing ----

    // Product name cleanup (unless --rawproduct)
    if !config.rawproduct && !config.noanno {
        let good_products = crate::postprocess::product::default_good_products();
        let mut num_cleaned = 0;
        for contig in &mut contigs {
            for feature in &mut contig.features {
                if feature.feature_type != crate::model::FeatureType::CDS {
                    continue;
                }
                if let Some(product) = feature.get_tag("product").map(|s| s.to_string()) {
                    let cleaned = crate::postprocess::product::cleanup_product(&product, &good_products);
                    if cleaned != product {
                        feature.remove_tag("product");
                        feature.add_tag("product", &cleaned);
                        // If product became HYPO, clear gene/EC_number
                        if cleaned == crate::postprocess::product::HYPO {
                            feature.add_tag("note", &product);
                            feature.remove_tag("gene");
                            feature.remove_tag("EC_number");
                        }
                        num_cleaned += 1;
                    }
                }
            }
        }
        if num_cleaned > 0 {
            log.push(format!("Cleaned {} /product names", num_cleaned));
        }
    }

    // Detect possible pseudo genes (adjacent CDS with same annotation)
    for contig in &contigs {
        let mut prev_product = String::new();
        for feature in contig.features.iter().filter(|f| f.feature_type == crate::model::FeatureType::CDS) {
            let product = feature.get_tag("product").unwrap_or("").to_string();
            if product == prev_product
                && product != crate::postprocess::product::HYPO
                && product != "unannotated protein"
            {
                log.push(format!(
                    "Possible /pseudo '{}' at {} position {}",
                    product, feature.seq_id, feature.start
                ));
            }
            prev_product = product;
        }
    }

    // Gene deduplication
    crate::postprocess::gene_dedup::deduplicate_genes(&mut contigs);

    // Locus tag assignment
    crate::postprocess::locus_assign::assign_locus_tags(&mut contigs, config);

    // Label unannotated CDS (after locus_tag assignment, to match Perl tag order)
    let empty_label = if config.noanno { "unannotated protein" } else { crate::postprocess::product::HYPO };
    let mut num_hypo = 0;
    for contig in &mut contigs {
        for feature in &mut contig.features {
            if feature.feature_type == crate::model::FeatureType::CDS && !feature.has_tag("product") {
                feature.add_tag("product", empty_label);
                num_hypo += 1;
            }
        }
    }
    if num_hypo > 0 {
        log.push(format!("Labelled {} proteins as '{}'", num_hypo, empty_label));
    }

    // ---- Build statistics ----
    let mut stats = AnnotationStats {
        num_contigs: contigs.len(),
        total_bp,
        ..Default::default()
    };

    for contig in &contigs {
        for f in &contig.features {
            *stats.feature_counts.entry(f.feature_type.as_str().to_string()).or_insert(0) += 1;
        }
    }

    Ok(AnnotationResult {
        contigs,
        stats,
        log,
    })
}

/// Run the full pipeline from an input FASTA file, writing all outputs.
///
/// This is the high-level function used by the CLI.
pub fn run(input: &Path, config: &mut ProkkaConfig) -> Result<(), ProkkaError> {
    let start_time = std::time::Instant::now();

    config.apply_compliant();
    config.validate()?;

    // Incompatible option checks (matching Perl Prokka lines 267-270)
    if config.noanno && config.proteins.is_some() {
        return Err(ProkkaError::Other(
            "In --noanno mode, the --proteins option will not be used!".into(),
        ));
    }
    if config.noanno && config.hmms.is_some() {
        return Err(ProkkaError::Other(
            "In --noanno mode, the --hmms option will not be used!".into(),
        ));
    }
    if config.fast && config.hmms.is_some() {
        return Err(ProkkaError::Other(
            "In --fast mode, the --hmms option will not be used".into(),
        ));
    }

    // Welcome banner
    if !config.quiet {
        eprintln!("This is prokka-rs {}", env!("CARGO_PKG_VERSION"));
        eprintln!("Homepage is {}", env!("CARGO_PKG_HOMEPAGE"));
    }

    // Resolve locus_tag
    if config.locustag.is_none() {
        config.locustag = Some(crate::locus_tag::generate_locus_tag(input)?);
    }

    // Resolve output directory
    let prefix = config.prefix.clone().unwrap_or_else(|| {
        // Generate date-based prefix without chrono dependency.
        // Use std::time to get seconds since epoch, then compute MMDDYYYY.
        let secs = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs();
        // Simple date calculation (UTC, matching Perl Prokka's US date format)
        let days = (secs / 86400) as i64;
        let (y, m, d) = days_to_date(days);
        format!("PROKKA_{:02}{:02}{:04}", m, d, y)
    });
    config.prefix = Some(prefix.clone());

    let outdir = config.outdir.clone().unwrap_or_else(|| {
        std::path::PathBuf::from(&prefix)
    });

    // Create output directory
    if outdir.exists() {
        if !config.force {
            return Err(ProkkaError::OutputDirExists(outdir));
        }
    } else {
        std::fs::create_dir_all(&outdir)?;
    }

    // Load and sanitize input
    let contigs = crate::input::load_and_sanitize_fasta(input, config)?;

    let total_bp: usize = contigs.iter().map(|c| c.len()).sum();
    if !config.quiet {
        eprintln!(
            "Wrote {} contigs totalling {} bp.",
            contigs.len(),
            total_bp
        );
    }

    // Write .fna
    let fna_path = outdir.join(format!("{}.fna", prefix));
    {
        let mut fna_file = std::fs::File::create(&fna_path)?;
        crate::output::fasta::write_fna(&mut fna_file, &contigs)?;
    }

    // Run annotation pipeline
    let result = annotate(contigs, config, Some(&fna_path))?;

    // Print log with timestamps
    if !config.quiet {
        for line in &result.log {
            eprintln!("[{}] {}", current_time_str(), line);
        }
    }

    // Write output files
    let gcode = config.effective_gcode();

    // .gff
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.gff", prefix)))?;
        crate::output::gff::write_gff3(&mut f, &result, config.gffver)?;
    }
    // .tbl
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.tbl", prefix)))?;
        crate::output::tbl::write_tbl(&mut f, &result)?;
    }
    // .tsv
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.tsv", prefix)))?;
        crate::output::tsv::write_tsv(&mut f, &result)?;
    }
    // .faa (protein FASTA)
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.faa", prefix)))?;
        crate::output::fasta::write_faa(&mut f, &result, gcode)?;
    }
    // .ffn (nucleotide feature FASTA)
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.ffn", prefix)))?;
        crate::output::fasta::write_ffn(&mut f, &result)?;
    }
    // .fsa (annotated input FASTA for tbl2asn)
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.fsa", prefix)))?;
        crate::output::fasta::write_fsa(&mut f, &result, config)?;
    }
    // .txt (statistics)
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.txt", prefix)))?;
        crate::output::txt::write_txt(&mut f, &result, config)?;
    }
    // .log
    {
        let mut f = std::fs::File::create(outdir.join(format!("{}.log", prefix)))?;
        for line in &result.log {
            writeln!(f, "{}", line)?;
        }
    }

    // Generate GenBank (.gbk, .sqn) via tbl2asn or fallback
    match crate::output::genbank::generate_genbank(&outdir, &prefix, config) {
        Ok(()) => {}
        Err(ProkkaError::ToolNotFound { .. }) => {
            if !config.quiet {
                eprintln!("tbl2asn not found, .gbk/.sqn files not generated.");
            }
        }
        Err(e) => {
            if !config.quiet {
                eprintln!("Warning: GenBank generation failed: {}", e);
            }
        }
    }

    if !config.quiet {
        eprintln!("Annotation finished successfully.");
        let walltime = start_time.elapsed();
        eprintln!("Walltime used: {:.2} minutes", walltime.as_secs_f64() / 60.0);
        eprintln!("Output files:");
        if let Ok(entries) = std::fs::read_dir(&outdir) {
            let mut paths: Vec<_> = entries.flatten().map(|e| e.path()).collect();
            paths.sort();
            for path in paths {
                eprintln!("  {}", path.display());
            }
        }
    }

    Ok(())
}
