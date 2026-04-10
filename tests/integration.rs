use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("prokka").join("test")
}

#[allow(dead_code)]
fn db_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("prokka").join("db")
}

/// Test basic pipeline with --noanno (no database annotation, fast).
#[test]
fn test_plasmid_noanno() {
    let dir = tempfile::tempdir().unwrap();
    let outdir = dir.path().join("output");

    let mut config = prokka_rs::ProkkaConfig {
        outdir: Some(outdir.clone()),
        prefix: Some("plasmid".into()),
        force: true,
        notrna: true,
        norrna: true,
        noanno: true,
        quiet: true,
        ..Default::default()
    };

    let input = test_data_dir().join("plasmid.fna");
    prokka_rs::pipeline::run(&input, &mut config).unwrap();

    // Check output files exist
    assert!(outdir.join("plasmid.fna").exists(), "Missing .fna");
    assert!(outdir.join("plasmid.gff").exists(), "Missing .gff");
    assert!(outdir.join("plasmid.tbl").exists(), "Missing .tbl");
    assert!(outdir.join("plasmid.tsv").exists(), "Missing .tsv");
    assert!(outdir.join("plasmid.faa").exists(), "Missing .faa");
    assert!(outdir.join("plasmid.ffn").exists(), "Missing .ffn");
    assert!(outdir.join("plasmid.fsa").exists(), "Missing .fsa");
    assert!(outdir.join("plasmid.txt").exists(), "Missing .txt");
    assert!(outdir.join("plasmid.log").exists(), "Missing .log");

    // Check .fna sequence content matches input
    let input_seq: String = std::fs::read_to_string(&input)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.chars())
        .collect();
    let output_seq: String = std::fs::read_to_string(outdir.join("plasmid.fna"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.chars())
        .collect();
    assert_eq!(input_seq, output_seq, "Sequences differ");

    // Check .txt has expected stats
    let txt = std::fs::read_to_string(outdir.join("plasmid.txt")).unwrap();
    assert!(txt.contains("contigs: 1"), "Wrong contig count");
    assert!(txt.contains("bases: 56520"), "Wrong base count");
    assert!(txt.contains("CDS:"), "Missing CDS count");

    // Check .tsv has header and data
    let tsv = std::fs::read_to_string(outdir.join("plasmid.tsv")).unwrap();
    let lines: Vec<&str> = tsv.lines().collect();
    assert!(lines[0].contains("locus_tag\tftype"), "Bad TSV header");
    assert!(lines.len() > 50, "Too few TSV lines: {}", lines.len());

    // Check .gff is valid GFF3
    let gff = std::fs::read_to_string(outdir.join("plasmid.gff")).unwrap();
    assert!(gff.starts_with("##gff-version 3"), "Not GFF3");
    assert!(gff.contains("##FASTA"), "Missing FASTA section");
    assert!(gff.contains("CDS"), "No CDS features in GFF");

    // In noanno mode, all CDS should be "unannotated protein"
    for line in tsv.lines().skip(1) {
        if line.contains("\tCDS\t") {
            assert!(
                line.contains("unannotated protein"),
                "CDS not labeled as unannotated: {}",
                line
            );
        }
    }
}

/// Test that CDS prediction finds a reasonable number of genes.
#[test]
fn test_plasmid_cds_count() {
    let dir = tempfile::tempdir().unwrap();
    let outdir = dir.path().join("output");

    let mut config = prokka_rs::ProkkaConfig {
        outdir: Some(outdir.clone()),
        prefix: Some("plasmid".into()),
        force: true,
        notrna: true,
        norrna: true,
        noanno: true,
        quiet: true,
        ..Default::default()
    };

    let input = test_data_dir().join("plasmid.fna");
    prokka_rs::pipeline::run(&input, &mut config).unwrap();

    let txt = std::fs::read_to_string(outdir.join("plasmid.txt")).unwrap();
    // Extract CDS count
    let cds_line = txt.lines().find(|l| l.starts_with("CDS:")).unwrap();
    let cds_count: usize = cds_line.split(':').nth(1).unwrap().trim().parse().unwrap();

    // Prokka typically finds ~60-65 CDS on this plasmid
    assert!(
        cds_count >= 55 && cds_count <= 75,
        "Unexpected CDS count: {} (expected 55-75)",
        cds_count
    );
}

/// Test locus_tag generation is deterministic.
#[test]
fn test_locus_tag_deterministic() {
    let input = test_data_dir().join("plasmid.fna");
    let tag1 = prokka_rs::locus_tag::generate_locus_tag(&input).unwrap();
    let tag2 = prokka_rs::locus_tag::generate_locus_tag(&input).unwrap();
    assert_eq!(tag1, tag2);
    assert_eq!(tag1.len(), 8);
    for c in tag1.chars() {
        assert!(c.is_ascii_uppercase(), "Non-uppercase char in tag: {}", c);
    }
}

/// Test input sanitization.
#[test]
fn test_input_sanitization() {
    let dir = tempfile::tempdir().unwrap();

    // Create a test FASTA with tricky content
    let fasta_path = dir.path().join("test.fna");
    std::fs::write(
        &fasta_path,
        ">seq1|pipe\nACGTacgtRYSW\n>seq2\nACGT*-NNNN\n",
    )
    .unwrap();

    let config = prokka_rs::ProkkaConfig {
        mincontiglen: 1,
        quiet: true,
        ..Default::default()
    };
    let contigs = prokka_rs::input::load_and_sanitize_fasta(&fasta_path, &config).unwrap();

    assert_eq!(contigs.len(), 2);
    // Pipe replaced with underscore
    assert_eq!(contigs[0].id, "seq1_pipe");
    // Uppercase, RYSW -> NNNN
    assert_eq!(contigs[0].sequence, b"ACGTACGTNNNN");
    // Stars and dashes removed, N kept
    assert_eq!(contigs[1].sequence, b"ACGTNNNN");
}

/// Test that protein translation works correctly.
#[test]
fn test_protein_translation() {
    use prokka_rs::annotate::blast::translate_dna;

    // Standard translation
    assert_eq!(translate_dna(b"ATGGCTAAA", 11), b"MAK");
    // GTG start codon should become M
    assert_eq!(translate_dna(b"GTGGCTAAA", 11), b"MAK");
    // Stop codon stripped
    assert_eq!(translate_dna(b"ATGGCTAAATAA", 11), b"MAK");
    // Too short
    assert_eq!(translate_dna(b"AT", 11), b"");
}

/// Test the annotation header parser.
#[test]
fn test_annotation_parser() {
    use prokka_rs::annotate::database::parse_annotation_header;

    let ann = parse_annotation_header("2.1.1.48~~~ermC~~~rRNA adenine N-6-methyltransferase~~~COG1234");
    assert_eq!(ann.ec_number, "2.1.1.48");
    assert_eq!(ann.gene, "ermC");
    assert_eq!(ann.product, "rRNA adenine N-6-methyltransferase");
    assert_eq!(ann.cog, "COG1234");

    // Plain product (no ~~~)
    let ann = parse_annotation_header("hypothetical protein");
    assert_eq!(ann.product, "hypothetical protein");
    assert!(ann.gene.is_empty());
}

/// Test with larger genome input (genome.fna, 6.7MB).
#[test]
fn test_genome_noanno() {
    let dir = tempfile::tempdir().unwrap();
    let outdir = dir.path().join("output");

    let mut config = prokka_rs::ProkkaConfig {
        outdir: Some(outdir.clone()),
        prefix: Some("genome".into()),
        force: true,
        notrna: true,
        norrna: true,
        noanno: true,
        quiet: true,
        ..Default::default()
    };

    let input = test_data_dir().join("genome.fna");
    if !input.exists() {
        return; // skip if test data not present
    }
    prokka_rs::pipeline::run(&input, &mut config).unwrap();

    let txt = std::fs::read_to_string(outdir.join("genome.txt")).unwrap();
    assert!(txt.contains("contigs: 1"));
    assert!(txt.contains("bases: 6601757"));

    let cds_line = txt.lines().find(|l| l.starts_with("CDS:")).unwrap();
    let cds_count: usize = cds_line.split(':').nth(1).unwrap().trim().parse().unwrap();
    // Genome should have thousands of CDS
    assert!(
        cds_count >= 5000 && cds_count <= 7000,
        "Unexpected CDS count: {} (expected 5000-7000)",
        cds_count
    );
}

/// Test the library API (annotate function directly).
#[test]
fn test_library_api() {
    let config = prokka_rs::ProkkaConfig {
        noanno: true,
        notrna: true,
        norrna: true,
        quiet: true,
        ..Default::default()
    };

    let contigs = vec![prokka_rs::Contig::new(
        "test_seq".to_string(),
        b"ATGAAACCCGGGTTTTAA".to_vec(),
    )];

    let result = prokka_rs::annotate(contigs, &config, None).unwrap();
    assert_eq!(result.contigs.len(), 1);
    assert_eq!(result.stats.num_contigs, 1);
    assert_eq!(result.stats.total_bp, 18);
}

/// Test product name cleanup.
#[test]
fn test_product_cleanup() {
    use prokka_rs::postprocess::product::{cleanup_product, default_good_products, HYPO};

    let good = default_good_products();

    assert_eq!(cleanup_product("DUF1234 domain protein", &good), HYPO);
    assert_eq!(
        cleanup_product("Probable methyltransferase", &good),
        "putative methyltransferase"
    );
    assert_eq!(
        cleanup_product("ABC-type domain", &good),
        "ABC-type domain protein"
    );
    assert_eq!(cleanup_product("rep", &good), "rep"); // whitelisted
}
