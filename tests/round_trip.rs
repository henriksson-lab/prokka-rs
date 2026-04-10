//! Round-trip tests: write output → parse it back → verify consistency.

use prokka_rs::model::*;
use prokka_rs::config::ProkkaConfig;

fn sample_result() -> AnnotationResult {
    let mut contig = Contig::new("NC_001.1".into(), b"ATGGCTAAACCCGGGTTTAAATAA".to_vec());

    let mut cds1 = SeqFeature::new(
        FeatureType::CDS, "NC_001.1".into(), "Prodigal:2.6".into(),
        1, 21, Strand::Forward,
    );
    cds1.add_tag("ID", "T_00001");
    cds1.add_tag("locus_tag", "T_00001");
    cds1.add_tag("product", "test kinase");
    cds1.add_tag("gene", "tknA");
    cds1.add_tag("EC_number", "2.7.1.1");
    cds1.add_tag("db_xref", "COG:COG0001");
    cds1.add_tag("inference", "similar to AA sequence:UniProtKB:P12345");

    let mut trna = SeqFeature::new(
        FeatureType::TRNA, "NC_001.1".into(), "Aragorn:1.2".into(),
        22, 23, Strand::Reverse,
    );
    trna.add_tag("ID", "T_00002");
    trna.add_tag("locus_tag", "T_00002");
    trna.add_tag("product", "tRNA-Ala(TGC)");

    contig.features = vec![cds1, trna];

    let mut stats = AnnotationStats::default();
    stats.num_contigs = 1;
    stats.total_bp = 24;
    stats.feature_counts.insert("CDS".into(), 1);
    stats.feature_counts.insert("tRNA".into(), 1);

    AnnotationResult {
        contigs: vec![contig],
        stats,
        log: vec!["test log".into()],
    }
}

/// Round-trip GFF3: write → parse back key properties
#[test]
fn test_gff3_round_trip() {
    let result = sample_result();
    let mut buf = Vec::new();
    prokka_rs::output::gff::write_gff3(&mut buf, &result, 3).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Parse GFF3 back manually
    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines[0], "##gff-version 3");
    assert!(lines[1].starts_with("##sequence-region NC_001.1 1 24"),
        "Got: '{}'", lines[1]);

    // Find CDS line
    let cds_line = lines.iter().find(|l| l.contains("\tCDS\t")).unwrap();
    let fields: Vec<&str> = cds_line.split('\t').collect();
    assert_eq!(fields[0], "NC_001.1"); // seqid
    assert_eq!(fields[2], "CDS");       // type
    assert_eq!(fields[3], "1");         // start
    assert_eq!(fields[4], "21");        // end
    assert_eq!(fields[6], "+");         // strand
    assert_eq!(fields[7], "0");         // phase for CDS

    // Attributes should contain key tags
    let attrs = fields[8];
    assert!(attrs.contains("locus_tag=T_00001"));
    assert!(attrs.contains("product=test kinase"));
    assert!(attrs.contains("gene=tknA"));
    assert!(attrs.contains("EC_number=2.7.1.1"));

    // Find tRNA line
    let trna_line = lines.iter().find(|l| l.contains("\ttRNA\t")).unwrap();
    let fields: Vec<&str> = trna_line.split('\t').collect();
    assert_eq!(fields[6], "-"); // reverse strand
    assert_eq!(fields[7], "."); // phase for non-CDS

    // FASTA section
    assert!(output.contains("##FASTA\n"));
    assert!(output.contains(">NC_001.1\n"));
}

/// Round-trip TBL: write → parse back key properties
#[test]
fn test_tbl_round_trip() {
    let result = sample_result();
    let mut buf = Vec::new();
    prokka_rs::output::tbl::write_tbl(&mut buf, &result).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Parse TBL back
    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines[0], ">Feature NC_001.1");

    // CDS on forward strand: start < end
    assert_eq!(lines[1], "1\t21\tCDS");

    // Should have product and locus_tag tags (lowercase)
    assert!(output.contains("\t\t\tproduct\ttest kinase\n"));
    assert!(output.contains("\t\t\tlocus_tag\tT_00001\n"));
    assert!(output.contains("\t\t\tEC_number\t2.7.1.1\n"));

    // Should NOT have uppercase tags (ID, Parent, Name)
    assert!(!output.lines().any(|l| l.contains("\tID\t")));

    // tRNA on reverse strand: end < start
    let trna_line = lines.iter().find(|l| l.contains("tRNA")).unwrap();
    assert_eq!(*trna_line, "23\t22\ttRNA"); // end, start for reverse
}

/// Round-trip TSV: write → parse back all fields
#[test]
fn test_tsv_round_trip() {
    let result = sample_result();
    let mut buf = Vec::new();
    prokka_rs::output::tsv::write_tsv(&mut buf, &result).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines[0], "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct");

    // CDS line
    let cds_fields: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(cds_fields[0], "T_00001");   // locus_tag
    assert_eq!(cds_fields[1], "CDS");       // ftype
    assert_eq!(cds_fields[2], "21");        // length
    assert_eq!(cds_fields[3], "tknA");      // gene
    assert_eq!(cds_fields[4], "2.7.1.1");   // EC_number
    assert_eq!(cds_fields[5], "COG0001");   // COG
    assert_eq!(cds_fields[6], "test kinase"); // product

    // tRNA line
    let trna_fields: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(trna_fields[0], "T_00002");
    assert_eq!(trna_fields[1], "tRNA");
    assert_eq!(trna_fields[6], "tRNA-Ala(TGC)");
}

/// Round-trip FASTA (.fna): write → read back → verify identical sequence
#[test]
fn test_fna_round_trip() {
    let result = sample_result();
    let mut buf = Vec::new();
    prokka_rs::output::fasta::write_fna(&mut buf, &result.contigs).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Parse back
    let mut seq = String::new();
    let mut id = String::new();
    for line in output.lines() {
        if let Some(h) = line.strip_prefix('>') {
            id = h.split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.push_str(line);
        }
    }
    assert_eq!(id, "NC_001.1");
    assert_eq!(seq.as_bytes(), result.contigs[0].sequence.as_slice());
}

/// Round-trip .faa: translate and verify protein sequence
#[test]
fn test_faa_round_trip() {
    let result = sample_result();
    let mut buf = Vec::new();
    prokka_rs::output::fasta::write_faa(&mut buf, &result, 11).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Should have 1 record (CDS only, not tRNA)
    assert_eq!(output.matches('>').count(), 1);
    assert!(output.contains(">T_00001 test kinase\n"));

    // Verify protein: ATGGCTAAACCCGGGTTT = ATG GCT AAA CCC GGG TTT = M A K P G F
    // (positions 1-21 = ATGGCTAAACCCGGGTTTAAA, 7 codons)
    let protein_line: String = output.lines()
        .filter(|l| !l.starts_with('>'))
        .collect();
    assert!(protein_line.starts_with("MAKPGF"));
}

/// Round-trip TXT: write → verify format
#[test]
fn test_txt_round_trip() {
    let result = sample_result();
    let config = ProkkaConfig {
        genus: "Test".into(),
        species: "species".into(),
        strain: "str1".into(),
        ..Default::default()
    };
    let mut buf = Vec::new();
    prokka_rs::output::txt::write_txt(&mut buf, &result, &config).unwrap();
    let output = String::from_utf8(buf).unwrap();

    assert!(output.contains("organism: Test species str1"));
    assert!(output.contains("contigs: 1"));
    assert!(output.contains("bases: 24"));
    assert!(output.contains("CDS: 1"));
    assert!(output.contains("tRNA: 1"));
}

/// Test GFF3 escaping for special characters in product names
#[test]
fn test_gff3_special_chars_round_trip() {
    let mut contig = Contig::new("ctg1".into(), b"ACGT".to_vec());
    let mut cds = SeqFeature::new(
        FeatureType::CDS, "ctg1".into(), "test".into(), 1, 4, Strand::Forward,
    );
    cds.add_tag("product", "RNA-dependent RNA polymerase; antiviral");
    cds.add_tag("Note", "score=100% identity");
    contig.features = vec![cds];

    let result = AnnotationResult {
        contigs: vec![contig],
        stats: AnnotationStats::default(),
        log: Vec::new(),
    };

    let mut buf = Vec::new();
    prokka_rs::output::gff::write_gff3(&mut buf, &result, 3).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Semicolons and percent should be escaped
    let cds_line = output.lines().find(|l| l.contains("\tCDS\t")).unwrap();
    assert!(cds_line.contains("RNA-dependent RNA polymerase%3B antiviral"));
    assert!(cds_line.contains("score%3D100%25 identity"));
}

/// End-to-end pipeline round trip: run full pipeline, verify all outputs exist and parse
#[test]
fn test_pipeline_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let outdir = dir.path().join("output");
    let input = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("prokka")
        .join("test")
        .join("plasmid.fna");

    let mut config = ProkkaConfig {
        outdir: Some(outdir.clone()),
        prefix: Some("rt".into()),
        force: true,
        notrna: true,
        norrna: true,
        noanno: true,
        quiet: true,
        ..Default::default()
    };

    prokka_rs::pipeline::run(&input, &mut config).unwrap();

    // Read back each file and verify basic structure

    // .gff starts with ##gff-version
    let gff = std::fs::read_to_string(outdir.join("rt.gff")).unwrap();
    assert!(gff.starts_with("##gff-version 3"));
    assert!(gff.contains("##FASTA"));
    let cds_count = gff.lines().filter(|l| l.contains("\tCDS\t")).count();
    assert!(cds_count > 50, "Expected >50 CDS in GFF, got {}", cds_count);

    // .tsv has correct number of data lines (header + features)
    let tsv = std::fs::read_to_string(outdir.join("rt.tsv")).unwrap();
    let tsv_data_lines = tsv.lines().count() - 1; // minus header
    assert_eq!(tsv_data_lines, cds_count, "TSV data lines should match GFF CDS count");

    // .fna sequence matches input
    let input_seq: String = std::fs::read_to_string(&input).unwrap()
        .lines().filter(|l| !l.starts_with('>')).collect();
    let output_seq: String = std::fs::read_to_string(outdir.join("rt.fna")).unwrap()
        .lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(input_seq, output_seq);

    // .faa has same number of records as CDS
    let faa = std::fs::read_to_string(outdir.join("rt.faa")).unwrap();
    let faa_count = faa.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(faa_count, cds_count, ".faa records should match CDS count");

    // .ffn has same number of records as CDS (when only CDS features exist)
    let ffn = std::fs::read_to_string(outdir.join("rt.ffn")).unwrap();
    let ffn_count = ffn.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(ffn_count, cds_count, ".ffn records should match CDS count");

    // .txt has expected stats
    let txt = std::fs::read_to_string(outdir.join("rt.txt")).unwrap();
    assert!(txt.contains("contigs: 1"));
    assert!(txt.contains("bases: 56520"));

    // .tbl starts with >Feature
    let tbl = std::fs::read_to_string(outdir.join("rt.tbl")).unwrap();
    assert!(tbl.starts_with(">Feature"));
}
