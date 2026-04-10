/// Verify BLOSUM62 matrix self-scores are positive for standard amino acids.
#[test]
fn test_blosum62_self_scores() {
    let m = &blast_rs::matrix::BLOSUM62;
    let m_idx = blast_rs::input::aminoacid_to_ncbistdaa(b'M') as usize; // 12
    let a_idx = blast_rs::input::aminoacid_to_ncbistdaa(b'A') as usize; // 1
    assert_eq!(m[a_idx][a_idx], 4, "A-A should be 4");
    assert_eq!(m[m_idx][m_idx], 5, "M-M should be 5");
    assert_eq!(m[20][20], 11, "W-W should be 11");
}

/// Diagnostic test: verify blast-rs can find a self-hit (protein against itself).
/// If this fails, something is fundamentally wrong with the BLAST search.
#[test]
fn test_blast_self_hit() {
    // A known protein sequence
    let protein = b"MGIFDGKVAIITGGGKAKSIGYGIAVAYAKEGANLVLTGRNEQKLLDAKEELERLYGIKV";

    let query: Vec<u8> = protein.iter().map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();
    let subject = query.clone();

    // First try ungapped scan to see if seeds are found
    let ungapped = blast_rs::protein_lookup::protein_scan(
        &query, &subject,
        &blast_rs::matrix::BLOSUM62,
        3,    // word_size
        11.0, // threshold
        30,   // ungap_x_dropoff
    );
    eprintln!("Ungapped hits: {}", ungapped.len());
    for h in &ungapped {
        eprintln!("  score={} q={}..{} s={}..{}", h.score, h.query_start, h.query_end, h.subject_start, h.subject_end);
    }

    let hits = blast_rs::protein_lookup::protein_gapped_scan(
        &query, &subject,
        &blast_rs::matrix::BLOSUM62,
        3,    // word_size
        11.0, // threshold
        30,   // ungap_x_dropoff
        11,   // gap_open
        1,    // gap_extend
        50,   // gap_x_dropoff
        0,    // ungap_cutoff
    );
    eprintln!("Gapped hits: {}", hits.len());
    for h in &hits {
        eprintln!("  score={} q={}..{} s={}..{}", h.score, h.query_start, h.query_end, h.subject_start, h.subject_end);
    }

    assert!(!hits.is_empty(), "BLAST should find a self-hit (ungapped={} gapped={})", ungapped.len(), hits.len());
}

/// Test that blast-rs can find a hit between similar but not identical proteins.
#[test]
fn test_blast_similar_proteins() {
    // Two proteins with ~80% identity
    let query_str = b"MGIFDGKVAIITGGGKAKSIGYGIAVAYAKEGANLVLTGRNEQKLLDAKEELERLYGIKV";
    let subject_str = b"MGIFDGKVAIITGGGKAKSIGYGIAVAYAKEGANXVXTGRNEQKLXDAKEELERLYGIKV"; // 3 substitutions

    let query: Vec<u8> = query_str.iter().map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();
    let subject: Vec<u8> = subject_str.iter().map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();

    let hits = blast_rs::protein_lookup::protein_gapped_scan(
        &query, &subject,
        &blast_rs::matrix::BLOSUM62,
        3, 11.0, 30, 11, 1, 50, 0,
    );

    assert!(!hits.is_empty(), "BLAST should find a hit between similar proteins");
}

/// Test that the full annotation pipeline finds at least some hits on the sprot database.
/// This is a slow test since it searches the entire database.
#[test]
#[ignore] // Run with: cargo test -- --ignored test_blast_sprot
fn test_blast_sprot_finds_hits() {
    use std::path::PathBuf;

    let db_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("prokka/db/kingdom/Bacteria/sprot");
    if !db_path.exists() {
        return;
    }

    // Load database
    let db_file = std::fs::File::open(&db_path).unwrap();
    let subjects = blast_rs::input::parse_fasta(db_file);
    assert!(subjects.len() > 100, "Should have many subjects");

    // A known bacterial protein (replication initiation protein)
    let query_str = b"MKQIKEYLEEFVHSRLNKNIILRAAGFEYAKENPNFSQYYGNTVVSLPHRGKYGGPVNRIAPEMFHQIVAKPGERTFEGMFAIFKHRFPDWRDAES";
    let query: Vec<u8> = query_str.iter().map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();

    // Search against first 1000 subjects
    let mut best_score = 0i32;
    for subj in subjects.iter().take(1000) {
        let subj_aa: Vec<u8> = subj.sequence.iter().map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();
        let hits = blast_rs::protein_lookup::protein_gapped_scan(
            &query, &subj_aa,
            &blast_rs::matrix::BLOSUM62,
            3, 11.0, 30, 11, 1, 50, 0,
        );
        for h in &hits {
            if h.score > best_score {
                best_score = h.score;
                eprintln!("New best: {} score={} desc={}", subj.id, h.score, &subj.defline[..60.min(subj.defline.len())]);
            }
        }
    }

    assert!(best_score > 0, "Should find at least one hit in sprot");
}
