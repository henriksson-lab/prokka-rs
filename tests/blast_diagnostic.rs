fn simple_matrix() -> [[i32; blast_rs::matrix::AA_SIZE]; blast_rs::matrix::AA_SIZE] {
    let mut m = [[-2i32; blast_rs::matrix::AA_SIZE]; blast_rs::matrix::AA_SIZE];
    for i in 1..21 { m[i][i] = 5; }
    m
}

/// Verify simple protein matrix self-scores are positive for standard amino acids.
#[test]
fn test_simple_matrix_self_scores() {
    // Simple matrix: +5 diagonal for indices 1-20, -2 elsewhere
    let mut m = [[-2i32; blast_rs::matrix::AA_SIZE]; blast_rs::matrix::AA_SIZE];
    for i in 1..21 { m[i][i] = 5; }

    let m_idx = blast_rs::input::aminoacid_to_ncbistdaa(b'M') as usize;
    let a_idx = blast_rs::input::aminoacid_to_ncbistdaa(b'A') as usize;
    assert!(m[m_idx][m_idx] > 0, "M self-score: {}", m[m_idx][m_idx]);
    assert!(m[a_idx][a_idx] > 0, "A self-score: {}", m[a_idx][a_idx]);

    // 3-mer self-score should be 15 (5+5+5), well above threshold=11
    assert_eq!(m[m_idx][m_idx] + m[a_idx][a_idx] + m[a_idx][a_idx], 15);
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
        &simple_matrix(),
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
        &simple_matrix(),
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
        &simple_matrix(),
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
            &simple_matrix(),
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
