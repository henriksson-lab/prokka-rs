/// NCBI genetic code translation tables.
///
/// Each table is represented by two strings:
/// - `amino_acids`: 64 amino acid translations for all codons in TCAG order
/// - `starts`: start codon markers ('M' = start, '-' = not start)
///
/// Codons are indexed as: T=0, C=1, A=2, G=3
/// Index = first*16 + second*4 + third
/// Order: TTT, TTC, TTA, TTG, TCT, TCC, TCA, TCG, TAT, TAC, TAA, TAG, TGT, TGC, TGA, TGG,
///        CTT, CTC, CTA, CTG, CCT, CCC, CCA, CCG, CAT, CAC, CAA, CAG, CGT, CGC, CGA, CGG,
///        ATT, ATC, ATA, ATG, ACT, ACC, ACA, ACG, AAT, AAC, AAA, AAG, AGT, AGC, AGA, AGG,
///        GTT, GTC, GTA, GTG, GCT, GCC, GCA, GCG, GAT, GAC, GAA, GAG, GGT, GGC, GGA, GGG
struct GeneticCode {
    amino_acids: &'static [u8; 64],
    starts: &'static [u8; 64],
}

// Table 1: Standard
const TABLE_1_AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_1_ST: &[u8; 64] = b"---M---------------M---------------M----------------------------";

// Table 2: Vertebrate Mitochondrial
const TABLE_2_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
const TABLE_2_ST: &[u8; 64] = b"--------------------------------MMMM---------------M------------";

// Table 3: Yeast Mitochondrial
const TABLE_3_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_3_ST: &[u8; 64] = b"----------------------------------MM----------------------------";

// Table 4: Mycoplasma/Spiroplasma (TGA=W)
const TABLE_4_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_4_ST: &[u8; 64] = b"--MM---------------M------------MMMM---------------M------------";

// Table 5: Invertebrate Mitochondrial
const TABLE_5_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
const TABLE_5_ST: &[u8; 64] = b"---M----------------------------MMMM---------------M------------";

// Table 6: Ciliate
const TABLE_6_AA: &[u8; 64] = b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_6_ST: &[u8; 64] = b"-----------------------------------M----------------------------";

// Table 9: Echinoderm/Flatworm Mitochondrial
const TABLE_9_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
const TABLE_9_ST: &[u8; 64] = b"-----------------------------------M---------------M------------";

// Table 10: Euplotid Nuclear
const TABLE_10_AA: &[u8; 64] = b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_10_ST: &[u8; 64] = b"-----------------------------------M----------------------------";

// Table 11: Bacterial/Archaeal/Plant Plastid
const TABLE_11_AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_11_ST: &[u8; 64] = b"---M---------------M------------MMMM---------------M------------";

// Table 12: Alternative Yeast Nuclear
const TABLE_12_AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_12_ST: &[u8; 64] = b"-------------------M---------------M----------------------------";

// Table 13: Ascidian Mitochondrial
const TABLE_13_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
const TABLE_13_ST: &[u8; 64] = b"---M------------------------------MM---------------M------------";

// Table 14: Alternative Flatworm Mitochondrial
const TABLE_14_AA: &[u8; 64] = b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
const TABLE_14_ST: &[u8; 64] = b"-----------------------------------M----------------------------";

// Table 15: Blepharisma Nuclear
const TABLE_15_AA: &[u8; 64] = b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_15_ST: &[u8; 64] = b"-----------------------------------M----------------------------";

// Table 16: Chlorophycean Mitochondrial
const TABLE_16_AA: &[u8; 64] = b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_16_ST: &[u8; 64] = b"-----------------------------------M----------------------------";

// Table 21: Trematode Mitochondrial
const TABLE_21_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
const TABLE_21_ST: &[u8; 64] = b"-----------------------------------M---------------M------------";

// Table 22: Scenedesmus obliquus Mitochondrial
const TABLE_22_AA: &[u8; 64] = b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_22_ST: &[u8; 64] = b"-----------------------------------M----------------------------";

// Table 23: Thraustochytrium Mitochondrial
const TABLE_23_AA: &[u8; 64] = b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_23_ST: &[u8; 64] = b"--------------------------------M--M---------------M------------";

// Table 24: Rhabdopleuridae Mitochondrial
const TABLE_24_AA: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG";
const TABLE_24_ST: &[u8; 64] = b"---M---------------M---------------M---------------M------------";

// Table 25: Candidate Division SR1 and Gracilibacteria
const TABLE_25_AA: &[u8; 64] = b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_25_ST: &[u8; 64] = b"---M-------------------------------M---------------M------------";

fn get_table(gcode: u8) -> &'static GeneticCode {
    // Use a static dispatch via match
    static TABLES: [GeneticCode; 26] = [
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 0 -> standard
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 1
        GeneticCode { amino_acids: TABLE_2_AA, starts: TABLE_2_ST }, // 2
        GeneticCode { amino_acids: TABLE_3_AA, starts: TABLE_3_ST }, // 3
        GeneticCode { amino_acids: TABLE_4_AA, starts: TABLE_4_ST }, // 4
        GeneticCode { amino_acids: TABLE_5_AA, starts: TABLE_5_ST }, // 5
        GeneticCode { amino_acids: TABLE_6_AA, starts: TABLE_6_ST }, // 6
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 7 (invalid, use standard)
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 8 (invalid, use standard)
        GeneticCode { amino_acids: TABLE_9_AA, starts: TABLE_9_ST }, // 9
        GeneticCode { amino_acids: TABLE_10_AA, starts: TABLE_10_ST }, // 10
        GeneticCode { amino_acids: TABLE_11_AA, starts: TABLE_11_ST }, // 11
        GeneticCode { amino_acids: TABLE_12_AA, starts: TABLE_12_ST }, // 12
        GeneticCode { amino_acids: TABLE_13_AA, starts: TABLE_13_ST }, // 13
        GeneticCode { amino_acids: TABLE_14_AA, starts: TABLE_14_ST }, // 14
        GeneticCode { amino_acids: TABLE_15_AA, starts: TABLE_15_ST }, // 15
        GeneticCode { amino_acids: TABLE_16_AA, starts: TABLE_16_ST }, // 16
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 17 (invalid)
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 18 (invalid)
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 19 (invalid)
        GeneticCode { amino_acids: TABLE_1_AA, starts: TABLE_1_ST }, // 20 (invalid)
        GeneticCode { amino_acids: TABLE_21_AA, starts: TABLE_21_ST }, // 21
        GeneticCode { amino_acids: TABLE_22_AA, starts: TABLE_22_ST }, // 22
        GeneticCode { amino_acids: TABLE_23_AA, starts: TABLE_23_ST }, // 23
        GeneticCode { amino_acids: TABLE_24_AA, starts: TABLE_24_ST }, // 24
        GeneticCode { amino_acids: TABLE_25_AA, starts: TABLE_25_ST }, // 25
    ];
    &TABLES[gcode.min(25) as usize]
}

/// Convert a nucleotide base to its index (T=0, C=1, A=2, G=3).
fn base_to_idx(b: u8) -> Option<usize> {
    match b {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Translate a single codon using the specified NCBI genetic code table.
pub fn translate_codon(codon: &[u8], gcode: u8) -> u8 {
    if codon.len() < 3 {
        return b'X';
    }
    match (base_to_idx(codon[0]), base_to_idx(codon[1]), base_to_idx(codon[2])) {
        (Some(a), Some(b), Some(c)) => {
            let idx = a * 16 + b * 4 + c;
            get_table(gcode).amino_acids[idx]
        }
        _ => b'X',
    }
}

/// Check if a codon is a valid start codon for the specified genetic code.
pub fn is_start_codon(codon: &[u8], gcode: u8) -> bool {
    if codon.len() < 3 {
        return false;
    }
    match (base_to_idx(codon[0]), base_to_idx(codon[1]), base_to_idx(codon[2])) {
        (Some(a), Some(b), Some(c)) => {
            let idx = a * 16 + b * 4 + c;
            get_table(gcode).starts[idx] == b'M'
        }
        _ => false,
    }
}

/// Check if a codon is a stop codon for the specified genetic code.
pub fn is_stop_codon(codon: &[u8], gcode: u8) -> bool {
    translate_codon(codon, gcode) == b'*'
}

/// Translate a DNA sequence to protein using the given genetic code.
///
/// Mimics BioPerl's translate(-codontable_id=>N, -complete=>1):
/// - First codon is translated to M if it's a valid start codon
/// - Stop codon (*) is stripped from the end
pub fn translate_dna(dna: &[u8], gcode: u8) -> Vec<u8> {
    if dna.len() < 3 {
        return Vec::new();
    }

    let mut protein = Vec::with_capacity(dna.len() / 3);
    let mut first = true;

    for i in (0..dna.len() - 2).step_by(3) {
        let codon = &dna[i..i + 3];
        let aa = translate_codon(codon, gcode);

        if first {
            first = false;
            if is_start_codon(codon, gcode) {
                protein.push(b'M');
            } else {
                protein.push(aa);
            }
        } else {
            protein.push(aa);
        }
    }

    // Strip trailing stop codon
    if protein.last() == Some(&b'*') {
        protein.pop();
    }

    protein
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_code() {
        assert_eq!(translate_codon(b"ATG", 1), b'M');
        assert_eq!(translate_codon(b"TAA", 1), b'*');
        assert_eq!(translate_codon(b"TAG", 1), b'*');
        assert_eq!(translate_codon(b"TGA", 1), b'*');
        assert_eq!(translate_codon(b"TTT", 1), b'F');
        assert_eq!(translate_codon(b"GCT", 1), b'A');
    }

    #[test]
    fn test_bacterial_code_11() {
        // Table 11 has same translations as table 1
        assert_eq!(translate_codon(b"ATG", 11), b'M');
        assert_eq!(translate_codon(b"TAA", 11), b'*');
        assert_eq!(translate_codon(b"TGA", 11), b'*');
        // But different start codons
        assert!(is_start_codon(b"ATG", 11));
        assert!(is_start_codon(b"GTG", 11));
        assert!(is_start_codon(b"TTG", 11));
        assert!(is_start_codon(b"CTG", 11));
        assert!(is_start_codon(b"ATT", 11));
    }

    #[test]
    fn test_mycoplasma_code_4() {
        // In table 4, TGA codes for Tryptophan (W), not stop
        assert_eq!(translate_codon(b"TGA", 4), b'W');
        assert!(!is_stop_codon(b"TGA", 4));
        // TAA and TAG are still stops
        assert_eq!(translate_codon(b"TAA", 4), b'*');
        assert_eq!(translate_codon(b"TAG", 4), b'*');
    }

    #[test]
    fn test_mitochondrial_code_5() {
        // Table 5 (Invertebrate Mitochondrial): AGA/AGG = Ser (not Arg)
        assert_eq!(translate_codon(b"AGA", 5), b'S');
        assert_eq!(translate_codon(b"AGG", 5), b'S');
        // Standard: AGA/AGG = Arg
        assert_eq!(translate_codon(b"AGA", 1), b'R');
        assert_eq!(translate_codon(b"AGG", 1), b'R');
    }

    #[test]
    fn test_translate_dna_table_4() {
        // TGA is W in table 4, not a stop codon
        let dna = b"ATGGCTTGA"; // ATG GCT TGA
        let protein_t1 = translate_dna(dna, 1);
        let protein_t4 = translate_dna(dna, 4);
        assert_eq!(protein_t1, b"MA");  // TGA is stop, stripped
        assert_eq!(protein_t4, b"MAW"); // TGA is tryptophan
    }

    #[test]
    fn test_start_codons() {
        // Standard (1): ATG, TTG, CTG are starts per NCBI
        assert!(is_start_codon(b"ATG", 1));
        assert!(is_start_codon(b"TTG", 1));
        assert!(is_start_codon(b"CTG", 1));
        assert!(!is_start_codon(b"GTG", 1));
        assert!(!is_start_codon(b"AAA", 1));

        // Bacterial (11): ATG, GTG, TTG, CTG, ATT are starts
        assert!(is_start_codon(b"ATG", 11));
        assert!(is_start_codon(b"GTG", 11));
        assert!(is_start_codon(b"TTG", 11));
        assert!(is_start_codon(b"CTG", 11));
        assert!(is_start_codon(b"ATT", 11));
    }

    #[test]
    fn test_unknown_bases() {
        assert_eq!(translate_codon(b"NNN", 11), b'X');
        assert_eq!(translate_codon(b"ANT", 11), b'X');
    }

    #[test]
    fn test_translate_dna_basic() {
        assert_eq!(translate_dna(b"ATGGCTAAA", 11), b"MAK");
        assert_eq!(translate_dna(b"ATGGCTAAATAA", 11), b"MAK"); // strip stop
        assert_eq!(translate_dna(b"GTGGCTAAA", 11), b"MAK"); // GTG start -> M
    }
}
