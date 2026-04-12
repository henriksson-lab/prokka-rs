use std::io::{BufRead, Write};
use std::path::Path;

use crate::error::ProkkaError;

/// Detect format of a protein file (FASTA, GenBank, or EMBL).
pub fn detect_format(path: &Path) -> Result<&'static str, ProkkaError> {
    let file = std::fs::File::open(path)
        .map_err(|_| ProkkaError::FileNotReadable(path.to_path_buf()))?;
    let reader = std::io::BufReader::new(file);

    for line in reader.lines().take(10) {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('>') {
            return Ok("fasta");
        }
        if trimmed.starts_with("LOCUS") {
            return Ok("genbank");
        }
        if trimmed.starts_with("ID ") {
            return Ok("embl");
        }
    }
    Ok("fasta") // default assumption
}

/// Convert a GenBank file to Prokka-format FASTA.
///
/// Extracts CDS translations with ~~~-delimited annotations:
/// `>accession EC_number~~~gene~~~product~~~`
///
/// Replicates the logic of `prokka-genbank_to_fasta_db`.
pub fn genbank_to_fasta(input: &Path, output: &Path) -> Result<usize, ProkkaError> {
    let content = std::fs::read_to_string(input)?;
    let mut out = std::fs::File::create(output)?;
    let mut count = 0;

    let mut in_features = false;
    let mut in_cds = false;
    let mut gene = String::new();
    let mut product = String::new();
    let mut ec_number = String::new();
    let mut locus_tag = String::new();
    let mut translation = String::new();
    let mut in_translation = false;
    let mut accession = String::new();

    for line in content.lines() {
        if line.starts_with("ACCESSION") {
            accession = line.split_whitespace().nth(1).unwrap_or("").to_string();
        }

        if line.starts_with("FEATURES") {
            in_features = true;
            continue;
        }

        if line.starts_with("ORIGIN") || line.starts_with("//") {
            // Flush current CDS if any
            if in_cds && !translation.is_empty() {
                write_protein_entry(
                    &mut out, &accession, &locus_tag, &ec_number, &gene, &product, &translation,
                )?;
                count += 1;
            }
            in_features = false;
            in_cds = false;
            in_translation = false;
            gene.clear();
            product.clear();
            ec_number.clear();
            locus_tag.clear();
            translation.clear();
            continue;
        }

        if !in_features {
            continue;
        }

        // Check for new feature (starts at column 6, 21-char key)
        if line.len() > 5 && line.as_bytes()[5] != b' ' {
            // Flush previous CDS
            if in_cds && !translation.is_empty() {
                write_protein_entry(
                    &mut out, &accession, &locus_tag, &ec_number, &gene, &product, &translation,
                )?;
                count += 1;
            }
            gene.clear();
            product.clear();
            ec_number.clear();
            locus_tag.clear();
            translation.clear();
            in_translation = false;

            let trimmed = line.trim();
            in_cds = trimmed.starts_with("CDS ");
            continue;
        }

        if !in_cds {
            continue;
        }

        let trimmed = line.trim();

        // Handle multi-line translation
        if in_translation {
            if trimmed.starts_with('/') || trimmed.starts_with("ORIGIN") {
                in_translation = false;
                // Remove trailing quote
                if translation.ends_with('"') {
                    translation.pop();
                }
            } else {
                let clean = trimmed.trim_end_matches('"');
                translation.push_str(clean);
                if trimmed.ends_with('"') {
                    in_translation = false;
                }
                continue;
            }
        }

        // Parse qualifiers
        if let Some(val) = extract_qualifier(trimmed, "/gene=") {
            gene = val;
        } else if let Some(val) = extract_qualifier(trimmed, "/product=") {
            product = val;
        } else if let Some(val) = extract_qualifier(trimmed, "/EC_number=") {
            ec_number = val;
        } else if let Some(val) = extract_qualifier(trimmed, "/locus_tag=") {
            locus_tag = val;
        } else if let Some(rest) = trimmed.strip_prefix("/translation=\"") {
            if let Some(stripped) = rest.strip_suffix('"') {
                translation = stripped.to_string();
            } else {
                translation = rest.to_string();
                in_translation = true;
            }
        }
    }

    Ok(count)
}

/// Convert an EMBL file to Prokka-format FASTA.
pub fn embl_to_fasta(input: &Path, output: &Path) -> Result<usize, ProkkaError> {
    let content = std::fs::read_to_string(input)?;
    let mut out = std::fs::File::create(output)?;
    let mut count = 0;

    let mut in_cds = false;
    let mut gene = String::new();
    let mut product = String::new();
    let mut ec_number = String::new();
    let mut locus_tag = String::new();
    let mut translation = String::new();
    let mut in_translation = false;
    let mut accession = String::new();

    for line in content.lines() {
        if line.starts_with("AC ") {
            accession = line[5..].trim().trim_end_matches(';').to_string();
        }

        if line.starts_with("FT ") {
            let rest = line[5..].trim();

            // New feature key (not indented qualifier)
            if !rest.starts_with('/') && !rest.is_empty() && !rest.starts_with("\"") {
                // Flush previous CDS
                if in_cds && !translation.is_empty() {
                    write_protein_entry(
                        &mut out, &accession, &locus_tag, &ec_number, &gene, &product, &translation,
                    )?;
                    count += 1;
                }
                gene.clear();
                product.clear();
                ec_number.clear();
                locus_tag.clear();
                translation.clear();
                in_translation = false;

                in_cds = rest.starts_with("CDS ");
                continue;
            }

            if !in_cds {
                continue;
            }

            if in_translation {
                let clean = rest.trim_end_matches('"');
                translation.push_str(clean);
                if rest.ends_with('"') {
                    in_translation = false;
                }
                continue;
            }

            if let Some(val) = extract_qualifier(rest, "/gene=") {
                gene = val;
            } else if let Some(val) = extract_qualifier(rest, "/product=") {
                product = val;
            } else if let Some(val) = extract_qualifier(rest, "/EC_number=") {
                ec_number = val;
            } else if let Some(val) = extract_qualifier(rest, "/locus_tag=") {
                locus_tag = val;
            } else if let Some(r) = rest.strip_prefix("/translation=\"") {
                if let Some(stripped) = r.strip_suffix('"') {
                    translation = stripped.to_string();
                } else {
                    translation = r.to_string();
                    in_translation = true;
                }
            }
        } else if line.starts_with("SQ") || line.starts_with("//") {
            if in_cds && !translation.is_empty() {
                write_protein_entry(
                    &mut out, &accession, &locus_tag, &ec_number, &gene, &product, &translation,
                )?;
                count += 1;
            }
            in_cds = false;
            in_translation = false;
            gene.clear();
            product.clear();
            ec_number.clear();
            locus_tag.clear();
            translation.clear();
        }
    }

    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_fasta() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("test.faa");
        std::fs::write(&p, ">seq1\nACGT\n").unwrap();
        assert_eq!(detect_format(&p).unwrap(), "fasta");
    }

    #[test]
    fn test_detect_genbank() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("test.gbk");
        std::fs::write(&p, "LOCUS       test\nFEATURES\n").unwrap();
        assert_eq!(detect_format(&p).unwrap(), "genbank");
    }

    #[test]
    fn test_detect_embl() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("test.embl");
        std::fs::write(&p, "ID   test; SV 1\nFT   CDS\n").unwrap();
        assert_eq!(detect_format(&p).unwrap(), "embl");
    }

    #[test]
    fn test_genbank_to_fasta() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("test.gbk");
        let output = dir.path().join("test.faa");
        std::fs::write(&input, "\
LOCUS       test
ACCESSION   TEST001
FEATURES             Location/Qualifiers
     CDS             1..300
                     /gene=\"abc\"
                     /product=\"test protein\"
                     /EC_number=\"1.2.3.4\"
                     /locus_tag=\"T_00001\"
                     /translation=\"MAKTPGF\"
ORIGIN
//
").unwrap();
        let count = genbank_to_fasta(&input, &output).unwrap();
        assert_eq!(count, 1);
        let fasta = std::fs::read_to_string(&output).unwrap();
        assert!(fasta.contains(">T_00001 1.2.3.4~~~abc~~~test protein~~~"));
        assert!(fasta.contains("MAKTPGF"));
    }
}

fn extract_qualifier(line: &str, prefix: &str) -> Option<String> {
    line.strip_prefix(prefix)
        .map(|val| val.trim_matches('"').to_string())
}

fn write_protein_entry(
    out: &mut impl Write,
    accession: &str,
    locus_tag: &str,
    ec_number: &str,
    gene: &str,
    product: &str,
    translation: &str,
) -> Result<(), ProkkaError> {
    let id = if !locus_tag.is_empty() {
        locus_tag
    } else {
        accession
    };
    // Prokka ~~~-delimited format: EC~~~gene~~~product~~~COG
    writeln!(out, ">{} {}~~~{}~~~{}~~~", id, ec_number, gene, product)?;
    for chunk in translation.as_bytes().chunks(60) {
        out.write_all(chunk)?;
        writeln!(out)?;
    }
    Ok(())
}
