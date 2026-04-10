use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::Contig;

const MAX_CONTIG_ID_LEN: usize = 37;

/// Load contigs from a FASTA file, sanitize sequences, and filter by length.
///
/// Replicates the input processing from Perl Prokka lines 417-469:
/// - Filter contigs shorter than mincontiglen
/// - Replace `|` with `_` in IDs
/// - Detect duplicate IDs
/// - Optionally rename contigs with --centre prefix
/// - Uppercase sequences, remove gaps (`*`, `-`), replace non-ACGT with N
pub fn load_and_sanitize_fasta(
    path: &Path,
    config: &ProkkaConfig,
) -> Result<Vec<Contig>, ProkkaError> {
    let file = File::open(path).map_err(|_| ProkkaError::FileNotReadable(path.to_path_buf()))?;
    let reader = BufReader::new(file);

    let mut contigs = Vec::new();
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = Vec::new();
    let mut ncontig: usize = 0;

    let contig_prefix = config.locustag.as_deref()
        .or(config.prefix.as_deref())
        .unwrap_or("");
    let contig_prefix = if contig_prefix.is_empty() {
        String::new()
    } else {
        format!("{}_", contig_prefix)
    };

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            // Process previous contig if any
            if let Some(id) = current_id.take() {
                let seq = sanitize_sequence(&current_seq);
                if seq.len() >= config.mincontiglen as usize {
                    contigs.push(Contig::new(id, seq));
                }
                current_seq.clear();
            }

            // Parse new header
            let raw_id = header.split_whitespace().next().unwrap_or("").to_string();
            let mut id = raw_id.replace('|', "_");

            // Check for duplicates
            if seen_ids.contains(&id) {
                return Err(ProkkaError::DuplicateContigId {
                    path: path.to_path_buf(),
                    id,
                });
            }
            seen_ids.insert(id.clone());

            ncontig += 1;

            // Rename if --centre is set
            if let Some(ref centre) = config.centre {
                id = format!("gnl|{}|{}{}", centre, contig_prefix, ncontig);
            }

            // Validate ID length
            if id.len() > MAX_CONTIG_ID_LEN {
                return Err(ProkkaError::ContigIdTooLong {
                    id: id.clone(),
                    len: id.len(),
                    max: MAX_CONTIG_ID_LEN,
                });
            }

            current_id = Some(id);
        } else {
            current_seq.extend_from_slice(line.as_bytes());
        }
    }

    // Process the last contig
    if let Some(id) = current_id {
        let seq = sanitize_sequence(&current_seq);
        if seq.len() >= config.mincontiglen as usize {
            contigs.push(Contig::new(id, seq));
        }
    }

    if contigs.is_empty() {
        return Err(ProkkaError::NoContigs {
            path: path.to_path_buf(),
        });
    }

    Ok(contigs)
}

/// Sanitize a DNA sequence:
/// - Uppercase
/// - Remove `*` and `-` (gaps/pads)
/// - Replace anything that's not ACGT with N
fn sanitize_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .filter_map(|&b| {
            let b = b.to_ascii_uppercase();
            match b {
                b'*' | b'-' => None,
                b'A' | b'C' | b'G' | b'T' => Some(b),
                _ => Some(b'N'),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sanitize_sequence() {
        assert_eq!(sanitize_sequence(b"ACGTacgt"), b"ACGTACGT");
        assert_eq!(sanitize_sequence(b"AC*GT"), b"ACGT");
        assert_eq!(sanitize_sequence(b"AC-GT"), b"ACGT");
        assert_eq!(sanitize_sequence(b"ACRTGT"), b"ACNTGT");
        assert_eq!(sanitize_sequence(b"ACWSGT"), b"ACNNGT");
    }
}
