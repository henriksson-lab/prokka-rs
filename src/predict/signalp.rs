use std::io::Write as IoWrite;
use std::path::Path;
use std::process::Command;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{Contig, FeatureType, SeqFeature, Strand};

const SIGNALP_MAXSEQ: usize = 10_000;

/// Detect signal peptides using SignalP (external process).
///
/// Replicates Perl Prokka lines 801-938.
/// Supports SignalP versions 3, 4, and 5.
pub fn predict_signalp(
    outdir: &Path,
    contigs: &[Contig],
    config: &ProkkaConfig,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    // SignalP requires --gram
    let gram = match &config.gram {
        Some(g) => g.clone(),
        None => return Ok(Vec::new()),
    };

    // Only for Bacteria
    if config.kingdom != crate::config::Kingdom::Bacteria {
        return Ok(Vec::new());
    }

    // Detect SignalP version
    let sigpver = detect_signalp_version()?;
    if sigpver == 0 {
        return Ok(Vec::new());
    }

    let gram_str = if gram.contains('+') || gram.to_lowercase().contains("pos") {
        "gram+"
    } else {
        "gram-"
    };

    let gcode = config.effective_gcode();

    // Collect CDS proteins, write to temp file
    let spout_path = outdir.join("signalp.faa");
    let mut cds_map: Vec<(usize, usize)> = Vec::new(); // (contig_idx, feature_idx)

    {
        let mut faa = std::fs::File::create(&spout_path)?;
        let mut count = 0usize;
        for (ci, contig) in contigs.iter().enumerate() {
            for (fi, feature) in contig.features.iter().enumerate() {
                if feature.feature_type != FeatureType::CDS {
                    continue;
                }
                count += 1;
                let dna = contig.extract(feature.start, feature.end, feature.strand);
                let protein = crate::codon_table::translate_dna(&dna, gcode);
                writeln!(faa, ">{}", count)?;
                faa.write_all(&protein)?;
                writeln!(faa)?;
                cds_map.push((ci, fi));
            }
        }
        if count > SIGNALP_MAXSEQ {
            let _ = std::fs::remove_file(&spout_path);
            eprintln!("Skipping SignalP: {} CDS exceeds limit of {}", count, SIGNALP_MAXSEQ);
            return Ok(Vec::new());
        }
    }

    // Run SignalP
    let output_text = run_signalp(sigpver, gram_str, &spout_path, outdir)?;

    // Parse output
    let mut features = Vec::new();
    let tool = format!("SignalP:{}", sigpver);

    for line in output_text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();

        let parsed = match sigpver {
            3 => parse_signalp_v3(&fields),
            4 => parse_signalp_v4(&fields),
            5 => parse_signalp_v5(&fields),
            _ => None,
        };

        if let Some((cds_idx, cleave_pos, note)) = parsed {
            let idx = cds_idx.saturating_sub(1); // 1-based to 0-based
            if idx >= cds_map.len() {
                continue;
            }
            let (ci, fi) = cds_map[idx];
            if ci >= contigs.len() || fi >= contigs[ci].features.len() {
                continue;
            }
            let parent = &contigs[ci].features[fi];

            // Convert cleavage position to DNA coordinates
            let (start, end) = if parent.strand == Strand::Forward {
                let s = parent.start;
                let e = s + cleave_pos * 3 - 1;
                (s, e)
            } else {
                let s = parent.end;
                let e = s - (cleave_pos * 3 - 1);
                (e.min(s), s.max(e))
            };

            let mut sigpep = SeqFeature::new(
                FeatureType::SigPeptide,
                parent.seq_id.clone(),
                tool.clone(),
                start,
                end,
                parent.strand,
            );
            sigpep.frame = 0; // compulsory for peptides
            sigpep.add_tag("product", "putative signal peptide");
            sigpep.add_tag("inference", &format!("ab initio prediction:{}", tool));
            sigpep.add_tag("note", &note);

            features.push(sigpep);
        }
    }

    let _ = std::fs::remove_file(&spout_path);

    Ok(features)
}

fn detect_signalp_version() -> Result<u8, ProkkaError> {
    // Try `signalp -version` first (v5), then `signalp -v` (v3/v4)
    if let Ok(out) = Command::new("signalp").arg("-version").output() {
        let s = String::from_utf8_lossy(&out.stdout);
        if let Some(ver) = s.split_whitespace().find_map(|w| w.parse::<f64>().ok()) {
            return Ok(ver as u8);
        }
    }
    if let Ok(out) = Command::new("signalp").arg("-v").stdin(std::process::Stdio::null()).output() {
        let s = String::from_utf8_lossy(&out.stderr);
        let combined = format!("{}{}", s, String::from_utf8_lossy(&out.stdout));
        if combined.contains("SignalP-3") || combined.contains("3.") {
            return Ok(3);
        }
        if combined.contains("SignalP-4") || combined.contains("4.") {
            return Ok(4);
        }
    }
    Ok(0) // Not found
}

fn run_signalp(ver: u8, gram: &str, faa_path: &Path, outdir: &Path) -> Result<String, ProkkaError> {
    match ver {
        3 => {
            let out = Command::new("signalp")
                .arg("-t").arg(gram)
                .arg("-f").arg("short")
                .arg("-m").arg("hmm")
                .arg(faa_path)
                .stderr(std::process::Stdio::null())
                .output()
                .map_err(|_| ProkkaError::ToolNotFound { tool: "signalp".into() })?;
            Ok(String::from_utf8_lossy(&out.stdout).to_string())
        }
        4 => {
            let out = Command::new("signalp")
                .arg("-t").arg(gram)
                .arg("-f").arg("short")
                .arg(faa_path)
                .stderr(std::process::Stdio::null())
                .output()
                .map_err(|_| ProkkaError::ToolNotFound { tool: "signalp".into() })?;
            Ok(String::from_utf8_lossy(&out.stdout).to_string())
        }
        5 => {
            let summary_path = outdir.join("signalp_summary.signalp5");
            let _ = Command::new("signalp")
                .arg("-tmp").arg(outdir)
                .arg("-prefix").arg(outdir.join("signalp"))
                .arg("-org").arg(gram)
                .arg("-format").arg("short")
                .arg("-fasta").arg(faa_path)
                .stderr(std::process::Stdio::null())
                .output();
            let text = std::fs::read_to_string(&summary_path).unwrap_or_default();
            let _ = std::fs::remove_file(&summary_path);
            Ok(text)
        }
        _ => Ok(String::new()),
    }
}

/// Parse SignalP v3 short output. Returns (cds_index, cleavage_pos, note).
fn parse_signalp_v3(fields: &[&str]) -> Option<(usize, usize, String)> {
    if fields.len() != 7 || fields[6] != "Y" {
        return None;
    }
    let idx: usize = fields[0].parse().ok()?;
    let cleave: usize = fields[3].parse().ok()?;
    let prob = fields[5];
    Some((idx, cleave, format!("predicted cleavage at residue {} with probability {}", fields[3], prob)))
}

/// Parse SignalP v4 short output.
fn parse_signalp_v4(fields: &[&str]) -> Option<(usize, usize, String)> {
    if fields.len() != 12 || fields[9] != "Y" {
        return None;
    }
    let idx: usize = fields[0].parse().ok()?;
    let cleave: usize = fields[2].parse().ok()?;
    Some((idx, cleave, format!("predicted cleavage at residue {}", fields[2])))
}

/// Parse SignalP v5 short output.
fn parse_signalp_v5(fields: &[&str]) -> Option<(usize, usize, String)> {
    if fields.len() != 12 {
        return None;
    }
    // Check if it's SP, TAT, or LIPO type
    if !fields[1].starts_with("SP") && !fields[1].starts_with("TAT") && !fields[1].starts_with("LIPO") {
        return None;
    }
    let idx: usize = fields[0].parse().ok()?;
    let tp_type = fields[1];
    let tp_prob = match tp_type {
        t if t.starts_with("SP") => fields[2],
        t if t.starts_with("TAT") => fields[3],
        t if t.starts_with("LIPO") => fields[4],
        _ => return None,
    };

    // Parse cleavage position from field 8 (format: "pos1-pos2.")
    let cleave_str = fields[8];
    let cleave: usize = cleave_str.split('-').next()?.trim_end_matches('.').parse().ok()?;
    let clprob = fields[11];

    Some((idx, cleave, format!(
        "{} (Probability: {}), predicted cleavage at residue {} with probability {}",
        tp_type, tp_prob, cleave, clprob
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_signalp_v3() {
        let fields = vec!["1", "Y", "0.900", "25", "Y", "0.95", "Y"];
        let result = parse_signalp_v3(&fields);
        assert!(result.is_some());
        let (idx, cleave, note) = result.unwrap();
        assert_eq!(idx, 1);
        assert_eq!(cleave, 25);
        assert!(note.contains("probability 0.95"));
    }

    #[test]
    fn test_parse_signalp_v3_no_sigpep() {
        let fields = vec!["1", "Y", "0.100", "0", "N", "0.10", "N"];
        assert!(parse_signalp_v3(&fields).is_none());
    }

    #[test]
    fn test_parse_signalp_v4() {
        let fields = vec!["1", "SP", "21", "0.9", "0.1", "0.0", "0.0", "0.0", "0.0", "Y", "0.8", "?"];
        let result = parse_signalp_v4(&fields);
        assert!(result.is_some());
        let (idx, cleave, _) = result.unwrap();
        assert_eq!(idx, 1);
        assert_eq!(cleave, 21);
    }

    #[test]
    fn test_parse_signalp_v5() {
        let fields = vec!["1", "SP(Sec/SPI)", "0.95", "0.01", "0.02", "0.01", "0.01", "0.0", "21-22.", "ANA-AT.", "0.0", "0.95"];
        let result = parse_signalp_v5(&fields);
        assert!(result.is_some());
        let (idx, cleave, note) = result.unwrap();
        assert_eq!(idx, 1);
        assert_eq!(cleave, 21);
        assert!(note.contains("SP(Sec/SPI)"));
    }
}
