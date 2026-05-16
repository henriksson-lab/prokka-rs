//! Ribosomal RNA prediction via Barrnap (default) or RNAmmer (`--rnammer`).
//!
//! Implements pipeline step 3. Skipped when `--kingdom Viruses` or
//! `--norrna` is set. Both backends produce `rRNA` features carrying a
//! `product` tag derived from the model name.

use std::path::Path;
use std::process::Command;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{FeatureType, SeqFeature, Strand};

/// Predict rRNA features using Barrnap (external process).
///
/// Replicates Perl Prokka lines 558-613:
/// - Command: `LC_ALL=C barrnap --kingdom {mode} --threads {cpus} --quiet {fna}`
/// - Parse GFF3 output
pub fn predict_rrna(
    fna_path: &Path,
    config: &ProkkaConfig,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    let barrnap_mode = match config.kingdom.barrnap_mode() {
        Some(m) => m,
        None => return Ok(Vec::new()), // Viruses: no rRNA search
    };

    // If --rnammer is set, try RNAmmer first
    if config.rnammer {
        match predict_rrna_rnammer(fna_path, barrnap_mode) {
            Ok(features) => return Ok(features),
            Err(ProkkaError::ToolNotFound { .. }) => {
                // Fall through to Barrnap
            }
            Err(e) => return Err(e),
        }
    }

    let output = Command::new("barrnap")
        .env("LC_ALL", "C")
        .arg("--kingdom")
        .arg(barrnap_mode)
        .arg("--threads")
        .arg(config.cpus.to_string())
        .arg("--quiet")
        .arg(fna_path)
        .output()
        .map_err(|_| ProkkaError::ToolNotFound {
            tool: "barrnap".to_string(),
        })?;

    if !output.status.success() {
        return Err(ProkkaError::ToolFailed {
            tool: "barrnap".to_string(),
            message: String::from_utf8_lossy(&output.stderr).to_string(),
        });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_gff3_rrna(&stdout)
}

/// Parse Barrnap's GFF3 output into a list of `rRNA` features.
///
/// The Perl pipeline removes the GFF `Name` attribute and keeps only
/// `product`; we read whichever of `product=` or `Name=` is present so
/// downstream code always has a product label.
pub fn parse_gff3_rrna(gff_output: &str) -> Result<Vec<SeqFeature>, ProkkaError> {
    let mut features = Vec::new();

    for line in gff_output.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let seq_id = fields[0].to_string();
        let source = fields[1].to_string();
        let start: usize = match fields[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: usize = match fields[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let strand = match fields[6] {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Forward,
        };

        // Parse attributes to get product (from Name= attribute)
        let attrs = fields[8];
        let mut product = None;
        for attr in attrs.split(';') {
            if let Some(val) = attr.strip_prefix("Name=") {
                // Prokka removes the Name tag and only keeps product
                product = Some(val.to_string());
            } else if let Some(val) = attr.strip_prefix("product=") {
                product = Some(val.to_string());
            }
        }

        let product = match product {
            Some(p) => p,
            None => continue,
        };

        let mut feature = SeqFeature::new(
            FeatureType::RRNA,
            seq_id,
            source.clone(),
            start,
            end,
            strand,
        );

        feature.add_tag("product", &product);
        feature.add_tag("inference", &format!("COORDINATES:profile:{}", source));

        features.push(feature);
    }

    Ok(features)
}

/// Predict rRNA using RNAmmer (alternative to Barrnap).
///
/// Runs `rnammer -S {mode} -xml {xmlfile} {fna}` and parses the XML output.
/// Matches Perl Prokka lines 573-604.
fn predict_rrna_rnammer(
    fna_path: &Path,
    mode: &str,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    let rnammer_mode = match mode {
        "bac" => "bac",
        "arc" => "arc",
        "mito" => "euk",
        _ => return Ok(Vec::new()),
    };

    // Write XML to temp file
    let xml_path = fna_path.with_extension("rnammer.xml");

    let status = Command::new("rnammer")
        .arg("-S")
        .arg(rnammer_mode)
        .arg("-xml")
        .arg(&xml_path)
        .arg(fna_path)
        .status()
        .map_err(|_| ProkkaError::ToolNotFound {
            tool: "rnammer".to_string(),
        })?;

    if !status.success() {
        let _ = std::fs::remove_file(&xml_path);
        return Err(ProkkaError::ToolFailed {
            tool: "rnammer".to_string(),
            message: "rnammer returned non-zero exit code".to_string(),
        });
    }

    let xml_content = std::fs::read_to_string(&xml_path).unwrap_or_default();
    let _ = std::fs::remove_file(&xml_path);

    parse_rnammer_xml(&xml_content)
}

/// Parse RNAmmer XML output into a list of `rRNA` features.
///
/// Walks each `<entry>` block extracting `sequenceEntry`, `mol`, `start`,
/// `stop`, `direction`. The `mol` description is rewritten via
/// `s_r` -> `S ribosomal `, matching Perl Prokka line 585 so labels read
/// e.g. `23S ribosomal RNA` rather than `23S_rRNA`.
fn parse_rnammer_xml(xml: &str) -> Result<Vec<SeqFeature>, ProkkaError> {
    let mut features = Vec::new();

    // Simple tag-based XML parsing (avoids XML crate dependency).
    // Each <entry> block contains the fields we need.
    for entry_block in xml.split("<entry>").skip(1) {
        let seq_id = extract_xml_tag(entry_block, "sequenceEntry").unwrap_or_default();
        let mut desc = extract_xml_tag(entry_block, "mol").unwrap_or_default();
        let start_str = extract_xml_tag(entry_block, "start").unwrap_or_default();
        let stop_str = extract_xml_tag(entry_block, "stop").unwrap_or_default();
        let direction_str = extract_xml_tag(entry_block, "direction").unwrap_or_default();

        if seq_id.is_empty() {
            continue;
        }

        let start: usize = match start_str.parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: usize = match stop_str.parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let strand = match direction_str.as_str() {
            "-1" | "-" => Strand::Reverse,
            _ => Strand::Forward,
        };

        // Perl: $desc =~ s/s_r/S ribosomal /i;
        desc = desc.replace("s_r", "S ribosomal ").replace("S_r", "S ribosomal ");

        let tool = "RNAmmer".to_string();
        let mut feature = SeqFeature::new(
            FeatureType::RRNA,
            seq_id,
            tool.clone(),
            start,
            end,
            strand,
        );
        feature.add_tag("product", &desc);
        feature.add_tag("inference", &format!("COORDINATES:profile:{}", tool));

        features.push(feature);
    }

    Ok(features)
}

/// Extract the trimmed text content of a simple XML element
/// (`<tag>...</tag>`) from `block`, or `None` if either delimiter is
/// missing. Used to avoid pulling in a full XML parser for the small,
/// well-known RNAmmer schema.
fn extract_xml_tag(block: &str, tag: &str) -> Option<String> {
    let open = format!("<{}>", tag);
    let close = format!("</{}>", tag);
    let start = block.find(&open)? + open.len();
    let end = block[start..].find(&close)? + start;
    Some(block[start..end].trim().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_barrnap_gff() {
        let gff = "\
##gff-version 3
contig1\tbarrnap:0.9\trRNA\t100\t1600\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA
contig1\tbarrnap:0.9\trRNA\t2000\t4900\t0\t-\t.\tName=23S_rRNA;product=23S ribosomal RNA
";
        let features = parse_gff3_rrna(gff).unwrap();
        assert_eq!(features.len(), 2);

        assert_eq!(features[0].feature_type, FeatureType::RRNA);
        assert_eq!(features[0].start, 100);
        assert_eq!(features[0].end, 1600);
        assert_eq!(features[0].strand, Strand::Forward);
        // When both Name= and product= exist, product= takes precedence (parsed last)
        assert_eq!(features[0].get_tag("product"), Some("16S ribosomal RNA"));

        assert_eq!(features[1].strand, Strand::Reverse);
    }
}
