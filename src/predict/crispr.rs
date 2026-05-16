//! CRISPR repeat detection via the external `minced` tool.
//!
//! Implements pipeline step 5. CRISPRs are emitted as `repeat_region`
//! features and (in the parent pipeline) included alongside RNAs in the
//! CDS overlap-filter input, because Prodigal will occasionally call ORFs
//! inside CRISPR arrays.

use std::path::Path;
use std::process::Command;

use crate::error::ProkkaError;
use crate::model::{FeatureType, SeqFeature, Strand};

/// Detect CRISPR repeats using minced (external process).
///
/// Replicates Perl Prokka lines 677-704:
/// - Command: `minced -gff {fna}`
/// - Parse GFF3 output for repeat_region features
pub fn predict_crispr(
    fna_path: &Path,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    let output = Command::new("minced")
        .arg("-gff")
        .arg(fna_path)
        .output()
        .map_err(|_| ProkkaError::ToolNotFound {
            tool: "minced".to_string(),
        })?;

    // minced may not be installed — treat as optional
    if !output.status.success() {
        return Ok(Vec::new());
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_minced_gff(&stdout)
}

/// Parse minced GFF3 output into a list of `repeat_region` features.
///
/// Reads only rows whose feature type is `repeat_region`. The GFF score
/// column carries the repeat-unit count, which becomes the `note` tag
/// (`"CRISPR with N repeat units"`), and the `rpt_unit_seq` attribute is
/// preserved. The Perl pipeline also strips the `ID` and `score` tags
/// (issue #328); we never add them in the first place.
pub fn parse_minced_gff(gff_output: &str) -> Result<Vec<SeqFeature>, ProkkaError> {
    let mut features = Vec::new();

    for line in gff_output.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        // Only process repeat_region features
        if fields[2] != "repeat_region" {
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
        let num_rpt: usize = fields[5].parse().unwrap_or(0);

        // Parse attributes for rpt_unit_seq
        let attrs = fields[8];
        let mut rpt_unit_seq = None;
        for attr in attrs.split(';') {
            if let Some(val) = attr.strip_prefix("rpt_unit_seq=") {
                rpt_unit_seq = Some(val.to_string());
            }
        }

        let mut feature = SeqFeature::new(
            FeatureType::RepeatRegion,
            seq_id,
            source,
            start,
            end,
            Strand::Forward,
        );

        feature.add_tag("rpt_family", "CRISPR");
        feature.add_tag("rpt_type", "direct");
        feature.add_tag("note", &format!("CRISPR with {} repeat units", num_rpt));
        if let Some(seq) = rpt_unit_seq {
            feature.add_tag("rpt_unit_seq", &seq);
        }

        features.push(feature);
    }

    Ok(features)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_minced_gff() {
        let gff = "\
##gff-version 3
contig1\tminced:0.3.0\trepeat_region\t11491\t12592\t16\t.\t.\tID=CRISPR1;rpt_family=CRISPR;rpt_unit_seq=GTCGAAAGACATTGCCCTGTTCCAAAGGGATTGAGAC
";
        let features = parse_minced_gff(gff).unwrap();
        assert_eq!(features.len(), 1);
        assert_eq!(features[0].feature_type, FeatureType::RepeatRegion);
        assert_eq!(features[0].start, 11491);
        assert_eq!(features[0].end, 12592);
        assert_eq!(features[0].get_tag("rpt_family"), Some("CRISPR"));
        assert_eq!(features[0].get_tag("note"), Some("CRISPR with 16 repeat units"));
    }
}
