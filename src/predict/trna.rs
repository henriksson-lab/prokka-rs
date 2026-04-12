use std::path::Path;
use std::process::Command;

use regex::Regex;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{FeatureType, SeqFeature, Strand};

/// Predict tRNA and tmRNA features using Aragorn (external process).
///
/// Replicates Perl Prokka lines 476-553:
/// - Command: `aragorn -l -gc{gcode} {aragorn_opt} -w {fna}`
/// - Parse output: `>seqid` headers, then 5-field data lines
/// - Location regex: `(c)?\[-?(\d+),(\d+)\]`
/// - Filter: skip pseudo (contains `?`), skip start > end, skip > 500bp
pub fn predict_trna(
    fna_path: &Path,
    config: &ProkkaConfig,
    contig_lengths: &[(String, usize)],
) -> Result<Vec<SeqFeature>, ProkkaError> {
    let gcode = config.effective_gcode();
    let aragorn_opt = config.kingdom.aragorn_opt();

    let mut cmd = Command::new("aragorn");
    cmd.arg("-l")
        .arg(format!("-gc{}", gcode))
        .arg("-w");
    if !aragorn_opt.is_empty() {
        cmd.arg(aragorn_opt);
    }
    cmd.arg(fna_path);

    let output = cmd.output().map_err(|_| ProkkaError::ToolNotFound {
        tool: "aragorn".to_string(),
    })?;

    if !output.status.success() {
        return Err(ProkkaError::ToolFailed {
            tool: "aragorn".to_string(),
            message: String::from_utf8_lossy(&output.stderr).to_string(),
        });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_aragorn_output(&stdout, contig_lengths)
}

/// Parse Aragorn text output into SeqFeature list.
pub fn parse_aragorn_output(
    output: &str,
    contig_lengths: &[(String, usize)],
) -> Result<Vec<SeqFeature>, ProkkaError> {
    static COORD_RE: std::sync::LazyLock<Regex> = std::sync::LazyLock::new(|| {
        Regex::new(r"(c)?\[-?(\d+),(\d+)\]").unwrap()
    });
    let coord_re = &*COORD_RE;

    // Build a lookup for contig lengths
    let contig_len_map: std::collections::HashMap<&str, usize> = contig_lengths
        .iter()
        .map(|(id, len)| (id.as_str(), *len))
        .collect();

    let mut features = Vec::new();
    let mut current_sid: Option<String> = None;

    for line in output.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            current_sid = rest.split_whitespace().next().map(|s| s.to_string());
            continue;
        }

        let sid = match &current_sid {
            Some(s) => s.clone(),
            None => continue,
        };

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 5 {
            continue;
        }

        // First field must be a number (index)
        if fields[0].parse::<usize>().is_err() {
            continue;
        }

        // Skip pseudo/wacky genes (contain '?')
        if fields[1].contains('?') {
            continue;
        }

        // Parse location
        let location = fields[2];
        let caps = match coord_re.captures(location) {
            Some(c) => c,
            None => continue,
        };

        let is_revcom = caps.get(1).is_some();
        let start: usize = match caps[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: usize = match caps[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // Skip if start > end
        if start > end {
            continue;
        }

        // Correct coordinates: clamp to sequence bounds
        let seq_len = contig_len_map.get(sid.as_str()).copied().unwrap_or(usize::MAX);
        let start = start.max(1);
        let end = end.min(seq_len);

        // Skip if too big (>500bp)
        if end.saturating_sub(start) > 500 {
            continue;
        }

        let strand = if is_revcom { Strand::Reverse } else { Strand::Forward };

        // Determine feature type and product
        let (ftype, product, gene) = if fields[1].starts_with("tmRNA") {
            (FeatureType::TMRNA, "transfer-messenger RNA, SsrA".to_string(), Some("ssrA"))
        } else {
            let product = format!("{}{}", fields[1], fields[4]);
            (FeatureType::TRNA, product, None)
        };

        let tool = "Aragorn".to_string();
        let mut feature = SeqFeature::new(
            ftype,
            sid.clone(),
            tool.clone(),
            start,
            end,
            strand,
        );

        feature.add_tag("product", &product);
        feature.add_tag("inference", &format!("COORDINATES:profile:{}", tool));
        if let Some(g) = gene {
            feature.add_tag("gene", g);
        }

        features.push(feature);
    }

    Ok(features)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_aragorn_trna() {
        let output = "\
>contig1
1      tRNA-Ala              c[100,172]      34      (TGC)
2      tRNA-Arg              [500,572]       34      (TCT)
3      tmRNA                 [1000,1350]     35      (*)
";
        let contig_lengths = vec![("contig1".to_string(), 5000)];
        let features = parse_aragorn_output(output, &contig_lengths).unwrap();
        assert_eq!(features.len(), 3);

        assert_eq!(features[0].feature_type, FeatureType::TRNA);
        assert_eq!(features[0].strand, Strand::Reverse);
        assert_eq!(features[0].start, 100);
        assert_eq!(features[0].end, 172);
        assert_eq!(features[0].get_tag("product"), Some("tRNA-Ala(TGC)"));

        assert_eq!(features[1].strand, Strand::Forward);

        assert_eq!(features[2].feature_type, FeatureType::TMRNA);
        assert_eq!(features[2].get_tag("product"), Some("transfer-messenger RNA, SsrA"));
        assert_eq!(features[2].get_tag("gene"), Some("ssrA"));
    }

    #[test]
    fn test_skip_pseudo_and_large() {
        let output = "\
>contig1
1      tRNA-Ala?             [100,172]       34      (TGC)
2      tRNA-Arg              [100,700]       34      (TCT)
";
        let contig_lengths = vec![("contig1".to_string(), 5000)];
        let features = parse_aragorn_output(output, &contig_lengths).unwrap();
        // First skipped due to '?', second skipped due to >500bp
        assert_eq!(features.len(), 0);
    }
}
