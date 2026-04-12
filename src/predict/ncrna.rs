use std::path::Path;
use std::process::Command;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{FeatureType, SeqFeature, Strand};

/// Predict ncRNA features using Infernal/cmscan (external process).
///
/// Replicates Perl Prokka lines 618-663:
/// - Command: `cmscan -Z $dbsize --cut_ga --rfam --nohmmonly --fmt 2 --cpu $cpus
///            --tblout /dev/stdout -o /dev/null --noali $cmdb $fna`
/// - Parse tabular output (fmt 2)
pub fn predict_ncrna(
    fna_path: &Path,
    config: &ProkkaConfig,
    total_bp: usize,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    if !config.rfam || total_bp == 0 {
        return Ok(Vec::new());
    }

    let cmdb = config.dbdir.join("cm").join(config.kingdom.as_str());
    let index_file = format!("{}.i1m", cmdb.display());
    if !Path::new(&index_file).exists() {
        return Ok(Vec::new());
    }

    let dbsize = (total_bp * 2) as f64 / 1_000_000.0;

    let output = Command::new("cmscan")
        .arg("-Z")
        .arg(format!("{:.6}", dbsize))
        .arg("--cut_ga")
        .arg("--rfam")
        .arg("--nohmmonly")
        .arg("--fmt")
        .arg("2")
        .arg("--cpu")
        .arg(config.cpus.max(1).to_string())
        .arg("--tblout")
        .arg("/dev/stdout")
        .arg("-o")
        .arg("/dev/null")
        .arg("--noali")
        .arg(&cmdb)
        .arg(fna_path)
        .output()
        .map_err(|_| ProkkaError::ToolNotFound {
            tool: "cmscan".to_string(),
        })?;

    if !output.status.success() {
        return Err(ProkkaError::ToolFailed {
            tool: "cmscan".to_string(),
            message: String::from_utf8_lossy(&output.stderr).to_string(),
        });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_cmscan_output(&stdout)
}

/// Parse cmscan tabular output (fmt 2) into SeqFeature list.
pub fn parse_cmscan_output(
    output: &str,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    let mut features = Vec::new();

    for line in output.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 27 {
            continue;
        }

        // Field 2 (0-indexed) must be RF accession
        if !fields[2].starts_with("RF") {
            continue;
        }

        let seq_id = fields[3].to_string();

        // Skip overlapping hits (field 19 == "=")
        if fields.get(19) == Some(&"=") {
            continue;
        }

        let start1: usize = match fields[9].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let start2: usize = match fields[10].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let start = start1.min(start2);
        let end = start1.max(start2);

        let strand = if fields[11] == "-" {
            Strand::Reverse
        } else {
            Strand::Forward
        };

        let score: f64 = fields[16].parse().unwrap_or(0.0);

        let product = fields[1].to_string();
        let accession = fields[2].to_string();
        let note: String = fields[26..].join(" ");

        let tool = "Infernal".to_string();
        let mut feature = SeqFeature::new(
            FeatureType::MiscRNA,
            seq_id,
            tool.clone(),
            start,
            end,
            strand,
        );
        feature.score = Some(score);
        feature.add_tag("product", &product);
        feature.add_tag("inference", &format!("COORDINATES:profile:{}", tool));
        feature.add_tag("accession", &accession);
        if !note.is_empty() {
            feature.add_tag("Note", &note);
        }

        features.push(feature);
    }

    Ok(features)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{FeatureType, Strand};

    #[test]
    fn test_parse_cmscan_basic() {
        // Real cmscan fmt2 output has 27+ whitespace-delimited fields:
        // 0:target 1:product 2:accession 3:query_id 4:mdl 5:mdl_from 6:mdl_to
        // 7:? 8:? 9:seq_from 10:seq_to 11:strand 12:trunc 13:pass 14:gc
        // 15:bias 16:score 17:evalue 18:inc 19:overlap? 20-25:more 26+:description
        let line = "tRNA tRNA RF00005 contig1 cm 1 71 - - 100 170 + no 1 0.40 0.0 65.3 1.2e-14 ! - - - - - - tRNA gene";
        let features = parse_cmscan_output(line).unwrap();
        assert_eq!(features.len(), 1);
        assert_eq!(features[0].feature_type, FeatureType::MiscRNA);
        assert_eq!(features[0].seq_id, "contig1");
        assert_eq!(features[0].start, 100);
        assert_eq!(features[0].end, 170);
        assert_eq!(features[0].strand, Strand::Forward);
        assert_eq!(features[0].get_tag("product"), Some("tRNA"));
        assert_eq!(features[0].get_tag("accession"), Some("RF00005"));
    }

    #[test]
    fn test_parse_cmscan_skips_comments() {
        let output = "# this is a comment\n# another comment\n";
        let features = parse_cmscan_output(output).unwrap();
        assert_eq!(features.len(), 0);
    }

    #[test]
    fn test_parse_cmscan_skips_overlap() {
        // Field 19 is "=" for overlapping hits — should be skipped
        // Need 27+ fields; field[2] starts with RF, field[19] = "="
        let mut fields = vec!["a"; 30];
        fields[2] = "RF00001";
        fields[3] = "contig1";
        fields[9] = "100";
        fields[10] = "200";
        fields[11] = "+";
        fields[16] = "50.0";
        fields[19] = "="; // overlap marker
        let line = fields.join(" ");
        let features = parse_cmscan_output(&line).unwrap();
        assert_eq!(features.len(), 0); // should be skipped
    }
}
