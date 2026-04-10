use std::io::Write;
use std::path::Path;
use std::process::Command;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;


/// Generate GenBank (.gbk) and Sequin (.sqn) files.
///
/// Currently uses tbl2asn as an external process (same approach as Perl Prokka).
/// Requires .fsa and .tbl files to already exist in the output directory.
///
/// Replicates Perl Prokka lines 1416-1437.
pub fn generate_genbank(
    outdir: &Path,
    prefix: &str,
    config: &ProkkaConfig,
) -> Result<(), ProkkaError> {
    let fsa_path = outdir.join(format!("{}.fsa", prefix));
    let err_path = outdir.join(format!("{}.err", prefix));

    if !fsa_path.exists() {
        return Err(ProkkaError::Other(format!(
            "Missing .fsa file: {}",
            fsa_path.display()
        )));
    }

    // tbl2asn options: -M b for >10000 contigs, -M n otherwise
    let num_contigs = count_fasta_records(&fsa_path);
    let m_opt = if num_contigs > 10_000 { "b" } else { "n" };

    let annotation_comment = format!(
        "Annotated using prokka-rs {} from https://github.com/henriksson-lab/prokka-rs",
        env!("CARGO_PKG_VERSION")
    );

    let status = Command::new("tbl2asn")
        .arg("-V").arg("b")
        .arg("-a").arg("r10k")
        .arg("-l").arg("paired-ends")
        .arg("-M").arg(m_opt)
        .arg("-N").arg(config.accver.to_string())
        .arg("-y").arg(&annotation_comment)
        .arg("-Z").arg(&err_path)
        .arg("-i").arg(&fsa_path)
        .stderr(std::process::Stdio::null())
        .status();

    match status {
        Ok(s) if s.success() => {}
        Ok(_) => {
            return Err(ProkkaError::ToolFailed {
                tool: "tbl2asn".to_string(),
                message: "tbl2asn returned non-zero exit code".to_string(),
            });
        }
        Err(_) => {
            // tbl2asn not available — try to write a basic GenBank file ourselves
            if !config.quiet {
                eprintln!("tbl2asn not found, writing basic GenBank file.");
            }
            let gbk_path = outdir.join(format!("{}.gbk", prefix));
            let mut f = std::fs::File::create(&gbk_path)?;
            write_basic_genbank(&mut f, outdir, prefix, config)?;
            return Ok(());
        }
    }

    // Clean up tbl2asn output files
    let _ = std::fs::remove_file(outdir.join("errorsummary.val"));
    for ext in &["dr", "fixedproducts", "ecn", "val"] {
        let _ = std::fs::remove_file(outdir.join(format!("{}.{}", prefix, ext)));
    }

    // Fix tbl2asn output: rename .gbf to .gbk and fix COORDINATES spacing
    let gbf_path = outdir.join(format!("{}.gbf", prefix));
    let gbk_path = outdir.join(format!("{}.gbk", prefix));
    if gbf_path.exists() {
        let content = std::fs::read_to_string(&gbf_path)?;
        let fixed = content.replace("COORDINATES: profile", "COORDINATES:profile");
        std::fs::write(&gbk_path, fixed)?;
        let _ = std::fs::remove_file(&gbf_path);
    }

    Ok(())
}

/// Write a basic GenBank flat file when tbl2asn is not available.
fn write_basic_genbank(
    writer: &mut impl Write,
    outdir: &Path,
    prefix: &str,
    config: &ProkkaConfig,
) -> Result<(), ProkkaError> {
    // Read the .fsa to get sequences
    let fsa_path = outdir.join(format!("{}.fsa", prefix));
    let fsa_content = std::fs::read_to_string(&fsa_path)?;

    // Parse simple FASTA
    let mut contigs: Vec<(String, String)> = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in fsa_content.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                contigs.push((current_id.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_id = header.split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_seq.push_str(line);
        }
    }
    if !current_id.is_empty() {
        contigs.push((current_id, current_seq));
    }

    let date = chrono::Local::now().format("%d-%b-%Y").to_string().to_uppercase();

    for (id, seq) in &contigs {
        let len = seq.len();

        // LOCUS line
        writeln!(
            writer,
            "LOCUS       {:<16} {:>10} bp    DNA     linear   UNK {}",
            id, len, date
        )?;
        writeln!(writer, "DEFINITION  {} {} {} {}.",
            config.genus, config.species, config.strain,
            config.plasmid.as_deref().unwrap_or(""))?;
        writeln!(writer, "ACCESSION   ")?;
        writeln!(writer, "VERSION")?;
        writeln!(writer, "KEYWORDS    .")?;
        writeln!(writer, "SOURCE      {} {}", config.genus, config.species)?;
        writeln!(writer, "  ORGANISM  {} {}", config.genus, config.species)?;
        writeln!(writer, "            Unclassified.")?;
        writeln!(writer, "COMMENT     Annotated using prokka-rs {} from",
            env!("CARGO_PKG_VERSION"))?;
        writeln!(writer, "            https://github.com/henriksson-lab/prokka-rs")?;
        writeln!(writer, "FEATURES             Location/Qualifiers")?;
        writeln!(writer, "     source          1..{}", len)?;
        writeln!(writer, "                     /organism=\"{} {}\"", config.genus, config.species)?;
        writeln!(writer, "                     /mol_type=\"genomic DNA\"")?;
        writeln!(writer, "                     /strain=\"{}\"", config.strain)?;

        // TODO: Read .tbl and add features here

        // ORIGIN section
        writeln!(writer, "ORIGIN")?;
        let seq_bytes = seq.as_bytes();
        let mut pos = 0;
        while pos < seq_bytes.len() {
            write!(writer, "{:>9}", pos + 1)?;
            for group in 0..6 {
                let start = pos + group * 10;
                let end = (start + 10).min(seq_bytes.len());
                if start >= seq_bytes.len() {
                    break;
                }
                write!(writer, " ")?;
                writer.write_all(&seq_bytes[start..end].to_ascii_lowercase())?;
            }
            writeln!(writer)?;
            pos += 60;
        }
        writeln!(writer, "//")?;
    }

    Ok(())
}

/// Count FASTA records in a file.
fn count_fasta_records(path: &Path) -> usize {
    std::fs::read_to_string(path)
        .map(|s| s.lines().filter(|l| l.starts_with('>')).count())
        .unwrap_or(0)
}
