//! FASTA-family output writers.
//!
//! Produces the four FASTA outputs Prokka emits:
//!
//! - `.fna` — raw input contigs ([`write_fna`])
//! - `.faa` — translated CDS protein sequences ([`write_faa`])
//! - `.ffn` — nucleotide sequence of each CDS/RNA feature ([`write_ffn`])
//! - `.fsa` — input contigs with `tbl2asn`-style modifier tags in the
//!   description line ([`write_fsa`])
//!
//! All records are wrapped at 60 columns to match BioPerl's `Bio::SeqIO`
//! FASTA writer, which the Perl Prokka pipeline uses.

use std::io::Write;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{AnnotationResult, Contig, FeatureType};

/// Write a single FASTA record with 60-character line wrapping.
///
/// Emits `>id` (plus `" desc"` if a non-empty description is supplied) on
/// the header line, then the sequence in 60-column blocks. This matches
/// the line wrapping BioPerl uses, which Perl Prokka relies on for all of
/// its FASTA output.
pub fn write_fasta(
    writer: &mut impl Write,
    id: &str,
    desc: Option<&str>,
    sequence: &[u8],
) -> Result<(), ProkkaError> {
    match desc {
        Some(d) if !d.is_empty() => writeln!(writer, ">{} {}", id, d)?,
        _ => writeln!(writer, ">{}", id)?,
    }
    for chunk in sequence.chunks(60) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Write all input contigs to a `.fna` FASTA file.
///
/// Each contig is emitted as `>id` followed by the contig sequence with no
/// description line. This is the unannotated "raw nucleotide" output.
pub fn write_fna(
    writer: &mut impl Write,
    contigs: &[Contig],
) -> Result<(), ProkkaError> {
    for contig in contigs {
        write_fasta(writer, &contig.id, None, &contig.sequence)?;
    }
    Ok(())
}

/// Write the protein FASTA (`.faa`) — translated CDS sequences.
///
/// Iterates every CDS feature, extracts its nucleotide subsequence
/// (reverse-complementing if on the minus strand), translates with the
/// configured genetic code, and writes the protein as `>locus_tag product`.
/// Replicates Perl Prokka line 1362 (`$p->translate(-codontable_id=>$gcode, -complete=>1)`).
pub fn write_faa(
    writer: &mut impl Write,
    result: &AnnotationResult,
    gcode: u8,
) -> Result<(), ProkkaError> {
    for contig in &result.contigs {
        for feature in &contig.features {
            if feature.feature_type != FeatureType::CDS {
                continue;
            }

            let dna = contig.extract(feature.start, feature.end, feature.strand);
            let protein = crate::annotate::blast::translate_dna(&dna, gcode);

            let locus_tag = feature.get_tag("locus_tag").unwrap_or("unknown");
            let product = feature.get_tag("product").unwrap_or("");

            write_fasta(writer, locus_tag, Some(product), &protein)?;
        }
    }
    Ok(())
}

/// Write the per-feature nucleotide FASTA (`.ffn`).
///
/// Emits the nucleotide subsequence of every CDS, rRNA, tRNA, tmRNA, and
/// misc_RNA feature, in their original strand orientation. Sister `gene`
/// and `mRNA` features are intentionally skipped. Replicates the feature-type
/// filter in Perl Prokka line 1364
/// (`if ($f->primary_tag =~ m/^(CDS|rRNA|tmRNA|tRNA|misc_RNA)$/)`).
pub fn write_ffn(
    writer: &mut impl Write,
    result: &AnnotationResult,
) -> Result<(), ProkkaError> {
    for contig in &result.contigs {
        for feature in &contig.features {
            let include = matches!(
                feature.feature_type,
                FeatureType::CDS
                    | FeatureType::RRNA
                    | FeatureType::TMRNA
                    | FeatureType::TRNA
                    | FeatureType::MiscRNA
            );
            if !include {
                continue;
            }

            let subseq = contig.extract(feature.start, feature.end, feature.strand);
            let locus_tag = feature.get_tag("locus_tag").unwrap_or("unknown");
            let product = feature.get_tag("product").unwrap_or("");

            write_fasta(writer, locus_tag, Some(product), &subseq)?;
        }
    }
    Ok(())
}

/// Write the annotated input FASTA (`.fsa`) for `tbl2asn`.
///
/// Each contig is written with a description containing the `tbl2asn`
/// modifier tags `[gcode=...] [organism=...] [strain=...]` (plus an optional
/// `[plasmid=...]`). Replicates Perl Prokka lines 1320-1326.
pub fn write_fsa(
    writer: &mut impl Write,
    result: &AnnotationResult,
    config: &ProkkaConfig,
) -> Result<(), ProkkaError> {
    let gcode = config.effective_gcode();
    let mut desc = format!(
        "[gcode={}] [organism={} {}] [strain={}]",
        gcode, config.genus, config.species, config.strain
    );
    if let Some(ref plasmid) = config.plasmid {
        desc.push_str(&format!(" [plasmid={}]", plasmid));
    }

    for contig in &result.contigs {
        write_fasta(writer, &contig.id, Some(&desc), &contig.sequence)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::*;

    #[test]
    fn test_fasta_60_char_wrapping() {
        let mut buf = Vec::new();
        let seq = b"A".repeat(150);
        write_fasta(&mut buf, "test", None, &seq).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[0], ">test");
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[2].len(), 60);
        assert_eq!(lines[3].len(), 30);
    }

    #[test]
    fn test_fasta_with_description() {
        let mut buf = Vec::new();
        write_fasta(&mut buf, "seq1", Some("my protein"), b"ACGT").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with(">seq1 my protein\n"));
    }

    #[test]
    fn test_fasta_no_description() {
        let mut buf = Vec::new();
        write_fasta(&mut buf, "seq1", None, b"ACGT").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with(">seq1\n"));
        assert!(!output.contains(">seq1 "));
    }

    #[test]
    fn test_faa_translates_cds() {
        let mut contig = Contig::new("ctg1".into(), b"NNATGGCTAAATAAANN".to_vec());
        let mut cds = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "t".into(), 3, 14, Strand::Forward,
        );
        cds.add_tag("locus_tag", "T_00001");
        cds.add_tag("product", "test protein");
        contig.features.push(cds);

        let result = AnnotationResult {
            contigs: vec![contig],
            stats: AnnotationStats::default(),
            log: Vec::new(),
        };

        let mut buf = Vec::new();
        write_faa(&mut buf, &result, 11).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains(">T_00001 test protein\n"));
        assert!(output.contains("MAK\n")); // ATG GCT AAA TAA -> MAK (stop stripped)
    }

    #[test]
    fn test_ffn_extracts_features() {
        let mut contig = Contig::new("ctg1".into(), b"AAACCCGGGTTT".to_vec());
        let mut cds = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "t".into(), 1, 6, Strand::Forward,
        );
        cds.add_tag("locus_tag", "T_00001");
        cds.add_tag("product", "p1");
        let mut trna = SeqFeature::new(
            FeatureType::TRNA, "ctg1".into(), "t".into(), 7, 12, Strand::Forward,
        );
        trna.add_tag("locus_tag", "T_00002");
        trna.add_tag("product", "p2");
        // Gene should NOT appear in ffn
        let gene = SeqFeature::new(
            FeatureType::Gene, "ctg1".into(), "t".into(), 1, 6, Strand::Forward,
        );
        contig.features = vec![cds, trna, gene];

        let result = AnnotationResult {
            contigs: vec![contig],
            stats: AnnotationStats::default(),
            log: Vec::new(),
        };

        let mut buf = Vec::new();
        write_ffn(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // Should have 2 records (CDS + tRNA), not 3 (no Gene)
        let count = output.matches('>').count();
        assert_eq!(count, 2);
        assert!(output.contains(">T_00001 p1\nAAACCC\n"));
        assert!(output.contains(">T_00002 p2\nGGGTTT\n"));
    }

    #[test]
    fn test_ffn_reverse_strand() {
        let mut contig = Contig::new("ctg1".into(), b"AACCGGTT".to_vec());
        let mut cds = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "t".into(), 3, 6, Strand::Reverse,
        );
        cds.add_tag("locus_tag", "T_00001");
        cds.add_tag("product", "p1");
        contig.features = vec![cds];

        let result = AnnotationResult {
            contigs: vec![contig],
            stats: AnnotationStats::default(),
            log: Vec::new(),
        };

        let mut buf = Vec::new();
        write_ffn(&mut buf, &result).unwrap();
        let output = String::from_utf8(buf).unwrap();
        // positions 3-6 = CCGG, reverse complement = CCGG
        assert!(output.contains("CCGG\n"));
    }

    #[test]
    fn test_fsa_has_metadata() {
        let contig = Contig::new("ctg1".into(), b"ACGT".to_vec());
        let result = AnnotationResult {
            contigs: vec![contig],
            stats: AnnotationStats::default(),
            log: Vec::new(),
        };
        let config = ProkkaConfig {
            genus: "Escherichia".into(),
            species: "coli".into(),
            strain: "K12".into(),
            plasmid: Some("pUC19".into()),
            ..Default::default()
        };

        let mut buf = Vec::new();
        write_fsa(&mut buf, &result, &config).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("[gcode=11]"));
        assert!(output.contains("[organism=Escherichia coli]"));
        assert!(output.contains("[strain=K12]"));
        assert!(output.contains("[plasmid=pUC19]"));
    }
}
