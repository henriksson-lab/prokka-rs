use std::io::Write;

use crate::error::ProkkaError;
use crate::model::{AnnotationResult, FeatureType};

/// Percent-encode a GFF3 attribute value per the GFF3 specification.
///
/// Characters that must be encoded: tab, newline, carriage return,
/// percent, semicolon, equals, ampersand, comma.
fn gff3_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.bytes() {
        match c {
            b'\t' => out.push_str("%09"),
            b'\n' => out.push_str("%0A"),
            b'\r' => out.push_str("%0D"),
            b'%' => out.push_str("%25"),
            b';' => out.push_str("%3B"),
            b'=' => out.push_str("%3D"),
            b'&' => out.push_str("%26"),
            b',' => out.push_str("%2C"),
            _ => out.push(c as char),
        }
    }
    out
}

/// Write GFF3 output.
///
/// Replicates Perl Prokka lines 1305-1389.
pub fn write_gff3(
    writer: &mut impl Write,
    result: &AnnotationResult,
    gff_version: u8,
) -> Result<(), ProkkaError> {
    writeln!(writer, "##gff-version {}", gff_version)?;

    // Sequence-region lines
    for contig in &result.contigs {
        writeln!(writer, "##sequence-region {} 1 {}", contig.id, contig.len())?;
    }

    // Feature lines
    for contig in &result.contigs {
        // Sort features: by start, then end descending, then Parent presence
        let mut sorted_features: Vec<_> = contig.features.iter().collect();
        sorted_features.sort_by(|a, b| {
            a.start.cmp(&b.start)
                .then(b.end.cmp(&a.end))
                .then(a.has_tag("Parent").cmp(&b.has_tag("Parent")))
        });

        for feature in &sorted_features {
            // GFF columns: seqid source type start end score strand phase attributes
            let score = match feature.score {
                Some(s) => format!("{}", s),
                None => ".".to_string(),
            };
            let strand = match feature.strand {
                crate::model::Strand::Forward => "+",
                crate::model::Strand::Reverse => "-",
            };
            let phase = if feature.feature_type == FeatureType::CDS {
                "0"
            } else {
                "."
            };

            // Build attributes with proper GFF3 escaping.
            // GFF3 convention: ID first, then other tags in insertion order.
            let mut attrs = Vec::new();

            // ID first (GFF3 convention, matches BioPerl behavior)
            if let Some(values) = feature.tags.get("ID") {
                for value in values {
                    attrs.push(format!("ID={}", gff3_escape(value)));
                }
            }

            // Add Name tag from gene (Perl Prokka adds this during GFF output)
            if let Some(gene) = feature.get_tag("gene") {
                attrs.push(format!("Name={}", gff3_escape(gene)));
            }

            // Remaining tags in insertion order
            for (key, values) in &feature.tags {
                if key == "ID" {
                    continue; // already emitted
                }
                for value in values {
                    attrs.push(format!("{}={}", gff3_escape(key), gff3_escape(value)));
                }
            }
            let attr_str = attrs.join(";");

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                feature.seq_id,
                feature.source,
                feature.feature_type.as_str(),
                feature.start,
                feature.end,
                score,
                strand,
                phase,
                attr_str
            )?;
        }
    }

    // FASTA section
    if !result.contigs.is_empty() {
        writeln!(writer, "##FASTA")?;
        for contig in &result.contigs {
            crate::output::fasta::write_fasta(writer, &contig.id, None, &contig.sequence)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gff3_escape() {
        assert_eq!(gff3_escape("simple text"), "simple text");
        assert_eq!(gff3_escape("a;b=c"), "a%3Bb%3Dc");
        assert_eq!(gff3_escape("100%"), "100%25");
        assert_eq!(gff3_escape("a\tb"), "a%09b");
    }
}
