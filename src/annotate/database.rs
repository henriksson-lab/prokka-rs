//! Annotation database descriptors and `~~~`-delimited defline parsing.
//!
//! Prokka stores its reference protein annotations as FASTA records whose
//! description line is split into four fields by `~~~`:
//! `EC_number~~~gene~~~product~~~COG`. Both the BLAST and HMMER annotators
//! consume these deflines through [`parse_annotation_header`].

/// Search backend used to query an annotation database.
///
/// Mirrors the `FMT` field on each entry of Perl Prokka's `@database` list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DbFormat {
    /// Protein BLAST database (indexed via `makeblastdb`).
    Blast,
    /// HMMER3 profile database (pressed via `hmmpress`).
    Hmmer3,
}

/// Descriptor for an annotation database in the search hierarchy.
#[derive(Debug, Clone)]
pub struct AnnotationDb {
    /// Path to the database file.
    pub path: String,
    /// Source prefix for /inference tag (e.g. "similar to AA sequence:UniProtKB:").
    pub source_prefix: String,
    /// Database format.
    pub format: DbFormat,
    /// E-value threshold (overrides global if set).
    pub evalue: Option<f64>,
    /// Minimum query coverage percent (overrides global if set).
    pub min_coverage: Option<f64>,
}

/// Parsed annotation from a ~~~-delimited FASTA description.
#[derive(Debug, Clone, Default)]
pub struct ParsedAnnotation {
    pub ec_number: String,
    pub gene: String,
    pub product: String,
    pub cog: String,
}

/// Parse a ~~~-delimited annotation header.
///
/// Format: `EC_number~~~gene~~~product~~~COG`
/// Fields may be empty. If no ~~~ delimiter, the entire string is the product.
pub fn parse_annotation_header(desc: &str) -> ParsedAnnotation {
    if desc.contains("~~~") {
        let parts: Vec<&str> = desc.splitn(4, "~~~").collect();
        ParsedAnnotation {
            ec_number: collapse_ec(parts.first().unwrap_or(&"")),
            gene: parts.get(1).unwrap_or(&"").to_string(),
            product: parts.get(2).unwrap_or(&"").to_string(),
            cog: parts.get(3).unwrap_or(&"").to_string(),
        }
    } else {
        ParsedAnnotation {
            product: desc.to_string(),
            ..Default::default()
        }
    }
}

/// Collapse transitionary EC numbers: replace `n\d+` with `-`.
fn collapse_ec(ec: &str) -> String {
    static RE: std::sync::LazyLock<regex::Regex> = std::sync::LazyLock::new(|| {
        regex::Regex::new(r"n\d+").unwrap()
    });
    RE.replace_all(ec, "-").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_full_annotation() {
        let ann = parse_annotation_header("2.1.1.48~~~ermC~~~rRNA adenine N-6-methyltransferase~~~COG1234");
        assert_eq!(ann.ec_number, "2.1.1.48");
        assert_eq!(ann.gene, "ermC");
        assert_eq!(ann.product, "rRNA adenine N-6-methyltransferase");
        assert_eq!(ann.cog, "COG1234");
    }

    #[test]
    fn test_parse_partial_annotation() {
        let ann = parse_annotation_header("~~~traB~~~transfer complex protein TraB~~~");
        assert_eq!(ann.ec_number, "");
        assert_eq!(ann.gene, "traB");
        assert_eq!(ann.product, "transfer complex protein TraB");
        assert_eq!(ann.cog, "");
    }

    #[test]
    fn test_parse_minimal_annotation() {
        let ann = parse_annotation_header("~~~~~~transposase~~~");
        assert_eq!(ann.ec_number, "");
        assert_eq!(ann.gene, "");
        assert_eq!(ann.product, "transposase");
    }

    #[test]
    fn test_parse_plain_product() {
        let ann = parse_annotation_header("hypothetical protein");
        assert_eq!(ann.product, "hypothetical protein");
        assert_eq!(ann.ec_number, "");
    }

    #[test]
    fn test_ec_collapse() {
        assert_eq!(collapse_ec("2.1.n1.48"), "2.1.-.48");
        assert_eq!(collapse_ec("1.2.3.4"), "1.2.3.4");
    }
}
