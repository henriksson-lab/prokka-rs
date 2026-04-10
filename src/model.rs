use indexmap::IndexMap;
use std::fmt;

/// Strand of a genomic feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Returns +1 for Forward, -1 for Reverse (matching BioPerl convention).
    pub fn to_int(self) -> i8 {
        match self {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// Feature type for genomic annotations.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FeatureType {
    CDS,
    Gene,
    MRNA,
    TRNA,
    TMRNA,
    RRNA,
    MiscRNA,
    RepeatRegion,
    SigPeptide,
}

impl FeatureType {
    /// Returns the string representation matching Prokka output.
    pub fn as_str(&self) -> &'static str {
        match self {
            FeatureType::CDS => "CDS",
            FeatureType::Gene => "gene",
            FeatureType::MRNA => "mRNA",
            FeatureType::TRNA => "tRNA",
            FeatureType::TMRNA => "tmRNA",
            FeatureType::RRNA => "rRNA",
            FeatureType::MiscRNA => "misc_RNA",
            FeatureType::RepeatRegion => "repeat_region",
            FeatureType::SigPeptide => "sig_peptide",
        }
    }

    /// Whether this feature type is an RNA type (for overlap checking).
    pub fn is_rna(&self) -> bool {
        matches!(self, FeatureType::TRNA | FeatureType::TMRNA | FeatureType::RRNA)
    }
}

impl fmt::Display for FeatureType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl PartialOrd for FeatureType {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FeatureType {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.as_str().cmp(other.as_str())
    }
}

/// A single genomic feature annotation.
///
/// Replaces BioPerl's `Bio::SeqFeature::Generic`.
#[derive(Debug, Clone)]
pub struct SeqFeature {
    /// Feature type (CDS, tRNA, rRNA, etc.).
    pub feature_type: FeatureType,
    /// Contig/sequence ID this feature belongs to.
    pub seq_id: String,
    /// Source tool that predicted this feature.
    pub source: String,
    /// 1-based start position (inclusive).
    pub start: usize,
    /// 1-based end position (inclusive).
    pub end: usize,
    /// Strand.
    pub strand: Strand,
    /// Score (e-value, probability, etc.) — None if not applicable.
    pub score: Option<f64>,
    /// GFF frame: 0 for CDS, -1 for others (written as '.').
    pub frame: i8,
    /// Key-value tags (preserves insertion order).
    /// A key can have multiple values.
    pub tags: IndexMap<String, Vec<String>>,
}

impl SeqFeature {
    /// Create a new feature with the given core attributes and no tags.
    pub fn new(
        feature_type: FeatureType,
        seq_id: String,
        source: String,
        start: usize,
        end: usize,
        strand: Strand,
    ) -> Self {
        let frame = if feature_type == FeatureType::CDS { 0 } else { -1 };
        SeqFeature {
            feature_type,
            seq_id,
            source,
            start,
            end,
            strand,
            score: None,
            frame,
            tags: IndexMap::new(),
        }
    }

    /// Add a tag value. Multiple values for the same key are allowed.
    pub fn add_tag(&mut self, key: &str, value: &str) {
        self.tags
            .entry(key.to_string())
            .or_default()
            .push(value.to_string());
    }

    /// Get the first value of a tag, or None.
    pub fn get_tag(&self, key: &str) -> Option<&str> {
        self.tags.get(key).and_then(|v| v.first()).map(|s| s.as_str())
    }

    /// Check if a tag exists.
    pub fn has_tag(&self, key: &str) -> bool {
        self.tags.contains_key(key)
    }

    /// Remove a tag entirely.
    pub fn remove_tag(&mut self, key: &str) {
        self.tags.swap_remove(key);
    }

    /// Feature length in base pairs.
    pub fn length(&self) -> usize {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }

    /// Check if this feature overlaps with another on the same contig.
    pub fn overlaps(&self, other: &SeqFeature) -> bool {
        self.seq_id == other.seq_id && self.start <= other.end && self.end >= other.start
    }
}

/// A contig (DNA sequence) with its annotations.
///
/// Replaces BioPerl's `Bio::Seq` + associated features.
#[derive(Debug, Clone)]
pub struct Contig {
    /// Sequence identifier.
    pub id: String,
    /// DNA sequence (uppercase ACGTN only).
    pub sequence: Vec<u8>,
    /// Annotated features on this contig.
    pub features: Vec<SeqFeature>,
}

impl Contig {
    /// Create a new contig with no features.
    pub fn new(id: String, sequence: Vec<u8>) -> Self {
        Contig {
            id,
            sequence,
            features: Vec::new(),
        }
    }

    /// Sequence length in base pairs.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Whether the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Extract a subsequence (1-based, inclusive start and end).
    /// Returns the subsequence on the given strand (reverse-complemented if Reverse).
    pub fn extract(&self, start: usize, end: usize, strand: Strand) -> Vec<u8> {
        if start < 1 || end > self.sequence.len() || start > end {
            return Vec::new();
        }
        let subseq = &self.sequence[start - 1..end];
        match strand {
            Strand::Forward => subseq.to_vec(),
            Strand::Reverse => reverse_complement(subseq),
        }
    }
}

/// Reverse complement a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'N' => b'N',
            _ => b'N',
        })
        .collect()
}

/// Statistics about an annotation run.
#[derive(Debug, Clone, Default)]
pub struct AnnotationStats {
    /// Count of features by type.
    pub feature_counts: IndexMap<String, usize>,
    /// Total contigs.
    pub num_contigs: usize,
    /// Total base pairs.
    pub total_bp: usize,
}

/// The complete result of an annotation pipeline run.
#[derive(Debug, Clone)]
pub struct AnnotationResult {
    /// Annotated contigs with their features.
    pub contigs: Vec<Contig>,
    /// Summary statistics.
    pub stats: AnnotationStats,
    /// Log messages generated during the run.
    pub log: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GATTACA"), b"TGTAATC");
        assert_eq!(reverse_complement(b"N"), b"N");
    }

    #[test]
    fn test_contig_extract() {
        let c = Contig::new("test".into(), b"ACGTACGT".to_vec());
        assert_eq!(c.extract(1, 4, Strand::Forward), b"ACGT");
        assert_eq!(c.extract(1, 4, Strand::Reverse), b"ACGT"); // palindrome
        assert_eq!(c.extract(1, 3, Strand::Forward), b"ACG");
        assert_eq!(c.extract(1, 3, Strand::Reverse), b"CGT");
    }

    #[test]
    fn test_feature_overlap() {
        let f1 = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "test".into(), 100, 200, Strand::Forward,
        );
        let f2 = SeqFeature::new(
            FeatureType::TRNA, "ctg1".into(), "test".into(), 150, 250, Strand::Forward,
        );
        let f3 = SeqFeature::new(
            FeatureType::RRNA, "ctg1".into(), "test".into(), 300, 400, Strand::Forward,
        );
        let f4 = SeqFeature::new(
            FeatureType::CDS, "ctg2".into(), "test".into(), 150, 250, Strand::Forward,
        );
        assert!(f1.overlaps(&f2));
        assert!(!f1.overlaps(&f3));
        assert!(!f1.overlaps(&f4)); // different contig
    }

    #[test]
    fn test_feature_tags() {
        let mut f = SeqFeature::new(
            FeatureType::CDS, "ctg1".into(), "test".into(), 1, 100, Strand::Forward,
        );
        f.add_tag("product", "hypothetical protein");
        f.add_tag("gene", "abc");
        assert_eq!(f.get_tag("product"), Some("hypothetical protein"));
        assert_eq!(f.get_tag("gene"), Some("abc"));
        assert!(f.has_tag("product"));
        assert!(!f.has_tag("EC_number"));
        f.remove_tag("gene");
        assert!(!f.has_tag("gene"));
    }
}
