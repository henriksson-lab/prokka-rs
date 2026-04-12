use std::collections::HashSet;
use std::sync::LazyLock;

use regex::Regex;

/// The product name assigned to unannotated proteins.
pub const HYPO: &str = "hypothetical protein";
/// The product name for --noanno mode.
pub const UNANN: &str = "unannotated protein";

// Pre-compiled regexes for cleanup_product (compiled once, used many times).
static RE_HYPO: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?i)DUF\d|UPF\d|conserved|domain of unknown|\b[CN].term|paralog").unwrap()
});
static RE_HOMOLOG: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\bhomolog( \d)?\b").unwrap()
});
static RE_ARCOG: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^arCOG\d+\s+").unwrap()
});
static RE_EC_COG: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\((EC|COG).*?\)").unwrap()
});
static RE_IS: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\bIS\d+\b").unwrap()
});
static RE_LOCUS: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\s+\w+\d{4,}c?").unwrap()
});
static RE_AND: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r" and (inactivated|related) \w+").unwrap()
});
static RE_FAMILY: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r",\s*family$").unwrap()
});
static RE_PUTATIVE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?i)^(potential|possible|probable|predicted|uncharacteri.ed)").unwrap()
});
static RE_SUFFIX: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?i)(domain|family|binding|fold|like|motif|repeat)\s*$").unwrap()
});
static RE_SPACE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\s+").unwrap()
});

/// Clean up a /product annotation to be Genbank/ENA compliant.
///
/// Exact port of Perl Prokka's `cleanup_product()` function (lines 1461-1490).
pub fn cleanup_product(product: &str, good_products: &HashSet<String>) -> String {
    // Check the whitelist first
    if good_products.contains(product) {
        return product.to_string();
    }

    if RE_HYPO.is_match(product) {
        return HYPO.to_string();
    }

    // Return HYPO if no lowercase letters at all
    if !product.chars().any(|c| c.is_ascii_lowercase()) {
        return HYPO.to_string();
    }

    let mut p = product.to_string();

    p = RE_HOMOLOG.replace_all(&p, "").to_string();
    p = RE_ARCOG.replace(&p, "").to_string();
    p = RE_EC_COG.replace_all(&p, "").to_string();

    if !RE_IS.is_match(&p) {
        p = RE_LOCUS.replace_all(&p, "").to_string();
    }

    p = RE_AND.replace_all(&p, "").to_string();
    p = RE_FAMILY.replace(&p, "").to_string();
    p = RE_PUTATIVE.replace(&p, "putative").to_string();

    if RE_SUFFIX.is_match(&p) && !p.contains(',') {
        p.push_str(" protein");
    }

    p = RE_SPACE.replace_all(&p, " ").to_string();

    p.trim().to_string()
}

/// Default whitelist of product names that should not be cleaned.
pub fn default_good_products() -> HashSet<String> {
    let mut set = HashSet::new();
    set.insert("rep".to_string());
    set.insert("Conserved virulence factor B".to_string());
    set.insert("Cypemycin N-terminal methyltransferase".to_string());
    set
}

#[cfg(test)]
mod tests {
    use super::*;

    fn good() -> HashSet<String> {
        default_good_products()
    }

    #[test]
    fn test_whitelist_passes_through() {
        assert_eq!(cleanup_product("rep", &good()), "rep");
        assert_eq!(
            cleanup_product("Conserved virulence factor B", &good()),
            "Conserved virulence factor B"
        );
    }

    #[test]
    fn test_duf_becomes_hypo() {
        assert_eq!(cleanup_product("DUF1234 domain protein", &good()), HYPO);
        assert_eq!(cleanup_product("UPF0123 family protein", &good()), HYPO);
    }

    #[test]
    fn test_no_lowercase_becomes_hypo() {
        assert_eq!(cleanup_product("ABCDEF", &good()), HYPO);
    }

    #[test]
    fn test_homolog_removed() {
        assert_eq!(
            cleanup_product("Leu-binding protein homolog", &good()),
            "Leu-binding protein"
        );
        assert_eq!(
            cleanup_product("ABC homolog 3 protein", &good()),
            "ABC protein"
        );
    }

    #[test]
    fn test_putative_normalization() {
        assert_eq!(
            cleanup_product("Probable methyltransferase", &good()),
            "putative methyltransferase"
        );
        assert_eq!(
            cleanup_product("Predicted kinase", &good()),
            "putative kinase"
        );
    }

    #[test]
    fn test_domain_gets_protein_suffix() {
        assert_eq!(
            cleanup_product("ABC-type domain", &good()),
            "ABC-type domain protein"
        );
    }

    #[test]
    fn test_conserved_becomes_hypo() {
        assert_eq!(cleanup_product("conserved hypothetical protein", &good()), HYPO);
    }
}
