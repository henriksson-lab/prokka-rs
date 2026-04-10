use std::collections::HashSet;

use regex::Regex;

/// The product name assigned to unannotated proteins.
pub const HYPO: &str = "hypothetical protein";
/// The product name for --noanno mode.
pub const UNANN: &str = "unannotated protein";

/// Clean up a /product annotation to be Genbank/ENA compliant.
///
/// Exact port of Perl Prokka's `cleanup_product()` function (lines 1461-1490).
pub fn cleanup_product(product: &str, good_products: &HashSet<String>) -> String {
    // Check the whitelist first
    if good_products.contains(product) {
        return product.to_string();
    }

    let p = product.to_string();

    // Return HYPO if matches specific patterns
    let hypo_re = Regex::new(r"(?i)DUF\d|UPF\d|conserved|domain of unknown|\b[CN].term|paralog").unwrap();
    if hypo_re.is_match(&p) {
        return HYPO.to_string();
    }

    // Return HYPO if no lowercase letters at all
    if !p.chars().any(|c| c.is_ascii_lowercase()) {
        return HYPO.to_string();
    }

    let mut p = p;

    // Remove "homolog" (with optional " N")
    let homolog_re = Regex::new(r"\bhomolog( \d)?\b").unwrap();
    p = homolog_re.replace_all(&p, "").to_string();

    // Remove arCOG prefix
    let arcog_re = Regex::new(r"^arCOG\d+\s+").unwrap();
    p = arcog_re.replace(&p, "").to_string();

    // Remove (EC...) and (COG...)
    let ec_cog_re = Regex::new(r"\((EC|COG).*?\)").unwrap();
    p = ec_cog_re.replace_all(&p, "").to_string();

    // Remove possible locus tags (word+4+ digits), but not IS names
    let is_re = Regex::new(r"\bIS\d+\b").unwrap();
    if !is_re.is_match(&p) {
        let locus_re = Regex::new(r"\s+\w+\d{4,}c?").unwrap();
        p = locus_re.replace_all(&p, "").to_string();
    }

    // Remove "and (inactivated|related) word"
    let and_re = Regex::new(r" and (inactivated|related) \w+").unwrap();
    p = and_re.replace_all(&p, "").to_string();

    // Remove trailing ", family"
    let family_re = Regex::new(r",\s*family$").unwrap();
    p = family_re.replace(&p, "").to_string();

    // Replace potential/possible/probable/predicted/uncharacterized with "putative"
    let putative_re = Regex::new(r"(?i)^(potential|possible|probable|predicted|uncharacteri.ed)").unwrap();
    p = putative_re.replace(&p, "putative").to_string();

    // If ends with domain/family/binding/fold/like/motif/repeat, append " protein"
    let suffix_re = Regex::new(r"(?i)(domain|family|binding|fold|like|motif|repeat)\s*$").unwrap();
    if suffix_re.is_match(&p) && !p.contains(',') {
        p.push_str(" protein");
    }

    // Collapse multiple spaces
    let space_re = Regex::new(r"\s+").unwrap();
    p = space_re.replace_all(&p, " ").to_string();

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
