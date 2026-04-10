use std::path::PathBuf;

/// Kingdom/annotation mode for Prokka.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Kingdom {
    Bacteria,
    Archaea,
    Viruses,
    Mitochondria,
}

impl Kingdom {
    /// Parse from string (case-insensitive prefix matching like Perl Prokka).
    pub fn parse(s: &str) -> Result<Self, crate::error::ProkkaError> {
        let lower = s.to_lowercase();
        if lower.starts_with("bac") || lower.starts_with("prok") {
            Ok(Kingdom::Bacteria)
        } else if lower.starts_with("arch") {
            Ok(Kingdom::Archaea)
        } else if lower.starts_with("vir") {
            Ok(Kingdom::Viruses)
        } else if lower.starts_with("mito") || lower.starts_with("mt") {
            Ok(Kingdom::Mitochondria)
        } else {
            Err(crate::error::ProkkaError::InvalidKingdom(s.to_string()))
        }
    }

    /// Display name as used in Prokka output.
    pub fn as_str(&self) -> &'static str {
        match self {
            Kingdom::Bacteria => "Bacteria",
            Kingdom::Archaea => "Archaea",
            Kingdom::Viruses => "Viruses",
            Kingdom::Mitochondria => "Mitochondria",
        }
    }

    /// Default genetic code for this kingdom.
    pub fn default_gcode(&self) -> u8 {
        match self {
            Kingdom::Bacteria => 11,
            Kingdom::Archaea => 11,
            Kingdom::Viruses => 1,
            Kingdom::Mitochondria => 5,
        }
    }

    /// Barrnap kingdom mode string, or None if not applicable.
    pub fn barrnap_mode(&self) -> Option<&'static str> {
        match self {
            Kingdom::Bacteria => Some("bac"),
            Kingdom::Archaea => Some("arc"),
            Kingdom::Mitochondria => Some("mito"),
            Kingdom::Viruses => None,
        }
    }

    /// Aragorn extra options for this kingdom.
    pub fn aragorn_opt(&self) -> &'static str {
        match self {
            Kingdom::Mitochondria => "-mt",
            _ => "",
        }
    }

    /// Whether CDS/RNA overlap is allowed (Mitochondria are highly packed).
    pub fn allow_cds_rna_overlap(&self) -> bool {
        matches!(self, Kingdom::Mitochondria)
    }
}

/// Full configuration for a Prokka annotation run.
///
/// Maps all CLI parameters from the original Perl Prokka.
#[derive(Debug, Clone)]
pub struct ProkkaConfig {
    // General
    pub quiet: bool,
    pub debug: bool,

    // Setup
    pub dbdir: PathBuf,
    pub listdb: bool,
    pub setupdb: bool,
    pub cleandb: bool,

    // Outputs
    pub outdir: Option<PathBuf>,
    pub force: bool,
    pub prefix: Option<String>,
    pub addgenes: bool,
    pub addmrna: bool,
    pub locustag: Option<String>,
    pub increment: u32,
    pub gffver: u8,
    pub compliant: bool,
    pub centre: Option<String>,
    pub accver: u32,

    // Organism details
    pub genus: String,
    pub species: String,
    pub strain: String,
    pub plasmid: Option<String>,

    // Annotations
    pub kingdom: Kingdom,
    pub gcode: u8,
    pub prodigaltf: Option<PathBuf>,
    pub gram: Option<String>,
    pub usegenus: bool,
    pub proteins: Option<PathBuf>,
    pub hmms: Option<PathBuf>,
    pub metagenome: bool,
    pub rawproduct: bool,
    pub cdsrnaolap: bool,

    // Matching
    pub evalue: f64,
    pub coverage: f64,

    // Computation
    pub cpus: u32,
    pub fast: bool,
    pub noanno: bool,
    pub mincontiglen: u32,
    pub rfam: bool,
    pub norrna: bool,
    pub notrna: bool,
    pub rnammer: bool,
}

impl Default for ProkkaConfig {
    fn default() -> Self {
        ProkkaConfig {
            quiet: false,
            debug: false,
            dbdir: PathBuf::from("db"),
            listdb: false,
            setupdb: false,
            cleandb: false,
            outdir: None,
            force: false,
            prefix: None,
            addgenes: false,
            addmrna: false,
            locustag: None,
            increment: 1,
            gffver: 3,
            compliant: false,
            centre: None,
            accver: 1,
            genus: "Genus".to_string(),
            species: "species".to_string(),
            strain: "strain".to_string(),
            plasmid: None,
            kingdom: Kingdom::Bacteria,
            gcode: 0, // 0 means "use kingdom default"
            prodigaltf: None,
            gram: None,
            usegenus: false,
            proteins: None,
            hmms: None,
            metagenome: false,
            rawproduct: false,
            cdsrnaolap: false,
            evalue: 1e-9,
            coverage: 80.0,
            cpus: 8,
            fast: false,
            noanno: false,
            mincontiglen: 1,
            rfam: false,
            norrna: false,
            notrna: false,
            rnammer: false,
        }
    }
}

impl ProkkaConfig {
    /// Resolve the effective genetic code (apply kingdom default if gcode == 0).
    pub fn effective_gcode(&self) -> u8 {
        if self.gcode == 0 {
            self.kingdom.default_gcode()
        } else {
            self.gcode
        }
    }

    /// Apply --compliant mode adjustments.
    pub fn apply_compliant(&mut self) {
        if self.compliant {
            self.addgenes = true;
            if self.centre.is_none() {
                self.centre = Some("Prokka".to_string());
            }
            self.mincontiglen = 200;
        }
    }

    /// Validate configuration, returning an error if invalid.
    pub fn validate(&self) -> Result<(), crate::error::ProkkaError> {
        let gcode = self.effective_gcode();
        if !(1..=25).contains(&gcode) {
            return Err(crate::error::ProkkaError::InvalidGeneticCode(gcode));
        }
        if self.evalue < 0.0 {
            return Err(crate::error::ProkkaError::InvalidEvalue(self.evalue));
        }
        if self.coverage < 0.0 || self.coverage > 100.0 {
            return Err(crate::error::ProkkaError::InvalidCoverage(self.coverage));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kingdom_parse() {
        assert_eq!(Kingdom::parse("Bacteria").unwrap(), Kingdom::Bacteria);
        assert_eq!(Kingdom::parse("bac").unwrap(), Kingdom::Bacteria);
        assert_eq!(Kingdom::parse("prok").unwrap(), Kingdom::Bacteria);
        assert_eq!(Kingdom::parse("Archaea").unwrap(), Kingdom::Archaea);
        assert_eq!(Kingdom::parse("arch").unwrap(), Kingdom::Archaea);
        assert_eq!(Kingdom::parse("Viruses").unwrap(), Kingdom::Viruses);
        assert_eq!(Kingdom::parse("vir").unwrap(), Kingdom::Viruses);
        assert_eq!(Kingdom::parse("Mitochondria").unwrap(), Kingdom::Mitochondria);
        assert_eq!(Kingdom::parse("mito").unwrap(), Kingdom::Mitochondria);
        assert_eq!(Kingdom::parse("mt").unwrap(), Kingdom::Mitochondria);
        assert!(Kingdom::parse("invalid").is_err());
    }

    #[test]
    fn test_kingdom_default_gcode() {
        assert_eq!(Kingdom::Bacteria.default_gcode(), 11);
        assert_eq!(Kingdom::Archaea.default_gcode(), 11);
        assert_eq!(Kingdom::Viruses.default_gcode(), 1);
        assert_eq!(Kingdom::Mitochondria.default_gcode(), 5);
    }

    #[test]
    fn test_effective_gcode() {
        let config = ProkkaConfig {
            gcode: 0,
            kingdom: Kingdom::Bacteria,
            ..Default::default()
        };
        assert_eq!(config.effective_gcode(), 11); // default for bacteria

        let config = ProkkaConfig {
            gcode: 4,
            kingdom: Kingdom::Bacteria,
            ..Default::default()
        };
        assert_eq!(config.effective_gcode(), 4); // explicit override
    }

    #[test]
    fn test_compliant_mode() {
        let mut config = ProkkaConfig::default();
        config.compliant = true;
        config.apply_compliant();
        assert!(config.addgenes);
        assert_eq!(config.centre, Some("Prokka".to_string()));
        assert_eq!(config.mincontiglen, 200);
    }

    #[test]
    fn test_compliant_preserves_existing_centre() {
        let mut config = ProkkaConfig {
            compliant: true,
            centre: Some("UoN".into()),
            ..Default::default()
        };
        config.apply_compliant();
        assert_eq!(config.centre, Some("UoN".to_string())); // not overwritten
    }

    #[test]
    fn test_validate_invalid_gcode() {
        let config = ProkkaConfig {
            gcode: 30,
            ..Default::default()
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn test_validate_invalid_evalue() {
        let config = ProkkaConfig {
            evalue: -1.0,
            ..Default::default()
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn test_validate_invalid_coverage() {
        let config = ProkkaConfig {
            coverage: 150.0,
            ..Default::default()
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn test_validate_ok() {
        let config = ProkkaConfig::default();
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_barrnap_mode() {
        assert_eq!(Kingdom::Bacteria.barrnap_mode(), Some("bac"));
        assert_eq!(Kingdom::Archaea.barrnap_mode(), Some("arc"));
        assert_eq!(Kingdom::Mitochondria.barrnap_mode(), Some("mito"));
        assert_eq!(Kingdom::Viruses.barrnap_mode(), None);
    }

    #[test]
    fn test_aragorn_opt() {
        assert_eq!(Kingdom::Bacteria.aragorn_opt(), "");
        assert_eq!(Kingdom::Mitochondria.aragorn_opt(), "-mt");
    }
}
