//! Runtime configuration for a Prokka annotation run.
//!
//! [`ProkkaConfig`] mirrors every CLI flag from Perl Prokka's `setOptions`
//! (see `prokka/bin/prokka` around line 1734) so the library and the CLI binary
//! share the same options surface. The [`Kingdom`] enum captures the
//! `--kingdom` annotation mode and its kingdom-specific defaults (genetic
//! code, Barrnap mode, Aragorn flags, allowed feature overlaps).

use std::path::PathBuf;

/// Kingdom/annotation mode (`--kingdom`).
///
/// Determines the default genetic code, which rRNA tool flags to use, and
/// whether CDS/RNA overlap is tolerated. Parsed case-insensitively from the
/// same prefixes Perl Prokka accepts.
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
/// One field per CLI flag (see `sub setOptions` and `sub usage` in
/// `prokka/bin/prokka`). Use [`ProkkaConfig::default`] for sensible defaults,
/// then mutate the fields you care about. [`ProkkaConfig::apply_compliant`]
/// applies the `--compliant` preset and [`ProkkaConfig::validate`] checks
/// invariants before the pipeline runs.
#[derive(Debug, Clone)]
pub struct ProkkaConfig {
    // General
    /// `--quiet`: suppress screen output.
    pub quiet: bool,
    /// `--debug`: keep temporary files for inspection.
    pub debug: bool,

    // Setup
    /// `--dbdir`: Prokka database root folder.
    pub dbdir: PathBuf,
    /// `--listdb`: list all configured databases and exit.
    pub listdb: bool,
    /// `--setupdb`: index all installed databases and exit.
    pub setupdb: bool,
    /// `--cleandb`: remove all database indices and exit.
    pub cleandb: bool,

    // Outputs
    /// `--outdir`: output folder (auto-derived from prefix if `None`).
    pub outdir: Option<PathBuf>,
    /// `--force`: overwrite an existing output folder.
    pub force: bool,
    /// `--prefix`: filename prefix for all outputs (auto-derived if `None`).
    pub prefix: Option<String>,
    /// `--addgenes`: add a `gene` feature for each `CDS` feature.
    pub addgenes: bool,
    /// `--addmrna`: add an `mRNA` feature for each `CDS` feature.
    pub addmrna: bool,
    /// `--locustag`: locus tag prefix (MD5-derived if `None`; see
    /// [`crate::locus_tag::generate_locus_tag`]).
    pub locustag: Option<String>,
    /// `--increment`: locus tag counter increment.
    pub increment: u32,
    /// `--gffver`: GFF version to write (2 or 3).
    pub gffver: u8,
    /// `--compliant`: enforce GenBank/ENA/DDBJ submission rules
    /// (see [`ProkkaConfig::apply_compliant`]).
    pub compliant: bool,
    /// `--centre`: sequencing centre ID for compliant contig renaming.
    pub centre: Option<String>,
    /// `--accver`: accession version to put in the GenBank file.
    pub accver: u32,

    // Organism details
    /// `--genus`: organism genus name.
    pub genus: String,
    /// `--species`: organism species name.
    pub species: String,
    /// `--strain`: organism strain name.
    pub strain: String,
    /// `--plasmid`: plasmid name or identifier.
    pub plasmid: Option<String>,

    // Annotations
    /// `--kingdom`: annotation mode.
    pub kingdom: Kingdom,
    /// `--gcode`: NCBI translation table (`0` = use kingdom default,
    /// resolved by [`ProkkaConfig::effective_gcode`]).
    pub gcode: u8,
    /// `--prodigaltf`: Prodigal training file (overrides built-in model).
    pub prodigaltf: Option<PathBuf>,
    /// `--gram`: Gram stain (`-`/`neg` or `+`/`pos`); enables SignalP.
    pub gram: Option<String>,
    /// `--usegenus`: also search the genus-specific BLAST database.
    pub usegenus: bool,
    /// `--proteins`: extra trusted protein FASTA/GenBank/EMBL searched first.
    pub proteins: Option<PathBuf>,
    /// `--hmms`: extra trusted HMM profile DB searched first.
    pub hmms: Option<PathBuf>,
    /// `--metagenome`: improve gene predictions for fragmented genomes
    /// (Prodigal "meta" mode).
    pub metagenome: bool,
    /// `--rawproduct`: skip the `/product` cleanup step.
    pub rawproduct: bool,
    /// `--cdsrnaolap`: allow tRNA/rRNA features to overlap CDS.
    pub cdsrnaolap: bool,

    // Matching
    /// `--evalue`: BLAST/HMMER similarity e-value cut-off.
    pub evalue: f64,
    /// `--coverage`: minimum percent coverage on query protein.
    pub coverage: f64,

    // Computation
    /// `--cpus`: thread count (`0` = use all cores).
    pub cpus: u32,
    /// `--fast`: skip HMM databases, BLASTP only.
    pub fast: bool,
    /// `--noanno`: skip CDS annotation; label all CDS as `unannotated protein`.
    pub noanno: bool,
    /// `--mincontiglen`: drop contigs shorter than this (NCBI requires 200).
    pub mincontiglen: u32,
    /// `--rfam`: enable ncRNA search with Infernal/Rfam (slow).
    pub rfam: bool,
    /// `--norrna`: disable rRNA prediction.
    pub norrna: bool,
    /// `--notrna`: disable tRNA prediction.
    pub notrna: bool,
    /// `--rnammer`: prefer RNAmmer over Barrnap for rRNA prediction.
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
        assert_eq!(
            Kingdom::parse("Mitochondria").unwrap(),
            Kingdom::Mitochondria
        );
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
