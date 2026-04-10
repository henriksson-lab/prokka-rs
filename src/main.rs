use std::path::PathBuf;
use clap::Parser;

/// Rapid prokaryotic genome annotation.
#[derive(Parser, Debug)]
#[command(name = "prokka-rs", version, about)]
struct Cli {
    /// Input contig FASTA file
    input: Option<PathBuf>,

    // General
    /// No screen output
    #[arg(long)]
    quiet: bool,
    /// Debug mode: keep all temporary files
    #[arg(long)]
    debug: bool,

    /// Print citation for referencing prokka-rs
    #[arg(long)]
    citation: bool,
    /// List all software dependencies
    #[arg(long)]
    depends: bool,

    // Setup
    /// Prokka database root folder
    #[arg(long)]
    dbdir: Option<PathBuf>,
    /// List all configured databases
    #[arg(long)]
    listdb: bool,
    /// Index all installed databases
    #[arg(long)]
    setupdb: bool,
    /// Remove all database indices
    #[arg(long)]
    cleandb: bool,

    // Outputs
    /// Output folder
    #[arg(long)]
    outdir: Option<PathBuf>,
    /// Force overwriting existing output folder
    #[arg(long)]
    force: bool,
    /// Filename output prefix
    #[arg(long)]
    prefix: Option<String>,
    /// Add 'gene' features for each 'CDS' feature
    #[arg(long)]
    addgenes: bool,
    /// Add 'mRNA' features for each 'CDS' feature
    #[arg(long)]
    addmrna: bool,
    /// Locus tag prefix
    #[arg(long)]
    locustag: Option<String>,
    /// Locus tag counter increment
    #[arg(long, default_value = "1")]
    increment: u32,
    /// GFF version
    #[arg(long, default_value = "3")]
    gffver: u8,
    /// Force Genbank/ENA/DDJB compliance
    #[arg(long)]
    compliant: bool,
    /// Sequencing centre ID
    #[arg(long)]
    centre: Option<String>,
    /// Version to put in Genbank file
    #[arg(long, default_value = "1")]
    accver: u32,

    // Organism details
    /// Genus name
    #[arg(long, default_value = "Genus")]
    genus: String,
    /// Species name
    #[arg(long, default_value = "species")]
    species: String,
    /// Strain name
    #[arg(long, default_value = "strain")]
    strain: String,
    /// Plasmid name or identifier
    #[arg(long)]
    plasmid: Option<String>,

    // Annotations
    /// Annotation mode: Bacteria|Archaea|Viruses|Mitochondria
    #[arg(long, default_value = "Bacteria")]
    kingdom: String,
    /// Genetic code / Translation table
    #[arg(long, default_value = "0")]
    gcode: u8,
    /// Prodigal training file
    #[arg(long)]
    prodigaltf: Option<PathBuf>,
    /// Gram: -/neg +/pos
    #[arg(long)]
    gram: Option<String>,
    /// Use genus-specific BLAST databases (needs --genus)
    #[arg(long)]
    usegenus: bool,
    /// FASTA or GBK file to use as 1st priority
    #[arg(long)]
    proteins: Option<PathBuf>,
    /// Trusted HMM to first annotate from
    #[arg(long)]
    hmms: Option<PathBuf>,
    /// Improve gene predictions for highly fragmented genomes
    #[arg(long)]
    metagenome: bool,
    /// Do not clean up /product annotation
    #[arg(long)]
    rawproduct: bool,
    /// Allow [tr]RNA to overlap CDS
    #[arg(long)]
    cdsrnaolap: bool,

    // Matching
    /// Similarity e-value cut-off
    #[arg(long, default_value = "1e-9")]
    evalue: f64,
    /// Minimum coverage on query protein
    #[arg(long, default_value = "80")]
    coverage: f64,

    // Computation
    /// Number of CPUs to use [0=all]
    #[arg(long, default_value = "8")]
    cpus: u32,
    /// Fast mode - only use basic BLASTP databases
    #[arg(long)]
    fast: bool,
    /// For CDS just set /product="unannotated protein"
    #[arg(long)]
    noanno: bool,
    /// Minimum contig size [NCBI needs 200]
    #[arg(long, default_value = "1")]
    mincontiglen: u32,
    /// Enable searching for ncRNAs with Infernal+Rfam (SLOW!)
    #[arg(long)]
    rfam: bool,
    /// Don't run rRNA search
    #[arg(long)]
    norrna: bool,
    /// Don't run tRNA search
    #[arg(long)]
    notrna: bool,
    /// Prefer RNAmmer over Barrnap for rRNA prediction
    #[arg(long)]
    rnammer: bool,
}

impl Cli {
    fn to_config(&self) -> Result<prokka_rs::ProkkaConfig, prokka_rs::ProkkaError> {
        let kingdom = prokka_rs::config::Kingdom::parse(&self.kingdom)?;
        let default_dbdir = std::path::PathBuf::from(
            std::env::var("PROKKA_DBDIR").unwrap_or_else(|_| "db".to_string()),
        );

        Ok(prokka_rs::ProkkaConfig {
            quiet: self.quiet,
            debug: self.debug,
            dbdir: self.dbdir.clone().unwrap_or(default_dbdir),
            listdb: self.listdb,
            setupdb: self.setupdb,
            cleandb: self.cleandb,
            outdir: self.outdir.clone(),
            force: self.force,
            prefix: self.prefix.clone(),
            addgenes: self.addgenes,
            addmrna: self.addmrna,
            locustag: self.locustag.clone(),
            increment: self.increment,
            gffver: self.gffver,
            compliant: self.compliant,
            centre: self.centre.clone(),
            accver: self.accver,
            genus: self.genus.clone(),
            species: self.species.clone(),
            strain: self.strain.clone(),
            plasmid: self.plasmid.clone(),
            kingdom,
            gcode: self.gcode,
            prodigaltf: self.prodigaltf.clone(),
            gram: self.gram.clone(),
            usegenus: self.usegenus,
            proteins: self.proteins.clone(),
            hmms: self.hmms.clone(),
            metagenome: self.metagenome,
            rawproduct: self.rawproduct,
            cdsrnaolap: self.cdsrnaolap,
            evalue: self.evalue,
            coverage: self.coverage,
            cpus: self.cpus,
            fast: self.fast,
            noanno: self.noanno,
            mincontiglen: self.mincontiglen,
            rfam: self.rfam,
            norrna: self.norrna,
            notrna: self.notrna,
            rnammer: self.rnammer,
        })
    }
}

fn main() {
    let cli = Cli::parse();

    // Handle info commands
    if cli.citation {
        eprintln!("\nIf you use prokka-rs in your work, please cite:\n");
        eprintln!("  Seemann T, \"Prokka: Rapid Prokaryotic Genome Annotation\",");
        eprintln!("  Bioinformatics, 2014 Jul 15;30(14):2068-9.\n");
        eprintln!("  PMID:24642063");
        eprintln!("  doi:10.1093/bioinformatics/btu153");
        eprintln!("  http://www.ncbi.nlm.nih.gov/pubmed/24642063\n");
        return;
    }

    if cli.depends {
        // Native Rust (no external binary needed)
        eprintln!("prodigal-rs (native, built-in)");
        eprintln!("blast-rs (native, built-in)");
        eprintln!("hmmer-pure-rs (native, built-in)");
        // External tools (optional)
        eprintln!("aragorn >= 1.2 (optional, tRNA prediction)");
        eprintln!("barrnap >= 0.4 (optional, rRNA prediction)");
        eprintln!("minced >= 2.0 (optional, CRISPR detection)");
        eprintln!("cmscan >= 1.1 (optional, ncRNA with --rfam)");
        eprintln!("signalp >= 3.0 (optional, signal peptides with --gram)");
        eprintln!("tbl2asn (optional, GenBank/Sequin output)");
        return;
    }

    // Handle setup commands that don't need input
    if cli.listdb || cli.setupdb || cli.cleandb {
        let config = match cli.to_config() {
            Ok(c) => c,
            Err(e) => {
                eprintln!("Error: {}", e);
                std::process::exit(2);
            }
        };
        if cli.cleandb {
            if let Err(e) = prokka_rs::database::clean::clean_db(&config.dbdir) {
                eprintln!("Error: {}", e);
                std::process::exit(2);
            }
        }
        if cli.setupdb {
            if let Err(e) = prokka_rs::database::setup::setup_db(&config.dbdir) {
                eprintln!("Error: {}", e);
                std::process::exit(2);
            }
        }
        if cli.listdb {
            let inv = prokka_rs::database::list::list_db(&config.dbdir);
            eprintln!("Kingdoms: {:?}", inv.kingdoms);
            eprintln!("Genera: {:?}", inv.genera);
            eprintln!("HMMs: {:?}", inv.hmms);
            eprintln!("CMs: {:?}", inv.cms);
        }
        return;
    }

    let input = match &cli.input {
        Some(p) => p.clone(),
        None => {
            eprintln!("Error: Please supply a contig FASTA file on the command line.");
            std::process::exit(2);
        }
    };

    let mut config = match cli.to_config() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(2);
        }
    };

    if let Err(e) = prokka_rs::pipeline::run(&input, &mut config) {
        eprintln!("Error: {}", e);
        std::process::exit(2);
    }
}
