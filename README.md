# prokka-rs

Rapid prokaryotic genome annotation — [Prokka](https://github.com/tseemann/prokka) rewritten in Rust.

prokka-rs annotates bacterial, archaeal, viral, and mitochondrial genomes by identifying features
(tRNA, rRNA, ncRNA, CRISPRs, CDS) and assigning functional labels via hierarchical database searches.

## Features

- **Native gene prediction** via [prodigal-rs](https://github.com/henriksson-lab/prodigal-rs) (no external Prodigal binary needed)
- **Native BLAST searches** via [blast-rs](https://github.com/henriksson-lab/blast-rs)
- **Native HMM searches** via [hmmer-pure-rs](https://github.com/henriksson-lab/hmmer-pure-rs)
- **Library API** — use prokka-rs as a crate in your Rust code with input/output via variables
- All output files match the original Perl Prokka format

## Installation

```bash
cargo install prokka-rs
```

Or build from source:
```bash
git clone https://github.com/henriksson-lab/prokka-rs
cd prokka-rs
cargo build --release
```

## Quick Start

### CLI Usage

```bash
# Basic annotation
prokka-rs contigs.fna

# With custom output directory and organism details
prokka-rs --outdir mydir --prefix mygenome --genus Escherichia --species coli contigs.fna

# Fast mode (BLAST only, skip HMM searches)
prokka-rs --fast --cpus 8 contigs.fna

# Specify database directory
prokka-rs --dbdir /path/to/prokka/db contigs.fna
```

### Library Usage

```rust
use prokka_rs::{ProkkaConfig, Contig, annotate};
use std::path::Path;

// Load your sequences
let contigs = vec![
    Contig::new("seq1".to_string(), b"ATGAAACCC...".to_vec()),
];

let config = ProkkaConfig {
    kingdom: prokka_rs::config::Kingdom::Bacteria,
    ..Default::default()
};

let result = annotate(contigs, &config, None).unwrap();

for contig in &result.contigs {
    for feature in &contig.features {
        println!("{}\t{}\t{}\t{}",
            feature.seq_id,
            feature.feature_type,
            feature.start,
            feature.get_tag("product").unwrap_or(""));
    }
}
```

## Output Files

| Extension | Description |
|-----------|-------------|
| `.gff`    | Master annotation in GFF3 format (features + sequences) |
| `.fna`    | Nucleotide FASTA of input contigs |
| `.faa`    | Protein FASTA of translated CDS |
| `.ffn`    | Nucleotide FASTA of predicted transcripts |
| `.tbl`    | NCBI feature table |
| `.tsv`    | Tab-separated feature summary |
| `.fsa`    | Annotated FASTA for GenBank submission |
| `.gbk`    | GenBank flat file (requires tbl2asn) |
| `.txt`    | Annotation statistics |
| `.log`    | Execution log |

## Benchmarks

Measured on Linux x86_64, compiled with `-C target-cpu=native`:

| Test input | Size | Mode | Time |
|-----------|------|------|------|
| plasmid.fna | 57 KB | --noanno | 0.4s |
| genome.fna | 6.7 MB | --noanno | 18s |

## Command Line Options

```
prokka-rs [OPTIONS] <INPUT>

General:
  --quiet       No screen output
  --debug       Debug mode: keep temporary files

Setup:
  --dbdir       Prokka database root folder
  --listdb      List all configured databases
  --setupdb     Index all installed databases
  --cleandb     Remove all database indices

Outputs:
  --outdir      Output folder
  --force       Force overwriting existing output folder
  --prefix      Filename output prefix
  --addgenes    Add 'gene' features for each 'CDS' feature
  --locustag    Locus tag prefix
  --compliant   Force Genbank/ENA/DDJB compliance
  --centre      Sequencing centre ID

Organism:
  --genus       Genus name [default: Genus]
  --species     Species name [default: species]
  --strain      Strain name [default: strain]
  --kingdom     Bacteria|Archaea|Viruses|Mitochondria [default: Bacteria]
  --gcode       Genetic code / Translation table

Annotations:
  --proteins    FASTA file of trusted proteins (1st priority)
  --hmms        Trusted HMM file (1st priority)
  --usegenus    Use genus-specific BLAST databases
  --metagenome  Improve predictions for fragmented genomes
  --rawproduct  Do not clean up /product annotation
  --evalue      Similarity e-value cut-off [default: 1e-9]
  --coverage    Minimum query coverage [default: 80]

Computation:
  --cpus        Number of CPUs [default: 8]
  --fast        Fast mode — skip HMM searches
  --noanno      Skip all annotation searches
  --rfam        Enable ncRNA search with Infernal+Rfam (slow)
  --norrna      Skip rRNA search
  --notrna      Skip tRNA search
```

## Citation

If you use prokka-rs, please cite the original Prokka paper:

> Seemann T. "Prokka: rapid prokaryotic genome annotation."
> Bioinformatics, 2014 Jul 15;30(14):2068-9.
> [DOI:10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153)

## License

GPL-3.0
