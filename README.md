# prokka-rs

Rapid prokaryotic genome annotation — [Prokka](https://github.com/tseemann/prokka) rewritten in Rust.

prokka-rs annotates bacterial, archaeal, viral, and mitochondrial genomes by identifying features
(tRNA, rRNA, ncRNA, CRISPRs, CDS) and assigning functional labels via hierarchical database searches.

This is a translation of the original code and not the authorative implementation. This code should generate bitwise
equal output to the original. Please report any deviations

The aim of this project is to increase performance, especially by providing this code through a type-safe library interface.
The code can also be compiled to be used for webassembly.


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

Compared against Prokka 1.15.6 (Perl) on Linux x86_64 (single-threaded, `--cpus 1`).
Rust binary compiled with `RUSTFLAGS="-C target-cpu=native" cargo build --release`.

### CDS prediction only (`--noanno`)

| Input | Perl Prokka | prokka-rs | Speedup | CDS found |
|-------|------------|-----------|---------|-----------|
| plasmid.fna (57 KB) | 4.6s | 0.7s | **6.6x** | 63 (both) |
| genome.fna (6.7 MB) | 58.7s | 39.7s | **1.5x** | identical |

Both tools find the same number of CDS. prokka-rs uses prodigal-rs natively
(no subprocess overhead), which accounts for most of the speedup on small inputs.

### With BLAST annotation (`--fast`, 1 CPU)

| Input | Perl Prokka | prokka-rs | Annotations found |
|-------|------------|-----------|-------------------|
| plasmid.fna (57 KB) | 16.1s | 2m 6s | Perl: 10, Rust: 7 |
| genome.fna (6.7 MB) | 20m 8s | TBD | |

blast-rs uses an indexed protein database (built at runtime) for searching.
Still ~8x slower than NCBI BLAST+ per query, but functional. The 3 missing
annotations are due to prodigal-rs gene boundary differences producing
slightly different protein sequences.

### With BLAST annotation (`--fast`, 8 CPUs)

| Input | Perl Prokka |
|-------|------------|
| plasmid.fna (57 KB) | 8.8s |

prokka-rs does not yet parallelize BLAST searches across CPUs.
This is planned using rayon.

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
