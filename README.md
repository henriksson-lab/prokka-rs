# prokka-rs

Rapid prokaryotic genome annotation — [Prokka](https://github.com/tseemann/prokka) rewritten in Rust.

prokka-rs annotates bacterial, archaeal, viral, and mitochondrial genomes by identifying features
(tRNA, rRNA, ncRNA, CRISPRs, CDS) and assigning functional labels via hierarchical database searches.


## This is an LLM-mediated faithful (hopefully) translation, not the original code!

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this page](https://rewrites.bio/)
but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back to the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information



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

Compared against Prokka 1.15.6 (Perl) on Linux x86_64.
Rust binary compiled with `RUSTFLAGS="-C target-cpu=native" cargo build --release`.

### CDS prediction only (`--noanno`, 1 CPU)

| Input | Perl Prokka | prokka-rs | Speedup | CDS found |
|-------|------------|-----------|---------|-----------|
| plasmid.fna (57 KB) | 3.6s | 0.35s | **10x** | 63 (both) |
| genome.fna (6.7 MB) | 41.3s | 15.9s | **2.6x** | 6048 vs 6059 |

prokka-rs uses prodigal-rs natively (no subprocess, no BioPerl overhead).

### With BLAST annotation (`--fast`)

| Input | CPUs | Perl Prokka | prokka-rs | Speedup | Annotations |
|-------|------|------------|-----------|---------|-------------|
| plasmid.fna (57 KB) | 1 | 12.8s | 13.8s | 0.9x | 10 vs 11 |
| plasmid.fna (57 KB) | 8 | 7.5s | **2.7s** | **2.8x** | 10 vs 11 |

blast-rs uses the NCBI BLAST+ two-hit algorithm with thick backbone lookup
tables and rolling hash. Per-query speed is ~1.7x NCBI on small queries.
With 8 CPUs, queries are parallelized across threads via rayon.

prokka-rs finds 1 extra annotation because blast-rs doesn't yet implement
NCBI's composition-based statistics adjustment (comp_based_stats mode 2).
This causes one false positive that NCBI would filter. See TODO.md.

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
