//! CDS annotation (pipeline step 7).
//!
//! Annotates predicted CDS features by similarity search against reference
//! databases. Replaces Perl Prokka's external `blastp`/`hmmscan` calls with
//! native Rust libraries (`blast-rs`, `hmmer-pure-rs`).
//!
//! Annotation databases use the Prokka `~~~`-delimited FASTA description
//! format: `>SeqID EC_number~~~gene~~~product~~~COG`. Parsing of these
//! deflines is shared between BLAST and HMMER searches and lives in
//! [`database`].

pub mod blast;
pub mod hmmer;
pub mod database;
