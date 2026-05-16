//! # prokka-rs
//!
//! A Rust rewrite of [Prokka](https://github.com/tseemann/prokka), the rapid
//! prokaryotic genome annotation pipeline by Torsten Seemann. The goal is
//! byte-identical output with the original Perl `prokka` script, while using
//! native Rust libraries (`prodigal-rs`, `blast-rs`, `hmmer-pure-rs`) in place
//! of external subprocess calls where possible.
//!
//! ## Library entry points
//!
//! Most library consumers only need a handful of items, which are re-exported
//! from the crate root:
//!
//! - [`annotate()`] / [`run()`] — pipeline orchestrators
//! - [`ProkkaConfig`] — configuration (every CLI flag has a field)
//! - [`Contig`], [`SeqFeature`], [`AnnotationResult`] — core data model
//! - [`ProkkaError`] — unified error type
//!
//! The CLI binary (`src/main.rs`) is a thin wrapper around [`run`], gated
//! behind the on-by-default `cli` Cargo feature.

pub mod codon_table;
pub mod config;
pub mod genbank_to_fasta;
pub mod error;
pub mod input;
pub mod locus_tag;
pub mod model;
pub mod pipeline;

pub mod predict;
pub mod annotate;
pub mod postprocess;
pub mod output;
pub mod database;

// Re-export key types for library users.
pub use config::ProkkaConfig;
pub use error::ProkkaError;
pub use model::{AnnotationResult, AnnotationStats, Contig, SeqFeature};
pub use pipeline::{annotate, run};
