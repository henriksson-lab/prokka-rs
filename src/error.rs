//! Unified error type for the crate.
//!
//! Every fallible function in prokka-rs returns `Result<T, ProkkaError>`.
//! The Perl pipeline uses [`sub err`](../../prokka/bin/prokka) (line 1557),
//! which simply prints a message and exits — here the variants carry enough
//! context for callers to act on or render a user-friendly message.

use std::path::PathBuf;

/// All errors that can occur in prokka-rs.
#[derive(Debug, thiserror::Error)]
pub enum ProkkaError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FASTA file '{path}' contains no suitable sequence entries")]
    NoContigs { path: PathBuf },

    #[error("FASTA file '{path}' contains duplicate sequence ID: {id}")]
    DuplicateContigId { path: PathBuf, id: String },

    #[error("contig ID too long ({len} > {max}): {id}")]
    ContigIdTooLong { id: String, len: usize, max: usize },

    #[error("invalid genetic code: {0} (must be 1..25)")]
    InvalidGeneticCode(u8),

    #[error("invalid kingdom: {0}")]
    InvalidKingdom(String),

    #[error("invalid e-value: {0}")]
    InvalidEvalue(f64),

    #[error("invalid coverage: {0} (must be 0..100)")]
    InvalidCoverage(f64),

    #[error("output folder '{0}' already exists (use --force to overwrite)")]
    OutputDirExists(PathBuf),

    #[error("database not indexed — run 'prokka-rs --setupdb' first")]
    DatabaseNotIndexed,

    #[error("external tool '{tool}' not found in PATH")]
    ToolNotFound { tool: String },

    #[error("external tool '{tool}' failed: {message}")]
    ToolFailed { tool: String, message: String },

    #[error("gene prediction error: {0}")]
    Prodigal(String),

    #[error("BLAST search error: {0}")]
    Blast(String),

    #[error("HMMER search error: {0}")]
    Hmmer(String),

    #[error("cannot read file: {0}")]
    FileNotReadable(PathBuf),

    #[error("{0}")]
    Other(String),
}
