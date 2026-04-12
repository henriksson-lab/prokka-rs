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
