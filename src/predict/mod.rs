//! Feature prediction modules.
//!
//! Implements pipeline steps 2-6 from the Prokka pipeline (see CLAUDE.md):
//! tRNA/tmRNA (Aragorn), rRNA (Barrnap/RNAmmer), ncRNA (Infernal/cmscan),
//! CRISPRs (minced), and CDS (Prodigal via the native `prodigal-rs` crate).
//! Also provides the post-prediction CDS/RNA overlap filter and the
//! SignalP signal-peptide detection step.

pub mod cds;
pub mod trna;
pub mod rrna;
pub mod ncrna;
pub mod crispr;
pub mod signalp;
pub mod overlap;
