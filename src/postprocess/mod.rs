//! Post-processing of predicted features before output (pipeline steps 8-9).
//!
//! The submodules here run in this order over the collected features:
//!
//! 1. [`product`] — cleanup of `/product` annotation strings so they are
//!    Genbank/ENA compliant (drops DUF/UPF references, normalises
//!    "putative"/"predicted", etc.).
//! 2. [`gene_dedup`] — renames colliding `/gene` names by appending `_1`,
//!    `_2`, ... so each CDS gene tag is unique across all contigs.
//! 3. [`locus_assign`] — sorts features by start position and assigns
//!    `locus_tag`, `ID`, and (for `--centre`) `protein_id` tags. Optionally
//!    spawns sister `gene` / `mRNA` features for `--addgenes` / `--addmrna`.
//!
//! The default locus-tag prefix is derived deterministically from the MD5
//! hash of the input FASTA file (see [`locus_assign`] / `crate::locus_tag`).

pub mod product;
pub mod gene_dedup;
pub mod locus_assign;
