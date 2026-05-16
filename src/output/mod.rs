//! Output file writers for the annotation pipeline (pipeline step 10).
//!
//! Each submodule writes one of the Prokka output formats:
//!
//! - [`gff`] — GFF3 annotation table with embedded FASTA section
//! - [`tbl`] — NCBI feature table (used as input to `tbl2asn`)
//! - [`tsv`] — tab-separated annotation summary (one row per feature)
//! - [`fasta`] — `.fna` (contigs), `.faa` (proteins), `.ffn` (feature nucleotide), `.fsa` (annotated FASTA)
//! - [`txt`] — annotation statistics
//! - [`genbank`] — `.gbk` / `.sqn` via `tbl2asn`, with a basic in-house fallback
//!
//! All writers consume an [`AnnotationResult`](crate::model::AnnotationResult)
//! produced by the upstream [`pipeline`](crate::pipeline) and reproduce the
//! exact output formatting of the Perl Prokka reference implementation
//! (`prokka/bin/prokka`, lines ~1305-1437).

pub mod gff;
pub mod tbl;
pub mod tsv;
pub mod fasta;
pub mod txt;
pub mod genbank;
