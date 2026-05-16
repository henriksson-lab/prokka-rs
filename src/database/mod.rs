//! Annotation database management for the `--setupdb`, `--cleandb` and
//! `--listdb` CLI subcommands.
//!
//! These mirror the equivalent Perl Prokka subroutines (`setup_db`,
//! `clean_db`, `list_db` around lines 1632-1694 of `prokka/bin/prokka`) and
//! operate on the standard `$dbdir/{kingdom,genus,hmm,cm}` layout.

pub mod setup;
pub mod clean;
pub mod list;
