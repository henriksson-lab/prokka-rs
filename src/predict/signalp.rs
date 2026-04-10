use std::path::Path;

use crate::config::ProkkaConfig;
use crate::error::ProkkaError;
use crate::model::{Contig, SeqFeature};

/// Detect signal peptides using SignalP (external process).
///
/// SignalP is proprietary and version-dependent (v3/v4/v5 have different interfaces).
/// This is only enabled when --gram is specified and signalp is in PATH.
///
/// Replicates Perl Prokka lines 801-938.
pub fn predict_signalp(
    outdir: &Path,
    contigs: &[Contig],
    config: &ProkkaConfig,
) -> Result<Vec<SeqFeature>, ProkkaError> {
    // SignalP is optional and requires --gram
    if config.gram.is_none() {
        return Ok(Vec::new());
    }

    // Check if signalp is available
    if std::process::Command::new("signalp")
        .arg("-version")
        .output()
        .is_err()
    {
        if !config.quiet {
            eprintln!("SignalP not found in PATH, skipping signal peptide detection.");
        }
        return Ok(Vec::new());
    }

    // TODO: Implement full SignalP v3/v4/v5 parsing
    // This is complex due to version-specific output formats and is deferred
    // to a later phase. For now, signal peptide detection is skipped.
    let _ = (outdir, contigs);
    Ok(Vec::new())
}
