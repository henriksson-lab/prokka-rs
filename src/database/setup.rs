//! Implementation of `--setupdb`: build BLAST/HMM/CM indices for every
//! installed annotation database under `$dbdir`.
//!
//! This shells out to the upstream NCBI/HMMER/Infernal tools
//! (`makeblastdb`, `hmmpress`, `cmpress`) — the same external commands the
//! Perl `setup_db()` routine drives near line 1643 of `prokka/bin/prokka`.

use std::path::Path;
use std::process::Command;

use crate::error::ProkkaError;

/// Index all installed databases (BLAST, HMM, CM).
///
/// Replicates Perl Prokka's `setup_db()` (lines 1643-1681).
/// Currently uses external tools (makeblastdb, hmmpress, cmpress).
pub fn setup_db(dbdir: &Path) -> Result<(), ProkkaError> {
    // First clean existing indices
    crate::database::clean::clean_db(dbdir)?;

    // Index kingdom BLAST databases
    for db_name in &["sprot", "IS", "AMR"] {
        let kingdom_dir = dbdir.join("kingdom");
        if let Ok(entries) = std::fs::read_dir(&kingdom_dir) {
            for entry in entries.flatten() {
                let fasta = entry.path().join(db_name);
                if fasta.exists() && fasta.is_file() {
                    eprintln!("Making kingdom BLASTP database: {}", fasta.display());
                    run_makeblastdb(&fasta)?;
                }
            }
        }
    }

    // Index genus databases
    let genus_dir = dbdir.join("genus");
    if genus_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&genus_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() && path.extension().is_none() {
                    eprintln!("Making genus BLASTP database: {}", path.display());
                    run_makeblastdb(&path)?;
                }
            }
        }
    }

    // Press HMM databases
    let hmm_dir = dbdir.join("hmm");
    if hmm_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&hmm_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let name = path.file_name().unwrap_or_default().to_string_lossy().to_string();

                // Decompress .gz first
                if name.ends_with(".hmm.gz") {
                    let uncompressed = path.with_extension(""); // remove .gz
                    eprintln!("Uncompressing: {}", path.display());
                    let output = Command::new("gzip")
                        .arg("-dc")
                        .arg(&path)
                        .output();
                    if let Ok(out) = output {
                        if out.status.success() {
                            let _ = std::fs::write(&uncompressed, &out.stdout);
                            run_hmmpress(&uncompressed)?;
                        }
                    }
                } else if name.ends_with(".hmm") {
                    eprintln!("Pressing HMM database: {}", path.display());
                    run_hmmpress(&path)?;
                }
            }
        }
    }

    // Press CM databases
    let cm_dir = dbdir.join("cm");
    if cm_dir.exists() {
        for kingdom in &["Viruses", "Bacteria", "Archaea"] {
            let cm_file = cm_dir.join(kingdom);
            if cm_file.exists() {
                eprintln!("Pressing CM database: {}", cm_file.display());
                run_cmpress(&cm_file)?;
            }
        }
    }

    Ok(())
}

/// Index a protein FASTA with `makeblastdb -hash_index -dbtype prot`.
///
/// Errors map missing executables to [`ProkkaError::ToolNotFound`] and
/// non-zero exits to [`ProkkaError::ToolFailed`]. Matches the command line
/// at line 1653 of `prokka/bin/prokka`.
fn run_makeblastdb(fasta: &Path) -> Result<(), ProkkaError> {
    let status = Command::new("makeblastdb")
        .arg("-hash_index")
        .arg("-dbtype")
        .arg("prot")
        .arg("-in")
        .arg(fasta)
        .arg("-logfile")
        .arg("/dev/null")
        .status();

    match status {
        Ok(s) if s.success() => Ok(()),
        Ok(_) => Err(ProkkaError::ToolFailed {
            tool: "makeblastdb".to_string(),
            message: format!("Failed to index {}", fasta.display()),
        }),
        Err(_) => Err(ProkkaError::ToolNotFound {
            tool: "makeblastdb".to_string(),
        }),
    }
}

/// Run `hmmpress` on an HMM database to build the `.h3{i,f,m,p}` sidecars.
///
/// Mirrors the `hmmpress \Q$hmm\E` call at line 1671 of `prokka/bin/prokka`.
fn run_hmmpress(hmm: &Path) -> Result<(), ProkkaError> {
    let status = Command::new("hmmpress")
        .arg(hmm)
        .status();

    match status {
        Ok(s) if s.success() => Ok(()),
        Ok(_) => Err(ProkkaError::ToolFailed {
            tool: "hmmpress".to_string(),
            message: format!("Failed to press {}", hmm.display()),
        }),
        Err(_) => Err(ProkkaError::ToolNotFound {
            tool: "hmmpress".to_string(),
        }),
    }
}

/// Run `cmpress` on an Infernal covariance model to build the
/// `.i1{i,f,m,p}` sidecars.
///
/// Mirrors the `cmpress \Q$cm\E` call at line 1677 of `prokka/bin/prokka`.
fn run_cmpress(cm: &Path) -> Result<(), ProkkaError> {
    let status = Command::new("cmpress")
        .arg(cm)
        .status();

    match status {
        Ok(s) if s.success() => Ok(()),
        Ok(_) => Err(ProkkaError::ToolFailed {
            tool: "cmpress".to_string(),
            message: format!("Failed to press {}", cm.display()),
        }),
        Err(_) => Err(ProkkaError::ToolNotFound {
            tool: "cmpress".to_string(),
        }),
    }
}
