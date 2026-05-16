//! Implementation of `--cleandb`: remove generated database index sidecars.

use std::path::Path;

use crate::error::ProkkaError;

/// Remove all database index files.
///
/// Replicates Perl Prokka's `clean_db()` (lines 1685-1694):
/// - Remove kingdom/*/*.p?? (BLAST indices)
/// - Remove genus/*.p?? (BLAST indices)
/// - Remove hmm/*.h3? (HMMER indices)
/// - Remove cm/*.i1? (Infernal indices)
pub fn clean_db(dbdir: &Path) -> Result<(), ProkkaError> {
    remove_glob(dbdir, "kingdom/*/", &["pin", "psq", "phr", "pdb", "pot", "ptf", "pto"]);
    remove_glob(dbdir, "genus/", &["pin", "psq", "phr", "pdb", "pot", "ptf", "pto"]);
    remove_glob(dbdir, "hmm/", &["h3i", "h3f", "h3m", "h3p"]);
    remove_glob(dbdir, "cm/", &["i1i", "i1f", "i1m", "i1p"]);
    Ok(())
}

/// Remove every file under `dbdir/subdir` whose extension appears in
/// `extensions`.
///
/// If `subdir` contains a `*`, the portion before the first `*` is treated
/// as a parent directory and every immediate subdirectory of it is swept
/// (so `"kingdom/*/"` walks each per-kingdom directory).
fn remove_glob(dbdir: &Path, subdir: &str, extensions: &[&str]) {
    // Handle wildcard in subdir (e.g. "kingdom/*/")
    if subdir.contains('*') {
        let parts: Vec<&str> = subdir.split('*').collect();
        let parent = dbdir.join(parts[0].trim_end_matches('/'));
        if let Ok(entries) = std::fs::read_dir(&parent) {
            for entry in entries.flatten() {
                if entry.path().is_dir() {
                    remove_index_files(&entry.path(), extensions);
                }
            }
        }
    } else {
        let dir = dbdir.join(subdir.trim_end_matches('/'));
        remove_index_files(&dir, extensions);
    }
}

/// Delete every regular file in `dir` whose extension matches one of
/// `extensions`. Missing directories and per-file removal errors are
/// silently ignored, matching Perl Prokka's `delfile(<glob>)` semantics.
fn remove_index_files(dir: &Path, extensions: &[&str]) {
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(ext) = path.extension() {
                let ext_str = ext.to_string_lossy();
                if extensions.iter().any(|e| ext_str == *e) {
                    let _ = std::fs::remove_file(&path);
                }
            }
        }
    }
}
