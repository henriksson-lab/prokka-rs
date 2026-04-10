use std::path::Path;

/// Inventory of available databases.
#[derive(Debug, Default)]
pub struct DbInventory {
    pub kingdoms: Vec<String>,
    pub genera: Vec<String>,
    pub hmms: Vec<String>,
    pub cms: Vec<String>,
}

/// List all configured databases.
///
/// Replicates Perl Prokka's `list_db()`, `kingdoms()`, `genera()`, `hmms()`, `cms()`.
pub fn list_db(dbdir: &Path) -> DbInventory {
    let mut inv = DbInventory::default();

    // Kingdoms: look for *.pin files in kingdom/*/
    let kingdom_dir = dbdir.join("kingdom");
    if kingdom_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&kingdom_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    // Check if sprot.pin exists (indexed)
                    if path.join("sprot.pin").exists() || path.join("sprot").exists() {
                        if let Some(name) = path.file_name() {
                            inv.kingdoms.push(name.to_string_lossy().to_string());
                        }
                    }
                }
            }
        }
    }
    inv.kingdoms.sort();

    // Genera: look for *.pin files in genus/
    let genus_dir = dbdir.join("genus");
    if genus_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&genus_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if let Some(ext) = path.extension() {
                    if ext == "pin" {
                        if let Some(stem) = path.file_stem() {
                            inv.genera.push(stem.to_string_lossy().to_string());
                        }
                    }
                }
                // Also detect unindexed genus files (no .pin but file exists)
                if path.extension().is_none() && path.is_file() {
                    if let Some(name) = path.file_name() {
                        let name = name.to_string_lossy().to_string();
                        if !inv.genera.contains(&name) {
                            inv.genera.push(name);
                        }
                    }
                }
            }
        }
    }
    inv.genera.sort();

    // HMMs: look for *.h3m files in hmm/
    let hmm_dir = dbdir.join("hmm");
    if hmm_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&hmm_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let name = path.file_name().unwrap_or_default().to_string_lossy().to_string();
                if name.ends_with(".h3m") {
                    inv.hmms.push(name.trim_end_matches(".h3m").trim_end_matches(".hmm").to_string());
                } else if name.ends_with(".hmm") || name.ends_with(".hmm.gz") {
                    let base = name.trim_end_matches(".gz").trim_end_matches(".hmm").to_string();
                    if !inv.hmms.contains(&base) {
                        inv.hmms.push(base);
                    }
                }
            }
        }
    }
    inv.hmms.sort();

    // CMs: look for *.i1m files in cm/
    let cm_dir = dbdir.join("cm");
    if cm_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&cm_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let name = path.file_name().unwrap_or_default().to_string_lossy().to_string();
                if name.ends_with(".i1m") {
                    inv.cms.push(name.trim_end_matches(".i1m").to_string());
                } else if path.is_file() && path.extension().is_none() {
                    let n = name.clone();
                    if !inv.cms.contains(&n) && n != "README" && !n.starts_with("__") {
                        inv.cms.push(n);
                    }
                }
            }
        }
    }
    inv.cms.sort();

    inv
}
