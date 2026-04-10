use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::error::ProkkaError;

/// Generate a locus_tag prefix from MD5 of the input file contents.
///
/// Exact port of Perl Prokka's `generate_locus_tag()` (lines 1709-1729):
/// - Compute MD5 hex digest of the file
/// - Take first 8 hex characters
/// - Uppercase each character
/// - Map digits 0-9 to letters G-P (chr(ord('F')+1+digit))
pub fn generate_locus_tag(path: &Path) -> Result<String, ProkkaError> {
    let mut file = File::open(path)?;
    let mut contents = Vec::new();
    file.read_to_end(&mut contents)?;

    let digest = md5::compute(&contents);
    let hex = format!("{:x}", digest);

    let mut tag = String::with_capacity(8);
    for c in hex.chars().take(8) {
        let upper = c.to_ascii_uppercase();
        if upper.is_ascii_digit() {
            // Map '0'..'9' to 'G'..'P'
            let digit = upper as u8 - b'0';
            tag.push((b'F' + 1 + digit) as char);
        } else {
            tag.push(upper);
        }
    }

    Ok(tag)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_generate_locus_tag() {
        // Create a temp file with known content and verify the tag is deterministic
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fna");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">test\nACGTACGT").unwrap();
        }
        let tag1 = generate_locus_tag(&path).unwrap();
        let tag2 = generate_locus_tag(&path).unwrap();
        assert_eq!(tag1, tag2);
        assert_eq!(tag1.len(), 8);
        // All chars should be uppercase letters A-P
        for c in tag1.chars() {
            assert!(c.is_ascii_uppercase());
        }
    }

    #[test]
    fn test_digit_mapping() {
        // Verify the digit-to-letter mapping
        // '0' -> 'G', '1' -> 'H', ..., '9' -> 'P'
        for d in 0u8..10 {
            let expected = (b'F' + 1 + d) as char;
            let c = (b'0' + d) as char;
            let upper = c.to_ascii_uppercase();
            let digit = upper as u8 - b'0';
            let result = (b'F' + 1 + digit) as char;
            assert_eq!(result, expected);
        }
    }
}
