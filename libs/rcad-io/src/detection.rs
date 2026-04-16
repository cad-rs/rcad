//! Format detection for layout files.

use std::path::Path;
use std::io::Read;

use crate::format::LayoutFormat;
use crate::error::IoError;

/// Magic bytes for GDSII files.
/// GDSII starts with a 2-byte record header followed by BGNSTR (0x0002)
/// or more commonly, the first record is HEADER (0x0000).
/// GDSII uses big-endian 2-byte length + 2-byte record type.
#[allow(dead_code)]
const GDS_MAGIC: &[u8] = &[0x00, 0x06]; // Typical first bytes: length=6 for HEADER record

/// Magic bytes for OASIS files.
/// OASIS files start with "%SEMI-OASIS\r\n" (13 bytes)
const OASIS_MAGIC: &[u8] = b"%SEMI-OASIS";

/// Detects the format of a layout file from its path.
///
/// # Arguments
/// * `path` - Path to the file
///
/// # Returns
/// * `Some(LayoutFormat)` if the format could be detected
/// * `None` if the format is unknown or cannot be read
pub fn detect_format<P: AsRef<Path>>(path: P) -> Option<LayoutFormat> {
    let path = path.as_ref();

    // Try extension first
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        if let Some(format) = LayoutFormat::from_extension(ext) {
            return Some(format);
        }
    }

    // Fall back to magic byte detection
    detect_format_from_content(path).ok().flatten()
}

/// Detects format from file content using magic bytes.
fn detect_format_from_content(path: &Path) -> Result<Option<LayoutFormat>, IoError> {
    let mut file = std::fs::File::open(path)?;
    let mut buffer = [0u8; 32];
    let bytes_read = file.read(&mut buffer)?;

    if bytes_read < 4 {
        return Ok(None);
    }

    Ok(detect_format_from_bytes(&buffer[..bytes_read]))
}

/// Detects format from the first few bytes of a file.
///
/// # Arguments
/// * `bytes` - First few bytes of the file (at least 4 bytes recommended)
///
/// # Returns
/// * `Some(LayoutFormat)` if the format could be detected
/// * `None` if the format is unknown
pub fn detect_format_from_bytes(bytes: &[u8]) -> Option<LayoutFormat> {
    if bytes.len() < 4 {
        return None;
    }

    // Check for OASIS magic bytes first (more specific)
    if bytes.starts_with(OASIS_MAGIC) {
        return Some(LayoutFormat::Oasis);
    }

    // Check for GDSII format
    // GDSII records: 2-byte length (big-endian) + 2-byte record type
    // First record is typically HEADER (0x0002) with length 6
    // So first bytes are typically: 0x00 0x06 0x00 0x02
    if is_gds_format(bytes) {
        return Some(LayoutFormat::Gds);
    }

    None
}

/// Checks if bytes represent a GDSII file.
fn is_gds_format(bytes: &[u8]) -> bool {
    if bytes.len() < 4 {
        return false;
    }

    // GDSII record structure: 2-byte length (big-endian) + 2-byte record type
    let len = u16::from_be_bytes([bytes[0], bytes[1]]) as usize;
    let record_type = u16::from_be_bytes([bytes[2], bytes[3]]);

    // Valid GDSII record lengths are 4-65532 bytes (including the 4-byte header)
    // First record should be HEADER (0x0002) or BGNLIB (0x0100)
    // Length should be reasonable (typically 6 for HEADER)
    let is_valid_length = len >= 4 && len <= 65532;
    let is_valid_first_record = matches!(record_type, 0x0002 | 0x0100);

    is_valid_length && is_valid_first_record
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_gds_from_extension() {
        assert_eq!(detect_format("test.gds"), Some(LayoutFormat::Gds));
        assert_eq!(detect_format("test.gds2"), Some(LayoutFormat::Gds));
        assert_eq!(detect_format("test.GDS"), Some(LayoutFormat::Gds));
    }

    #[test]
    fn test_detect_oasis_from_extension() {
        assert_eq!(detect_format("test.oas"), Some(LayoutFormat::Oasis));
        assert_eq!(detect_format("test.oasis"), Some(LayoutFormat::Oasis));
        assert_eq!(detect_format("test.OAS"), Some(LayoutFormat::Oasis));
    }

    #[test]
    fn test_detect_from_bytes_gds() {
        // Typical GDSII header record: length=6, type=HEADER(0x0002)
        let gds_bytes = [0x00, 0x06, 0x00, 0x02, 0x00, 0x01, 0x00, 0x00];
        assert_eq!(detect_format_from_bytes(&gds_bytes), Some(LayoutFormat::Gds));
    }

    #[test]
    fn test_detect_from_bytes_oasis() {
        let oasis_bytes = b"%SEMI-OASIS\r\n";
        assert_eq!(detect_format_from_bytes(oasis_bytes), Some(LayoutFormat::Oasis));
    }

    #[test]
    fn test_detect_from_bytes_unknown() {
        let unknown_bytes = [0x00, 0x00, 0x00, 0x00];
        assert_eq!(detect_format_from_bytes(&unknown_bytes), None);
    }

    #[test]
    fn test_format_properties() {
        assert_eq!(LayoutFormat::Gds.extension(), "gds");
        assert_eq!(LayoutFormat::Oasis.extension(), "oas");
        assert_eq!(LayoutFormat::Gds.name(), "GDSII");
        assert_eq!(LayoutFormat::Oasis.name(), "OASIS");
    }
}
