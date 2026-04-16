//! GDS reader module.

use crate::{GdsLibrary, GdsError};

/// GDS file reader.
pub struct GdsReader;

impl GdsReader {
    /// Create a new GDS reader.
    pub fn new() -> Self {
        Self
    }

    /// Read a GDS file.
    pub fn read(&self, _path: &std::path::Path) -> Result<GdsLibrary, GdsError> {
        // TODO: Implement GDS reading
        Err(GdsError::Io(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "GDS reading not yet implemented",
        )))
    }
}

impl Default for GdsReader {
    fn default() -> Self {
        Self::new()
    }
}
