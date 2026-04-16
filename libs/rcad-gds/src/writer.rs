//! GDS writer module.

use crate::{GdsLibrary, GdsError};

/// GDS file writer.
pub struct GdsWriter;

impl GdsWriter {
    /// Create a new GDS writer.
    pub fn new() -> Self {
        Self
    }

    /// Write a GDS library to file.
    pub fn write(&self, _library: &GdsLibrary, _path: &std::path::Path) -> Result<(), GdsError> {
        // TODO: Implement GDS writing
        Err(GdsError::Io(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "GDS writing not yet implemented",
        )))
    }
}

impl Default for GdsWriter {
    fn default() -> Self {
        Self::new()
    }
}
