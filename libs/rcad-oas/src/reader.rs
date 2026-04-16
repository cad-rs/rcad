//! OASIS file reader.

use crate::{OasError, OasLibrary};
use std::path::Path;

/// OASIS file reader.
pub struct OasReader {
    // Placeholder for future implementation
}

impl OasReader {
    /// Create a new OASIS reader.
    pub fn new() -> Self {
        Self {}
    }

    /// Read an OASIS file.
    pub fn read_file<P: AsRef<Path>>(&self, _path: P) -> Result<OasLibrary, OasError> {
        // TODO: Implement using laykit
        Err(OasError::InvalidFormat("OASIS reading not yet implemented".to_string()))
    }
}

impl Default for OasReader {
    fn default() -> Self {
        Self::new()
    }
}
