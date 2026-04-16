//! Traits for unified layout reading and writing.

use std::path::Path;

use crate::error::Result;
use crate::format::LayoutFormat;

/// Trait for reading layout files in a format-agnostic way.
pub trait LayoutReader {
    /// The library type produced by this reader.
    type Library;

    /// Reads a layout file from the given path.
    fn read_file<P: AsRef<Path>>(&self, path: P) -> Result<Self::Library>;

    /// Reads a layout from bytes.
    fn read_bytes(&self, bytes: &[u8]) -> Result<Self::Library>;
}

/// Trait for writing layout files in a format-agnostic way.
pub trait LayoutWriter {
    /// The library type to write.
    type Library;

    /// Writes a layout library to a file.
    fn write_file<P: AsRef<Path>>(&self, library: &Self::Library, path: P) -> Result<()>;

    /// Writes a layout library to a byte vector.
    fn write_bytes(&self, library: &Self::Library) -> Result<Vec<u8>>;
}

/// Factory for creating format-specific readers and writers.
pub struct LayoutIoFactory;

impl LayoutIoFactory {
    /// Creates a reader for the given format.
    pub fn reader(format: LayoutFormat) -> Box<dyn std::any::Any> {
        match format {
            LayoutFormat::Gds => Box::new(rcad_gds::GdsReader),
            LayoutFormat::Oasis => Box::new(rcad_oas::OasReader),
        }
    }

    /// Creates a writer for the given format.
    pub fn writer(format: LayoutFormat) -> Box<dyn std::any::Any> {
        match format {
            LayoutFormat::Gds => Box::new(rcad_gds::GdsWriter::new()),
            LayoutFormat::Oasis => Box::new(rcad_oas::OasWriter::new()),
        }
    }
}

/// Helper functions for reading layout files.
pub mod read {
    use super::*;
    use crate::IoError;

    /// Reads a GDS file and returns a GdsLibrary.
    pub fn gds<P: AsRef<Path>>(path: P) -> Result<rcad_gds::GdsLibrary> {
        rcad_gds::GdsReader::read_file(path.as_ref()).map_err(IoError::from)
    }

    /// Reads an OASIS file and returns an OasLibrary.
    pub fn oasis<P: AsRef<Path>>(path: P) -> Result<rcad_oas::OasLibrary> {
        rcad_oas::OasReader::read_file(path.as_ref()).map_err(IoError::from)
    }
}

/// Helper functions for writing layout files.
pub mod write {
    use super::*;
    use crate::IoError;

    /// Writes a GdsLibrary to a GDS file.
    pub fn gds<P: AsRef<Path>>(library: &rcad_gds::GdsLibrary, path: P) -> Result<()> {
        let writer = rcad_gds::GdsWriter::new();
        writer.write_file(library, path.as_ref()).map_err(IoError::from)
    }

    /// Writes an OasLibrary to an OASIS file.
    pub fn oasis<P: AsRef<Path>>(library: &rcad_oas::OasLibrary, path: P) -> Result<()> {
        let writer = rcad_oas::OasWriter::new();
        writer.write_file(library, path.as_ref()).map_err(IoError::from)
    }
}
