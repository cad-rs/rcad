//! Error types for rcad-io.

use std::io;
use thiserror::Error;

/// I/O error type for layout file operations.
#[derive(Debug, Error)]
pub enum IoError {
    /// GDS format error.
    #[error("GDS error: {0}")]
    Gds(#[from] rcad_gds::GdsError),

    /// OASIS format error.
    #[error("OASIS error: {0}")]
    Oasis(#[from] rcad_oas::OasError),

    /// I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    /// Format detection error.
    #[error("Unknown format: {0}")]
    UnknownFormat(String),

    /// Unsupported format.
    #[error("Unsupported format: {0}")]
    UnsupportedFormat(String),

    /// Invalid layer configuration.
    #[error("Invalid layer configuration: {0}")]
    InvalidLayerConfig(String),

    /// Serialization error.
    #[error("Serialization error: {0}")]
    Serialization(#[from] serde_json::Error),

    /// File not found.
    #[error("File not found: {0}")]
    FileNotFound(String),

    /// Invalid data.
    #[error("Invalid data: {0}")]
    InvalidData(String),
}

/// Result type for I/O operations.
pub type Result<T> = std::result::Result<T, IoError>;
