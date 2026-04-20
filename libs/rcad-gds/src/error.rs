use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum GdsError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),

    #[error("Invalid GDS format: {0}")]
    InvalidFormat(String),

    #[error("Cell not found: {0}")]
    CellNotFound(String),

    #[error("Layer {0} not configured")]
    LayerNotConfigured(i32),

    #[error("Geometry conversion failed: {0}")]
    GeometryError(String),

    #[error("Empty structure: {0}")]
    EmptyStructure(String),

    #[error("laykit parsing error: {0}")]
    Laykit(String),
}

pub type Result<T> = std::result::Result<T, GdsError>;
