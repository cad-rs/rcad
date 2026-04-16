//! Conversion between OASIS geometry and RCAD kernel types.

use crate::{OasError, OasCell, OasPolygon, OasPath};
use crate::layer_config::LayerConfig;

/// Convert OASIS geometry to RCAD kernel shapes.
pub struct OasConverter {
    config: LayerConfig,
}

impl OasConverter {
    /// Create a new converter with layer configuration.
    pub fn new(config: LayerConfig) -> Self {
        Self { config }
    }

    /// Convert an OASIS polygon to an RCAD face.
    pub fn polygon_to_face(&self, _polygon: &OasPolygon) -> Result<(), OasError> {
        // TODO: Implement using rcad-kernel
        Err(OasError::GeometryError("Polygon conversion not yet implemented".to_string()))
    }

    /// Convert an OASIS path to an RCAD face.
    pub fn path_to_face(&self, _path: &OasPath) -> Result<(), OasError> {
        // TODO: Implement using rcad-kernel
        Err(OasError::GeometryError("Path conversion not yet implemented".to_string()))
    }

    /// Convert an OASIS cell to an RCAD compound.
    pub fn cell_to_compound(&self, _cell: &OasCell) -> Result<(), OasError> {
        // TODO: Implement using rcad-kernel
        Err(OasError::GeometryError("Cell conversion not yet implemented".to_string()))
    }
}

impl Default for OasConverter {
    fn default() -> Self {
        Self::new(LayerConfig::default())
    }
}
