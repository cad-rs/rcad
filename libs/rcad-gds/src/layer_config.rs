//! Layer configuration for GDS import.

use serde::{Deserialize, Serialize};

/// Layer configuration for extrusion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerConfig {
    pub layer: i16,
    pub settings: LayerSettings,
}

/// Settings for layer extrusion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerSettings {
    pub height: f64,
    pub z_offset: f64,
}

impl Default for LayerSettings {
    fn default() -> Self {
        Self {
            height: 1.0,
            z_offset: 0.0,
        }
    }
}
