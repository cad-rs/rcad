//! Unified layer configuration for layout formats.

use serde::{Deserialize, Serialize};

/// Unified layer configuration that works across GDS and OASIS formats.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnifiedLayerConfig {
    /// Layer number (0-65535 for GDS, 0-65535 for OASIS)
    pub layer: u16,
    /// Data type (GDS) / purpose (OASIS)
    pub datatype: u16,
    /// Optional layer name (OASIS only, ignored in GDS)
    #[serde(default)]
    pub name: Option<String>,
    /// Optional thickness for 3D extrusion
    #[serde(default)]
    pub thickness: Option<f64>,
    /// Optional height offset (z-position)
    #[serde(default)]
    pub z_offset: Option<f64>,
    /// Optional color for visualization (RGBA, 0.0-1.0 range)
    #[serde(default)]
    pub color: Option<[f32; 4]>,
    /// Optional transparency (0.0-1.0)
    #[serde(default)]
    pub transparency: Option<f32>,
}

impl UnifiedLayerConfig {
    /// Creates a new layer config with the given layer and datatype.
    pub fn new(layer: u16, datatype: u16) -> Self {
        Self {
            layer,
            datatype,
            name: None,
            thickness: None,
            z_offset: None,
            color: None,
            transparency: None,
        }
    }

    /// Sets the layer name.
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Sets the thickness for 3D extrusion.
    pub fn with_thickness(mut self, thickness: f64) -> Self {
        self.thickness = Some(thickness);
        self
    }

    /// Sets the z-offset (height position).
    pub fn with_z_offset(mut self, z_offset: f64) -> Self {
        self.z_offset = Some(z_offset);
        self
    }

    /// Sets the visualization color (RGB values 0-255).
    pub fn with_color(mut self, r: u8, g: u8, b: u8, a: u8) -> Self {
        self.color = Some([
            r as f32 / 255.0,
            g as f32 / 255.0,
            b as f32 / 255.0,
            a as f32 / 255.0,
        ]);
        self
    }

    /// Sets the color from f32 values (0.0-1.0 range).
    pub fn with_color_f32(mut self, r: f32, g: f32, b: f32, a: f32) -> Self {
        self.color = Some([r, g, b, a]);
        self
    }

    /// Sets the transparency (0.0 = opaque, 1.0 = fully transparent).
    pub fn with_transparency(mut self, transparency: f32) -> Self {
        self.transparency = Some(transparency.clamp(0.0, 1.0));
        self
    }

    /// Converts to GDS layer settings.
    pub fn to_gds_settings(&self) -> rcad_gds::LayerSettings {
        rcad_gds::LayerSettings {
            thickness: self.thickness.unwrap_or(1.0),
            z_offset: self.z_offset.unwrap_or(0.0),
            color: self.color,
            name: self.name.clone(),
        }
    }

    /// Converts to OASIS layer settings.
    pub fn to_oas_settings(&self) -> rcad_oas::LayerSettings {
        rcad_oas::LayerSettings {
            thickness: self.thickness.unwrap_or(1.0),
            z_offset: self.z_offset.unwrap_or(0.0),
            color: self.color,
            name: self.name.clone(),
        }
    }
}

impl From<(u16, u16)> for UnifiedLayerConfig {
    fn from((layer, datatype): (u16, u16)) -> Self {
        Self::new(layer, datatype)
    }
}

impl Default for UnifiedLayerConfig {
    fn default() -> Self {
        Self::new(0, 0)
    }
}

/// A collection of layer configurations.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct LayerConfigMap {
    /// Map from (layer, datatype) to configuration.
    #[serde(default)]
    pub layers: std::collections::BTreeMap<(u16, u16), UnifiedLayerConfig>,
}

impl LayerConfigMap {
    /// Creates an empty layer config map.
    pub fn new() -> Self {
        Self {
            layers: std::collections::BTreeMap::new(),
        }
    }

    /// Adds a layer configuration.
    pub fn add(&mut self, config: UnifiedLayerConfig) {
        self.layers.insert((config.layer, config.datatype), config);
    }

    /// Gets a layer configuration.
    pub fn get(&self, layer: u16, datatype: u16) -> Option<&UnifiedLayerConfig> {
        self.layers.get(&(layer, datatype))
    }

    /// Returns the number of configured layers.
    pub fn len(&self) -> usize {
        self.layers.len()
    }

    /// Returns true if there are no configured layers.
    pub fn is_empty(&self) -> bool {
        self.layers.is_empty()
    }

    /// Loads layer configurations from a JSON file.
    pub fn from_json_file(path: &std::path::Path) -> Result<Self, crate::error::IoError> {
        let content = std::fs::read_to_string(path)?;
        let config: Self = serde_json::from_str(&content)?;
        Ok(config)
    }

    /// Saves layer configurations to a JSON file.
    pub fn to_json_file(&self, path: &std::path::Path) -> Result<(), crate::error::IoError> {
        let content = serde_json::to_string_pretty(self)?;
        std::fs::write(path, content)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_layer_config_new() {
        let config = UnifiedLayerConfig::new(1, 0);
        assert_eq!(config.layer, 1);
        assert_eq!(config.datatype, 0);
        assert!(config.name.is_none());
    }

    #[test]
    fn test_layer_config_builder() {
        let config = UnifiedLayerConfig::new(1, 0)
            .with_name("metal1")
            .with_thickness(0.5)
            .with_z_offset(1.0)
            .with_color(255, 0, 0, 255)
            .with_transparency(0.5);

        assert_eq!(config.name, Some("metal1".to_string()));
        assert_eq!(config.thickness, Some(0.5));
        assert_eq!(config.z_offset, Some(1.0));
        assert_eq!(config.color, Some([1.0, 0.0, 0.0, 1.0]));
        assert_eq!(config.transparency, Some(0.5));
    }

    #[test]
    fn test_layer_config_from_tuple() {
        let config: UnifiedLayerConfig = (1, 2).into();
        assert_eq!(config.layer, 1);
        assert_eq!(config.datatype, 2);
    }

    #[test]
    fn test_layer_config_map() {
        let mut map = LayerConfigMap::new();
        map.add(UnifiedLayerConfig::new(1, 0));
        map.add(UnifiedLayerConfig::new(2, 0));

        assert_eq!(map.len(), 2);
        assert!(map.get(1, 0).is_some());
        assert!(map.get(3, 0).is_none());
    }

    #[test]
    fn test_layer_config_serialization() {
        let config = UnifiedLayerConfig::new(1, 0)
            .with_name("metal1")
            .with_thickness(0.5);

        let json = serde_json::to_string(&config).unwrap();
        let decoded: UnifiedLayerConfig = serde_json::from_str(&json).unwrap();

        assert_eq!(config, decoded);
    }

    #[test]
    fn test_to_gds_settings() {
        let config = UnifiedLayerConfig::new(1, 0)
            .with_name("metal1")
            .with_thickness(0.5)
            .with_color_f32(1.0, 0.0, 0.0, 1.0);

        let settings = config.to_gds_settings();
        assert_eq!(settings.thickness, 0.5);
        assert_eq!(settings.name, Some("metal1".to_string()));
        assert_eq!(settings.color, Some([1.0, 0.0, 0.0, 1.0]));
    }

    #[test]
    fn test_to_oas_settings() {
        let config = UnifiedLayerConfig::new(1, 0)
            .with_name("metal1")
            .with_thickness(0.5);

        let settings = config.to_oas_settings();
        assert_eq!(settings.thickness, 0.5);
        assert_eq!(settings.name, Some("metal1".to_string()));
    }
}
