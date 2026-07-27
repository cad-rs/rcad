use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Settings for a single layer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerSettings {
    /// Extrusion thickness in user units.
    pub thickness: f64,
    /// Z offset for this layer.
    #[serde(default)]
    pub z_offset: f64,
    /// Optional RGBA color (0.0-1.0 range).
    #[serde(default)]
    pub color: Option<[f32; 4]>,
    /// Optional layer name.
    #[serde(default)]
    pub name: Option<String>,
}

impl Default for LayerSettings {
    fn default() -> Self {
        Self {
            thickness: 1.0,
            z_offset: 0.0,
            color: None,
            name: None,
        }
    }
}

impl LayerSettings {
    pub fn new(thickness: f64) -> Self {
        Self {
            thickness,
            ..Default::default()
        }
    }
}

/// Layer configuration for OAS-to-3D conversion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerConfig {
    /// Mapping from layer number to settings.
    pub layers: HashMap<i32, LayerSettings>,
    /// Default thickness for unconfigured layers.
    #[serde(default = "default_default_thickness")]
    pub default_thickness: f64,
}

fn default_default_thickness() -> f64 {
    1.0
}

impl Default for LayerConfig {
    fn default() -> Self {
        Self {
            layers: HashMap::new(),
            default_thickness: 1.0,
        }
    }
}

impl LayerConfig {
    /// Create an empty config with default thickness.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create config from a JSON file.
    pub fn from_json_file<P: AsRef<std::path::Path>>(
        path: P,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(path)?;
        Self::from_json(&content).map_err(|e| e.into())
    }

    /// Create config from a JSON string.
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    /// Export config to JSON string.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Add a layer setting.
    pub fn with_layer(mut self, layer: i32, settings: LayerSettings) -> Self {
        self.layers.insert(layer, settings);
        self
    }

    /// Get settings for a layer, using defaults if not configured.
    pub fn get(&self, layer: i32) -> LayerSettings {
        self.layers
            .get(&layer)
            .cloned()
            .unwrap_or_else(|| LayerSettings {
                thickness: self.default_thickness,
                ..Default::default()
            })
    }

    /// Get the thickness for a layer.
    pub fn thickness(&self, layer: i32) -> f64 {
        self.get(layer).thickness
    }

    /// Check if a layer is explicitly configured.
    pub fn is_configured(&self, layer: i32) -> bool {
        self.layers.contains_key(&layer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_layer_config_default() {
        let config = LayerConfig::default();
        assert_eq!(config.thickness(1), 1.0);
        assert_eq!(config.thickness(999), 1.0);
    }

    #[test]
    fn test_layer_config_with_layer() {
        let config = LayerConfig::new()
            .with_layer(1, LayerSettings::new(2.0))
            .with_layer(2, LayerSettings::new(3.0));

        assert_eq!(config.thickness(1), 2.0);
        assert_eq!(config.thickness(2), 3.0);
        assert_eq!(config.thickness(3), 1.0); // default
        assert!(config.is_configured(1));
        assert!(!config.is_configured(3));
    }

    #[test]
    fn test_layer_config_json() {
        let json = r#"{
            "layers": {
                "1": { "thickness": 0.5, "z_offset": 0.0 },
                "2": { "thickness": 1.0, "name": "metal" }
            },
            "default_thickness": 0.2
        }"#;

        let config = LayerConfig::from_json(json).unwrap();
        assert_eq!(config.thickness(1), 0.5);
        assert_eq!(config.thickness(2), 1.0);
        assert_eq!(config.thickness(3), 0.2);
        assert_eq!(
            config.layers.get(&2).unwrap().name,
            Some("metal".to_string())
        );
    }
}
