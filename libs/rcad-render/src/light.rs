//! Light source definitions.

use glam::Vec3;
use serde::{Deserialize, Serialize};

/// Light identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct LightId(pub u64);

/// Light source.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Light {
    pub id: LightId,
    pub name: String,
    pub light_type: LightType,
    pub color: Vec3,
    pub intensity: f32,
    pub enabled: bool,
}

/// Types of light sources.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LightType {
    /// Directional light (sun).
    Directional { direction: Vec3 },
    /// Point light.
    Point { position: Vec3, radius: f32 },
    /// Spot light.
    Spot {
        position: Vec3,
        direction: Vec3,
        inner_angle: f32,
        outer_angle: f32,
    },
    /// Area light.
    Area {
        position: Vec3,
        normal: Vec3,
        width: f32,
        height: f32,
    },
}

impl Light {
    /// Create a directional light.
    pub fn directional(id: LightId, name: &str, direction: Vec3) -> Self {
        Self {
            id,
            name: name.to_string(),
            light_type: LightType::Directional {
                direction: direction.normalize_or_zero(),
            },
            color: Vec3::new(1.0, 1.0, 1.0),
            intensity: 1.0,
            enabled: true,
        }
    }

    /// Create a point light.
    pub fn point(id: LightId, name: &str, position: Vec3) -> Self {
        Self {
            id,
            name: name.to_string(),
            light_type: LightType::Point { position, radius: 0.1 },
            color: Vec3::new(1.0, 1.0, 1.0),
            intensity: 1.0,
            enabled: true,
        }
    }

    /// Create a spot light.
    pub fn spot(
        id: LightId,
        name: &str,
        position: Vec3,
        direction: Vec3,
        inner_angle_deg: f32,
        outer_angle_deg: f32,
    ) -> Self {
        Self {
            id,
            name: name.to_string(),
            light_type: LightType::Spot {
                position,
                direction: direction.normalize_or_zero(),
                inner_angle: inner_angle_deg.to_radians(),
                outer_angle: outer_angle_deg.to_radians(),
            },
            color: Vec3::new(1.0, 1.0, 1.0),
            intensity: 1.0,
            enabled: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_directional_light() {
        let light = Light::directional(LightId(1), "Sun", Vec3::new(0.0, -1.0, 0.0));
        assert!(matches!(light.light_type, LightType::Directional { .. }));
    }

    #[test]
    fn test_point_light() {
        let light = Light::point(LightId(1), "Bulb", Vec3::new(0.0, 2.0, 0.0));
        assert!(matches!(light.light_type, LightType::Point { .. }));
    }
}
