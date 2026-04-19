//! Environment and lighting environment.

use glam::Vec3;
use serde::{Deserialize, Serialize};

/// Environment for scene lighting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Environment {
    pub env_type: EnvironmentType,
    pub intensity: f32,
    pub rotation: f32,
}

/// Types of environment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EnvironmentType {
    Constant { color: Vec3 },
    Hdri(HdriEnvironment),
    Sky(SkyEnvironment),
}

/// HDRI environment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HdriEnvironment {
    pub file_path: String,
}

/// Procedural sky environment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SkyEnvironment {
    pub sun_direction: Vec3,
    pub turbidity: f32,
    pub ground_albedo: f32,
}

impl Default for Environment {
    fn default() -> Self {
        Self {
            env_type: EnvironmentType::Constant {
                color: Vec3::new(0.5, 0.5, 0.5),
            },
            intensity: 1.0,
            rotation: 0.0,
        }
    }
}

impl Environment {
    pub fn constant(color: Vec3) -> Self {
        Self {
            env_type: EnvironmentType::Constant { color },
            ..Default::default()
        }
    }

    pub fn hdri(file_path: &str) -> Self {
        Self {
            env_type: EnvironmentType::Hdri(HdriEnvironment {
                file_path: file_path.to_string(),
            }),
            ..Default::default()
        }
    }

    pub fn sky(sun_direction: Vec3) -> Self {
        Self {
            env_type: EnvironmentType::Sky(SkyEnvironment {
                sun_direction: sun_direction.normalize_or_zero(),
                turbidity: 2.0,
                ground_albedo: 0.1,
            }),
            ..Default::default()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_environment() {
        let env = Environment::default();
        assert!(matches!(env.env_type, EnvironmentType::Constant { .. }));
    }

    #[test]
    fn test_sky_environment() {
        let env = Environment::sky(Vec3::new(0.0, -1.0, 0.5));
        assert!(matches!(env.env_type, EnvironmentType::Sky(_)));
    }
}
