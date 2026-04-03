//! User-facing geometry construction helpers.
//!
//! This module is the public modeling entry layer for RCAD.
//! The API intentionally prefers OCCT-style direct constructor functions
//! over fluent builder structs.

mod curve;
mod solid;
mod surface;
pub mod brep_builder;
pub mod fillet;
pub mod ops;

pub use curve::*;
pub use solid::*;
pub use surface::*;
pub use brep_builder::*;
pub use fillet::{chamfer_edge, fillet_edge};
pub use ops::*;

use glam::DVec3;
use rcad_kernel::BRep;
use std::error::Error;
use std::fmt;

const EPS: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildError {
    NonFiniteValue(&'static str),
    NonPositiveValue(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    DegenerateGeometry(&'static str),
    InvalidIndex(usize),
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteValue(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveValue(name) => write!(f, "{name} must be > 0"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::DegenerateGeometry(msg) => write!(f, "degenerate geometry: {msg}"),
            Self::InvalidIndex(idx) => write!(f, "invalid index: {idx}"),
        }
    }
}

impl Error for BuildError {}

fn validate_point(name: &'static str, point: DVec3) -> Result<DVec3, BuildError> {
    if point.is_finite() {
        Ok(point)
    } else {
        Err(BuildError::NonFiniteValue(name))
    }
}

fn validate_positive(name: &'static str, value: f64) -> Result<f64, BuildError> {
    if !value.is_finite() {
        Err(BuildError::NonFiniteValue(name))
    } else if value <= 0.0 {
        Err(BuildError::NonPositiveValue(name))
    } else {
        Ok(value)
    }
}

fn normalize_vector(name: &'static str, vector: DVec3) -> Result<DVec3, BuildError> {
    validate_point(name, vector)?;
    if vector.length_squared() <= EPS {
        Err(BuildError::ZeroVector(name))
    } else {
        Ok(vector.normalize())
    }
}

fn normalize_rejection(
    name: &'static str,
    vector: DVec3,
    reference_name: &'static str,
    reference: DVec3,
) -> Result<DVec3, BuildError> {
    let vector = normalize_vector(name, vector)?;
    let rejected = vector - reference * vector.dot(reference);
    if rejected.length_squared() <= EPS {
        Err(BuildError::ParallelVectors(name, reference_name))
    } else {
        Ok(rejected.normalize())
    }
}

fn basis_from_x_y(x_dir: DVec3, y_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BuildError> {
    let x_axis = normalize_vector("x_dir", x_dir)?;
    let y_axis = normalize_rejection("y_dir", y_dir, "x_dir", x_axis)?;
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

fn basis_from_axis_ref(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BuildError> {
    let y_axis = normalize_vector("axis", axis)?;
    let x_axis = normalize_rejection("ref_dir", ref_dir, "axis", y_axis)?;
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

fn translate_brep(brep: &mut BRep, offset: DVec3) {
    for vertex in &mut brep.vertices {
        vertex.point += offset;
    }
    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                face.normal = face.normal.normalize_or_zero();
            }
        }
    }
}

fn transform_brep(brep: &mut BRep, origin: DVec3, x_axis: DVec3, y_axis: DVec3, z_axis: DVec3) {
    for vertex in &mut brep.vertices {
        vertex.point = origin
            + x_axis * vertex.point.x
            + y_axis * vertex.point.y
            + z_axis * vertex.point.z;
    }

    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                let transformed = x_axis * face.normal.x + y_axis * face.normal.y + z_axis * face.normal.z;
                face.normal = transformed.normalize_or_zero();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{PrimitiveSolid, Surface3};

    #[test]
    fn line_rejects_zero_direction() {
        let err = line(DVec3::ZERO, DVec3::ZERO).unwrap_err();
        assert_eq!(err, BuildError::ZeroVector("direction"));
    }

    #[test]
    fn ellipse_rejects_parallel_major_direction() {
        let err = ellipse(DVec3::ZERO, DVec3::Z, DVec3::Z, 2.0, 1.0).unwrap_err();
        assert_eq!(err, BuildError::ParallelVectors("major_dir", "normal"));
    }

    #[test]
    fn box_brep_builds_transformed_vertices() {
        let brep = box_brep(DVec3::new(1.0, 2.0, 3.0), DVec3::Y, DVec3::Z, 2.0, 3.0, 4.0).unwrap();

        assert_eq!(brep.vertices.len(), 8);
        assert!(brep.vertices.iter().any(|v| v.point == DVec3::new(1.0, 2.0, 3.0)));
        assert!(brep.vertices.iter().any(|v| v.point == DVec3::new(5.0, 4.0, 6.0)));
    }

    #[test]
    fn sphere_brep_translates_bounds() {
        let brep = sphere_brep(DVec3::new(10.0, -2.0, 4.0), 2.0).unwrap();

        let min_y = brep
            .vertices
            .iter()
            .map(|v| v.point.y)
            .fold(f64::INFINITY, f64::min);
        let max_y = brep
            .vertices
            .iter()
            .map(|v| v.point.y)
            .fold(f64::NEG_INFINITY, f64::max);

        assert!((min_y - (-4.0)).abs() < 1e-6);
        assert!((max_y - 0.0).abs() < 1e-6);
    }

    #[test]
    fn cylinder_primitive_returns_expected_shape() {
        let primitive = cylinder_primitive(3.0, 5.0).unwrap();

        match primitive {
            PrimitiveSolid::Cylinder { radius, height } => {
                assert_eq!(radius, 3.0);
                assert_eq!(height, 5.0);
            }
            other => panic!("expected cylinder primitive, got {other:?}"),
        }
    }

    #[test]
    fn make_plane_alias_matches_plane_constructor() {
        let surface = make_plane(DVec3::new(1.0, 2.0, 3.0), DVec3::Z).unwrap();

        match surface {
            Surface3::Plane(plane) => {
                assert_eq!(plane.origin, DVec3::new(1.0, 2.0, 3.0));
                assert_eq!(plane.normal, DVec3::Z);
            }
            other => panic!("expected plane surface, got {other:?}"),
        }
    }
}