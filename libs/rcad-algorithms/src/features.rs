//! First-stage feature operations (TKFeat-like APIs).
//!
//! This module builds practical feature workflows on top of the existing
//! boolean kernel. The first shipped feature is a cylindrical hole.

use glam::{DAffine3, DMat3, DVec3};
use rcad_kernel::{BRep, PrimitiveSolid};

use crate::{BooleanError, BooleanOpType, boolean_op};

/// Errors returned by feature operations.
#[derive(Debug)]
pub enum FeatureError {
    NonFiniteInput(&'static str),
    NonPositiveInput(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    Boolean(BooleanError),
}

impl std::fmt::Display for FeatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFiniteInput(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveInput(name) => write!(f, "{name} must be > 0"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::Boolean(err) => write!(f, "boolean operation failed: {err}"),
        }
    }
}

impl std::error::Error for FeatureError {}

impl From<BooleanError> for FeatureError {
    fn from(value: BooleanError) -> Self {
        Self::Boolean(value)
    }
}

const EPS: f64 = 1e-12;

fn validate_finite(name: &'static str, v: f64) -> Result<f64, FeatureError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(FeatureError::NonFiniteInput(name))
    }
}

fn validate_positive(name: &'static str, v: f64) -> Result<f64, FeatureError> {
    let v = validate_finite(name, v)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(FeatureError::NonPositiveInput(name))
    }
}

fn normalize(name: &'static str, v: DVec3) -> Result<DVec3, FeatureError> {
    if !v.is_finite() {
        return Err(FeatureError::NonFiniteInput(name));
    }
    if v.length_squared() <= EPS {
        return Err(FeatureError::ZeroVector(name));
    }
    Ok(v.normalize())
}

fn axis_ref_basis(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), FeatureError> {
    let y_axis = normalize("axis", axis)?;
    let ref_dir = normalize("ref_dir", ref_dir)?;
    let x_reject = ref_dir - y_axis * ref_dir.dot(y_axis);
    if x_reject.length_squared() <= EPS {
        return Err(FeatureError::ParallelVectors("ref_dir", "axis"));
    }
    let x_axis = x_reject.normalize();
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

/// Create a cylindrical through/blind hole by subtracting an oriented cylinder
/// from `target`.
///
/// - `center`: center of the tool cylinder.
/// - `axis`: cylinder axis direction.
/// - `ref_dir`: reference direction used to build local orientation.
/// - `radius`: hole radius.
/// - `depth`: tool cylinder height.
///
/// For through holes, pass a `depth` larger than the part thickness along
/// `axis`.
pub fn make_cylindrical_hole(
    target: &BRep,
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    depth: f64,
) -> Result<BRep, FeatureError> {
    if !center.is_finite() {
        return Err(FeatureError::NonFiniteInput("center"));
    }
    let radius = validate_positive("radius", radius)?;
    let depth = validate_positive("depth", depth)?;

    let (x_axis, y_axis, z_axis) = axis_ref_basis(axis, ref_dir)?;

    let mut tool = BRep::from_primitive(PrimitiveSolid::Cylinder { radius, height: depth });
    let rot = DMat3::from_cols(x_axis, y_axis, z_axis);
    tool.apply_transform(DAffine3::from_mat3_translation(rot, center));

    Ok(boolean_op(BooleanOpType::Difference, target, &tool)?)
}

#[cfg(test)]
mod tests {
    use glam::DVec3;
    use rcad_kernel::{BRep, PrimitiveSolid};

    use super::*;
    #[test]
    fn cylindrical_hole_subtracts_from_box() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        });

        let result = make_cylindrical_hole(
            &target,
            DVec3::ZERO,
            DVec3::Y,
            DVec3::X,
            0.6,
            6.0,
        )
        .expect("cylindrical hole should succeed");

        assert!(
            result.solids[0].shells[0].faces.len() >= target.solids[0].shells[0].faces.len(),
            "hole operation should keep or increase face count"
        );
        assert!(!result.edges.is_empty(), "hole result should keep edge topology");
    }

    #[test]
    fn cylindrical_hole_rejects_parallel_axis_ref_dir() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        let err = make_cylindrical_hole(
            &target,
            DVec3::ZERO,
            DVec3::Y,
            DVec3::Y,
            0.3,
            3.0,
        )
        .expect_err("parallel axis/ref_dir must be rejected");

        assert!(matches!(err, FeatureError::ParallelVectors(_, _)));
    }
}
