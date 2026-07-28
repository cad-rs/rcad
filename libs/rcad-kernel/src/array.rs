//! Array (pattern) operations  ?linear and circular repetition of BRep solids.
//!
//! Analogous to OCCT `BRepOffsetAPI_MakeThickSolid`-style patterns and
//! `BRepFeat_MakeLinearForm` / `BRepFeat_MakeRevol` for feature repetition.
//!
//! # Operations
//!
//! - **Linear pattern**: repeat along a direction with uniform spacing
//! - **Circular pattern**: rotate around an axis with uniform angular spacing
//! - **Mirror pattern**: mirror across a plane, optionally including original
//! - **Rectangular grid pattern**: 2D array with staggered options
//! - **Variable spacing pattern**: non-uniform spacing along a direction
//! - **Path pattern**: pattern along a curve with optional alignment

use crate::geom::{
    Circle3, ConicalSurface, Curve3, CylindricalSurface, Ellipse3, Hyperbola3, Line3,
    LinearExtrusionSurface, OffsetSurface, Plane, RevolutionSurface, SphericalSurface, Surface3,
    ToroidalSurface, TrimmedSurface,
};
use crate::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
use crate::{BRep, any_perpendicular, topods};
use glam::{DMat4, DVec3};

/// Parameters for a linear pattern.
#[derive(Debug, Clone)]
pub struct LinearPatternParams {
    /// Direction of the pattern.
    pub direction: DVec3,
    /// Number of copies (including the original). Must be >= 1.
    pub count: usize,
    /// Spacing between consecutive copies.
    pub spacing: f64,
}

/// Parameters for a circular pattern.
#[derive(Debug, Clone)]
pub struct CircularPatternParams {
    /// A point on the rotation axis.
    pub axis_origin: DVec3,
    /// Normalized rotation axis direction.
    pub axis_direction: DVec3,
    /// Number of copies (including the original). Must be >= 1.
    pub count: usize,
    /// Total angle in radians for the full pattern (copies are evenly spaced).
    pub total_angle: f64,
}

/// Parameters for a mirror pattern.
#[derive(Debug, Clone)]
pub struct MirrorPatternParams {
    /// Origin point on the mirror plane.
    pub plane_origin: DVec3,
    /// Normal vector of the mirror plane (will be normalized).
    pub plane_normal: DVec3,
    /// Whether to include the original shape in the result.
    pub include_original: bool,
}

/// Parameters for a rectangular grid pattern.
#[derive(Debug, Clone)]
pub struct RectangularPatternParams {
    /// Direction for the first (X) axis of the pattern.
    pub direction1: DVec3,
    /// Number of copies along direction1 (including original).
    pub count1: usize,
    /// Spacing between copies along direction1.
    pub spacing1: f64,
    /// Direction for the second (Y) axis of the pattern.
    pub direction2: DVec3,
    /// Number of copies along direction2 (including original).
    pub count2: usize,
    /// Spacing between copies along direction2.
    pub spacing2: f64,
    /// Stagger pattern configuration.
    pub stagger: StaggerConfig,
}

/// Stagger configuration for rectangular patterns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StaggerConfig {
    /// No staggering - regular grid.
    #[default]
    None,
    /// Offset odd rows by half spacing1.
    OddRows,
    /// Offset even rows by half spacing1.
    EvenRows,
}

/// Parameters for variable spacing patterns.
#[derive(Debug, Clone)]
pub struct VariableSpacingPatternParams {
    /// Direction of the pattern.
    pub direction: DVec3,
    /// List of spacings between consecutive copies.
    /// Number of copies = spacings.len() + 1 (includes original).
    pub spacings: Vec<f64>,
}

/// Parameters for distance-based spacing pattern.
#[derive(Debug, Clone)]
pub struct DistanceSpacingPatternParams {
    /// Direction of the pattern.
    pub direction: DVec3,
    /// Total distance to cover.
    pub total_distance: f64,
    /// Number of copies (including the original).
    pub count: usize,
}

/// Parameters for a path-based pattern.
#[derive(Debug, Clone)]
pub struct PathPatternParams {
    /// List of parameter values (0.0 to 1.0) along the path where copies are placed.
    /// 0.0 = start of path, 1.0 = end of path.
    pub parameters: Vec<f64>,
    /// Whether to align instances with the path tangent.
    pub align_to_path: bool,
    /// Up vector for alignment when align_to_path is true.
    pub up_vector: DVec3,
}

/// Parameters for a pattern with suppression support.
#[derive(Debug, Clone)]
pub struct PatternWithSuppressionParams {
    /// Indices of instances to suppress (0-indexed, 0 = original).
    pub suppressed_indices: Vec<usize>,
    /// Base pattern transformation matrices for each instance.
    pub transforms: Vec<DMat4>,
}

/// Parameters for pattern within a boundary.
#[derive(Debug, Clone)]
pub struct BoundaryPatternParams {
    /// Grid direction 1.
    pub direction1: DVec3,
    /// Grid count 1.
    pub count1: usize,
    /// Grid spacing 1.
    pub spacing1: f64,
    /// Grid direction 2.
    pub direction2: DVec3,
    /// Grid count 2.
    pub spacing2: f64,
    /// Function to test if a point is within the boundary.
    /// Returns true if the point (grid position) should include an instance.
    pub boundary_test: fn(DVec3) -> bool,
}

/// Instance-specific transformation for patterns.
#[derive(Debug, Clone)]
pub struct InstanceTransform {
    /// Index of the instance (0 = original).
    pub index: usize,
    /// Additional transformation to apply to this instance.
    pub transform: DMat4,
}

/// Error type for pattern operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PatternError {
    /// Count must be at least 1.
    InvalidCount,
    /// Spacing must be positive.
    InvalidSpacing,
    /// Direction vector must be non-zero.
    ZeroDirection,
    /// Axis direction must be non-zero.
    ZeroAxis,
    /// Total angle must be non-zero and <= 2*pi.
    InvalidAngle,
    /// Input BRep has no solids.
    NoSolids,
    /// Plane normal must be non-zero.
    ZeroPlaneNormal,
    /// Spacings list must not be empty.
    EmptySpacings,
    /// Negative spacing in list.
    NegativeSpacing,
    /// Total distance must be positive.
    InvalidDistance,
    /// Parameters list must not be empty.
    EmptyParameters,
    /// Parameter must be in range [0, 1].
    InvalidParameter,
    /// Transform list must not be empty.
    EmptyTransforms,
    /// Suppressed index out of range.
    SuppressedIndexOutOfRange,
    /// Path curve evaluation failed.
    PathEvaluationFailed,
}

impl std::fmt::Display for PatternError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCount => write!(f, "pattern count must be >= 1"),
            Self::InvalidSpacing => write!(f, "pattern spacing must be > 0"),
            Self::ZeroDirection => write!(f, "pattern direction must be non-zero"),
            Self::ZeroAxis => write!(f, "pattern axis must be non-zero"),
            Self::InvalidAngle => write!(f, "pattern angle must be > 0 and <= 2*pi"),
            Self::NoSolids => write!(f, "input BRep has no solids"),
            Self::ZeroPlaneNormal => write!(f, "mirror plane normal must be non-zero"),
            Self::EmptySpacings => write!(f, "spacings list must not be empty"),
            Self::NegativeSpacing => write!(f, "all spacings must be >= 0"),
            Self::InvalidDistance => write!(f, "total distance must be > 0"),
            Self::EmptyParameters => write!(f, "parameters list must not be empty"),
            Self::InvalidParameter => write!(f, "parameters must be in range [0, 1]"),
            Self::EmptyTransforms => write!(f, "transform list must not be empty"),
            Self::SuppressedIndexOutOfRange => write!(f, "suppressed index out of range"),
            Self::PathEvaluationFailed => write!(f, "path curve evaluation failed"),
        }
    }
}

/// Apply a linear pattern to a BRep  ?repeat copies along a direction.
///
/// Returns a new BRep containing all copies merged into a single solid.
/// The original is included as the first copy (offset 0).
pub fn linear_pattern(brep: &BRep, params: &LinearPatternParams) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }
    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    for i in 0..params.count {
        let offset = dir * (i as f64 * params.spacing);
        append_transformed_brep(&mut out, brep, &translation_matrix(offset))?;
    }

    Ok(out)
}

/// Apply a circular pattern to a BRep  ?rotate copies around an axis.
///
/// Returns a new BRep containing all copies merged into a single solid.
/// The original is included as the first copy (angle 0).
pub fn circular_pattern(brep: &BRep, params: &CircularPatternParams) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.total_angle <= 0.0 || params.total_angle > std::f64::consts::TAU {
        return Err(PatternError::InvalidAngle);
    }
    let axis = params
        .axis_direction
        .try_normalize()
        .ok_or(PatternError::ZeroAxis)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();
    let angle_step = params.total_angle / params.count as f64;

    for i in 0..params.count {
        let angle = i as f64 * angle_step;
        let mat = rotation_matrix(params.axis_origin, axis, angle);
        append_transformed_brep(&mut out, brep, &mat)?;
    }

    Ok(out)
}

//  € € Mirror Pattern  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Apply a mirror pattern to a BRep  ?mirror across a plane.
///
/// Returns a new BRep containing the mirrored copy and optionally the original.
pub fn mirror_pattern(brep: &BRep, params: &MirrorPatternParams) -> Result<BRep, PatternError> {
    let normal = params
        .plane_normal
        .try_normalize()
        .ok_or(PatternError::ZeroPlaneNormal)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // Include original if requested
    if params.include_original {
        append_transformed_brep(&mut out, brep, &DMat4::IDENTITY)?;
    }

    // Mirror across the plane
    let mirror_mat = mirror_matrix(params.plane_origin, normal);
    append_transformed_brep(&mut out, brep, &mirror_mat)?;

    Ok(out)
}

/// Apply a compound mirror pattern  ?mirror and linear pattern combined.
///
/// Mirrors the shape across a plane, then applies a linear pattern to both
/// the original and mirrored copies.
pub fn mirror_linear_pattern(
    brep: &BRep,
    mirror_params: &MirrorPatternParams,
    linear_params: &LinearPatternParams,
) -> Result<BRep, PatternError> {
    // First apply mirror
    let mirrored = mirror_pattern(brep, mirror_params)?;

    // Then apply linear pattern to the mirrored result
    linear_pattern(&mirrored, linear_params)
}

/// Apply a compound mirror and circular pattern.
///
/// Mirrors the shape across a plane, then applies a circular pattern.
pub fn mirror_circular_pattern(
    brep: &BRep,
    mirror_params: &MirrorPatternParams,
    circular_params: &CircularPatternParams,
) -> Result<BRep, PatternError> {
    let mirrored = mirror_pattern(brep, mirror_params)?;
    circular_pattern(&mirrored, circular_params)
}

//  € € Rectangular Grid Pattern  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Apply a rectangular grid pattern to a BRep.
///
/// Creates a 2D grid of copies along two orthogonal (or non-orthogonal) directions.
/// Supports staggered patterns for alternating row offsets.
pub fn rectangular_pattern(
    brep: &BRep,
    params: &RectangularPatternParams,
) -> Result<BRep, PatternError> {
    if params.count1 < 1 || params.count2 < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing1 <= 0.0 || params.spacing2 <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir1 = params
        .direction1
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;
    let dir2 = params
        .direction2
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    for j in 0..params.count2 {
        let row_offset = dir2 * (j as f64 * params.spacing2);
        let stagger_offset = match params.stagger {
            StaggerConfig::None => DVec3::ZERO,
            StaggerConfig::OddRows if j % 2 == 1 => dir1 * (params.spacing1 * 0.5),
            StaggerConfig::EvenRows if j % 2 == 0 && j > 0 => dir1 * (params.spacing1 * 0.5),
            _ => DVec3::ZERO,
        };

        for i in 0..params.count1 {
            let col_offset = dir1 * (i as f64 * params.spacing1);
            let total_offset = row_offset + col_offset + stagger_offset;
            append_transformed_brep(&mut out, brep, &translation_matrix(total_offset))?;
        }
    }

    Ok(out)
}

/// Compute the transformation matrix for a specific position in a rectangular grid.
///
/// Returns the transformation matrix for the instance at (i, j) in the grid.
pub fn rectangular_pattern_transform(
    params: &RectangularPatternParams,
    i: usize,
    j: usize,
) -> Result<DMat4, PatternError> {
    if params.spacing1 <= 0.0 || params.spacing2 <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir1 = params
        .direction1
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;
    let dir2 = params
        .direction2
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    let row_offset = dir2 * (j as f64 * params.spacing2);
    let stagger_offset = match params.stagger {
        StaggerConfig::None => DVec3::ZERO,
        StaggerConfig::OddRows if j % 2 == 1 => dir1 * (params.spacing1 * 0.5),
        StaggerConfig::EvenRows if j.is_multiple_of(2) && j > 0 => dir1 * (params.spacing1 * 0.5),
        _ => DVec3::ZERO,
    };
    let col_offset = dir1 * (i as f64 * params.spacing1);
    let total_offset = row_offset + col_offset + stagger_offset;

    Ok(translation_matrix(total_offset))
}

//  € € Variable Spacing Pattern  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Apply a pattern with non-uniform spacing.
///
/// Creates copies at positions determined by cumulative spacings.
pub fn variable_spacing_pattern(
    brep: &BRep,
    params: &VariableSpacingPatternParams,
) -> Result<BRep, PatternError> {
    if params.spacings.is_empty() {
        return Err(PatternError::EmptySpacings);
    }
    if params.spacings.iter().any(|&s| s < 0.0) {
        return Err(PatternError::NegativeSpacing);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // First copy is at offset 0 (original)
    append_transformed_brep(&mut out, brep, &DMat4::IDENTITY)?;

    // Subsequent copies at cumulative offsets
    let mut cumulative = 0.0;
    for &spacing in &params.spacings {
        cumulative += spacing;
        let offset = dir * cumulative;
        append_transformed_brep(&mut out, brep, &translation_matrix(offset))?;
    }

    Ok(out)
}

/// Apply a pattern with distance-based spacing.
///
/// Distributes copies evenly along a total distance.
pub fn distance_spacing_pattern(
    brep: &BRep,
    params: &DistanceSpacingPatternParams,
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.total_distance <= 0.0 {
        return Err(PatternError::InvalidDistance);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // Spacing between copies
    let spacing = if params.count > 1 {
        params.total_distance / (params.count - 1) as f64
    } else {
        0.0
    };

    for i in 0..params.count {
        let offset = dir * (i as f64 * spacing);
        append_transformed_brep(&mut out, brep, &translation_matrix(offset))?;
    }

    Ok(out)
}

//  € € Path Pattern  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Trait for path evaluation used in path patterns.
pub trait PathEvaluator {
    /// Evaluate the path at parameter t (0.0 to 1.0).
    /// Returns (position, tangent) tuple.
    fn evaluate(&self, t: f64) -> Option<(DVec3, DVec3)>;
}

/// Apply a pattern along a path.
///
/// Places copies at specified parameter values along a path curve.
/// Optionally aligns instances with the path tangent.
pub fn path_pattern(
    brep: &BRep,
    params: &PathPatternParams,
    path: &dyn PathEvaluator,
) -> Result<BRep, PatternError> {
    if params.parameters.is_empty() {
        return Err(PatternError::EmptyParameters);
    }
    if params.parameters.iter().any(|&t| !(0.0..=1.0).contains(&t)) {
        return Err(PatternError::InvalidParameter);
    }

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();
    let up = params.up_vector.normalize_or(DVec3::Z);

    for &t in &params.parameters {
        let (position, tangent) = path.evaluate(t).ok_or(PatternError::PathEvaluationFailed)?;

        let mat = if params.align_to_path {
            let tangent = tangent.normalize_or(DVec3::X);
            // Create a coordinate frame
            let side = up.cross(tangent).normalize_or(DVec3::Y);
            let real_up = tangent.cross(side).normalize_or(up);

            // Build rotation matrix from basis vectors
            DMat4::from_cols(
                glam::DVec4::new(tangent.x, tangent.y, tangent.z, 0.0),
                glam::DVec4::new(side.x, side.y, side.z, 0.0),
                glam::DVec4::new(real_up.x, real_up.y, real_up.z, 0.0),
                glam::DVec4::new(position.x, position.y, position.z, 1.0),
            )
        } else {
            translation_matrix(position)
        };

        append_transformed_brep(&mut out, brep, &mat)?;
    }

    Ok(out)
}

/// Apply a pattern along a path with equal spacing.
///
/// Creates count copies evenly distributed along the path.
pub fn path_pattern_equal_spacing(
    brep: &BRep,
    path: &dyn PathEvaluator,
    count: usize,
    align_to_path: bool,
    up_vector: DVec3,
) -> Result<BRep, PatternError> {
    if count < 1 {
        return Err(PatternError::InvalidCount);
    }

    let params = PathPatternParams {
        parameters: (0..count)
            .map(|i| i as f64 / (count - 1).max(1) as f64)
            .collect(),
        align_to_path,
        up_vector,
    };

    path_pattern(brep, &params, path)
}

//  € € Pattern with Suppression  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Apply a pattern with instance suppression.
///
/// Creates copies but excludes instances at the specified indices.
pub fn pattern_with_suppression(
    brep: &BRep,
    transforms: &[DMat4],
    suppressed_indices: &[usize],
) -> Result<BRep, PatternError> {
    if transforms.is_empty() {
        return Err(PatternError::EmptyTransforms);
    }

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    // Check suppressed indices are valid
    let max_index = transforms.len() - 1;
    if suppressed_indices.iter().any(|&i| i > max_index) {
        return Err(PatternError::SuppressedIndexOutOfRange);
    }

    let mut out = BRep::new();
    let suppressed_set: std::collections::HashSet<usize> =
        suppressed_indices.iter().copied().collect();

    for (i, mat) in transforms.iter().enumerate() {
        if !suppressed_set.contains(&i) {
            append_transformed_brep(&mut out, brep, mat)?;
        }
    }

    Ok(out)
}

/// Apply a linear pattern with suppression support.
pub fn linear_pattern_with_suppression(
    brep: &BRep,
    params: &LinearPatternParams,
    suppressed_indices: &[usize],
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let transforms: Vec<DMat4> = (0..params.count)
        .map(|i| translation_matrix(dir * (i as f64 * params.spacing)))
        .collect();

    pattern_with_suppression(brep, &transforms, suppressed_indices)
}

/// Apply a circular pattern with suppression support.
pub fn circular_pattern_with_suppression(
    brep: &BRep,
    params: &CircularPatternParams,
    suppressed_indices: &[usize],
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.total_angle <= 0.0 || params.total_angle > std::f64::consts::TAU {
        return Err(PatternError::InvalidAngle);
    }

    let axis = params
        .axis_direction
        .try_normalize()
        .ok_or(PatternError::ZeroAxis)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let angle_step = params.total_angle / params.count as f64;
    let transforms: Vec<DMat4> = (0..params.count)
        .map(|i| rotation_matrix(params.axis_origin, axis, i as f64 * angle_step))
        .collect();

    pattern_with_suppression(brep, &transforms, suppressed_indices)
}

//  € € Pattern with Instance Transforms  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Apply a pattern with instance-specific additional transformations.
///
/// Each instance can have its own additional transformation applied
/// on top of the base pattern transformation.
pub fn pattern_with_instance_transforms(
    brep: &BRep,
    base_transforms: &[DMat4],
    instance_transforms: &[InstanceTransform],
) -> Result<BRep, PatternError> {
    if base_transforms.is_empty() {
        return Err(PatternError::EmptyTransforms);
    }

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    // Build a map from index to instance transform
    let instance_map: std::collections::HashMap<usize, DMat4> = instance_transforms
        .iter()
        .map(|it| (it.index, it.transform))
        .collect();

    let mut out = BRep::new();

    for (i, base_mat) in base_transforms.iter().enumerate() {
        let final_mat = if let Some(instance_mat) = instance_map.get(&i) {
            *base_mat * *instance_mat
        } else {
            *base_mat
        };
        append_transformed_brep(&mut out, brep, &final_mat)?;
    }

    Ok(out)
}

//  € € Pattern within Boundary  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Apply a rectangular grid pattern constrained within a boundary.
///
/// Only creates instances whose grid positions fall within the boundary.
pub fn pattern_within_boundary(
    brep: &BRep,
    params: &BoundaryPatternParams,
) -> Result<BRep, PatternError> {
    if params.count1 < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing1 <= 0.0 || params.spacing2 <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir1 = params
        .direction1
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;
    let dir2 = params
        .direction2
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if !brep.has_solids() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // Estimate count2 from boundary test
    // We need to iterate until boundary test fails consistently
    let mut j = 0;
    let mut found_in_row;

    loop {
        found_in_row = false;
        let row_offset = dir2 * (j as f64 * params.spacing2);

        for i in 0..params.count1 {
            let col_offset = dir1 * (i as f64 * params.spacing1);
            let grid_pos = row_offset + col_offset;

            if (params.boundary_test)(grid_pos) {
                found_in_row = true;
                append_transformed_brep(&mut out, brep, &translation_matrix(grid_pos))?;
            }
        }

        // Stop if no instances were found in this row and we've gone past reasonable bounds
        if !found_in_row && j > 0 {
            break;
        }
        j += 1;

        // Safety limit to prevent infinite loops
        if j > 1000 {
            break;
        }
    }

    Ok(out)
}

//  € € Utility Functions  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Create a mirror transformation matrix across a plane.
fn mirror_matrix(plane_origin: DVec3, plane_normal: DVec3) -> DMat4 {
    // Mirror matrix: reflect across plane
    // M = I - 2 * n * n^T for reflection through origin
    // For arbitrary plane: translate to origin, reflect, translate back
    let n = plane_normal;
    let reflect = DMat4::from_cols(
        glam::DVec4::new(
            1.0 - 2.0 * n.x * n.x,
            -2.0 * n.y * n.x,
            -2.0 * n.z * n.x,
            0.0,
        ),
        glam::DVec4::new(
            -2.0 * n.x * n.y,
            1.0 - 2.0 * n.y * n.y,
            -2.0 * n.z * n.y,
            0.0,
        ),
        glam::DVec4::new(
            -2.0 * n.x * n.z,
            -2.0 * n.y * n.z,
            1.0 - 2.0 * n.z * n.z,
            0.0,
        ),
        glam::DVec4::new(0.0, 0.0, 0.0, 1.0),
    );

    DMat4::from_translation(plane_origin) * reflect * DMat4::from_translation(-plane_origin)
}

/// Generate transformation matrices for a linear pattern.
pub fn generate_linear_transforms(
    params: &LinearPatternParams,
) -> Result<Vec<DMat4>, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    Ok((0..params.count)
        .map(|i| translation_matrix(dir * (i as f64 * params.spacing)))
        .collect())
}

/// Generate transformation matrices for a circular pattern.
pub fn generate_circular_transforms(
    params: &CircularPatternParams,
) -> Result<Vec<DMat4>, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }

    let axis = params
        .axis_direction
        .try_normalize()
        .ok_or(PatternError::ZeroAxis)?;

    let angle_step = params.total_angle / params.count as f64;

    Ok((0..params.count)
        .map(|i| rotation_matrix(params.axis_origin, axis, i as f64 * angle_step))
        .collect())
}

/// Scale a pattern by uniformly scaling all spacing values.
pub fn scale_pattern_params(params: &LinearPatternParams, scale: f64) -> LinearPatternParams {
    LinearPatternParams {
        direction: params.direction,
        count: params.count,
        spacing: params.spacing * scale,
    }
}

/// Scale a rectangular pattern.
pub fn scale_rectangular_params(
    params: &RectangularPatternParams,
    scale: f64,
) -> RectangularPatternParams {
    RectangularPatternParams {
        direction1: params.direction1,
        count1: params.count1,
        spacing1: params.spacing1 * scale,
        direction2: params.direction2,
        count2: params.count2,
        spacing2: params.spacing2 * scale,
        stagger: params.stagger,
    }
}

//  € € Internal helpers  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

fn translation_matrix(offset: DVec3) -> DMat4 {
    DMat4::from_translation(offset)
}

fn rotation_matrix(origin: DVec3, axis: DVec3, angle: f64) -> DMat4 {
    DMat4::from_translation(origin)
        * DMat4::from_axis_angle(axis, angle)
        * DMat4::from_translation(-origin)
}

fn append_transformed_brep(
    target: &mut BRep,
    source: &BRep,
    mat: &DMat4,
) -> Result<(), PatternError> {
    use crate::geom::{transform_curve, transform_surface};
    let da = glam::DAffine3::from_mat4(*mat);
    // Build vertex/edge index maps: old index -> new Shape
    let mut v_map: Vec<topods::Shape> = Vec::new();
    let mut e_map: Vec<topods::Shape> = Vec::new();

    // Pass 1: copy vertices
    for ts in &source.tshapes {
        if let topods::TShape::Vertex(vd) = &**ts {
            let new_pt = da.transform_point3(vd.point);
            let sr = target.add_tvertex(new_pt);
            v_map.push(sr);
        }
    }

    // Pass 2: copy edges
    for ts in &source.tshapes {
        if let topods::TShape::Edge(ed) = &**ts {
            let first = v_map.get(ed.first.index).cloned().unwrap_or(topods::Shape::null());
            let last = v_map.get(ed.last.index).cloned().unwrap_or(topods::Shape::null());
            let curve = ed.curve.as_ref().map(|c| transform_curve(c, &da));
            let sr = target.add_tedge(curve, first.clone(), last.clone(), ed.range);
            e_map.push(sr);
        }
    }

    // Pass 3: copy solids with shells/faces/wires
    for ts in &source.tshapes {
        let topods::TShape::Solid(sd) = &**ts else {
            continue;
        };
        let mut shell_refs = Vec::new();
        for sr in &sd.shells {
            if let topods::TShape::Shell(shd) = &*source.tshapes[sr.index] {
                let mut face_refs: Vec<topods::Shape> = Vec::new();
                for fsr in &shd.faces {
                    if let topods::TShape::Face(fd) = &*source.tshapes[fsr.index] {
                        // Outer wire
                        let outer_edges: Vec<topods::Shape> = {
                            if let topods::TShape::Wire(wd) = &*source.tshapes[fd.outer_wire.index]
                            {
                                wd.edges
                                    .iter()
                                    .map(|esr| {
                                        let ne = e_map
                                            .get(esr.index)
                                            .cloned().unwrap_or(topods::Shape::null());
                                        topods::Shape::synthetic(
                                            ne.index,
                                            esr.orientation,
                                        )
                                    })
                                    .collect()
                            } else {
                                Vec::new()
                            }
                        };
                        let outer_wire = target.add_twire(outer_edges);
                        // Inner wires
                        let mut inner_wires = Vec::new();
                        for iwsr in &fd.inner_wires {
                            if let topods::TShape::Wire(iwd) = &*source.tshapes[iwsr.index] {
                                let ies: Vec<topods::Shape> = iwd
                                    .edges
                                    .iter()
                                    .map(|esr| {
                                        let ne = e_map
                                            .get(esr.index)
                                            .cloned().unwrap_or(topods::Shape::null());
                                        topods::Shape::synthetic(
                                            ne.index,
                                            esr.orientation,
                                        )
                                    })
                                    .collect();
                                inner_wires.push(target.add_twire(ies));
                            }
                        }
                        let surface = fd.surface.as_ref().map(|s| transform_surface(s, &da));
                        target.add_tface(
                            surface,
                            outer_wire,
                            inner_wires,
                            fd.sample_point,
                            fd.uv_domain,
                            Vec::new(),
                            true,
                        );
                    }
                }
            }
            // TODO: push face_refs to shell_refs when topods Shell/face tracking is available
        }
        target.add_tsolid(shell_refs);
    }

    Ok(())
}

fn transform_curve(curve: &Curve3, mat: &DMat4) -> Curve3 {
    let transform_point = |p: DVec3| {
        let r = mat.transform_point3(p);
        DVec3::new(r.x, r.y, r.z)
    };
    let transform_direction = |v: DVec3| {
        let r = mat.transform_vector3(v);
        DVec3::new(r.x, r.y, r.z).normalize_or(v)
    };

    match curve {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: transform_point(l.origin),
            direction: transform_direction(l.direction),
        }),
        Curve3::Circle(c) => Curve3::Circle(Circle3::new(
            transform_point(c.center),
            transform_direction(c.normal),
            c.radius,
        )),
        Curve3::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            center: transform_point(e.center),
            normal: transform_direction(e.normal),
            major_dir: transform_direction(e.major_dir),
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        }),
        Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
            center: transform_point(h.center),
            normal: transform_direction(h.normal),
            major_dir: transform_direction(h.major_dir),
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        }),
        Curve3::BSpline(b) => {
            let mut nb = b.clone();
            for cp in &mut nb.control_points {
                *cp = transform_point(*cp);
            }
            Curve3::BSpline(nb)
        }
        Curve3::Bezier(b) => {
            let mut nb = b.clone();
            for cp in &mut nb.control_points {
                *cp = transform_point(*cp);
            }
            Curve3::Bezier(nb)
        }
        _ => curve.clone(),
    }
}

fn transform_surface(surface: &Surface3, mat: &DMat4) -> Surface3 {
    let transform_point = |p: DVec3| {
        let r = mat.transform_point3(p);
        DVec3::new(r.x, r.y, r.z)
    };
    let transform_direction = |v: DVec3| {
        let r = mat.transform_vector3(v);
        DVec3::new(r.x, r.y, r.z).normalize_or(v)
    };

    match surface {
        Surface3::Plane(p) => Surface3::Plane(Plane::new(
            transform_point(p.origin),
            transform_direction(p.normal),
        )),
        Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
            origin: transform_point(c.origin),
            axis: transform_direction(c.axis),
            radius: c.radius,
            ref_dir: transform_direction(c.ref_dir),
        }),
        Surface3::Sphere(s) => Surface3::Sphere(SphericalSurface {
            center: transform_point(s.center),
            axis: transform_direction(s.axis),
            radius: s.radius,
            ref_dir: any_perpendicular(transform_direction(s.axis)),
        }),
        Surface3::Cone(c) => Surface3::Cone(ConicalSurface {
            apex: transform_point(c.apex),
            axis: transform_direction(c.axis),
            radius: c.radius,
            half_angle_rad: c.half_angle_rad,
        }),
        Surface3::Torus(t) => Surface3::Torus(ToroidalSurface {
            center: transform_point(t.center),
            axis: transform_direction(t.axis),
            major_radius: t.major_radius,
            minor_radius: t.minor_radius,
        }),
        Surface3::BSpline(b) => {
            let mut nb = b.clone();
            for row in &mut nb.control_points {
                for cp in row {
                    *cp = transform_point(*cp);
                }
            }
            Surface3::BSpline(nb)
        }
        Surface3::LinearExtrusion(le) => Surface3::LinearExtrusion(LinearExtrusionSurface {
            profile: Box::new(transform_curve(&le.profile, mat)),
            direction: le.direction,
        }),
        Surface3::Revolution(r) => Surface3::Revolution(RevolutionSurface {
            profile: Box::new(transform_curve(&r.profile, mat)),
            axis_origin: transform_point(r.axis_origin),
            axis_dir: transform_direction(r.axis_dir),
        }),
        Surface3::Offset(o) => Surface3::Offset(OffsetSurface {
            basis: Box::new(transform_surface(&o.basis, mat)),
            offset_distance: o.offset_distance,
        }),
        Surface3::Trimmed(t) => Surface3::Trimmed(TrimmedSurface {
            basis: Box::new(transform_surface(&t.basis, mat)),
            trim: t.trim,
        }),
        _ => surface.clone(),
    }
}

//  € € Tests  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
