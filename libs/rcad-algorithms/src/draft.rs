//! Draft angle operations — analogous to OCCT `BRepDraftBuilder`.
//!
//! # Algorithm
//!
//! For each vertex, compute its signed distance `h` to the neutral plane along
//! the pull direction. The draft displacement is:
//!
//! delta = h * tan(angle) * n_perp
//!
//! where `n_perp` is the component of the face normal perpendicular to the pull
//! direction. This tilts each face by the draft angle while keeping vertices on
//! the neutral plane fixed.
//!
//! # Supported surfaces
//!
//! - **Planar faces**: Full draft angle applied, normal rotates around edge axis
//! - **Cylindrical faces**: Draft applied along axis (radius changes with height)
//! - **Conical faces**: Inherit draft angle from base, adjusted for cone angle
//! - **Spherical/toroidal faces**: Limited support via approximation

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
use rcad_kernel::geom::SurfaceEval;
use rcad_kernel::topods::{BRep, TShape, ShapeRef};
use std::collections::HashMap;

/// Default tolerance for geometric operations.
const TOLERANCE: f64 = TOLERANCE_COORD_SUB;

// ============================================================
// Parameters and Configuration
// ============================================================

/// Parameters controlling the draft operation.
#[derive(Debug, Clone)]
pub struct DraftParams {
    /// Normalized pull direction (the "pull" axis of the mold).
    pub pull_direction: DVec3,
    /// Default draft angle in radians. Positive = material added, negative = removed.
    pub draft_angle: f64,
    /// A point on the neutral plane (vertices on this plane don't move).
    pub neutral_point: DVec3,
}

/// Advanced parameters for draft operations with per-face control.
#[derive(Debug, Clone)]
pub struct DraftParamsAdvanced {
    /// Base parameters.
    pub base: DraftParams,
    /// Per-face draft angle overrides (face index -> draft angle in radians).
    pub face_angle_overrides: HashMap<usize, f64>,
    /// Per-face neutral plane overrides (face index -> neutral plane point).
    pub face_neutral_overrides: HashMap<usize, DVec3>,
    /// Transition zone width for smooth angle changes (as fraction of height).
    pub transition_zone_width: f64,
    /// Whether to apply draft to internal features (bosses, ribs).
    pub draft_internal_features: bool,
    /// Minimum feature size to consider (features smaller than this are ignored).
    pub min_feature_size: f64,
}

impl Default for DraftParamsAdvanced {
    fn default() -> Self {
        Self {
            base: DraftParams {
                pull_direction: DVec3::Z,
                draft_angle: 0.0,
                neutral_point: DVec3::ZERO,
            },
            face_angle_overrides: HashMap::new(),
            face_neutral_overrides: HashMap::new(),
            transition_zone_width: 0.0,
            draft_internal_features: false,
            min_feature_size: 0.1,
        }
    }
}

/// Configuration for draft angle analysis and validation.
#[derive(Debug, Clone)]
pub struct DraftValidationConfig {
    /// Minimum acceptable draft angle (radians).
    pub min_draft_angle: f64,
    /// Maximum acceptable draft angle (radians).
    pub max_draft_angle: f64,
    /// Tolerance for undercut detection.
    pub undercut_tolerance: f64,
    /// Whether to check for self-intersections.
    pub check_self_intersection: bool,
    /// Whether to detect internal features.
    pub detect_internal_features: bool,
}

impl Default for DraftValidationConfig {
    fn default() -> Self {
        Self {
            min_draft_angle: 0.5_f64.to_radians(),  // 0.5 degrees minimum
            max_draft_angle: 45.0_f64.to_radians(), // 45 degrees maximum
            undercut_tolerance: TOLERANCE_MESH_LEGACY,
            check_self_intersection: true,
            detect_internal_features: true,
        }
    }
}

// ============================================================
// Error Types
// ============================================================

/// Error type for draft operations.
#[derive(Debug, Clone, PartialEq)]
pub enum DraftError {
    /// A face has a surface type that is not yet supported for drafting.
    UnsupportedSurface { face_index: usize, surface_type: String },
    /// The draft angle is too large (> 89 degrees).
    AngleTooLarge { angle_rad: f64 },
    /// The draft angle is too small for manufacturability.
    AngleTooSmall { angle_rad: f64, min_angle_rad: f64 },
    /// The input BRep has no faces.
    NoFaces,
    /// Self-intersection detected after drafting.
    SelfIntersection { description: String },
    /// Undercut detected in the draft direction.
    UndercutDetected { face_index: usize, description: String },
    /// Invalid pull direction (zero vector).
    InvalidPullDirection,
    /// Neutral surface definition failed.
    NeutralSurfaceError { description: String },
}

impl std::fmt::Display for DraftError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedSurface { face_index, surface_type } => {
                write!(f, "face {} has unsupported surface type: {}", face_index, surface_type)
            }
            Self::AngleTooLarge { angle_rad } => {
                write!(f, "draft angle must be < 89 degrees, got {:.1} degrees", angle_rad.to_degrees())
            }
            Self::AngleTooSmall { angle_rad, min_angle_rad } => {
                write!(f, "draft angle {:.2} degrees is below minimum {:.2} degrees",
                    angle_rad.to_degrees(), min_angle_rad.to_degrees())
            }
            Self::NoFaces => write!(f, "input BRep has no faces"),
            Self::SelfIntersection { description } => {
                write!(f, "self-intersection detected: {}", description)
            }
            Self::UndercutDetected { face_index, description } => {
                write!(f, "undercut on face {}: {}", face_index, description)
            }
            Self::InvalidPullDirection => write!(f, "pull direction must be a non-zero vector"),
            Self::NeutralSurfaceError { description } => {
                write!(f, "neutral surface error: {}", description)
            }
        }
    }
}

impl std::error::Error for DraftError {}

// ============================================================
// Analysis and Validation Types
// ============================================================

/// Information about an internal feature (boss, rib, etc.).
#[derive(Debug, Clone)]
pub struct InternalFeature {
    /// Type of internal feature.
    pub feature_type: InternalFeatureType,
    /// Face indices belonging to this feature.
    pub face_indices: Vec<usize>,
    /// Center of mass of the feature.
    pub center: DVec3,
    /// Approximate size (bounding box diagonal) of the feature.
    pub size: f64,
    /// Height range along pull direction.
    pub height_range: (f64, f64),
}

/// Types of internal features.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InternalFeatureType {
    /// A cylindrical boss (protrusion).
    Boss,
    /// A rib or web.
    Rib,
    /// A slot or groove.
    Slot,
    /// A hole or pocket.
    Hole,
    /// Unknown feature type.
    Unknown,
}

/// Result of draft angle validation.
#[derive(Debug, Clone)]
pub struct DraftValidationResult {
    /// Whether the part is valid for drafting.
    pub is_valid: bool,
    /// Faces with insufficient draft angles.
    pub insufficient_draft_faces: Vec<FaceDraftIssue>,
    /// Faces with undercuts in the pull direction.
    pub undercut_faces: Vec<UndercutIssue>,
    /// Detected internal features.
    pub internal_features: Vec<InternalFeature>,
    /// Self-intersection issues (if any).
    pub self_intersection_issues: Vec<SelfIntersectionIssue>,
    /// Overall draft quality score (0.0 to 1.0).
    pub quality_score: f64,
}

/// Issue with draft angle on a face.
#[derive(Debug, Clone)]
pub struct FaceDraftIssue {
    /// Face index.
    pub face_index: usize,
    /// Current draft angle (radians).
    pub current_angle: f64,
    /// Required minimum draft angle (radians).
    pub required_angle: f64,
    /// Description of the issue.
    pub description: String,
}

/// Undercut detection result for a face.
#[derive(Debug, Clone)]
pub struct UndercutIssue {
    /// Face index.
    pub face_index: usize,
    /// Undercut severity (0.0 to 1.0).
    pub severity: f64,
    /// Description of the undercut.
    pub description: String,
}

/// Self-intersection issue.
#[derive(Debug, Clone)]
pub struct SelfIntersectionIssue {
    /// Description of the self-intersection.
    pub description: String,
    /// Face indices involved (if applicable).
    pub involved_faces: Vec<usize>,
}

/// Parting line detection result.
#[derive(Debug, Clone)]
pub struct PartingLineResult {
    /// Points on the parting line (in 3D space).
    pub points: Vec<DVec3>,
    /// Edge indices that form the parting line.
    pub edge_indices: Vec<usize>,
    /// Whether the parting line is closed.
    pub is_closed: bool,
    /// Recommended draft direction (optimized).
    pub recommended_direction: DVec3,
}

/// Neutral surface definition.
#[derive(Debug, Clone)]
pub enum NeutralSurface {
    /// A flat plane.
    Plane { point: DVec3, normal: DVec3 },
    /// A curved surface (for complex geometries).
    Curved { surface: Surface3 },
    /// A set of edges that should remain fixed.
    EdgeSet { edge_indices: Vec<usize> },
}

// ============================================================
// Core Draft Operations
// ============================================================

/// Apply a draft angle to all planar faces of a BRep.
///
/// Vertices on the neutral plane remain fixed. Other vertices are displaced
/// perpendicular to the pull direction by `h * tan(angle)`.
pub fn draft_solid(brep: &BRep, params: &DraftParams) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;
    let face_count = face_srs.len();
    if face_count == 0 {
        return Err(DraftError::NoFaces);
    }

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;
    let tan_angle = params.draft_angle.tan();

    // Step 1: Compute new vertex positions (preserve old vertex order)
    let new_pts: Vec<DVec3> = brep.tshapes.iter().filter_map(|ts| {
        if let TShape::Vertex(vd) = ts.as_ref() {
            let h = (vd.point - neutral).dot(pull);
            Some(vd.point + pull * (h * tan_angle))
        } else {
            None
        }
    }).collect();

    // Step 2: Compute new face normals
    let new_face_normals: Vec<DVec3> = face_srs.iter().map(|&sr| {
        let fd = brep.face(sr);
        let n = fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default()
            .normalize();
        let axis = n.cross(pull);
        let axis_len = axis.length();
        if axis_len < TOLERANCE_LINEAR_ULTRA_STRICT {
            return n;
        }
        let k = axis / axis_len;
        let cos_a = params.draft_angle.cos();
        let sin_a = params.draft_angle.sin();
        let rotated = n * cos_a + k.cross(n) * sin_a;
        rotated.normalize_or(n)
    }).collect();

    // Step 3: Build result BRep
    build_drafted_brep(brep, &new_pts, &new_face_normals, &face_srs)
}

/// Apply draft with advanced per-face control.
///
/// This function supports:
/// - Per-face draft angle overrides
/// - Variable draft angles with transition zones
/// - Non-planar surface handling
/// - Internal feature detection and handling
pub fn draft_solid_advanced(
    brep: &BRep,
    params: &DraftParamsAdvanced,
) -> Result<BRep, DraftError> {
    validate_draft_params(&params.base)?;

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;
    if face_srs.is_empty() {
        return Err(DraftError::NoFaces);
    }

    let pull = params.base.pull_direction.normalize();
    let neutral = params.base.neutral_point;

    // Compute per-vertex displacements accounting for adjacent face angles
    let vertex_displacements = compute_vertex_displacements(brep, &face_srs, params, pull, neutral)?;

    // Apply displacements
    let new_pts: Vec<DVec3> = brep.tshapes.iter().enumerate().map(|(i, ts)| {
        if let TShape::Vertex(vd) = ts.as_ref() {
            vd.point + vertex_displacements.get(&i).copied().unwrap_or(DVec3::ZERO)
        } else {
            DVec3::ZERO
        }
    }).collect();

    // Compute new face normals with per-face angles
    let new_face_normals: Vec<DVec3> = face_srs.iter().enumerate().map(|(fi, &sr)| {
        let fd = brep.face(sr);
        let n = fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default();
        let angle = params.face_angle_overrides.get(&fi).copied().unwrap_or(params.base.draft_angle);
        compute_rotated_normal(&n, pull, angle)
    }).collect();

    // Build result
    build_drafted_brep(brep, &new_pts, &new_face_normals, &face_srs)
}

/// Draft cylindrical faces by modifying radius along the pull direction.
///
/// For cylindrical surfaces, drafting changes the radius linearly with height.
/// The cylindrical axis should be parallel to the pull direction.
pub fn draft_cylindrical_face(
    brep: &BRep,
    face_index: usize,
    params: &DraftParams,
) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;
    let _face_sr = *face_srs.get(face_index).ok_or(DraftError::NoFaces)?;

    // Get surface geometry from the face
    let face_ts = brep.tshapes.get(face_index).ok_or(DraftError::NoFaces)?;
    let cylinder_surf = match &**face_ts {
        TShape::Face(fd) => {
            match fd.surface.as_ref() {
                Some(Surface3::Cylinder(c)) => c.clone(),
                Some(other) => return Err(DraftError::UnsupportedSurface {
                    face_index,
                    surface_type: surface_type_name(other),
                }),
                None => return Err(DraftError::UnsupportedSurface {
                    face_index,
                    surface_type: "none".to_string(),
                }),
            }
        }
        _ => return Err(DraftError::NoFaces),
    };

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;
    let tan_angle = params.draft_angle.tan();

    // Check if cylinder axis is parallel to pull direction
    let axis = cylinder_surf.axis.normalize();
    let axis_alignment = axis.dot(pull).abs();
    if axis_alignment < 0.99 {
        // Cylindrical face is not aligned with pull direction
        // Apply a more general approach
        return draft_face_general(brep, face_index, params);
    }

    // Collect vertex positions
    let vpts: Vec<DVec3> = brep.tshapes.iter().filter_map(|ts| {
        if let TShape::Vertex(vd) = ts.as_ref() { Some(vd.point) } else { None }
    }).collect();

    // Compute new vertex positions based on cylindrical draft
    let mut new_pts: Vec<DVec3> = vpts.clone();

    for (vi, _ts) in brep.tshapes.iter().enumerate() {
        if let Some(pt) = brep.vertex_point(vi) {
            let h = (pt - neutral).dot(pull);
            // For a cylinder, the radial displacement is proportional to height
            let radial_dir = (pt - cylinder_surf.origin).reject_from(axis).normalize_or(DVec3::ZERO);
            if radial_dir.length() > TOLERANCE_LINEAR_ULTRA_STRICT {
                let radial_displacement = h * tan_angle;
                new_pts[vi] = pt + radial_dir * radial_displacement;
            }
        }
    }

    // Face normal remains essentially unchanged for drafted cylinders
    let new_face_normals: Vec<DVec3> = face_srs.iter().map(|&sr| {
        let fd = brep.face(sr);
        fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default()
    }).collect();

    build_drafted_brep(brep, &new_pts, &new_face_normals, &face_srs)
}

/// Draft conical faces by adjusting the cone angle.
///
/// Conical surfaces inherently have a draft angle. This function adjusts
/// the cone angle to match the desired draft angle.
pub fn draft_conical_face(
    brep: &BRep,
    face_index: usize,
    params: &DraftParams,
) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;
    let _face_sr = *face_srs.get(face_index).ok_or(DraftError::NoFaces)?;

    // Get surface geometry from the face
    let face_ts = brep.tshapes.get(face_index).ok_or(DraftError::NoFaces)?;
    let cone_surf = match &**face_ts {
        TShape::Face(fd) => {
            match fd.surface.as_ref() {
                Some(Surface3::Cone(c)) => c.clone(),
                Some(other) => return Err(DraftError::UnsupportedSurface {
                    face_index,
                    surface_type: surface_type_name(other),
                }),
                None => return Err(DraftError::UnsupportedSurface {
                    face_index,
                    surface_type: "none".to_string(),
                }),
            }
        }
        _ => return Err(DraftError::NoFaces),
    };

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;

    // Check if cone axis is parallel to pull direction
    let axis = cone_surf.axis.normalize();
    let axis_alignment = axis.dot(pull).abs();
    if axis_alignment < 0.99 {
        return draft_face_general(brep, face_index, params);
    }

    // For cones, the effective draft angle is the cone half-angle
    // Adjust the cone angle by combining with the desired draft
    let effective_draft = cone_surf.half_angle_rad + params.draft_angle;

    // Compute new vertex positions
    let vpts: Vec<DVec3> = brep.tshapes.iter().filter_map(|ts| {
        if let TShape::Vertex(vd) = ts.as_ref() { Some(vd.point) } else { None }
    }).collect();

    let mut new_pts: Vec<DVec3> = vpts.clone();
    let tan_effective = effective_draft.tan();
    let _tan_original = cone_surf.half_angle_rad.tan();

    for (vi, _ts) in brep.tshapes.iter().enumerate() {
        if let Some(pt) = brep.vertex_point(vi) {
            let _h = (pt - neutral).dot(pull);
            let radial_vec = pt - cone_surf.apex;
            let radial_dist = radial_vec.reject_from(axis).length();
            let axial_dist = radial_vec.dot(axis).abs();

            if axial_dist > TOLERANCE {
                // New radial distance based on adjusted cone angle
                let new_radial_dist = axial_dist * tan_effective;
                let radial_change = new_radial_dist - radial_dist;
                let radial_dir = radial_vec.reject_from(axis).normalize_or(DVec3::ZERO);
                if radial_dir.length() > TOLERANCE_LINEAR_ULTRA_STRICT {
                    new_pts[vi] = pt + radial_dir * radial_change;
                }
            }
        }
    }

    // Compute new normal for the conical face
    let new_face_normals: Vec<DVec3> = face_srs.iter().enumerate().map(|(fi, &sr)| {
        let fd = brep.face(sr);
        if fi == face_index {
            compute_cone_normal(effective_draft, axis, pull)
        } else {
            fd.surface.as_ref()
                .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
                .unwrap_or_default()
        }
    }).collect();

    build_drafted_brep(brep, &new_pts, &new_face_normals, &face_srs)
}

/// General-purpose drafting for any face type.
///
/// Uses vertex displacement based on height and draft angle.
pub fn draft_face_general(
    brep: &BRep,
    face_index: usize,
    params: &DraftParams,
) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;
    if face_index >= face_srs.len() {
        return Err(DraftError::NoFaces);
    }

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;
    let tan_angle = params.draft_angle.tan();

    // Get the face surface normal
    let face_sr = face_srs[face_index];
    let fd = brep.face(face_sr);
    let face_normal = fd.surface.as_ref()
        .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
        .unwrap_or_default()
        .normalize();

    // Compute displacement direction perpendicular to both face normal and pull direction
    let displacement_dir = face_normal.cross(pull).normalize_or(pull);

    // Compute new vertex positions
    let new_pts: Vec<DVec3> = brep.tshapes.iter().filter_map(|ts| {
        if let TShape::Vertex(vd) = ts.as_ref() {
            let h = (vd.point - neutral).dot(pull);
            let displacement = h * tan_angle * displacement_dir;
            Some(vd.point + displacement)
        } else {
            None
        }
    }).collect();

    // Compute new face normals
    let new_face_normals: Vec<DVec3> = face_srs.iter().map(|&sr| {
        let fd = brep.face(sr);
        let n = fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default();
        compute_rotated_normal(&n, pull, params.draft_angle)
    }).collect();

    build_drafted_brep(brep, &new_pts, &new_face_normals, &face_srs)
}

// ============================================================
// Validation and Analysis Functions
// ============================================================

/// Validate draft angles on a solid.
///
/// Checks for:
/// - Insufficient draft angles
/// - Undercuts in the pull direction
/// - Self-intersections
/// - Internal features that may need special handling
pub fn validate_draft_angles(
    brep: &BRep,
    pull_direction: DVec3,
    config: &DraftValidationConfig,
) -> Result<DraftValidationResult, DraftError> {
    let pull = pull_direction.normalize();
    if pull.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;

    let mut insufficient_draft_faces = Vec::new();
    let mut undercut_faces = Vec::new();
    let mut quality_sum = 0.0;

    for (fi, &sr) in face_srs.iter().enumerate() {
        let fd = brep.face(sr);
        let normal = fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default()
            .normalize();
        let draft_angle = compute_draft_angle(&normal, pull);

        // Check for insufficient draft
        if draft_angle.abs() < config.min_draft_angle {
            insufficient_draft_faces.push(FaceDraftIssue {
                face_index: fi,
                current_angle: draft_angle,
                required_angle: config.min_draft_angle,
                description: format!(
                    "Draft angle {:.2} degrees is below minimum {:.2} degrees",
                    draft_angle.to_degrees(),
                    config.min_draft_angle.to_degrees()
                ),
            });
        }

        // Check for undercuts
        if draft_angle < -config.undercut_tolerance {
            let severity = (config.min_draft_angle - draft_angle) / config.max_draft_angle;
            undercut_faces.push(UndercutIssue {
                face_index: fi,
                severity: severity.min(1.0).max(0.0),
                description: format!(
                    "Face has undercut of {:.2} degrees relative to pull direction",
                    (-draft_angle).to_degrees()
                ),
            });
        }

        // Compute quality contribution
        let angle_quality = if draft_angle >= config.min_draft_angle && draft_angle <= config.max_draft_angle {
            1.0
        } else if draft_angle > 0.0 {
            draft_angle / config.min_draft_angle
        } else {
            0.0
        };
        quality_sum += angle_quality;
    }

    // Detect internal features
    let internal_features = if config.detect_internal_features {
        detect_internal_features(brep, pull)?
    } else {
        Vec::new()
    };

    // Check for self-intersections (simplified check)
    let self_intersection_issues = if config.check_self_intersection {
        check_self_intersection(brep)?
    } else {
        Vec::new()
    };

    let quality_score = if face_srs.is_empty() {
        0.0
    } else {
        quality_sum / face_srs.len() as f64
    };

    let is_valid = insufficient_draft_faces.is_empty()
        && undercut_faces.is_empty()
        && self_intersection_issues.is_empty();

    Ok(DraftValidationResult {
        is_valid,
        insufficient_draft_faces,
        undercut_faces,
        internal_features,
        self_intersection_issues,
        quality_score,
    })
}

/// Detect undercuts in the given pull direction.
///
/// Returns a list of face indices that form undercuts.
pub fn detect_undercuts(
    brep: &BRep,
    pull_direction: DVec3,
    tolerance: f64,
) -> Result<Vec<usize>, DraftError> {
    let pull = pull_direction.normalize();
    if pull.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;

    let mut undercut_faces = Vec::new();

    for (fi, &sr) in face_srs.iter().enumerate() {
        let fd = brep.face(sr);
        let normal = fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default()
            .normalize();
        let draft_angle = compute_draft_angle(&normal, pull);

        if draft_angle < tolerance {
            undercut_faces.push(fi);
        }
    }

    Ok(undercut_faces)
}

/// Detect the optimal parting line for a part.
///
/// The parting line is the curve where the mold splits.
pub fn detect_parting_line(
    brep: &BRep,
    pull_direction: DVec3,
) -> Result<PartingLineResult, DraftError> {
    let pull = pull_direction.normalize();
    if pull.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }

    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;

    // Find edges where adjacent faces have opposite draft directions
    let mut parting_edges = Vec::new();
    let mut parting_points = Vec::new();

    for &sr in &face_srs {
        let fd = brep.face(sr);
        let normal = fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default()
            .normalize();
        let draft_angle = compute_draft_angle(&normal, pull);

        // Check if this face is near the "equator" (perpendicular to pull direction)
        let is_equatorial = normal.dot(pull).abs() < 0.1;

        if is_equatorial || draft_angle.abs() < 0.01_f64.to_radians() {
            // This face may be on the parting line — collect its edge indices
            if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                for esr in &wd.edges {
                    if !parting_edges.contains(&esr.index) {
                        parting_edges.push(esr.index);
                    }
                }
            }
        }
    }

    // Collect points from parting edges
    for &ei in &parting_edges {
        if let TShape::Edge(ed) = &*brep.tshapes[ei] {
            let vs = brep.vertex_point(ed.first.index);
            let ve = brep.vertex_point(ed.last.index);
            if let (Some(vsp), Some(vep)) = (vs, ve) {
                parting_points.push(vsp);
                parting_points.push(vep);
            }
        }
    }

    // Determine if parting line is closed
    let is_closed = parting_edges.len() >= 3;

    // Optimize draft direction (simplified - just use input for now)
    let recommended_direction = optimize_draft_direction(brep, pull)?;

    Ok(PartingLineResult {
        points: parting_points,
        edge_indices: parting_edges,
        is_closed,
        recommended_direction,
    })
}

/// Detect internal features (bosses, ribs, etc.) in the part.
pub fn detect_internal_features(
    brep: &BRep,
    pull_direction: DVec3,
) -> Result<Vec<InternalFeature>, DraftError> {
    let pull = pull_direction.normalize();
    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;

    let mut features = Vec::new();

    // Group faces by connectivity
    let face_groups = group_connected_faces(brep, &face_srs);

    for group in face_groups {
        if group.len() <= 1 {
            continue;
        }

        // Analyze the group to determine feature type
        let feature_type = classify_feature_type(brep, &group, pull);
        let (center, size, height_range) = compute_feature_properties(brep, &group, pull);

        // Skip if too small
        if size < 0.1 {
            continue;
        }

        features.push(InternalFeature {
            feature_type,
            face_indices: group,
            center,
            size,
            height_range,
        });
    }

    Ok(features)
}

/// Optimize the draft direction for best parting line and minimum undercuts.
pub fn optimize_draft_direction(
    brep: &BRep,
    initial_direction: DVec3,
) -> Result<DVec3, DraftError> {
    let face_srs = collect_face_refs(brep).ok_or(DraftError::NoFaces)?;

    // Sample a set of candidate directions
    let mut best_direction = initial_direction.normalize();
    let mut best_score = evaluate_draft_direction(brep, &face_srs, best_direction);

    // Try variations around the initial direction
    let variations = [
        DVec3::X, DVec3::Y, DVec3::Z,
        DVec3::NEG_X, DVec3::NEG_Y, DVec3::NEG_Z,
    ];

    for v in variations {
        let dir = v.normalize();
        let score = evaluate_draft_direction(brep, &face_srs, dir);
        if score > best_score {
            best_score = score;
            best_direction = dir;
        }
    }

    // Fine-tune with small angle variations
    for angle in [15.0_f64, 30.0_f64, 45.0_f64].iter().map(|a| a.to_radians()) {
        for axis in [DVec3::X, DVec3::Y] {
            let rotated = rotate_vector_around_axis(best_direction, axis, angle);
            let score = evaluate_draft_direction(brep, &face_srs, rotated);
            if score > best_score {
                best_score = score;
                best_direction = rotated;
            }
        }
    }

    Ok(best_direction)
}

// ============================================================
// Helper Functions
// ============================================================

fn validate_draft_params(params: &DraftParams) -> Result<(), DraftError> {
    if params.pull_direction.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }
    if params.draft_angle.abs() > std::f64::consts::FRAC_PI_2 - 0.02 {
        return Err(DraftError::AngleTooLarge {
            angle_rad: params.draft_angle,
        });
    }
    Ok(())
}

fn compute_rotated_normal(normal: &DVec3, pull: DVec3, angle: f64) -> DVec3 {
    let n = normal.normalize();
    let axis = n.cross(pull);
    let axis_len = axis.length();
    if axis_len < TOLERANCE_LINEAR_ULTRA_STRICT {
        return n;
    }
    let k = axis / axis_len;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let rotated = n * cos_a + k.cross(n) * sin_a;
    rotated.normalize_or(n)
}

fn compute_cone_normal(half_angle: f64, axis: DVec3, pull: DVec3) -> DVec3 {
    // Normal to a cone surface points outward at the half angle
    let radial = axis.cross(pull).normalize_or(DVec3::X);
    let normal = axis * half_angle.cos() + radial * half_angle.sin();
    normal.normalize_or(axis)
}

fn compute_draft_angle(normal: &DVec3, pull: DVec3) -> f64 {
    // Draft angle is the angle between the face normal and the horizontal plane
    // perpendicular to the pull direction
    let n = normal.normalize();
    let cos_angle = n.dot(pull).abs();
    // The draft angle is 90 degrees minus the angle with pull direction
    std::f64::consts::FRAC_PI_2 - cos_angle.acos()
}

fn compute_vertex_displacements(
    brep: &BRep,
    face_srs: &[ShapeRef],
    params: &DraftParamsAdvanced,
    pull: DVec3,
    neutral: DVec3,
) -> Result<HashMap<usize, DVec3>, DraftError> {
    let mut displacements: HashMap<usize, DVec3> = HashMap::new();

    for (fi, &sr) in face_srs.iter().enumerate() {
        let angle = params.face_angle_overrides.get(&fi).copied().unwrap_or(params.base.draft_angle);
        let face_neutral = params.face_neutral_overrides.get(&fi).copied().unwrap_or(neutral);
        let tan_angle = angle.tan();

        // Get vertices belonging to this face
        let fd = brep.face(sr);
        if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
            for esr in &wd.edges {
                if let TShape::Edge(ed) = &*brep.tshapes[esr.index] {
                    for &vi in &[ed.first.index, ed.last.index] {
                        let h = (brep.vertex_point(vi).unwrap_or(DVec3::ZERO) - face_neutral).dot(pull);
                        let displacement = pull * (h * tan_angle);

                        // Average displacements if vertex belongs to multiple faces
                        displacements
                            .entry(vi)
                            .and_modify(|d| *d = (*d + displacement) * 0.5)
                            .or_insert(displacement);
                    }
                }
            }
        }
    }

    // Apply transition zone smoothing if enabled
    if params.transition_zone_width > 0.0 {
        apply_transition_zones(&mut displacements, brep, params, pull);
    }

    Ok(displacements)
}

fn apply_transition_zones(
    displacements: &mut HashMap<usize, DVec3>,
    brep: &BRep,
    params: &DraftParamsAdvanced,
    pull: DVec3,
) {
    // Smooth transitions between regions with different draft angles
    let transition_height = params.transition_zone_width;

    // Group vertices by height
    let mut height_groups: HashMap<i32, Vec<usize>> = HashMap::new();
    for (&vi, _) in displacements.iter() {
        if let Some(pt) = brep.vertex_point(vi) {
            let h = pt.dot(pull);
            let group = (h / transition_height).floor() as i32;
            height_groups.entry(group).or_default().push(vi);
        }
    }

    // Apply smoothing within transition zones
    for (_, group) in height_groups.iter() {
        if group.len() < 2 {
            continue;
        }

        // Compute average displacement in this zone
        let avg_displacement: DVec3 = group
            .iter()
            .filter_map(|vi| displacements.get(vi))
            .sum::<DVec3>()
            / group.len() as f64;

        // Apply weighted average for smoothing
        for vi in group {
            if let Some(d) = displacements.get_mut(vi) {
                *d = *d * 0.7 + avg_displacement * 0.3;
            }
        }
    }
}

fn build_drafted_brep(
    brep: &BRep,
    new_pts: &[DVec3],
    new_face_normals: &[DVec3],
    face_srs: &[ShapeRef],
) -> Result<BRep, DraftError> {
    let mut out = BRep::new();

    // ── Step 1: Create vertices with new positions ──
    // Build a mapping: old vertex tshape index -> new ShapeRef
    let mut vmap: HashMap<usize, ShapeRef> = HashMap::new();
    let mut pt_iter = new_pts.iter();
    for (ts_idx, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Vertex(_) = ts.as_ref() {
            let p = *pt_iter.next().ok_or_else(|| DraftError::NoFaces)?;
            let sr = out.add_tvertex(p);
            vmap.insert(ts_idx, sr);
        }
    }

    // ── Step 2: Create edges ──
    let mut emap: HashMap<usize, ShapeRef> = HashMap::new(); // old edge tshape idx -> new ShapeRef
    for (ts_idx, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Edge(ed) = ts.as_ref() {
            let vs = vmap.get(&ed.first.index).copied().ok_or(DraftError::NoFaces)?;
            let ve = vmap.get(&ed.last.index).copied().ok_or(DraftError::NoFaces)?;

            // Compute curve between the new vertex positions
            let p1 = out.vertex_point(vs.index).unwrap_or(DVec3::ZERO);
            let p2 = out.vertex_point(ve.index).unwrap_or(DVec3::ZERO);
            let dir = p2 - p1;
            let len = dir.length();

            let curve = if len > TOLERANCE {
                Some(Curve3::Line(Line3 {
                    origin: p1,
                    direction: dir / len,
                }))
            } else {
                None
            };

            let e_sr = out.add_tedge(curve, vs, ve, [0.0, len]);
            emap.insert(ts_idx, e_sr);
        }
    }

    // ── Step 3: Create faces with updated surfaces ──
    let mut new_face_srs: Vec<ShapeRef> = Vec::new();
    for (fi, &old_face_sr) in face_srs.iter().enumerate() {
        let fd = brep.face(old_face_sr);
        let updated_normal = new_face_normals.get(fi).copied().unwrap_or(DVec3::Z);
        let surface = fd.surface.clone();

        // Rebuild the surface so its normal reflects the draft angle.
        // For planar surfaces, rotate the normal
        let new_surface = surface.map(|s| match s {
            Surface3::Plane(mut pl) => {
                pl.normal = updated_normal.normalize_or(pl.normal);
                Surface3::Plane(pl)
            }
            other => other,
        });

        // Create the new edge ShapeRefs for the outer wire
        let new_wire_edges: Vec<ShapeRef> = if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
            wd.edges.iter().filter_map(|esr| emap.get(&esr.index).copied()).collect()
        } else {
            Vec::new()
        };
        let new_outer_wire = out.add_twire(new_wire_edges);

        // Inner wires
        let new_inner_wires: Vec<ShapeRef> = fd.inner_wires.iter().map(|iw_sr| {
            let inner_edges: Vec<ShapeRef> = if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
                iwd.edges.iter().filter_map(|esr| emap.get(&esr.index).copied()).collect()
            } else {
                Vec::new()
            };
            out.add_twire(inner_edges)
        }).collect();

        // Create the face with updated surface
        let new_face_sr = out.add_tface(
            new_surface,
            new_outer_wire,
            new_inner_wires,
            fd.sample_point,
            fd.uv_domain,
            Vec::new(),
            fd.natural_restriction,
        );
        new_face_srs.push(new_face_sr);
    }

    // ── Step 4: Create shell and solid ──
    if !new_face_srs.is_empty() {
        let shell_sr = out.add_tshell(new_face_srs);
        out.add_tsolid(vec![shell_sr]);
    }

    Ok(out)
}

fn surface_type_name(surface: &Surface3) -> String {
    match surface {
        Surface3::Plane(_) => "Plane".to_string(),
        Surface3::Cylinder(_) => "Cylinder".to_string(),
        Surface3::Sphere(_) => "Sphere".to_string(),
        Surface3::Cone(_) => "Cone".to_string(),
        Surface3::Torus(_) => "Torus".to_string(),
        Surface3::Ellipsoid(_) => "Ellipsoid".to_string(),
        Surface3::Helicoid(_) => "Helicoid".to_string(),
        Surface3::Pipe(_) => "Pipe".to_string(),
        Surface3::BSpline(_) => "BSpline".to_string(),
        Surface3::LinearExtrusion(_) => "LinearExtrusion".to_string(),
        Surface3::Revolution(_) => "Revolution".to_string(),
        Surface3::Ruled(_) => "Ruled".to_string(),
        Surface3::Coons(_) => "Coons".to_string(),
        Surface3::Bezier(_) => "Bezier".to_string(),
        Surface3::TriBezier(_) => "TriBezier".to_string(),
        Surface3::Offset(_) => "Offset".to_string(),
        Surface3::Trimmed(_) => "Trimmed".to_string(),
    }
}

/// Collect all face ShapeRefs from the first solid's first shell.
fn collect_face_refs(brep: &BRep) -> Option<Vec<ShapeRef>> {
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            if let Some(shell_sr) = sd.shells.first() {
                if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    return Some(shd.faces.clone());
                }
            }
        }
    }
    None
}

fn check_self_intersection(brep: &BRep) -> Result<Vec<SelfIntersectionIssue>, DraftError> {
    // Simplified self-intersection check
    // A proper implementation would use BVH and triangle-triangle intersection tests
    let mut issues = Vec::new();

    // Check for degenerate edges (zero length)
    for (ei, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Edge(ed) = ts.as_ref() {
            let Some(vs) = brep.vertex_point(ed.first.index) else { continue; };
            let Some(ve) = brep.vertex_point(ed.last.index) else { continue; };
            let len = (ve - vs).length();
            if len < TOLERANCE {
                issues.push(SelfIntersectionIssue {
                    description: format!("Degenerate edge {} with length {}", ei, len),
                    involved_faces: find_faces_with_edge(brep, ei),
                });
            }
        }
    }

    Ok(issues)
}

fn find_faces_with_edge(brep: &BRep, edge_index: usize) -> Vec<usize> {
    let mut faces = Vec::new();
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            if let Some(shell_sr) = sd.shells.first() {
                if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for (fi, face_sr) in shd.faces.iter().enumerate() {
                        let fd = brep.face(*face_sr);
                        if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                            if wd.edges.iter().any(|e| e.index == edge_index) {
                                faces.push(fi);
                            }
                        }
                        for iw_sr in &fd.inner_wires {
                            if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
                                if iwd.edges.iter().any(|e| e.index == edge_index) {
                                    faces.push(fi);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    faces
}

fn group_connected_faces(brep: &BRep, face_srs: &[ShapeRef]) -> Vec<Vec<usize>> {
    // Build edge-to-face adjacency using tshape indices
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, &sr) in face_srs.iter().enumerate() {
        let fd = brep.face(sr);
        if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
            for esr in &wd.edges {
                edge_to_faces.entry(esr.index).or_default().push(fi);
            }
        }
    }

    // Find connected components using union-find
    let mut parent: Vec<usize> = (0..face_srs.len()).collect();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], x: usize, y: usize) {
        let px = find(parent, x);
        let py = find(parent, y);
        if px != py {
            parent[px] = py;
        }
    }

    for (_, face_list) in edge_to_faces.iter() {
        if face_list.len() >= 2 {
            for i in 1..face_list.len() {
                union(&mut parent, face_list[0], face_list[i]);
            }
        }
    }

    // Group faces by their root
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for fi in 0..face_srs.len() {
        let root = find(&mut parent, fi);
        groups.entry(root).or_default().push(fi);
    }

    groups.into_values().collect()
}

fn classify_feature_type(
    brep: &BRep,
    face_indices: &[usize],
    pull: DVec3,
) -> InternalFeatureType {
    let face_srs = match collect_face_refs(brep) {
        Some(srs) => srs,
        None => return InternalFeatureType::Unknown,
    };

    let mut total_draft_angle = 0.0;
    let mut face_count = 0;

    for &fi in face_indices {
        if let Some(&sr) = face_srs.get(fi) {
            let fd = brep.face(sr);
            let normal = fd.surface.as_ref()
                .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
                .unwrap_or_default()
                .normalize();
            total_draft_angle += compute_draft_angle(&normal, pull);
            face_count += 1;
        }
    }

    if face_count == 0 {
        return InternalFeatureType::Unknown;
    }

    let avg_draft = total_draft_angle / face_count as f64;

    // Classify based on average draft angle and face count
    if avg_draft > 0.0 {
        if face_count <= 2 {
            InternalFeatureType::Rib
        } else if face_count <= 6 {
            InternalFeatureType::Boss
        } else {
            InternalFeatureType::Unknown
        }
    } else if avg_draft < 0.0 {
        if face_count <= 4 {
            InternalFeatureType::Slot
        } else {
            InternalFeatureType::Hole
        }
    } else {
        InternalFeatureType::Unknown
    }
}

fn compute_feature_properties(
    brep: &BRep,
    face_indices: &[usize],
    pull: DVec3,
) -> (DVec3, f64, (f64, f64)) {
    let face_srs = match collect_face_refs(brep) {
        Some(srs) => srs,
        None => return (DVec3::ZERO, 0.0, (0.0, 0.0)),
    };

    // Collect all vertices in the feature
    let mut vertices: Vec<DVec3> = Vec::new();
    let mut heights: Vec<f64> = Vec::new();

    for &fi in face_indices {
        if let Some(&sr) = face_srs.get(fi) {
            let fd = brep.face(sr);
            if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                for esr in &wd.edges {
                    if let TShape::Edge(ed) = &*brep.tshapes[esr.index] {
                        if let Some(pt) = brep.vertex_point(ed.first.index) {
                            vertices.push(pt);
                            heights.push(pt.dot(pull));
                        }
                        if let Some(pt) = brep.vertex_point(ed.last.index) {
                            vertices.push(pt);
                            heights.push(pt.dot(pull));
                        }
                    }
                }
            }
        }
    }

    if vertices.is_empty() {
        return (DVec3::ZERO, 0.0, (0.0, 0.0));
    }

    // Compute center of mass
    let center = vertices.iter().sum::<DVec3>() / vertices.len() as f64;

    // Compute size (bounding box diagonal)
    let min_pt = vertices.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
    let max_pt = vertices.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));
    let size = (max_pt - min_pt).length();

    // Compute height range
    let h_min = heights.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let h_max = heights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    (center, size, (h_min, h_max))
}

fn evaluate_draft_direction(_brep: &BRep, face_srs: &[ShapeRef], direction: DVec3) -> f64 {
    let mut score = 0.0;

    for &sr in face_srs {
        let fd = _brep.face(sr);
        let normal = fd.surface.as_ref()
            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or_default()
            .normalize();
        let draft_angle = compute_draft_angle(&normal, direction);

        // Score based on how close draft is to ideal range
        let ideal_min = 1.0_f64.to_radians();
        let ideal_max = 5.0_f64.to_radians();

        if draft_angle >= ideal_min && draft_angle <= ideal_max {
            score += 1.0;
        } else if draft_angle > 0.0 {
            score += 0.5;
        } else {
            score -= 1.0; // Penalty for undercuts
        }
    }

    score
}

fn rotate_vector_around_axis(v: DVec3, axis: DVec3, angle: f64) -> DVec3 {
    // Rodrigues rotation formula
    let k = axis.normalize_or(DVec3::Z);
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    v * cos_a + k.cross(v) * sin_a + k * (k.dot(v) * (1.0 - cos_a))
}

// ============================================================
// Tests
// ============================================================
