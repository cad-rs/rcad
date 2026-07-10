//! BRepFeat-style feature-based modeling operations.
//!
//! This module provides feature-based modeling operations analogous to OCCT's TKFeat
//! (BRepFeat package). Features are operations that add or remove material from a
//! base shape while maintaining design intent.
//!
//! # Feature Types
//!
//! - **Rib**: Add reinforcing rib features from a wire profile
//! - **Groove**: Create slot/groove features by removing material
//! - **Prism**: Feature-based prism (extrusion) with fuse modes
//! - **Revol**: Feature-based revolution with fuse modes
//! - **Pipe**: Feature-based pipe along a spine curve
//! - **Draft**: Apply draft angle to faces for moldability


use std::sync::Arc;
use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::topods::{self, ShapeRef, TShape, Orientation};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};

use crate::{BooleanError, BooleanOpType};
use crate::bop_occt_union::boolean_op_generic as boolean_op;

// ===========================================================?
// Error Types
// ===========================================================?

/// Errors returned by BRepFeat operations.
#[derive(Debug)]
pub enum BRepFeatError {
 /// Input contains non-finite values.
 NonFiniteInput(&'static str),
 /// Input value must be positive.
 NonPositiveInput(&'static str),
 /// Invalid input geometry or parameters.
 InvalidInput(String),
 /// Zero-length vector where non-zero is required.
 ZeroVector(&'static str),
 /// Vectors are parallel when they should not be.
 ParallelVectors(&'static str, &'static str),
 /// Boolean operation failed.
 BooleanFailed(BooleanError),
 /// Modeling operation failed.
 ModelingFailed(String),
 /// Profile wire is invalid (not closed, too few edges, etc.).
 InvalidProfile(String),
 /// Feature does not intersect with target shape.
 NoIntersection,
 /// Draft angle is out of valid range.
 InvalidDraftAngle { angle_rad: f64 },
 /// Face index out of range.
 FaceNotFound { face_index: usize },
 /// Neutral plane is invalid.
 InvalidNeutralPlane,
}

impl std::fmt::Display for BRepFeatError {
 fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
 match self {
 Self::NonFiniteInput(name) => write!(f, "{name} must be finite"),
 Self::NonPositiveInput(name) => write!(f, "{name} must be > 0"),
 Self::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
 Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
 Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
 Self::BooleanFailed(err) => write!(f, "boolean operation failed: {err}"),
 Self::ModelingFailed(msg) => write!(f, "modeling operation failed: {msg}"),
 Self::InvalidProfile(msg) => write!(f, "invalid profile: {msg}"),
 Self::NoIntersection => write!(f, "feature does not intersect with target shape"),
 Self::InvalidDraftAngle { angle_rad } => {
 write!(f, "draft angle {:.1} degrees is out of valid range (-89, 89)", angle_rad.to_degrees())
 }
 Self::FaceNotFound { face_index } => write!(f, "face index {face_index} not found"),
 Self::InvalidNeutralPlane => write!(f, "neutral plane definition is invalid"),
 }
 }
}

impl std::error::Error for BRepFeatError {}

impl From<BooleanError> for BRepFeatError {
 fn from(value: BooleanError) -> Self {
 Self::BooleanFailed(value)
 }
}

impl From<rcad_modeling::BuildError> for BRepFeatError {
 fn from(value: rcad_modeling::BuildError) -> Self {
 Self::ModelingFailed(value.to_string())
 }
}

// ===========================================================?
// Fuse Mode and Parameters
// ===========================================================?

/// Defines how a feature interacts with the base shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FuseMode {
 /// Add material to the base shape (boolean union).
 Add,
 /// Remove material from the base shape (boolean difference).
 Remove,
 /// Compute the intersection with the base shape.
 Common,
}

impl From<FuseMode> for BooleanOpType {
 fn from(mode: FuseMode) -> Self {
 match mode {
 FuseMode::Add => BooleanOpType::Union,
 FuseMode::Remove => BooleanOpType::Difference,
 FuseMode::Common => BooleanOpType::Intersection,
 }
 }
}

/// Parameters for feature operations.
#[derive(Debug, Clone)]
pub struct FeatureParams {
 /// Tolerance for merging vertices and edges after the operation.
 ///
 /// [`Default`] uses [`TOLERANCE_MESH_LEGACY`] for backward compatibility.
 /// For phase-C alignment with pairwise BRep stitching, derive a floor from
 /// [`Self::merge_tolerance_for_operands`].
 pub merge_tolerance: f64,
 /// Whether to perform validity checks after the operation.
 pub validate_result: bool,
 /// Whether to simplify the result (merge coplanar faces, etc.).
 pub simplify_result: bool,
}

impl Default for FeatureParams {
 fn default() -> Self {
 Self {
 merge_tolerance: TOLERANCE_MESH_LEGACY,
 validate_result: true,
 simplify_result: true,
 }
 }
}

impl FeatureParams {
 /// Linear merge tolerance aligned with phase-C mesh chaining: **`tessellation_merge_linear_from_two_breps`**
 /// (Relaxed adaptive, [`TOLERANCE_MESH_LEGACY`] minimum, topological fold).
 ///
 /// Use when post eature-result stitching should match [`crate::section::intersect_triangle_soups_adaptive`]-style pipelines.
 pub fn merge_tolerance_for_operands(base: &rcad_kernel::BRep, tool: &rcad_kernel::BRep) -> f64 {
 tessellation_merge_linear_from_two_breps(base, tool)
 }
}

/// Parameters specific to rib features.
#[derive(Debug, Clone)]
pub struct RibParams {
 /// Thickness of the rib (perpendicular to the profile).
 pub thickness: f64,
 /// Height of the rib (along the extrusion direction).
 pub height: f64,
 /// Draft angle for the rib sides (radians).
 pub draft_angle: f64,
 /// Whether to merge the rib with the target shape.
 pub fuse: bool,
}

impl Default for RibParams {
 fn default() -> Self {
 Self {
 thickness: 1.0,
 height: 1.0,
 draft_angle: 0.0,
 fuse: true,
 }
 }
}

/// Parameters specific to groove features.
#[derive(Debug, Clone)]
pub struct GrooveParams {
 /// Depth of the groove.
 pub depth: f64,
 /// Width of the groove (if applicable).
 pub width: Option<f64>,
 /// Whether the groove should go through the entire shape.
 pub through_all: bool,
 /// Draft angle for the groove sides (radians).
 pub draft_angle: f64,
}

impl Default for GrooveParams {
 fn default() -> Self {
 Self {
 depth: 1.0,
 width: None,
 through_all: false,
 draft_angle: 0.0,
 }
 }
}

// ===========================================================?
// Helper Functions
// ===========================================================?

const EPS: f64 = TOLERANCE_LEN_MIN;

fn validate_finite(name: &'static str, v: f64) -> Result<f64, BRepFeatError> {
 if v.is_finite() {
 Ok(v)
 } else {
 Err(BRepFeatError::NonFiniteInput(name))
 }
}

fn validate_positive(name: &'static str, v: f64) -> Result<f64, BRepFeatError> {
 let v = validate_finite(name, v)?;
 if v > 0.0 {
 Ok(v)
 } else {
 Err(BRepFeatError::NonPositiveInput(name))
 }
}

fn normalize(name: &'static str, v: DVec3) -> Result<DVec3, BRepFeatError> {
 if !v.is_finite() {
 return Err(BRepFeatError::NonFiniteInput(name));
 }
 if v.length_squared() <= EPS {
 return Err(BRepFeatError::ZeroVector(name));
 }
 Ok(v.normalize())
}

fn axis_ref_basis(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BRepFeatError> {
 let y_axis = normalize("axis", axis)?;
 let ref_dir = normalize("ref_dir", ref_dir)?;
 let x_reject = ref_dir - y_axis * ref_dir.dot(y_axis);
 if x_reject.length_squared() <= EPS {
 return Err(BRepFeatError::ParallelVectors("ref_dir", "axis"));
 }
 let x_axis = x_reject.normalize();
 let z_axis = x_axis.cross(y_axis).normalize();
 Ok((x_axis, y_axis, z_axis))
}

/// Build a prism solid from bottom and top polygon sections.
fn build_prism_from_sections(bot: &[DVec3], top: &[DVec3], dir: DVec3) -> Result<rcad_kernel::BRep, BRepFeatError> {
 let n = bot.len();
 if n < 3 || top.len() != n {
 return Err(BRepFeatError::InvalidInput("section vertex count mismatch or too few vertices".to_string()));
 }

 let mut brep = rcad_kernel::BRep::new();

 // Add vertices: bot[0..n] then top[0..n]
 let mut verts: Vec<ShapeRef> = Vec::with_capacity(2 * n);
 for &p in bot {
 verts.push(brep.add_tvertex(p));
 }
 for &p in top {
 verts.push(brep.add_tvertex(p));
 }

 fn add_line_edge(brep: &mut rcad_kernel::BRep, verts: &[ShapeRef], start: usize, end: usize) -> ShapeRef {
 let p0 = brep.vertex_point(verts[start].index).unwrap();
 let p1 = brep.vertex_point(verts[end].index).unwrap();
 let d = p1 - p0;
 let len = d.length();
 let dir_vec = if len > 0.0 { d / len } else { DVec3::X };
 brep.add_tedge(
 Some(Curve3::Line(Line3 { origin: p0, direction: dir_vec })),
 verts[start],
 ShapeRef::synthetic_with_orientation(verts[end].index, Orientation::Reversed),
 [0.0, len],
 )
 }

 // Bottom-cap edges
 let bot_edges: Vec<ShapeRef> = (0..n).map(|i| add_line_edge(&mut brep, &verts, i, (i + 1) % n)).collect();
 // Top-cap edges
 let top_edges: Vec<ShapeRef> = (0..n).map(|i| add_line_edge(&mut brep, &verts, n + i, n + (i + 1) % n)).collect();
 // Vertical edges
 let vert_edges: Vec<ShapeRef> = (0..n).map(|i| add_line_edge(&mut brep, &verts, i, n + i)).collect();

 let mut face_refs: Vec<ShapeRef> = Vec::with_capacity(n + 2);

 // Bottom cap (outward normal = -dir)
 {
 let wire_edge_refs: Vec<ShapeRef> = (0..n)
 .map(|i| ShapeRef::synthetic_with_orientation(bot_edges[n - 1 - i].index, Orientation::Reversed))
 .collect();
 let wire_sr = brep.add_twire(wire_edge_refs);
 let face_sr = brep.add_tface(
 Some(Surface3::Plane(Plane { origin: bot[0], normal: -dir })),
 wire_sr,
 vec![],
 None,
 None,
 vec![],
 true,
 );
 face_refs.push(face_sr);
 }

 // Top cap (outward normal = +dir)
 {
 let wire_edge_refs: Vec<ShapeRef> = (0..n)
 .map(|i| ShapeRef::synthetic_with_orientation(top_edges[i].index, Orientation::Forward))
 .collect();
 let wire_sr = brep.add_twire(wire_edge_refs);
 let face_sr = brep.add_tface(
 Some(Surface3::Plane(Plane { origin: top[0], normal: dir })),
 wire_sr,
 vec![],
 None,
 None,
 vec![],
 true,
 );
 face_refs.push(face_sr);
 }

 // Lateral quad faces
 for i in 0..n {
 let j = (i + 1) % n;
 let a = bot[i];
 let b = bot[j];
 let c = top[j];
 let face_normal = {
 let ab = b - a;
 let ac = c - a;
 let nv = ab.cross(ac);
 if nv.length_squared() > TOLERANCE_VEC_SQ_MIN { nv.normalize() } else { -dir.cross(ab).normalize_or(DVec3::X) }
 };
 let wire_edge_refs = vec![
 ShapeRef::synthetic_with_orientation(bot_edges[i].index, Orientation::Forward),
 ShapeRef::synthetic_with_orientation(vert_edges[j].index, Orientation::Forward),
 ShapeRef::synthetic_with_orientation(top_edges[i].index, Orientation::Reversed),
 ShapeRef::synthetic_with_orientation(vert_edges[i].index, Orientation::Reversed),
 ];
 let wire_sr = brep.add_twire(wire_edge_refs);
 let face_sr = brep.add_tface(
 Some(Surface3::Plane(Plane { origin: a, normal: face_normal })),
 wire_sr,
 vec![],
 None,
 None,
 vec![],
 true,
 );
 face_refs.push(face_sr);
 }

 let shell_sr = brep.add_tshell(face_refs);
 brep.add_tsolid(vec![shell_sr]);
 Ok(brep)
}

/// Build a polygon face BRep from vertices.
fn build_polygon_face_brep(profile_verts: &[DVec3]) -> Result<topods::BRep, BRepFeatError> {
 if profile_verts.len() < 3 {
 return Err(BRepFeatError::InvalidInput("profile needs >= 3 vertices".to_string()));
 }

 let n = profile_verts.len();
 let mut brep = rcad_kernel::BRep::new();

 let verts: Vec<ShapeRef> = profile_verts.iter().map(|&p| brep.add_tvertex(p)).collect();

 let normal = {
 let a = profile_verts[0];
 let b = profile_verts[1];
 let c = profile_verts[2];
 let n_vec = (b - a).cross(c - a);
 if n_vec.length_squared() <= EPS {
 return Err(BRepFeatError::InvalidInput("profile vertices are degenerate".to_string()));
 }
 n_vec.normalize()
 };

 // Build edges for the polygon loop
 let edge_refs: Vec<ShapeRef> = (0..n)
 .map(|i| {
 let j = (i + 1) % n;
 let p0 = profile_verts[i];
 let p1 = profile_verts[j];
 let d = p1 - p0;
 let len = d.length();
 let dir_vec = if len > 0.0 { d / len } else { DVec3::X };
 brep.add_tedge(
 Some(Curve3::Line(Line3 { origin: p0, direction: dir_vec })),
 verts[i],
 ShapeRef::synthetic_with_orientation(verts[j].index, Orientation::Reversed),
 [0.0, len],
 )
 })
 .collect();

 // Wire with all edges Forward
 let wire_edge_refs: Vec<ShapeRef> = edge_refs.iter()
 .map(|e| ShapeRef::synthetic_with_orientation(e.index, Orientation::Forward))
 .collect();
 let wire_sr = brep.add_twire(wire_edge_refs);
 let face_sr = brep.add_tface(
 Some(Surface3::Plane(Plane { origin: profile_verts[0], normal })),
 wire_sr,
 vec![],
 None,
 None,
 vec![],
 true,
 );

 let shell_sr = brep.add_tshell(vec![face_sr]);
 brep.add_tsolid(vec![shell_sr]);
 Ok(brep)
}

// ===========================================================?
// Rib Operations
// ===========================================================?

/// Create a rib feature from a wire profile.
///
/// A rib is a thin-wall feature that reinforces a part. The profile wire defines
/// the cross-section of the rib, and it is extruded in the given direction with
/// the specified thickness.
///
/// # Arguments
///
/// * `target` - The base shape to add the rib to.
/// * `profile_wire` - Vertices defining the rib profile (closed polygon).
/// * `direction` - Direction of rib extrusion.
/// * `thickness` - Thickness of the rib perpendicular to the profile.
///
/// # Returns
///
/// The resulting shape with the rib added.
///
/// # Example
///
/// ```ignore
/// let result = make_rib(&base_shape, &profile, DVec3::Y, 2.0)?;
/// ```
pub fn make_rib(
 target: &rcad_kernel::BRep,
 profile_wire: &[DVec3],
 direction: DVec3,
 thickness: f64,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profile_wire.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }

 let dir = normalize("direction", direction)?;
 let thickness = validate_positive("thickness", thickness)?;

 // Extrude the profile in both directions for thickness
 let half_thickness = thickness / 2.0;

 // Find the profile normal
 let a = profile_wire[0];
 let b = profile_wire[1];
 let c = profile_wire[2];
 let profile_normal = (b - a).cross(c - a).normalize();

 // Create two offset copies of the profile
 let profile_offset_neg: Vec<DVec3> = profile_wire.iter()
 .map(|&p| p - profile_normal * half_thickness)
 .collect();
 let profile_offset_pos: Vec<DVec3> = profile_wire.iter()
 .map(|&p| p + profile_normal * half_thickness)
 .collect();

 // Build a thick prism for the rib
 let rib_solid = build_rib_solid(&profile_offset_neg, &profile_offset_pos, dir, thickness)?;

 // Fuse with target — boolean_op already returns topods::BRep
 Ok(boolean_op(BooleanOpType::Union, target, &rib_solid)?)
}

/// Create a linear rib from a profile with specified height.
///
/// Similar to `make_rib` but with explicit control over the rib height.
///
/// # Arguments
///
/// * `target` - The base shape to add the rib to.
/// * `profile` - Vertices defining the rib profile.
/// * `direction` - Direction of rib extrusion.
///
/// # Returns
///
/// The resulting shape with the linear rib added.
pub fn make_linear_rib(
 target: &rcad_kernel::BRep,
 profile: &[DVec3],
 direction: DVec3,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profile.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }

 let dir = normalize("direction", direction)?;

 // Compute profile centroid and height
 let centroid: DVec3 = profile.iter().sum::<DVec3>() / profile.len() as f64;

 // Find the profile extents along the direction
 let heights: Vec<f64> = profile.iter().map(|&p| (p - centroid).dot(dir)).collect();
 let min_h = heights.iter().fold(f64::INFINITY, |a, &b| a.min(b));
 let max_h = heights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
 let profile_height = max_h - min_h;

 // Create a prism from the profile
 let bottom: Vec<DVec3> = profile.iter().map(|&p| p - dir * min_h).collect();
 let top: Vec<DVec3> = bottom.iter().map(|&p| p + dir * profile_height).collect();

 let rib_solid = build_prism_from_sections(&bottom, &top, dir)?;

 // Fuse with target
 Ok(boolean_op(BooleanOpType::Union, target, &rib_solid)?)
}

/// Build a solid rib from two offset profile sections.
fn build_rib_solid(
 bot: &[DVec3],
 top: &[DVec3],
 dir: DVec3,
 _thickness: f64,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 build_prism_from_sections(bot, top, dir)
}

// ===========================================================?
// Groove Operations
// ===========================================================?

/// Create a groove (slot) feature from a profile wire.
///
/// A groove is a depression cut into a part. The profile wire defines the
/// cross-section of the groove.
///
/// # Arguments
///
/// * `target` - The base shape to cut the groove into.
/// * `profile_wire` - Vertices defining the groove profile.
/// * `direction` - Direction of groove extrusion.
/// * `depth` - Depth of the groove.
///
/// # Returns
///
/// The resulting shape with the groove cut.
pub fn make_groove(
 target: &rcad_kernel::BRep,
 profile_wire: &[DVec3],
 direction: DVec3,
 depth: f64,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profile_wire.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }

 let dir = normalize("direction", direction)?;
 let depth = validate_positive("depth", depth)?;

 // Build the groove tool
 let bottom: Vec<DVec3> = profile_wire.to_vec();
 let top: Vec<DVec3> = bottom.iter().map(|&p| p + dir * depth).collect();

 let groove_tool = build_prism_from_sections(&bottom, &top, dir)?;

 // Subtract from target
 Ok(boolean_op(BooleanOpType::Difference, target, &groove_tool)?)
}

/// Create a through groove (slot) that goes through the entire shape.
///
/// Similar to `make_groove` but the groove extends through the entire target shape.
///
/// # Arguments
///
/// * `target` - The base shape to cut the groove into.
/// * `profile` - Vertices defining the groove profile.
/// * `direction` - Direction of groove extrusion.
///
/// # Returns
///
/// The resulting shape with the through groove cut.
pub fn make_through_groove(
 target: &rcad_kernel::BRep,
 profile: &[DVec3],
 direction: DVec3,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profile.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }

 let dir = normalize("direction", direction)?;

 // Compute the bounding box of the target to determine the through distance
 let (min_pt, max_pt) = compute_bounding_box(target);
 let extent = (max_pt - min_pt).length() + 1.0;

 // Build the groove tool that extends beyond the target
 let bottom: Vec<DVec3> = profile.iter().map(|&p| p - dir * extent).collect();
 let top: Vec<DVec3> = profile.iter().map(|&p| p + dir * extent).collect();

 let groove_tool = build_prism_from_sections(&bottom, &top, dir)?;

 // Subtract from target
 Ok(boolean_op(BooleanOpType::Difference, target, &groove_tool)?)
}

/// Compute the axis-aligned bounding box of a BRep.
fn compute_bounding_box(brep: &rcad_kernel::BRep) -> (DVec3, DVec3) {
 let vcount = brep.vertex_count();
 if vcount == 0 {
 return (DVec3::ZERO, DVec3::ZERO);
 }

 let mut min_pt = DVec3::INFINITY;
 let mut max_pt = DVec3::NEG_INFINITY;
 for ts in &brep.tshapes {
 if let TShape::Vertex(vd) = ts.as_ref() {
 min_pt = min_pt.min(vd.point);
 max_pt = max_pt.max(vd.point);
 }
 }
 (min_pt, max_pt)
}

// ===========================================================?
// Prism Feature
// ===========================================================?

/// Create a feature-based prism from a profile.
///
/// Creates a prism by extruding a profile in the given direction and combining
/// it with the target shape using the specified fuse mode.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the prism profile.
/// * `direction` - Direction of extrusion.
/// * `fuse_mode` - How to combine with the target (Add, Remove, Common).
///
/// # Returns
///
/// The resulting shape after the prism operation.
pub fn make_prism_feature(
 target: &rcad_kernel::BRep,
 profile: &[DVec3],
 direction: DVec3,
 fuse_mode: FuseMode,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profile.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }

 let dir = normalize("direction", direction)?;

 // Compute extrusion depth based on profile extents
 let centroid: DVec3 = profile.iter().sum::<DVec3>() / profile.len() as f64;
 let heights: Vec<f64> = profile.iter().map(|&p| (p - centroid).dot(dir)).collect();
 let max_h = heights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

 // Extrude the profile
 let bottom: Vec<DVec3> = profile.to_vec();
 let top: Vec<DVec3> = bottom.iter().map(|&p| p + dir * max_h.abs().max(1.0)).collect();

 let prism_tool = build_prism_from_sections(&bottom, &top, dir)?;

 // Apply boolean operation based on fuse mode
 let op = BooleanOpType::from(fuse_mode);
 Ok(boolean_op(op, target, &prism_tool)?)
}

// ===========================================================?
// Revolution Feature
// ===========================================================?

/// Create a feature-based revolution from a profile.
///
/// Creates a revolution by rotating a profile around an axis and combining
/// it with the target shape using the specified fuse mode.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the revolution profile.
/// * `axis` - Origin point of the revolution axis.
/// * `axis_dir` - Direction of the revolution axis.
/// * `angle` - Revolution angle in radians (full circle = 2*PI).
/// * `fuse_mode` - How to combine with the target (Add, Remove, Common).
///
/// # Returns
///
/// The resulting shape after the revolution operation.
pub fn make_revol_feature(
 target: &topods::BRep,
 profile: &[DVec3],
 axis: DVec3,
 axis_dir: DVec3,
 angle: f64,
 fuse_mode: FuseMode,
) -> Result<topods::BRep, BRepFeatError> {
 if profile.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }

 if !axis.is_finite() {
 return Err(BRepFeatError::NonFiniteInput("axis"));
 }

 let axis_dir = normalize("axis_dir", axis_dir)?;
 let angle = validate_finite("angle", angle)?;

 // Build the profile face
 let profile_brep = build_polygon_face_brep(profile)?;

 // Create the revolution
 let revol_tool = rcad_modeling::revolve(&profile_brep, 0, axis, axis_dir, angle)?;

 // Apply boolean operation based on fuse mode — both target and revol_tool are topods::BRep
 let op = BooleanOpType::from(fuse_mode);
 Ok(boolean_op(op, target, &revol_tool)?)
}

// ===========================================================?
// Pipe Feature
// ===========================================================?

/// Create a feature-based pipe along a spine curve.
///
/// Creates a pipe by sweeping a profile along a spine path and combining
/// it with the target shape using the specified fuse mode.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the pipe cross-section.
/// * `spine` - Vertices defining the spine path.
/// * `fuse_mode` - How to combine with the target (Add, Remove, Common).
///
/// # Returns
///
/// The resulting shape after the pipe operation.
pub fn make_pipe_feature(
 target: &rcad_kernel::BRep,
 profile: &[DVec3],
 spine: &[DVec3],
 fuse_mode: FuseMode,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profile.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }
 if spine.len() < 2 {
 return Err(BRepFeatError::InvalidInput("spine needs >= 2 vertices".to_string()));
 }

 // Build the pipe by sweeping the profile along the spine
 let pipe_tool = build_pipe_solid(profile, spine)?;

 // Apply boolean operation based on fuse mode
 let op = BooleanOpType::from(fuse_mode);
 Ok(boolean_op(op, target, &pipe_tool)?)
}

/// Build a pipe solid by sweeping a profile along a spine.
fn build_pipe_solid(profile: &[DVec3], spine: &[DVec3]) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if spine.len() < 2 {
 return Err(BRepFeatError::InvalidInput("spine needs at least 2 points".to_string()));
 }

 // Compute the spine direction at each point
 let mut frames: Vec<(DVec3, DVec3, DVec3)> = Vec::with_capacity(spine.len());

 for i in 0..spine.len() {
 let tangent = if i == 0 {
 (spine[1] - spine[0]).normalize_or(DVec3::Z)
 } else if i == spine.len() - 1 {
 (spine[i] - spine[i - 1]).normalize_or(DVec3::Z)
 } else {
 (spine[i + 1] - spine[i - 1]).normalize_or(DVec3::Z)
 };

 // Build a local coordinate frame
 let up = if tangent.cross(DVec3::Y).length() > 0.1 {
 tangent.cross(DVec3::Y).normalize()
 } else {
 tangent.cross(DVec3::X).normalize()
 };
 let right = tangent.cross(up).normalize();

 frames.push((spine[i], right, up));
 }

 // Build cross-sections at each spine point
 let sections: Vec<Vec<DVec3>> = frames.iter().map(|(origin, right, up)| {
 profile.iter().map(|&p| {
 origin + right * p.x + up * p.y
 }).collect()
 }).collect();

 // Build the pipe by lofting through sections
 build_loft_solid(&sections)
}

/// Build a loft solid through multiple sections.
pub(crate) fn build_loft_solid(sections: &[Vec<DVec3>]) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if sections.len() < 2 {
 return Err(BRepFeatError::InvalidInput("need at least 2 sections for loft".to_string()));
 }

 let n = sections[0].len();
 if n < 3 {
 return Err(BRepFeatError::InvalidInput("each section needs at least 3 vertices".to_string()));
 }

 // Check all sections have the same number of vertices
 for (i, s) in sections.iter().enumerate() {
 if s.len() != n {
 return Err(BRepFeatError::InvalidInput(format!(
 "section {} has {} vertices, expected {}", i, s.len(), n
 )));
 }
 }

 let mut brep = rcad_kernel::BRep::new();

 let num_sections = sections.len();

 // Add all vertices
 let mut verts: Vec<ShapeRef> = Vec::with_capacity(num_sections * n);
 for si in 0..num_sections {
 for &p in &sections[si] {
 verts.push(brep.add_tvertex(p));
 }
 }

 // Helper: add a line edge between two vertex indices
 let add_line_edge = |brep: &mut rcad_kernel::BRep, start: usize, end: usize| -> ShapeRef {
 let p0 = brep.vertex_point(verts[start].index).unwrap();
 let p1 = brep.vertex_point(verts[end].index).unwrap();
 let d = p1 - p0;
 let len = d.length();
 let dir_vec = if len > 0.0 { d / len } else { DVec3::X };
 brep.add_tedge(
 Some(Curve3::Line(Line3 { origin: p0, direction: dir_vec })),
 verts[start],
 ShapeRef::synthetic_with_orientation(verts[end].index, Orientation::Reversed),
 [0.0, len],
 )
 };

 // Build edge tables
 // Cap edges for each section
 let mut cap_edges: Vec<Vec<ShapeRef>> = Vec::with_capacity(num_sections);
 for si in 0..num_sections {
 let base = si * n;
 let mut edges: Vec<ShapeRef> = Vec::with_capacity(n);
 for i in 0..n {
 let ed = add_line_edge(&mut brep, base + i, base + (i + 1) % n);
 edges.push(ed);
 }
 cap_edges.push(edges);
 }

 // Lateral edges
 let mut lateral_edges: Vec<Vec<ShapeRef>> = Vec::with_capacity(num_sections - 1);
 for si in 0..num_sections - 1 {
 let base0 = si * n;
 let base1 = (si + 1) * n;
 let mut edges: Vec<ShapeRef> = Vec::with_capacity(n);
 for i in 0..n {
 let ed = add_line_edge(&mut brep, base0 + i, base1 + i);
 edges.push(ed);
 }
 lateral_edges.push(edges);
 }

 // Build faces
 let mut face_refs = Vec::new();

 // Bottom cap
 let bottom_normal = {
 let a = sections[0][0];
 let b = sections[0][1];
 let c = sections[0][2];
 (b - a).cross(c - a).normalize_or(DVec3::Z)
 };
 {
 let wire_edge_refs: Vec<ShapeRef> = (0..n)
 .map(|i| ShapeRef::synthetic_with_orientation(cap_edges[0][n - 1 - i].index, Orientation::Reversed))
 .collect();
 let wire_sr = brep.add_twire(wire_edge_refs);
 let face_sr = brep.add_tface(
 Some(Surface3::Plane(Plane { origin: sections[0][0], normal: bottom_normal })),
 wire_sr,
 vec![],
 None,
 None,
 vec![],
 true,
 );
 face_refs.push(face_sr);
 }

 // Top cap
 let top_normal = {
 let a = sections[num_sections - 1][0];
 let b = sections[num_sections - 1][1];
 let c = sections[num_sections - 1][2];
 (b - a).cross(c - a).normalize_or(DVec3::Z)
 };
 {
 let wire_edge_refs: Vec<ShapeRef> = (0..n)
 .map(|i| ShapeRef::synthetic_with_orientation(cap_edges[num_sections - 1][i].index, Orientation::Forward))
 .collect();
 let wire_sr = brep.add_twire(wire_edge_refs);
 let face_sr = brep.add_tface(
 Some(Surface3::Plane(Plane { origin: sections[num_sections - 1][0], normal: top_normal })),
 wire_sr,
 vec![],
 None,
 None,
 vec![],
 true,
 );
 face_refs.push(face_sr);
 }

 // Lateral faces (quads between sections)
 for si in 0..num_sections - 1 {
 for i in 0..n {
 let j = (i + 1) % n;

 // Compute face normal
 let p0 = sections[si][i];
 let p1 = sections[si][j];
 let p2 = sections[si + 1][j];
 let normal = (p1 - p0).cross(p2 - p0).normalize_or(top_normal);

 let wire_edge_refs = vec![
 ShapeRef::synthetic_with_orientation(cap_edges[si][i].index, Orientation::Forward),
 ShapeRef::synthetic_with_orientation(lateral_edges[si][j].index, Orientation::Forward),
 ShapeRef::synthetic_with_orientation(cap_edges[si + 1][i].index, Orientation::Reversed),
 ShapeRef::synthetic_with_orientation(lateral_edges[si][i].index, Orientation::Reversed),
 ];

 let wire_sr = brep.add_twire(wire_edge_refs);
 let face_sr = brep.add_tface(
 Some(Surface3::Plane(Plane { origin: p0, normal })),
 wire_sr,
 vec![],
 None,
 None,
 vec![],
 true,
 );
 face_refs.push(face_sr);
 }
 }

 let shell_sr = brep.add_tshell(face_refs);
 brep.add_tsolid(vec![shell_sr]);
 Ok(brep)
}

// ===========================================================?
// Draft Feature
// ===========================================================?

/// Parameters for draft application.
#[derive(Debug, Clone)]
pub struct DraftFeatureParams {
 /// Draft angle in radians.
 pub angle: f64,
 /// Neutral plane origin point.
 pub neutral_point: DVec3,
 /// Neutral plane normal (pull direction).
 pub pull_direction: DVec3,
}

impl Default for DraftFeatureParams {
 fn default() -> Self {
 Self {
 angle: 2.0_f64.to_radians(),
 neutral_point: DVec3::ZERO,
 pull_direction: DVec3::Z,
 }
 }
}

/// Collect all vertex points from a BRep into a Vec<DVec3> indexed by vertex ordinal.
fn collect_vertex_points(brep: &rcad_kernel::BRep) -> Vec<DVec3> {
 brep.tshapes.iter().filter_map(|ts| {
 if let TShape::Vertex(vd) = ts.as_ref() { Some(vd.point) } else { None }
 }).collect()
}

/// Get the first solid's first shell's face list, along with the face count.
/// This mirrors the old `brep.solids[0].shells[0].faces` pattern.
fn first_shell_faces(brep: &rcad_kernel::BRep) -> Option<Vec<ShapeRef>> {
 for ts in &brep.tshapes {
 if let TShape::Solid(sd) = ts.as_ref() {
 let shell_sr = sd.shells.first()?;
 if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
 return Some(shd.faces.clone());
 }
 }
 }
 None
}

/// Extract edge data (first vertex index, last vertex index) from a BRep edge at a given index.
fn edge_vertex_indices(brep: &rcad_kernel::BRep, idx: usize) -> Option<(usize, usize)> {
 let ts = brep.tshapes.get(idx)?;
 if let TShape::Edge(ed) = ts.as_ref() {
 Some((ed.first.index, ed.last.index))
 } else {
 None
 }
}

/// Get the wire edge refs for a given face's outer wire.
fn face_outer_wire_edges(brep: &rcad_kernel::BRep, face_sr: &ShapeRef) -> Option<Vec<ShapeRef>> {
 let ts = brep.tshapes.get(face_sr.index)?;
 if let TShape::Face(fd) = ts.as_ref() {
 let wire_sr = fd.outer_wire;
 let wts = brep.tshapes.get(wire_sr.index)?;
 if let TShape::Wire(wd) = wts.as_ref() {
 return Some(wd.edges.clone());
 }
 }
 None
}

/// Apply draft angle to specified faces.
///
/// Draft angle is applied to allow parts to be removed from molds. The draft
/// tilts the faces so that the part can be pulled out in the pull direction.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `face_indices` - Indices of faces to apply draft to.
/// * `angle` - Draft angle in radians. Positive values add draft for easy removal.
/// * `neutral_plane` - A point on the neutral plane (vertices here stay fixed).
///
/// # Returns
///
/// The shape with draft applied.
pub fn apply_draft_feature(
 target: &rcad_kernel::BRep,
 face_indices: &[usize],
 angle: f64,
 neutral_plane: DVec3,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 // Validate angle
 if angle.abs() >= std::f64::consts::FRAC_PI_2 - TOLERANCE_MESH_LEGACY {
 return Err(BRepFeatError::InvalidDraftAngle { angle_rad: angle });
 }

 let faces = first_shell_faces(target)
 .ok_or_else(|| BRepFeatError::InvalidInput("target has no solids".to_string()))?;

 // Validate face indices
 for &fi in face_indices {
 if fi >= faces.len() {
 return Err(BRepFeatError::FaceNotFound { face_index: fi });
 }
 }

 // If no faces specified, apply draft to all non-horizontal faces
 let faces_to_draft: Vec<usize> = if face_indices.is_empty() {
 // Find all faces that appear non-horizontal
 // Use surface normal heuristic from target data
 faces.iter().enumerate()
 .filter(|(_, face_sr)| {
 if let Some(wire_edges) = face_outer_wire_edges(target, face_sr) {
 // Check first edge to guess face orientation
 if let Some(edge_idx) = wire_edges.first().map(|e| e.index) {
 if let Some(TShape::Edge(ed)) = target.tshapes.get(edge_idx).map(|ts| ts.as_ref()) {
 // Rough check: if the first edge's start and end vertices have different Z,
 // the face is likely non-horizontal
 let p0 = target.vertex_point(ed.first.index).unwrap_or(DVec3::ZERO);
 let p1 = target.vertex_point(ed.last.index).unwrap_or(DVec3::ZERO);
 let edge_dir = (p1 - p0).normalize_or(DVec3::Z);
 return edge_dir.dot(DVec3::Z).abs() < 0.99;
 }
 }
 }
 false
 })
 .map(|(i, _)| i)
 .collect()
 } else {
 face_indices.to_vec()
 };

 if faces_to_draft.is_empty() {
 // No faces to draft, return a clone
 return Ok(target.clone());
 }

 // Apply draft transformation
 let pull_dir = DVec3::Z; // Default pull direction
 let tan_angle = angle.tan();

 // Compute new vertex positions
 let vcount = target.vertex_count();
 let mut new_positions: Vec<DVec3> = collect_vertex_points(target);

 for &fi in &faces_to_draft {
 if let Some(face_sr) = faces.get(fi) {
 if let Some(wire_edges) = face_outer_wire_edges(target, face_sr) {
 for we in &wire_edges {
 if let Some((vi_start, vi_end)) = edge_vertex_indices(target, we.index) {
 for &vi in &[vi_start, vi_end] {
 if vi < vcount {
 if let Some(v) = target.vertex_point(vi) {
 let h = (v - neutral_plane).dot(pull_dir);
 // Apply draft displacement
 let radial_dir = (v - neutral_plane).reject_from(pull_dir).normalize_or(DVec3::ZERO);
 if radial_dir.length() > TOLERANCE_LINEAR_ULTRA_STRICT {
 new_positions[vi] = v + radial_dir * (h * tan_angle);
 }
 }
 }
 }
 }
 }
 }
 }
 }

 // Build the new BRep with modified vertex positions
 build_drafted_brep(target, &new_positions)
}

/// OCCT DRAW `depouille` / `BRepOffsetAPI_DraftAngle` approximation: one pull direction and a
/// sequence of per-face blocks `(face_index, angle_rad, neutral_point, neutral_plane_normal)`.
///
/// Face indices match `rcad-kernel` box face order (see `occt-test-gen` OCCT `bx_1` `bx_6` map).
///
/// Multi-face blocks: each vertex displacement is the **average** of the per-face draft
/// suggestions (from the original vertex position). This avoids double-counting edge vertices
/// when several faces are drafted with the same pull direction (OCCT `depouille` scripts).
pub fn apply_depouille(
 target: &rcad_kernel::BRep,
 pull_direction: DVec3,
 blocks: &[(usize, f64, DVec3, DVec3)],
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 let pull = pull_direction.normalize_or(DVec3::Z);
 if pull.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
 return Err(BRepFeatError::ZeroVector("pull_direction"));
 }

 let faces = first_shell_faces(target)
 .ok_or_else(|| BRepFeatError::InvalidInput("target has no solids".to_string()))?;

 let vcount = target.vertex_count();
 let mut sum_disp = vec![DVec3::ZERO; vcount];
 let mut contrib = vec![0u32; vcount];

 for &(face_index, angle, neutral_point, neutral_normal) in blocks {
 if angle.abs() >= std::f64::consts::FRAC_PI_2 - TOLERANCE_MESH_LEGACY {
 return Err(BRepFeatError::InvalidDraftAngle { angle_rad: angle });
 }
 let n = neutral_normal.normalize_or(DVec3::Z);
 if n.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
 return Err(BRepFeatError::ZeroVector("neutral_normal"));
 }
 if face_index >= faces.len() {
 return Err(BRepFeatError::FaceNotFound { face_index });
 }

 let tan_a = angle.tan();
 if let Some(face_sr) = faces.get(face_index) {
 if let Some(wire_edges) = face_outer_wire_edges(target, face_sr) {
 for we in &wire_edges {
 if let Some((vi_start, vi_end)) = edge_vertex_indices(target, we.index) {
 for &vi in &[vi_start, vi_end] {
 if vi < vcount {
 if let Some(v0) = target.vertex_point(vi) {
 let disp = draft_vertex_displacement_occt(v0, neutral_point, n, pull, tan_a);
 sum_disp[vi] += disp;
 contrib[vi] += 1;
 }
 }
 }
 }
 }
 }
 }
 }

 let new_positions: Vec<DVec3> = collect_vertex_points(target)
 .into_iter()
 .enumerate()
 .map(|(i, v)| {
 let d = if contrib[i] > 0 {
 sum_disp[i] / contrib[i] as f64
 } else {
 DVec3::ZERO
 };
 v + d
 })
 .collect();
 build_drafted_brep(target, &new_positions)
}

/// Draft displacement for one vertex: neutral plane \((P,\hat n)\), pull \(\hat D\), angle \(\theta\).
fn draft_vertex_displacement_occt(
 v: DVec3,
 neutral_point: DVec3,
 n: DVec3,
 pull: DVec3,
 tan_a: f64,
) -> DVec3 {
 let w = v - neutral_point;
 let h = w.dot(pull);
 let radial = w - pull * h;
 const EPS_PLANE: f64 = TOLERANCE_ABS;
 const EPS_PAR: f64 = TOLERANCE_MESH_LEGACY;

 // Face lies in (or parallel to) neutral plane and pull is along plane normal: in-plane taper
 // (same sign convention as OCCT `depouille` / `checkprops -s` on `tests/draft/angle/A1`).
 if pull.dot(n).abs() >= 1.0 - EPS_PAR && w.dot(n).abs() < EPS_PLANE {
 let tang = w - n * w.dot(n);
 if tang.length_squared() > TOLERANCE_METRIC_SQ_NEAR_ZERO {
 return tang * tan_a;
 }
 return DVec3::ZERO;
 }

 if radial.length_squared() > TOLERANCE_METRIC_SQ_NEAR_ZERO {
 return radial.normalize() * (h * tan_a);
 }
 DVec3::ZERO
}

/// Build a new BRep with drafted vertex positions.
/// Clones the original BRep (preserving all surfaces, curves, and topology)
/// and updates only the vertex positions.
fn build_drafted_brep(
 original: &rcad_kernel::BRep,
 new_positions: &[DVec3],
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 let mut brep = original.clone();
 // Find vertex TShape indices and update their points
 let mut vi = 0;
 for ts in &mut brep.tshapes {
 let shape = Arc::make_mut(ts);
 match shape {
 TShape::Vertex(vd) => {
 if vi < new_positions.len() {
 vd.point = new_positions[vi];
 }
 vi += 1;
 }
 _ => {}
 }
 }
 Ok(brep)
}

// ===========================================================?
// Advanced Feature Operations
// ===========================================================?

/// Create a drafted prism feature (tapered extrusion).
///
/// Creates a prism with draft angle, useful for molded parts.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the prism profile.
/// * `direction` - Direction of extrusion.
/// * `depth` - Extrusion depth.
/// * `draft_angle` - Draft angle in radians.
/// * `fuse_mode` - How to combine with the target.
///
/// # Returns
///
/// The resulting shape with the drafted prism.
pub fn make_drafted_prism(
 target: &rcad_kernel::BRep,
 profile: &[DVec3],
 direction: DVec3,
 depth: f64,
 draft_angle: f64,
 fuse_mode: FuseMode,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profile.len() < 3 {
 return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
 }

 let dir = normalize("direction", direction)?;
 let depth = validate_positive("depth", depth)?;

 if draft_angle.abs() >= std::f64::consts::FRAC_PI_2 - TOLERANCE_MESH_LEGACY {
 return Err(BRepFeatError::InvalidDraftAngle { angle_rad: draft_angle });
 }

 // Compute centroid
 let centroid: DVec3 = profile.iter().sum::<DVec3>() / profile.len() as f64;

 // Apply draft by scaling the top profile
 let taper = depth * draft_angle.tan();

 let bottom: Vec<DVec3> = profile.to_vec();
 let top: Vec<DVec3> = profile.iter().map(|&p| {
 let radial = p - centroid;
 let radial_2d = radial - dir * radial.dot(dir);
 let radial_dir = if radial_2d.length() > EPS {
 radial_2d.normalize()
 } else {
 DVec3::ZERO
 };
 p + dir * depth + radial_dir * taper
 }).collect();

 let prism_tool = build_prism_from_sections(&bottom, &top, dir)?;

 let op = BooleanOpType::from(fuse_mode);
 Ok(boolean_op(op, target, &prism_tool)?)
}

/// Create a multi-profile pipe (loft) feature.
///
/// Creates a solid by lofting through multiple profiles.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profiles` - Vector of profiles, each being a vector of vertices.
/// * `fuse_mode` - How to combine with the target.
///
/// # Returns
///
/// The resulting shape with the loft feature.
pub fn make_loft_feature(
 target: &rcad_kernel::BRep,
 profiles: &[Vec<DVec3>],
 fuse_mode: FuseMode,
) -> Result<rcad_kernel::BRep, BRepFeatError> {
 if profiles.len() < 2 {
 return Err(BRepFeatError::InvalidInput("need at least 2 profiles for loft".to_string()));
 }

 for (i, profile) in profiles.iter().enumerate() {
 if profile.len() < 3 {
 return Err(BRepFeatError::InvalidProfile(format!(
 "profile {} needs >= 3 vertices", i
 )));
 }
 }

 let loft_tool = build_loft_solid(profiles)?;

 let op = BooleanOpType::from(fuse_mode);
 Ok(boolean_op(op, target, &loft_tool)?)
}

// ===========================================================?
// Tests
// ===========================================================?
