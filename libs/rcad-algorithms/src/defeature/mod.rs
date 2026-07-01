//! Defeaturing pass: suppress small cylindrical holes, bosses, and very small faces.
//!
//! Analogous to `BRepAlgoAPI_Defeaturing` in OCCT 8.0.
//!
//! # Overview
//!
//! The defeaturing pass identifies and removes small geometric features from a
//! B-Rep solid that are irrelevant to downstream analysis (meshing, simulation,
//! manufacturing tolerancing).  The baseline implementation handles:
//!
//! - **Cylindrical holes** (through-holes and blind holes): detected by finding groups
//!   of connected cylindrical faces whose radius is below `max_hole_radius`.  The hole
//!   is filled by boolean-unioning a capped cylinder solid into the host body.
//!
//! - **Cylindrical bosses** (protruding cylinders): same detection, opposite normal
//!   direction, filled by boolean-differencing the boss cylinder from the host body.
//!
//! - **Conical holes/bosses**: similar detection for conical features.
//!
//! - **Small-face identification**: faces whose approximate polygon area is below
//!   `max_small_face_area` are reported (see [`identify_small_faces`]).  Removal is
//!   left to the caller because patching isolated small faces without topology
//!   information is highly geometry-specific.
//!
//! # Enhanced Features
//!
//! The enhanced implementation also supports:
//! - **Retry mechanism**: Boolean failures trigger fuzzy tolerance escalation.
//! - **Topology healing**: Post-defeature connectivity repair.
//! - **Feature group detection**: Connected compound features are handled together.
//! - **Slot/pocket detection**: Rectangular and circular slot features.
//! - **Blend/chamfer detection**: Fillets and chamfers below a size threshold.
//!
//! # Usage
//!
//! ```rust
//! use glam::DVec3;
//! use rcad_algorithms::{DefeaturingOptions, defeature_brep};
//! use rcad_modeling::make_box_brep;
//!
//! let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 8.0, 6.0).unwrap();
//! let opts = DefeaturingOptions {
//!     max_hole_radius: 5.0,  // fill holes <= 5 mm radius
//!     ..Default::default()
//! };
//!
//! let (_defeatured, report) = defeature_brep(&brep, &opts).unwrap();
//! assert_eq!(report.holes_removed, 0);
//! ```

use std::collections::{HashMap, HashSet, VecDeque};

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{ConicalSurface, CylindricalSurface, Plane, SphericalSurface, Surface3, ToroidalSurface, any_perpendicular};
use rcad_kernel::topology::{Face, Wire};
use rcad_modeling::make_cylinder_brep;

use crate::tolerance::*;
use crate::{BooleanOpType, BooleanOptions, boolean_op, boolean_op_robust, BooleanRobustOptions, BooleanRetryPolicy};
use crate::brep_repair::make_connected_enhanced;

// -- Tolerances --------------------------------------------------------------

/// Maximum cross-product magnitude for two normalized axis vectors to be
/// considered parallel.
const AXIS_PARALLEL_TOL: f64 = TOLERANCE_RETRY_LADDER_MID;

/// Maximum allowable difference in cylinder radii to be grouped together.
const RADIUS_TOL: f64 = TOLERANCE_RETRY_LADDER_MID;

/// Default fill margin along the axis: how much the fill solid extends beyond
/// the detected hole extent to ensure a clean boolean union.
const DEFAULT_FILL_MARGIN: f64 = TOLERANCE_ABS * 4.0;

// -- Public types ------------------------------------------------------------

/// A detected cylindrical feature (hole or boss) in a B-Rep solid.
///
/// Produced by [`detect_cylindrical_features`].
#[derive(Debug, Clone)]
pub struct CylindricalFeature {
    /// Local face indices *within `solids[0].shells[0].faces`* that make up
    /// the cylindrical wall of this feature.
    pub face_indices: Vec<usize>,

    /// `true` if this is a hole (the material surrounds the cylinder from the
    /// outside; the cylindrical face normal points toward the axis).
    /// `false` if this is a boss (the material is inside the cylinder; the
    /// normal points away from the axis).
    pub is_hole: bool,

    /// A point on the cylinder axis (taken from the underlying surface origin).
    pub origin: DVec3,

    /// Normalized cylinder axis direction.
    pub axis: DVec3,

    /// Cylinder radius.
    pub radius: f64,

    /// Minimum parametric extent along `axis` from `origin` (in model units).
    /// Computed by projecting all wall-face vertex positions onto the axis.
    pub t_min: f64,

    /// Maximum parametric extent along `axis` from `origin`.
    pub t_max: f64,
}

impl CylindricalFeature {
    /// Height of the feature along the cylinder axis.
    pub fn height(&self) -> f64 {
        (self.t_max - self.t_min).max(0.0)
    }
}

/// A detected conical feature (hole or boss) in a B-Rep solid.
#[derive(Debug, Clone)]
pub struct ConicalFeature {
    /// Local face indices within the shell that make up the conical wall.
    pub face_indices: Vec<usize>,
    /// True if this is a hole (material surrounds the cone from outside).
    pub is_hole: bool,
    /// Apex point of the cone.
    pub apex: DVec3,
    /// Normalized axis direction.
    pub axis: DVec3,
    /// Reference radius at a specific height.
    pub reference_radius: f64,
    /// Half angle in radians.
    pub half_angle: f64,
    /// Minimum parametric extent along axis from apex.
    pub t_min: f64,
    /// Maximum parametric extent along axis from apex.
    pub t_max: f64,
}

/// A detected slot feature in a B-Rep solid.
///
/// Slots are elongated recesses or protrusions, typically with rectangular
/// or rounded cross-sections.
#[derive(Debug, Clone)]
pub struct SlotFeature {
    /// Local face indices within the shell that make up the slot.
    pub face_indices: Vec<usize>,
    /// True if this is a recess (slot), false if protrusion.
    pub is_recess: bool,
    /// Slot length along the major direction.
    pub length: f64,
    /// Slot width.
    pub width: f64,
    /// Slot depth (for recesses) or height (for protrusions).
    pub depth: f64,
    /// Origin point at the center of the slot bottom.
    pub origin: DVec3,
    /// Direction along the slot length.
    pub length_dir: DVec3,
    /// Direction along the slot width.
    pub width_dir: DVec3,
    /// Direction along the slot depth.
    pub depth_dir: DVec3,
    /// Whether the slot has rounded ends (cylindrical end caps).
    pub has_rounded_ends: bool,
}

/// A detected pocket feature in a B-Rep solid.
///
/// Pockets are enclosed recesses, typically with flat bottoms and
/// vertical or drafted side walls.
#[derive(Debug, Clone)]
pub struct PocketFeature {
    /// Local face indices within the shell that make up the pocket.
    pub face_indices: Vec<usize>,
    /// True if this is a pocket (recess), false if a pad (protrusion).
    pub is_recess: bool,
    /// Pocket diameter for circular pockets, or max dimension for rectangular.
    pub diameter: f64,
    /// Pocket depth.
    pub depth: f64,
    /// Center point of the pocket opening.
    pub center: DVec3,
    /// Normal direction pointing out of the pocket.
    pub normal: DVec3,
    /// Whether the pocket is circular (true) or rectangular (false).
    pub is_circular: bool,
    /// Width for rectangular pockets (0.0 for circular).
    pub width: f64,
    /// Length for rectangular pockets (0.0 for circular).
    pub length: f64,
    /// Whether the pocket is a through-pocket (passes through the solid).
    pub is_through: bool,
    /// Bottom face index if available (for blind pockets).
    pub bottom_face_index: Option<usize>,
    /// Face indices of the side walls.
    pub wall_face_indices: Vec<usize>,
}

/// A detected blend (fillet) or chamfer feature in a B-Rep solid.
#[derive(Debug, Clone)]
pub struct BlendFeature {
    /// Local face indices within the shell that make up the blend.
    pub face_indices: Vec<usize>,
    /// True if this is a fillet (curved), false if a chamfer (flat).
    pub is_fillet: bool,
    /// Radius for fillets, 0.0 for chamfers.
    pub radius: f64,
    /// Chamfer distance (both sides equal) for chamfers.
    pub chamfer_distance: f64,
    /// Representative point on the blend.
    pub sample_point: DVec3,
    /// Approximate normal direction at the sample point.
    pub normal: DVec3,
}

/// Configuration for pocket detection.
#[derive(Debug, Clone)]
pub struct PocketDetectionConfig {
    /// Maximum pocket diameter (or max dimension) to consider.
    pub max_diameter: f64,
    /// Maximum pocket depth to consider.
    pub max_depth: f64,
    /// Minimum pocket depth (to filter out shallow recesses).
    pub min_depth: f64,
    /// Tolerance for determining if a pocket is through.
    pub through_tolerance: f64,
    /// Whether to detect rectangular pockets.
    pub detect_rectangular: bool,
    /// Whether to detect circular pockets.
    pub detect_circular: bool,
    /// Minimum aspect ratio (depth/width) for pocket detection.
    pub min_aspect_ratio: f64,
}

impl Default for PocketDetectionConfig {
    fn default() -> Self {
        Self {
            max_diameter: 50.0,
            max_depth: 100.0,
            min_depth: 0.1,
            through_tolerance: TOLERANCE_ABS * 10.0,
            detect_rectangular: true,
            detect_circular: true,
            min_aspect_ratio: 0.01,
        }
    }
}

impl PocketDetectionConfig {
    /// Create config for small features only.
    pub fn small_features() -> Self {
        Self {
            max_diameter: 10.0,
            max_depth: 20.0,
            min_depth: 0.05,
            through_tolerance: TOLERANCE_ABS * 5.0,
            detect_rectangular: true,
            detect_circular: true,
            min_aspect_ratio: 0.01,
        }
    }

    /// Create config for large features.
    pub fn large_features() -> Self {
        Self {
            max_diameter: 200.0,
            max_depth: 500.0,
            min_depth: 1.0,
            through_tolerance: TOLERANCE_ABS * 20.0,
            detect_rectangular: true,
            detect_circular: true,
            min_aspect_ratio: 0.005,
        }
    }
}

/// A detected boss feature in a B-Rep solid.
///
/// Bosses are protruding features, typically cylindrical or rectangular pads.
#[derive(Debug, Clone)]
pub struct BossFeature {
    /// Local face indices within the shell that make up the boss.
    pub face_indices: Vec<usize>,
    /// Boss diameter for circular bosses, or max dimension for rectangular.
    pub diameter: f64,
    /// Boss height (protrusion from base surface).
    pub height: f64,
    /// Center point of the boss base.
    pub base_center: DVec3,
    /// Normal direction of the boss (pointing away from base).
    pub normal: DVec3,
    /// Whether the boss is circular (true) or rectangular (false).
    pub is_circular: bool,
    /// Width for rectangular bosses (0.0 for circular).
    pub width: f64,
    /// Length for rectangular bosses (0.0 for circular).
    pub length: f64,
    /// Face indices of the side walls.
    pub wall_face_indices: Vec<usize>,
    /// Face index of the top face (if available).
    pub top_face_index: Option<usize>,
}

/// A detected fillet feature in a B-Rep solid.
#[derive(Debug, Clone)]
pub struct FilletFeature {
    /// Local face indices within the shell that make up the fillet.
    pub face_indices: Vec<usize>,
    /// Fillet radius.
    pub radius: f64,
    /// Representative point on the fillet.
    pub sample_point: DVec3,
    /// Approximate axis direction for the fillet (for edge fillets).
    pub axis: DVec3,
    /// Whether this is a variable-radius fillet.
    pub is_variable: bool,
    /// Min radius for variable fillets.
    pub min_radius: f64,
    /// Max radius for variable fillets.
    pub max_radius: f64,
    /// Adjacent face indices that the fillet connects.
    pub adjacent_faces: Vec<usize>,
}

/// A detected chamfer feature in a B-Rep solid.
#[derive(Debug, Clone)]
pub struct ChamferFeature {
    /// Local face indices within the shell that make up the chamfer.
    pub face_indices: Vec<usize>,
    /// Chamfer distance (equal on both sides for 45-degree chamfers).
    pub distance: f64,
    /// Second distance for asymmetric chamfers (equal to distance for symmetric).
    pub distance2: f64,
    /// Chamfer angle in radians (PI/4 for 45-degree).
    pub angle: f64,
    /// Representative point on the chamfer.
    pub sample_point: DVec3,
    /// Normal direction of the chamfer face.
    pub normal: DVec3,
    /// Adjacent face indices that the chamfer connects.
    pub adjacent_faces: Vec<usize>,
}

/// Feature type enumeration for unified feature handling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FeatureType {
    /// Cylindrical hole or boss.
    Cylindrical,
    /// Conical feature.
    Conical,
    /// Slot feature.
    Slot,
    /// Pocket feature.
    Pocket,
    /// Boss feature.
    Boss,
    /// Fillet feature.
    Fillet,
    /// Chamfer feature.
    Chamfer,
    /// Blend feature (fillet or chamfer).
    Blend,
}

/// A detected hole pattern (array of similar holes).
///
/// Hole patterns represent groups of similar holes that may be processed
/// together for efficiency or that share geometric relationships.
#[derive(Debug, Clone)]
pub struct HolePattern {
    /// Indices of cylindrical features that form this pattern.
    pub feature_indices: Vec<usize>,
    /// Pattern type: linear, circular, rectangular grid, or irregular.
    pub pattern_type: HolePatternType,
    /// Number of holes in the pattern.
    pub count: usize,
    /// Spacing between holes (for linear/circular patterns).
    pub spacing: f64,
    /// Pattern origin (first hole center or pattern center).
    pub origin: DVec3,
    /// Pattern direction (for linear patterns) or axis (for circular patterns).
    pub direction: DVec3,
    /// Common radius for all holes in the pattern.
    pub common_radius: f64,
    /// Common depth for all holes (0.0 for through-holes).
    pub common_depth: f64,
}

/// Type of hole pattern arrangement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HolePatternType {
    /// Holes arranged in a single line.
    Linear,
    /// Holes arranged in a circle.
    Circular,
    /// Holes arranged in a rectangular grid.
    RectangularGrid,
    /// Holes that don't fit a regular pattern.
    Irregular,
}

/// Feature group representing connected features that should be processed together.
#[derive(Debug, Clone)]
pub struct FeatureGroup {
    /// Group ID.
    pub id: usize,
    /// Cylindrical feature indices in this group.
    pub cylindrical_indices: Vec<usize>,
    /// Conical feature indices in this group.
    pub conical_indices: Vec<usize>,
    /// Slot feature indices in this group.
    pub slot_indices: Vec<usize>,
    /// Pocket feature indices in this group.
    pub pocket_indices: Vec<usize>,
    /// Blend feature indices in this group.
    pub blend_indices: Vec<usize>,
    /// Total number of faces in this group.
    pub total_faces: usize,
}

/// Options controlling the defeaturing pass.
#[derive(Debug, Clone, Copy)]
pub struct DefeaturingOptions {
    /// Maximum radius of cylindrical **holes** to fill.  Set to `0.0` (or any
    /// non-positive value) to skip hole removal.
    pub max_hole_radius: f64,

    /// Maximum radius of cylindrical **bosses** to remove.  Set to `0.0` to
    /// skip boss removal.
    pub max_boss_radius: f64,

    /// Maximum approximate polygon area for a face to be flagged as "small"
    /// by [`identify_small_faces`].  Set to `0.0` to disable.
    pub max_small_face_area: f64,

    /// Safety margin (in model units) added on each side of the fill solid
    /// along the cylinder axis to prevent numerical slivers.
    pub fill_margin: f64,

    /// Enable conical feature detection and removal.
    pub enable_conical_features: bool,

    /// Maximum reference radius for conical holes.
    pub max_conical_hole_radius: f64,

    /// Enable retry mechanism for failed boolean operations.
    pub enable_retry: bool,

    /// Fuzzy tolerance multiplier for retry attempts.
    pub retry_fuzzy_multiplier: f64,

    /// Maximum number of retry attempts per feature.
    pub max_retries: usize,

    /// Run post-defeature connectivity healing.
    pub run_post_healing: bool,

    /// Tolerance for post-defeature healing.
    pub healing_tolerance: f64,

    // -- Slot/Pocket feature options --
    /// Enable slot feature detection and removal.
    pub enable_slot_features: bool,

    /// Maximum slot width to consider for removal.
    pub max_slot_width: f64,

    /// Maximum slot depth to consider for removal.
    pub max_slot_depth: f64,

    /// Enable pocket feature detection and removal.
    pub enable_pocket_features: bool,

    /// Maximum pocket diameter (or max dimension) for removal.
    pub max_pocket_diameter: f64,

    /// Maximum pocket depth for removal.
    pub max_pocket_depth: f64,

    // -- Blend/Chamfer feature options --
    /// Enable blend (fillet/chamfer) feature detection.
    pub enable_blend_features: bool,

    /// Maximum blend radius to consider for removal.
    /// Fillets with radius <= this value will be targeted.
    pub max_blend_radius: f64,

    /// Maximum chamfer distance to consider for removal.
    pub max_chamfer_distance: f64,
}

impl Default for DefeaturingOptions {
    fn default() -> Self {
        Self {
            max_hole_radius: 0.0,
            max_boss_radius: 0.0,
            max_small_face_area: 0.0,
            fill_margin: DEFAULT_FILL_MARGIN,
            enable_conical_features: false,
            max_conical_hole_radius: 0.0,
            enable_retry: false,
            retry_fuzzy_multiplier: 10.0,
            max_retries: 3,
            run_post_healing: false,
            healing_tolerance: TOLERANCE_ABS * 10.0,
            // Slot/Pocket defaults
            enable_slot_features: false,
            max_slot_width: 0.0,
            max_slot_depth: 0.0,
            enable_pocket_features: false,
            max_pocket_diameter: 0.0,
            max_pocket_depth: 0.0,
            // Blend defaults
            enable_blend_features: false,
            max_blend_radius: 0.0,
            max_chamfer_distance: 0.0,
        }
    }
}

/// Report produced by [`defeature_brep`].
#[derive(Debug, Clone, Default)]
pub struct DefeaturingReport {
    /// Number of cylindrical holes successfully filled.
    pub holes_removed: usize,

    /// Number of cylindrical bosses successfully removed.
    pub bosses_removed: usize,

    /// Number of conical features removed.
    pub conical_features_removed: usize,

    /// Number of features that were detected but could not be suppressed
    /// (e.g. due to a boolean failure).
    pub failed_features: usize,

    /// Number of retry attempts made.
    pub retry_attempts: usize,

    /// Number of features that succeeded after retry.
    pub succeeded_after_retry: usize,

    /// Number of faces identified as "small" (area <= `max_small_face_area`).
    /// These are *not* removed automatically; use the returned face indices
    /// from [`identify_small_faces`] for targeted treatment.
    pub small_faces_identified: usize,

    /// Whether post-defeature healing was performed.
    pub healing_performed: bool,

    /// Number of vertices merged during healing.
    pub healing_vertices_merged: usize,

    /// Number of small edges removed during healing.
    pub healing_small_edges_removed: usize,

    // -- Slot/Pocket statistics --
    /// Number of slot features removed.
    pub slots_removed: usize,

    /// Number of pocket features removed.
    pub pockets_removed: usize,

    // -- Blend statistics --
    /// Number of blend (fillet/chamfer) features removed.
    pub blends_removed: usize,

    // -- Feature group statistics --
    /// Number of feature groups processed.
    pub feature_groups_processed: usize,

    /// Number of faces that are part of detected feature groups.
    pub grouped_faces: usize,
}

/// Errors returned by the defeaturing pass.
#[derive(Debug)]
pub enum DefeaturingError {
    /// The input BRep has no solids or shells.
    EmptyInput,
    /// Every detected feature failed to be suppressed.
    AllFeaturesFailed,
}

impl std::fmt::Display for DefeaturingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "input BRep has no geometry"),
            Self::AllFeaturesFailed => write!(f, "all detected features failed to be suppressed"),
        }
    }
}

impl std::error::Error for DefeaturingError {}

// -- Internal helpers --------------------------------------------------------

/// Compute the flat face index for a face in `solids[si].shells[shi].faces[fi]`.
fn flat_face_index(brep: &BRep, si: usize, shi: usize, fi: usize) -> usize {
    let mut idx = 0usize;
    for s in 0..si {
        for sh in &brep.solids[s].shells {
            idx += sh.faces.len();
        }
    }
    for sh in 0..shi {
        idx += brep.solids[si].shells[sh].faces.len();
    }
    idx + fi
}

/// Return the `CylindricalSurface` backing a face, or `None` if the face has
/// no surface data or is not a cylinder.
fn face_cylinder(brep: &BRep, si: usize, shi: usize, fi: usize) -> Option<CylindricalSurface> {
    let ffi = flat_face_index(brep, si, shi, fi);
    let sid = brep.geom.face_surface.get(ffi)?.as_ref().copied()?;
    match brep.geom.surfaces.get(sid)? {
        Surface3::Cylinder(c) => Some(*c),
        _ => None,
    }
}

/// Return the `ConicalSurface` backing a face, or `None` if the face has
/// no surface data or is not a cone.
fn face_cone(brep: &BRep, si: usize, shi: usize, fi: usize) -> Option<ConicalSurface> {
    let ffi = flat_face_index(brep, si, shi, fi);
    let sid = brep.geom.face_surface.get(ffi)?.as_ref().copied()?;
    match brep.geom.surfaces.get(sid)? {
        Surface3::Cone(c) => Some(*c),
        _ => None,
    }
}

/// Return the `Plane` backing a face, or `None` if the face has
/// no surface data or is not a plane.
fn face_plane(brep: &BRep, si: usize, shi: usize, fi: usize) -> Option<Plane> {
    let ffi = flat_face_index(brep, si, shi, fi);
    let sid = brep.geom.face_surface.get(ffi)?.as_ref().copied()?;
    match brep.geom.surfaces.get(sid)? {
        Surface3::Plane(p) => Some(*p),
        _ => None,
    }
}

/// Return the `ToroidalSurface` backing a face, or `None` if the face has
/// no surface data or is not a torus.
fn face_torus(brep: &BRep, si: usize, shi: usize, fi: usize) -> Option<ToroidalSurface> {
    let ffi = flat_face_index(brep, si, shi, fi);
    let sid = brep.geom.face_surface.get(ffi)?.as_ref().copied()?;
    match brep.geom.surfaces.get(sid)? {
        Surface3::Torus(t) => Some(*t),
        _ => None,
    }
}

/// Return the `SphericalSurface` backing a face, or `None` if the face has
/// no surface data or is not a sphere.
fn face_sphere(brep: &BRep, si: usize, shi: usize, fi: usize) -> Option<SphericalSurface> {
    let ffi = flat_face_index(brep, si, shi, fi);
    let sid = brep.geom.face_surface.get(ffi)?.as_ref().copied()?;
    match brep.geom.surfaces.get(sid)? {
        Surface3::Sphere(s) => Some(*s),
        _ => None,
    }
}

/// Check if a face is likely a planar blend/chamfer face.
/// Returns Some((is_fillet, radius_or_chamfer_dist)) if detected.
fn detect_blend_face(brep: &BRep, si: usize, shi: usize, fi: usize, max_blend_radius: f64, max_chamfer_distance: f64) -> Option<BlendFeature> {
    // Check for torus (fillet)
    if max_blend_radius > 0.0 {
        if let Some(torus) = face_torus(brep, si, shi, fi) {
            // Torus minor radius indicates fillet radius
            if torus.minor_radius > 0.0 && torus.minor_radius <= max_blend_radius {
                let face = &brep.solids[si].shells[shi].faces[fi];
                let sample_point = get_face_sample_point(brep, si, shi, fi).unwrap_or(torus.center);
                return Some(BlendFeature {
                    face_indices: vec![fi],
                    is_fillet: true,
                    radius: torus.minor_radius,
                    chamfer_distance: 0.0,
                    sample_point,
                    normal: face.normal.normalize_or_zero(),
                });
            }
        }

        // Check for sphere (ball-end fillet)
        if let Some(sphere) = face_sphere(brep, si, shi, fi)
            && sphere.radius > 0.0 && sphere.radius <= max_blend_radius {
                let face = &brep.solids[si].shells[shi].faces[fi];
                return Some(BlendFeature {
                    face_indices: vec![fi],
                    is_fillet: true,
                    radius: sphere.radius,
                    chamfer_distance: 0.0,
                    sample_point: sphere.center,
                    normal: face.normal.normalize_or_zero(),
                });
            }

        // Check for cylinder with small radius (edge fillet)
        if let Some(cyl) = face_cylinder(brep, si, shi, fi)
            && cyl.radius > 0.0 && cyl.radius <= max_blend_radius {
                let face = &brep.solids[si].shells[shi].faces[fi];
                let sample_point = cyl.origin;
                return Some(BlendFeature {
                    face_indices: vec![fi],
                    is_fillet: true,
                    radius: cyl.radius,
                    chamfer_distance: 0.0,
                    sample_point,
                    normal: face.normal.normalize_or_zero(),
                });
            }
    }

    // Check for chamfer (small planar face connecting two other faces at an angle)
    if max_chamfer_distance > 0.0
        && let Some(_plane) = face_plane(brep, si, shi, fi) {
            // Heuristic: small planar face that connects two non-parallel faces
            // This is a simplified check; a full implementation would analyze
            // the adjacent faces and check if they meet at an angle
            let face = &brep.solids[si].shells[shi].faces[fi];

            // Estimate chamfer size from face dimensions
            let face_area = estimate_face_area(brep, si, shi, fi);
            let chamfer_estimate = face_area.sqrt() / 1.414; // Approximate for 45-degree chamfer

            if chamfer_estimate > 0.0 && chamfer_estimate <= max_chamfer_distance {
                let sample_point = get_face_sample_point(brep, si, shi, fi).unwrap_or_default();
                return Some(BlendFeature {
                    face_indices: vec![fi],
                    is_fillet: false,
                    radius: 0.0,
                    chamfer_distance: chamfer_estimate,
                    sample_point,
                    normal: face.normal.normalize_or_zero(),
                });
            }
        }

    None
}

/// Get a sample point from a face (first vertex).
fn get_face_sample_point(brep: &BRep, si: usize, shi: usize, fi: usize) -> Option<DVec3> {
    let face = &brep.solids[si].shells[shi].faces[fi];
    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx)
            && let Some(v) = brep.vertices.get(edge.start) {
                return Some(v.point);
            }
    }
    None
}

/// Estimate the area of a face using fan triangulation.
fn estimate_face_area(brep: &BRep, si: usize, shi: usize, fi: usize) -> f64 {
    let face = &brep.solids[si].shells[shi].faces[fi];

    // Collect vertex positions in order
    let mut pts: Vec<DVec3> = Vec::new();
    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let vi = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vi) {
                pts.push(v.point);
            }
        }
    }

    if pts.len() < 3 {
        return 0.0;
    }

    // Fan triangulation from first point
    let p0 = pts[0];
    let mut area = 0.0f64;
    for i in 1..pts.len() - 1 {
        area += (pts[i] - p0).cross(pts[i + 1] - p0).length() * 0.5;
    }

    area
}

/// Return `true` if two normalized axis vectors are parallel (or antiparallel).
fn axes_parallel(a1: DVec3, a2: DVec3) -> bool {
    a1.normalize_or_zero()
        .cross(a2.normalize_or_zero())
        .length()
        < AXIS_PARALLEL_TOL
}

/// Return `true` if two infinite axis lines (origin + direction) are the same
/// line in 3-D space.
fn axes_same_line(o1: DVec3, ax1: DVec3, o2: DVec3, ax2: DVec3) -> bool {
    if !axes_parallel(ax1, ax2) {
        return false;
    }
    let ax = ax1.normalize_or_zero();
    let d = o2 - o1;
    let dist_sq = (d - d.dot(ax) * ax).length_squared();
    dist_sq < AXIS_PARALLEL_TOL * AXIS_PARALLEL_TOL
}

/// Determine whether a cylindrical face is likely a hole wall by checking
/// the majority voting of `face.normal` against the radial outward directions
/// at each boundary vertex.
///
/// **Limitation**: after a boolean operation the stored `face.normal` may be
/// the cylinder's seam direction rather than the true outward-from-solid normal
/// (this is a known limitation of the legacy curved-face split path).  We use a
/// majority vote across ALL boundary vertices to reduce sensitivity to any single
/// seam-direction artifact.  Falls back to `true` (hole) on tie or missing data.
fn is_hole_face(face: &Face, brep: &BRep, cyl: &CylindricalSurface) -> bool {
    let ax = cyl.axis.normalize_or_zero();
    let face_n = face.normal.normalize_or_zero();
    if face_n.length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
        return true; // no normal stored -> assume hole
    }

    // Collect unique vertex indices to avoid biasing the vote on seam
    // vertices that appear as both `edge.end` and `edge.start` on adjacent
    // edges in the outer wire.
    let mut seen: HashSet<usize> = HashSet::new();
    let mut collect_verts = |wire: &Wire| {
        for we in &wire.edges {
            let Some(edge) = brep.edges.get(we.idx) else { continue; };
            seen.insert(edge.start);
            seen.insert(edge.end);
        }
    };
    collect_verts(&face.outer_wire);
    for iw in &face.inner_wires {
        collect_verts(iw);
    }

    let mut hole_votes: i32 = 0;
    let mut boss_votes: i32 = 0;
    for &vi in &seen {
        let Some(v) = brep.vertices.get(vi) else { continue; };
        let to_pt = v.point - cyl.origin;
        let radial = to_pt - to_pt.dot(ax) * ax;
        if radial.length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
            continue;
        }
        let radial_dir = radial.normalize();
        // For a cylindrical HOLE wall (drill removed from solid), the
        // boolean builder stores face.normal pointing OUTWARD from the
        // cylinder axis (i.e. in the +radial direction), because the
        // cylinder's seam normal is the ref_dir == the outward direction
        // at the seam.  dot > 0 -> face_n agrees with outward radial -> hole.
        // For a BOSS, the face is part of the exterior of the added cylinder,
        // so the same seam normal convention still holds: dot > 0 -> hole wall.
        // We therefore use dot > 0 as the "outward (hole)" signal.
        let dot = face_n.dot(radial_dir);
        if dot > TOLERANCE_MESH_LEGACY {
            hole_votes += 1;
        } else if dot < -TOLERANCE_MESH_LEGACY {
            boss_votes += 1;
        }
    }

    // Tie or majority hole votes -> assume hole.
    hole_votes >= boss_votes
}

/// Compute the min/max projection of all wall-face vertices onto the cylinder axis.
fn axis_extent_of_group(
    brep: &BRep,
    si: usize,
    shi: usize,
    face_indices: &[usize],
    cyl: &CylindricalSurface,
) -> (f64, f64) {
    let ax = cyl.axis.normalize_or_zero();
    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;

    for &fi in face_indices {
        let face = &brep.solids[si].shells[shi].faces[fi];
        for we in &face.outer_wire.edges {
            let Some(edge) = brep.edges.get(we.idx) else {
                continue;
            };
            for &vi in &[edge.start, edge.end] {
                let Some(v) = brep.vertices.get(vi) else {
                    continue;
                };
                let t = (v.point - cyl.origin).dot(ax);
                if t < t_min {
                    t_min = t;
                }
                if t > t_max {
                    t_max = t;
                }
            }
        }
    }

    if t_min.is_infinite() {
        (0.0, 0.0)
    } else {
        (t_min, t_max)
    }
}

// -- Public detection functions ---------------------------------------------

/// Detect all cylindrical features (holes and bosses) in `solids[0].shells[0]`
/// whose radius falls within the specified bounds.
///
/// Pass `max_hole_radius = 0.0` to skip hole detection, and similarly for
/// `max_boss_radius`.
///
/// Returns a list of [`CylindricalFeature`] objects, one per connected group.
pub fn detect_cylindrical_features(
    brep: &BRep,
    max_hole_radius: f64,
    max_boss_radius: f64,
) -> Vec<CylindricalFeature> {
    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> [local face_idx] adjacency.
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face is a cylinder surface.
        let Some(cyl) = face_cylinder(brep, si, shi, start) else {
            continue;
        };

        // Use the larger of the two thresholds so we collect the full
        // group without pre-filtering on the (unreliable) is_hole flag.
        let effective_max = max_hole_radius.max(max_boss_radius);
        if effective_max <= 0.0 || cyl.radius > effective_max {
            // Radius out of range -> skip but don't mark visited so other
            // features sharing an edge can still be explored.
            continue;
        }

        // BFS: collect all connected cylindrical faces on the same axis/radius.
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(fi) = queue.pop_front() {
            group.push(fi);

            let face_edges: Vec<usize> = {
                let f = &shell.faces[fi];
                let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                for iw in &f.inner_wires {
                    es.extend(iw.edges.iter().map(|we| we.idx));
                }
                es
            };

            for ei in face_edges {
                let Some(neighbours) = edge_to_faces.get(&ei) else {
                    continue;
                };
                for &nfi in neighbours {
                    if visited[nfi] {
                        continue;
                    }
                    let Some(ncyl) = face_cylinder(brep, si, shi, nfi) else {
                        continue;
                    };
                    if (ncyl.radius - cyl.radius).abs() > RADIUS_TOL {
                        continue;
                    }
                    if !axes_same_line(cyl.origin, cyl.axis, ncyl.origin, ncyl.axis) {
                        continue;
                    }
                    visited[nfi] = true;
                    queue.push_back(nfi);
                }
            }
        }

        // Determine is_hole by group-level majority vote.  Aggregating across
        // all faces in the group avoids sensitivity to the seam-direction
        // artefact that makes per-face voting unreliable after a boolean op.
        // Tie-breaks towards hole (the more common defeaturing target).
        let group_hole_count = group
            .iter()
            .filter(|&&fi| is_hole_face(&shell.faces[fi], brep, &cyl))
            .count();
        let is_hole = group_hole_count * 2 >= group.len();

        let (t_min, t_max) = axis_extent_of_group(brep, si, shi, &group, &cyl);

        features.push(CylindricalFeature {
            face_indices: group,
            is_hole,
            origin: cyl.origin,
            axis: cyl.axis.normalize_or_zero(),
            radius: cyl.radius,
            t_min,
            t_max,
        });
    }

    features
}

/// Detect all conical features (holes and bosses) in `solids[0].shells[0]`
/// whose reference radius falls within the specified bounds.
///
/// Pass `max_hole_radius = 0.0` to skip hole detection.
///
/// Returns a list of [`ConicalFeature`] objects, one per connected group.
pub fn detect_conical_features(
    brep: &BRep,
    max_hole_radius: f64,
) -> Vec<ConicalFeature> {
    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> [local face_idx] adjacency.
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face is a cone surface.
        let Some(cone) = face_cone(brep, si, shi, start) else {
            continue;
        };

        // Calculate reference radius at mid-height for size filtering.
        // Use a point on the cone surface to estimate the reference radius.
        let face = &shell.faces[start];
        let mut sample_point: Option<DVec3> = None;
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx)
                && let Some(v) = brep.vertices.get(edge.start) {
                    sample_point = Some(v.point);
                    break;
                }
        }

        let reference_radius = if let Some(pt) = sample_point {
            let ax = cone.axis.normalize_or_zero();
            let to_pt = pt - cone.apex;
            let t = to_pt.dot(ax);
            // Radius at height t from apex: r = t * tan(half_angle)
            t.abs() * cone.half_angle_rad.tan()
        } else {
            // Fallback: use the cone's stored reference radius if available
            cone.radius
        };

        if max_hole_radius <= 0.0 || reference_radius > max_hole_radius {
            continue;
        }

        // BFS: collect all connected conical faces on the same axis/apex.
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(fi) = queue.pop_front() {
            group.push(fi);

            let face_edges: Vec<usize> = {
                let f = &shell.faces[fi];
                let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                for iw in &f.inner_wires {
                    es.extend(iw.edges.iter().map(|we| we.idx));
                }
                es
            };

            for ei in face_edges {
                let Some(neighbours) = edge_to_faces.get(&ei) else {
                    continue;
                };
                for &nfi in neighbours {
                    if visited[nfi] {
                        continue;
                    }
                    let Some(ncone) = face_cone(brep, si, shi, nfi) else {
                        continue;
                    };
                    // Check same axis line and similar half angle.
                    if !axes_same_line(cone.apex, cone.axis, ncone.apex, ncone.axis) {
                        continue;
                    }
                    if (ncone.half_angle_rad - cone.half_angle_rad).abs() > TOLERANCE_RETRY_LADDER_COARSE {
                        continue;
                    }
                    visited[nfi] = true;
                    queue.push_back(nfi);
                }
            }
        }

        // Determine is_hole by checking if the cone widens away from apex
        // and the face normal points inward (toward axis).
        let ax = cone.axis.normalize_or_zero();
        let is_hole = {
            // Sample the first face's normal to determine hole vs boss.
            let f = &shell.faces[group[0]];
            let fnormal = f.normal.normalize_or_zero();
            // For a conical hole, the normal points outward from the solid,
            // which for an inward-facing cone wall means pointing toward the axis.
            // Check by computing radial direction at a sample point.
            let mut toward_axis_votes = 0i32;
            let mut away_axis_votes = 0i32;
            for we in &f.outer_wire.edges {
                if let Some(edge) = brep.edges.get(we.idx) {
                    for &vi in &[edge.start, edge.end] {
                        if let Some(v) = brep.vertices.get(vi) {
                            let to_pt = v.point - cone.apex;
                            let radial = to_pt - to_pt.dot(ax) * ax;
                            if radial.length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
                                continue;
                            }
                            let radial_dir = radial.normalize();
                            // Dot < 0 means normal points toward axis (hole).
                            let dot = fnormal.dot(radial_dir);
                            if dot < -TOLERANCE_MESH_LEGACY {
                                toward_axis_votes += 1;
                            } else if dot > TOLERANCE_MESH_LEGACY {
                                away_axis_votes += 1;
                            }
                        }
                    }
                }
            }
            toward_axis_votes >= away_axis_votes
        };

        // Compute axis extents from vertices.
        let mut t_min = f64::INFINITY;
        let mut t_max = f64::NEG_INFINITY;
        for &fi in &group {
            let f = &shell.faces[fi];
            for we in &f.outer_wire.edges {
                if let Some(edge) = brep.edges.get(we.idx) {
                    for &vi in &[edge.start, edge.end] {
                        if let Some(v) = brep.vertices.get(vi) {
                            let t = (v.point - cone.apex).dot(ax);
                            t_min = t_min.min(t);
                            t_max = t_max.max(t);
                        }
                    }
                }
            }
        }

        if t_min.is_infinite() {
            t_min = 0.0;
            t_max = 0.0;
        }

        features.push(ConicalFeature {
            face_indices: group,
            is_hole,
            apex: cone.apex,
            axis: ax,
            reference_radius,
            half_angle: cone.half_angle_rad,
            t_min,
            t_max,
        });
    }

    features
}

/// Detect all slot features in `solids[0].shells[0]`.
///
/// Slots are elongated features with rectangular or rounded cross-sections,
/// typically formed by a combination of planar and cylindrical faces.
///
/// Parameters:
/// - `max_width`: Maximum slot width to consider
/// - `max_depth`: Maximum slot depth to consider
///
/// Returns a list of [`SlotFeature`] objects.
pub fn detect_slot_features(
    brep: &BRep,
    max_width: f64,
    max_depth: f64,
) -> Vec<SlotFeature> {
    if max_width <= 0.0 || max_depth <= 0.0 {
        return Vec::new();
    }

    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    // Strategy: Find groups of connected planar faces that form a slot-like shape
    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face could be part of a slot (planar or small-radius cylinder)
        let is_slot_candidate = face_plane(brep, si, shi, start).is_some()
            || face_cylinder(brep, si, shi, start)
                .map(|c| c.radius <= max_width)
                .unwrap_or(false);

        if !is_slot_candidate {
            continue;
        }

        // BFS to find connected slot-like faces
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        // Collect geometry information
        let mut planes: Vec<Plane> = Vec::new();
        let mut cylinders: Vec<(CylindricalSurface, usize)> = Vec::new();

        while let Some(fi) = queue.pop_front() {
            group.push(fi);

            if let Some(plane) = face_plane(brep, si, shi, fi) {
                planes.push(plane);
            }
            if let Some(cyl) = face_cylinder(brep, si, shi, fi)
                && cyl.radius <= max_width {
                    cylinders.push((cyl, fi));
                }

            let face_edges: Vec<usize> = {
                let f = &shell.faces[fi];
                let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                for iw in &f.inner_wires {
                    es.extend(iw.edges.iter().map(|we| we.idx));
                }
                es
            };

            for ei in face_edges {
                let Some(neighbours) = edge_to_faces.get(&ei) else {
                    continue;
                };
                for &nfi in neighbours {
                    if visited[nfi] {
                        continue;
                    }
                    // Check if neighbor is also slot-like
                    let is_neighbor_candidate = face_plane(brep, si, shi, nfi).is_some()
                        || face_cylinder(brep, si, shi, nfi)
                            .map(|c| c.radius <= max_width)
                            .unwrap_or(false);

                    if is_neighbor_candidate {
                        visited[nfi] = true;
                        queue.push_back(nfi);
                    }
                }
            }
        }

        // Analyze the group to determine if it's a slot
        if group.len() < 3 {
            // A slot needs at least a bottom and two sides
            continue;
        }

        // Try to identify slot geometry
        if let Some(slot) = analyze_slot_group(brep, si, shi, &group, &planes, &cylinders, max_width, max_depth) {
            features.push(slot);
        }
    }

    features
}

/// Analyze a group of faces to determine if they form a slot.
fn analyze_slot_group(
    brep: &BRep,
    si: usize,
    shi: usize,
    group: &[usize],
    planes: &[Plane],
    cylinders: &[(CylindricalSurface, usize)],
    max_width: f64,
    max_depth: f64,
) -> Option<SlotFeature> {
    // Need at least one planar face (bottom)
    if planes.is_empty() {
        return None;
    }

    // Collect all vertices
    let mut vertices: Vec<DVec3> = Vec::new();
    for &fi in group {
        let face = &brep.solids[si].shells[shi].faces[fi];
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                if let Some(v) = brep.vertices.get(edge.start) {
                    vertices.push(v.point);
                }
                if let Some(v) = brep.vertices.get(edge.end) {
                    vertices.push(v.point);
                }
            }
        }
    }

    if vertices.len() < 4 {
        return None;
    }

    // Compute bounding box
    let mut min_pt = vertices[0];
    let mut max_pt = vertices[0];
    for pt in &vertices[1..] {
        min_pt = min_pt.min(*pt);
        max_pt = max_pt.max(*pt);
    }

    let dimensions = max_pt - min_pt;
    let mut dims = [dimensions.x, dimensions.y, dimensions.z];
    dims.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

    // Slot should be elongated (one dimension significantly larger than width)
    let length = dims[0];
    let width = dims[1];
    let depth = dims[2];

    if width > max_width || depth > max_depth {
        return None;
    }

    // Determine slot orientation
    let length_dir = if (dimensions.x - length).abs() < TOLERANCE_MESH_LEGACY {
        DVec3::X
    } else if (dimensions.y - length).abs() < TOLERANCE_MESH_LEGACY {
        DVec3::Y
    } else {
        DVec3::Z
    };

    let width_dir = if (dimensions.x - width).abs() < TOLERANCE_MESH_LEGACY {
        DVec3::X
    } else if (dimensions.y - width).abs() < TOLERANCE_MESH_LEGACY {
        DVec3::Y
    } else {
        DVec3::Z
    };

    let depth_dir = if (dimensions.x - depth).abs() < TOLERANCE_MESH_LEGACY {
        DVec3::X
    } else if (dimensions.y - depth).abs() < TOLERANCE_MESH_LEGACY {
        DVec3::Y
    } else {
        DVec3::Z
    };

    // Check for rounded ends (cylindrical faces at slot ends)
    let has_rounded_ends = !cylinders.is_empty();

    let center = (min_pt + max_pt) * 0.5;
    let origin = center - depth_dir * depth * 0.5; // Bottom center

    Some(SlotFeature {
        face_indices: group.to_vec(),
        is_recess: true, // Assume recess by default
        length,
        width,
        depth,
        origin,
        length_dir,
        width_dir,
        depth_dir,
        has_rounded_ends,
    })
}

/// Detect all pocket features in `solids[0].shells[0]`.
///
/// Pockets are enclosed recesses with flat bottoms and side walls.
/// Both circular and rectangular pockets are detected.
///
/// Parameters:
/// - `max_diameter`: Maximum pocket diameter (or max dimension) to consider
/// - `max_depth`: Maximum pocket depth to consider
///
/// Returns a list of [`PocketFeature`] objects.
pub fn detect_pocket_features(
    brep: &BRep,
    max_diameter: f64,
    max_depth: f64,
) -> Vec<PocketFeature> {
    if max_diameter <= 0.0 || max_depth <= 0.0 {
        return Vec::new();
    }

    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    // Find connected groups that could be pockets
    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face could be part of a pocket
        let is_pocket_candidate = face_plane(brep, si, shi, start).is_some()
            || face_cylinder(brep, si, shi, start)
                .map(|c| c.radius <= max_diameter)
                .unwrap_or(false);

        if !is_pocket_candidate {
            continue;
        }

        // BFS to find connected pocket-like faces
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        let mut has_cylindrical_walls = false;
        let mut cylindrical_radius = 0.0f64;
        let mut wall_planes: Vec<Plane> = Vec::new();

        while let Some(fi) = queue.pop_front() {
            group.push(fi);

            if let Some(plane) = face_plane(brep, si, shi, fi) {
                wall_planes.push(plane);
            }
            if let Some(cyl) = face_cylinder(brep, si, shi, fi)
                && cyl.radius <= max_diameter {
                    has_cylindrical_walls = true;
                    cylindrical_radius = cyl.radius;
                }

            let face_edges: Vec<usize> = {
                let f = &shell.faces[fi];
                let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                for iw in &f.inner_wires {
                    es.extend(iw.edges.iter().map(|we| we.idx));
                }
                es
            };

            for ei in face_edges {
                let Some(neighbours) = edge_to_faces.get(&ei) else {
                    continue;
                };
                for &nfi in neighbours {
                    if visited[nfi] {
                        continue;
                    }
                    let is_neighbor_candidate = face_plane(brep, si, shi, nfi).is_some()
                        || face_cylinder(brep, si, shi, nfi)
                            .map(|c| c.radius <= max_diameter)
                            .unwrap_or(false);

                    if is_neighbor_candidate {
                        visited[nfi] = true;
                        queue.push_back(nfi);
                    }
                }
            }
        }

        // Analyze the group
        if let Some(pocket) = analyze_pocket_group(
            brep,
            si,
            shi,
            &group,
            has_cylindrical_walls,
            cylindrical_radius,
            &wall_planes,
            max_diameter,
            max_depth,
        ) {
            features.push(pocket);
        }
    }

    features
}

/// Analyze a group of faces to determine if they form a pocket.
fn analyze_pocket_group(
    brep: &BRep,
    si: usize,
    shi: usize,
    group: &[usize],
    has_cylindrical_walls: bool,
    cylindrical_radius: f64,
    wall_planes: &[Plane],
    max_diameter: f64,
    max_depth: f64,
) -> Option<PocketFeature> {
    if group.is_empty() {
        return None;
    }

    // Collect all vertices
    let mut vertices: Vec<DVec3> = Vec::new();
    for &fi in group {
        let face = &brep.solids[si].shells[shi].faces[fi];
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                if let Some(v) = brep.vertices.get(edge.start) {
                    vertices.push(v.point);
                }
                if let Some(v) = brep.vertices.get(edge.end) {
                    vertices.push(v.point);
                }
            }
        }
    }

    if vertices.len() < 3 {
        return None;
    }

    // Compute bounding box
    let mut min_pt = vertices[0];
    let mut max_pt = vertices[0];
    for pt in &vertices[1..] {
        min_pt = min_pt.min(*pt);
        max_pt = max_pt.max(*pt);
    }

    let dimensions = max_pt - min_pt;
    let mut dims = [dimensions.x, dimensions.y, dimensions.z];
    dims.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

    let depth = dims[2]; // Smallest dimension is likely depth

    if depth > max_depth {
        return None;
    }

    let center = (min_pt + max_pt) * 0.5;

    // Determine if circular or rectangular
    let is_circular = has_cylindrical_walls && cylindrical_radius > 0.0;

    let (diameter, width, length) = if is_circular {
        (cylindrical_radius * 2.0, 0.0, 0.0)
    } else {
        (dims[0], dims[1], dims[0])
    };

    if diameter > max_diameter {
        return None;
    }

    // Compute approximate normal from wall planes
    let normal = wall_planes
        .first()
        .map(|p| p.normal.normalize_or_zero())
        .unwrap_or(DVec3::Z);

    Some(PocketFeature {
        face_indices: group.to_vec(),
        is_recess: true,
        diameter,
        depth,
        center,
        normal,
        is_circular,
        width,
        length,
        is_through: false, // Will be determined by enhanced detection
        bottom_face_index: None,
        wall_face_indices: Vec::new(),
    })
}
include!("e1.rs");
include!("e2.rs");
include!("tests_inc.rs");include!("test2_inc.rs");

