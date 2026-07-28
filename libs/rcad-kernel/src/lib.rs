use glam::DVec3;
use serde::{Deserialize, Serialize};

// ============================================================================
// Module tree — aligned with OCCT toolkit hierarchy
// ============================================================================

// Module-level re-exports for backward compatibility.
// Downstream crates import these as `rcad_kernel::topods`, `rcad_kernel::topology`, etc.
pub use topo::topods;
pub use topo::topology;
pub use topo::topo_query;
pub use topo::topo_shape;
pub use topo::brep_graph;
pub use core::precision;
pub use core::precision as tolerance;
pub use ocaf::persistent_naming;
pub use ocaf::naming;
pub use ocaf::appearance;
pub use ocaf::annotation;
pub use ocaf::dim_tol;
pub use core::units;
pub use core::message;
pub use math::projection;
pub use math::fit;
pub use math::math_utils;
pub use math::curvature;
pub use math::arc_length;
pub use base::extend;
pub use base::nurbs_convert;

/// Geometric (analytic) model types: position, curve, surface, primitive descriptors.
///
/// This module describes *what shape is*.
/// Corresponds to OCCT TKG2d + TKG3d.
pub mod geom;

/// Core TKernel-level types: precision, collections, naming, appearance, GD&T, annotations.
pub mod core;

/// Topology model types (TopoDS/ModelingData): vertex/edge/face/shell/solid incidence.
pub mod topo;

/// Math algorithms (TKMath): arc length, curvature, distance, properties, projection, fitting, extrema.
pub mod math;

/// Geometry foundation (TKGeomBase): NURBS conversion, curve/surface extension/trimming.
pub mod base;

/// Assemblies: hierarchy, instancing, and world-transform flattening.
///
/// Analogous to the shape hierarchy managed by OCCT `XCAFDoc_ShapeTool`.
pub mod ocaf;

// ============================================================================
// Re-exports from core (TKernel)
// ============================================================================

pub use core::precision::{
    ANGULAR, APPROXIMATION, COMPUTATIONAL, CONFUSION, INFINITE_VALUE, INTERSECTION,
    PCONFUSION, SQUARE_COMPUTATIONAL, SQUARE_CONFUSION, SQUARE_INTERSECTION,
    brep_same_parameter, edge_same_parameter,
    edge_same_range, edge_tolerance, face_domain, face_tolerance, finalize_tolerance_hierarchy,
    is_infinite_value, is_negative_infinite_value, is_positive_infinite_value, model_tolerance,
    p_approximation, p_approximation_with_tangent,
    p_confusion, p_confusion_with_tangent,
    p_intersection, p_intersection_with_tangent,
    parametric, parametric_default,
    resize_tolerance_arrays, set_edge_tolerance, set_face_tolerance, set_vertex_tolerance,
    square_p_confusion,
    step_export_uncertainty, update_edge_tolerance, update_face_tolerance, update_vertex_tolerance,
    vertex_tolerance,
};

pub use ocaf::appearance::{Color, FaceColor, StepColor};

pub use ocaf::annotation::{
    Annotation, AnnotationKind, AnnotationNote, AnnotationStore, ArrowType, BalloonAnnotation,
    LeaderLine, Note, NoteCategory, NoteTarget, NoteType, SurfaceTextureSymbol, TextAnnotation,
    View, ViewProjection, WeldSymbol, WeldType,
};

pub use ocaf::dim_tol::{
    DatumReference, DatumSystem, DimTolStore, DimensionType, DimensionalTolerance,
    GeometricToleranceObject, GeometricToleranceType, ToleranceModifier,
};

pub use ocaf::naming::{PersistentNamingHooks, TopoEntityRef};

pub use ocaf::persistent_naming::{
    ConflictResolution, CrossOperationHistory, CrossOperationStabilityReport, EntityGenealogy,
    EntityType, EntityTypeStability, IssueSeverity, NamePropagationPolicy,
    NamingConflictResolution, NamingContext, NamingEvent, NamingHistory, NamingIssue, NamingRule,
    NamingStabilityReport, OperationId, OperationRecord, OperationStats, OperationType,
    PersistentId, PersistentNamingEngine, PersistentNamingHooksExt,
};

// ============================================================================
// Re-exports from topo (TopoDS / ModelingData)
// ============================================================================

pub use topo::topology::{CompSolid, Compound, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

pub use topo::topo_query::{
    edge_adjacent_faces, edge_count, face_count, face_edges, is_degenerate_edge,
    periodic_seam_edge_indices, salient_vertex_indices, topological_vertex_count,
    vertex_adjacent_edges, vertex_indices, vertex_storage_len, wire_edges_unique_by_index,
};

pub use topo::topo_simplify::{merge_collinear_brep_edges, merge_collinear_edges_in_wires};

pub use topo::brep_graph::{
    BRepGraph, BRepGraphBuilder, BRepGraphCheckpointData, BRepGraphTool, BfsFaces,
    DfsEdgesFromVertex, DfsFaces, ManifoldRepairHints, NonManifoldSummary, RepairHint,
};

// ============================================================================
// Re-exports from math (TKMath)
// ============================================================================

pub use math::arc_length::arc_length;

pub use math::curvature::{gaussian_curvature, mean_curvature, principal_curvatures};

pub use math::distance::{ShapeDistance, min_distance, point_to_shape_distance};

pub use math::extrema::{CurveCurveExtrema, ExtremaPair, extrema_curve_curve};

pub use math::fit::{FitError, approximate_points, interpolate_points, interpolate_points_2d};

pub use math::projection::{
    CurveProjection, SurfaceProjection, closest_point_on_curve, closest_point_on_curve_range,
    closest_point_on_surface, closest_point_on_surface_near, make_pcurve_on_surface,
};

pub use math::properties::{
    InertiaTensor, centroid, face_surface_area, face_triangles_pub, inertia_tensor,
    point_in_spherical_polygon_3d_pub, signed_volume, surface_area,
    try_analytic_face_surface_area_pub, volume,
};
pub use math::properties::face_flat_iter;

// ============================================================================
// Re-exports from base (TKGeomBase)
// ============================================================================

pub use base::nurbs_convert::{
    bezier_curve_to_bspline, bezier_surface_to_bspline, circle_to_bspline, curve_to_bspline,
    cylinder_to_bspline, ellipse_to_bspline, line_to_bspline, line_to_bspline_range,
    plane_to_bspline, sphere_to_bspline, surface_to_bspline,
};

pub use base::extend::{
    CurveEnd, SurfaceBoundary, extend_bspline_surface, extend_curve_by_length,
    extend_curve_to_point, insert_knot_to_multiplicity, insert_knot_u_once, insert_knot_v_once,
    refine_bspline_surface_isoparametric_spans, trim_curve, trim_surface,
};

// ============================================================================
// Re-exports from geom (TKG2d + TKG3d) — kept alongside crate root for backward compat
// ============================================================================

pub use geom::PrimitiveSolid;
pub use geom::TrimmedSurface;
pub use geom::{
    ArchimedeanSpiral2d, BSplineCurve2, CircleInvolute2d, Ellipse2d, LogarithmicSpiral2d,
    SineWave2d,
};
pub use geom::{
    BSplineSurface, CoonsSurface, LinearExtrusionSurface, RevolutionSurface, RuledSurface,
};
pub use geom::{BezierCurve2, BezierCurve3, BezierSurface, TriBezierSurface};
pub use geom::{CircularHelix3, Hyperbola3, Parabola3, SineWave3};
pub use geom::{Curve2d, Curve3, Surface3};
pub use geom::{Curve2dEval, CurveEval, SurfaceEval, any_perpendicular};
pub use geom::{EllipsoidalSurface, HelicoidSurface, PipeSurface};
pub use geom::{OffsetCurve3, OffsetSurface};
pub use geom::{Point2, Point3, Vec2, Vec3};

// ============================================================================
// Top-level types defined in this crate (no OCCT toolkit mapping)
// ============================================================================

/// A parameter-space curve binding that ties a 3D edge to an adjacent face's
/// surface parameter domain (u, v).  Analogous to OCCT `BRep_CurveOnSurface`.
///
/// `surface_idx` indexes into `GeomStore.surfaces`.
/// `curve2d_idx` indexes into `GeomStore.curve2ds`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PCurve {
    pub surface_idx: usize,
    pub curve2d_idx: usize,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GeomStore {
    /// Pool of 3D analytic curves.
    pub curves: Vec<Curve3>,
    /// Pool of analytic surfaces.
    pub surfaces: Vec<Surface3>,
    /// Pool of 2D parameter-space curves used by PCurves.
    pub curve2ds: Vec<Curve2d>,
    /// Indexed by `BRep.edges` index; value is index into `curves`.
    pub edge_curve: Vec<Option<usize>>,
    /// Flattened face order across solids/shells; value is index into `surfaces`.
    pub face_surface: Vec<Option<usize>>,
    /// Indexed by `BRep.edges` index; each entry is the list of PCurves for
    /// that edge on its adjacent faces (usually 1, seam edges have 2).
    pub edge_pcurves: Vec<Vec<PCurve>>,
    /// Parallel to `edge_curve`: the parameter range [t1, t2] of the edge on its
    /// 3D curve. `None` = unknown (algorithms fall back to `CurveEval::default_domain`).
    /// Analogous to `BRep_Edge::Range()` in OCCT.
    #[serde(default)]
    pub edge_curve_range: Vec<Option<[f64; 2]>>,
    /// Parallel to `BRep.edges`: `true` if this is a degenerate edge (zero-length,
    /// e.g. a polar singularity). Analogous to `BRep_Edge::Degenerated()` in OCCT.
    #[serde(default)]
    pub edge_degenerated: Vec<bool>,
    /// Per-vertex tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to `BRep.vertices`. Analogous to `BRep_Tool::Tolerance(vertex)` in OCCT.
    #[serde(default)]
    pub vertex_tolerance: Vec<f64>,
    /// Per-edge tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to `BRep.edges`. Analogous to `BRep_Tool::Tolerance(edge)` in OCCT.
    #[serde(default)]
    pub edge_tolerance: Vec<f64>,
    /// Per-face tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to the flattened face order (same indexing as `face_surface`).
    /// Analogous to `BRep_Tool::Tolerance(face)` in OCCT.
    #[serde(default)]
    pub face_tolerance: Vec<f64>,
    /// Per-curve2d parameter range [t1, t2].
    ///
    /// Used when the PCurve originates from a STEP `TRIMMED_CURVE` entity in
    /// 2D parameter space. `None` means the natural domain of the curve is used.
    /// Parallel to `GeomStore.curve2ds`. Analogous to `edge_curve_range` for 3D.
    #[serde(default)]
    pub curve2d_range: Vec<Option<[f64; 2]>>,
    /// Per-face surface parameter domain override [u1, u2, v1, v2].
    ///
    /// When populated (e.g. from a STEP `RECTANGULAR_TRIMMED_SURFACE`), the face
    /// is restricted to this subdomain of its underlying surface. `None` means
    /// `SurfaceEval::default_domain()` is used. Parallel to `face_surface`.
    /// Analogous to `edge_curve_range` for 3D curves.
    #[serde(default)]
    pub face_surface_range: Vec<Option<[f64; 4]>>,
    /// Per-edge SameParameter flag.
    ///
    /// `true` if the 3D curve and all PCurves share the same parameterization
    /// (i.e. the parameter `t` on the 3D curve maps directly to the same `t`
    /// on every PCurve). Analogous to `BRep_Edge::SameParameter()` in OCCT.
    /// When absent or empty, assumed `true` for analytic primitives we generate.
    #[serde(default)]
    pub edge_same_parameter: Vec<bool>,
    /// Per-edge SameRange flag.
    ///
    /// `true` if all PCurves on this edge share the same `[t1, t2]` parameter
    /// range as the 3D curve's `edge_curve_range`. Analogous to
    /// `BRep_Edge::SameRange()` in OCCT.
    /// When absent or empty, assumed `true` for analytic primitives we generate.
    #[serde(default)]
    pub edge_same_range: Vec<bool>,
    /// ✅ OCCT : (FillInternalVertices  )
    /// face_surface (  flat face index)。
    #[serde(default)]
    pub face_internal_vertices: Vec<Vec<usize>>,
    /// Per-edge start/end vertex parameters on the edge's 3D curve.
    /// [start_param, end_param] — parallel to BRep.edges.
    #[serde(default)]
    pub edge_vertex_params: Vec<Option<[f64; 2]>>,
}

impl GeomStore {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Main BRep type — alias for the OCCT-aligned `topo::topods::BRep`.
pub use topo::topods::BRep;

/// Conservative bounding-box contribution from an analytic curve.
pub fn curve_bounding_box(curve: &geom::Curve3) -> Option<[DVec3; 2]> {
    match curve {
        geom::Curve3::Circle(c) => {
            // A circle only extends in its own plane, not in the normal
            // direction.  For each axis i the max extent is
            // radius * sqrt(1 - n_i^2) where n is the unit normal.
            let n = c.normal.normalize();
            let extent = DVec3::new(
                c.radius * (1.0 - n.x * n.x).sqrt(),
                c.radius * (1.0 - n.y * n.y).sqrt(),
                c.radius * (1.0 - n.z * n.z).sqrt(),
            );
            Some([c.center - extent, c.center + extent])
        }
        geom::Curve3::Ellipse(e) => {
            // Same plane-restricted extent for ellipses.
            let n = e.normal.normalize();
            let max_r = e.major_radius.max(e.minor_radius);
            let extent = DVec3::new(
                max_r * (1.0 - n.x * n.x).sqrt(),
                max_r * (1.0 - n.y * n.y).sqrt(),
                max_r * (1.0 - n.z * n.z).sqrt(),
            );
            Some([e.center - extent, e.center + extent])
        }
        geom::Curve3::BSpline(b) => {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for &p in &b.control_points {
                mn = mn.min(p);
                mx = mx.max(p);
            }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        geom::Curve3::Bezier(b) => {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for &p in &b.control_points {
                mn = mn.min(p);
                mx = mx.max(p);
            }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        _ => None,
    }
}

/// Conservative bounding-box contribution from an analytic surface,
/// expanding based on vertex positions projected onto the surface frame.
pub fn surface_bounding_box(
    surface: &geom::Surface3,
    vertices: &[crate::topo::topology::Vertex],
) -> Option<[DVec3; 2]> {
    match surface {
        geom::Surface3::Cylinder(cyl) => {
            let r = cyl.radius;
            let axis = cyl.axis.normalize_or_zero();
            if axis.length_squared() < 0.5 {
                return None;
            }

            // Project all vertices onto the cylinder axis to find the
            // axial extent of the trimmed surface.
            let mut min_axial = f64::INFINITY;
            let mut max_axial = f64::NEG_INFINITY;
            for v in vertices {
                let proj = (v.point - cyl.origin).dot(axis);
                min_axial = min_axial.min(proj);
                max_axial = max_axial.max(proj);
            }
            if !min_axial.is_finite() {
                return None;
            }

            let p_lo = cyl.origin + axis * min_axial;
            let p_hi = cyl.origin + axis * max_axial;
            let rv = DVec3::splat(r);
            Some([p_lo.min(p_hi) - rv, p_lo.max(p_hi) + rv])
        }
        geom::Surface3::Sphere(sph) => {
            let r = DVec3::splat(sph.radius);
            Some([sph.center - r, sph.center + r])
        }
        geom::Surface3::Torus(tor) => {
            let r = tor.major_radius + tor.minor_radius;
            let rv = DVec3::splat(r);
            Some([tor.center - rv, tor.center + rv])
        }
        geom::Surface3::Cone(cone) => {
            let axis = cone.axis_dir();
            let apex = cone.apex_point();

            // Project vertices onto the cone axis.
            let mut min_axial = f64::INFINITY;
            let mut max_axial = f64::NEG_INFINITY;
            for v in vertices {
                let proj = (v.point - apex).dot(axis);
                min_axial = min_axial.min(proj);
                max_axial = max_axial.max(proj);
            }
            if !min_axial.is_finite() {
                return None;
            }

            let max_r = cone
                .radius_at_axial(min_axial)
                .max(cone.radius_at_axial(max_axial));
            let rv = DVec3::splat(max_r.max(cone.radius));
            let p_lo = apex + axis * min_axial;
            let p_hi = apex + axis * max_axial;
            Some([p_lo.min(p_hi) - rv, p_lo.max(p_hi) + rv])
        }
        geom::Surface3::Ellipsoid(e) => {
            let max_r = e.radius_x.max(e.radius_y).max(e.radius_z);
            let rv = DVec3::splat(max_r);
            Some([e.center - rv, e.center + rv])
        }
        geom::Surface3::BSpline(b) => {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for row in &b.control_points {
                for p in row {
                    mn = mn.min(*p);
                    mx = mx.max(*p);
                }
            }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        geom::Surface3::Bezier(b) => {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for row in &b.control_points {
                for &p in row {
                    mn = mn.min(p);
                    mx = mx.max(p);
                }
            }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        _ => None,
    }
}
