use glam::DVec3;
use serde::{Deserialize, Serialize};

/// Geometric (analytic) model types: position, curve, surface, primitive descriptors.
///
/// This module describes *what shape is*.
pub mod geom;

/// Assemblies: hierarchy, instancing, and world-transform flattening.
///
/// Analogous to the shape hierarchy managed by OCCT `XCAFDoc_ShapeTool`.
pub mod assembly;

/// Topology model types: vertex/edge/face/shell/solid incidence relationships.
///
/// This module describes *how things are connected*.
pub mod topology;

/// Shape properties: surface area, volume, centroid.
///
/// Analogous to OCCT `GProp_GProps` + `BRepGProp`.
pub mod properties;
pub use crate::properties::face_flat_iter;

/// Topology query helpers: edge adjacency, vertex adjacency, shape counts.
///
/// Analogous to OCCT `TopExp_Explorer` and `TopExp::MapShapesAndAncestors`.
pub mod topo_query;
/// Topology simplification helpers: wire/edge cleanup for fragmented topology.
pub mod topo_simplify;

/// Cached graph-topology wrapper with O(1) adjacency, DFS/BFS traversal,
/// and mutation-dirty tracking.
///
/// Analogous to OCCT `BRepGraph` module (new in OCCT 7.7+).
pub mod brep_graph;

/// Persistent naming hooks for stable user-level topology labels.
///
/// Analogous to OCCT OCAF/TopoNaming-style name tables.
pub mod naming;

/// Persistent naming semantics for BRepGraph topology entities.
///
/// Provides stable, operation-surviving identifiers for topology entities.
pub mod persistent_naming;

/// Differential geometry: principal curvatures, Gaussian curvature, mean curvature.
///
/// Analogous to OCCT `GeomLProp_SLProps`.
pub mod curvature;

/// Curve arc-length computation.
///
/// Analogous to OCCT `GCPnts_AbscissaPoint` / `CPnts_AbscissaPoint::Length`.
pub mod arc_length;

/// Visual appearance: per-face/solid RGB color and basic material.
///
/// Analogous to OCCT `XCAFDoc_ColorTool`.
pub mod appearance;

/// Dimension and tolerance object model (GD&T).
///
/// Analogous to OCCT `XCAFDimTolObjects`.
pub mod dim_tol;

/// Annotation object model for CAD annotations (PMI).
///
/// Analogous to OCCT `XCAFNoteObjects`.
pub mod annotation;

/// Precision constants and per-entity tolerance query helpers.
///
/// Analogous to OCCT `Precision` class and `BRep_Tool::Tolerance`.
pub mod topods;
pub mod tolerance;

/// Curve fitting: B-spline interpolation and approximation through point sets.
///
/// Analogous to OCCT `GeomAPI_Interpolate` and `GeomAPI_PointsToBSpline`.
pub mod fit;

/// Closest-point projection from a 3D point onto a curve or surface.
///
/// Analogous to OCCT `GeomAPI_ProjectPointOnCurve` and
/// `GeomAPI_ProjectPointOnSurf`.
pub mod projection;

/// Shape-to-shape and point-to-shape minimum distance.
///
/// Analogous to OCCT `BRepExtrema_DistShapeShape`.
pub mod distance;

/// Curve-curve extrema: find (s,t) minimising |C1(s) − C2(t)|.
///
/// Analogous to OCCT `GeomAPI_ExtremaCurveCurve`.
pub mod extrema;

/// NURBS interoperability: convert analytic curves/surfaces to BSpline.
///
/// Analogous to OCCT `GeomConvert::CurveToBSplineCurve` /
/// `GeomConvert::SurfaceToBSplineSurface`.
pub mod nurbs_convert;

/// Standard collections analogous to OCCT TColStd package.
pub mod tcol_std;

/// Mathematical utilities analogous to OCCT TKMath package.
pub mod math_utils;

/// Design-feature array/pattern creation utilities.
pub mod array;

/// Curve and surface trimming and extension.
///
/// Analogous to OCCT `Geom_TrimmedCurve` construction helpers,
/// `GeomAPI_ExtendCurveToPoint`, and `Geom_RectangularTrimmedSurface`.
pub mod extend;

pub use distance::{ShapeDistance, min_distance, point_to_shape_distance};
pub use extend::{
 CurveEnd, SurfaceBoundary, extend_bspline_surface, extend_curve_by_length,
 extend_curve_to_point, insert_knot_to_multiplicity, insert_knot_u_once, insert_knot_v_once,
 refine_bspline_surface_isoparametric_spans, trim_curve, trim_surface,
};
pub use extrema::{CurveCurveExtrema, ExtremaPair, extrema_curve_curve};
pub use fit::{FitError, approximate_points, interpolate_points, interpolate_points_2d};
pub use nurbs_convert::{
 bezier_curve_to_bspline, bezier_surface_to_bspline, circle_to_bspline, curve_to_bspline,
 cylinder_to_bspline, ellipse_to_bspline, line_to_bspline, line_to_bspline_range,
 plane_to_bspline, sphere_to_bspline, surface_to_bspline,
};
pub use projection::{
 CurveProjection, SurfaceProjection, closest_point_on_curve, closest_point_on_surface,
 closest_point_on_surface_near, make_pcurve_on_surface,
};

pub use appearance::{Color, FaceColor, StepColor};
pub use dim_tol::{
 DimensionType, GeometricToleranceType, ToleranceModifier,
 DimensionalTolerance, GeometricToleranceObject, DatumReference, DatumSystem, DimTolStore,
};
pub use annotation::{
 NoteType, ArrowType, WeldType,
 AnnotationNote, TextAnnotation, LeaderLine, SurfaceTextureSymbol,
 WeldSymbol, BalloonAnnotation, AnnotationStore,
 Annotation, AnnotationKind, Note, NoteCategory, NoteTarget, View, ViewProjection,
};
pub use arc_length::arc_length;
pub use curvature::{gaussian_curvature, mean_curvature, principal_curvatures};
pub use geom::{Point3, Vec3, Point2, Vec2};
pub use geom::PrimitiveSolid;
pub use geom::TrimmedSurface;
pub use geom::{
 ArchimedeanSpiral2d, BSplineCurve2, CircleInvolute2d, Ellipse2d, LogarithmicSpiral2d,
 SineWave2d,
};
pub use geom::{BSplineSurface, CoonsSurface, LinearExtrusionSurface, RevolutionSurface, RuledSurface};
pub use geom::{BezierCurve2, BezierCurve3, BezierSurface, TriBezierSurface};
pub use geom::{Curve2d, Curve3, Surface3};
pub use geom::{Curve2dEval, CurveEval, SurfaceEval, any_perpendicular};
pub use geom::{EllipsoidalSurface, HelicoidSurface, PipeSurface};
pub use geom::{CircularHelix3, Hyperbola3, Parabola3, SineWave3};
pub use geom::{OffsetCurve3, OffsetSurface};
pub use properties::{
 InertiaTensor, centroid, face_surface_area, face_triangles_pub, try_analytic_face_surface_area_pub, point_in_spherical_polygon_3d_pub, inertia_tensor, signed_volume, surface_area, volume,
};
pub use tolerance::{
 ANGULAR, APPROXIMATION, brep_same_parameter, CONFUSION, edge_same_parameter, edge_same_range, edge_tolerance,
 face_domain, face_tolerance, model_tolerance, vertex_tolerance,
 resize_tolerance_arrays,
 set_vertex_tolerance, update_vertex_tolerance,
 set_edge_tolerance, update_edge_tolerance,
 set_face_tolerance, update_face_tolerance,
 finalize_tolerance_hierarchy,
 step_export_uncertainty,
};
pub use topo_query::{
 vertex_storage_len,
 salient_vertex_indices,
 edge_adjacent_faces, edge_count, face_count, face_edges, is_degenerate_edge,
 periodic_seam_edge_indices, wire_edges_unique_by_index,
 topological_vertex_count, vertex_adjacent_edges, vertex_indices,
};
pub use topo_simplify::{merge_collinear_edges_in_wires, merge_collinear_brep_edges};
pub use brep_graph::{
 BfsFaces, BRepGraph, BRepGraphBuilder, BRepGraphCheckpointData,
 BRepGraphTool, DfsEdgesFromVertex, DfsFaces,
 ManifoldRepairHints, NonManifoldSummary, RepairHint,
};
pub use naming::{PersistentNamingHooks, TopoEntityRef};
pub use persistent_naming::{
 ConflictResolution, CrossOperationHistory, CrossOperationStabilityReport, EntityGenealogy,
 EntityType, EntityTypeStability, IssueSeverity, NamingConflictResolution, NamingContext,
 NamingEvent, NamingHistory, NamingIssue, NamingRule, NamingStabilityReport,
 NamePropagationPolicy, OperationId, OperationRecord, OperationStats, OperationType,
 PersistentId, PersistentNamingEngine, PersistentNamingHooksExt,
};
pub use topology::{Compound, CompSolid, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

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

/// Main BRep type — alias for the OCCT-aligned `topods::BRep`.
pub use topods::BRep;

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
 pub fn surface_bounding_box(surface: &geom::Surface3, vertices: &[crate::Vertex]) -> Option<[DVec3; 2]> {
 match surface {
 geom::Surface3::Cylinder(cyl) => {
 let r = cyl.radius;
 let axis = cyl.axis.normalize_or_zero();
 if axis.length_squared() < 0.5 { return None; }

 // Project all vertices onto the cylinder axis to find the
 // axial extent of the trimmed surface.
 let mut min_axial = f64::INFINITY;
 let mut max_axial = f64::NEG_INFINITY;
 for v in vertices {
 let proj = (v.point - cyl.origin).dot(axis);
 min_axial = min_axial.min(proj);
 max_axial = max_axial.max(proj);
 }
 if !min_axial.is_finite() { return None; }

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
 if !min_axial.is_finite() { return None; }

 let max_r = cone.radius_at_axial(min_axial)
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

#[cfg(test)]
pub mod tkmath_gtests;

#[cfg(test)]
pub mod math_gtests;

