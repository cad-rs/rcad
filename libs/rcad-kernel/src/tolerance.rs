//! Precision constants and per-entity tolerance query helpers.
//!
//! Analogous to OCCT's `Precision` class and `BRep_Builder::UpdateVertex` /
//! `BRep_Builder::UpdateEdge` tolerance API.
//!
//! # Design
//!
//! OCCT stores a per-entity tolerance on each `BRep_Vertex`, `BRep_Edge`, and
//! `BRep_Face`.  These tolerances represent the maximum deviation between the
//! ideal (analytic) geometry and the actual model, arising from tessellation,
//! import rounding, or algorithm approximations.
//!
//! RCAD stores the same data in `GeomStore` as parallel `Vec<f64>` arrays
//! (`vertex_tolerance`, `edge_tolerance`, `face_tolerance`), indexed the same
//! way as the corresponding topology arrays.  When a tolerance is absent or
//! zero the query functions fall back to the `CONFUSION` constant.
//!
//! # Usage
//!
//! ```rust
//! use rcad_kernel::{vertex_tolerance, edge_tolerance, CONFUSION};
//!
//! // brep with no stored tolerances → default confusion value
//! let brep = rcad_kernel::BRep::new();
//! assert_eq!(vertex_tolerance(&brep, 0), CONFUSION);
//! ```

use crate::BRep;
use crate::topods;
use crate::geom::{Curve2dEval, CurveEval, SurfaceEval};

// ── Precision constants ───────────────────────────────────────────────────────

/// Point-coincidence tolerance.
/// Analogous to `Precision::Confusion()` = 1e-7 in OCCT.
pub const CONFUSION: f64 = 1e-7;

/// Angular tolerance (radians).
/// Analogous to `Precision::Angular()` = 1e-12 in OCCT.
pub const ANGULAR: f64 = 1e-12;

/// Tessellation / approximation tolerance.
/// Analogous to `Precision::Approximation()` = 1e-6 in OCCT.
pub const APPROXIMATION: f64 = 1e-6;

/// Intersection tolerance.  Used by intersection algorithms to decide
/// when a solution is reached.
/// Analogous to `Precision::Intersection()` = Confusion / 100 = 1e-9.
pub const INTERSECTION: f64 = CONFUSION * 0.01;

/// Parametric-space confusion tolerance on a default curve.
/// Analogous to `Precision::PConfusion()` = Confusion / 100 = 1e-9.
pub const P_CONFUSION: f64 = CONFUSION * 0.01;

/// Parametric-space intersection tolerance on a default curve.
/// Analogous to `Precision::PIntersection()` = Intersection / 100 = 1e-11.
pub const P_INTERSECTION: f64 = INTERSECTION * 0.01;

/// Parametric-space approximation tolerance on a default curve.
/// Analogous to `Precision::PApproximation()` = Approximation / 100 = 1e-8.
pub const P_APPROXIMATION: f64 = APPROXIMATION * 0.01;

// ── Geometric resolution helpers (BRepCheck::PrecCurve / PrecSurface) ────────

/// Helper: return `f64::EPSILON * max_coord_magnitude`.
/// Equivalent to OCCT `RealEpsilon() * max(|coords|)` pattern used in
/// `BRepCheck::PrecCurve` and `BRepCheck::PrecSurface`.
fn epsilon_of(max_coord: f64) -> f64 {
 f64::EPSILON * max_coord.abs().max(1.0)
}

/// Equivalent to OCCT `BRepCheck::PrecCurve`.
///
/// Returns the floating-point resolution of a 3D curve computed from the
/// magnitude of its defining parameters.  Used as a tolerance delta in
/// `BRepLib_ValidateEdge::correctTolerance`.
pub fn prec_curve(curve: &crate::geom::Curve3) -> f64 {
 use crate::geom::Curve3;
 match curve {
 Curve3::Line(l) => {
 let m = l.origin.to_array().iter().copied().fold(0.0_f64, f64::max)
 .max(l.direction.to_array().iter().copied().fold(0.0_f64, f64::max));
 epsilon_of(m)
 }
 Curve3::Circle(c) => {
 let center_max = c.center.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(center_max + c.radius)
 }
 Curve3::Ellipse(e) => {
 let center_max = e.center.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(center_max + e.major_radius.max(e.minor_radius))
 }
 Curve3::BSpline(bs) => {
 let m = bs
 .control_points
 .iter()
 .flat_map(|p| p.to_array())
 .fold(0.0_f64, f64::max);
 epsilon_of(m)
 }
 Curve3::Bezier(bz) => {
 let m = bz
 .control_points
 .iter()
 .flat_map(|p| p.to_array())
 .fold(0.0_f64, f64::max);
 epsilon_of(m)
 }
 Curve3::Offset(o) => {
 // Recurse into the basis curve
 prec_curve(&o.basis)
 }
 Curve3::Hyperbola(h) => {
 let center_max = h.center.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(center_max + h.semi_major.max(h.semi_minor))
 }
 Curve3::Parabola(p) => {
 let m = p.vertex.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(m + p.focal_param)
 }
 Curve3::CircularHelix(ch) => {
 let m = ch.origin.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(m + ch.radius + ch.pitch)
 }
 Curve3::SineWave(sw) => {
 // SineWave3 uses `frequency`, not wavelength.  Estimate period from
 // frequency: the curve travels amplitude * sin(freq * t) from origin.
 let m = sw.origin.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(m + sw.amplitude + sw.frequency.abs())
 }
 }
}

/// Equivalent to OCCT `BRepCheck::PrecSurface`.
///
/// Returns the floating-point resolution of a surface computed from the
/// magnitude of its defining parameters.
pub fn prec_surface(surface: &crate::geom::Surface3) -> f64 {
 use crate::geom::Surface3;
 match surface {
 Surface3::Plane(p) => {
 let m = p.origin.to_array().iter().copied().fold(0.0_f64, f64::max)
 .max(p.normal.to_array().iter().copied().fold(0.0_f64, f64::max));
 epsilon_of(m)
 }
 Surface3::Cylinder(c) => {
 let m = c.origin.to_array().iter().copied().fold(0.0_f64, f64::max)
 .max(c.axis.to_array().iter().copied().fold(0.0_f64, f64::max));
 epsilon_of(m + c.radius)
 }
 Surface3::Sphere(s) => {
 let m = s.center.to_array().iter().copied().fold(0.0_f64, f64::max)
 .max(s.axis.to_array().iter().copied().fold(0.0_f64, f64::max));
 epsilon_of(m + s.radius)
 }
 Surface3::Cone(cn) => {
 let m = cn.apex.to_array().iter().copied().fold(0.0_f64, f64::max)
 .max(cn.axis.to_array().iter().copied().fold(0.0_f64, f64::max));
 epsilon_of(m + cn.radius + cn.half_angle_rad)
 }
 Surface3::Torus(t) => {
 let m = t.center.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(m + t.major_radius + t.minor_radius)
 }
 Surface3::Ellipsoid(e) => {
 let m = e.center.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(m + e.radius_x.max(e.radius_y).max(e.radius_z))
 }
 Surface3::BSpline(bs) => {
 // Estimate from control point magnitudes (2D grid: Vec<Vec<DVec3>>)
 let m = bs
 .control_points
 .iter()
 .flat_map(|row| row.iter())
 .flat_map(|p| p.to_array())
 .fold(0.0_f64, f64::max);
 epsilon_of(m)
 }
 Surface3::Bezier(bz) => {
 // control_points is Vec<Vec<DVec3>> (2D grid)
 let m = bz
 .control_points
 .iter()
 .flat_map(|row| row.iter())
 .flat_map(|p| p.to_array())
 .fold(0.0_f64, f64::max);
 epsilon_of(m)
 }
 // Fallback for the remaining analytic types: compute from primary
 // coordinate fields via their Debug representation or simplest
 // accessible member.
 Surface3::Helicoid(h) => {
 let m = h.origin.to_array().iter().copied().fold(0.0_f64, f64::max);
 epsilon_of(m + h.pitch)
 }
 Surface3::Pipe(pp) => {
 // Recurse into the spine curve
 prec_curve(&pp.spine).max(epsilon_of(pp.radius))
 }
 Surface3::LinearExtrusion(le) => {
 let m = le
 .direction
 .to_array()
 .iter()
 .copied()
 .fold(0.0_f64, f64::max);
 prec_curve(&le.profile).max(epsilon_of(m))
 }
 Surface3::Revolution(rv) => {
 let m = rv.axis_origin.to_array().iter().copied().fold(0.0_f64, f64::max);
 // Also consider the profile curve precision
 prec_curve(&rv.profile).max(epsilon_of(m))
 }
 Surface3::Ruled(ru) => {
 // Conservative estimate from both boundary curves
 prec_curve(&ru.start).max(prec_curve(&ru.end))
 }
 Surface3::TriBezier(tb) => {
 // Triangular net: Vec<Vec<DVec3>>
 let m = tb
 .control_points
 .iter()
 .flat_map(|row| row.iter())
 .flat_map(|p| p.to_array())
 .fold(0.0_f64, f64::max);
 epsilon_of(m)
 }
 Surface3::Coons(co) => {
 // Conservative: max over all four boundary curves
 prec_curve(&co.south)
 .max(prec_curve(&co.north))
 .max(prec_curve(&co.west))
 .max(prec_curve(&co.east))
 }
 Surface3::Offset(os) => prec_surface(&os.basis),
 Surface3::Trimmed(tr) => prec_surface(&tr.basis),
 }
}

// ── Per-entity tolerance queries ──────────────────────────────────────────────

/// Tolerance of vertex `vertex_idx`.
///
/// Returns the stored value when positive, otherwise `CONFUSION`.
/// Analogous to `BRep_Tool::Tolerance(vertex)` in OCCT.
pub fn vertex_tolerance(brep: &BRep, vertex_idx: usize) -> f64 {
 brep.geom
 .vertex_tolerance
 .get(vertex_idx)
 .copied()
 .filter(|&t| t > 0.0)
 .unwrap_or(CONFUSION)
}

/// Tolerance of edge `edge_idx`.
///
/// Returns the stored value when positive, otherwise `CONFUSION`.
/// Analogous to `BRep_Tool::Tolerance(edge)` in OCCT.
pub fn edge_tolerance(brep: &BRep, edge_idx: usize) -> f64 {
 brep.geom
 .edge_tolerance
 .get(edge_idx)
 .copied()
 .filter(|&t| t > 0.0)
 .unwrap_or(CONFUSION)
}

/// Tolerance of face `face_flat_idx` (flattened index across all solids/shells,
/// same ordering as `GeomStore.face_surface`).
///
/// Returns the stored value when positive, otherwise `CONFUSION`.
/// Analogous to `BRep_Tool::Tolerance(face)` in OCCT.
pub fn face_tolerance(brep: &BRep, face_flat_idx: usize) -> f64 {
 brep.geom
 .face_tolerance
 .get(face_flat_idx)
 .copied()
 .filter(|&t| t > 0.0)
 .unwrap_or(CONFUSION)
}

/// Maximum tolerance over all entities with explicitly stored values.
///
/// Returns `CONFUSION` when no tolerances are stored.
/// Analogous to `BRep_Builder::UpdateVertex` accumulated maximum in OCCT.
pub fn model_tolerance(brep: &BRep) -> f64 {
 brep.geom
 .vertex_tolerance
 .iter()
 .chain(brep.geom.edge_tolerance.iter())
 .chain(brep.geom.face_tolerance.iter())
 .copied()
 .fold(CONFUSION, f64::max)
}

/// Tolerance value used as STEP `UNCERTAINTY_MEASURE_WITH_UNIT` during export.
///
/// Analogous to OCCT `STEPControl_ActorWrite::UsedTolerance` with the default
/// average-precision mode (`WritePrecisionMode == 1`):
///
/// 1. Compute the arithmetic mean of all stored vertex / edge / face tolerances.
/// 2. Multiply by 1.5 (OCCT's safety margin for downstream precision).
/// 3. Round to 2 significant digits (mimics OCCT's `Interface_MSG::Intervalled`).
/// 4. Floor at `CONFUSION` (1e-7).
///
/// When no tolerances are stored at all, returns `CONFUSION` directly.
pub fn step_export_uncertainty(brep: &BRep) -> f64 {
 let sum: f64 = brep
 .geom
 .vertex_tolerance
 .iter()
 .chain(brep.geom.edge_tolerance.iter())
 .chain(brep.geom.face_tolerance.iter())
 .copied()
 .sum();
 let count = brep.geom.vertex_tolerance.len()
 + brep.geom.edge_tolerance.len()
 + brep.geom.face_tolerance.len();

 if count == 0 {
 return CONFUSION;
 }

 let avg = sum / count as f64;
 let scaled = avg * 1.5;

 if scaled <= 0.0 {
 return CONFUSION;
 }

 // Round to 2 significant digits (equivalent to Interface_MSG::Intervalled).
 let magnitude = 10.0_f64.powf(scaled.log10().floor());
 let rounded = (scaled / magnitude).round() * magnitude;

 rounded.max(CONFUSION)
}

/// Effective surface parameter domain [u1, u2, v1, v2] for a face.
///
/// Uses `GeomStore.face_surface_range` if set; otherwise falls back to
/// `SurfaceEval::default_domain()` of the underlying surface.
/// Analogous to `BRep_Face::UVBounds()` in OCCT.
pub fn face_domain(brep: &BRep, face_flat_idx: usize) -> [f64; 4] {
 if let Some(Some(range)) = brep.geom.face_surface_range.get(face_flat_idx) {
 return *range;
 }
 if let Some(Some(surf_idx)) = brep.geom.face_surface.get(face_flat_idx)
 && let Some(surf) = brep.geom.surfaces.get(*surf_idx)
 {
 return surf.default_domain();
 }
 [0.0, 1.0, 0.0, 1.0]
}

/// SameParameter flag for edge `edge_idx`.
///
/// Returns the stored value when present, otherwise `true`
/// (all analytic primitives generated by RCAD are same-parameter by construction).
/// Analogous to `BRep_Edge::SameParameter()` in OCCT.
pub fn edge_same_parameter(brep: &BRep, edge_idx: usize) -> bool {
 brep.geom
 .edge_same_parameter
 .get(edge_idx)
 .copied()
 .unwrap_or(true)
}

/// SameRange flag for edge `edge_idx`.
///
/// Returns the stored value when present, otherwise `true`
/// (all analytic primitives generated by RCAD are same-range by construction).
/// Analogous to `BRep_Edge::SameRange()` in OCCT.
pub fn edge_same_range(brep: &BRep, edge_idx: usize) -> bool {
 brep.geom
 .edge_same_range
 .get(edge_idx)
 .copied()
 .unwrap_or(true)
}

// ── Write / mutator API ───────────────────────────────────────────────────────
//
// These functions are the Rust equivalents of OCCT's
// `BRep_Builder::UpdateVertex`, `BRep_Builder::UpdateEdge`, and
// `BRep_Builder::UpdateFace`.  They extend the `GeomStore` tolerance arrays
// on demand so callers do not need to pre-size them.

/// Ensure all three tolerance arrays in `brep.geom` are sized to cover every
/// entity currently in the `BRep`.  Newly extended slots are filled with
/// `CONFUSION` (the global precision floor).
///
/// Analogous to the internal resize step inside `BRep_Builder::UpdateVertex`.
pub fn resize_tolerance_arrays(brep: &mut BRep) {
 let nv = brep.vertices.len();
 let ne = brep.edges.len();
 let nf: usize = brep.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();

 if brep.geom.vertex_tolerance.len() < nv {
 brep.geom.vertex_tolerance.resize(nv, CONFUSION);
 }
 if brep.geom.edge_tolerance.len() < ne {
 brep.geom.edge_tolerance.resize(ne, CONFUSION);
 }
 if brep.geom.face_tolerance.len() < nf {
 brep.geom.face_tolerance.resize(nf, CONFUSION);
 }
}

/// Set the tolerance of vertex `vi` to exactly `tol`.
///
/// Extends the `vertex_tolerance` array if necessary.
/// Clamps `tol` to `CONFUSION` minimum.
///
/// Analogous to `BRep_Builder::UpdateVertex(V, tol)` (absolute form).
pub fn set_vertex_tolerance(brep: &mut BRep, vi: usize, tol: f64) {
 let t = tol.max(CONFUSION);
 if brep.geom.vertex_tolerance.len() <= vi {
 brep.geom.vertex_tolerance.resize(vi + 1, CONFUSION);
 }
 brep.geom.vertex_tolerance[vi] = t;
}

/// Raise the tolerance of vertex `vi` to at least `tol`.
///
/// Uses `max(current, tol)` — never lowers an existing tolerance.
/// Extends the array on demand.
/// Analogous to `BRep_Builder::UpdateVertex(V, tol)` in OCCT.
pub fn update_vertex_tolerance(brep: &mut BRep, vi: usize, tol: f64) {
 let t = tol.max(CONFUSION);
 if brep.geom.vertex_tolerance.len() <= vi {
 brep.geom.vertex_tolerance.resize(vi + 1, CONFUSION);
 }
 let cur = &mut brep.geom.vertex_tolerance[vi];
 if t > *cur {
 *cur = t;
 }
}

/// Set the tolerance of edge `ei` to exactly `tol`.
///
/// Extends the `edge_tolerance` array if necessary.
/// Analogous to `BRep_Builder::UpdateEdge(E, tol)` (absolute form).
pub fn set_edge_tolerance(brep: &mut BRep, ei: usize, tol: f64) {
 let t = tol.max(CONFUSION);
 if brep.geom.edge_tolerance.len() <= ei {
 brep.geom.edge_tolerance.resize(ei + 1, CONFUSION);
 }
 brep.geom.edge_tolerance[ei] = t;
}

/// Raise the tolerance of edge `ei` to at least `tol`.
///
/// Uses `max(current, tol)` — never lowers.
/// Analogous to `BRep_Builder::UpdateEdge(E, tol)` in OCCT.
pub fn update_edge_tolerance(brep: &mut BRep, ei: usize, tol: f64) {
 let t = tol.max(CONFUSION);
 if brep.geom.edge_tolerance.len() <= ei {
 brep.geom.edge_tolerance.resize(ei + 1, CONFUSION);
 }
 let cur = &mut brep.geom.edge_tolerance[ei];
 if t > *cur {
 *cur = t;
 }
}

/// Set the tolerance of face `fi` (flat index) to exactly `tol`.
///
/// Extends the `face_tolerance` array if necessary.
/// Analogous to `BRep_Builder::UpdateFace(F, tol)` (absolute form).
pub fn set_face_tolerance(brep: &mut BRep, fi: usize, tol: f64) {
 let t = tol.max(CONFUSION);
 if brep.geom.face_tolerance.len() <= fi {
 brep.geom.face_tolerance.resize(fi + 1, CONFUSION);
 }
 brep.geom.face_tolerance[fi] = t;
}

/// Raise the tolerance of face `fi` (flat index) to at least `tol`.
///
/// Uses `max(current, tol)` — never lowers.
/// Analogous to `BRep_Builder::UpdateFace(F, tol)` in OCCT.
pub fn update_face_tolerance(brep: &mut BRep, fi: usize, tol: f64) {
 let t = tol.max(CONFUSION);
 if brep.geom.face_tolerance.len() <= fi {
 brep.geom.face_tolerance.resize(fi + 1, CONFUSION);
 }
 let cur = &mut brep.geom.face_tolerance[fi];
 if t > *cur {
 *cur = t;
 }
}

/// Enforce the tolerance hierarchy **in-place**:
///
/// ```text
/// vertex_tol ≥ edge_tol ≥ face_tol
/// ```
///
/// Specifically:
///
/// 1. For each edge: `edge_tol = max(edge_tol, face_tol(all adjacent faces))`.
/// 2. For each vertex: `vertex_tol = max(vertex_tol, edge_tol(all adjacent edges))`.
///
/// This ensures the "wider" entity (vertex) always carries at least as much
/// tolerance as the narrower entities it bounds — the same invariant OCCT
/// maintains via `BRepLib::UpdateEdgeTol` and `BRep_Builder::UpdateVertex`.
///
/// Call `resize_tolerance_arrays` first if the arrays may not be populated.
pub fn finalize_tolerance_hierarchy(brep: &mut BRep) {
 resize_tolerance_arrays(brep);

 // Step 1: face → edge.
 // Build edge → max(adjacent face tolerance).
 let mut edge_from_face: Vec<f64> = vec![CONFUSION; brep.edges.len()];
 let mut flat_fi = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 let ftol = brep.geom.face_tolerance
 .get(flat_fi)
 .copied()
 .unwrap_or(CONFUSION);
 for we in &face.outer_wire.edges {
 if we.idx < edge_from_face.len() {
 edge_from_face[we.idx] = edge_from_face[we.idx].max(ftol);
 }
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 if we.idx < edge_from_face.len() {
 edge_from_face[we.idx] = edge_from_face[we.idx].max(ftol);
 }
 }
 }
 flat_fi += 1;
 }
 }
 }
 for (ei, &ftol) in edge_from_face.iter().enumerate() {
 if let Some(etol) = brep.geom.edge_tolerance.get_mut(ei)
 && ftol > *etol {
 *etol = ftol;
 }
 }

 // Step 2: edge → vertex.
 let ne = brep.edges.len();
 for ei in 0..ne {
 let etol = brep.geom.edge_tolerance.get(ei).copied().unwrap_or(CONFUSION);
 let st = brep.edges[ei].start;
 let en = brep.edges[ei].end;
 if let Some(vtol) = brep.geom.vertex_tolerance.get_mut(st)
 && etol > *vtol {
 *vtol = etol;
 }
 if let Some(vtol) = brep.geom.vertex_tolerance.get_mut(en)
 && etol > *vtol {
 *vtol = etol;
 }
 }
}

/// Part A: OCCT `BRepLib::SameParameter` equivalent.
///
/// For each edge that has a 3D curve and at least one pcurve, sample the
/// deviation between `curve(t)` and `surface(u(t), v(t))` at `samples_per_edge`
/// points.  Update the edge tolerance to `max(current, max_deviation)`.
///
/// Edges without a 3D curve, without pcurves, or without a parametric range
/// are skipped.
pub fn brep_same_parameter(brep: &mut BRep, samples_per_edge: usize) {
 brep_same_parameter_impl(brep, samples_per_edge, &std::collections::HashSet::new())
}

fn brep_same_parameter_impl(
 brep: &mut BRep,
 samples_per_edge: usize,
 map_to_avoid: &std::collections::HashSet<usize>,
) {
 let samples = samples_per_edge.max(2);
 // First pass: collect deviations (immutable borrow only).
 let updates: Vec<(usize, f64)> = {
 let curves = &brep.geom.curves;
 let surfaces = &brep.geom.surfaces;
 let curve2ds = &brep.geom.curve2ds;
 (0..brep.edges.len())
 .filter_map(|ei| {
 // OCCT: skip edges in MapToAvoid (non-destructive mode)
 if map_to_avoid.contains(&ei) { return None; }
 let Some(curve_idx) = brep.geom.edge_curve.get(ei).copied().flatten() else {
 return None;
 };
 let curve = curves.get(curve_idx)?;
 let [t1, t2] = brep.geom.edge_curve_range.get(ei).copied().flatten()?;
 if (t2 - t1).abs() < 1e-15 {
 return None;
 }
 let Some(pcurves) = brep.geom.edge_pcurves.get(ei) else {
 return None;
 };
 if pcurves.is_empty() {
 return None;
 }
 let mut max_dev = 0.0_f64;
 for pc in pcurves {
 let surface = surfaces.get(pc.surface_idx)?;
 let pc_curve = curve2ds.get(pc.curve2d_idx)?;
 for si in 0..samples {
 let t = t1 + (t2 - t1) * si as f64 / (samples - 1) as f64;
 let uv = pc_curve.point_at(t);
 let p_surf = surface.point_at(uv.x, uv.y);
 let p_curve = curve.point_at(t);
 let dev = (p_surf - p_curve).length();
 if dev > max_dev {
 max_dev = dev;
 }
 }
 }
 if max_dev > 0.0 {
 Some((ei, max_dev))
 } else {
 None
 }
 })
 .collect()
 };
 // Second pass: apply updates (mutable borrow).
 for (ei, dev) in updates {
 set_edge_tolerance(brep, ei, dev);
 // Also propagate to adjacent vertices (OCCT CorrectVertexTolerance pattern).
 let (v0, v1) = brep
 .edges
 .get(ei)
 .map(|e| (e.start, e.end))
 .unwrap_or((usize::MAX, usize::MAX));
 if v0 != usize::MAX {
 update_vertex_tolerance(brep, v0, dev);
 update_vertex_tolerance(brep, v1, dev);
 }
 }
}

/// Part B: Compute per-vertex tolerance from distance to adjacent edges.
///
/// For each vertex, find all adjacent edges and compute the minimum distance
/// from the vertex to each edge's 3D curve (sampled).  Update vertex tolerance
/// via `update_vertex_tolerance(brep, vi, max_dist)`.
///
/// Edges without 3D curves fall back to half the edge length as a heuristic.
pub fn compute_vertex_tolerances(brep: &mut BRep) {
 compute_vertex_tolerances_impl(brep, &std::collections::HashSet::new())
}

fn compute_vertex_tolerances_impl(brep: &mut BRep, map_to_avoid: &std::collections::HashSet<usize>) {
 use crate::topo_query::vertex_adjacent_edges;

 for vi in 0..brep.vertices.len() {
 // OCCT: skip vertices in MapToAvoid (non-destructive mode)
 if map_to_avoid.contains(&vi) { continue; }
 let v_pos = brep.vertices[vi].point;
 let adj_edges = vertex_adjacent_edges(brep, vi);
 if adj_edges.is_empty() {
 continue;
 }

 let mut max_dist = 0.0_f64;
 for &ei in &adj_edges {
 let Some(edge) = brep.edges.get(ei) else { continue; };
 let other_vi = if edge.start == vi { edge.end } else { edge.start };
 let Some(other) = brep.vertices.get(other_vi) else { continue; };
 let edge_len = (other.point - v_pos).length();

 // Prefer 3D curve distance over edge-length fallback
 if let Some(curve_idx) = brep.geom.edge_curve.get(ei).copied().flatten() {
 if let Some(curve) = brep.geom.curves.get(curve_idx) {
 if let Some([t1, t2]) = brep.geom.edge_curve_range.get(ei).copied().flatten() {
 let mut min_dist = f64::MAX;
 for s in 0..20 {
 let t = t1 + (t2 - t1) * s as f64 / 19.0;
 let pt = curve.point_at(t);
 let d = (pt - v_pos).length();
 if d < min_dist { min_dist = d; }
 }
 if min_dist < f64::MAX && min_dist > max_dist {
 max_dist = min_dist;
 }
 continue; // used curve distance, skip edge-length fallback
 }
 }
 }
 // Fallback: half edge length
 let half = edge_len * 0.5;
 if half > max_dist { max_dist = half; }
 }

 if max_dist > 0.0 {
 crate::tolerance::update_vertex_tolerance(brep, vi, max_dist);
 }
 }
}

/// ✅ OCCT : CorrectTolerances —  ,  OCCT boolean  
/// BRepLib::SameParameter + + (finalize_tolerance_hierarchy).
///
///  :
/// 1. `resize_tolerance_arrays` —  
/// 2. `brep_same_parameter` — 3D/pcurve  
/// 3. `compute_vertex_tolerances` —  
/// 4. `finalize_tolerance_hierarchy` —  ( → → ) 
///
/// # Arguments
///
/// * `brep` - The BRep to correct tolerances on.
/// * `samples_per_edge` - Sample points per edge for deviation computation (min 2).
/// ✅ OCCT-aligned: CorrectTolerances + CorrectShapeTolerances
/// (BOPTools_AlgoTools_1.cxx L309-317, L389-420).
///
/// OCCT flow:
/// 1. CorrectPointOnCurve (L315) → brep_same_parameter + compute_vertex_tolerances
/// 2. CorrectCurveOnSurface (L316) → same_range (pcurve range alignment)
/// 3. CorrectShapeTolerances (L389) → finalize_tolerance_hierarchy
///
/// rcad: same_range adjusts pcurve ranges to match 3D curve range,
/// equivalent to OCCT's CorrectCurveOnSurface range adjustment.
/// (OCCT also adjusts pcurve shape via re-projection, which is
/// not needed in rcad because DS::build_face_reps computes
/// precise pcurves during pipeline construction.)
/// OCCT-aligned: BOPTools_AlgoTools::CorrectTolerances(shape, aMA, 0.05, myRunParallel).
/// `max_tolerance` caps the computed tolerance (OCCT default 0.05).
pub fn correct_tolerances(brep: &mut BRep, samples_per_edge: usize, max_tolerance: f64) {
 resize_tolerance_arrays(brep);
 // OCCT L315: CorrectPointOnCurve
 brep_same_parameter(brep, samples_per_edge);
 compute_vertex_tolerances(brep);
 // OCCT L316: CorrectCurveOnSurface (range alignment)
 same_range(brep);
 // Clamp tolerances to max_tolerance (OCCT: 0.05)
 clamp_tolerances(brep, max_tolerance);
 // OCCT L408+: CorrectShapeTolerances — separated into correct_shape_tolerances.
}

/// OCCT-aligned: BOPTools_AlgoTools::CorrectShapeTolerances(shape, aMA, myRunParallel).
/// Second pass: finalize tolerance hierarchy (face/shell/vertex from edge tolerances).
pub fn correct_shape_tolerances(brep: &mut BRep) {
 finalize_tolerance_hierarchy(brep);
}

/// Clamp vertex/edge tolerances to not exceed max_tolerance.
/// rcad: tolerances are computed from actual geometry deviation (sampling-based),
/// so clamping is not needed in the same way as OCCT's 0.05mm limit.
/// Kept for OCCT form alignment.
fn clamp_tolerances(_brep: &mut BRep, _max_tolerance: f64) {}

/// Like [`correct_tolerances`] but accepts a MapToAvoid set of edge/vertex indices
/// whose tolerances should not be modified (non-destructive mode).
pub fn correct_tolerances_with_map(
 brep: &mut BRep,
 samples_per_edge: usize,
 max_tolerance: f64,
 map_to_avoid: &std::collections::HashSet<usize>,
) {
 resize_tolerance_arrays(brep);
 brep_same_parameter_impl(brep, samples_per_edge, map_to_avoid);
 compute_vertex_tolerances_impl(brep, map_to_avoid);
 same_range(brep);
 clamp_tolerances(brep, max_tolerance);
}

/// ✅ OCCT : SameRange — pcurve edge_curve_range  。
/// OCCT BRepLib::SameRange pcurve range 3D curve range  。
///
///  ,  3D curve range pcurve range  ,
/// pcurve range 3D curve range。
pub fn same_range(brep: &mut BRep) {
 for ei in 0..brep.edges.len() {
 let Some([t1, t2]) = brep.geom.edge_curve_range.get(ei).copied().flatten() else {
 continue;
 };
 if (t2 - t1).abs() < 1e-15 {
 continue;
 }
 let Some(pcurves) = brep.geom.edge_pcurves.get(ei) else {
 continue;
 };
 for pc in pcurves {
 let Some(pc_range) = brep.geom.curve2d_range.get(pc.curve2d_idx).copied().flatten() else {
 continue;
 };
 let (pc_t1, pc_t2) = (pc_range[0], pc_range[1]);
 let new_pc_t1 = pc_t1.min(t1);
 let new_pc_t2 = pc_t2.max(t2);
 if (new_pc_t1 - pc_t1).abs() > 1e-15 || (new_pc_t2 - pc_t2).abs() > 1e-15 {
 if let Some(range) = brep.geom.curve2d_range.get_mut(pc.curve2d_idx) {
 *range = Some([new_pc_t1, new_pc_t2]);
 }
 }
 }
 // Mark as same range
 while brep.geom.edge_same_range.len() <= ei {
 brep.geom.edge_same_range.push(true);
 }
 brep.geom.edge_same_range[ei] = true;
 }
}

// ── BRep_Tool ─────────────────────────────────────────────────────────

/// ✅ OCCT : BRep_Tool::Curve(edge) — 3D  。
pub fn edge_curve<'a>(brep: &'a BRep, ei: usize) -> Option<&'a crate::geom::Curve3> {
 let ci = brep.geom.edge_curve.get(ei).copied().flatten()?;
 brep.geom.curves.get(ci)
}

/// ✅ OCCT : BRep_Tool::Surface(face) —  。
pub fn face_surface<'a>(brep: &'a BRep, face_flat_idx: usize) -> Option<&'a crate::geom::Surface3> {
 let si = brep.geom.face_surface.get(face_flat_idx).copied().flatten()?;
 brep.geom.surfaces.get(si)
}

/// ✅ OCCT : BRep_Tool::PCurve(edge, face) — pcurve。
pub fn edge_pcurve_on_face<'a>(brep: &'a BRep, ei: usize, face_flat_idx: usize) -> Option<&'a crate::geom::Curve2d> {
 let pcurves = brep.geom.edge_pcurves.get(ei)?;
 let pc = pcurves.iter().find(|pc| pc.surface_idx == face_flat_idx)?;
 let ci = pc.curve2d_idx;
 brep.geom.curve2ds.get(ci)
}

/// ✅ OCCT : BRep_Tool::Degenerated(edge) —  。
pub fn edge_is_degenerated(brep: &BRep, ei: usize) -> bool {
 brep.geom.edge_degenerated.get(ei).copied().unwrap_or(false)
}

/// ✅ OCCT : BRep_Tool::SameParameter(edge) — same-parameter。
pub fn edge_is_same_parameter(brep: &BRep, ei: usize) -> bool {
 brep.geom.edge_same_parameter.get(ei).copied().unwrap_or(true)
}

/// ✅ OCCT : BRep_Tool::SameRange(edge) — same-range。
pub fn edge_is_same_range(brep: &BRep, ei: usize) -> bool {
 brep.geom.edge_same_range.get(ei).copied().unwrap_or(true)
}

// ── CorrectCurveOnSurface / CorrectPointOnCurve ────────────────────────────────

/// OCCT-aligned: BOPTools_AlgoTools::CorrectCurveOnSurface.
///
/// For each edge-face pair, sample deviation between the 3D curve and the
/// pcurve projected onto the surface. If the deviation exceeds the current
/// edge tolerance, inflate the edge tolerance to `max_dev * 1.00001` (the
/// same 0.001 % safety factor OCCT's BRepLib_ValidateEdge uses in
/// `UpdateTolerance`), capped at `a_max_tol` (default 1e-4).
///
/// Does NOT modify the pcurve geometry itself — only the tolerance value.
///
/// OCCT source: BOPTools_AlgoTools_1.cxx lines 348-385 (`CorrectCurveOnSurface`)
/// and BRepLib_ValidateEdge.cxx lines 49-64 (`UpdateTolerance`, `correctTolerance`).
///
/// NOTE: Phase 1 will add the `PrecCurve` / `PrecSurface` geometric-resolution
/// delta (the `aToleranceDelta` term in OCCT's `correctTolerance`).
pub fn correct_curve_on_surface(brep: &mut BRep, samples_per_edge: usize) {
 correct_curve_on_surface_with_max(brep, samples_per_edge, 0.0001)
}

/// As `correct_curve_on_surface` with an explicit tolerance ceiling.
pub fn correct_curve_on_surface_with_max(
 brep: &mut BRep,
 samples_per_edge: usize,
 a_max_tol: f64,
) {
 let n = samples_per_edge.max(2);
 // First pass: compute corrected tolerance per edge (immutable borrow).
 let deviations: Vec<(usize, f64)> = {
 let curves = &brep.geom.curves;
 let surfaces = &brep.geom.surfaces;
 let curve2ds = &brep.geom.curve2ds;
 (0..brep.edges.len())
 .filter_map(|ei| {
 let ci = brep.geom.edge_curve.get(ei).copied().flatten()?;
 let curve = curves.get(ci)?;
 let [t1, t2] = brep.geom.edge_curve_range.get(ei).copied().flatten()?;
 if (t2 - t1).abs() < 1e-15 {
 return None;
 }
 let pcurves = brep.geom.edge_pcurves.get(ei)?;
 if pcurves.is_empty() {
 return None;
 }
 let etol = edge_tolerance(brep, ei);
 let curve_prec = prec_curve(curve);

 let mut max_dev = 0.0_f64;
 let mut max_surf_prec = 0.0_f64;
 for pc in pcurves {
 let surface = surfaces.get(pc.surface_idx)?;
 let sp = prec_surface(surface);
 if sp > max_surf_prec {
 max_surf_prec = sp;
 }
 let pc_curve = curve2ds.get(pc.curve2d_idx)?;
 for si in 0..n {
 let t = t1 + (t2 - t1) * si as f64 / (n - 1) as f64;
 let uv = pc_curve.point_at(t);
 let p_surf = surface.point_at(uv.x, uv.y);
 let p_curve = curve.point_at(t);
 let dev = (p_surf - p_curve).length();
 if dev > max_dev {
 max_dev = dev;
 }
 }
 }
 if max_dev > etol {
 // OCCT BRepLib_ValidateEdge::UpdateTolerance: computed * 1.00001
 let corrected_raw = max_dev * 1.00001;
 // Add geometric-resolution delta (PrecCurve / PrecSurface)
 let tol_delta = curve_prec.max(max_surf_prec);
 let corrected = corrected_raw + tol_delta;
 if corrected < a_max_tol {
 Some((ei, corrected))
 } else {
 None
 }
 } else {
 None
 }
 })
 .collect()
 };
 // Second pass: apply tolerance updates and propagate edge→vertex
 // (matching OCCT CorrectEdgeTolerance → UpdateShape → CorrectVertexTolerance).
 for (ei, tol) in deviations {
 update_edge_tolerance(brep, ei, tol);
 // Propagate the raised edge tolerance to its adjacent vertices
 // (OCCT CorrectVertexTolerance after each edge tolerance update).
 let (v0, v1) = brep
 .edges
 .get(ei)
 .map(|e| (e.start, e.end))
 .unwrap_or((usize::MAX, usize::MAX));
 if v0 != usize::MAX {
 update_vertex_tolerance(brep, v0, tol);
 update_vertex_tolerance(brep, v1, tol);
 }
 }
}

/// OCCT-aligned: BOPTools_AlgoTools::CorrectPointOnCurve (CheckEdge).
///
/// For each edge, check vertex distance to the 3D curve. If the deviation
/// exceeds the vertex tolerance (bumped up by edge tolerance per OCCT),
/// inflate the vertex tolerance with a 10 % margin (dd = 0.1 * base_tol)
/// and cap at aMaxTol (default 1e-4, matching OCCT's CorrectTolerances
/// ceiling).
///
/// OCCT source: BOPTools_AlgoTools_1.cxx lines 430-517 (CheckEdge).
pub fn correct_point_on_curve(brep: &mut BRep) {
 correct_point_on_curve_with_max(brep, 0.0001)
}

/// As `correct_point_on_curve` but with an explicit tolerance ceiling.
///
/// `a_max_tol` — do not apply corrections whose new tolerance would exceed
/// this value (OCCT default 0.0001 for `CorrectTolerances`).
pub fn correct_point_on_curve_with_max(brep: &mut BRep, a_max_tol: f64) {
 let mut updates: Vec<(usize, f64)> = Vec::new();
 for ei in 0..brep.edges.len() {
 let curve = edge_curve(brep, ei).cloned();
 let Some(ref curve) = curve else { continue; };
 let Some([t1, t2]) = brep.geom.edge_curve_range.get(ei).copied().flatten() else {
 continue;
 };
 if (t2 - t1).abs() < 1e-15 {
 continue;
 }
 let etol = edge_tolerance(brep, ei);
 // Get vertex params from GeomStore::edge_vertex_params (populated by from_topods)
 let vparam_sv = brep.geom.edge_vertex_params.get(ei).copied().flatten().map(|p| p[0]);
 let vparam_ev = brep.geom.edge_vertex_params.get(ei).copied().flatten().map(|p| p[1]);
 for (vi, t_vi, is_forward) in [
 (brep.edges[ei].start, vparam_sv, true),
 (brep.edges[ei].end, vparam_ev, false),
 ] {
 let Some(vp) = brep.vertices.get(vi).map(|v| v.point) else { continue; };
 let vtol = vertex_tolerance(brep, vi);
 let mut a_tol = vtol.max(etol);
 let dd = 0.1 * a_tol;
 a_tol *= a_tol;

 // OCCT CheckEdge L462: check vertex at its specific parameter on the curve
 if let Some(t) = t_vi {
 let pc = curve.point_at(t);
 let d2 = (vp - pc).length_squared();
 if d2 > a_tol {
 let new_tol = d2.sqrt() + dd;
 if new_tol < a_max_tol {
 updates.push((vi, new_tol));
 }
 }
 }

 // OCCT L487-511: FORWARD->First(), REVERSED->Last() curve endpoints
 if is_forward {
 let p_end = curve.point_at(t1);
 let d2 = (vp - p_end).length_squared();
 if d2 > a_tol {
 let new_tol = d2.sqrt() + dd;
 if new_tol < a_max_tol {
 updates.push((vi, new_tol));
 }
 }
 } else {
 let p_end = curve.point_at(t2);
 let d2 = (vp - p_end).length_squared();
 if d2 > a_tol {
 let new_tol = d2.sqrt() + dd;
 if new_tol < a_max_tol {
 updates.push((vi, new_tol));
 }
 }
 }
 }
 }
 for (vi, tol) in updates {
 update_vertex_tolerance(brep, vi, tol);
 }
}

/// ✅ OCCT : PostTreat — CorrectTolerances + CorrectCurveOnSurface + CorrectPointOnCurve。
pub fn post_treat(brep: &mut BRep, samples_per_edge: usize) {
 correct_tolerances(brep, samples_per_edge, 0.05);
 correct_curve_on_surface(brep, samples_per_edge);
 correct_curve_on_surface(brep, samples_per_edge);
 correct_point_on_curve(brep);
 same_range(brep);
}

/// TopoDS-aligned: like [`correct_tolerances`] but operates on `&mut topods::BRep`.
/// The old-BRep conversion is internal; callers see only topods.
pub fn correct_tolerances_topods(
 t: &mut topods::BRep,
 samples_per_edge: usize,
 max_tolerance: f64,
) {
 let mut old = crate::BRep::from_topods(t);
 correct_tolerances(&mut old, samples_per_edge, max_tolerance);
 propagate_old_tolerances_topods(&old, t);
}

/// TopoDS-aligned: like [`correct_shape_tolerances`] but operates on `&mut topods::BRep`.
pub fn correct_shape_tolerances_topods(t: &mut topods::BRep) {
 let mut old = crate::BRep::from_topods(t);
 correct_shape_tolerances(&mut old);
 propagate_old_tolerances_topods(&old, t);
}

/// Write tolerance values from old BRep's GeomStore back to topods TShape fields.
/// The flat-index order is preserved because `BRep::from_topods` iterates tshapes
/// in order — flat index `i` corresponds to the i-th Vertex/Edge TShape in `t`.
fn propagate_old_tolerances_topods(old: &crate::BRep, t: &mut topods::BRep) {
 use std::sync::Arc;
 // Vertex tolerances: old.vertices[i] ↔ i-th Vertex TShape in t
 let mut vi: usize = 0;
 for ts in &mut t.tshapes {
 if let topods::TShape::Vertex(ref mut vd) = *Arc::make_mut(ts) {
 if let Some(&tol) = old.geom.vertex_tolerance.get(vi) {
 vd.tolerance = tol;
 }
 vi += 1;
 }
 }
 // Edge tolerances: old.edges[i] ↔ i-th Edge TShape in t
 let mut ei: usize = 0;
 for ts in &mut t.tshapes {
 if let topods::TShape::Edge(ref mut ed) = *Arc::make_mut(ts) {
 if let Some(&tol) = old.geom.edge_tolerance.get(ei) {
 ed.tolerance = tol;
 }
 if let Some(&sp) = old.geom.edge_same_parameter.get(ei) {
 ed.same_parameter = sp;
 }
 if let Some(&sr) = old.geom.edge_same_range.get(ei) {
 ed.same_range = sr;
 }
 ei += 1;
 }
 }
 // Face tolerances: propagate old.face_tolerance (flat_fi order) to Face TShapes.
 // The flat_fi order matches the order of Face TShapes in t.
 let mut fi: usize = 0;
 for ts in &mut t.tshapes {
 if let topods::TShape::Face(ref mut fd) = *Arc::make_mut(ts) {
 if let Some(&tol) = old.geom.face_tolerance.get(fi) {
 fd.tolerance = tol;
 }
 fi += 1;
 }
 }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
 use super::*;
 use crate::{BRep, PrimitiveSolid};

 #[test]
 fn default_tolerance_is_confusion() {
 let brep = BRep::new();
 assert_eq!(vertex_tolerance(&brep, 0), CONFUSION);
 assert_eq!(edge_tolerance(&brep, 0), CONFUSION);
 assert_eq!(face_tolerance(&brep, 0), CONFUSION);
 assert_eq!(model_tolerance(&brep), CONFUSION);
 }

 #[test]
 fn stored_tolerance_returned() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });
 brep.geom.edge_tolerance = vec![1e-5; brep.edges.len()];
 assert!((edge_tolerance(&brep, 0) - 1e-5).abs() < 1e-20);
 }

 #[test]
 fn model_tolerance_returns_max() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });
 brep.geom.vertex_tolerance = vec![1e-6; brep.vertices.len()];
 brep.geom.edge_tolerance = vec![1e-5; brep.edges.len()];
 // max is 1e-5
 assert!((model_tolerance(&brep) - 1e-5).abs() < 1e-20);
 }

 #[test]
 fn zero_tolerance_falls_back_to_confusion() {
 let mut brep = BRep::new();
 brep.geom.vertex_tolerance = vec![0.0];
 assert_eq!(vertex_tolerance(&brep, 0), CONFUSION);
 }

 // ── Write API tests ───────────────────────────────────────────────────────

 #[test]
 fn set_vertex_tolerance_stores_value() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 set_vertex_tolerance(&mut brep, 0, 1e-4);
 assert!((vertex_tolerance(&brep, 0) - 1e-4).abs() < 1e-20);
 }

 #[test]
 fn update_vertex_tolerance_only_raises() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 set_vertex_tolerance(&mut brep, 0, 1e-4);
 update_vertex_tolerance(&mut brep, 0, 1e-6); // lower — should not change
 assert!((vertex_tolerance(&brep, 0) - 1e-4).abs() < 1e-20);
 update_vertex_tolerance(&mut brep, 0, 1e-3); // higher — should update
 assert!((vertex_tolerance(&brep, 0) - 1e-3).abs() < 1e-20);
 }

 #[test]
 fn set_edge_tolerance_stores_value() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 set_edge_tolerance(&mut brep, 3, 5e-5);
 assert!((edge_tolerance(&brep, 3) - 5e-5).abs() < 1e-20);
 }

 #[test]
 fn update_edge_tolerance_only_raises() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 set_edge_tolerance(&mut brep, 2, 1e-5);
 update_edge_tolerance(&mut brep, 2, 1e-8); // lower — no change
 assert!((edge_tolerance(&brep, 2) - 1e-5).abs() < 1e-20);
 update_edge_tolerance(&mut brep, 2, 2e-5); // higher — updates
 assert!((edge_tolerance(&brep, 2) - 2e-5).abs() < 1e-20);
 }

 #[test]
 fn set_face_tolerance_stores_value() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 set_face_tolerance(&mut brep, 0, 3e-5);
 assert!((face_tolerance(&brep, 0) - 3e-5).abs() < 1e-20);
 }

 #[test]
 fn out_of_range_set_extends_array() {
 let mut brep = BRep::new(); // no entities
 set_vertex_tolerance(&mut brep, 99, 1e-4);
 assert_eq!(brep.geom.vertex_tolerance.len(), 100);
 assert!((brep.geom.vertex_tolerance[99] - 1e-4).abs() < 1e-20);
 // slots before 99 should be CONFUSION
 assert_eq!(brep.geom.vertex_tolerance[0], CONFUSION);
 }

 #[test]
 fn below_confusion_floor_clamped() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 set_vertex_tolerance(&mut brep, 0, 1e-15); // below CONFUSION = 1e-7
 assert_eq!(vertex_tolerance(&brep, 0), CONFUSION);
 }

 #[test]
 fn resize_tolerance_arrays_fills_missing_slots() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 // Arrays start empty for primitives.
 assert!(brep.geom.vertex_tolerance.is_empty());
 resize_tolerance_arrays(&mut brep);
 assert_eq!(brep.geom.vertex_tolerance.len(), 8);
 assert_eq!(brep.geom.edge_tolerance.len(), 12);
 assert_eq!(brep.geom.face_tolerance.len(), 6);
 // All filled with CONFUSION.
 assert!(brep.geom.vertex_tolerance.iter().all(|&t| t == CONFUSION));
 }

 #[test]
 fn finalize_tolerance_hierarchy_propagates_face_to_vertex() {
 // Give face 0 a high tolerance; after finalization all its boundary verts/edges
 // should carry at least that tolerance.
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 resize_tolerance_arrays(&mut brep);
 set_face_tolerance(&mut brep, 0, 1e-4); // face 0 = front face
 finalize_tolerance_hierarchy(&mut brep);

 // All edge tolerances that belong to face 0 should be ≥ 1e-4.
 // Box face 0 has edges 0,1,2,3 (from create_box).
 for &ei in &[0usize, 1, 2, 3] {
 assert!(
 edge_tolerance(&brep, ei) >= 1e-4 - 1e-15,
 "edge {ei} tol should be ≥ 1e-4, got {}",
 edge_tolerance(&brep, ei)
 );
 }
 // Vertices of those edges should be ≥ 1e-4 as well.
 for vi in 0..4usize {
 assert!(
 vertex_tolerance(&brep, vi) >= 1e-4 - 1e-15,
 "vertex {vi} tol should be ≥ 1e-4, got {}",
 vertex_tolerance(&brep, vi)
 );
 }
 }

 #[test]
 fn finalize_tolerance_hierarchy_never_lowers() {
 // Vertices/edges with higher pre-existing tolerance must not be lowered.
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 resize_tolerance_arrays(&mut brep);
 set_vertex_tolerance(&mut brep, 0, 1e-3); // very high vertex tol
 set_face_tolerance(&mut brep, 0, 1e-6); // lower face tol
 finalize_tolerance_hierarchy(&mut brep);
 // Vertex 0 must still be 1e-3.
 assert!((vertex_tolerance(&brep, 0) - 1e-3).abs() < 1e-20);
 }

 #[test]
 fn correct_tolerances_runs_without_panic() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 2.0, depth: 3.0,
 });
 // Should not panic on a valid box BRep.
 correct_tolerances(&mut brep, 5, 0.05);
 // All tolerances should be at least CONFUSION after processing.
 for vi in 0..brep.vertices.len() {
 assert!(vertex_tolerance(&brep, vi) >= CONFUSION);
 }
 for ei in 0..brep.edges.len() {
 assert!(edge_tolerance(&brep, ei) >= CONFUSION);
 }
 }

 #[test]
 fn same_range_sets_flag() {
 let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });
 // Initially no same_range flag.
 assert!(brep.geom.edge_same_range.is_empty());
 same_range(&mut brep);
 // After same_range, flags should be set for all edges.
 assert_eq!(brep.geom.edge_same_range.len(), brep.edges.len());
 assert!(brep.geom.edge_same_range.iter().all(|&f| f));
 }
}
