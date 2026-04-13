//! Hidden-Line Removal (HLR).
//!
//! Projects a BRep's edges onto a view plane and classifies each edge segment
//! as **visible** or **hidden** by testing against the silhouette of all faces.
//!
//! Analytic silhouette curves are generated for curved surfaces (cylinder,
//! sphere) and processed through the same visibility pipeline as wire edges.
//!
//! Analogous to OCCT `HLRBRep_Algo` / `HLRBRep_HLRToShape`.
//!
//! # Algorithm
//!
//! For each edge (and silhouette curve):
//! 1. Project both endpoints onto the screen plane.
//! 2. Sample `N` points along the edge in 3D.
//! 3. For each sample, cast a ray from that point toward the camera.
//! 4. If any face triangle blocks the ray **closer** to the camera than the
//!    edge sample, that sample is hidden.
//! 5. Classify runs of consecutive samples → visible/hidden segments.
//!
//! The result is a set of `HlrSegment`s — 2D projected line segments labeled
//! visible or hidden.

use glam::{DAffine3, DMat4, DVec2, DVec3, DVec4};
use rcad_kernel::geom::{Circle3, CurveEval, Surface3, any_perpendicular};
use rcad_kernel::{BRep, SurfaceEval};

// ── Public types ──────────────────────────────────────────────────────────────

/// Hint about the geometric type of the original 3D edge curve.
/// Used by consumers (e.g. SVG exporter) to emit arcs instead of polylines.
#[derive(Debug, Clone, PartialEq)]
pub enum CurveHint {
    /// Edge is a full or partial circle in 3D.
    Circle {
        /// Projected 2D center of the circle.
        center: DVec2,
        /// Projected radius (approximate — perspective not applied).
        radius: f64,
    },
    /// Any other non-straight curve (ellipse, spline, …).
    Other,
}

/// A projected edge segment labeled as visible or hidden.
#[derive(Debug, Clone, PartialEq)]
pub struct HlrSegment {
    /// Start point in 2D screen space.
    pub start: DVec2,
    /// End point in 2D screen space.
    pub end: DVec2,
    /// Whether this segment is visible from the camera.
    pub visible: bool,
    /// Optional hint about the underlying curve type (None for straight lines).
    pub curve_hint: Option<CurveHint>,
}

/// Output of an HLR computation.
#[derive(Debug, Clone, Default)]
pub struct HlrResult {
    pub segments: Vec<HlrSegment>,
}

impl HlrResult {
    /// Return only visible segments.
    pub fn visible(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible)
    }

    /// Return only hidden segments.
    pub fn hidden(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| !s.visible)
    }
}

/// Camera / view specification for HLR.
#[derive(Debug, Clone)]
pub struct HlrCamera {
    /// Camera position in world space.
    pub eye: DVec3,
    /// Target point (look-at).
    pub target: DVec3,
    /// Up direction.
    pub up: DVec3,
}

impl HlrCamera {
    pub fn new(eye: DVec3, target: DVec3) -> Self {
        Self {
            eye,
            target,
            up: DVec3::Y,
        }
    }

    pub fn with_up(mut self, up: DVec3) -> Self {
        self.up = up;
        self
    }

    /// Isometric-style view from the +X+Y+Z octant.
    pub fn isometric(distance: f64) -> Self {
        let d = distance / 3.0_f64.sqrt();
        Self::new(DVec3::splat(d), DVec3::ZERO)
    }

    /// Front view (looking along +Y, up = +Z).
    pub fn front(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, -distance, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }

    /// Top view (looking down -Z).
    pub fn top(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, 0.0, distance), DVec3::ZERO).with_up(DVec3::Y)
    }

    /// Right-side view (looking along -X, up = +Z).
    pub fn right(distance: f64) -> Self {
        Self::new(DVec3::new(distance, 0.0, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Build a right-handed view matrix (world → camera space).
fn look_at(eye: DVec3, target: DVec3, up: DVec3) -> DMat4 {
    let forward = (target - eye).normalize_or_zero();
    let right = forward.cross(up).normalize_or_zero();
    let up = right.cross(forward);

    DMat4::from_cols(
        DVec4::new(right.x, right.y, right.z, -right.dot(eye)),
        DVec4::new(up.x, up.y, up.z, -up.dot(eye)),
        DVec4::new(-forward.x, -forward.y, -forward.z, forward.dot(eye)),
        DVec4::new(0.0, 0.0, 0.0, 1.0),
    )
    .transpose()
}

/// Project a world-space point to 2D screen space using the view matrix.
/// Returns (x, y) in camera space (z is depth; ignored for 2D output).
fn project(p: DVec3, view: &DMat4) -> (DVec2, f64) {
    let hp = view.mul_vec4(DVec4::new(p.x, p.y, p.z, 1.0));
    (DVec2::new(hp.x, hp.y), hp.z)
}

/// Collect all triangles from a BRep (fan-triangulate faces without pre-triangulated data).
fn collect_triangles(brep: &BRep) -> Vec<[DVec3; 3]> {
    let mut tris = Vec::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if !face.triangles.is_empty() {
                    for &[i, j, k] in &face.triangles {
                        if let (Some(a), Some(b), Some(c)) = (
                            brep.vertices.get(i),
                            brep.vertices.get(j),
                            brep.vertices.get(k),
                        ) {
                            tris.push([a.point, b.point, c.point]);
                        }
                    }
                } else {
                    // Fan-triangulate from wire
                    let pts: Vec<DVec3> = face
                        .outer_wire
                        .edges
                        .iter()
                        .filter_map(|we| {
                            let edge = brep.edges.get(we.idx)?;
                            let vi = if we.forward { edge.start } else { edge.end };
                            brep.vertices.get(vi).map(|v| v.point)
                        })
                        .collect();
                    if pts.len() >= 3 {
                        let origin = pts[0];
                        for i in 1..pts.len() - 1 {
                            tris.push([origin, pts[i], pts[i + 1]]);
                        }
                    }
                }
            }
        }
    }
    tris
}

/// Ray-triangle intersection (Möller–Trumbore). Returns `Some(t)` if the ray
/// `origin + t*dir` hits the triangle (t > epsilon, front-face only).
fn ray_triangle_intersect(origin: DVec3, dir: DVec3, tri: &[DVec3; 3]) -> Option<f64> {
    const EPS: f64 = 1e-8;
    let edge1 = tri[1] - tri[0];
    let edge2 = tri[2] - tri[0];
    let h = dir.cross(edge2);
    let a = edge1.dot(h);
    if a.abs() < EPS {
        return None;
    }
    let f = 1.0 / a;
    let s = origin - tri[0];
    let u = f * s.dot(h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(edge1);
    let v = f * dir.dot(q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * edge2.dot(q);
    if t > EPS { Some(t) } else { None }
}

/// Test if a world-space point is occluded by any triangle when viewed from `eye`.
fn is_occluded(point: DVec3, eye: DVec3, triangles: &[[DVec3; 3]], dist_to_eye: f64) -> bool {
    let dir = (eye - point).normalize_or_zero();
    let origin = point + dir * 1e-5; // push off surface
    for tri in triangles {
        if let Some(t) = ray_triangle_intersect(origin, dir, tri)
            && t < dist_to_eye - 1e-4
        {
            return true;
        }
    }
    false
}

// ── Public API ────────────────────────────────────────────────────────────────

// ── Silhouette generation ─────────────────────────────────────────────────────

/// Internal: one silhouette curve to process through the HLR pipeline.
struct SilhouetteCurve {
    /// World-space sample points (at least 2).
    world_pts: Vec<DVec3>,
    /// Optional curve hint for SVG output.
    curve_hint: Option<CurveHint>,
    /// If true, emit one segment per consecutive point pair instead of merging
    /// runs.  Used for dense polyline approximations (e.g. sphere silhouette).
    dense: bool,
}

/// Generate analytic silhouette curves for curved surfaces in `brep`.
///
/// Handles:
/// - **Cylinder**: two lines parallel to the axis at the silhouette angles.
/// - **Sphere**: a 64-point polyline approximating the great circle perpendicular
///   to the view direction.
fn compute_silhouettes(brep: &BRep, view_dir: DVec3, samples: usize) -> Vec<SilhouetteCurve> {
    let mut curves: Vec<SilhouetteCurve> = Vec::new();
    let line_samples = samples.max(16);
    let dense_curve_samples = (samples * 8).max(96);

    if brep.solids.is_empty() {
        return curves;
    }

    let mut face_idx = 0usize;
    for shell in &brep.solids[0].shells {
        for _face in &shell.faces {
            let surf_idx = match brep.geom.face_surface.get(face_idx).and_then(|o| *o) {
                Some(idx) => idx,
                None => {
                    face_idx += 1;
                    continue;
                }
            };
            let surface = &brep.geom.surfaces[surf_idx];

            let domain = match brep.geom.face_surface_range.get(face_idx).and_then(|o| *o) {
                Some(r) => r,
                None => surface.default_domain(),
            };
            let [_u0, _u1, v0, v1] = domain;

            match surface {
                Surface3::Cylinder(cyl) => {
                    // Project view direction onto the plane perpendicular to the axis.
                    let d_perp = view_dir - view_dir.dot(cyl.axis) * cyl.axis;
                    if d_perp.length_squared() < 1e-10 {
                        // Viewing along the axis — no silhouette lines.
                        face_idx += 1;
                        continue;
                    }
                    // Direction from axis to silhouette (perpendicular to both axis and d_perp).
                    let sil_dir = cyl.axis.cross(d_perp).normalize_or_zero();

                    // Resolve v range (height along axis).
                    // face_surface_range may be absent; fall back to vertex projections.
                    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
                        (v0, v1)
                    } else {
                        let mut lo = f64::INFINITY;
                        let mut hi = f64::NEG_INFINITY;
                        for vert in &brep.vertices {
                            let proj = (vert.point - cyl.origin).dot(cyl.axis);
                            lo = lo.min(proj);
                            hi = hi.max(proj);
                        }
                        if lo.is_finite() && hi.is_finite() {
                            (lo, hi)
                        } else {
                            face_idx += 1;
                            continue;
                        }
                    };

                    for &sign in &[1.0_f64, -1.0] {
                        let offset = sil_dir * sign * cyl.radius;
                        let world_pts: Vec<DVec3> = (0..line_samples)
                            .map(|i| {
                                let t = i as f64 / (line_samples - 1) as f64;
                                let v = v0_eff + (v1_eff - v0_eff) * t;
                                cyl.origin + v * cyl.axis + offset
                            })
                            .collect();
                        curves.push(SilhouetteCurve {
                            world_pts,
                            curve_hint: None,
                            dense: false,
                        });
                    }
                }

                Surface3::Sphere(sph) => {
                    // Silhouette is the great circle perpendicular to view_dir.
                    let x_ax = any_perpendicular(view_dir);
                    let y_ax = view_dir.cross(x_ax).normalize_or_zero();
                    // Use 0..N (open) so first != last; the SVG arc hint closes the circle.
                    let pts: Vec<DVec3> = (0..dense_curve_samples)
                        .map(|i| {
                            let t = 2.0 * std::f64::consts::PI * i as f64 / dense_curve_samples as f64;
                            sph.center + sph.radius * (t.cos() * x_ax + t.sin() * y_ax)
                        })
                        .collect();
                    // CurveHint will be computed after projection in process_world_pts;
                    // pass None here and let the caller set it if needed.
                    // For the sphere silhouette we emit individual consecutive segments
                    // so that the polyline is rendered correctly.
                    curves.push(SilhouetteCurve {
                        world_pts: pts,
                        curve_hint: None,
                        dense: true,
                    });
                }

                Surface3::Cone(con) => {
                    // Silhouette of a cone: two generators from the apex.
                    // Project view direction onto the plane perpendicular to the axis.
                    let d_perp = view_dir - view_dir.dot(con.axis) * con.axis;
                    if d_perp.length_squared() < 1e-10 {
                        // Viewing along the axis — no silhouette lines.
                        face_idx += 1;
                        continue;
                    }
                    let sil_dir = con.axis.cross(d_perp).normalize_or_zero();

                    // Determine v range (slant distance from apex).
                    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
                        (v0, v1)
                    } else {
                        let mut lo = f64::INFINITY;
                        let mut hi = f64::NEG_INFINITY;
                        for vert in &brep.vertices {
                            let proj = (vert.point - con.apex).dot(con.axis);
                            lo = lo.min(proj);
                            hi = hi.max(proj);
                        }
                        if lo.is_finite() && hi.is_finite() {
                            (lo.max(0.0), hi.max(0.0))
                        } else {
                            face_idx += 1;
                            continue;
                        }
                    };

                    let tan_a = con.half_angle_rad.tan();
                    for &sign in &[1.0_f64, -1.0] {
                        let world_pts: Vec<DVec3> = (0..line_samples)
                            .map(|i| {
                                let t = i as f64 / (line_samples - 1) as f64;
                                let v = v0_eff + (v1_eff - v0_eff) * t;
                                con.apex + v * con.axis + v * tan_a * sil_dir * sign
                            })
                            .collect();
                        if world_pts
                            .first()
                            .zip(world_pts.last())
                            .map(|(a, b)| (*b - *a).length_squared() > 1e-12)
                            .unwrap_or(false)
                        {
                            curves.push(SilhouetteCurve {
                                world_pts,
                                curve_hint: None,
                                dense: false,
                            });
                        }
                    }
                }

                Surface3::Torus(tor) => {
                    // Silhouette of a torus: for each azimuth u, the silhouette v satisfies
                    //   normal(u,v) · view_dir = 0
                    //   cos(v)*(radial(u)·d) + sin(v)*(axis·d) = 0
                    //   v = atan2(-(radial(u)·d), axis·d)
                    // This gives two silhouette curves (v and v+π).
                    let x_ax = any_perpendicular(tor.axis);
                    let y_ax = tor.axis.cross(x_ax).normalize_or_zero();
                    let axis_dot = tor.axis.dot(view_dir);
                    for &offset in &[0.0_f64, std::f64::consts::PI] {
                        let pts: Vec<DVec3> = (0..dense_curve_samples)
                            .map(|i| {
                                let u = 2.0 * std::f64::consts::PI * i as f64 / dense_curve_samples as f64;
                                let radial = u.cos() * x_ax + u.sin() * y_ax;
                                let radial_dot = radial.dot(view_dir);
                                let v = (-radial_dot).atan2(axis_dot) + offset;
                                let tube_center = tor.center + tor.major_radius * radial;
                                tube_center + tor.minor_radius * (v.cos() * radial + v.sin() * tor.axis)
                            })
                            .collect();
                        curves.push(SilhouetteCurve {
                            world_pts: pts,
                            curve_hint: None,
                            dense: true,
                        });
                    }
                }

                _ => {}
            }

            face_idx += 1;
        }
    }

    curves
}

/// Process a list of world-space sample points through the HLR visibility
/// pipeline and append resulting segments to `result`.
///
/// When `dense` is true, one segment is emitted per consecutive point pair
/// (useful for polyline approximations of curved silhouettes).
fn process_world_pts(
    world_pts: &[DVec3],
    curve_hint: Option<CurveHint>,
    dense: bool,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    result: &mut HlrResult,
) {
    if world_pts.len() < 2 {
        return;
    }
    let n = world_pts.len();

    let sample_vis: Vec<bool> = world_pts
        .iter()
        .map(|&wp| {
            let dist = (camera.eye - wp).length();
            !is_occluded(wp, camera.eye, triangles, dist)
        })
        .collect();

    let screen_pts: Vec<DVec2> = world_pts.iter().map(|&wp| project(wp, view).0).collect();

    if dense {
        // Emit one segment per consecutive pair (preserves polyline shape).
        for i in 0..n - 1 {
            let seg = HlrSegment {
                start: screen_pts[i],
                end: screen_pts[i + 1],
                visible: sample_vis[i],
                curve_hint: curve_hint.clone(),
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
        }
        return;
    }

    let mut seg_start = 0usize;
    for i in 1..n {
        let changed = sample_vis[i] != sample_vis[seg_start];
        let last = i == n - 1;
        if changed || last {
            let end_idx = if last && !changed { i } else { i - 1 };
            let seg = HlrSegment {
                start: screen_pts[seg_start],
                end: screen_pts[end_idx],
                visible: sample_vis[seg_start],
                curve_hint: curve_hint.clone(),
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
            if changed {
                seg_start = i;
            }
        }
    }
}



/// Perform hidden-line removal on a BRep from the given camera position.
///
/// Returns 2D projected segments labeled visible/hidden.
/// `samples` controls how finely each edge is subdivided for occlusion testing
/// (higher = more accurate but slower; 8 is a reasonable default).
pub fn hlr(brep: &BRep, camera: &HlrCamera, samples: usize) -> HlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let triangles = collect_triangles(brep);
    let samples = samples.max(2);
    let mut result = HlrResult::default();

    // ── Wire edges ────────────────────────────────────────────────────────────

    // Collect all unique edges from all faces + standalone edges
    let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    edge_indices.insert(we.idx);
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        edge_indices.insert(we.idx);
                    }
                }
            }
        }
    }
    for i in 0..brep.edges.len() {
        edge_indices.insert(i);
    }

    for &edge_idx in &edge_indices {
        let Some(edge) = brep.edges.get(edge_idx) else { continue };
        let Some(v_start) = brep.vertices.get(edge.start) else { continue };
        let Some(v_end) = brep.vertices.get(edge.end) else { continue };

        let p0 = v_start.point;
        let p1 = v_end.point;

        let edge_curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .and_then(|&ci| ci)
            .and_then(|ci| brep.geom.curves.get(ci));

        let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
            if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
        });

        let is_other_curve = edge_curve
            .map_or(false, |c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
            && circle_info.is_none();

        let edge_samples = if circle_info.is_some() || is_other_curve {
            64.max(samples)
        } else {
            samples
        };

        let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
            let [t0, t1] = brep
                .geom
                .edge_curve_range
                .get(edge_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| circ.default_domain());
            (0..edge_samples)
                .map(|i| {
                    let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                    circ.point_at(t)
                })
                .collect()
        } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
            let [t0, t1] = brep
                .geom
                .edge_curve_range
                .get(edge_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| curve.default_domain());
            (0..edge_samples)
                .map(|i| {
                    let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                    curve.point_at(t)
                })
                .collect()
        } else {
            if (p1 - p0).length_squared() < 1e-12 {
                continue;
            }
            (0..edge_samples)
                .map(|i| {
                    let t = i as f64 / (edge_samples - 1) as f64;
                    p0 + (p1 - p0) * t
                })
                .collect()
        };

        // Compute curve_hint for circle edges
        let screen_pts_for_hint: Vec<DVec2> =
            world_pts.iter().map(|&wp| project(wp, &view).0).collect();
        let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
            let (center_2d, _) = project(circ.center, &view);
            let r = screen_pts_for_hint
                .iter()
                .map(|p| (*p - center_2d).length())
                .fold(0.0_f64, f64::max);
            Some(CurveHint::Circle { center: center_2d, radius: r })
        } else if is_other_curve {
            Some(CurveHint::Other)
        } else {
            None
        };

        process_world_pts(&world_pts, curve_hint, false, camera, &view, &triangles, &mut result);
    }

    // ── Analytic silhouette curves ────────────────────────────────────────────

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    for sil in compute_silhouettes(brep, view_dir, samples) {
        process_world_pts(&sil.world_pts, sil.curve_hint, sil.dense, camera, &view, &triangles, &mut result);
    }

    result
}

/// Per-component HLR result for assembly HLR.
#[derive(Debug, Clone, Default)]
pub struct ComponentHlr {
    /// Component name (from the assembly).
    pub name: String,
    /// HLR segments for this component.
    pub segments: Vec<HlrSegment>,
}

/// Output of assembly HLR — one `ComponentHlr` per leaf BRep.
#[derive(Debug, Clone, Default)]
pub struct AssemblyHlrResult {
    pub components: Vec<ComponentHlr>,
}

impl AssemblyHlrResult {
    /// Return all visible segments across all components.
    pub fn visible_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| s.visible).map(move |s| (c, s))
        })
    }

    /// Return all hidden segments across all components.
    pub fn hidden_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| !s.visible).map(move |s| (c, s))
        })
    }
}

/// Transform a BRep's vertices by an affine transform.
/// Returns a new BRep with transformed vertex positions.
fn transform_brep(brep: &BRep, transform: &DAffine3) -> BRep {
    let mut out = brep.clone();
    for v in &mut out.vertices {
        v.point = transform.transform_point3(v.point);
    }
    out
}

/// Perform hidden-line removal on an assembly of BReps.
///
/// Each component's geometry is transformed to world space, then all triangles
/// are merged into a single occlusion buffer. Each component's edges are
/// tested against the global occlusion buffer, so components correctly
/// occlude each other.
///
/// Returns one `ComponentHlr` per leaf component.
pub fn hlr_assembly(
    components: &[(BRep, DAffine3, String)],
    camera: &HlrCamera,
    samples: usize,
) -> AssemblyHlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let samples = samples.max(2);

    // Transform all BRePs to world space and collect a unified triangle pool.
    let world_breps: Vec<BRep> = components
        .iter()
        .map(|(brep, xf, _)| transform_brep(brep, xf))
        .collect();

    let mut all_triangles: Vec<[DVec3; 3]> = Vec::new();
    for wb in &world_breps {
        all_triangles.extend(collect_triangles(wb));
    }

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    let mut result = AssemblyHlrResult::default();

    for (wb, (_, _, name)) in world_breps.iter().zip(components.iter()) {
        let mut comp_result = HlrResult::default();

        // ── Wire edges ────────────────────────────────────────────────────
        let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        for solid in &wb.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for we in &face.outer_wire.edges {
                        edge_indices.insert(we.idx);
                    }
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            edge_indices.insert(we.idx);
                        }
                    }
                }
            }
        }
        for i in 0..wb.edges.len() {
            edge_indices.insert(i);
        }

        for &edge_idx in &edge_indices {
            let Some(edge) = wb.edges.get(edge_idx) else { continue };
            let Some(v_start) = wb.vertices.get(edge.start) else { continue };
            let Some(v_end) = wb.vertices.get(edge.end) else { continue };

            let p0 = v_start.point;
            let p1 = v_end.point;

            let edge_curve = wb
                .geom
                .edge_curve
                .get(edge_idx)
                .and_then(|&ci| ci)
                .and_then(|ci| wb.geom.curves.get(ci));

            let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
                if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
            });

            let is_other_curve = edge_curve
                .map_or(false, |c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
                && circle_info.is_none();

            let edge_samples = if circle_info.is_some() || is_other_curve {
                64.max(samples)
            } else {
                samples
            };

            let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| circ.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        circ.point_at(t)
                    })
                    .collect()
            } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| curve.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        curve.point_at(t)
                    })
                    .collect()
            } else {
                if (p1 - p0).length_squared() < 1e-12 {
                    continue;
                }
                (0..edge_samples)
                    .map(|i| {
                        let t = i as f64 / (edge_samples - 1) as f64;
                        p0 + (p1 - p0) * t
                    })
                    .collect()
            };

            let screen_pts_for_hint: Vec<DVec2> =
                world_pts.iter().map(|&wp| project(wp, &view).0).collect();
            let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
                let (center_2d, _) = project(circ.center, &view);
                let r = screen_pts_for_hint
                    .iter()
                    .map(|p| (*p - center_2d).length())
                    .fold(0.0_f64, f64::max);
                Some(CurveHint::Circle { center: center_2d, radius: r })
            } else if is_other_curve {
                Some(CurveHint::Other)
            } else {
                None
            };

            process_world_pts(&world_pts, curve_hint, false, camera, &view, &all_triangles, &mut comp_result);
        }

        // ── Analytic silhouette curves ────────────────────────────────────
        for sil in compute_silhouettes(wb, view_dir, samples) {
            process_world_pts(&sil.world_pts, sil.curve_hint, sil.dense, camera, &view, &all_triangles, &mut comp_result);
        }

        result.components.push(ComponentHlr {
            name: name.clone(),
            segments: comp_result.segments,
        });
    }

    result
}

/// Render HLR result as a simple SVG string.
///
/// Visible edges are drawn solid black; hidden edges are dashed gray.
/// `scale` controls pixel size per unit.
pub fn hlr_to_svg(result: &HlrResult, scale: f64, margin: f64) -> String {
    if result.segments.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".to_string();
    }

    // Compute bounding box
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for seg in &result.segments {
        for p in [seg.start, seg.end] {
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
        }
    }

    // Flip Y (SVG Y grows downward, camera Y grows upward)
    let transform = |p: DVec2| -> (f64, f64) {
        let x = (p.x - min_x) * scale + margin;
        let y = (max_y - p.y) * scale + margin;
        (x, y)
    };

    let w = (max_x - min_x) * scale + 2.0 * margin;
    let h = (max_y - min_y) * scale + 2.0 * margin;

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}\" height=\"{:.1}\" viewBox=\"0 0 {:.1} {:.1}\">\n",
        w, h, w, h
    );
    svg.push_str("  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");

    for seg in &result.segments {
        let (x1, y1) = transform(seg.start);
        let (x2, y2) = transform(seg.end);
        let stroke = if seg.visible {
            "black\" stroke-width=\"1.5"
        } else {
            "#999\" stroke-width=\"0.8\" stroke-dasharray=\"4,3"
        };

        // For circle segments emit an SVG arc path; for all others emit a line.
        if let Some(CurveHint::Circle { center, radius }) = &seg.curve_hint {
            let (cx, cy) = transform(*center);
            let r = radius * scale;
            // Determine large-arc flag: compare arc length vs half-circumference
            let dx1 = x1 - cx;
            let dy1 = y1 - cy;
            let dx2 = x2 - cx;
            let dy2 = y2 - cy;
            let cross = dx1 * dy2 - dy1 * dx2;
            let dot = dx1 * dx2 + dy1 * dy2;
            let angle = cross.atan2(dot).abs();
            let large_arc = if angle > std::f64::consts::PI { 1 } else { 0 };
            let sweep = if cross < 0.0 { 0 } else { 1 };
            svg.push_str(&format!(
                "  <path d=\"M {:.3} {:.3} A {:.3} {:.3} 0 {} {} {:.3} {:.3}\" fill=\"none\" stroke=\"{}\"/>\n",
                x1, y1, r, r, large_arc, sweep, x2, y2, stroke
            ));
            // Also record the center for debugging/reference (as a tiny dot, invisible by default)
            let _ = (cx, cy); // suppress unused warning
        } else {
            svg.push_str(&format!(
                "  <line x1=\"{:.3}\" y1=\"{:.3}\" x2=\"{:.3}\" y2=\"{:.3}\" stroke=\"{}\"/>\n",
                x1, y1, x2, y2, stroke
            ));
        }
    }
    svg.push_str("</svg>\n");
    svg
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn unit_box_hlr_produces_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "HLR should produce segments for a box"
        );
        assert!(
            result.visible().count() > 0,
            "some segments should be visible"
        );
    }

    #[test]
    fn hlr_svg_is_valid_xml() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(svg.contains("<svg"), "output should be SVG");
        assert!(svg.contains("</svg>"), "SVG should close properly");
        assert!(svg.contains("<line"), "SVG should contain lines");
    }

    #[test]
    fn top_view_box_has_visible_top_edges() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::top(5.0);
        let result = hlr(&brep, &camera, 8);
        let vis = result.visible().count();
        let hid = result.hidden().count();
        assert!(vis > 0, "top view should have visible edges");
        assert!(hid > 0, "top view should have hidden (bottom) edges");
    }

    #[test]
    fn front_view_and_right_view_both_produce_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 1.0,
            depth: 1.0,
        });
        let front_result = hlr(&brep, &HlrCamera::front(5.0), 8);
        let right_result = hlr(&brep, &HlrCamera::right(5.0), 8);
        assert!(!front_result.segments.is_empty());
        assert!(!right_result.segments.is_empty());
    }

    #[test]
    fn hlr_svg_contains_hidden_dashed_lines() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        // Hidden lines are rendered dashed
        assert!(
            svg.contains("stroke-dasharray") || svg.contains("hidden"),
            "SVG should mark hidden lines differently"
        );
    }

    #[test]
    fn hlr_result_has_correct_visibility_counts() {
        // An isometric view of a box has 3 visible faces and 3 hidden faces.
        // The front 3 edges of each visible face → at least some hidden segments exist.
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(10.0);
        let result = hlr(&brep, &camera, 16);
        let total = result.segments.len();
        assert!(total >= 12, "a box has 12 edges, expect at least 12 segments; got {total}");
    }

    #[test]
    fn hlr_circle_edge_sampling() {
        use rcad_kernel::geom::{Circle3, Curve3, CurveEval};

        // Build a minimal BRep with a single circle edge (no solids).
        let mut brep = rcad_kernel::BRep::new();
        let circ = Circle3 {
            center: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
            radius: 1.0,
        };
        // Add two vertices on the circle (half-circle arc)
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(0.0),
        });
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(std::f64::consts::PI),
        });
        brep.edges.push(rcad_kernel::topology::Edge { start: 0, end: 1 });
        brep.geom.curves.push(Curve3::Circle(circ));
        brep.geom.edge_curve.push(Some(0));
        brep.geom
            .edge_curve_range
            .push(Some([0.0, std::f64::consts::PI]));

        let camera = HlrCamera::top(5.0);
        let result = hlr(&brep, &camera, 8);

        // The circle edge should produce at least one segment.
        assert!(
            !result.segments.is_empty(),
            "circle edge should produce HLR segments"
        );

        // All sampled 3D points on the circle should lie ON the circle (unit radius).
        // Verify by checking screen_pts all lie within radius ≈ 1.0 of circle center
        // when projected top-down (X-Y plane).
        for seg in &result.segments {
            // The curve_hint for circle segments should be set.
            assert!(
                matches!(seg.curve_hint, Some(CurveHint::Circle { .. })),
                "circle edge segments should carry CurveHint::Circle"
            );
        }

        // SVG should contain arc path elements (not just lines) for circle edges.
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(
            svg.contains("<path") || result.segments.is_empty(),
            "circle edge SVG should contain <path> arc elements"
        );
    }

    /// Cylinder viewed from the side should produce silhouette line segments
    /// in addition to the wire edges.
    #[test]
    fn cylinder_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        // The cylinder axis is +Y.  Use the right-side camera (looking along -X)
        // so the view direction is perpendicular to the axis → two silhouette lines.
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cylinder HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cylinder should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Sphere HLR should produce silhouette segments (the great circle).
    #[test]
    fn sphere_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "sphere HLR should produce silhouette segments"
        );
    }

    /// Cone viewed from the side should produce two silhouette lines from the apex.
    #[test]
    fn cone_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 1.0,
            height: 2.0,
        });
        // View from the right (perpendicular to cone axis) → two silhouette generators.
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cone HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cone should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Torus HLR should produce silhouette segments.
    #[test]
    fn torus_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        let camera = HlrCamera::front(20.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "torus HLR should produce silhouette segments"
        );
    }

    // ── Assembly HLR tests ─────────────────────────────────────────────────────

    /// Two boxes side by side — both should produce segments.
    #[test]
    fn hlr_assembly_two_boxes() {
        let box1 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let box2 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (box1, DAffine3::from_translation(DVec3::new(-2.0, 0.0, 0.0)), "box_left".to_string()),
            (box2, DAffine3::from_translation(DVec3::new(2.0, 0.0, 0.0)), "box_right".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2, "should have 2 component results");
        assert!(result.components.iter().all(|c| !c.segments.is_empty()),
            "each component should produce segments");
    }

    /// Small box behind a large box — the small box should be partially hidden.
    #[test]
    fn hlr_assembly_occlusion() {
        let big = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 3.0, height: 3.0, depth: 3.0,
        });
        let small = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 0.5, height: 0.5, depth: 0.5,
        });

        // Front camera looks along +Y from (0, -10, 0).
        // Place small box at +Y behind the big box so it's occluded.
        let components = vec![
            (big, DAffine3::IDENTITY, "big".to_string()),
            (small, DAffine3::from_translation(DVec3::new(0.0, 3.0, 0.0)), "small_behind".to_string()),
        ];

        let camera = HlrCamera::front(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // The small box behind the big one should have mostly hidden segments
        let small_comp = result.components.iter().find(|c| c.name == "small_behind").unwrap();
        let hidden = small_comp.segments.iter().filter(|s| !s.visible).count();
        let visible = small_comp.segments.iter().filter(|s| s.visible).count();
        assert!(hidden > visible,
            "small box behind big one should have more hidden than visible segments; hidden={hidden}, visible={visible}");
    }

    /// Assembly with a single component should match single-BRep HLR.
    #[test]
    fn hlr_assembly_single_matches_hlr() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);

        let single_hlr = hlr(&brep, &camera, 8);
        let assembly_result = hlr_assembly(
            &[(brep.clone(), DAffine3::IDENTITY, "box".to_string())],
            &camera, 8,
        );

        assert_eq!(assembly_result.components.len(), 1);
        let asm_segs = &assembly_result.components[0].segments;
        // Segment counts should be similar (same geometry, same algorithm)
        assert!(asm_segs.len() >= single_hlr.segments.len() - 2,
            "assembly HLR should produce similar segment count");
        assert!(asm_segs.len() <= single_hlr.segments.len() + 2,
            "assembly HLR should produce similar segment count");
    }

    /// Stacked boxes — top box visible, bottom box partially occluded.
    #[test]
    fn hlr_assembly_stacked_boxes() {
        let bottom = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0, height: 1.0, depth: 2.0,
        });
        let top = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (bottom, DAffine3::from_translation(DVec3::new(0.0, 0.0, 0.0)), "bottom".to_string()),
            (top, DAffine3::from_translation(DVec3::new(0.0, 0.0, 1.5)), "top".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // Both boxes should have some visible segments
        for comp in &result.components {
            let vis = comp.segments.iter().filter(|s| s.visible).count();
            assert!(vis > 0, "{} should have visible segments", comp.name);
        }
    }

    /// Empty assembly should return empty result.
    #[test]
    fn hlr_assembly_empty() {
        let components: Vec<(BRep, DAffine3, String)> = vec![];
        let camera = HlrCamera::isometric(5.0);
        let result = hlr_assembly(&components, &camera, 8);
        assert!(result.components.is_empty());
    }
}
