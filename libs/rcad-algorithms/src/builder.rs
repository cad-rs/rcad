use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval, *};
use rcad_kernel::topology::*;

use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::history::{BooleanHistory, FaceOrigin};
use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;
use crate::triangulate::triangulate_polygon;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    Union,
    Intersection,
    Difference,
}

#[derive(Debug)]
pub enum BooleanError {
    EmptyInput,
    MissingGeometry(&'static str),
    DegenerateResult,
    /// A numeric operation produced a non-finite or NaN value.
    NumericalFailure(&'static str),
    /// An expected non-empty collection was empty (e.g. polyline with no points).
    EmptyCollection(&'static str),
}

impl std::fmt::Display for BooleanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty input"),
            Self::MissingGeometry(msg) => write!(f, "missing geometry: {msg}"),
            Self::DegenerateResult => write!(f, "degenerate result"),
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {msg}"),
            Self::EmptyCollection(msg) => write!(f, "unexpected empty collection: {msg}"),
        }
    }
}

impl std::error::Error for BooleanError {}

/// A sub-region of an original face after splitting by intersection curves.
#[derive(Debug, Clone)]
pub struct SubFace {
    /// Boundary vertex positions in 3D (ordered polygon).
    pub boundary: Vec<DVec3>,
    /// The surface this lies on.
    pub surface: Surface3,
    /// Normal direction.
    pub normal: DVec3,
    /// UV centroid of this sub-face's parameter-space polygon (for curved surfaces).
    /// Used by `sample_point` to produce a geometrically representative interior point.
    pub uv_centroid: Option<DVec2>,
    /// Explicit override for the sample point. When set, `sample_point()` uses this
    /// instead of computing it from the boundary centroid. Used when the centroid would
    /// fall in a different classification region (e.g. the outer annular region around
    /// an embedded circle, whose centroid falls inside the circle).
    pub sample_override: Option<DVec3>,
}

impl SubFace {
    fn sample_point(&self) -> DVec3 {
        // Returns a point slightly INSIDE the surface (toward the interior of the solid),
        // so classify_point can tell whether this sub-face is inside or outside
        // the other solid.
        //
        // For sphere sub-faces the outward normal points AWAY from the sphere center,
        // so we must offset toward the center to stay inside the sphere's volume.
        // We use the UV centroid to get a point in the middle of the spherical cap.
        if let Some(pt) = self.sample_override {
            return pt;
        }
        match &self.surface {
            Surface3::Sphere(s) => {
                // Use UV centroid to pick a point in the CENTER of the spherical cap
                let surface_pt = if let Some(uv) = self.uv_centroid {
                    s.point_at(uv.x, uv.y)
                } else if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else {
                    s.center + s.radius * DVec3::X
                };
                // Offset inward toward sphere center
                let to_center = (s.center - surface_pt).normalize_or_zero();
                let inward = if to_center.length_squared() > 0.5 {
                    to_center
                } else {
                    -self.normal
                };
                surface_pt + inward * (TOLERANCE_ABS * 10.0)
            }
            _ => {
                let centroid = if self.boundary.is_empty() {
                    DVec3::ZERO
                } else {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                };
                centroid + self.normal * TOLERANCE_ABS * 10.0
            }
        }
    }
}

type FaceEntry = (Vec<usize>, Vec<[usize; 3]>, DVec3, Surface3);

/// Builds result BRep, deduplicating vertices and edges.
struct ResultBuilder {
    vertices: Vec<DVec3>,
    vertex_map: HashMap<u64, usize>, // hash of position → index
    edges: Vec<(usize, usize)>,
    faces: Vec<FaceEntry>, // (boundary vertex indices, triangles, normal, surface)
    face_origins: Vec<FaceOrigin>,
}

impl ResultBuilder {
    fn new() -> Self {
        Self {
            vertices: Vec::new(),
            vertex_map: HashMap::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            face_origins: Vec::new(),
        }
    }

    fn add_vertex(&mut self, point: DVec3) -> usize {
        let key = hash_point(point);
        if let Some(&idx) = self.vertex_map.get(&key) {
            // Double-check actual coincidence (hash collision protection)
            if points_coincide(self.vertices[idx], point) {
                return idx;
            }
        }
        // Linear scan fallback for hash collisions
        for (i, v) in self.vertices.iter().enumerate() {
            if points_coincide(*v, point) {
                return i;
            }
        }
        let idx = self.vertices.len();
        self.vertices.push(point);
        self.vertex_map.insert(key, idx);
        idx
    }

    fn add_edge(&mut self, v1: usize, v2: usize) -> usize {
        let key = (v1.min(v2), v1.max(v2));
        for (i, e) in self.edges.iter().enumerate() {
            if (e.0.min(e.1), e.0.max(e.1)) == key {
                return i;
            }
        }
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        idx
    }

    fn emit_face_with_origin(&mut self, sub: &SubFace, flip: bool, origin: FaceOrigin) {
        let normal = if flip { -sub.normal } else { sub.normal };

        // Add vertices
        let vert_indices: Vec<usize> = sub.boundary.iter().map(|&p| self.add_vertex(p)).collect();

        // Add edges
        let mut edge_indices = Vec::new();
        for i in 0..vert_indices.len() {
            let j = (i + 1) % vert_indices.len();
            let ei = self.add_edge(vert_indices[i], vert_indices[j]);
            edge_indices.push(ei);
        }

        // Triangulate
        let mut tris = triangulate_polygon(&sub.boundary, normal);
        // Remap triangle indices from local (0..n) to result vertex indices
        for tri in &mut tris {
            for idx in tri.iter_mut() {
                *idx = vert_indices[*idx];
            }
        }

        self.faces
            .push((edge_indices, tris, normal, sub.surface.clone()));
        self.face_origins.push(origin);
    }

    fn build(self) -> (BRep, BooleanHistory) {
        let vertices = self
            .vertices
            .into_iter()
            .map(|point| Vertex { point })
            .collect();

        let edges = self
            .edges
            .into_iter()
            .map(|(start, end)| Edge { start, end })
            .collect();

        let mut geom = rcad_kernel::GeomStore::default();
        let mut faces = Vec::new();

        for (edge_indices, triangles, normal, surface) in self.faces {
            let wire = Wire {
                edges: edge_indices.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
            };
            faces.push(Face {
                outer_wire: wire,
                inner_wires: vec![],
                normal,
                triangles,
            });

            let surf_idx = geom.surfaces.len();
            geom.surfaces.push(surface);
            geom.face_surface.push(Some(surf_idx));
        }

        let history = BooleanHistory {
            face_origins: self.face_origins,
        };

        let brep = BRep {
            vertices,
            edges,
            solids: vec![Solid {
                shells: vec![Shell { faces }],
            }],
            geom,
        };
        (brep, history)
    }
}

fn hash_point(p: DVec3) -> u64 {
    // Quantize to tolerance grid for spatial hashing
    let scale = 1.0 / TOLERANCE_ABS;
    let ix = (p.x * scale).round() as i64;
    let iy = (p.y * scale).round() as i64;
    let iz = (p.z * scale).round() as i64;
    // FNV-1a style hash
    let mut h: u64 = 14695981039346656037;
    for v in [ix, iy, iz] {
        h ^= v as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h
}

/// Boolean result builder (OCCT: BOPAlgo_BOP).
pub struct BooleanBuilder<'a> {
    ds: &'a DS,
    op: BooleanOpType,
}

impl<'a> BooleanBuilder<'a> {
    pub fn new(ds: &'a DS, op: BooleanOpType) -> Self {
        Self { ds, op }
    }

    pub fn build(&self) -> Result<BRep, BooleanError> {
        let (brep, _) = self.build_with_history()?;
        Ok(brep)
    }

    pub fn build_with_history(&self) -> Result<(BRep, BooleanHistory), BooleanError> {
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);

        if a_faces.is_empty() || b_faces.is_empty() {
            return Err(BooleanError::EmptyInput);
        }

        let mut result = ResultBuilder::new();

        // Process A faces against B solid
        for &fi in &a_faces {
            let sub_faces = self.split_face(fi);
            for sub in &sub_faces {
                let sample = sub.sample_point();
                let class = classify_point(sample, &b_faces, self.ds);

                let keep = match self.op {
                    BooleanOpType::Union => {
                        class == Classification::Out || class == Classification::On
                    }
                    BooleanOpType::Intersection => {
                        class == Classification::In || class == Classification::On
                    }
                    BooleanOpType::Difference => class == Classification::Out,
                };

                if keep {
                    result.emit_face_with_origin(sub, false, FaceOrigin::FromA(fi));
                }
            }
        }

        // Process B faces against A solid
        for &fi in &b_faces {
            let sub_faces = self.split_face(fi);
            for sub in &sub_faces {
                let sample = sub.sample_point();
                let class = classify_point(sample, &a_faces, self.ds);

                let keep = match self.op {
                    BooleanOpType::Union => class == Classification::Out,
                    BooleanOpType::Intersection => class == Classification::In,
                    BooleanOpType::Difference => class == Classification::In,
                };

                if keep {
                    let flip = self.op == BooleanOpType::Difference;
                    result.emit_face_with_origin(sub, flip, FaceOrigin::FromB(fi));
                }
            }
        }

        let (brep, history) = result.build();
        if brep.solids[0].shells[0].faces.is_empty() {
            return Err(BooleanError::DegenerateResult);
        }

        // Debug-mode geometry integrity check.
        // Verifies that every face in the result has a non-zero normal vector.
        // This catches the most common class of geometry regression (degenerate faces
        // produced by a wrong normal computation) without requiring a full wire-closure
        // check (which the current builder doesn't yet guarantee for all curve types).
        #[cfg(debug_assertions)]
        for (fi, face) in brep.solids[0].shells[0].faces.iter().enumerate() {
            debug_assert!(
                face.normal != glam::DVec3::ZERO,
                "boolean_op result face {fi} has zero normal"
            );
        }

        Ok((brep, history))
    }

    /// Split a face by intersection curves. If no intersection curves cross this
    /// face, returns the whole face as a single SubFace.
    fn split_face(&self, face_idx: usize) -> Vec<SubFace> {
        let face = &self.ds.faces[face_idx];
        let fi = &face.face_info;

        if fi.curves_in.is_empty() {
            // No intersections — return whole face
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| self.ds.vertices[vi].point)
                .collect();
            return vec![SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
                uv_centroid: None,
                sample_override: None,
            }];
        }

        // For planar faces: project to 2D, split by intersection segments
        match &face.surface.clone() {
            Surface3::Plane(plane) => self.split_planar_face(face_idx, plane),
            Surface3::Cylinder(_)
            | Surface3::Sphere(_)
            | Surface3::Cone(_)
            | Surface3::Torus(_) => self.split_curved_face_parametric(face_idx),
            _ => {
                // Other curved surfaces — return whole face for now
                let boundary = face
                    .boundary_verts
                    .iter()
                    .map(|&vi| self.ds.vertices[vi].point)
                    .collect();
                vec![SubFace {
                    boundary,
                    surface: face.surface.clone(),
                    normal: face.normal,
                    uv_centroid: None,
                sample_override: None,
                }]
            }
        }
    }

    /// Split a planar face by intersection line segments.
    ///
    /// Algorithm:
    /// 1. Project boundary + intersection segment endpoints to 2D
    /// 2. Find where intersection segment endpoints lie on boundary edges
    /// 3. Insert intersection points into boundary at correct positions
    /// 4. Walk augmented boundary to extract sub-polygons on each side
    fn split_planar_face(&self, face_idx: usize, plane: &Plane) -> Vec<SubFace> {
        let face = &self.ds.faces[face_idx];

        // Collect 3D boundary points
        let boundary_3d: Vec<DVec3> = face
            .boundary_verts
            .iter()
            .map(|&vi| self.ds.vertices[vi].point)
            .collect();

        if boundary_3d.len() < 3 {
            return vec![];
        }

        // Project boundary to 2D in the plane
        let (u_axis, v_axis) = plane_local_basis(plane);
        let project_to_2d = |p: DVec3| -> DVec2 {
            let d = p - plane.origin;
            DVec2::new(d.dot(u_axis), d.dot(v_axis))
        };
        let lift_to_3d = |uv: DVec2| -> DVec3 { plane.origin + u_axis * uv.x + v_axis * uv.y };

        let boundary_2d: Vec<DVec2> = boundary_3d.iter().map(|&p| project_to_2d(p)).collect();

        // Process each intersection curve to split the polygon
        let mut polygons_2d: Vec<Vec<DVec2>> = vec![boundary_2d];
        // Track circles that were embedded inside polygons (center_2d, radius).
        // When such a circle is fully inside a polygon, that polygon's centroid
        // may fall inside the circle — we must use a vertex-based sample instead.
        let mut embedded_circles: Vec<(DVec2, f64)> = Vec::new();

        for &ci in &face.face_info.curves_in {
            let ic = &self.ds.intersection_curves[ci];

            let curve_halfspace_split: Option<Vec<Vec<DVec2>>> = match &ic.curve {
                Curve3::Circle(circle) => {
                    // Plane-sphere intersection produces a circle lying in the plane.
                    // Project the circle center to 2D and split by the circle boundary.
                    let center_2d = project_to_2d(circle.center);
                    let radius = circle.radius;
                    let mut next: Vec<Vec<DVec2>> = Vec::new();
                    for poly in &polygons_2d {
                        let halves = split_polygon_by_circle_2d(poly, center_2d, radius);
                        next.extend(halves);
                    }
                    // Track this circle so we can compute correct sample points later
                    embedded_circles.push((center_2d, radius));
                    Some(next)
                }
                Curve3::Line(line) => {
                    // Use segment from start to end vertex
                    let p_start = self.ds.vertices[ic.start_vertex].point;
                    let p_end = self.ds.vertices[ic.end_vertex].point;
                    if points_coincide(p_start, p_end) {
                        None
                    } else {
                        let seg_s2d = project_to_2d(p_start);
                        let _seg_e2d = project_to_2d(p_end);
                        let mut next: Vec<Vec<DVec2>> = Vec::new();
                        for poly in &polygons_2d {
                            // Use line direction to split
                            let dir = DVec2::new(
                                (line.direction - plane.normal * line.direction.dot(plane.normal))
                                    .dot(u_axis),
                                (line.direction - plane.normal * line.direction.dot(plane.normal))
                                    .dot(v_axis),
                            );
                            let halves = split_polygon_2d_by_line(poly, seg_s2d, dir);
                            next.extend(halves);
                        }
                        Some(next)
                    }
                }
                _ => {
                    // For other curves, fall back to segment approach
                    let p_start = self.ds.vertices[ic.start_vertex].point;
                    let p_end = self.ds.vertices[ic.end_vertex].point;
                    if !points_coincide(p_start, p_end) {
                        let seg_s2d = project_to_2d(p_start);
                        let seg_e2d = project_to_2d(p_end);
                        let mut next: Vec<Vec<DVec2>> = Vec::new();
                        for poly in &polygons_2d {
                            let halves = split_polygon_2d_by_segment(poly, seg_s2d, seg_e2d);
                            next.extend(halves);
                        }
                        Some(next)
                    } else {
                        None
                    }
                }
            };

            if let Some(new_polys) = curve_halfspace_split
                && !new_polys.is_empty()
            {
                polygons_2d = new_polys;
            }
        }

        polygons_2d
            .into_iter()
            .filter(|p| p.len() >= 3)
            .map(|poly_2d| {
                let boundary: Vec<DVec3> = poly_2d.iter().map(|&uv| lift_to_3d(uv)).collect();
                // If there are embedded circles and this polygon's centroid falls inside
                // one of them, use the first boundary vertex (offset by normal) as the
                // sample point instead. All polygon vertices of the outer region are
                // outside all embedded circles, so the first vertex is a valid sample.
                let sample_override = if !embedded_circles.is_empty() {
                    let centroid_2d = {
                        let sum = poly_2d.iter().fold(DVec2::ZERO, |acc, &p| acc + p);
                        sum / poly_2d.len() as f64
                    };
                    let centroid_in_circle = embedded_circles.iter().any(|&(c, r)| {
                        (centroid_2d - c).length() < r
                    });
                    if centroid_in_circle && !boundary.is_empty() {
                        // Pick first vertex (outside the circle) + normal offset
                        Some(boundary[0] + face.normal * crate::tolerance::TOLERANCE_ABS * 10.0)
                    } else {
                        None
                    }
                } else {
                    None
                };
                SubFace {
                    boundary,
                    surface: face.surface.clone(),
                    normal: face.normal,
                    uv_centroid: None,
                    sample_override,
                }
            })
            .collect()
    }

    fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }

    /// Split a curved face (Cylinder, Sphere, Cone, Torus) by intersection polylines.
    ///
    /// Legacy approximate method: for each intersection polyline that crosses the face,
    /// we split the boundary point list into two halves at the points closest to the
    /// polyline endpoints. Kept as fallback when UV data or PCurves are unavailable.
    fn split_curved_face_legacy(&self, face_idx: usize) -> Vec<SubFace> {
        let face = &self.ds.faces[face_idx];
        let surface = face.surface.clone();
        let normal = face.normal;

        // Collect all intersection polylines for this face
        let mut all_polylines: Vec<Vec<DVec3>> = Vec::new();
        for &ci in &face.face_info.curves_in {
            let ic = &self.ds.intersection_curves[ci];
            if ic.polyline.len() >= 2 {
                all_polylines.push(ic.polyline.clone());
            } else {
                // Analytic curve — sample it into a polyline (e.g. circle)
                let pts: Vec<DVec3> = (0..=16)
                    .map(|i| {
                        let t = ic.t_range[0] + (ic.t_range[1] - ic.t_range[0]) * i as f64 / 16.0;
                        use rcad_kernel::CurveEval;
                        ic.curve.point_at(t)
                    })
                    .collect();
                all_polylines.push(pts);
            }
        }

        if all_polylines.is_empty() {
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| self.ds.vertices[vi].point)
                .collect();
            return vec![SubFace {
                boundary,
                surface,
                normal,
                uv_centroid: None,
                sample_override: None,
            }];
        }

        // Collect boundary vertices
        let boundary_pts: Vec<DVec3> = face
            .boundary_verts
            .iter()
            .map(|&vi| self.ds.vertices[vi].point)
            .collect();

        if boundary_pts.len() < 3 {
            return vec![SubFace {
                boundary: boundary_pts,
                surface,
                normal,
                uv_centroid: None,
                sample_override: None,
            }];
        }

        // For each intersection polyline, split the boundary into two sub-faces
        // by finding the boundary points closest to each polyline endpoint.
        let mut result_boundaries: Vec<Vec<DVec3>> = vec![boundary_pts];

        for polyline in &all_polylines {
            let (Some(&seg_start), Some(&seg_end)) = (polyline.first(), polyline.last()) else {
                continue;
            };

            let mut next_result: Vec<Vec<DVec3>> = Vec::new();
            for bnd in result_boundaries.drain(..) {
                let n = bnd.len();
                if n < 3 {
                    next_result.push(bnd);
                    continue;
                }

                // Find indices of boundary points closest to the two polyline endpoints
                let Some((i_start, _)) = bnd
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        a.distance_squared(seg_start)
                            .partial_cmp(&b.distance_squared(seg_start))
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                else {
                    next_result.push(bnd);
                    continue;
                };
                let Some((i_end, _)) = bnd
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        a.distance_squared(seg_end)
                            .partial_cmp(&b.distance_squared(seg_end))
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                else {
                    next_result.push(bnd);
                    continue;
                };

                // Ensure i_start < i_end for consistent splitting
                let (ia, ib, p_a, p_b) = if i_start <= i_end {
                    (i_start, i_end, seg_start, seg_end)
                } else {
                    (i_end, i_start, seg_end, seg_start)
                };

                if ia == ib {
                    // Degenerate: can't split, keep as is
                    next_result.push(bnd);
                    continue;
                }

                // Sub-face A: bnd[0..=ia] + polyline + bnd[ib..=n-1]
                let mut sub_a: Vec<DVec3> = bnd[..=ia].to_vec();
                sub_a.push(p_a);
                for &p in polyline.iter().skip(1).rev().skip(1) {
                    sub_a.push(p);
                }
                sub_a.push(p_b);
                sub_a.extend_from_slice(&bnd[ib..]);

                // Sub-face B: bnd[ia..=ib] + reverse polyline
                let mut sub_b: Vec<DVec3> = bnd[ia..=ib].to_vec();
                sub_b.push(p_b);
                for &p in polyline.iter().skip(1).rev().skip(1) {
                    sub_b.push(p);
                }
                sub_b.push(p_a);

                if sub_a.len() >= 3 {
                    next_result.push(sub_a);
                }
                if sub_b.len() >= 3 {
                    next_result.push(sub_b);
                }
            }
            result_boundaries = next_result;
        }

        result_boundaries
            .into_iter()
            .filter(|b| b.len() >= 3)
            .map(|boundary| SubFace {
                boundary,
                surface: surface.clone(),
                normal,
                uv_centroid: None,
                sample_override: None,
            })
            .collect()
    }

    /// Split a curved face using parameter-space (UV) 2D clipping.
    ///
    /// For each intersection curve on this face, samples the associated PCurve
    /// into a 2D trim polyline in UV space, then splits the UV boundary polygon.
    /// Maps resulting sub-polygons back to 3D via surface evaluation.
    ///
    /// Falls back to `split_curved_face_legacy` when UV data or PCurves are missing.
    fn split_curved_face_parametric(&self, face_idx: usize) -> Vec<SubFace> {
        let face = &self.ds.faces[face_idx];

        // Need UV boundary to operate in parameter space
        let uv_boundary = match &face.uv_boundary {
            Some(b) if b.len() >= 3 => b.clone(),
            _ => return self.split_curved_face_legacy(face_idx),
        };

        let surface = face.surface.clone();
        let normal = face.normal;

        // Collect 2D trim polylines from PCurves for each intersection curve
        let mut trim_polylines: Vec<Vec<DVec2>> = Vec::new();
        for &ci in &face.face_info.curves_in {
            if let Some(pcurve) = self.find_pcurve_for_face(ci, face_idx) {
                let ic = &self.ds.intersection_curves[ci];
                let [t0, t1] = ic.t_range;
                const N: usize = 32;
                let pts: Vec<DVec2> = match &pcurve {
                    // BSpline PCurves from polyline_pcurve_by_projection use
                    // chord-length parameterization normalized to [0,1].
                    // The 3D arc-length t_range is unrelated to the BSpline domain.
                    rcad_kernel::geom::Curve2d::BSpline(_) => (0..=N)
                        .map(|i| {
                            let t = i as f64 / N as f64;
                            pcurve.point_at(t)
                        })
                        .collect(),
                    // Analytic curves (Line2d, Circle2d, Ellipse2d) use the same
                    // t parameterization as the 3D intersection curve.
                    _ => (0..=N)
                        .map(|i| {
                            let t = t0 + (t1 - t0) * i as f64 / N as f64;
                            pcurve.point_at(t)
                        })
                        .collect(),
                };
                if pts.len() >= 2 {
                    trim_polylines.push(pts);
                }
            }
        }

        // If no PCurves available, fall back to legacy method
        if trim_polylines.is_empty() {
            return self.split_curved_face_legacy(face_idx);
        }

        // Split UV polygon by each trim polyline
        let mut uv_polygons: Vec<Vec<DVec2>> = vec![uv_boundary];

        for trim in &trim_polylines {
            let mut next: Vec<Vec<DVec2>> = Vec::new();
            for poly in uv_polygons.drain(..) {
                let halves = split_uv_polygon_by_trim(&poly, trim);
                next.extend(halves);
            }
            uv_polygons = next;
        }

        // Map each UV sub-polygon back to 3D
        uv_polygons
            .into_iter()
            .filter(|p| p.len() >= 3)
            .map(|uv_poly| {
                let n = uv_poly.len() as f64;
                let centroid_uv = uv_poly.iter().copied().sum::<DVec2>() / n;

                let boundary: Vec<DVec3> = match &surface {
                    Surface3::Sphere(_) => {
                        // Sphere poles: v=0 and v=π are degenerate (all u map to one point).
                        // Skip UV corners at v < ε or v > π-ε and supplement with trim points.
                        sphere_subface_boundary_3d(&uv_poly, &trim_polylines, &surface)
                    }
                    _ => uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect(),
                };
                // For curved surfaces, compute the actual surface normal at the centroid UV
                let sub_normal = {
                    let computed = surface.normal_at(centroid_uv.x, centroid_uv.y);
                    // If normal computation failed, fall back to face normal
                    if computed.length_squared() > 0.5 {
                        computed
                    } else {
                        normal
                    }
                };
                SubFace {
                    boundary,
                    surface: surface.clone(),
                    normal: sub_normal,
                    uv_centroid: Some(centroid_uv),
                    sample_override: None,
                }
            })
            .collect()
    }

    /// Find the PCurve (2D parametric curve) for the given intersection curve
    /// as it lies on the given face. Searches FaceFace interferences to determine
    /// whether this face is f1 (use pcurve_on_a) or f2 (use pcurve_on_b).
    fn find_pcurve_for_face(
        &self,
        curve_idx: usize,
        face_idx: usize,
    ) -> Option<rcad_kernel::geom::Curve2d> {
        for interference in &self.ds.interferences {
            if let Interference::FaceFace { f1, f2, curves, .. } = interference
                && curves.contains(&curve_idx)
            {
                let ic = &self.ds.intersection_curves[curve_idx];
                if *f1 == face_idx {
                    return ic.pcurve_on_a.clone();
                } else if *f2 == face_idx {
                    return ic.pcurve_on_b.clone();
                }
            }
        }
        None
    }
}

/// Build a 3D boundary polygon for a sphere sub-face given its UV polygon.
///
/// The sphere has singularities at v=0 (north pole) and v=π (south pole):
/// `point_at(u, 0)` and `point_at(u, π)` both return a single point regardless
/// of u, so UV polygon corners near the poles map to degenerate 3D vertices.
///
/// Strategy:
/// 1. Map each UV corner that is NOT near a pole to a 3D point normally.
/// 2. Supplement with intersection-curve (trim) 3D points that lie within
///    or near the boundary of this UV sub-polygon.
/// 3. Deduplicate and return at least 3 distinct 3D points.
fn sphere_subface_boundary_3d(
    uv_poly: &[DVec2],
    trim_polylines_uv: &[Vec<DVec2>],
    surface: &Surface3,
) -> Vec<DVec3> {
    use std::f64::consts::PI;
    const POLE_EPSILON: f64 = 0.05; // ~3° from pole

    let mut pts: Vec<DVec3> = Vec::new();

    // 1. Map non-degenerate UV corners to 3D
    for uv in uv_poly {
        let v = uv.y;
        if v > POLE_EPSILON && v < PI - POLE_EPSILON {
            pts.push(surface.point_at(uv.x, uv.y));
        }
        // At poles, collapse to a single well-defined 3D point
        else if v <= POLE_EPSILON {
            pts.push(surface.point_at(0.0, 0.0)); // north pole
        } else {
            pts.push(surface.point_at(0.0, PI)); // south pole
        }
    }

    // 2. Add trim polyline 3D points.
    //    Include trim UV points that are inside the sub-polygon OR very close to its boundary.
    //    This captures the intersection circle that forms the actual edge of the sub-face.
    for trim_uv in trim_polylines_uv {
        // Sample up to every-other trim point to avoid excessive density
        for (idx, uv) in trim_uv.iter().enumerate() {
            if idx % 2 != 0 { continue; } // sparse sampling
            let v = uv.y;
            if v > POLE_EPSILON && v < PI - POLE_EPSILON {
                let p3 = surface.point_at(uv.x, uv.y);
                // Include if inside or very near the boundary of the sub-polygon
                if point_in_polygon_2d(uv_poly, *uv) || point_near_polygon_2d(uv_poly, *uv, 0.1) {
                    pts.push(p3);
                }
            }
        }
    }

    // 3. Deduplicate (merge points within TOLERANCE_ABS)
    let mut deduped: Vec<DVec3> = Vec::new();
    for p in &pts {
        if deduped.iter().all(|q: &DVec3| (*p - *q).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS) {
            deduped.push(*p);
        }
    }

    deduped
}

/// Check if a 2D point is within `margin` of any edge of a polygon.
fn point_near_polygon_2d(poly: &[DVec2], pt: DVec2, margin: f64) -> bool {
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = poly[i];
        let b = poly[j];
        let ab = b - a;
        let len_sq = ab.length_squared();
        let t = if len_sq < 1e-14 { 0.0 } else { ((pt - a).dot(ab) / len_sq).clamp(0.0, 1.0) };
        let closest = a + t * ab;
        if (pt - closest).length() < margin {
            return true;
        }
    }
    false
}

/// Split a 2D UV polygon by a 2D trim polyline.
///
/// Algorithm:
/// 1. Find trim start/end's closest edge on the polygon boundary.
/// 2. Project trim endpoints onto boundary edges to find exact split points.
/// 3. Split polygon into two halves at those points, inserting the trim polyline
///    between them.
///
/// For closed trim polylines (start ≈ end), uses a closed-curve splitting
/// algorithm: the trim forms an interior polygon that divides the outer polygon
/// into "inside trim" and "outside trim" regions.
///
/// Returns 1 polygon if splitting is degenerate, or 2 sub-polygons otherwise.
fn split_uv_polygon_by_trim(poly: &[DVec2], trim: &[DVec2]) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 || trim.len() < 2 {
        return vec![poly.to_vec()];
    }

    let trim_start = trim[0];
    let trim_end = trim[trim.len() - 1];

    // Detect closed trim: start ≈ end in UV space
    let is_closed_trim = (trim_start - trim_end).length_squared() < 1e-6;
    if is_closed_trim {
        // The trim is a closed loop. Check if it's inside the polygon.
        // Use the trim as an interior boundary for "inside" and return
        // [trim_polygon, outer_polygon].
        let trim_centroid = trim.iter().copied().sum::<DVec2>() / trim.len() as f64;
        let is_inside = point_in_polygon_2d(poly, trim_centroid);
        if is_inside {
            // Return: [trim polygon (inner cap), outer polygon (without cap)]
            // Deduplicate closing point from trim if needed
            let trim_dedup: Vec<DVec2> = {
                let mut v = trim.to_vec();
                if v.len() > 1 && (v[0] - v[v.len() - 1]).length_squared() < 1e-12 {
                    v.pop();
                }
                v
            };
            if trim_dedup.len() >= 3 {
                return vec![trim_dedup, poly.to_vec()];
            }
        }
        // If trim centroid is not inside the polygon, can't split
        return vec![poly.to_vec()];
    }

    // Find closest point on each polygon edge for trim_start and trim_end,
    // returning (edge_index, parameter t in [0,1], projected point).
    let closest_on_boundary = |q: DVec2| -> (usize, f64, DVec2) {
        let mut best_edge = 0usize;
        let mut best_t = 0.0f64;
        let mut best_pt = poly[0];
        let mut best_dist = f64::INFINITY;

        for i in 0..n {
            let j = (i + 1) % n;
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            let len_sq = ab.dot(ab);
            let t = if len_sq < 1e-14 {
                0.0
            } else {
                ((q - a).dot(ab) / len_sq).clamp(0.0, 1.0)
            };
            let proj = a + t * ab;
            let dist = (q - proj).length_squared();
            if dist < best_dist {
                best_dist = dist;
                best_edge = i;
                best_t = t;
                best_pt = proj;
            }
        }
        (best_edge, best_t, best_pt)
    };

    let (edge_s, _t_s, pt_s) = closest_on_boundary(trim_start);
    let (edge_e, _t_e, pt_e) = closest_on_boundary(trim_end);

    // Ensure ia <= ib for consistent polygon walking
    let (ia, ib, p_a, p_b, trim_forward) = if edge_s <= edge_e {
        (edge_s, edge_e, pt_s, pt_e, true)
    } else {
        (edge_e, edge_s, pt_e, pt_s, false)
    };

    if ia == ib {
        // Both endpoints project to the same edge — degenerate, can't split
        return vec![poly.to_vec()];
    }

    // Build the trim points in the correct order for each half
    let trim_pts: Vec<DVec2> = if trim_forward {
        trim.to_vec()
    } else {
        trim.iter().copied().rev().collect()
    };

    // Sub-polygon A: poly[0..=ia] + p_a + trim_pts (interior only) + p_b + poly[ib+1..]
    let mut sub_a: Vec<DVec2> = poly[..=ia].to_vec();
    sub_a.push(p_a);
    // Interior trim points (skip first and last which are endpoints)
    for &p in trim_pts.iter().skip(1).rev().skip(1) {
        sub_a.push(p);
    }
    sub_a.push(p_b);
    sub_a.extend_from_slice(&poly[ib + 1..]);

    // Sub-polygon B: p_a + poly[ia+1..=ib] + p_b + trim_pts reversed (interior only)
    let mut sub_b: Vec<DVec2> = vec![p_a];
    sub_b.extend_from_slice(&poly[ia + 1..=ib]);
    sub_b.push(p_b);
    for &p in trim_pts.iter().skip(1).rev().skip(1).rev() {
        sub_b.push(p);
    }

    // Deduplicate consecutive near-equal points
    let dedup_2d = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < 1e-18 {
            result.pop();
        }
        result
    };

    let sub_a = dedup_2d(sub_a);
    let sub_b = dedup_2d(sub_b);

    let mut out = Vec::new();
    if sub_a.len() >= 3 {
        out.push(sub_a);
    }
    if sub_b.len() >= 3 {
        out.push(sub_b);
    }

    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}

/// Split a 2D polygon by a circle boundary.
///
/// Vertices inside the circle (distance < radius) are on the "inside" group,
/// vertices outside (distance > radius) are on the "outside" group.
/// Returns up to 2 sub-polygons: the part inside and the part outside.
///
/// When the circle is fully inside the polygon (all vertices outside),
/// samples the circle at N_CIRCLE_SAMPLES points and returns both
/// the approximate circular cap and the annular region.
fn split_polygon_by_circle_2d(poly: &[DVec2], center: DVec2, radius: f64) -> Vec<Vec<DVec2>> {
    const N_CIRCLE_SAMPLES: usize = 24;
    let n = poly.len();
    if n < 3 {
        return vec![poly.to_vec()];
    }

    let tol = TOLERANCE_ABS;

    // Signed distance: negative = inside circle, positive = outside
    let signed_dist = |p: DVec2| -> f64 { (p - center).length() - radius };

    let dists: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();

    // Check if all vertices are on the same side
    let all_inside = dists.iter().all(|&d| d <= tol);
    let all_outside = dists.iter().all(|&d| d >= -tol);

    if all_inside {
        // All polygon vertices inside circle — keep whole polygon
        return vec![poly.to_vec()];
    }

    if all_outside {
        // Circle is fully inside the polygon OR polygon is fully outside circle.
        // Check if circle center is inside the polygon:
        let center_inside = point_in_polygon_2d(poly, center);
        if !center_inside {
            // Circle doesn't overlap with this polygon — keep as-is
            return vec![poly.to_vec()];
        }
        // Circle is fully inside the polygon — produce circular cap + annular region
        // Sample the circle at N points
        let circle_poly: Vec<DVec2> = (0..N_CIRCLE_SAMPLES)
            .map(|i| {
                let theta = std::f64::consts::TAU * i as f64 / N_CIRCLE_SAMPLES as f64;
                center + DVec2::new(theta.cos(), theta.sin()) * radius
            })
            .collect();
        // Return: inside = circle polygon, outside = original polygon (with circle as hole)
        // For simplicity, return just the circle as the "inside" part
        // and the original polygon as the "outside" part (approximate)
        return vec![circle_poly, poly.to_vec()];
    }

    // Find crossings: edges where signed distance changes sign
    let mut crossings: Vec<(usize, DVec2)> = Vec::new();

    for i in 0..n {
        let j = (i + 1) % n;
        let di = dists[i];
        let dj = dists[j];

        if di.abs() < tol {
            continue; // vertex i is on the circle
        }
        if dj.abs() < tol {
            continue; // vertex j is on the circle (handled when edge starting at j is processed)
        }

        if di * dj < 0.0 {
            // Edge crosses the circle boundary
            // Find exact crossing: solve |a + t*(b-a) - center|² = r²
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            let ac = a - center;
            let qa = ab.dot(ab);
            let qb = 2.0 * ab.dot(ac);
            let qc = ac.dot(ac) - radius * radius;
            let disc = qb * qb - 4.0 * qa * qc;
            if disc < 0.0 {
                continue;
            }
            let sq = disc.sqrt();
            for &sign in &[-1.0_f64, 1.0_f64] {
                let t = (-qb + sign * sq) / (2.0 * qa);
                if t > -tol && t < 1.0 + tol {
                    let t = t.clamp(0.0, 1.0);
                    let pt = a + t * ab;
                    crossings.push((i, pt));
                    break; // take the first valid crossing on this edge
                }
            }
        }
    }

    if crossings.len() < 2 {
        // Can't split — keep as-is
        return vec![poly.to_vec()];
    }

    // Sort crossings by edge index
    crossings.sort_by_key(|(idx, _)| *idx);

    // Take the first two crossings
    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];

    if idx1 == idx2 {
        return vec![poly.to_vec()];
    }

    // Sample the arc between pt1 and pt2 (going through the inside of the polygon)
    // Determine which arc (minor or major) connects pt1 to pt2 and stays inside the polygon
    let theta1 = (pt1 - center).to_angle();
    let theta2 = (pt2 - center).to_angle();

    // For the "inside" sub-polygon, we need the arc that passes through the inside of the polygon.
    // Try both arcs and pick the one whose midpoint is inside the polygon.
    let mid_theta_cw = (theta1 + theta2) * 0.5;
    let mid_theta_ccw = mid_theta_cw + std::f64::consts::PI;
    let mid_cw = center + DVec2::new(mid_theta_cw.cos(), mid_theta_cw.sin()) * radius;
    let _mid_ccw = center + DVec2::new(mid_theta_ccw.cos(), mid_theta_ccw.sin()) * radius;

    // The arc midpoint that is inside the polygon corresponds to the "inside" portion
    let arc_goes_cw_inside = point_in_polygon_2d(poly, mid_cw);
    let inner_mid_theta = if arc_goes_cw_inside {
        mid_theta_cw
    } else {
        mid_theta_ccw
    };

    // Determine angular span and direction for the inner arc
    let arc_n = ((N_CIRCLE_SAMPLES as f64 * (theta2 - theta1).abs() / std::f64::consts::TAU)
        as usize)
        .max(3);

    // Build arc points from pt1 to pt2 going through inner_mid_theta
    let inner_arc: Vec<DVec2> = {
        // Compute proper arc from theta1 through inner_mid_theta to theta2
        let delta = {
            let mut d = theta2 - theta1;
            // Adjust delta to go through inner_mid_theta
            let going_ccw = inner_mid_theta > theta1 || inner_mid_theta < theta2;
            if going_ccw {
                while d < 0.0 {
                    d += std::f64::consts::TAU;
                }
                if d > std::f64::consts::TAU {
                    d -= std::f64::consts::TAU;
                }
            } else {
                while d > 0.0 {
                    d -= std::f64::consts::TAU;
                }
                if d < -std::f64::consts::TAU {
                    d += std::f64::consts::TAU;
                }
            }
            d
        };
        (0..=arc_n)
            .map(|i| {
                let t = i as f64 / arc_n as f64;
                let theta = theta1 + delta * t;
                center + DVec2::new(theta.cos(), theta.sin()) * radius
            })
            .collect()
    };

    // Sub-polygon "inside" (circle side): pt1 → arc → pt2 + polygon walk from idx2 to idx1
    // Actually: vertices of polygon that are INSIDE the circle + arc from pt1 to pt2
    let poly_inside_verts: Vec<DVec2> = poly[idx1 + 1..=idx2].to_vec();

    let mut sub_inside: Vec<DVec2> = vec![pt1];
    sub_inside.extend_from_slice(&poly_inside_verts);
    sub_inside.push(pt2);
    // Add arc back (reversed, so the boundary goes: inside polygon verts, then arc back to pt1)
    for &p in inner_arc.iter().skip(1).rev().skip(1) {
        sub_inside.push(p);
    }

    // Sub-polygon "outside" (non-circle side): pt2 → arc → pt1 + polygon walk
    let poly_outside_verts_a: Vec<DVec2> = poly[..=idx1].to_vec();
    let poly_outside_verts_b: Vec<DVec2> = poly[idx2 + 1..].to_vec();

    let mut sub_outside: Vec<DVec2> = poly_outside_verts_a;
    sub_outside.push(pt1);
    // Add inner arc (forward) as the "hole" boundary
    for &p in inner_arc.iter().skip(1).rev().skip(1) {
        sub_outside.push(p);
    }
    sub_outside.push(pt2);
    sub_outside.extend(poly_outside_verts_b);

    let dedup = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < 1e-18 {
            result.pop();
        }
        result
    };

    let sub_inside = dedup(sub_inside);
    let sub_outside = dedup(sub_outside);

    let mut out = Vec::new();
    if sub_inside.len() >= 3 {
        out.push(sub_inside);
    }
    if sub_outside.len() >= 3 {
        out.push(sub_outside);
    }

    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}

/// Check if a 2D point is inside a 2D polygon using ray casting.
fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i];
        let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Split a 2D polygon by an infinite line through `point` with direction `dir`.
///
/// Vertices on the positive side (cross product > 0) form one group, negative side the other.
fn split_polygon_2d_by_line(poly: &[DVec2], point: DVec2, dir: DVec2) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 {
        return vec![poly.to_vec()];
    }
    let tol = TOLERANCE_ABS;

    // Signed distance from line
    let signed_dist = |p: DVec2| -> f64 {
        let d = p - point;
        dir.x * d.y - dir.y * d.x // perpendicular component
    };

    let dists: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();
    let all_pos = dists.iter().all(|&d| d >= -tol);
    let all_neg = dists.iter().all(|&d| d <= tol);

    if all_pos || all_neg {
        return vec![poly.to_vec()];
    }

    let mut crossings: Vec<(usize, DVec2)> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let di = dists[i];
        let dj = dists[j];
        if di.abs() < tol || dj.abs() < tol {
            continue;
        }
        if di * dj < 0.0 {
            let t = di / (di - dj);
            let pt = poly[i] + t * (poly[j] - poly[i]);
            crossings.push((i, pt));
        }
    }

    if crossings.len() < 2 {
        return vec![poly.to_vec()];
    }

    crossings.sort_by_key(|(idx, _)| *idx);

    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];
    if idx1 == idx2 {
        return vec![poly.to_vec()];
    }

    let mut sub_a: Vec<DVec2> = poly[..=idx1].to_vec();
    sub_a.push(pt1);
    sub_a.push(pt2);
    sub_a.extend_from_slice(&poly[idx2 + 1..]);

    let mut sub_b: Vec<DVec2> = vec![pt1];
    sub_b.extend_from_slice(&poly[idx1 + 1..=idx2]);
    sub_b.push(pt2);

    let dedup = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < 1e-18 {
            result.pop();
        }
        result
    };

    let sub_a = dedup(sub_a);
    let sub_b = dedup(sub_b);
    let mut out = Vec::new();
    if sub_a.len() >= 3 {
        out.push(sub_a);
    }
    if sub_b.len() >= 3 {
        out.push(sub_b);
    }
    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}

/// Split a 2D polygon by a segment from `seg_start` to `seg_end`.
fn split_polygon_2d_by_segment(
    poly: &[DVec2],
    seg_start: DVec2,
    seg_end: DVec2,
) -> Vec<Vec<DVec2>> {
    let dir = seg_end - seg_start;
    if dir.length_squared() < 1e-18 {
        return vec![poly.to_vec()];
    }
    split_polygon_2d_by_line(poly, seg_start, dir.normalize())
}
