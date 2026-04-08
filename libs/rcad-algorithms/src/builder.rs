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
}

impl std::fmt::Display for BooleanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty input"),
            Self::MissingGeometry(msg) => write!(f, "missing geometry: {msg}"),
            Self::DegenerateResult => write!(f, "degenerate result"),
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
}

impl SubFace {
    fn sample_point(&self) -> DVec3 {
        // Use centroid offset slightly inward along the normal.
        // The offset avoids the sample sitting exactly on a face of the other
        // solid (which causes ambiguous IN/OUT classification).
        let centroid = self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64;
        centroid + self.normal * TOLERANCE_ABS * 10.0
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

        // Collect all intersection curve segments that cross this face
        let mut segments: Vec<(DVec3, DVec3)> = Vec::new();
        for &ci in &face.face_info.curves_in {
            let ic = &self.ds.intersection_curves[ci];
            let p_start = self.ds.vertices[ic.start_vertex].point;
            let p_end = self.ds.vertices[ic.end_vertex].point;
            segments.push((p_start, p_end));
        }

        if segments.is_empty() {
            return vec![SubFace {
                boundary: boundary_3d,
                surface: face.surface.clone(),
                normal: face.normal,
            }];
        }

        // For each intersection segment, split the polygon into two halves.
        // For box-box cases, each face gets at most one intersection segment,
        // producing at most two sub-faces.
        let mut polygons = vec![boundary_3d];

        for (seg_start, seg_end) in &segments {
            let mut new_polygons = Vec::new();

            for poly in &polygons {
                let split = split_polygon_by_segment(poly, *seg_start, *seg_end, plane);
                new_polygons.extend(split);
            }

            polygons = new_polygons;
        }

        polygons
            .into_iter()
            .filter(|p| p.len() >= 3)
            .map(|boundary| SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
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
            }];
        }

        // For each intersection polyline, split the boundary into two sub-faces
        // by finding the boundary points closest to each polyline endpoint.
        let mut result_boundaries: Vec<Vec<DVec3>> = vec![boundary_pts];

        for polyline in &all_polylines {
            let seg_start = *polyline.first().unwrap();
            let seg_end = *polyline.last().unwrap();

            let mut next_result: Vec<Vec<DVec3>> = Vec::new();
            for bnd in result_boundaries.drain(..) {
                let n = bnd.len();
                if n < 3 {
                    next_result.push(bnd);
                    continue;
                }

                // Find indices of boundary points closest to the two polyline endpoints
                let (i_start, _) = bnd
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        a.distance_squared(seg_start)
                            .partial_cmp(&b.distance_squared(seg_start))
                            .unwrap()
                    })
                    .unwrap();
                let (i_end, _) = bnd
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        a.distance_squared(seg_end)
                            .partial_cmp(&b.distance_squared(seg_end))
                            .unwrap()
                    })
                    .unwrap();

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
                let [t0, t1] = {
                    let ic = &self.ds.intersection_curves[ci];
                    ic.t_range
                };
                const N: usize = 32;
                let pts: Vec<DVec2> = (0..=N)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * i as f64 / N as f64;
                        pcurve.point_at(t)
                    })
                    .collect();
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
                let boundary: Vec<DVec3> = uv_poly
                    .iter()
                    .map(|uv| surface.point_at(uv.x, uv.y))
                    .collect();
                SubFace {
                    boundary,
                    surface: surface.clone(),
                    normal,
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

/// Split a 3D polygon by a line segment. Returns the resulting sub-polygons
/// (1 if segment doesn't cross, 2 if it splits the polygon).
fn split_polygon_by_segment(
    poly: &[DVec3],
    seg_start: DVec3,
    seg_end: DVec3,
    plane: &Plane,
) -> Vec<Vec<DVec3>> {
    let n = poly.len();
    let (u_axis, v_axis) = plane_local_basis(plane);

    let project = |p: DVec3| -> [f64; 2] {
        let d = p - plane.origin;
        [d.dot(u_axis), d.dot(v_axis)]
    };

    let seg_s = project(seg_start);
    let seg_e = project(seg_end);
    let seg_dir_2d = [seg_e[0] - seg_s[0], seg_e[1] - seg_s[1]];

    // Classify each polygon vertex by which side of the segment it's on
    let signed_dist = |p: [f64; 2]| -> f64 {
        // Cross product of seg_dir with (p - seg_s)
        seg_dir_2d[0] * (p[1] - seg_s[1]) - seg_dir_2d[1] * (p[0] - seg_s[0])
    };

    let poly_2d: Vec<[f64; 2]> = poly.iter().map(|&p| project(p)).collect();
    let sides: Vec<f64> = poly_2d.iter().map(|p| signed_dist(*p)).collect();

    // Find edges that cross the segment line
    let mut crossings: Vec<(usize, DVec3)> = Vec::new(); // (edge_index, intersection_point)

    for i in 0..n {
        let j = (i + 1) % n;
        let si = sides[i];
        let sj = sides[j];

        if si.abs() < TOLERANCE_ABS {
            // Vertex i is on the line — don't double-count
            continue;
        }
        if sj.abs() < TOLERANCE_ABS {
            // Vertex j is on the line — will be handled when processing edge starting at j
            continue;
        }

        if si * sj < 0.0 {
            // Edge crosses the line
            let t = si / (si - sj);
            let p = poly[i] + (poly[j] - poly[i]) * t;
            crossings.push((i, p));
        }
    }

    if crossings.len() < 2 {
        // Segment doesn't properly split this polygon
        return vec![poly.to_vec()];
    }

    // Sort crossings by position along the polygon boundary
    crossings.sort_by_key(|(idx, _)| *idx);

    // Take the first two crossings to split the polygon
    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];

    // Build two sub-polygons by walking the boundary
    let mut poly_a = Vec::new();
    let mut poly_b = Vec::new();

    // Walk from start to first crossing
    for &p in poly.iter().take(idx1 + 1) {
        poly_a.push(p);
    }
    poly_a.push(pt1);
    poly_a.push(pt2);
    for &p in poly.iter().skip(idx2 + 1) {
        poly_a.push(p);
    }

    // Walk from first crossing to second crossing
    poly_b.push(pt1);
    for &p in poly.iter().skip(idx1 + 1).take(idx2 - idx1) {
        poly_b.push(p);
    }
    poly_b.push(pt2);

    // Remove near-duplicate consecutive vertices
    let dedup = |v: Vec<DVec3>| -> Vec<DVec3> {
        let mut result: Vec<DVec3> = Vec::new();
        for p in v {
            if result.is_empty() || !points_coincide(*result.last().unwrap(), p) {
                result.push(p);
            }
        }
        // Check wrap-around
        if result.len() > 1 && points_coincide(*result.first().unwrap(), *result.last().unwrap()) {
            result.pop();
        }
        result
    };

    let poly_a = dedup(poly_a);
    let poly_b = dedup(poly_b);

    let mut result = Vec::new();
    if poly_a.len() >= 3 {
        result.push(poly_a);
    }
    if poly_b.len() >= 3 {
        result.push(poly_b);
    }

    if result.is_empty() {
        vec![poly.to_vec()]
    } else {
        result
    }
}

/// Split a 2D UV polygon by a 2D trim polyline.
///
/// Algorithm:
/// 1. Find trim start/end's closest edge on the polygon boundary.
/// 2. Project trim endpoints onto boundary edges to find exact split points.
/// 3. Split polygon into two halves at those points, inserting the trim polyline
///    between them.
///
/// Returns 1 polygon if splitting is degenerate, or 2 sub-polygons otherwise.
fn split_uv_polygon_by_trim(poly: &[DVec2], trim: &[DVec2]) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 || trim.len() < 2 {
        return vec![poly.to_vec()];
    }

    let trim_start = *trim.first().unwrap();
    let trim_end = *trim.last().unwrap();

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
            if result.is_empty() || (p - *result.last().unwrap()).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1
            && (result[0] - *result.last().unwrap()).length_squared() < 1e-18
        {
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
