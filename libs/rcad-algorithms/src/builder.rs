use std::collections::HashMap;

use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topology::*;
use rcad_kernel::BRep;

use crate::bopds::ds::*;
use crate::classify::{classify_point, Classification};
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
struct SubFace {
    /// Boundary vertex positions in 3D (ordered polygon).
    boundary: Vec<DVec3>,
    /// The surface this lies on.
    surface: Surface3,
    /// Normal direction.
    normal: DVec3,
}

impl SubFace {
    fn sample_point(&self) -> DVec3 {
        // Use centroid offset slightly inward along the normal.
        // The offset avoids the sample sitting exactly on a face of the other
        // solid (which causes ambiguous IN/OUT classification).
        let centroid =
            self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64;
        centroid + self.normal * TOLERANCE_ABS * 10.0
    }
}

/// Builds result BRep, deduplicating vertices and edges.
struct ResultBuilder {
    vertices: Vec<DVec3>,
    vertex_map: HashMap<u64, usize>, // hash of position → index
    edges: Vec<(usize, usize)>,
    faces: Vec<(Vec<usize>, Vec<[usize; 3]>, DVec3, Surface3)>, // (boundary vertex indices, triangles, normal, surface)
}

impl ResultBuilder {
    fn new() -> Self {
        Self {
            vertices: Vec::new(),
            vertex_map: HashMap::new(),
            edges: Vec::new(),
            faces: Vec::new(),
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

    fn emit_face(&mut self, sub: &SubFace, flip: bool) {
        let normal = if flip { -sub.normal } else { sub.normal };

        // Add vertices
        let vert_indices: Vec<usize> = sub
            .boundary
            .iter()
            .map(|&p| self.add_vertex(p))
            .collect();

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
            .push((edge_indices, tris, normal, sub.surface));
    }

    fn build(self) -> BRep {
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
                edges: edge_indices,
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

        BRep {
            vertices,
            edges,
            solids: vec![Solid {
                shells: vec![Shell { faces }],
            }],
            geom,
        }
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
                    result.emit_face(sub, false);
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
                    result.emit_face(sub, flip);
                }
            }
        }

        let brep = result.build();
        if brep.solids[0].shells[0].faces.is_empty() {
            return Err(BooleanError::DegenerateResult);
        }

        Ok(brep)
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
                surface: face.surface,
                normal: face.normal,
            }];
        }

        // For planar faces: project to 2D, split by intersection segments
        match &face.surface {
            Surface3::Plane(plane) => self.split_planar_face(face_idx, plane),
            _ => {
                // Curved surfaces — return whole face for now
                let boundary = face
                    .boundary_verts
                    .iter()
                    .map(|&vi| self.ds.vertices[vi].point)
                    .collect();
                vec![SubFace {
                    boundary,
                    surface: face.surface,
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
                surface: face.surface,
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
                surface: face.surface,
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
    let (idx1, pt1) = crossings[0].clone();
    let (idx2, pt2) = crossings[1].clone();

    // Build two sub-polygons by walking the boundary
    let mut poly_a = Vec::new();
    let mut poly_b = Vec::new();

    // Walk from start to first crossing
    for i in 0..=idx1 {
        poly_a.push(poly[i]);
    }
    poly_a.push(pt1);
    poly_a.push(pt2);
    for i in (idx2 + 1)..n {
        poly_a.push(poly[i]);
    }

    // Walk from first crossing to second crossing
    poly_b.push(pt1);
    for i in (idx1 + 1)..=idx2 {
        poly_b.push(poly[i]);
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
