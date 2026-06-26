use std::collections::{HashMap, HashSet};
use glam::DVec2; use glam::DVec3;
use rcad_kernel::geom::*; use rcad_kernel::BRep;
use rcad_kernel::topods;
use crate::history::{BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin};
use crate::bopds::ds::*; use crate::tolerance::*;
use crate::builder::types::{WireFace, WireSegment, WireEdgeSource, FaceEntry, FaceSampleData};
use crate::builder::SourceSide;
use crate::builder::{hash_point, curve_eq};
use crate::triangulate::{triangulate_polygon, triangulate_polygon_with_holes};
use rcad_kernel::topology::*;

/// Builds result BRep from accumulated DS face data.
///
/// OCCT-aligned: pure conversion — BuildResult does no dedup/merge/cull.
pub(crate) struct ResultBuilder {
    pub(crate) vertices: Vec<DVec3>,
    pub(crate) vertex_map: HashMap<u64, usize>,
    pub(crate) ds_vertex_map: HashMap<usize, usize>,
    pub(crate) edges: Vec<(usize, usize)>,
    pub(crate) faces: Vec<FaceEntry>,
    pub(crate) face_origins: Vec<FaceOrigin>,
    pub(crate) co_face_origins: Vec<(usize, FaceOrigin)>,
    pub(crate) shells: Vec<Vec<usize>>,
    pub(crate) solids: Vec<Vec<usize>>,
    pub(crate) custom_edge_curves: Vec<Option<Curve3>>,
    pub(crate) face_internal_vtx: Vec<Vec<usize>>,
    pub(crate) deg_edge_indices: HashSet<usize>,
    pub(crate) ic_edge_map: HashMap<usize, usize>,
    pub(crate) source_has_compound: bool,
    /// CompSolid solid groups — each entry is solid indices forming one CompSolid.
    /// Populated immediately by fill_images_containers_compsolid.
    pub(crate) compsolid_groups: Vec<Vec<usize>>,
}

impl ResultBuilder {
    fn estimate_boundary_normal(poly: &[DVec3]) -> DVec3 {
        if poly.len() < 3 {
            return DVec3::ZERO;
        }

        // Newell's method gives a stable polygon normal for arbitrary winding.
        let mut n = DVec3::ZERO;
        for i in 0..poly.len() {
            let p = poly[i];
            let q = poly[(i + 1) % poly.len()];
            n.x += (p.y - q.y) * (p.z + q.z);
            n.y += (p.z - q.z) * (p.x + q.x);
            n.z += (p.x - q.x) * (p.y + q.y);
        }
        let len = n.length();
        if len > TOLERANCE_LEN_MIN { n / len } else { DVec3::ZERO }
    }

    /// OCCT-aligned: emit BRep face from WireFace (replaces emit_face_with_origin).
    ///     Builds edges directly from WireSegments: seam edges use add_seam_edge /
    ///     add_edge_seam_degenerate; IC edges use add_circle_edge for Circle3 curves.
/// ✅ OCCT-aligned: emit_wire_face — builds BRep edges/face from WireSegments.
    pub(crate) fn emit_wire_face(
        &mut self,
        face_idx: usize,
        wf: &WireFace,
        segments: &[WireSegment],
        ds: &DS,
        flip: bool,
        origin: FaceOrigin,
        vertex_positions: &HashMap<usize, DVec3>,
    ) {
        let face = &ds.faces[face_idx];
        let mut normal = if flip { -face.normal } else { face.normal };
        if normal.length_squared() <= TOLERANCE_METRIC_SQ_NEAR_ZERO {
            normal = Self::estimate_boundary_normal_from_segments(&wf.outer_wire, segments, ds);
        }
        if normal.length_squared() <= TOLERANCE_METRIC_SQ_NEAR_ZERO {
            return;
        }

        // Outer wire: vertices + edges from WireSegments
        let mut vert_indices = Vec::new();
        let mut edge_indices = Vec::new();
        let ow: Vec<&usize> = wf.outer_wire.iter().filter(|&&si| segments[si].start_vertex != segments[si].end_vertex).collect();
        for &&si in &ow {
            let seg = &segments[si];
            // ✅ OCCT-aligned: canonical vertices use stored positions
            let get_pos = |vi: usize| -> DVec3 {
                vertex_positions.get(&vi).copied().unwrap_or(ds.vertices[vi].point)
            };
            let v1 = if seg.start_vertex < ds.vertices.len() {
                self.add_ds_vertex(seg.start_vertex, ds.vertices[seg.start_vertex].point)
            } else {
                self.add_vertex(vertex_positions[&seg.start_vertex])
            };
            let v2 = if seg.end_vertex < ds.vertices.len() {
                self.add_ds_vertex(seg.end_vertex, ds.vertices[seg.end_vertex].point)
            } else {
                self.add_vertex(vertex_positions[&seg.end_vertex])
            };
            if vert_indices.is_empty() || vert_indices.last() != Some(&v1) {
                vert_indices.push(v1);
            }
            let (ei, forward) = if seg.is_seam {
                let seam_deg = (get_pos(seg.start_vertex)
                    - get_pos(seg.end_vertex)).length_squared() < TOLERANCE_ABS_SQ;
                let sphere_surf = match &ds.faces[face_idx].surface {
                    Surface3::Sphere(s) => s,
                    _ => &SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, ref_dir: DVec3::X },
                };
                // ✅ OCCT-aligned: canonical deg edges (vertex >= ds.vertices.len())
                let is_canon_deg = seg.start_vertex >= ds.vertices.len() || seg.end_vertex >= ds.vertices.len();
                let ei = if seam_deg || is_canon_deg {
                    self.add_edge_seam_degenerate(v1, v2, sphere_surf)
                } else {
                    let seam_normal = any_perpendicular(sphere_surf.axis).normalize();
                    let seam_circle = Curve3::Circle(Circle3 {
                        center: sphere_surf.center,
                        normal: seam_normal,
                        radius: sphere_surf.radius,
                    });                    self.add_seam_edge(v1, v2, seam_circle)
                };
                (ei, true)
            } else {
                let ei = match &seg.source {
                    // ✅ OCCT-aligned: IC edge identity (section edges shared).
                    WireEdgeSource::IntersectionCurve(ci) => {
                        let crv = &ds.intersection_curves[*ci].curve;
                        self.add_ic_edge(*ci, v1, v2, crv.clone())
                    }
                    WireEdgeSource::DsEdge(_) => self.add_edge(v1, v2),
                    _ => self.add_edge(v1, v2),
                };
                let forward = self.edges[ei].0 == v1;
                (ei, forward)
            };
            edge_indices.push((ei, forward));
        }

        let mut inner_wire_edges: Vec<Vec<(usize, bool)>> = Vec::new();
        let mut iw_vert_indices_all: Vec<usize> = Vec::new();
        for iw in &wf.inner_wires {
            let mut iw_verts = Vec::new();
            let mut iw_edges = Vec::new();
            for &si in iw {
                let seg = &segments[si];
                let getp = |vi: usize| -> DVec3 { if vi < ds.vertices.len() { ds.vertices[vi].point } else { *vertex_positions.get(&vi).unwrap_or(&DVec3::ZERO) } };
                let v1 = if seg.start_vertex < ds.vertices.len() { self.add_ds_vertex(seg.start_vertex, ds.vertices[seg.start_vertex].point) } else { let p = getp(seg.start_vertex); self.add_vertex(p) };
                let v2 = if seg.end_vertex < ds.vertices.len() { self.add_ds_vertex(seg.end_vertex, ds.vertices[seg.end_vertex].point) } else { let p = getp(seg.end_vertex); self.add_vertex(p) };
                if iw_verts.is_empty() || iw_verts.last() != Some(&v1) {
                    iw_verts.push(v1);
                }
                let ei = match &seg.source {
                    // ✅ OCCT-aligned: IC edge identity (inner/internal wires).
                    WireEdgeSource::IntersectionCurve(ci) => {
                        let crv = &ds.intersection_curves[*ci].curve;
                        self.add_ic_edge(*ci, v1, v2, crv.clone())
                    }
                    WireEdgeSource::DsEdge(_) | WireEdgeSource::SeamEdge => self.add_edge(v1, v2),
                };
                let forward = self.edges[ei].0 == v1;
                iw_edges.push((ei, forward));
            }
            inner_wire_edges.push(iw_edges);
            iw_vert_indices_all.extend(iw_verts);
        }

        // ✅ OCCT-aligned: Internal wire edges (TopAbs_INTERNAL).
        //    Seam edges use add_seam_edge for curve-aware unique identity.
        let mut internal_wire_edges: Vec<Vec<(usize, bool)>> = Vec::new();
        for iw in &wf.internal_wires {
            let mut iw_edges = Vec::new();
            for &si in iw {
                let seg = &segments[si];
                let getp = |vi: usize| -> DVec3 { if vi < ds.vertices.len() { ds.vertices[vi].point } else { *vertex_positions.get(&vi).unwrap_or(&DVec3::ZERO) } };
                let v1 = if seg.start_vertex < ds.vertices.len() { self.add_ds_vertex(seg.start_vertex, ds.vertices[seg.start_vertex].point) } else { let p = getp(seg.start_vertex); self.add_vertex(p) };
                let v2 = if seg.end_vertex < ds.vertices.len() { self.add_ds_vertex(seg.end_vertex, ds.vertices[seg.end_vertex].point) } else { let p = getp(seg.end_vertex); self.add_vertex(p) };
                let ei = match &seg.source {
                    WireEdgeSource::IntersectionCurve(ci) => {
                        let crv = &ds.intersection_curves[*ci].curve;
                        if let Curve3::Circle(_) = crv { self.add_circle_edge(v1, v2, crv.clone()) }
                        else { self.add_edge(v1, v2) }
                    }
                    WireEdgeSource::DsEdge(_) | WireEdgeSource::SeamEdge if seg.is_seam => {
                        let s = match &ds.faces[face_idx].surface {
                            Surface3::Sphere(sph) => sph,
                            _ => &SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, ref_dir: DVec3::X },
                        };
                        let c = Curve3::Circle(Circle3 { center: s.center, normal: any_perpendicular(s.axis).normalize(), radius: s.radius });
                        self.add_seam_edge(v1, v2, c)
                    }
                    _ => self.add_edge(v1, v2),
                };
                let forward = self.edges[ei].0 == v1;
                iw_edges.push((ei, forward));
            }
            internal_wire_edges.push(iw_edges);
        }

        // Triangulation
        let outer_boundary: Vec<DVec3> = vert_indices.iter().map(|&vi| self.vertices[vi]).collect();
        let iw_boundaries: Vec<Vec<DVec3>> = inner_wire_edges.iter().map(|iw_es| {
            let mut pts = Vec::new();
            for &(ei, _) in iw_es {
                let (a, b) = self.edges[ei];
                if pts.is_empty() || pts.last() != Some(&a) {
                    pts.push(a);
                }
            }
            pts.iter().map(|&vi| self.vertices[vi]).collect()
        }).collect();
        let all_vert_indices: Vec<usize> = [vert_indices.as_slice(), iw_vert_indices_all.as_slice()].concat();
        let mut tris = if iw_boundaries.is_empty() {
            triangulate_polygon(&outer_boundary, normal)
        } else {
            triangulate_polygon_with_holes(&outer_boundary, &iw_boundaries, normal)
        };
        for tri in &mut tris {
            for idx in tri.iter_mut() {
                *idx = all_vert_indices[*idx];
            }
        }

        // Coincident face dedup
        let centroid = outer_boundary.iter().copied().sum::<DVec3>() / outer_boundary.len().max(1) as f64;
        let area = Self::polygon_signed_area_on_normal(&outer_boundary, normal);
        let mut outer_sig: Vec<usize> = edge_indices.iter().map(|&(eid, _)| eid).collect();
        outer_sig.sort_unstable();
        let nlen = normal.length();
        let nunit = if nlen > TOLERANCE_LEN_MIN { normal / nlen } else { normal };
        for (existing_idx, (existing_outer, existing_inner, _existing_tris, existing_normal,
             _surf, _uv, existing_centroid, existing_area, _existing_sp, _existing_iw))
            in self.faces.iter().enumerate()
        {
            let mut ex_sig: Vec<usize> = existing_outer.iter().map(|&(eid, _)| eid).collect();
            for iw_edges in existing_inner {
                ex_sig.extend(iw_edges.iter().map(|&(eid, _)| eid));
            }
            ex_sig.sort_unstable();
            let elen = existing_normal.length();
            if elen <= TOLERANCE_LEN_MIN { continue; }
            let eunit = *existing_normal / elen;
            let sig_match = ex_sig == outer_sig;
            let geo_match = nunit.dot(eunit).abs() >= 0.99
                && (*existing_centroid - centroid).length() <= TOLERANCE_LINEAR_RELAX_8
                && (existing_area - area).abs() <= TOLERANCE_LINEAR_RELAX_8 * existing_area.max(area).max(1.0);
            if sig_match || geo_match {
                self.co_face_origins.push((existing_idx, origin));
                return;
            }
        }

        // ✅ OCCT-aligned: No extra internal vertices needed — wire pipeline handles
        //    seam edges via WireSegment virtual edges; BuilderFace does not add
        //    degenerate vertices to the result face.

        // Compute UV domain for sphere faces
        let sphere_uv = if matches!(face.surface, Surface3::Sphere(_)) {
            let uvs: Vec<DVec2> = if !wf.outer_wire.is_empty() {
                wf.outer_wire.iter().map(|&si| {
                    let seg = &segments[si];
                    let sph = match &face.surface {
                        Surface3::Sphere(s) => s,
                        _ => unreachable!(),
                    };
                    sph.world_to_uv(ds.vertices[seg.start_vertex].point)
                }).collect()
            } else { vec![] };
            if !uvs.is_empty() {
                let u_min = uvs.iter().map(|uv| uv.x).fold(f64::INFINITY, f64::min);
                let u_max = uvs.iter().map(|uv| uv.x).fold(f64::NEG_INFINITY, f64::max);
                let v_min = uvs.iter().map(|uv| uv.y).fold(f64::INFINITY, f64::min);
                let v_max = uvs.iter().map(|uv| uv.y).fold(f64::NEG_INFINITY, f64::max);
                if (u_max - u_min).abs() > TOLERANCE_FLOAT_LOOSE && (v_max - v_min).abs() > TOLERANCE_FLOAT_LOOSE {
                    Some([u_min, u_max, v_min, v_max])
                } else { None }
            } else { None }
        } else { None };

        self.face_internal_vtx.push(Vec::new());
        let sample_pt = if !wf.outer_wire.is_empty() {
            let si = wf.outer_wire[0];
            let seg = &segments[si];
            ds.vertices[seg.start_vertex].point
        } else {
            ds.vertices.get(0).map(|v| v.point).unwrap_or(DVec3::ZERO)
        };
        self.faces.push((
            edge_indices,
            inner_wire_edges,
            tris,
            normal,
            face.surface.clone(),
            sphere_uv,
            centroid,
            area,
            sample_pt,
            internal_wire_edges,
        ));
        self.face_origins.push(origin);
    }

    /// ✅ OCCT-aligned: estimate face normal from wire segments.
    ///     Uses Newell's method on the outer wire boundary vertices.
    fn estimate_boundary_normal_from_segments(
        outer_wire: &[usize],
        segments: &[WireSegment],
        ds: &DS,
    ) -> DVec3 {
        if outer_wire.len() < 3 { return DVec3::ZERO; }
        let pts: Vec<DVec3> = outer_wire.iter().map(|&si| {
            let seg = &segments[si];
            ds.vertices[seg.start_vertex].point
        }).collect();
        Self::estimate_boundary_normal(&pts)
    }

    fn polygon_signed_area_on_normal(poly: &[DVec3], normal: DVec3) -> f64 {
        if poly.len() < 3 {
            return 0.0;
        }
        let n = normal.normalize_or_zero();
        let ax = n.x.abs();
        let ay = n.y.abs();
        let az = n.z.abs();
        let axis = if ax >= ay && ax >= az {
            0usize
        } else if ay >= az {
            1usize
        } else {
            2usize
        };

        let mut area2 = 0.0;
        for i in 0..poly.len() {
            let p = poly[i];
            let q = poly[(i + 1) % poly.len()];
            area2 += match axis {
                0 => p.y * q.z - q.y * p.z,
                1 => p.x * q.z - q.x * p.z,
                _ => p.x * q.y - q.x * p.y,
            };
        }
        0.5 * area2.abs()
    }

    pub(crate) fn new() -> Self {
        Self {
            vertices: Vec::new(),
            vertex_map: HashMap::new(),
            ds_vertex_map: HashMap::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            face_origins: Vec::new(),
            co_face_origins: Vec::new(),
            custom_edge_curves: Vec::new(),
            face_internal_vtx: Vec::new(),
            deg_edge_indices: std::collections::HashSet::new(),
            ic_edge_map: HashMap::new(),
            shells: Vec::new(),
            solids: Vec::new(),
            source_has_compound: false,
            compsolid_groups: Vec::new(),
        }
    }

    /// ✅ OCCT-aligned: BuildResult(EDGE) — build edges from split_edges.
    /// ✅ OCCT-aligned: BuildResult(FACE) — build faces from accumulated face data.
    ///   OCCT Builder_1.cxx L130-168: iterate myImages for TopAbs_FACE, add to myShape.
    ///   rcad: build faces from self.faces, referencing already-built self.edges.
    ///   Maps each face's per-vertex-pair edges to the BRep edge indices from build_edges.
    /// ✅ OCCT-aligned: BuildResult(FACE) — build faces from accumulated face data.
    ///   OCCT Builder_1.cxx L130-168: iterate myImages for TopAbs_FACE, add to myShape.
    ///   rcad: validate face edge refs against built edges, prepare for shell/solid assembly.
    pub(crate) fn build_faces(&mut self) {
        // Validate that all face edge references are within bounds of built edges
        let n_edges = self.edges.len();
        for (fi, face) in self.faces.iter().enumerate() {
            for &(ei, _) in &face.0 {
                assert!(ei < n_edges,
                    "face[{}] edge ref {} out of range ({} edges)", fi, ei, n_edges);
            }
            for iw in &face.1 {
                for &(ei, _) in iw {
                    assert!(ei < n_edges,
                        "face[{}] inner edge ref {} out of range", fi, ei);
                }
            }
        }
    }

    /// ✅ OCCT-aligned: BuildResult(FACE) — add unmodified source face.
    /// ✅ OCCT-aligned: BuildResult(FACE) — add original source face (Builder_1.cxx L146-152).
    ///   OCCT adds the original TopoDS_Face regardless of surface type.
    ///   rcad: builds FaceEntry from DS boundary_edges + inner_boundary_edges.
    ///   Handles all surface types (Plane, Cylinder, Sphere, Cone, Torus).
    /// ✅ OCCT-aligned: BuildResult(FACE) — add original faces without images
    /// (Builder_1.cxx L146-152).  When myImages does NOT contain the face,
    /// OCCT adds the original TopoDS_Face as-is.  rcad: reconstruct the face
    /// from DS boundary edges, surface, and centroid.
    pub(crate) fn build_original_face(&mut self, ds: &DS, fi: usize, origin: FaceOrigin) {
        let face = &ds.faces[fi];

        // --- Outer wire from boundary_edges ---
        let mut edge_indices: Vec<(usize, bool)> = Vec::new();
        let mut prev_end: Option<usize> = None;
        for &ei in &face.boundary_edges {
            if ei >= ds.edges.len() { continue; }
            let e = &ds.edges[ei];
            let (sv, ev) = match prev_end {
                Some(pe) if e.start_vertex == pe => (e.start_vertex, e.end_vertex),
                Some(pe) if e.end_vertex == pe => (e.end_vertex, e.start_vertex),
                _ => (e.start_vertex, e.end_vertex),
            };
            let brep_sv = self.add_ds_vertex(sv, ds.vertices[sv].point);
            let brep_ev = self.add_ds_vertex(ev, ds.vertices[ev].point);
            let bei = self.add_edge(brep_sv, brep_ev);
            let fwd = (self.edges[bei].0, self.edges[bei].1) == (brep_sv, brep_ev);
            edge_indices.push((bei, fwd));
            prev_end = Some(ev);
        }
        if edge_indices.len() < 3 { return; }

        // --- Inner wires (holes) from inner_boundary_edges ---
        let mut inner_wires: Vec<Vec<(usize, bool)>> = Vec::new();
        for iw_edges in &face.inner_boundary_edges {
            let mut wire: Vec<(usize, bool)> = Vec::new();
            for &(ei, forward_in_ds) in iw_edges {
                if ei >= ds.edges.len() { continue; }
                let e = &ds.edges[ei];
                let (sv, ev) = if forward_in_ds {
                    (e.start_vertex, e.end_vertex)
                } else {
                    (e.end_vertex, e.start_vertex)
                };
                let brep_sv = self.add_ds_vertex(sv, ds.vertices[sv].point);
                let brep_ev = self.add_ds_vertex(ev, ds.vertices[ev].point);
                let bei = self.add_edge(brep_sv, brep_ev);
                let fwd = (self.edges[bei].0, self.edges[bei].1) == (brep_sv, brep_ev);
                wire.push((bei, fwd));
            }
            if wire.len() >= 2 {
                inner_wires.push(wire);
            }
        }

        let normal = face.normal;
        let surface = face.surface.clone();
        let centroid = edge_indices.iter()
            .filter_map(|&(ei, fwd)| {
                let e = self.edges.get(ei)?;
                self.vertices.get(if fwd { e.1 } else { e.0 }).copied()
            })
            .sum::<DVec3>() / edge_indices.len() as f64;
        self.faces.push((
            edge_indices, inner_wires, vec![], normal, surface, None, centroid, 0.0, centroid, vec![],
        ));
        self.face_origins.push(origin);
    }

    /// ✅ OCCT-aligned: BRep_Builder::MakeVertex — dedup by position.
    ///    OCCT: each TopoDS_Vertex is unique by TShape identity (may share position).
    ///    rcad: dedup by hash of position + linear scan (geometric, not identity-based).
    ///    Equivalent behavior: same-position vertices return same index.
    pub(crate) fn add_vertex(&mut self, point: DVec3) -> usize {
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

    /// ✅ OCCT-aligned: add vertex by DS index identity (TopoDS_Vertex TShape).
    pub(crate) fn add_ds_vertex(&mut self, ds_vi: usize, point: DVec3) -> usize {
        if let Some(&idx) = self.ds_vertex_map.get(&ds_vi) {
            return idx;
        }
        let idx = self.add_vertex(point);
        self.ds_vertex_map.insert(ds_vi, idx);
        idx
    }

    /// Geometric edge key for OCCT-aligned edge-set matching (BOPTools_Set analog).
    /// Returns a hash of the two quantized vertex positions, sorted for direction
    /// independence, so geometrically identical edges from different operands
    /// produce the same key regardless of traversal direction or edge index.
    fn edge_geo_key(&self, ei: usize) -> u64 {
        let (v1, v2) = self.edges[ei];
        let p1 = self.vertices[v1];
        let p2 = self.vertices[v2];
        // Quantize to 1e-4 grid (building-level tolerance, per OCCT Precision)
        let q = |v: f64| (v / 1e-4).round() as i64;
        let k1 = (q(p1.x), q(p1.y), q(p1.z));
        let k2 = (q(p2.x), q(p2.y), q(p2.z));
        // Sort for direction independence
        let (ka, kb) = if k1 < k2 { (k1, k2) } else { (k2, k1) };
        // FNV-1a hash of the two quantized tuples
        let mut h: u64 = 14695981039346656037;
        h ^= ka.0 as u64; h = h.wrapping_mul(1099511628211);
        h ^= ka.1 as u64; h = h.wrapping_mul(1099511628211);
        h ^= ka.2 as u64; h = h.wrapping_mul(1099511628211);
        h ^= kb.0 as u64; h = h.wrapping_mul(1099511628211);
        h ^= kb.1 as u64; h = h.wrapping_mul(1099511628211);
        h ^= kb.2 as u64; h = h.wrapping_mul(1099511628211);
        h
    }

    /// ✅ OCCT-aligned: BRep_Builder::MakeEdge — creates new unique edge.
    ///    OCCT: each TopoDS_Edge is a distinct entity (per TShape identity).
    ///    Even edges connecting the same vertices are distinct TopoDS_Edges.
    ///    rcad: always appends a new edge, same semantics.
    pub(crate) fn add_edge_occt(&mut self, v1: usize, v2: usize) -> usize {
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        idx
    }

    /// ✅ OCCT-aligned: BOPTools_AlgoTools::MakeSectEdge — shared section edge.
    ///    OCCT: MakeSectEdge creates ONE TopoDS_Edge that both intersecting faces
    ///    reference via BRep_Builder::Add (shared TShape identity).
    ///    rcad: maps intersection curve index → result edge index so both faces
    ///    emit_wire_face calls get the same edge index for the same IC curve.
    ///    OCCT: each TopoDS_Edge is a distinct handle — no post-hoc merge needed.
    pub(crate) fn add_ic_edge(&mut self, ici: usize, v1: usize, v2: usize, curve: Curve3) -> usize {
        if let Some(&idx) = self.ic_edge_map.get(&ici) {
            // OCCT-aligned: the edge must have same vertices for both faces.
            // If remap_ic_v produced different vertices for the same IC on
            // different faces, log a warning (indicates remap inconsistency).
            let existing = self.edges[idx];
            if (existing.0 != v1 || existing.1 != v2) && (existing.0 != v2 || existing.1 != v1) {
                eprintln!("[IC_VTX] ci={} existing=({}, {}) called=({}, {})", ici, existing.0, existing.1, v1, v2);
            }
            return idx;
        }
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        self.custom_edge_curves[idx] = Some(curve);
        self.ic_edge_map.insert(ici, idx);
        // no_merge_edges removed — edges are inherently unique by index
        idx
    }

    /// ✅ OCCT-aligned: BRep_Builder::Add edge sharing — dedup by (v1,v2) pair.
    ///    OCCT: BRep_Builder::Add(theSameEdge, faceA) then Add(theSameEdge, faceB)
    ///    shares the same TopoDS_Edge between faces (TShape identity).
    ///    rcad: add_edge(v1,v2) returns the same index for the same vertex pair,
    ///    achieving the same sharing without requiring TopoDS shape handles.
    pub(crate) fn add_edge(&mut self, v1: usize, v2: usize) -> usize {
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

    /// ✅ OCCT-aligned: create degenerate seam edge with hemisphere circle curve.
    ///    OCCT sphere face outer wire always has a degenerate seam edge (same vertex at both ends).
    ///    Adds a sphere horizontal circle curve to make the edge recognizable in STEP export.
    pub(crate) fn add_edge_seam_degenerate(&mut self, v1: usize, v2: usize, sphere_surf: &SphericalSurface) -> usize {
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        // Store seam circle curve for STEP writer
        // ✅ OCCT-aligned: seam edge = sphere meridian through pole (not IC circle).
        //    If normal = axis, it would coincide with plane-sphere IC causing curve merge errors.
        let seam_normal = any_perpendicular(sphere_surf.axis).normalize();
        let seam_circle = Curve3::Circle(Circle3 {
            center: sphere_surf.center,
            normal: seam_normal,
            radius: sphere_surf.radius,
        });
        self.custom_edge_curves[idx] = Some(seam_circle);
        self.deg_edge_indices.insert(idx);
        idx
    }

    /// ✅ OCCT-aligned: circle edge with curve-aware dedup.
    ///    OCCT: TopoDS_Edge identity is per-TShape, not per vertex pair.
    ///    Two edges sharing vertices but with different curves are distinct.
    ///    rcad: dedup by both (v1,v2) AND curve identity (Circle3 geometry).
    pub(crate) fn add_circle_edge(&mut self, v1: usize, v2: usize, circle: Curve3) -> usize {
        let key = (v1.min(v2), v1.max(v2));
        for (i, e) in self.edges.iter().enumerate() {
            if (e.0.min(e.1), e.0.max(e.1)) == key {
                if let Some(ref existing) = self.custom_edge_curves.get(i).and_then(|c| c.as_ref()) {
                    // Different curve at same vertex pair → distinct TopoDS_Edge
                    if !curve_eq(existing, &circle) {
                        let idx = self.add_edge_occt(v1, v2);
                        while self.custom_edge_curves.len() <= idx {
                            self.custom_edge_curves.push(None);
                        }
                        self.custom_edge_curves[idx] = Some(circle);
                        return idx;
                    }
                }
                // Same curve or no existing curve → reuse
                while self.custom_edge_curves.len() <= i {
                    self.custom_edge_curves.push(None);
                }
                self.custom_edge_curves[i] = Some(circle);
                return i;
            }
        }
        let idx = self.add_edge_occt(v1, v2);
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        self.custom_edge_curves[idx] = Some(circle);
        idx
    }
    /// ✅ OCCT-aligned: BOPTools_AlgoTools::MakeEdge 等价 -- 始终创建新边,不进行顶点去重。
    ///    使用 add_edge_occt,确保不被其他面的边合并。
    ///    适用于 seam 子段与 IC 弧在 OCCT 中是不同的 TopoDS_Edge。
    pub(crate) fn add_circle_edge_occt(&mut self, v1: usize, v2: usize, circle: Curve3) -> usize {
        let idx = self.add_edge_occt(v1, v2);
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        self.custom_edge_curves[idx] = Some(circle);
        // no_merge_edges removed — edges are inherently unique by index
        idx
    }


    /// ✅ OCCT-aligned: MakeEdge for seam edges (BRep_Builder::MakeEdge pattern).
    ///    OCCT: BRep_Builder::MakeEdge creates a TopoDS_Edge with the 3D curve.
    ///    Seam edges and IC arcs at the same vertex pair are distinct TopoDS_Edges
    ///    (different TShapes).  rcad: same vertex pair + same curve → reuse (shared
    ///    TShape); same vertex pair + different curve → create new via add_edge_occt
    ///    (distinct TShape).  This matches OCCT's per-TShape edge identity.
    pub(crate) fn add_seam_edge(&mut self, v1: usize, v2: usize, circle: Curve3) -> usize {
        // Same logic as add_circle_edge: check for existing edge with same
        // vertex pair but different curve → create new; same curve → reuse.
        let key = (v1.min(v2), v1.max(v2));
        for (i, e) in self.edges.iter().enumerate() {
            if (e.0.min(e.1), e.0.max(e.1)) == key {
                if let Some(ref existing) = self.custom_edge_curves.get(i).and_then(|c| c.as_ref()) {
                    if !curve_eq(existing, &circle) {
                        let idx = self.add_edge_occt(v1, v2);
                        while self.custom_edge_curves.len() <= idx {
                            self.custom_edge_curves.push(None);
                        }
                        self.custom_edge_curves[idx] = Some(circle);
                        // no_merge_edges removed — edges are inherently unique by index
                        return idx;
                    }
                }
                while self.custom_edge_curves.len() <= i {
                    self.custom_edge_curves.push(None);
                }
                self.custom_edge_curves[i] = Some(circle);
                return i;
            }
        }
        let idx = self.add_edge_occt(v1, v2);
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        self.custom_edge_curves[idx] = Some(circle);
        // no_merge_edges removed — edges are inherently unique by index
        idx
    }

    /// DEPRECATED (FaceSampleData 鍐呴儴): 鍦嗗姬鍐呰竟鐣屾娴?浠呭湪 split_planar_face 璺緞浣跨敤銆?
    ///    OCCT: MakeBlocks 鈫?BOPTools_AlgoTools::MakeEdge(aIC,...)
    ///    split_planar_face 鐢熸垚鐨勫唴杈圭晫鏈?28+鐐?绠€鍖栦负2绔偣(arc_simplify),
    ///    鐒跺悗 emit_face_with_origin 鐢?add_circle_edge 鍒涘缓绮剧‘ Circle3 杈广€?
    /// DEPRECATED (FaceSampleData 鍐呴儴): 鍦嗗姬澶栬竟鐣屸啋鍐呰竟鐣岃浆鎹€俉ireFace 涓嶉渶瑕佹姝ラ銆?
// SubFace removed: convert

// SubFace removed: find_inner


    /// ✅ OCCT-aligned: BuildResult — pure conversion from ResultBuilder arrays to BRep.
    ///   OCCT BuildResult (Builder_1.cxx L130-168) iterates myImages and adds shapes
    ///   to myShape.  rcad converts internal arrays (vertices, edges, faces) to BRep
    ///   topology (Vertex, Edge, Face).  Both do NO merge/cull/post-processing.

    /// Build result as a TopoDS-aligned BRep with shared TShape references.
    /// Each vertex/edge is a unique TShape; edges reference vertices via ShapeRef with Orientation.
    pub(crate) fn build_topods(mut self) -> (topods::BRep, BooleanHistory) {
        use topods::{Orientation, ShapeRef, TShape, TVertexData, TEdgeData, TWireData, TFaceData, TShellData, TSolidData};
        use std::sync::Arc;
        let mut t = topods::BRep::new();

        // 1. Vertices (each becomes a TShape::Vertex)
        for v in &self.vertices {
            t.add_tvertex(*v);
        }

        // 2. Edges (reference vertices by their TShape index)
        let mut e_map: Vec<ShapeRef> = Vec::with_capacity(self.edges.len());
        for (ei, &(start, end)) in self.edges.iter().enumerate() {
            let first = ShapeRef::new(start);
            let last = ShapeRef::new(end);
            // Store curve if available (needed by surface area sampler)
            let curve_idx = self.custom_edge_curves.get(ei).and_then(|c| c.as_ref()).map(|crv| {
                let ci = t.curves.len();
                t.curves.push(crv.clone());
                ci
            });
            let sr = t.add_tedge(curve_idx, first, last, [0.0, 1.0]);
            e_map.push(sr);
        }

        // Map from old flat face index to face ShapeRef
        let mut face_refs: Vec<ShapeRef> = Vec::with_capacity(self.faces.len());

        // 4. Faces — build wires, then face TShapes
        for (edge_indices, inner_wire_edges, _triangles, _normal, surface, _uv_domain, _centroid, _area, sample_point, internal_wire_edges) in self.faces {
            // Build outer wire
            let outer_edges: Vec<ShapeRef> = edge_indices.iter().map(|&(idx, forward)| {
                let orient = if forward { Orientation::Forward } else { Orientation::Reversed };
                if idx < e_map.len() {
                    ShapeRef::with_orientation(e_map[idx].index, orient)
                } else {
                    ShapeRef::with_orientation(idx, orient)
                }
            }).collect();
            let outer_wire = t.add_twire(outer_edges);

            // Build inner wires
            let mut inner_wires = Vec::new();
            for wire_idxs in inner_wire_edges {
                let iw_edges: Vec<ShapeRef> = wire_idxs.iter().map(|&(idx, forward)| {
                    let orient = if forward { Orientation::Forward } else { Orientation::Reversed };
                    if idx < e_map.len() {
                        ShapeRef::with_orientation(e_map[idx].index, orient)
                    } else {
                        ShapeRef::with_orientation(idx, orient)
                    }
                }).collect();
                if !iw_edges.is_empty() {
                    inner_wires.push(t.add_twire(iw_edges));
                }
            }

            // Add internal wire edges as inner wires
            for iw_edges in internal_wire_edges {
                let iw: Vec<ShapeRef> = iw_edges.iter().map(|&(idx, forward)| {
                    let orient = if forward { Orientation::Forward } else { Orientation::Reversed };
                    if idx < e_map.len() {
                        ShapeRef::with_orientation(e_map[idx].index, orient)
                    } else {
                        ShapeRef::with_orientation(idx, orient)
                    }
                }).collect();
                if iw.len() >= 2 {
                    inner_wires.push(t.add_twire(iw));
                }
            }

            let surf_idx = t.surfaces.len();
            t.surfaces.push(surface);
            let face_sr = t.add_tface(Some(surf_idx), outer_wire, inner_wires, Some(sample_point), _uv_domain);
            face_refs.push(face_sr);
        }

        if self.solids.is_empty() && self.shells.is_empty() {
            // Legacy path: single shell, single solid
            let shell = t.add_tshell(face_refs);
            t.add_tsolid(vec![shell]);
        } else if !self.compsolid_groups.is_empty() {
            // OCCT-aligned: CompSolid wraps multiple solids sharing boundary faces.
            for cs_group in &self.compsolid_groups {
                let mut solid_refs = Vec::new();
                for &si in cs_group {
                    if si >= self.solids.len() { continue; }
                    let shell_refs: Vec<ShapeRef> = self.solids[si].iter().map(|&shi| {
                        let sf = self.shells.get(shi).map_or(vec![], |sf| {
                            sf.iter().map(|&fi| face_refs.get(fi).copied().unwrap_or(ShapeRef::new(fi))).collect()
                        });
                        t.add_tshell(sf)
                    }).collect();
                    solid_refs.push(t.add_tsolid(shell_refs));
                }
                if !solid_refs.is_empty() {
                    t.add_tcompsolid(solid_refs);
                }
            }
        } else if !self.solids.is_empty() {
            for solid_shells in &self.solids {
                let shell_refs: Vec<ShapeRef> = solid_shells.iter().map(|&si| {
                    let shell_faces = self.shells.get(si).map_or(vec![], |sf| {
                        sf.iter().map(|&fi| face_refs.get(fi).copied().unwrap_or(ShapeRef::new(fi))).collect()
                    });
                    t.add_tshell(shell_faces)
                }).collect();
                t.add_tsolid(shell_refs);
            }
        } else {
            let shell_refs: Vec<ShapeRef> = self.shells.iter().map(|shell_faces| {
                let sfr: Vec<ShapeRef> = shell_faces.iter().map(|&fi| face_refs.get(fi).copied().unwrap_or(ShapeRef::new(fi))).collect();
                t.add_tshell(sfr)
            }).collect();
            t.add_tsolid(shell_refs);
        }

        let history = BooleanHistory {
            face_origins: self.face_origins,
            co_face_origins: self.co_face_origins,
            edge_origins: Vec::new(),
            vertex_origins: Vec::new(),
            shell_origins: Vec::new(),
            solid_origins: Vec::new(),
            tracker: HistoryTracker::new(),
            deleted_from_a: Vec::new(),
            deleted_from_b: Vec::new(),
            deletion_reasons: std::collections::HashMap::new(),
        };

        (t, history)
    }

    pub(crate) fn build(mut self) -> (BRep, BooleanHistory) {
        eprintln!("ResultBuilder::build: {} vertices, {} edges, {} faces", self.vertices.len(), self.edges.len(), self.faces.len());
        for (vi, p) in self.vertices.iter().enumerate() {
            eprintln!("  V[{}] = ({:.12}, {:.12}, {:.12})", vi, p.x, p.y, p.z);
        }
        // ✅ OCCT-aligned: pure conversion (BuildResult, Builder_1.cxx L130-168).
        // OCCT does NO vertex/edge merge, NO orphan edge removal, NO face culling.
        let vertices = self
            .vertices
            .into_iter()
            .map(|point| Vertex { point })
            .collect();

        let mut edges: Vec<Edge> = self
            .edges
            .into_iter()
            .map(|(start, end)| Edge { start, end })
            .collect();

        let mut geom = rcad_kernel::GeomStore::default();
        let mut faces = Vec::new();

        for (edge_indices, inner_wire_edges, triangles, normal, surface, uv_domain, _centroid, _area, sample_point, internal_wire_edges) in self.faces {
            let wire = Wire {
                edges: edge_indices.iter().map(|&(idx, forward)| {
                    if forward { WireEdge::fwd(idx) } else { WireEdge::rev(idx) }
                }).collect(),
            };
            let inner_wires: Vec<Wire> = inner_wire_edges
                .into_iter()
                .map(|wire_edge_idxs| Wire {
                    edges: wire_edge_idxs.iter().map(|&(idx, forward)| {
                        if forward { WireEdge::fwd(idx) } else { WireEdge::rev(idx) }
                    }).collect(),
                })
                .collect();
            // OCCT-aligned: Add internal wire edges to inner_wires for edge ref counting
            let mut inner_wires = inner_wires;
            for iw_edges in internal_wire_edges {
                let iw: Vec<WireEdge> = iw_edges.iter().map(|&(idx, forward)| {
                    if forward { WireEdge::fwd(idx) } else { WireEdge::rev(idx) }
                }).collect();
                if iw.len() >= 2 {
                    inner_wires.push(Wire { edges: iw });
                }
            }
            // BooleanBuilder faces inherit or accumulate triangles during
            // splitting, but those tessellations are not guaranteed to remain
            // valid for exact property evaluation after trimming/rewiring.
            // Keep them as fallback display meshes only; exact consumers should
            // regenerate or use analytic surface integration.
            let mesh_dirty = true;
            let surf_idx = geom.surfaces.len();
            faces.push(Face {
                outer_wire: wire,
                inner_wires,
                normal,
                triangles,
                sample_point: Some(sample_point),
                mesh_dirty,
                surface_idx: Some(surf_idx),
            });            geom.surfaces.push(surface);
            geom.face_surface.push(Some(surf_idx));
            geom.face_surface_range.push(uv_domain);
        }
        geom.face_internal_vertices = self.face_internal_vtx;

        // ✅ OCCT-aligned: set section edge curves from custom_edge_curves.
        //    OCCT BuildResult (Builder_1.cxx L130-168) does NOT:
        //      - remove orphan edges (every edge created by MakeSplitEdges is valid)
        //      - cull faces with <3 outer edges (BuilderFace produces valid wires)
        //      - check face outer-wire edge count
        //    OCCT simply iterates argument shapes of matching type and adds their
        //    images to myShape.  No post-processing needed because the DS ensures
        //    correct topology.
        //    rcad: custom_edge_curves store Circle3/BSpline curves for section edges.
        //    OCCT: MakeEdge(aIC, ...) creates BRep edge with exact analytic curve.
        //    rcad defaults to recompute_plane_surfaces (Line3), override here.
        if !edges.is_empty() {
            for (ei, curve_opt) in self.custom_edge_curves.iter().enumerate() {
                if ei >= edges.len() { break; }
                if let Some(crv) = curve_opt {
                    while geom.edge_curve.len() <= ei {
                        geom.edge_curve.push(None);
                    }
                    let curve_idx = geom.curves.len();
                    geom.curves.push(crv.clone());
                    geom.edge_curve[ei] = Some(curve_idx);
                }
            }
        }

        let history = BooleanHistory {
            face_origins: self.face_origins,
            co_face_origins: self.co_face_origins,
            edge_origins: Vec::new(),
            vertex_origins: Vec::new(),
            shell_origins: Vec::new(),
            solid_origins: Vec::new(),
            tracker: HistoryTracker::new(),
            deleted_from_a: Vec::new(),
            deleted_from_b: Vec::new(),
            deletion_reasons: std::collections::HashMap::new(),
        };

        // OCCT-aligned: set edge_degenerated flag for degenerated seam edges
        for &ei in &self.deg_edge_indices {
            while geom.edge_degenerated.len() <= ei {
                geom.edge_degenerated.push(false);
            }
            geom.edge_degenerated[ei] = true;
        }
        // ✅ OCCT-aligned: build shell/solid structure (Phase 5+6 or fallback).
        let brep_solids = if self.solids.is_empty() && self.shells.is_empty() {
            // Legacy path: single shell, single solid.
            vec![Solid { shells: vec![Shell { faces }] }]
        } else if !self.solids.is_empty() {
            // Phase 5+6: explicit shell/solid groups.
            let faceref = &faces;
            self.solids.iter().map(|solid_shells| Solid {
                shells: solid_shells.iter().map(|&si| Shell {
                    faces: self.shells.get(si).map_or(vec![], |shell_faces| {
                        shell_faces.iter().map(|&fi| faceref[fi].clone()).collect()
                    }),
                }).collect(),
            }).collect()
        } else {
            // Phase 5 only: shells exist but not grouped into solids.
            let faceref = &faces;
            vec![Solid {
                shells: self.shells.iter().map(|shell_faces| Shell {
                    faces: shell_faces.iter().map(|&fi| faceref[fi].clone()).collect(),
                }).collect(),
            }]
        };
        // ✅ OCCT-aligned: wrap result solids in CompSolid when source had one.
        let (brep_solids_out, compsolid) = if !self.compsolid_groups.is_empty() && !brep_solids.is_empty() {
            let cs = CompSolid { solids: brep_solids, label: None };
            (vec![], Some(cs))
        } else {
            (brep_solids, None)
        };
        let brep = BRep {
            vertices,
            edges,
            solids: brep_solids_out,
            geom,
            compound: None,
            compsolid,
        };
        eprintln!("BRep built: {} faces", brep.solids[0].shells[0].faces.len());
        (brep, history)
    }
}
