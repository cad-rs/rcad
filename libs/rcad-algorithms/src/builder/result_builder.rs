use std::collections::{HashMap, HashSet};
use glam::DVec2; use glam::DVec3;
use rcad_kernel::geom::*; use rcad_kernel::BRep;
use rcad_kernel::topods;
use crate::history::{BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin};
use crate::bopds::ds::*; use crate::tolerance::*;
use crate::builder::types::{WireFace, WireSegment, WireEdgeSource, WireSegmentTopoDS, WireEdgeSourceTopoDS, FaceEntry, FaceSampleData};
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
    pub(crate) tmp_shells: Vec<Vec<usize>>,
    pub(crate) tmp_solids: Vec<Vec<usize>>,
    pub(crate) custom_edge_curves: Vec<Option<Curve3>>,
    pub(crate) face_internal_vtx: Vec<Vec<usize>>,
    pub(crate) deg_edge_indices: HashSet<usize>,
    pub(crate) ic_edge_map: HashMap<usize, usize>,
    pub(crate) source_has_compound: bool,
    pub(crate) tmp_compsolid_groups: Vec<Vec<usize>>,
    /// OCCT-aligned: per-source side tracking for solids (0=ShapeA/Args, 1=ShapeB/Tools).
    ///   Parallel to tmp_solids: solid_side_origin[i] = side that produced tmp_solids[i].
    pub(crate) solid_side_origin: Vec<usize>,
    /// OCCT-aligned: compound groups (FillImagesCompound output).
    ///   Each entry is a Vec of solid indices (into result.solids) forming one compound.
    pub(crate) compound_groups: Vec<Vec<usize>>,
    /// OCCT-aligned: natural_restriction for each face in self.faces.
    ///   Parallel to face_origins.  Populated by emit_wire_face / build_original_face.
    pub(crate) face_natural_restriction: Vec<bool>,
    /// Final topods shell references (populated by build_topods from tmp_shells).
    pub(crate) shells: Vec<topods::ShapeRef>,
    /// Final topods solid references (populated by build_topods from tmp_solids).
    pub(crate) solids: Vec<topods::ShapeRef>,
    /// Final topods compsolid references (populated by build_topods from tmp_compsolid_groups).
    pub(crate) compsolid_groups: Vec<topods::ShapeRef>,
    /// Topods face ShapeRefs (populated by build_topods_faces after fill_images_faces).
    pub(crate) face_refs: Vec<topods::ShapeRef>,
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
            let (ei, forward) = if seg.is_closed_on_face {
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
                    WireEdgeSource::DsEdge(_) | WireEdgeSource::SeamEdge if seg.is_closed_on_face => {
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
        self.face_natural_restriction.push(ds.faces[face_idx].natural_restriction);
    }

    /// ✅ OCCT-aligned: emit_wire_face using WireSegmentTopoDS.
    ///   Same logic as emit_wire_face but reads edge/vertex data through BRepTool.
    pub(crate) fn emit_wire_face_topods(
        &mut self,
        face_idx: usize,
        wf: &WireFace,
        segments: &[WireSegmentTopoDS],
        tool: &dyn rcad_kernel::topods::BRepTool,
        ic_curves: &HashMap<usize, Curve3>,
        flip: bool,
        origin: FaceOrigin,
        vertex_positions: &HashMap<usize, DVec3>,
        face_ref: rcad_kernel::topods::ShapeRef,
        natural_restriction: bool,
    ) {
        let normal = if let Some(surf) = tool.face_surface(face_ref) {
            match surf {
                Surface3::Plane(p) => if flip { -p.normal } else { p.normal },
                _ => {
                    let n = Self::estimate_boundary_normal_from_segments_topo(
                        &wf.outer_wire, segments, tool, vertex_positions);
                    if n.length_squared() > TOLERANCE_METRIC_SQ_NEAR_ZERO { n } else { return; }
                }
            }
        } else { return; };

        let get_pos = |vi: usize| -> DVec3 {
            vertex_positions.get(&vi).copied()
                .unwrap_or_else(|| tool.vertex_position(rcad_kernel::topods::ShapeRef::new(vi)))
        };

        let mut vert_indices = Vec::new();
        let mut edge_indices = Vec::new();
        let ow: Vec<&usize> = wf.outer_wire.iter()
            .filter(|&&si| segments[si].start_vertex.index != segments[si].end_vertex.index)
            .collect();
        for &&si in &ow {
            let seg = &segments[si];
            let v1 = if vertex_positions.contains_key(&seg.start_vertex.index) {
                self.add_vertex(vertex_positions[&seg.start_vertex.index])
            } else {
                self.add_ds_vertex(seg.start_vertex.index, get_pos(seg.start_vertex.index))
            };
            let v2 = if vertex_positions.contains_key(&seg.end_vertex.index) {
                self.add_vertex(vertex_positions[&seg.end_vertex.index])
            } else {
                self.add_ds_vertex(seg.end_vertex.index, get_pos(seg.end_vertex.index))
            };
            if vert_indices.is_empty() || vert_indices.last() != Some(&v1) {
                vert_indices.push(v1);
            }
            let (ei, forward) = if seg.is_closed_on_face {
                let seam_deg = (get_pos(seg.start_vertex.index)
                    - get_pos(seg.end_vertex.index)).length_squared() < TOLERANCE_ABS_SQ;
                let sphere_surf = match tool.face_surface(face_ref) {
                    Some(Surface3::Sphere(s)) => s.clone(),
                    _ => SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, ref_dir: DVec3::X },
                };
                let is_canon_deg = vertex_positions.contains_key(&seg.start_vertex.index)
                    || vertex_positions.contains_key(&seg.end_vertex.index);
                let ei = if seam_deg || is_canon_deg {
                    self.add_edge_seam_degenerate(v1, v2, &sphere_surf)
                } else {
                    let seam_normal = any_perpendicular(sphere_surf.axis).normalize();
                    let seam_circle = Curve3::Circle(Circle3 {
                        center: sphere_surf.center, normal: seam_normal, radius: sphere_surf.radius,
                    });
                    self.add_seam_edge(v1, v2, seam_circle)
                };
                (ei, true)
            } else {
                let ei = match &seg.source {
                    super::types::WireEdgeSourceTopoDS::IntersectionCurve(ci) => {
                        let crv = ic_curves.get(&ci.index)
                            .cloned().unwrap_or(Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X }));
                        self.add_ic_edge(ci.index, v1, v2, crv)
                    }
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
                let v1 = if vertex_positions.contains_key(&seg.start_vertex.index) {
                    self.add_vertex(vertex_positions[&seg.start_vertex.index])
                } else {
                    self.add_ds_vertex(seg.start_vertex.index, get_pos(seg.start_vertex.index))
                };
                let v2 = if vertex_positions.contains_key(&seg.end_vertex.index) {
                    self.add_vertex(vertex_positions[&seg.end_vertex.index])
                } else {
                    self.add_ds_vertex(seg.end_vertex.index, get_pos(seg.end_vertex.index))
                };
                if iw_verts.is_empty() || iw_verts.last() != Some(&v1) {
                    iw_verts.push(v1);
                }
                let ei = match &seg.source {
                    super::types::WireEdgeSourceTopoDS::IntersectionCurve(ci) => {
                        let crv = ic_curves.get(&ci.index)
                            .cloned().unwrap_or(Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X }));
                        self.add_ic_edge(ci.index, v1, v2, crv)
                    }
                    _ => self.add_edge(v1, v2),
                };
                let forward = self.edges[ei].0 == v1;
                iw_edges.push((ei, forward));
            }
            inner_wire_edges.push(iw_edges);
            iw_vert_indices_all.extend(iw_verts);
        }

        // Internal wire edges
        let mut internal_wire_edges: Vec<Vec<(usize, bool)>> = Vec::new();
        for iw in &wf.internal_wires {
            let mut iw_edges = Vec::new();
            for &si in iw {
                let seg = &segments[si];
                let v1 = if vertex_positions.contains_key(&seg.start_vertex.index) {
                    self.add_vertex(vertex_positions[&seg.start_vertex.index])
                } else {
                    self.add_ds_vertex(seg.start_vertex.index, get_pos(seg.start_vertex.index))
                };
                let v2 = if vertex_positions.contains_key(&seg.end_vertex.index) {
                    self.add_vertex(vertex_positions[&seg.end_vertex.index])
                } else {
                    self.add_ds_vertex(seg.end_vertex.index, get_pos(seg.end_vertex.index))
                };
                let ei = match &seg.source {
                    super::types::WireEdgeSourceTopoDS::IntersectionCurve(ci) => {
                        let crv = ic_curves.get(&ci.index)
                            .cloned().unwrap_or(Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X }));
                        if let Curve3::Circle(_) = &crv { self.add_circle_edge(v1, v2, crv) }
                        else { self.add_edge(v1, v2) }
                    }
                    _ if seg.is_closed_on_face => {
                        let sphere_surf = match tool.face_surface(face_ref) {
                            Some(Surface3::Sphere(s)) => s.clone(),
                            _ => SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, ref_dir: DVec3::X },
                        };
                        let c = Curve3::Circle(Circle3 {
                            center: sphere_surf.center,
                            normal: any_perpendicular(sphere_surf.axis).normalize(),
                            radius: sphere_surf.radius,
                        });
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
                if pts.is_empty() || pts.last() != Some(&a) { pts.push(a); }
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

        let centroid = outer_boundary.iter().copied().sum::<DVec3>() / outer_boundary.len().max(1) as f64;
        let area = Self::polygon_signed_area_on_normal(&outer_boundary, normal);
        let mut outer_sig: Vec<usize> = edge_indices.iter().map(|&(eid, _)| eid).collect();
        outer_sig.sort_unstable();
        let nlen = normal.length();
        let nunit = if nlen > TOLERANCE_LEN_MIN { normal / nlen } else { normal };
        // Same coincident face dedup as emit_wire_face
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

        // UV domain for sphere faces
        let sphere_uv = if let Some(Surface3::Sphere(sph)) = tool.face_surface(face_ref) {
            let uvs: Vec<DVec2> = if !wf.outer_wire.is_empty() {
                wf.outer_wire.iter().map(|&si| {
                    let pos = get_pos(segments[si].start_vertex.index);
                    sph.world_to_uv(pos)
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

        let surface = tool.face_surface(face_ref)
            .cloned().unwrap_or(Surface3::Plane(rcad_kernel::geom::Plane { origin: glam::DVec3::ZERO, normal: glam::DVec3::Z }));

        let sample_pt = if !wf.outer_wire.is_empty() {
            get_pos(segments[wf.outer_wire[0]].start_vertex.index)
        } else { DVec3::ZERO };

        self.face_internal_vtx.push(Vec::new());
        self.faces.push((
            edge_indices, inner_wire_edges, tris, normal, surface,
            sphere_uv, centroid, area, sample_pt, internal_wire_edges,
        ));
        self.face_origins.push(origin);
        self.face_natural_restriction.push(natural_restriction);
    }

    /// ✅ OCCT-aligned: estimate face normal from wire segments (TopoDS variant).
    fn estimate_boundary_normal_from_segments_topo(
        outer_wire: &[usize],
        segments: &[super::types::WireSegmentTopoDS],
        tool: &dyn rcad_kernel::topods::BRepTool,
        vertex_positions: &HashMap<usize, DVec3>,
    ) -> DVec3 {
        if outer_wire.len() < 3 { return DVec3::ZERO; }
        let pts: Vec<DVec3> = outer_wire.iter().map(|&si| {
            let seg = &segments[si];
            let vi = seg.start_vertex.index;
            vertex_positions.get(&vi).copied()
                .unwrap_or_else(|| tool.vertex_position(rcad_kernel::topods::ShapeRef::new(vi)))
        }).collect();
        Self::estimate_boundary_normal(&pts)
    }

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
            tmp_shells: Vec::new(),
            tmp_solids: Vec::new(),
            source_has_compound: false,
            tmp_compsolid_groups: Vec::new(),
            solid_side_origin: Vec::new(),
            compound_groups: Vec::new(),
            face_natural_restriction: Vec::new(),
            shells: Vec::new(),
            solids: Vec::new(),
            compsolid_groups: Vec::new(),
            face_refs: Vec::new(),
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
        self.face_natural_restriction.push(ds.faces[fi].natural_restriction);
    }

    /// ✅ OCCT-aligned: BuildResult(COMPSOLID) — build compsolids via BRepBuilder.
    ///   OCCT: BOPAlgo_Builder::BuildResult (Builder_1.cxx L130-168) iterates
    ///   source COMPSOLID shapes and adds their split images to myShape via
    ///   BRep_Builder::Add.  rcad: processes tmp_compsolid_groups (groups of
    ///   solid indices from fill_images_containers_compsolid) and creates
    ///   topods compsolids using BRepBuilder::make_compsolid.
    pub(crate) fn build_compsolids(&mut self, t: &mut topods::BRep, groups: Vec<Vec<usize>>) {
        let mut bb = topods::BRepBuilder::new();
        for cs_group in &groups {
            let solid_refs: Vec<topods::ShapeRef> = cs_group.iter()
                .filter_map(|&si| self.solids.get(si).copied())
                .collect();
            if !solid_refs.is_empty() {
                self.compsolid_groups.push(bb.make_compsolid(t, solid_refs));
            }
        }
    }

    /// ✅ OCCT-aligned: BuildResult(COMPOUND) — build compounds via BRepBuilder.
    ///   OCCT: BOPAlgo_Builder::BuildResult (Builder_1.cxx L130-168) iterates
    ///   source COMPOUND shapes and adds their split images to myShape via
    ///   BRep_Builder::Add.  rcad: processes compound_groups (groups of solid
    ///   indices from fill_images_compounds) and creates topods compounds
    ///   using BRepBuilder::make_compound.
    pub(crate) fn build_compounds(&mut self, t: &mut topods::BRep, groups: &[Vec<usize>]) {
        let mut bb = topods::BRepBuilder::new();
        for group in groups {
            let solid_refs: Vec<topods::ShapeRef> = group.iter()
                .filter_map(|&si| self.solids.get(si).copied())
                .collect();
            if !solid_refs.is_empty() {
                bb.make_compound(t, solid_refs);
            }
        }
    }

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



    /// ✅ OCCT-aligned: BuildResult(FACE) — create topods vertices, edges, wires, faces
    ///   in t_brep from the flat arrays.  Called after fill_images_faces so that
    ///   later BuildResult(SHELL) / BuildResult(SOLID) can reference these ShapeRefs.
    ///   OCCT: BuildResult(FACE) (Builder_1.cxx L130-168) adds face images to myShape.
    ///   rcad: previously deferred to build_topods; now aligned to happen per-phase.
    pub(crate) fn build_topods_faces(&mut self, t: &mut topods::BRep) {
        use topods::{Orientation, ShapeRef};

        // 1. Vertices → TShape::Vertex
        for v in &self.vertices {
            t.add_tvertex(*v);
        }

        // 2. Edges → TShape::Edge
        let mut e_map: Vec<ShapeRef> = Vec::with_capacity(self.edges.len());
        for (ei, &(start, end)) in self.edges.iter().enumerate() {
            let first = ShapeRef::new(start);
            let last = ShapeRef::new(end);
            let curve_idx = self.custom_edge_curves.get(ei).and_then(|c| c.as_ref()).map(|crv| {
                let ci = t.curves.len();
                t.curves.push(crv.clone());
                ci
            });
            e_map.push(t.add_tedge(curve_idx, first, last, [0.0, 1.0]));
        }

        // 3. Faces → TShape::Face (with wires)
        self.face_refs.clear();
        let mut flat_fi = 0usize;
        for (edge_indices, inner_wire_edges, _triangles, _normal, surface, _uv_domain, _centroid, _area, sample_point, internal_wire_edges) in &self.faces {
            // Outer wire
            let outer_edges: Vec<ShapeRef> = edge_indices.iter().map(|&(idx, forward)| {
                let orient = if forward { Orientation::Forward } else { Orientation::Reversed };
                if idx < e_map.len() { ShapeRef::with_orientation(e_map[idx].index, orient) }
                else { ShapeRef::with_orientation(idx, orient) }
            }).collect();
            let outer_wire = t.add_twire(outer_edges);

            // Inner wires
            let mut inner_wires = Vec::new();
            for wire_idxs in inner_wire_edges {
                let iw_edges: Vec<ShapeRef> = wire_idxs.iter().map(|&(idx, forward)| {
                    let orient = if forward { Orientation::Forward } else { Orientation::Reversed };
                    if idx < e_map.len() { ShapeRef::with_orientation(e_map[idx].index, orient) }
                    else { ShapeRef::with_orientation(idx, orient) }
                }).collect();
                if !iw_edges.is_empty() {
                    let w = t.add_twire(iw_edges);
                    t.wire_mut(w).closed = true;
                    inner_wires.push(w);
                }
            }
            // Internal wire edges
            for iw_edges in internal_wire_edges {
                let iw: Vec<ShapeRef> = iw_edges.iter().map(|&(idx, forward)| {
                    let orient = if forward { Orientation::Forward } else { Orientation::Reversed };
                    if idx < e_map.len() { ShapeRef::with_orientation(e_map[idx].index, orient) }
                    else { ShapeRef::with_orientation(idx, orient) }
                }).collect();
                if iw.len() >= 2 {
                    let w = t.add_twire(iw);
                    t.wire_mut(w).closed = true;
                    inner_wires.push(w);
                }
            }

            let surf_idx = t.surfaces.len();
            t.surfaces.push(surface.clone());
            let internal_vtx: Vec<ShapeRef> = self.face_internal_vtx.get(flat_fi)
                .map_or(vec![], |v| v.iter().map(|&vi| ShapeRef::new(vi)).collect());
            let nr = self.face_natural_restriction.get(flat_fi).copied().unwrap_or(true);
            self.face_refs.push(t.add_tface(Some(surf_idx), outer_wire, inner_wires, Some(*sample_point), *_uv_domain, internal_vtx, nr));
            flat_fi += 1;
        }
    }

    /// ✅ OCCT-aligned: Final assembly — return history (PIOperation_FillHistory).
    ///   Per-dimension BuildResult calls (Face/Shell/Solid/CompSolid) have already
    ///   created the corresponding topods TShapes in t_brep. This method returns
    ///   the BooleanHistory from the accumulated result data.
    ///   When `fill_history` is false (OCCT: !HasHistory → !myFillHistory),
    ///   returns an empty history with no origins tracking.
    pub(crate) fn build_topods(&mut self, t: &mut topods::BRep, fill_history: bool) -> BooleanHistory {
        // Fallback: if no solids were created by BuildResult but faces exist,
        // create a default shell + solid (OCCT: defaults to single shell/solid).
        if self.shells.is_empty() && !self.face_refs.is_empty() {
            let shell = t.add_tshell(std::mem::take(&mut self.face_refs));
            t.add_tsolid(vec![shell]);
        }

        BooleanHistory {
            face_origins: if fill_history { std::mem::take(&mut self.face_origins) } else { vec![] },
            co_face_origins: if fill_history { std::mem::take(&mut self.co_face_origins) } else { vec![] },
            edge_origins: Vec::new(),
            vertex_origins: Vec::new(),
            shell_origins: Vec::new(),
            solid_origins: Vec::new(),
            tracker: HistoryTracker::new(),
            deleted_from_a: Vec::new(),
            deleted_from_b: Vec::new(),
            deletion_reasons: std::collections::HashMap::new(),
            source_history: Vec::new(),
        }
    }

    pub(crate) fn build(mut self) -> (BRep, BooleanHistory) {
        eprintln!("ResultBuilder::build: {} vertices, {} edges, {} faces", self.vertices.len(), self.edges.len(), self.faces.len());
        for (vi, p) in self.vertices.iter().enumerate() {
            eprintln!("  V[{}] = ({:.12}, {:.12}, {:.12})", vi, p.x, p.y, p.z);
        }
        // ✅ OCCT-aligned: pure conversion (BuildResult, Builder_1.cxx L130-168).
        // OCCT does NO vertex/edge merge, NO orphan edge removal, NO face culling.
        let vertices: Vec<rcad_kernel::topology::Vertex> = self
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
        let brep = BRep::new(); // dead code path — real builds use build_topods

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
            source_history: Vec::new(),
        };
        eprintln!("BRep built: {} faces", brep.solids[0].shells[0].faces.len());
        (brep, history)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create a simple ResultBuilder with 4 vertices, 4 edges forming a square,
    /// 1 planar face, and optionally populated tmp_shells / tmp_solids.
    fn make_test_builder(with_shells: bool, with_solids: bool) -> ResultBuilder {
        let mut rb = ResultBuilder::new();
        // 4 vertices: square at z=0
        let v0 = rb.add_vertex(DVec3::new(0.0, 0.0, 0.0));
        let v1 = rb.add_vertex(DVec3::new(1.0, 0.0, 0.0));
        let v2 = rb.add_vertex(DVec3::new(1.0, 1.0, 0.0));
        let v3 = rb.add_vertex(DVec3::new(0.0, 1.0, 0.0));
        // 4 edges: square perimeter (outer wire)
        let e0 = rb.add_edge_occt(v0, v1);
        let e1 = rb.add_edge_occt(v1, v2);
        let e2 = rb.add_edge_occt(v2, v3);
        let e3 = rb.add_edge_occt(v3, v0);
        // 1 planar face
        let outer = vec![(e0, true), (e1, true), (e2, true), (e3, true)];
        let centroid = DVec3::new(0.5, 0.5, 0.0);
        rb.faces.push((
            outer,
            vec![],                    // inner wires
            vec![],
            DVec3::Z,
            Surface3::Plane(rcad_kernel::geom::Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
            }),
            None,
            centroid,
            1.0,
            centroid,
            vec![],
        ));
        rb.face_origins.push(FaceOrigin::FromA(0));
        if with_shells {
            rb.tmp_shells.push(vec![0]); // shell 0 contains face 0
        }
        if with_solids {
            rb.tmp_solids.push(vec![0]); // solid 0 contains shell 0
        }
        rb
    }

    #[test]
    fn build_topods_legacy_fallback() {
        // empty tmp_shells + empty tmp_solids → legacy path:
        //   one shell from all faces, one solid wrapping that shell.
        let mut rb = make_test_builder(false, false);
        let mut t = topods::BRep::new();
        rb.build_topods_faces(&mut t);
        let _history = rb.build_topods(&mut t, true);

        // tshapes: 4 vertices + 4 edges + 2 wires + 1 face + 1 shell + 1 solid
        assert!(t.tshapes.len() >= 12, "expected >= 12 tshapes, got {}", t.tshapes.len());

        // Verify Solid exists and references the shell
        let solid_count = t.tshapes.iter().filter(|ts| matches!(&***ts, topods::TShape::Solid(_))).count();
        assert_eq!(solid_count, 1, "legacy path should produce exactly 1 solid");

        // The result BRep (via from_topods) should have 1 solid
        let brep = rcad_kernel::BRep::from_topods(&t);
        assert_eq!(brep.solids.len(), 1);
        assert_eq!(brep.solids[0].shells.len(), 1);
        assert_eq!(brep.solids[0].shells[0].faces.len(), 1);
    }

    #[test]
    fn round_trip_single_solid_preserves_topology() {
        let orig = rcad_modeling::make_box_brep(
            DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0,
        ).unwrap();

        let t_in = orig.to_topods();
        let (mut rb, _) = builder_from_topods(&t_in);
        let mut t = topods::BRep::new();
        rb.build_topods_faces(&mut t);
        let _history = rb.build_topods(&mut t, true);
        let rebuilt = rcad_kernel::BRep::from_topods(&t);

        assert_eq!(rebuilt.solids.len(), orig.solids.len(),
            "solid count mismatch");
        let ns = |b: &rcad_kernel::BRep| -> usize {
            b.solids.iter().flat_map(|s| &s.shells).count()
        };
        let nf = |b: &rcad_kernel::BRep| -> usize {
            b.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count()
        };
        assert_eq!(ns(&rebuilt), ns(&orig), "shell count mismatch");
        assert_eq!(nf(&rebuilt), nf(&orig), "face count mismatch");
    }

    /// Round-trip with two solids sharing a face (adjacent boxes).
    /// Verifies that shared-edge identity is preserved across conversion.
    #[test]
    /// Square face with a triangular hole (inner wire).
    /// The inner wire must survive the round-trip.
    #[test]
    fn round_trip_face_with_inner_wire() {
        let mut rb = ResultBuilder::new();
        // Outer square: 4 vertices
        let o0 = rb.add_vertex(DVec3::new(0.0, 0.0, 0.0));
        let o1 = rb.add_vertex(DVec3::new(3.0, 0.0, 0.0));
        let o2 = rb.add_vertex(DVec3::new(3.0, 3.0, 0.0));
        let o3 = rb.add_vertex(DVec3::new(0.0, 3.0, 0.0));
        let oe0 = rb.add_edge_occt(o0, o1);
        let oe1 = rb.add_edge_occt(o1, o2);
        let oe2 = rb.add_edge_occt(o2, o3);
        let oe3 = rb.add_edge_occt(o3, o0);
        // Inner triangle: 3 vertices
        let i0 = rb.add_vertex(DVec3::new(1.0, 1.0, 0.0));
        let i1 = rb.add_vertex(DVec3::new(2.0, 1.0, 0.0));
        let i2 = rb.add_vertex(DVec3::new(1.5, 2.0, 0.0));
        let ie0 = rb.add_edge_occt(i0, i1);
        let ie1 = rb.add_edge_occt(i1, i2);
        let ie2 = rb.add_edge_occt(i2, i0);

        let outer = vec![(oe0, true), (oe1, true), (oe2, true), (oe3, true)];
        let inner = vec![vec![(ie0, true), (ie1, true), (ie2, true)]];
        let c = DVec3::new(1.5, 1.5, 0.0);
        rb.faces.push((
            outer, inner, vec![], DVec3::Z,
            Surface3::Plane(rcad_kernel::geom::Plane { origin: DVec3::ZERO, normal: DVec3::Z }),
            None, c, 1.0, c, vec![],
        ));
        rb.face_origins.push(FaceOrigin::FromA(0));

        let mut t = topods::BRep::new();
        rb.build_topods_faces(&mut t);
        let _history = rb.build_topods(&mut t, true);
        let rebuilt = rcad_kernel::BRep::from_topods(&t);

        assert_eq!(rebuilt.solids.len(), 1);
        let face = &rebuilt.solids[0].shells[0].faces[0];
        assert_eq!(face.inner_wires.len(), 1,
            "should have 1 inner wire (triangle hole)");
        assert_eq!(face.inner_wires[0].edges.len(), 3,
            "inner triangle has 3 edges");
        assert_eq!(face.outer_wire.edges.len(), 4,
            "outer square has 4 edges");
    }

    /// Direct topods → builder → topods round-trip without going through BRep.
    /// Verifies that build_topods correctly reconstructs shells and solids
    /// given only tmp_shells/tmp_solids populated from any source.
    #[test]
    fn direct_topods_round_trip_preserves_tshape_count() {
        let box_brep = rcad_modeling::make_box_brep(
            DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0,
        ).unwrap();
        let t0 = box_brep.to_topods();
        // Collect face/edge/vertex info from t0
        let (mut rb, _) = builder_from_topods(&t0);

        let mut t1 = topods::BRep::new();
        rb.build_topods_faces(&mut t1);
        let _history = rb.build_topods(&mut t1, true);

        // Same number of Vertex/Edge/Face/Shell/Solid TShapes
        let count_by_type = |t: &topods::BRep| -> (usize, usize, usize, usize, usize) {
            let (mut v, mut e, mut f, mut sh, mut so) = (0, 0, 0, 0, 0);
            for ts in &t.tshapes {
                match &**ts {
                    topods::TShape::Vertex(_) => v += 1,
                    topods::TShape::Edge(_) => e += 1,
                    topods::TShape::Face(_) => f += 1,
                    topods::TShape::Shell(_) => sh += 1,
                    topods::TShape::Solid(_) => so += 1,
                    _ => {}
                }
            }
            (v, e, f, sh, so)
        };
        assert_eq!(count_by_type(&t0), count_by_type(&t1),
            "TShape type counts must match after round-trip");
    }

    /// Internal wire edges (TopAbs_INTERNAL) survive conversion.
    #[test]
    fn round_trip_internal_wire_edges() {
        let mut rb = ResultBuilder::new();
        let v0 = rb.add_vertex(DVec3::new(0.0, 0.0, 0.0));
        let v1 = rb.add_vertex(DVec3::new(1.0, 0.0, 0.0));
        let v2 = rb.add_vertex(DVec3::new(1.0, 1.0, 0.0));
        let v3 = rb.add_vertex(DVec3::new(0.0, 1.0, 0.0));
        let e0 = rb.add_edge_occt(v0, v1);
        let e1 = rb.add_edge_occt(v1, v2);
        let e2 = rb.add_edge_occt(v2, v3);
        let e3 = rb.add_edge_occt(v3, v0);
        // Internal edge: a seam line inside the face (needs >= 2 edges for a valid wire)
        let vi0 = rb.add_vertex(DVec3::new(0.2, 0.2, 0.0));
        let vi1 = rb.add_vertex(DVec3::new(0.5, 0.8, 0.0));
        let vi2 = rb.add_vertex(DVec3::new(0.8, 0.2, 0.0));
        let ei0 = rb.add_edge(vi0, vi1);
        let ei1 = rb.add_edge(vi1, vi2);
        let internal_edges = vec![vec![(ei0, true), (ei1, true)]];

        let outer = vec![(e0, true), (e1, true), (e2, true), (e3, true)];
        let c = DVec3::new(0.5, 0.5, 0.0);
        rb.faces.push((
            outer, vec![], vec![], DVec3::Z,
            Surface3::Plane(rcad_kernel::geom::Plane { origin: DVec3::ZERO, normal: DVec3::Z }),
            None, c, 1.0, c, internal_edges,
        ));
        rb.face_origins.push(FaceOrigin::FromA(0));

        let mut t = topods::BRep::new();
        rb.build_topods_faces(&mut t);
        let _history = rb.build_topods(&mut t, true);
        let rebuilt = rcad_kernel::BRep::from_topods(&t);

        assert_eq!(rebuilt.solids.len(), 1);
        let face = &rebuilt.solids[0].shells[0].faces[0];
        // Internal wire edges become inner wires in the rebuilt BRep
        let total_inner_edges: usize = face.inner_wires.iter()
            .map(|w| w.edges.len()).sum();
        assert_eq!(total_inner_edges, 2,
            "internal seam edges should survive round-trip");
    }

    /// merge two BReps into one topods::BRep by concatenating their vertices/edges/faces
    /// after re-indexing, creating separate solids.
    fn merge_two_breps_into_topods(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> topods::BRep {
        use topods::{Orientation, ShapeRef};
        let mut t = topods::BRep::new();
        t.curves = a.geom.curves.clone();
        t.curves.extend(b.geom.curves.clone());
        t.surfaces = a.geom.surfaces.clone();
        t.surfaces.extend(b.geom.surfaces.clone());

        // Helper to add all vertices/edges/faces/shells from one BRep into t,
        // returning (shell_sr) for that solid.
        fn push_brep_solid(t: &mut topods::BRep, brep: &rcad_kernel::BRep) -> ShapeRef {
            let mut v_map: Vec<ShapeRef> = Vec::new();
            for v in &brep.vertices {
                v_map.push(t.add_tvertex(v.point));
            }
            let mut e_map: Vec<ShapeRef> = Vec::new();
            for (ei, e) in brep.edges.iter().enumerate() {
                let first = v_map[e.start];
                let last = v_map[e.end];
                let curve = brep.geom.edge_curve.get(ei).copied().flatten();
                let range = brep.geom.edge_curve_range.get(ei).copied().flatten().unwrap_or([0.0, 0.0]);
                e_map.push(t.add_tedge(curve, first, last, range));
            }
            let mut face_refs = Vec::new();
            let mut fi = 0usize;
            for solid in &brep.solids {
                for shell in &solid.shells {
                    for face in &shell.faces {
                        let outer_edges: Vec<ShapeRef> = face.outer_wire.edges.iter().map(|we| {
                            let orient = if we.forward { Orientation::Forward } else { Orientation::Reversed };
                            ShapeRef::with_orientation(e_map[we.idx].index, orient)
                        }).collect();
            let outer_wire = t.add_twire(outer_edges);
            t.wire_mut(outer_wire).closed = true;
                        let inner_wires: Vec<ShapeRef> = face.inner_wires.iter().map(|w| {
                            let iwe: Vec<ShapeRef> = w.edges.iter().map(|we| {
                                let orient = if we.forward { Orientation::Forward } else { Orientation::Reversed };
                                ShapeRef::with_orientation(e_map[we.idx].index, orient)
                            }).collect();
                            t.add_twire(iwe)
                        }).collect();
                        let internal_vtx: Vec<ShapeRef> = brep.geom.face_internal_vertices
                            .get(fi)
                            .map(|v| v.iter().map(|&vi| ShapeRef::new(vi)).collect())
                            .unwrap_or_default();
                        face_refs.push(t.add_tface(
                            face.surface_idx, outer_wire, inner_wires,
                            face.sample_point, None, internal_vtx, true,
                        ));
                        fi += 1;
                    }
                }
            }
            let shell = t.add_tshell(face_refs);
            t.add_tsolid(vec![shell])
        }

        push_brep_solid(&mut t, a);
        push_brep_solid(&mut t, b);
        t
    }

    /// Extract a ResultBuilder from an existing topods::BRep (inverse of build_topods).
    /// Reads all TShape entries and populates ResultBuilder fields accordingly.
    fn builder_from_topods(t: &topods::BRep) -> (ResultBuilder, topods::BRep) {
        let mut rb = ResultBuilder::new();
        let mut v_map: Vec<usize> = Vec::new();
        let mut e_map: Vec<Option<usize>> = Vec::new();

        for (ti, ts) in t.tshapes.iter().enumerate() {
            match &**ts {
                topods::TShape::Vertex(vd) => {
                    let vi = rb.add_vertex(vd.point);
                    while v_map.len() <= ti { v_map.push(0); }
                    v_map[ti] = vi;
                }
                topods::TShape::Edge(ed) => {
                    let v1 = v_map[ed.first.index];
                    let v2 = v_map[ed.last.index];
                    let ei = rb.add_edge(v1, v2);
                    while e_map.len() <= ti { e_map.push(None); }
                    e_map[ti] = Some(ei);
                    if let Some(ci) = ed.curve {
                        let crv = t.curves[ci].clone();
                        while rb.custom_edge_curves.len() <= ei {
                            rb.custom_edge_curves.push(None);
                        }
                        rb.custom_edge_curves[ei] = Some(crv);
                    }
                }
                _ => {}
            }
        }

        let mut shell_to_rb: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        for (ti, ts) in t.tshapes.iter().enumerate() {
            if let topods::TShape::Solid(sd) = &**ts {
                for shell_sr in &sd.shells {
                    if let topods::TShape::Shell(shd) = &*t.tshapes[shell_sr.index] {
                        let mut shell_faces: Vec<usize> = Vec::new();
                        for face_sr in &shd.faces {
                            if let topods::TShape::Face(fd) = &*t.tshapes[face_sr.index] {
                                let outer: Vec<(usize, bool)> = t.wire(fd.outer_wire).edges.iter().map(|e_sr| {
                                    let ei = e_map.get(e_sr.index).copied().flatten().unwrap_or(e_sr.index);
                                    (ei, e_sr.orientation.is_forward())
                                }).collect();
                                let inner: Vec<Vec<(usize, bool)>> = fd.inner_wires.iter().map(|w_sr| {
                                    t.wire(*w_sr).edges.iter().map(|e_sr| {
                                        let ei = e_map.get(e_sr.index).copied().flatten().unwrap_or(e_sr.index);
                                        (ei, e_sr.orientation.is_forward())
                                    }).collect()
                                }).collect();
                                let surf = fd.surface.map(|si| t.surfaces[si].clone())
                                    .unwrap_or(Surface3::Plane(rcad_kernel::geom::Plane {
                                        origin: DVec3::ZERO, normal: DVec3::Z,
                                    }));
                                let sample = fd.sample_point.unwrap_or(DVec3::ZERO);
                                rb.faces.push((
                                    outer, inner, vec![], DVec3::Z, surf,
                                    fd.uv_domain, sample, 1.0, sample, vec![],
                                ));
                                rb.face_origins.push(FaceOrigin::FromA(0));
                                shell_faces.push(rb.faces.len() - 1);
                            }
                        }
                        if !shell_faces.is_empty() {
                            rb.tmp_shells.push(shell_faces);
                            shell_to_rb.insert(shell_sr.index, rb.tmp_shells.len() - 1);
                        }
                    }
                }
                let solid_shell_indices: Vec<usize> = sd.shells.iter()
                    .filter_map(|sr| shell_to_rb.get(&sr.index).copied())
                    .collect();
                if !solid_shell_indices.is_empty() {
                    rb.tmp_solids.push(solid_shell_indices);
                }
            }
        }

        (rb, t.clone())
    }

    #[test]
    fn test_face_natural_restriction_tracked_in_builder() {
        let mut rb = ResultBuilder::new();
        // Add a face directly (simulates emit_wire_face)
        let v0 = rb.add_vertex(DVec3::ZERO);
        let v1 = rb.add_vertex(DVec3::X);
        let v2 = rb.add_vertex(DVec3::new(0.5, 0.866, 0.0));
        let e0 = rb.add_edge(v0, v1);
        let e1 = rb.add_edge(v1, v2);
        let e2 = rb.add_edge(v2, v0);
        let outer = vec![(e0, true), (e1, true), (e2, true)];
        let centroid = DVec3::new(0.5, 0.289, 0.0);
        rb.faces.push((
            outer, vec![], vec![], DVec3::Z,
            Surface3::Plane(rcad_kernel::geom::Plane { origin: DVec3::ZERO, normal: DVec3::Z }),
            None, centroid, 1.0, centroid, vec![],
        ));
        rb.face_origins.push(FaceOrigin::FromA(0));
        rb.face_natural_restriction.push(false);

        assert_eq!(rb.face_natural_restriction.len(), 1);
        assert!(!rb.face_natural_restriction[0]);
    }

    #[test]
    fn test_face_natural_restriction_propagated_to_topods_via_build_topods_faces() {
        let mut rb = ResultBuilder::new();
        let v0 = rb.add_vertex(DVec3::ZERO);
        let v1 = rb.add_vertex(DVec3::X);
        let v2 = rb.add_vertex(DVec3::new(0.5, 0.866, 0.0));
        let e0 = rb.add_edge(v0, v1);
        let e1 = rb.add_edge(v1, v2);
        let e2 = rb.add_edge(v2, v0);
        let outer = vec![(e0, true), (e1, true), (e2, true)];
        let centroid = DVec3::new(0.5, 0.289, 0.0);
        rb.faces.push((
            outer, vec![], vec![], DVec3::Z,
            Surface3::Plane(rcad_kernel::geom::Plane { origin: DVec3::ZERO, normal: DVec3::Z }),
            None, centroid, 1.0, centroid, vec![],
        ));
        rb.face_origins.push(FaceOrigin::FromA(0));
        rb.face_natural_restriction.push(false);

        let mut t = topods::BRep::new();
        rb.build_topods_faces(&mut t);

        // Find the face TShape and check natural_restriction
        let face_count = t.tshapes.iter().filter(|ts| matches!(&***ts, topods::TShape::Face(_))).count();
        assert_eq!(face_count, 1);
        let face_sr = rb.face_refs[0];
        let fd = t.face(face_sr);
        assert!(!fd.natural_restriction);
    }

    #[test]
    fn test_build_compsolids_creates_compsolid() {
        let mut rb = ResultBuilder::new();
        let mut t = topods::BRep::new();
        // Create two solids in t
        let v0 = t.add_tvertex(DVec3::ZERO);
        let v1 = t.add_tvertex(DVec3::X);
        let e = t.add_tedge(None, v0, v1, [0.0, 1.0]);
        let w = t.add_twire(vec![e]);
        let f = t.add_tface(None, w, vec![], None, None, vec![], true);
        let sh = t.add_tshell(vec![f]);
        let s0 = t.add_tsolid(vec![sh]);

        let v2 = t.add_tvertex(DVec3::Y);
        let e2 = t.add_tedge(None, v1, v2, [0.0, 1.0]);
        let w2 = t.add_twire(vec![e2]);
        let f2 = t.add_tface(None, w2, vec![], None, None, vec![], true);
        let sh2 = t.add_tshell(vec![f2]);
        let s1 = t.add_tsolid(vec![sh2]);

        // Register solids in result.solids
        rb.solids.push(s0);
        rb.solids.push(s1);

        // Build compsolid from [0, 1]
        rb.build_compsolids(&mut t, vec![vec![0, 1]]);

        assert_eq!(rb.compsolid_groups.len(), 1);
        // t should contain the compsolid TShape (in addition to previous shapes)
        let n_compsolid = t.tshapes.iter().filter(|ts| matches!(&***ts, topods::TShape::CompSolid(_))).count();
        assert_eq!(n_compsolid, 1);
    }

    #[test]
    fn test_build_compounds_creates_compound() {
        let mut rb = ResultBuilder::new();
        let mut t = topods::BRep::new();
        // Create two solids
        let v0 = t.add_tvertex(DVec3::ZERO);
        let v1 = t.add_tvertex(DVec3::X);
        let e = t.add_tedge(None, v0, v1, [0.0, 1.0]);
        let w = t.add_twire(vec![e]);
        let f = t.add_tface(None, w, vec![], None, None, vec![], true);
        let sh = t.add_tshell(vec![f]);
        let s0 = t.add_tsolid(vec![sh]);

        let v2 = t.add_tvertex(DVec3::Y);
        let e2 = t.add_tedge(None, v1, v2, [0.0, 1.0]);
        let w2 = t.add_twire(vec![e2]);
        let f2 = t.add_tface(None, w2, vec![], None, None, vec![], true);
        let sh2 = t.add_tshell(vec![f2]);
        let s1 = t.add_tsolid(vec![sh2]);

        rb.solids.push(s0);
        rb.solids.push(s1);

        rb.build_compounds(&mut t, &[vec![0, 1]]);

        let n_compound = t.tshapes.iter().filter(|ts| matches!(&***ts, topods::TShape::Compound(_))).count();
        assert_eq!(n_compound, 1);
    }
}
