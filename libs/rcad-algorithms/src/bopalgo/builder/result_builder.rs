use crate::bopalgo::builder::SourceSide;
use crate::bopalgo::builder::types::{
    FaceEntry, WireEdgeSource, WireEdgeSourceTopoDS, WireFace, WireSegment, WireSegmentTopoDS,
};
use crate::bopalgo::builder::{curve_eq, hash_point};
use crate::bopds::ds::*;
use crate::history::{
    BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin,
};
use crate::tolerance::*;
use crate::triangulate::{triangulate_polygon, triangulate_polygon_with_holes};
use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topods;
use rcad_kernel::topology::*;
use std::collections::{HashMap, HashSet};

/// Builds result BRep from accumulated DS face data.
///
/// pure conversion --BuildResult does no dedup/merge/cull.
#[allow(dead_code)]
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
    pub(crate) custom_edge_ranges: Vec<Option<[f64; 2]>>,
    /// Reference to DS edge array for looking up curve data of section edges
    /// (architecture diff A6).  Set by split_face_and_emit_topo_ds before use.
    pub(crate) ds_edges: Option<std::sync::Arc<Vec<crate::bopds::ds::DSEdge>>>,
    pub(crate) face_internal_vtx: Vec<Vec<usize>>,
    pub(crate) deg_edge_indices: HashSet<usize>,
    pub(crate) ic_edge_map: HashMap<usize, usize>,
    /// DS edge index --flat edge index.
    /// Two faces sharing the same DS edge get the same flat edge index,
    /// matching OCCT's TopoDS_Edge identity sharing (same TShape* pointer).
    /// Populated by emit_wire_face_topods for DsEdge-sourced wire segments.
    pub(crate) ds_edge_to_flat: HashMap<usize, usize>,
    pub(crate) source_has_compound: bool,
    pub(crate) tmp_compsolid_groups: Vec<Vec<usize>>,
    /// per-source side tracking for solids (0=ShapeA/Args, 1=ShapeB/Tools).
    /// Parallel to tmp_solids: solid_side_origin[i] = side that produced tmp_solids[i].
    pub(crate) solid_side_origin: Vec<usize>,
    /// compound groups (FillImagesCompound output).
    /// Each entry is a Vec of solid indices (into result.solids) forming one compound.
    pub(crate) compound_groups: Vec<Vec<usize>>,
    /// natural_restriction for each face in self.faces.
    /// Parallel to face_origins.  Populated by emit_wire_face / build_original_face.
    pub(crate) face_natural_restriction: Vec<bool>,
    /// maps result face index --array of DSWire indices for all its wires
    /// (outer first, then inners).  Parallel to face_origins.  Used by
    /// build_topods_faces to reference pre-built TShape::Wire from wire_refs.
    pub(crate) face_all_wire_idxs: Vec<Vec<usize>>,
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
        if len > TOLERANCE_LEN_MIN {
            n / len
        } else {
            DVec3::ZERO
        }
    }

    /// emit BRep face from WireFace (replaces emit_face_with_origin).
    /// Builds edges directly from WireSegments: seam edges use add_seam_edge /
    /// add_edge_seam_degenerate; IC edges use add_circle_edge for Circle3 curves.
    /// --emit_wire_face --builds BRep edges/face from WireSegments.
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
        let ow: Vec<&usize> = wf
            .outer_wire
            .iter()
            .filter(|&&si| segments[si].start_vertex != segments[si].end_vertex)
            .collect();
        for &&si in &ow {
            let seg = &segments[si];
            // --canonical vertices use stored positions
            let get_pos = |vi: usize| -> DVec3 {
                vertex_positions
                    .get(&vi)
                    .copied()
                    .unwrap_or(ds.vertex_point(vi))
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
                let seam_deg = (get_pos(seg.start_vertex) - get_pos(seg.end_vertex))
                    .length_squared()
                    < TOLERANCE_ABS_SQ;
                let sphere_surf = match &ds.faces[face_idx].surface {
                    Surface3::Sphere(s) => s,
                    _ => &SphericalSurface {
                        center: DVec3::ZERO,
                        axis: DVec3::Z,
                        radius: 1.0,
                        ref_dir: DVec3::X,
                    },
                };
                // --canonical deg edges (vertex >= ds.vertices.len())
                let is_canon_deg =
                    seg.start_vertex >= ds.vertices.len() || seg.end_vertex >= ds.vertices.len();
                let ei = if seam_deg || is_canon_deg {
                    self.add_edge_seam_degenerate(v1, v2, sphere_surf)
                } else {
                    let seam_normal = any_perpendicular(sphere_surf.axis).normalize();
                    let seam_circle = Curve3::Circle(Circle3::new(
                        sphere_surf.center,
                        seam_normal,
                        sphere_surf.radius,
                    ));
                    self.add_seam_edge(v1, v2, seam_circle)
                };
                (ei, true)
            } else {
                let ei = match &seg.source {
                    // --IC edge identity (section edges shared).
                    WireEdgeSource::IntersectionCurve(ci) => {
                        let crv = &ds.intersection_curves[*ci].curve;
                        self.add_ic_edge(*ci, v1, v2, crv.clone(), Some(seg.t_range))
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
                let getp = |vi: usize| -> DVec3 {
                    if vi < ds.vertices.len() {
                        ds.vertex_point(vi)
                    } else {
                        *vertex_positions.get(&vi).unwrap_or(&DVec3::ZERO)
                    }
                };
                let v1 = if seg.start_vertex < ds.vertices.len() {
                    self.add_ds_vertex(seg.start_vertex, ds.vertices[seg.start_vertex].point)
                } else {
                    let p = getp(seg.start_vertex);
                    self.add_vertex(p)
                };
                let v2 = if seg.end_vertex < ds.vertices.len() {
                    self.add_ds_vertex(seg.end_vertex, ds.vertices[seg.end_vertex].point)
                } else {
                    let p = getp(seg.end_vertex);
                    self.add_vertex(p)
                };
                if iw_verts.is_empty() || iw_verts.last() != Some(&v1) {
                    iw_verts.push(v1);
                }
                let ei = match &seg.source {
                    // --IC edge identity (inner/internal wires).
                    WireEdgeSource::IntersectionCurve(ci) => {
                        let crv = &ds.intersection_curves[*ci].curve;
                        self.add_ic_edge(*ci, v1, v2, crv.clone(), Some(seg.t_range))
                    }
                    WireEdgeSource::DsEdge(_) | WireEdgeSource::SeamEdge => self.add_edge(v1, v2),
                };
                let forward = self.edges[ei].0 == v1;
                iw_edges.push((ei, forward));
            }
            inner_wire_edges.push(iw_edges);
            iw_vert_indices_all.extend(iw_verts);
        }

        // --Internal wire edges (TopAbs_INTERNAL).
        // Seam edges use add_seam_edge for curve-aware unique identity.
        let mut internal_wire_edges: Vec<Vec<(usize, bool)>> = Vec::new();
        for iw in &wf.internal_wires {
            let mut iw_edges = Vec::new();
            for &si in iw {
                let seg = &segments[si];
                let getp = |vi: usize| -> DVec3 {
                    if vi < ds.vertices.len() {
                        ds.vertex_point(vi)
                    } else {
                        *vertex_positions.get(&vi).unwrap_or(&DVec3::ZERO)
                    }
                };
                let v1 = if seg.start_vertex < ds.vertices.len() {
                    self.add_ds_vertex(seg.start_vertex, ds.vertices[seg.start_vertex].point)
                } else {
                    let p = getp(seg.start_vertex);
                    self.add_vertex(p)
                };
                let v2 = if seg.end_vertex < ds.vertices.len() {
                    self.add_ds_vertex(seg.end_vertex, ds.vertices[seg.end_vertex].point)
                } else {
                    let p = getp(seg.end_vertex);
                    self.add_vertex(p)
                };
                let ei = match &seg.source {
                    WireEdgeSource::IntersectionCurve(ci) => {
                        let crv = &ds.intersection_curves[*ci].curve;
                        if let Curve3::Circle(_) = crv {
                            self.add_circle_edge(v1, v2, crv.clone())
                        } else {
                            self.add_edge(v1, v2)
                        }
                    }
                    WireEdgeSource::DsEdge(_) | WireEdgeSource::SeamEdge
                        if seg.is_closed_on_face =>
                    {
                        let s = match &ds.faces[face_idx].surface {
                            Surface3::Sphere(sph) => sph,
                            _ => &SphericalSurface {
                                center: DVec3::ZERO,
                                axis: DVec3::Z,
                                radius: 1.0,
                                ref_dir: DVec3::X,
                            },
                        };
                        let c = Curve3::Circle(Circle3::new(
                            s.center,
                            any_perpendicular(s.axis).normalize(),
                            s.radius,
                        ));
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
        let iw_boundaries: Vec<Vec<DVec3>> = inner_wire_edges
            .iter()
            .map(|iw_es| {
                let mut pts = Vec::new();
                for &(ei, _) in iw_es {
                    let (a, b) = self.edges[ei];
                    if pts.is_empty() || pts.last() != Some(&a) {
                        pts.push(a);
                    }
                }
                pts.iter().map(|&vi| self.vertices[vi]).collect()
            })
            .collect();
        let all_vert_indices: Vec<usize> =
            [vert_indices.as_slice(), iw_vert_indices_all.as_slice()].concat();
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
        let centroid =
            outer_boundary.iter().copied().sum::<DVec3>() / outer_boundary.len().max(1) as f64;
        let area = Self::polygon_signed_area_on_normal(&outer_boundary, normal);
        let mut outer_sig: Vec<usize> = edge_indices.iter().map(|&(eid, _)| eid).collect();
        outer_sig.sort_unstable();
        let nlen = normal.length();
        let nunit = if nlen > TOLERANCE_LEN_MIN {
            normal / nlen
        } else {
            normal
        };
        for (
            existing_idx,
            (
                existing_outer,
                existing_inner,
                _existing_tris,
                existing_normal,
                _surf,
                _uv,
                existing_centroid,
                existing_area,
                _existing_sp,
                _existing_iw,
            ),
        ) in self.faces.iter().enumerate()
        {
            let mut ex_sig: Vec<usize> = existing_outer.iter().map(|&(eid, _)| eid).collect();
            for iw_edges in existing_inner {
                ex_sig.extend(iw_edges.iter().map(|&(eid, _)| eid));
            }
            ex_sig.sort_unstable();
            let elen = existing_normal.length();
            if elen <= TOLERANCE_LEN_MIN {
                continue;
            }
            let eunit = *existing_normal / elen;
            let sig_match = ex_sig == outer_sig;
            let geo_match = nunit.dot(eunit).abs() >= 0.99
                && (*existing_centroid - centroid).length() <= TOLERANCE_LINEAR_RELAX_8
                && (existing_area - area).abs()
                    <= TOLERANCE_LINEAR_RELAX_8 * existing_area.max(area).max(1.0);
            if sig_match || geo_match {
                self.co_face_origins.push((existing_idx, origin));
                return;
            }
        }

        // --No extra internal vertices needed --wire pipeline handles
        // seam edges via WireSegment virtual edges; BuilderFace does not add
        // degenerate vertices to the result face.

        // Compute UV domain for sphere faces
        let sphere_uv = if matches!(face.surface, Surface3::Sphere(_)) {
            let uvs: Vec<DVec2> = if !wf.outer_wire.is_empty() {
                wf.outer_wire
                    .iter()
                    .map(|&si| {
                        let seg = &segments[si];
                        let sph = match &face.surface {
                            Surface3::Sphere(s) => s,
                            _ => unreachable!(),
                        };
                        sph.world_to_uv(ds.vertices[seg.start_vertex].point)
                    })
                    .collect()
            } else {
                vec![]
            };
            if !uvs.is_empty() {
                let u_min = uvs.iter().map(|uv| uv.x).fold(f64::INFINITY, f64::min);
                let u_max = uvs.iter().map(|uv| uv.x).fold(f64::NEG_INFINITY, f64::max);
                let v_min = uvs.iter().map(|uv| uv.y).fold(f64::INFINITY, f64::min);
                let v_max = uvs.iter().map(|uv| uv.y).fold(f64::NEG_INFINITY, f64::max);
                if (u_max - u_min).abs() > TOLERANCE_FLOAT_LOOSE
                    && (v_max - v_min).abs() > TOLERANCE_FLOAT_LOOSE
                {
                    Some([u_min, u_max, v_min, v_max])
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };

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
        self.face_natural_restriction
            .push(ds.face_natural_restriction(face_idx));
        {
            let mut widxs = Vec::new();
            if let Some(owi) = ds.face_outer_wire_idx(face_idx) {
                widxs.push(owi);
            }
            widxs.extend(&ds.faces[face_idx].inner_wire_idxs);
            self.face_all_wire_idxs.push(widxs);
        }
    }

    /// --emit_wire_face using WireSegmentTopoDS.
    /// Same logic as emit_wire_face but reads edge/vertex data through BRepTool.
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
        ds_ei_to_sr: &HashMap<usize, topods::ShapeRef>,
        sr_index_to_ds_ei: &HashMap<usize, usize>,
        ds: &crate::bopds::ds::DS,
    ) {
        let get_pos = |vi: usize| -> DVec3 {
            // use DS vertex world coordinate (TopLoc_Location already baked
            // by from_topods_with_location / load_brep).  The BRepTool fallback is only
            // used when vi is NOT a DS vertex index.
            if vi < ds.vertices.len() {
                ds.vertex_point(vi)
            } else {
                vertex_positions.get(&vi).copied().unwrap_or_else(|| {
                    tool.vertex_position(rcad_kernel::topods::ShapeRef::synthetic(vi))
                })
            }
        };

        let mut vert_indices = Vec::new();
        let mut edge_indices = Vec::new();
        let ow: Vec<&usize> = wf
            .outer_wire
            .iter()
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
                let seam_deg = (get_pos(seg.start_vertex.index) - get_pos(seg.end_vertex.index))
                    .length_squared()
                    < TOLERANCE_ABS_SQ;
                let sphere_surf = match tool.face_surface(face_ref) {
                    Some(Surface3::Sphere(s)) => s.clone(),
                    _ => SphericalSurface {
                        center: DVec3::ZERO,
                        axis: DVec3::Z,
                        radius: 1.0,
                        ref_dir: DVec3::X,
                    },
                };
                let is_canon_deg = vertex_positions.contains_key(&seg.start_vertex.index)
                    || vertex_positions.contains_key(&seg.end_vertex.index);
                let ei = if seam_deg || is_canon_deg {
                    self.add_edge_seam_degenerate(v1, v2, &sphere_surf)
                } else {
                    let seam_normal = any_perpendicular(sphere_surf.axis).normalize();
                    let seam_circle = Curve3::Circle(Circle3::new(
                        sphere_surf.center,
                        seam_normal,
                        sphere_surf.radius,
                    ));
                    self.add_seam_edge(v1, v2, seam_circle)
                };
                (ei, true)
            } else {
                let ei = match &seg.source {
                    super::types::WireEdgeSourceTopoDS::IntersectionCurve(ci) => {
                        let crv =
                            ic_curves
                                .get(&ci.index)
                                .cloned()
                                .unwrap_or(Curve3::Line(Line3 {
                                    origin: DVec3::ZERO,
                                    direction: DVec3::X,
                                }));
                        self.add_ic_edge(ci.index, v1, v2, crv, Some(seg.t_range))
                    }
                    super::types::WireEdgeSourceTopoDS::DsEdge(sr) => {
                        // Architecture diff A6: section edges created by MakeBlocks
                        // have their curve stored in the DS edge.
                        // sr.index = e_base + ds_ei; use sr_index_to_ds_ei reverse
                        // lookup to get the original DS edge index.
                        if let Some(&ds_ei) = sr_index_to_ds_ei.get(&sr.index) {
                            // same DS edge --same flat edge index
                            // (matching TopoDS_Edge identity sharing).
                            if let Some(&existing) = self.ds_edge_to_flat.get(&ds_ei) {
                                existing
                            } else if let Some(ref ds_edges) = self.ds_edges {
                                if ds_ei < ds_edges.len() {
                                    let crv = ds_edges[ds_ei].curve.clone();
                                    let range = ds_edges[ds_ei].t_range;
                                    let ei = self.add_edge_with_curve(v1, v2, crv, range);
                                    self.ds_edge_to_flat.insert(ds_ei, ei);
                                    ei
                                } else {
                                    let ei = self.add_edge(v1, v2);
                                    self.ds_edge_to_flat.insert(ds_ei, ei);
                                    ei
                                }
                            } else {
                                let ei = self.add_edge(v1, v2);
                                self.ds_edge_to_flat.insert(ds_ei, ei);
                                ei
                            }
                        } else {
                            self.add_edge(v1, v2)
                        }
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
                        let crv =
                            ic_curves
                                .get(&ci.index)
                                .cloned()
                                .unwrap_or(Curve3::Line(Line3 {
                                    origin: DVec3::ZERO,
                                    direction: DVec3::X,
                                }));
                        self.add_ic_edge(ci.index, v1, v2, crv, Some(seg.t_range))
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
                        let crv =
                            ic_curves
                                .get(&ci.index)
                                .cloned()
                                .unwrap_or(Curve3::Line(Line3 {
                                    origin: DVec3::ZERO,
                                    direction: DVec3::X,
                                }));
                        if let Curve3::Circle(_) = &crv {
                            self.add_circle_edge(v1, v2, crv)
                        } else {
                            self.add_edge(v1, v2)
                        }
                    }
                    _ if seg.is_closed_on_face => {
                        let sphere_surf = match tool.face_surface(face_ref) {
                            Some(Surface3::Sphere(s)) => s.clone(),
                            _ => SphericalSurface {
                                center: DVec3::ZERO,
                                axis: DVec3::Z,
                                radius: 1.0,
                                ref_dir: DVec3::X,
                            },
                        };
                        let c = Curve3::Circle(Circle3::new(
                            sphere_surf.center,
                            any_perpendicular(sphere_surf.axis).normalize(),
                            sphere_surf.radius,
                        ));
                        self.add_seam_edge(v1, v2, c)
                    }
                    _ => self.add_edge(v1, v2),
                };
                let forward = self.edges[ei].0 == v1;
                iw_edges.push((ei, forward));
            }
            internal_wire_edges.push(iw_edges);
        }

        // OCCT L611-612: store face (no triangulation, no dedup).
        let surface = tool
            .face_surface(face_ref)
            .cloned()
            .unwrap_or(Surface3::Plane(rcad_kernel::geom::Plane::new(
                DVec3::ZERO,
                DVec3::Z,
            )));
        let sample_pt = if !wf.outer_wire.is_empty() {
            get_pos(segments[wf.outer_wire[0]].start_vertex.index)
        } else {
            DVec3::ZERO
        };
        self.face_internal_vtx.push(Vec::new());
        self.faces.push((
            edge_indices,
            inner_wire_edges,
            vec![],
            DVec3::Z,
            surface,
            None,
            DVec3::ZERO,
            0.0,
            sample_pt,
            internal_wire_edges,
        ));
        self.face_origins.push(origin);
        self.face_natural_restriction.push(natural_restriction);
        {
            let mut widxs = Vec::new();
            if let Some(owi) = ds.face_outer_wire_idx(face_idx) {
                widxs.push(owi);
            }
            widxs.extend(&ds.faces[face_idx].inner_wire_idxs);
            self.face_all_wire_idxs.push(widxs);
        }
    }

    /// --estimate face normal from wire segments (TopoDS variant).
    fn estimate_boundary_normal_from_segments_topo(
        outer_wire: &[usize],
        segments: &[super::types::WireSegmentTopoDS],
        tool: &dyn rcad_kernel::topods::BRepTool,
        vertex_positions: &HashMap<usize, DVec3>,
    ) -> DVec3 {
        if outer_wire.len() < 3 {
            return DVec3::ZERO;
        }
        let pts: Vec<DVec3> = outer_wire
            .iter()
            .map(|&si| {
                let seg = &segments[si];
                let vi = seg.start_vertex.index;
                vertex_positions.get(&vi).copied().unwrap_or_else(|| {
                    tool.vertex_position(rcad_kernel::topods::ShapeRef::synthetic(vi))
                })
            })
            .collect();
        Self::estimate_boundary_normal(&pts)
    }

    fn estimate_boundary_normal_from_segments(
        outer_wire: &[usize],
        segments: &[WireSegment],
        ds: &DS,
    ) -> DVec3 {
        if outer_wire.len() < 3 {
            return DVec3::ZERO;
        }
        let pts: Vec<DVec3> = outer_wire
            .iter()
            .map(|&si| {
                let seg = &segments[si];
                ds.vertices[seg.start_vertex].point
            })
            .collect();
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
            custom_edge_ranges: Vec::new(),
            ds_edges: None,
            face_internal_vtx: Vec::new(),
            deg_edge_indices: std::collections::HashSet::new(),
            ic_edge_map: HashMap::new(),
            ds_edge_to_flat: HashMap::new(),
            tmp_shells: Vec::new(),
            tmp_solids: Vec::new(),
            source_has_compound: false,
            tmp_compsolid_groups: Vec::new(),
            solid_side_origin: Vec::new(),
            compound_groups: Vec::new(),
            face_natural_restriction: Vec::new(),
            face_all_wire_idxs: Vec::new(),
        }
    }

    /// --BuildResult(EDGE) --build edges from split_edges.
    /// --BuildResult(FACE) --build faces from accumulated face data.
    /// OCCT Builder_1.cxx L130-168: iterate myImages for TopAbs_FACE, add to myShape.
    /// rcad: build faces from self.faces, referencing already-built self.edges.
    /// Maps each face's per-vertex-pair edges to the BRep edge indices from build_edges.
    /// --BuildResult(FACE) --build faces from accumulated face data.
    /// OCCT Builder_1.cxx L130-168: iterate myImages for TopAbs_FACE, add to myShape.
    /// rcad: validate face edge refs against built edges, prepare for shell/solid assembly.
    pub(crate) fn build_faces(&mut self) {
        // Validate that all face edge references are within bounds of built edges
        let n_edges = self.edges.len();
        for (fi, face) in self.faces.iter().enumerate() {
            for &(ei, _) in &face.0 {
                assert!(
                    ei < n_edges,
                    "face[{}] edge ref {} out of range ({} edges)",
                    fi,
                    ei,
                    n_edges
                );
            }
            for iw in &face.1 {
                for &(ei, _) in iw {
                    assert!(
                        ei < n_edges,
                        "face[{}] inner edge ref {} out of range",
                        fi,
                        ei
                    );
                }
            }
        }
    }

    /// --BuildResult(FACE) --add unmodified source face.
    /// --BuildResult(FACE) --add original source face (Builder_1.cxx L146-152).
    /// OCCT adds the original TopoDS_Face regardless of surface type.
    /// rcad: builds FaceEntry from DS boundary_edges + inner_boundary_edges.
    /// Handles all surface types (Plane, Cylinder, Sphere, Cone, Torus).
    /// --BuildResult(FACE) --add original faces without images.
    /// Now creates TShape::Face directly (OCCT: adds existing TopoDS_Face to myShape).
    pub(crate) fn build_original_face(
        &mut self,
        ds: &DS,
        fi: usize,
        origin: FaceOrigin,
        t: &mut topods::BRep,
        face_refs: &mut Vec<topods::ShapeRef>,
    ) {
        let face = &ds.faces[fi];

        // --- Outer wire from boundary_edges, creating TShape vertices/edges directly ---
        let e_base = ds.vertices.len();
        let mut outer_edges: Vec<topods::ShapeRef> = Vec::new();
        let mut prev_end: Option<usize> = None;
        for &ei in &face.boundary_edges {
            if ei >= ds.edges.len() {
                continue;
            }
            let e = &ds.edges[ei];
            let (sv, ev) = match prev_end {
                Some(pe) if e.start_vertex == pe => (e.start_vertex, e.end_vertex),
                Some(pe) if e.end_vertex == pe => (e.end_vertex, e.start_vertex),
                _ => (e.start_vertex, e.end_vertex),
            };
            let sv_sr = t.add_tvertex(ds.vertex_point(sv)).with_location(e.location);
            let ev_sr = t.add_tvertex(ds.vertex_point(ev)).with_location(e.location);
            let e_sr = t
                .add_tedge(Some(e.curve.clone()), sv_sr, ev_sr, e.t_range)
                .with_location(e.location);
            outer_edges.push(e_sr);
            prev_end = Some(ev);
        }
        if outer_edges.len() < 3 {
            return;
        }
        let outer_wire = t.add_twire(outer_edges);

        // --- Inner wires ---
        let mut inner_wires: Vec<topods::ShapeRef> = Vec::new();
        for iw_edges in &face.inner_boundary_edges {
            let mut wire_edges: Vec<topods::ShapeRef> = Vec::new();
            for &(ei, forward_in_ds) in iw_edges {
                if ei >= ds.edges.len() {
                    continue;
                }
                let e = &ds.edges[ei];
                let (sv, ev) = if forward_in_ds {
                    (e.start_vertex, e.end_vertex)
                } else {
                    (e.end_vertex, e.start_vertex)
                };
                let sv_sr = t.add_tvertex(ds.vertex_point(sv)).with_location(e.location);
                let ev_sr = t.add_tvertex(ds.vertex_point(ev)).with_location(e.location);
                let e_sr = t
                    .add_tedge(Some(e.curve.clone()), sv_sr, ev_sr, e.t_range)
                    .with_location(e.location);
                wire_edges.push(e_sr);
            }
            if wire_edges.len() >= 2 {
                inner_wires.push(t.add_twire(wire_edges));
            }
        }

        let sample_pt = ds.vertices[face.boundary_verts.first().copied().unwrap_or(0)].point;
        let face_sr = t
            .add_tface(
                Some(face.surface.clone()),
                outer_wire,
                inner_wires,
                Some(sample_pt),
                None,
                vec![],
                ds.face_natural_restriction(fi),
            )
            .with_location(face.location);
        face_refs.push(face_sr);
        self.face_origins.push(origin);

        // Legacy flat-index path (kept for downstream consumers).
        // TODO: remove when downstream is migrated.
        {
            self.face_natural_restriction
                .push(ds.face_natural_restriction(fi));
            {
                let mut widxs = Vec::new();
                if let Some(owi) = ds.face_outer_wire_idx(fi) {
                    widxs.push(owi);
                }
                widxs.extend(&ds.faces[fi].inner_wire_idxs);
                self.face_all_wire_idxs.push(widxs);
            }
        }
    }

    /// --BuildResult(COMPSOLID) --build compsolids via BRepBuilder.
    /// OCCT: BOPAlgo_Builder::BuildResult (Builder_1.cxx L130-168) iterates
    /// source COMPSOLID shapes and adds their split images to myShape via
    /// BRep_Builder::Add.  rcad: processes tmp_compsolid_groups (groups of
    /// solid indices from fill_images_container_compsolid) and creates
    /// topods compsolids using BRepBuilder::make_compsolid.
    pub(crate) fn build_compsolids(
        &mut self,
        t: &mut topods::BRep,
        groups: Vec<Vec<usize>>,
        solids: &[topods::ShapeRef],
        compsolid_groups: &mut Vec<topods::ShapeRef>,
    ) {
        let mut bb = topods::BRepBuilder::new();
        for cs_group in &groups {
            let solid_refs: Vec<topods::ShapeRef> = cs_group
                .iter()
                .filter_map(|&si| solids.get(si).copied())
                .collect();
            if !solid_refs.is_empty() {
                compsolid_groups.push(bb.make_compsolid(t, solid_refs));
            }
        }
    }

    /// --BuildResult(COMPOUND) --build compounds via BRepBuilder.
    /// OCCT: BOPAlgo_Builder::BuildResult (Builder_1.cxx L130-168) iterates
    /// source COMPOUND shapes and adds their split images to myShape via
    /// BRep_Builder::Add.  rcad: processes compound_groups (groups of solid
    /// indices from fill_images_compounds) and creates topods compounds
    /// using BRepBuilder::make_compound.
    pub(crate) fn build_compounds(
        &mut self,
        t: &mut topods::BRep,
        groups: &[Vec<usize>],
        solids: &[topods::ShapeRef],
    ) {
        let mut bb = topods::BRepBuilder::new();
        for group in groups {
            let solid_refs: Vec<topods::ShapeRef> = group
                .iter()
                .filter_map(|&si| solids.get(si).copied())
                .collect();
            if !solid_refs.is_empty() {
                bb.make_compound(t, solid_refs);
            }
        }
    }

    pub(crate) fn add_vertex(&mut self, point: DVec3) -> usize {
        let key = hash_point(point);
        if let Some(&idx) = self.vertex_map.get(&key) {
            if points_coincide(self.vertices[idx], point) {
                return idx;
            }
        }
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

    /// --add vertex by DS index identity (TopoDS_Vertex TShape).
    pub(crate) fn add_ds_vertex(&mut self, ds_vi: usize, point: DVec3) -> usize {
        if let Some(&idx) = self.ds_vertex_map.get(&ds_vi) {
            return idx;
        }
        let idx = self.add_vertex(point);
        self.ds_vertex_map.insert(ds_vi, idx);
        idx
    }

    /// Geometric edge key for edge-set matching (BOPTools_Set analog).
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
        h ^= ka.0 as u64;
        h = h.wrapping_mul(1099511628211);
        h ^= ka.1 as u64;
        h = h.wrapping_mul(1099511628211);
        h ^= ka.2 as u64;
        h = h.wrapping_mul(1099511628211);
        h ^= kb.0 as u64;
        h = h.wrapping_mul(1099511628211);
        h ^= kb.1 as u64;
        h = h.wrapping_mul(1099511628211);
        h ^= kb.2 as u64;
        h = h.wrapping_mul(1099511628211);
        h
    }

    /// --BRep_Builder::MakeEdge --creates new unique edge.
    /// OCCT: each TopoDS_Edge is a distinct entity (per TShape identity).
    /// Even edges connecting the same vertices are distinct TopoDS_Edges.
    /// rcad: always appends a new edge, same semantics.
    pub(crate) fn add_edge_occt(&mut self, v1: usize, v2: usize) -> usize {
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        idx
    }

    /// --BOPTools_AlgoTools::MakeSectEdge --shared section edge.
    /// OCCT: MakeSectEdge creates ONE TopoDS_Edge that both intersecting faces
    /// reference via BRep_Builder::Add (shared TShape identity).
    /// rcad: maps intersection curve index --result edge index so both faces
    /// emit_wire_face calls get the same edge index for the same IC curve.
    /// OCCT: each TopoDS_Edge is a distinct handle --no post-hoc merge needed.
    pub(crate) fn add_ic_edge(
        &mut self,
        ici: usize,
        v1: usize,
        v2: usize,
        curve: Curve3,
        range: Option<[f64; 2]>,
    ) -> usize {
        if let Some(&idx) = self.ic_edge_map.get(&ici) {
            let existing = self.edges[idx];
            if (existing.0 != v1 || existing.1 != v2) && (existing.0 != v2 || existing.1 != v1) {
                eprintln!(
                    "[IC_VTX] ci={} existing=({}, {}) called=({}, {})",
                    ici, existing.0, existing.1, v1, v2
                );
            }
            return idx;
        }
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        self.custom_edge_curves[idx] = Some(curve);
        while self.custom_edge_ranges.len() <= idx {
            self.custom_edge_ranges.push(None);
        }
        self.custom_edge_ranges[idx] = range;
        self.ic_edge_map.insert(ici, idx);
        idx
    }

    /// --BRep_Builder::Add edge sharing --dedup by (v1,v2) pair.
    /// OCCT: BRep_Builder::Add(theSameEdge, faceA) then Add(theSameEdge, faceB)
    /// shares the same TopoDS_Edge between faces (TShape identity).
    /// rcad: add_edge(v1,v2) returns the same index for the same vertex pair,
    /// achieving the same sharing without requiring TopoDS shape handles.
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

    /// Add edge with known curve and parametric range (OCCT BRep_Builder::Add).
    /// Like add_edge but also stores the 3D curve geometry and its parameter
    /// range so that downstream consumers (area computation, STEP export) can
    /// evaluate the curve correctly without relying on hardcoded [0,1].
    pub(crate) fn add_edge_with_curve(
        &mut self,
        v1: usize,
        v2: usize,
        curve: Curve3,
        range: [f64; 2],
    ) -> usize {
        // No dedup by (v1,v2): two different curves may share endpoints.
        // Sharing of the same DSEdge between faces is handled upstream
        // (same segment generates the same edge index via add_ic_edge or
        //  the fallback add_edge dedup for uncurved edges).
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        self.custom_edge_curves[idx] = Some(curve);
        while self.custom_edge_ranges.len() <= idx {
            self.custom_edge_ranges.push(None);
        }
        self.custom_edge_ranges[idx] = Some(range);
        idx
    }

    /// --create degenerate seam edge with hemisphere circle curve.
    /// OCCT sphere face outer wire always has a degenerate seam edge (same vertex at both ends).
    /// Adds a sphere horizontal circle curve to make the edge recognizable in STEP export.
    pub(crate) fn add_edge_seam_degenerate(
        &mut self,
        v1: usize,
        v2: usize,
        sphere_surf: &SphericalSurface,
    ) -> usize {
        let idx = self.edges.len();
        self.edges.push((v1, v2));
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        // Store seam circle curve for STEP writer
        // --seam edge = sphere meridian through pole (not IC circle).
        // If normal = axis, it would coincide with plane-sphere IC causing curve merge errors.
        let seam_normal = any_perpendicular(sphere_surf.axis).normalize();
        let seam_circle = Curve3::Circle(Circle3::new(
            sphere_surf.center,
            seam_normal,
            sphere_surf.radius,
        ));
        self.custom_edge_curves[idx] = Some(seam_circle);
        self.deg_edge_indices.insert(idx);
        idx
    }

    /// --circle edge with curve-aware dedup.
    /// OCCT: TopoDS_Edge identity is per-TShape, not per vertex pair.
    /// Two edges sharing vertices but with different curves are distinct.
    /// rcad: dedup by both (v1,v2) AND curve identity (Circle3 geometry).
    pub(crate) fn add_circle_edge(&mut self, v1: usize, v2: usize, circle: Curve3) -> usize {
        let key = (v1.min(v2), v1.max(v2));
        for (i, e) in self.edges.iter().enumerate() {
            if (e.0.min(e.1), e.0.max(e.1)) == key {
                if let Some(ref existing) = self.custom_edge_curves.get(i).and_then(|c| c.as_ref())
                {
                    // Different curve at same vertex pair --distinct TopoDS_Edge
                    if !curve_eq(existing, &circle) {
                        let idx = self.add_edge_occt(v1, v2);
                        while self.custom_edge_curves.len() <= idx {
                            self.custom_edge_curves.push(None);
                        }
                        self.custom_edge_curves[idx] = Some(circle);
                        return idx;
                    }
                }
                // Same curve or no existing curve --reuse
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
    /// --BOPTools_AlgoTools::MakeEdge --  , --
    /// add_edge_occt, --
    /// seam IC OCCT TopoDS_Edge--
    pub(crate) fn add_circle_edge_occt(&mut self, v1: usize, v2: usize, circle: Curve3) -> usize {
        let idx = self.add_edge_occt(v1, v2);
        while self.custom_edge_curves.len() <= idx {
            self.custom_edge_curves.push(None);
        }
        self.custom_edge_curves[idx] = Some(circle);
        // no_merge_edges removed --edges are inherently unique by index
        idx
    }

    /// --MakeEdge for seam edges (BRep_Builder::MakeEdge pattern).
    /// OCCT: BRep_Builder::MakeEdge creates a TopoDS_Edge with the 3D curve.
    /// Seam edges and IC arcs at the same vertex pair are distinct TopoDS_Edges
    /// (different TShapes).  rcad: same vertex pair + same curve --reuse (shared
    /// TShape); same vertex pair + different curve --create new via add_edge_occt
    /// (distinct TShape).  This matches OCCT's per-TShape edge identity.
    pub(crate) fn add_seam_edge(&mut self, v1: usize, v2: usize, circle: Curve3) -> usize {
        // Same logic as add_circle_edge: check for existing edge with same
        // vertex pair but different curve --create new; same curve --reuse.
        let key = (v1.min(v2), v1.max(v2));
        for (i, e) in self.edges.iter().enumerate() {
            if (e.0.min(e.1), e.0.max(e.1)) == key {
                if let Some(ref existing) = self.custom_edge_curves.get(i).and_then(|c| c.as_ref())
                {
                    if !curve_eq(existing, &circle) {
                        let idx = self.add_edge_occt(v1, v2);
                        while self.custom_edge_curves.len() <= idx {
                            self.custom_edge_curves.push(None);
                        }
                        self.custom_edge_curves[idx] = Some(circle);
                        // no_merge_edges removed --edges are inherently unique by index
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
        // no_merge_edges removed --edges are inherently unique by index
        idx
    }

    /// DEPRECATED (FaceSampleData  ):  = --?  split_planar_face  -- ?
    /// OCCT: MakeBlocks  ?BOPTools_AlgoTools::MakeEdge(aIC,...)
    /// split_planar_face  ?28+ ? --2 --(arc_simplify),
    /// emit_face_with_origin  ?add_circle_edge  --Circle3  --

    /// --Architecture A1 --emit TShape for the face just added to self.faces.
    /// OCCT BRep_Builder creates edges/wires/faces incrementally during BuildSplitFaces.
    /// rcad previously deferred this to build_topods_faces; now creates TShapes per-face.
    pub(crate) fn emit_face_topods(
        &mut self,
        t: &mut topods::BRep,
        face_refs: &mut Vec<topods::ShapeRef>,
    ) {
        use topods::{Orientation, ShapeRef};
        let fi = self.faces.len().wrapping_sub(1);
        if fi >= self.faces.len() {
            return; // no face data
        }
        if fi < face_refs.len() && !face_refs[fi].is_null() {
            return; // already emitted
        }
        let (
            edge_indices,
            inner_wire_edges,
            _tris,
            _normal,
            surface,
            _uv_domain,
            _centroid,
            _area,
            sample_point,
            internal_wire_edges,
        ) = &self.faces[fi];

        // Edge --TShape::Edge (read curve/range from custom arrays)
        // Vertex identity is handled by BRep::add_tvertex (vert_by_pos cache).
        let mut e_map: Vec<ShapeRef> = Vec::with_capacity(edge_indices.len());
        for &(ei, _forward) in edge_indices.iter() {
            if ei >= self.edges.len() {
                continue;
            }
            let (v1, v2) = self.edges[ei];
            let first = if v1 < self.vertices.len() {
                t.add_tvertex(self.vertices[v1])
            } else {
                ShapeRef::NULL
            };
            let last = if v2 < self.vertices.len() {
                t.add_tvertex(self.vertices[v2])
            } else {
                ShapeRef::NULL
            };
            let curve = self
                .custom_edge_curves
                .get(ei)
                .and_then(|c| c.as_ref())
                .cloned();
            let curve_range = self
                .custom_edge_ranges
                .get(ei)
                .and_then(|r| *r)
                .or_else(|| {
                    self.custom_edge_curves
                        .get(ei)
                        .and_then(|c| c.as_ref())
                        .map(|crv| {
                            use rcad_kernel::geom::CurveEval;
                            crv.default_domain()
                        })
                })
                .unwrap_or([0.0, 1.0]);
            let e_ref = t.add_tedge(curve, first, last, curve_range);
            while e_map.len() <= ei {
                e_map.push(ShapeRef::NULL);
            }
            e_map[ei] = e_ref;
        }

        // Outer wire
        let outer_edges: Vec<ShapeRef> = edge_indices
            .iter()
            .filter_map(|&(idx, forward)| {
                let orient = if forward {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                if idx < e_map.len() {
                    Some(ShapeRef::synthetic_with_orientation(
                        e_map[idx].index,
                        orient,
                    ))
                } else {
                    None
                }
            })
            .collect();
        let outer_wire = t.add_twire(outer_edges);
        t.wire_mut(outer_wire).flags |= rcad_kernel::topods::tshape_flags::CLOSED;

        // Inner wires
        let mut inner_wires = Vec::new();
        for wire_idxs in inner_wire_edges {
            let iw_edges: Vec<ShapeRef> = wire_idxs
                .iter()
                .filter_map(|&(idx, forward)| {
                    let orient = if forward {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                    if idx < e_map.len() {
                        Some(ShapeRef::synthetic_with_orientation(
                            e_map[idx].index,
                            orient,
                        ))
                    } else {
                        None
                    }
                })
                .collect();
            if !iw_edges.is_empty() {
                let w = t.add_twire(iw_edges);
                t.wire_mut(w).flags |= rcad_kernel::topods::tshape_flags::CLOSED;
                inner_wires.push(w);
            }
        }

        // Internal wire edges
        for iw_edges in internal_wire_edges {
            let iw: Vec<ShapeRef> = iw_edges
                .iter()
                .filter_map(|&(idx, forward)| {
                    let orient = if forward {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                    if idx < e_map.len() {
                        Some(ShapeRef::synthetic_with_orientation(
                            e_map[idx].index,
                            orient,
                        ))
                    } else {
                        None
                    }
                })
                .collect();
            if iw.len() >= 2 {
                let w = t.add_twire(iw);
                t.wire_mut(w).flags |= rcad_kernel::topods::tshape_flags::CLOSED;
                inner_wires.push(w);
            }
        }

        // Face
        // surface passed directly to add_tface (geometry on TShape)
        let internal_vtx: Vec<ShapeRef> = Vec::new();
        let nr = self
            .face_natural_restriction
            .get(fi)
            .copied()
            .unwrap_or(false);
        let face_sr = t.add_tface(
            Some(surface.clone()),
            outer_wire,
            inner_wires,
            Some(*sample_point),
            *_uv_domain,
            internal_vtx,
            nr,
        );
        face_refs.push(face_sr);
    }

    /// --BuildResult(FACE) --create topods vertices, edges, wires, faces
    /// in t_brep from the flat arrays.  Called after fill_images_faces so that
    /// split faces have already been emitted as TShapes via emit_face_topods.
    pub(crate) fn build_topods_faces(
        &mut self,
        t: &mut topods::BRep,
        wire_refs: &[topods::ShapeRef],
        face_refs: &mut Vec<topods::ShapeRef>,
    ) {
        use topods::{Orientation, ShapeRef};

        // Architecture A1: skip faces already emitted as TShapes (via emit_face_topods).
        // face_refs already contains ShapeRefs for incrementally-emitted faces.
        let start_fi = face_refs.len();
        if start_fi >= self.faces.len() {
            return; // all faces already have TShapes
        }

        // 1. Vertices --TShape::Vertex (identity-based dedup via BRep::add_tvertex).
        let n_verts = self.vertices.len();
        let mut vi_to_ti: Vec<usize> = Vec::with_capacity(n_verts);
        for v in &self.vertices {
            let sr = t.add_tvertex(*v);
            vi_to_ti.push(sr.index);
        }

        // 2. Edges --TShape::Edge (use vi_to_ti to map vertex indices).
        let mut e_map: Vec<ShapeRef> = Vec::with_capacity(self.edges.len());
        for (ei, &(start, end)) in self.edges.iter().enumerate() {
            let first = ShapeRef::synthetic(vi_to_ti[start]);
            let last = ShapeRef::synthetic(vi_to_ti[end]);
            let curve = self
                .custom_edge_curves
                .get(ei)
                .and_then(|c| c.as_ref())
                .cloned();
            let curve_range = self
                .custom_edge_ranges
                .get(ei)
                .and_then(|r| *r)
                .or_else(|| {
                    self.custom_edge_curves
                        .get(ei)
                        .and_then(|c| c.as_ref())
                        .map(|crv| {
                            use rcad_kernel::geom::CurveEval;
                            crv.default_domain()
                        })
                })
                .unwrap_or([0.0, 1.0]);
            e_map.push(t.add_tedge(curve, first, last, curve_range));
        }

        // 3. Faces --TShape::Face (with wires) --only for faces NOT yet in face_refs.
        for (
            flat_fi,
            (
                edge_indices,
                inner_wire_edges,
                _triangles,
                _normal,
                surface,
                _uv_domain,
                _centroid,
                _area,
                sample_point,
                internal_wire_edges,
            ),
        ) in self.faces.iter().enumerate().skip(start_fi)
        {
            // use pre-built wires from wire_refs when available.
            let wire_idxs: Option<&Vec<usize>> = self.face_all_wire_idxs.get(flat_fi);

            // Outer wire --use pre-built if available
            let outer_wire = wire_idxs
                .and_then(|idxs| idxs.first().copied())
                .and_then(|wi| wire_refs.get(wi).filter(|sr| !sr.is_null()).copied())
                .unwrap_or_else(|| {
                    // Fallback: create wire inline from edge_indices
                    let outer_edges: Vec<ShapeRef> = edge_indices
                        .iter()
                        .filter_map(|&(idx, forward)| {
                            let orient = if forward {
                                Orientation::Forward
                            } else {
                                Orientation::Reversed
                            };
                            if idx < e_map.len() {
                                Some(ShapeRef::synthetic_with_orientation(
                                    e_map[idx].index,
                                    orient,
                                ))
                            } else {
                                None
                            }
                        })
                        .collect();
                    t.add_twire(outer_edges)
                });

            // Inner wires
            let mut inner_wires = Vec::new();

            // Pre-built inner wires from wire_refs
            if let Some(ref idxs) = wire_idxs {
                for &wi in idxs.iter().skip(1) {
                    if let Some(&sr) = wire_refs.get(wi) {
                        if !sr.is_null() {
                            inner_wires.push(sr);
                        }
                    }
                }
            }

            // Fallback: create inner wires inline for any extra edges not covered by wire_refs
            let covered_count = inner_wires.len();
            for wire_idxs in inner_wire_edges.iter().skip(covered_count) {
                let iw_edges: Vec<ShapeRef> = wire_idxs
                    .iter()
                    .filter_map(|&(idx, forward)| {
                        let orient = if forward {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                        if idx < e_map.len() {
                            Some(ShapeRef::synthetic_with_orientation(
                                e_map[idx].index,
                                orient,
                            ))
                        } else {
                            None
                        }
                    })
                    .collect();
                if !iw_edges.is_empty() {
                    let w = t.add_twire(iw_edges);
                    t.wire_mut(w).flags |= rcad_kernel::topods::tshape_flags::CLOSED;
                    inner_wires.push(w);
                }
            }
            // Internal wire edges
            for iw_edges in internal_wire_edges {
                let iw: Vec<ShapeRef> = iw_edges
                    .iter()
                    .filter_map(|&(idx, forward)| {
                        let orient = if forward {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                        if idx < e_map.len() {
                            Some(ShapeRef::synthetic_with_orientation(
                                e_map[idx].index,
                                orient,
                            ))
                        } else {
                            None
                        }
                    })
                    .collect();
                if iw.len() >= 2 {
                    let w = t.add_twire(iw);
                    t.wire_mut(w).flags |= rcad_kernel::topods::tshape_flags::CLOSED;
                    inner_wires.push(w);
                }
            }

            let internal_vtx: Vec<ShapeRef> =
                self.face_internal_vtx.get(flat_fi).map_or(vec![], |v| {
                    v.iter()
                        .map(|&vi| ShapeRef::synthetic(vi_to_ti.get(vi).copied().unwrap_or(vi)))
                        .collect()
                });
            let nr = self
                .face_natural_restriction
                .get(flat_fi)
                .copied()
                .unwrap_or(true);
            face_refs.push(t.add_tface(
                Some(surface.clone()),
                outer_wire,
                inner_wires,
                Some(*sample_point),
                *_uv_domain,
                internal_vtx,
                nr,
            ));
        }
    }

    /// --Final assembly --return history (PIOperation_FillHistory).
    /// Per-dimension BuildResult calls (Face/Shell/Solid/CompSolid) have already
    /// created the corresponding topods TShapes in t_brep. This method returns
    /// the BooleanHistory from the accumulated result data.
    /// When `fill_history` is false (OCCT: !HasHistory --!myFillHistory),
    /// returns an empty history with no origins tracking.
    pub(crate) fn build_topods(
        &mut self,
        t: &mut topods::BRep,
        fill_history: bool,
        shells: &[topods::ShapeRef],
        face_refs: &mut Vec<topods::ShapeRef>,
        solids: &[topods::ShapeRef],
        compsolid_groups: &[topods::ShapeRef],
    ) -> BooleanHistory {
        // Use tmp_solids when BuildResult did not create solids.
        // tmp_solids entries contain shell indices (into tmp_shells).
        // OCCT L1203 fallback: BuildResult(SOLID) already created tshapes when
        // BuildRC/ComputeState ran, but if skipped (no images), create from scratch.
        if shells.is_empty() && !self.tmp_solids.is_empty() {
            for tmp_solid in std::mem::take(&mut self.tmp_solids) {
                let faces: Vec<topods::ShapeRef> = tmp_solid
                    .iter()
                    .filter_map(|&shi| self.tmp_shells.get(shi))
                    .flat_map(|v| v.iter())
                    .filter_map(|&fi| face_refs.get(fi).copied())
                    .filter(|sr| !sr.is_null())
                    .collect();
                if !faces.is_empty() {
                    let shell = t.add_tshell(faces);
                    t.add_tsolid(vec![shell]);
                }
            }
        } else if !face_refs.is_empty() {
            let shell = t.add_tshell(std::mem::take(face_refs));
            t.add_tsolid(vec![shell]);
        }

        BooleanHistory {
            face_origins: if fill_history {
                std::mem::take(&mut self.face_origins)
            } else {
                vec![]
            },
            co_face_origins: if fill_history {
                std::mem::take(&mut self.co_face_origins)
            } else {
                vec![]
            },
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
}
