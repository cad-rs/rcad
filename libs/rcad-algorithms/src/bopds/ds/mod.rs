pub mod face_aabb;
pub mod types;
pub use types::*;
pub mod iterator;
pub use iterator::BOPDS_Iterator;
pub mod topods_builder;

/// Phase 4 transition: empty HashMap for fallback edge_vertex_params.
static EMPTY_VERTEX_PARAMS: LazyLock<HashMap<usize, f64>> = LazyLock::new(HashMap::new);

#[cfg(test)]
mod integration_tests;

use super::pave::{Pave, PaveBlock, SharedPB, NO_EDGE};
use super::common_block::CommonBlock;
use super::face_info::FaceInfo;
use crate::tolerance::*;
use std::collections::HashMap;
use std::sync::LazyLock;
use glam::{DVec2, DVec3};
use rcad_kernel::topods;
use rcad_kernel::PCurve;
use rcad_kernel::{CurveEval, SurfaceEval};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Line2d, Line3, Plane, Surface3, any_perpendicular};

impl DS {
    /// Create an empty DS (equivalent to OCCT BOPDS_DS default constructor).
    pub fn new_empty() -> Self {
        Self {
            shapes: Vec::new(),
            vertices: Vec::new(), edges: Vec::new(), wires: Vec::new(),
            shells: Vec::new(), solids: Vec::new(), comp_solids: Vec::new(),
            faces: Vec::new(),
            vertex_origins: Vec::new(), vertex_is_internal: Vec::new(), vertex_locations: Vec::new(), vertex_shape_idx: Vec::new(),
            edge_start_vertex: Vec::new(), edge_end_vertex: Vec::new(), edge_origins: Vec::new(),
            edge_paves: Vec::new(), edge_pave_blocks: Vec::new(), edge_face_reps: Vec::new(),
            edge_is_internal: Vec::new(), edge_face_tols: Vec::new(), edge_locations: Vec::new(), edge_shape_idx: Vec::new(),
            face_boundary_verts: Vec::new(), face_boundary_edges: Vec::new(),
            face_boundary_forwards: Vec::new(), face_inner_boundary: Vec::new(),
            face_outer_wire_idxs: Vec::new(), face_inner_wire_idxs: Vec::new(),
            face_normals: Vec::new(), face_origins: Vec::new(),
            source_face_idxs: Vec::new(), face_locations: Vec::new(), face_uv_boundary: Vec::new(),
            source_shell_idxs: Vec::new(), source_solid_idxs: Vec::new(), source_compsolid_idxs: Vec::new(), face_shape_idx: Vec::new(),
            wire_shape_idx: Vec::new(), shell_shape_idx: Vec::new(), solid_shape_idx: Vec::new(), compsolid_shape_idx: Vec::new(),
            interf_vv: Vec::new(), interf_ve: Vec::new(), interf_vf: Vec::new(),
            interf_ee: Vec::new(), interf_ef: Vec::new(), interf_ff: Vec::new(),
            interf_vz: Vec::new(), interf_ez: Vec::new(), interf_fz: Vec::new(),
            interf_zz: Vec::new(),
            intersection_curves: Vec::new(), ff_points: Vec::new(),
            section_edge_refs: Vec::new(),
            fuzzy_tol: crate::tolerance::TOLERANCE_ABS,
            a_vertex_count: 0, a_edge_count: 0, a_face_count: 0,
            shared_topology: Default::default(),
            shape_sd: ShapeSD::new(0, &SharedTopologyInfo::default()),
            same_domain_overlaps: Vec::new(), common_blocks: Vec::new(),
            my_images: Vec::new(), my_origins: Vec::new(),
            wire_images: Vec::new(), shell_images: Vec::new(),
            solid_images: Vec::new(), locations: Vec::new(),
            pave_blocks: Vec::new(),
            increased_ss: std::collections::HashSet::new(),
            interf_tb: std::collections::HashSet::new(),
            map_ve: std::collections::HashMap::new(),
            shape_info: Vec::new(), nb_source_shapes: 0,
        }
    }

    // ----- Phase 4: Helper methods accessing parallel arrays or ds.shape(n) -----

    /// Mutable access to TVertexData for a vertex (for write operations).
    pub fn vertex_data_mut(&mut self, vi: usize) -> &mut topods::TVertexData {
        assert!(vi < self.vertices.len(),
            "vertex_data_mut: {} >= vertices.len() {}", vi, self.vertices.len());
        let si = if vi < self.vertex_shape_idx.len() { self.vertex_shape_idx[vi] } else { vi };
        match std::sync::Arc::make_mut(&mut self.shapes[si]) {
            topods::TShape::Vertex(vd) => vd,
            _ => panic!("vertex_data_mut: {} (shape[{}]) is not a Vertex", vi, si),
        }
    }

    /// Mutable access to TEdgeData for an edge (for write operations).
    pub fn edge_data_mut(&mut self, ei: usize) -> &mut topods::TEdgeData {
        let si = if ei < self.edge_shape_idx.len() { self.edge_shape_idx[ei] } else { self.vertices.len() + ei };
        match std::sync::Arc::make_mut(&mut self.shapes[si]) {
            topods::TShape::Edge(ed) => ed,
            _ => panic!("edge_data_mut: {} (shape[{}]) is not an Edge", ei, si),
        }
    }

    /// Mutable access to TFaceData for a face (for write operations).
    pub fn face_data_mut(&mut self, fi: usize) -> &mut topods::TFaceData {
        let si = if fi < self.face_shape_idx.len() { self.face_shape_idx[fi] } else { fi };
        match std::sync::Arc::make_mut(&mut self.shapes[si]) {
            topods::TShape::Face(fd) => fd,
            _ => panic!("face_data_mut: {} (shape[{}]) is not a Face", fi, si),
        }
    }

    /// Vertex point from TShape via vertex_shape_idx.
    pub fn vertex_point(&self, vi: usize) -> glam::DVec3 {
        let si = if vi < self.vertex_shape_idx.len() { self.vertex_shape_idx[vi] } else { vi };
        match self.shape(si) {
            topods::TShape::Vertex(d) => d.point,
            _ => glam::DVec3::ZERO,
        }
    }

    /// Vertex tolerance from TShape via vertex_shape_idx.
    pub fn vertex_tolerance(&self, vi: usize) -> f64 {
        let si = if vi < self.vertex_shape_idx.len() { self.vertex_shape_idx[vi] } else { vi };
        match self.shape(si) {
            topods::TShape::Vertex(d) => d.tolerance,
            _ => 0.0,
        }
    }

    /// Vertex origin from parallel array (intersection-only data).
    pub fn vertex_origin(&self, vi: usize) -> Option<ShapeOrigin> {
        if vi < self.vertex_origins.len() { self.vertex_origins[vi] }
        else { self.vertices.get(vi).and_then(|v| v.origin) }
    }

    /// Vertex is_internal from parallel array.
    pub fn vertex_is_internal(&self, vi: usize) -> bool {
        if vi < self.vertex_is_internal.len() { self.vertex_is_internal[vi] }
        else { self.vertices.get(vi).map_or(false, |v| v.is_internal) }
    }

    /// Vertex location from parallel array.
    pub fn vertex_location(&self, vi: usize) -> u32 {
        if vi < self.vertex_locations.len() { self.vertex_locations[vi] }
        else { self.vertices.get(vi).map_or(0, |v| v.location) }
    }

    /// Edge curve from TShape.
    pub fn edge_curve(&self, ei: usize) -> Option<&rcad_kernel::geom::Curve3> {
        let si = if ei < self.edge_shape_idx.len() { self.edge_shape_idx[ei] } else { self.vertices.len() + ei };
        match self.shape(si) {
            topods::TShape::Edge(d) => d.curve.as_ref(),
            _ => None,
        }
    }

    /// Edge tolerance from TShape.
    pub fn edge_tolerance(&self, ei: usize) -> f64 {
        let si = if ei < self.edge_shape_idx.len() { self.edge_shape_idx[ei] } else { self.vertices.len() + ei };
        match self.shape(si) {
            topods::TShape::Edge(d) => d.tolerance,
            _ => 0.0,
        }
    }

    /// Edge parametric range from TShape.
    pub fn edge_range(&self, ei: usize) -> [f64; 2] {
        let si = if ei < self.edge_shape_idx.len() { self.edge_shape_idx[ei] } else { self.vertices.len() + ei };
        match self.shape(si) {
            topods::TShape::Edge(d) => d.range,
            _ => [0.0, 0.0],
        }
    }

    /// Edge start vertex DS index from parallel array.
    pub fn edge_start_vertex_ds(&self, ei: usize) -> usize {
        if ei < self.edge_start_vertex.len() { self.edge_start_vertex[ei] }
        else { self.edges.get(ei).map_or(0, |e| e.start_vertex) }
    }

    /// Edge end vertex DS index from parallel array.
    pub fn edge_end_vertex_ds(&self, ei: usize) -> usize {
        if ei < self.edge_end_vertex.len() { self.edge_end_vertex[ei] }
        else { self.edges.get(ei).map_or(0, |e| e.end_vertex) }
    }

    /// Edge origin from parallel array.
    pub fn edge_origin(&self, ei: usize) -> ShapeOrigin {
        if ei < self.edge_origins.len() { self.edge_origins[ei] }
        else { self.edges.get(ei).map_or(ShapeOrigin::ShapeA, |e| e.origin) }
    }

    /// Edge paves from parallel array.
    pub fn edge_paves(&self, ei: usize) -> &[crate::bopds::pave::Pave] {
        if ei < self.edge_paves.len() { &self.edge_paves[ei] }
        else if let Some(e) = self.edges.get(ei) { &e.paves }
        else { &[] }
    }

    /// Edge pave_blocks.
    /// OCCT: myDS->ChangePaveBlocks(nE) — single path.
    /// rcad: prefer self.edges[ei].pave_blocks over parallel array.
    pub fn edge_pave_blocks(&self, ei: usize) -> &[crate::bopds::pave::SharedPB] {
        if let Some(e) = self.edges.get(ei) { &e.pave_blocks }
        else if ei < self.edge_pave_blocks.len() { &self.edge_pave_blocks[ei] }
        else { &[] }
    }

    /// OCCT: myDS->ChangePaveBlocks(nE) — lazy init via HasReference check.
    /// OCCT BOPDS_DS.cxx L425-431: if (!HasReference()) InitPaveBlocks(theIndex).
    pub fn edge_pave_blocks_mut(&mut self, ei: usize) -> &mut Vec<crate::bopds::pave::SharedPB> {
        if ei < self.edges.len() {
            let si = *self.edge_shape_idx.get(ei).unwrap_or(&usize::MAX);
            let has_ref = si < self.shape_info.len() && self.shape_info[si].has_reference();
            if !has_ref {
                self.init_pave_blocks_for_edge(ei);
            }
            &mut self.edges[ei].pave_blocks
        } else if ei < self.edge_pave_blocks.len() {
            &mut self.edge_pave_blocks[ei]
        } else {
            panic!("edge_pave_blocks_mut: index {} out of bounds", ei);
        }
    }

    /// Edge face_reps from parallel array.
    pub fn edge_face_reps(&self, ei: usize) -> &[DSCurveRepOnFace] {
        if ei < self.edge_face_reps.len() { &self.edge_face_reps[ei] }
        else if let Some(e) = self.edges.get(ei) { &e.face_reps }
        else { &[] }
    }

    /// Edge is_internal from parallel array.
    pub fn edge_is_internal(&self, ei: usize) -> bool {
        if ei < self.edge_is_internal.len() { self.edge_is_internal[ei] }
        else { self.edges.get(ei).map_or(false, |e| e.is_internal) }
    }

    /// Edge per-face tolerances from parallel array.
    pub fn edge_face_tols(&self, ei: usize) -> &[(usize, f64)] {
        if ei < self.edge_face_tols.len() { &self.edge_face_tols[ei] }
        else if let Some(e) = self.edges.get(ei) { &e.face_tolerances }
        else { &[] }
    }

    /// Edge location from parallel array.
    pub fn edge_location(&self, ei: usize) -> u32 {
        if ei < self.edge_locations.len() { self.edge_locations[ei] }
        else { self.edges.get(ei).map_or(0, |e| e.location) }
    }

    /// Edge vertex_params from TShape.
    pub fn edge_vertex_params(&self, ei: usize) -> &std::collections::HashMap<usize, f64> {
        let si = if ei < self.edge_shape_idx.len() { self.edge_shape_idx[ei] } else { self.vertices.len() + ei };
        match self.shape(si) {
            topods::TShape::Edge(d) => &d.vertex_params,
            _ => {
                // Fallback to old DSEdge during Phase 4 transition
                self.edges.get(ei).map_or(&*EMPTY_VERTEX_PARAMS, |e| &e.vertex_params)
            }
        }
    }

    /// Edge is_geometric from TShape (curve.is_some()).
    pub fn edge_is_geometric(&self, ei: usize) -> bool {
        let si = if ei < self.edge_shape_idx.len() { self.edge_shape_idx[ei] } else { self.vertices.len() + ei };
        match self.shape(si) {
            topods::TShape::Edge(d) => d.curve.is_some(),
            _ => false,
        }
    }

    /// Face surface from TShape via face_shape_idx.
    pub fn face_surface(&self, fi: usize) -> Option<&rcad_kernel::geom::Surface3> {
        let si = if fi < self.face_shape_idx.len() { self.face_shape_idx[fi] } else { fi };
        match self.shape(si) {
            topods::TShape::Face(d) => d.surface.as_ref(),
            _ => None,
        }
    }

    /// Face tolerance from TShape via face_shape_idx.
    pub fn face_tolerance(&self, fi: usize) -> f64 {
        let si = if fi < self.face_shape_idx.len() { self.face_shape_idx[fi] } else { fi };
        match self.shape(si) {
            topods::TShape::Face(d) => d.tolerance,
            _ => 0.0,
        }
    }

    /// Face natural_restriction from TShape.
    pub fn face_natural_restriction(&self, fi: usize) -> bool {
        match self.shape(fi) {
            topods::TShape::Face(d) => d.natural_restriction,
            _ => false,
        }
    }

    /// Face origin from parallel array.
    pub fn face_origin(&self, fi: usize) -> ShapeOrigin {
        if fi < self.face_origins.len() { self.face_origins[fi] }
        else { self.faces.get(fi).map_or(ShapeOrigin::ShapeA, |f| f.origin) }
    }

    /// FaceInfo from faces array (OCCT BOPDS_ShapeInfo equivalent).
    pub fn face_info(&self, fi: usize) -> &FaceInfo {
        &self.faces[fi].face_info
    }

    /// Mutable FaceInfo — OCCT: ChangeFaceInfo → InitFaceInfo → UpdateFaceInfoOn.
    /// lazy-init: on first access, add boundary vertices to vertices_on.
    pub fn face_info_mut(&mut self, fi: usize) -> &mut FaceInfo {
        if fi < self.faces.len() && self.faces[fi].face_info.vertices_on.is_empty() {
            let bverts: Vec<usize> = self.faces[fi].boundary_verts.clone();
            if !bverts.is_empty() {
                for &vi in &bverts {
                    self.faces[fi].face_info.vertices_on.insert(vi);
                }
            }
        }
        &mut self.faces[fi].face_info
    }

    /// Face boundary edges from parallel array.
    pub fn face_boundary_edges(&self, fi: usize) -> &[usize] {
        if fi < self.face_boundary_edges.len() { &self.face_boundary_edges[fi] }
        else if let Some(f) = self.faces.get(fi) { &f.boundary_edges }
        else { &[] }
    }

    /// Face boundary verts from parallel array.
    pub fn face_boundary_verts(&self, fi: usize) -> &[usize] {
        if fi < self.face_boundary_verts.len() { &self.face_boundary_verts[fi] }
        else if let Some(f) = self.faces.get(fi) { &f.boundary_verts }
        else { &[] }
    }

    /// Face inner boundary edges from parallel array.
    pub fn face_inner_boundary(&self, fi: usize) -> &[Vec<(usize, bool)>] {
        if fi < self.face_inner_boundary.len() { &self.face_inner_boundary[fi] }
        else if let Some(f) = self.faces.get(fi) { &f.inner_boundary_edges }
        else { &[] }
    }

    /// Face outer wire index from parallel array.
    pub fn face_outer_wire_idx(&self, fi: usize) -> Option<usize> {
        if fi < self.face_outer_wire_idxs.len() { self.face_outer_wire_idxs[fi] }
        else { self.faces.get(fi).and_then(|f| f.outer_wire_idx) }
    }

    /// Face inner wire indices from parallel array.
    pub fn face_inner_wire_idxs(&self, fi: usize) -> &[usize] {
        if fi < self.face_inner_wire_idxs.len() { &self.face_inner_wire_idxs[fi] }
        else if let Some(f) = self.faces.get(fi) { &f.inner_wire_idxs }
        else { &[] }
    }

    /// Face normal from parallel array.
    pub fn face_normal(&self, fi: usize) -> glam::DVec3 {
        if fi < self.face_normals.len() { self.face_normals[fi] }
        else { self.faces.get(fi).map_or(glam::DVec3::Z, |f| f.normal) }
    }

    /// Face location from parallel array.
    pub fn face_location(&self, fi: usize) -> u32 {
        if fi < self.face_locations.len() { self.face_locations[fi] }
        else { self.faces.get(fi).map_or(0, |f| f.location) }
    }

    /// Face uv_boundary from parallel array.
    pub fn face_uv_boundary(&self, fi: usize) -> Option<&[glam::DVec2]> {
        if fi < self.face_uv_boundary.len() {
            self.face_uv_boundary[fi].as_ref().map(|v| v.as_slice())
        } else {
            self.faces.get(fi).and_then(|f| f.uv_boundary.as_ref()).map(|v| v.as_slice())
        }
    }

    /// Source face index from parallel array.
    pub fn source_face_idx(&self, fi: usize) -> usize {
        if fi < self.source_face_idxs.len() { self.source_face_idxs[fi] }
        else { self.faces.get(fi).map_or(0, |f| f.source_face_idx) }
    }

    /// Source shell index from parallel array.
    pub fn source_shell_idx(&self, fi: usize) -> Option<usize> {
        if fi < self.source_shell_idxs.len() { self.source_shell_idxs[fi] }
        else { self.faces.get(fi).and_then(|f| f.source_shell_idx) }
    }

    /// Source solid index from parallel array.
    pub fn source_solid_idx(&self, fi: usize) -> Option<usize> {
        if fi < self.source_solid_idxs.len() { self.source_solid_idxs[fi] }
        else { self.faces.get(fi).and_then(|f| f.source_solid_idx) }
    }

    /// Source compsolid index from parallel array.
    pub fn source_compsolid_idx(&self, fi: usize) -> Option<usize> {
        if fi < self.source_compsolid_idxs.len() { self.source_compsolid_idxs[fi] }
        else { self.faces.get(fi).and_then(|f| f.source_compsolid_idx) }
    }

    /// myDS->Shape(n). Returns &TShape at the given flat index.
    pub fn shape(&self, idx: usize) -> &topods::TShape {
        &self.shapes[idx]
    }

    /// mutable shape access for PaveFiller tolerance updates.
    pub fn shape_mut(&mut self, idx: usize) -> &mut topods::TShape {
        std::sync::Arc::make_mut(&mut self.shapes[idx])
    }

    /// myDS->Append 閳?add a new TShape and return its index.
    pub fn append_shape(&mut self, ts: topods::TShape) -> usize {
        let idx = self.shapes.len();
        self.shapes.push(std::sync::Arc::new(ts));
        idx
    }

    /// BOPDS_DS::NbShapes 閳?total count of all shapes.
    pub fn nb_shapes(&self) -> usize {
        self.shapes.len()
    }

    /// myDS->ShapeType(n) == TopAbs_VERTEX.
    pub fn is_vertex(&self, idx: usize) -> bool {
        self.shapes.get(idx).map_or(false, |s| matches!(&**s, topods::TShape::Vertex(_)))
    }

    /// myDS->ShapeType(n) == TopAbs_EDGE.
    pub fn is_edge(&self, idx: usize) -> bool {
        self.shapes.get(idx).map_or(false, |s| matches!(&**s, topods::TShape::Edge(_)))
    }

    /// myDS->ShapeType(n) == TopAbs_FACE.
    pub fn is_face(&self, idx: usize) -> bool {
        self.shapes.get(idx).map_or(false, |s| matches!(&**s, topods::TShape::Face(_)))
    }

    /// push a vertex (myDS->Append) + track in flat array.
    /// When `tshape` is Some, uses that TShape for ds.shapes instead of creating a synthetic one.
    /// Returns the vertex index.
    pub fn push_vertex(&mut self, dv: DSVertex, tshape: Option<std::sync::Arc<topods::TShape>>) -> usize {
        let vi = self.vertices.len();
        let origin = dv.origin;
        let is_internal = dv.is_internal;
        let location = dv.location;
        let point = dv.point;
        let geom_tol = dv.geom_tol;
        self.vertices.push(dv);
        // Phase 4: populate parallel arrays
        self.vertex_origins.push(origin);
        self.vertex_is_internal.push(is_internal);
        self.vertex_locations.push(location);
        use topods::tshape_flags;
        self.shapes.push(tshape.unwrap_or_else(|| std::sync::Arc::new(topods::TShape::Vertex(
            topods::TVertexData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                point,
                tolerance: geom_tol,
                points: Vec::new(),
            },
        ))));
        self.vertex_shape_idx.push(self.shapes.len() - 1);
        // OCCT BOPDS_DS::Append: push ShapeInfo for this vertex (keep 1:1 with shapes[]).
        self.shape_info.push(types::ShapeInfo {
            shape_type: rcad_kernel::topods::ShapeType::Vertex,
            sub_shapes: Vec::new(),
            flag: -1,
            reference: -1,
            has_brep: true,
            box_min: Some(point),
            box_max: Some(point),
            box_gap: geom_tol + self.fuzzy_tol * 0.5,
            is_new: true,
            rank: 0,
            source_idx: usize::MAX,
        });
        vi
    }

    /// push an edge (myDS->Append) + track in flat array.
    /// When `tshape` is Some, uses that TShape for ds.shapes instead of creating a synthetic one.
    /// Returns the edge index.
    pub fn push_edge(&mut self, de: DSEdge, tshape: Option<std::sync::Arc<topods::TShape>>) -> usize {
        let ei = self.edges.len();
        let origin = de.origin;
        let is_internal = de.is_internal;
        let e_tol = de.geom_tol;
        let location = de.location;
        let start_vertex = de.start_vertex;
        let end_vertex = de.end_vertex;
        let curve = de.curve.clone();
        let t_range = de.t_range;
        let face_reps = de.face_reps.clone();
        let vertex_params = de.vertex_params.clone();
        let face_tols = de.face_tolerances.clone();
        let paves = de.paves.clone();
        let pave_blocks = de.pave_blocks.clone();
        self.edges.push(de);
        // Phase 4: populate parallel arrays
        self.edge_start_vertex.push(start_vertex);
        self.edge_end_vertex.push(end_vertex);
        self.edge_origins.push(origin);
        self.edge_paves.push(paves);
        self.edge_pave_blocks.push(pave_blocks);
        self.edge_face_reps.push(face_reps);
        self.edge_is_internal.push(is_internal);
        self.edge_face_tols.push(face_tols);
        self.edge_locations.push(location);
        // OCCT-aligned: register new edge in ShapeInfo array
        // (OCCT BOPDS_DS::Append creates ShapeInfo for each new shape).
        let sv_si = *self.vertex_shape_idx.get(start_vertex).unwrap_or(&usize::MAX);
        let ev_si = *self.vertex_shape_idx.get(end_vertex).unwrap_or(&usize::MAX);
        let [t0, t1] = t_range;
        let mut bb_min = curve.point_at(t0);
        let mut bb_max = curve.point_at(t1);
        // OCCT BndLib_Add3dCurve::Add + GeomBndLib_Curve:
        // Per-type optimal AABB computation.
        match &curve {
            Curve3::Line(_) => {
                // OCCT GeomBndLib_Line: endpoints only.
            }
            Curve3::Circle(circ) => {
                // OCCT GeomBndLib_Circle::Box L23-43, L48-113 (1:1)
                // For each axis k (0,1,2):
                //   aXk = XAxis.Direction().Coord(k+1) = x_dir dot axis_k
                //   aYk = YAxis.Direction().Coord(k+1) = y_dir dot axis_k
                // Full circle: aAmp = sqrt(R²*aXk² + R²*aYk²), box = center ± aAmp
                // Arc: per-axis analytical extrema via atan(aYk/aXk)
                use std::f64::consts::PI;
                let aR = circ.radius;
                let aO = circ.center;
                let aXd = circ.x_dir;
                let aYd = circ.y_dir;
                let t_range_span = (t1 - t0).abs();
                let a_period = 2.0 * PI - 1e-12; // Precision::PConfusion approximation
                if t_range_span >= a_period {
                    // Full circle: OCCT L32-43 — aBox.Update(aMin[0..2], aMax[0..2])
                    // Compute per-axis extrema, then set all 3 dims at once.
                    let axes = [DVec3::X, DVec3::Y, DVec3::Z];
                    let a_min = axes.map(|ax| aO.dot(ax) - (aR * aR * aXd.dot(ax) * aXd.dot(ax) + aR * aR * aYd.dot(ax) * aYd.dot(ax)).sqrt());
                    let a_max = axes.map(|ax| aO.dot(ax) + (aR * aR * aXd.dot(ax) * aXd.dot(ax) + aR * aR * aYd.dot(ax) * aYd.dot(ax)).sqrt());
                    // Replace bb_min/bb_max with analytical extrema, then endpoints below are added after match.
                    bb_min = DVec3::new(a_min[0], a_min[1], a_min[2]);
                    bb_max = DVec3::new(a_max[0], a_max[1], a_max[2]);
                } else {
                    // Arc: OCCT L63-109
                    let a_u1 = t0;
                    let a_u2 = t1;
                    // OCCT L64-65: ElCLib::AdjustPeriodic(0., 2.*M_PI, Epsilon(1.), aU1, aU2)
                    let tau = 2.0 * PI;
                    let a_period_inner = tau - 0.0; // ULast - UFirst = 2π
                    let preci = f64::EPSILON;
                    // OCCT L128-147: AdjustPeriodic
                    let adj_a_u1 = a_u1 - (a_u1 / a_period_inner).floor() * a_period_inner;
                    let adj_a_u1 = if tau - adj_a_u1 < preci { adj_a_u1 - a_period_inner } else { adj_a_u1 };
                    let adj_a_u2 = a_u2 - ((a_u2 - adj_a_u1) / a_period_inner).floor() * a_period_inner;
                    let adj_a_u2 = if adj_a_u2 - adj_a_u1 < preci { adj_a_u2 + a_period_inner } else { adj_a_u2 };
                    // OCCT L95-111: ElCLib::InPeriod helper (theU → [base, base+period))
                    let in_period = |the_u: f64, base: f64| -> f64 {
                        let period = tau;
                        if period < f64::EPSILON { return the_u; }
                        let shifted = the_u + period * ((base - the_u) / period).ceil();
                        if shifted >= base { shifted } else { base }
                    };
                    // Add arc endpoints (OCCT L68-70) — already done above via start_pt/end_pt.
                    let axes = [DVec3::X, DVec3::Y, DVec3::Z];
                    for (k, axis) in axes.iter().enumerate() {
                        let aXk = aXd.dot(*axis);
                        let aYk = aYd.dot(*axis);
                        // OCCT L79-88: extremal parameter for min
                        let a_t_extr_min;
                        if aXk.abs() > 1e-15 { // gp::Resolution()
                            a_t_extr_min = in_period((aYk / aXk).atan(), 0.0);
                        } else {
                            a_t_extr_min = PI / 2.0;
                        }
                        // OCCT L89: aTExtrMax = aTExtrMin <= PI ? +PI : -PI
                        let a_t_extr_max = if a_t_extr_min <= PI {
                            a_t_extr_min + PI
                        } else {
                            a_t_extr_min - PI
                        };
                        // OCCT L91-92: compute coordinate values
                        let a_val_min = aR * a_t_extr_min.cos() * aXk
                            + aR * a_t_extr_min.sin() * aYk
                            + aO.dot(*axis);
                        let a_val_max = aR * a_t_extr_max.cos() * aXk
                            + aR * a_t_extr_max.sin() * aYk
                            + aO.dot(*axis);
                        // OCCT L93-97: swap both values AND extremal params
                        let mut a_t_extr_min_mut = a_t_extr_min;
                        let mut a_t_extr_max_mut = a_t_extr_max;
                        if a_val_min > a_val_max {
                            // std::swap(aValMin, aValMax) — values implicitly swapped via usage
                            // std::swap(aTExtrMin, aTExtrMax)
                            std::mem::swap(&mut a_t_extr_min_mut, &mut a_t_extr_max_mut);
                        }
                        // OCCT L99-108: check if in arc range via InPeriod
                        let t_k = in_period(a_t_extr_min_mut, adj_a_u1);
                        if t_k >= adj_a_u1 && t_k <= adj_a_u2 {
                            let p = curve.point_at(a_t_extr_min_mut);
                            bb_min = bb_min.min(p);
                            bb_max = bb_max.max(p);
                        }
                        let t_k2 = in_period(a_t_extr_max_mut, adj_a_u1);
                        if t_k2 >= adj_a_u1 && t_k2 <= adj_a_u2 {
                            let p = curve.point_at(a_t_extr_max_mut);
                            bb_min = bb_min.min(p);
                            bb_max = bb_max.max(p);
                        }
                    }
                }
            }
            Curve3::Ellipse(_) => {
                // OCCT GeomBndLib_Ellipse: sample at 32 points.
                let ns = 32usize;
                for k in 1..ns {
                    let t = t0 + (t1 - t0) * (k as f64) / (ns as f64);
                    let p = curve.point_at(t);
                    bb_min = bb_min.min(p);
                    bb_max = bb_max.max(p);
                }
            }
            _ => {
                // OCCT GeomBndLib_BSplineCurve/BezierCurve/OtherCurve:
                // sample at 32 points + tolerance.
                let ns = 32usize;
                for k in 1..ns {
                    let t = t0 + (t1 - t0) * (k as f64) / (ns as f64);
                    let p = curve.point_at(t);
                    bb_min = bb_min.min(p);
                    bb_max = bb_max.max(p);
                }
            }
        }
        // OCCT L1678-1685: anEdgeBoundBox.Add(aVertexInfo.Box()) —
        // include endpoint vertex boxes (with tolerance) in edge AABB.
        if let Some(sv_p) = self.vertices.get(start_vertex).map(|v| v.point) {
            let sv_tol = self.vertex_tolerance(start_vertex);
            bb_min = bb_min.min(sv_p - DVec3::splat(sv_tol + self.fuzzy_tol));
            bb_max = bb_max.max(sv_p + DVec3::splat(sv_tol + self.fuzzy_tol));
        }
        if end_vertex != start_vertex {
            if let Some(ev_p) = self.vertices.get(end_vertex).map(|v| v.point) {
                let ev_tol = self.vertex_tolerance(end_vertex);
                bb_min = bb_min.min(ev_p - DVec3::splat(ev_tol + self.fuzzy_tol));
                bb_max = bb_max.max(ev_p + DVec3::splat(ev_tol + self.fuzzy_tol));
            }
        }
        let rank: usize = if origin == types::ShapeOrigin::ShapeB { 1usize } else { 0usize };
        self.shape_info.push(types::ShapeInfo {
            shape_type: rcad_kernel::topods::ShapeType::Edge,
            sub_shapes: vec![sv_si, ev_si],
            flag: -1,
            reference: -1, // OCCT BOPDS_DS: initialized without reference; set by InitPaveBlocks
            has_brep: true,
            box_min: Some(bb_min - DVec3::splat(e_tol)),
            box_max: Some(bb_max + DVec3::splat(e_tol)),
            box_gap: e_tol + self.fuzzy_tol * 0.5,
            is_new: false,
            rank,
            source_idx: ei,
        });
        self.shapes.push(tshape.unwrap_or_else(|| {
            let mut pcurves: std::collections::HashMap<usize, (Curve2d, f64, f64)> =
                std::collections::HashMap::new();
            for rep in &self.edge_face_reps[ei] {
                pcurves.insert(rep.face_idx, (rep.pcurve.clone(), rep.start_param, rep.end_param));
            }
            let first = topods::ShapeRef::synthetic(start_vertex);
            let last = topods::ShapeRef::synthetic(end_vertex);
            use topods::tshape_flags;
            std::sync::Arc::new(topods::TShape::Edge(
                topods::TEdgeData {
                    my_shapes: vec![first, last],
                    flags: tshape_flags::DEFAULT,
                    curve: Some(curve),
                    first,
                    last,
                    range: t_range,
                    degenerated: start_vertex == end_vertex,
                    pcurves,
                    representations: Vec::new(),
                    vertex_params,
                    tolerance: e_tol,
                    same_parameter: true,
                    same_range: true,
                },
            ))
        }));
        self.edge_shape_idx.push(self.shapes.len() - 1);
        ei
    }

    /// Phase 4: push a face, populating both old DSFace and parallel arrays.
    /// When `tshape` is Some, uses that TShape for ds.shapes instead of creating a minimal one.
    /// Returns the face index.
    pub fn push_face(&mut self, df: DSFace, tshape: Option<std::sync::Arc<topods::TShape>>) -> usize {
        let fi = self.faces.len();
        let surface = df.surface.clone();
        let natural_restriction = df.natural_restriction;
        let geom_tol = df.geom_tol;
        let location = df.location;
        let boundary_verts = df.boundary_verts.clone();
        let boundary_edges = df.boundary_edges.clone();
        let boundary_forwards = df.boundary_edge_forwards.clone();
        let inner_boundary = df.inner_boundary_edges.clone();
        let outer_wire_idx = df.outer_wire_idx;
        let inner_wire_idxs = df.inner_wire_idxs.clone();
        let normal = df.normal;
        let origin = df.origin;
        let source_face_idx = df.source_face_idx;
        let uv_boundary = df.uv_boundary.clone();
        let source_shell_idx = df.source_shell_idx;
        let source_solid_idx = df.source_solid_idx;
        let source_compsolid_idx = df.source_compsolid_idx;
        self.faces.push(df);
        // Phase 4: populate parallel arrays
        self.face_boundary_verts.push(boundary_verts);
        self.face_boundary_edges.push(boundary_edges);
        self.face_boundary_forwards.push(boundary_forwards);
        self.face_inner_boundary.push(inner_boundary);
        self.face_outer_wire_idxs.push(outer_wire_idx);
        self.face_inner_wire_idxs.push(inner_wire_idxs);
        self.face_normals.push(normal);
        self.face_origins.push(origin);
        self.source_face_idxs.push(source_face_idx);
        self.face_locations.push(location);
        self.face_uv_boundary.push(uv_boundary);
        self.source_shell_idxs.push(source_shell_idx);
        self.source_solid_idxs.push(source_solid_idx);
        self.source_compsolid_idxs.push(source_compsolid_idx);
        self.shapes.push(tshape.unwrap_or_else(|| std::sync::Arc::new(topods::TShape::Face(
            topods::TFaceData {
                my_shapes: Vec::new(),
                flags: topods::tshape_flags::DEFAULT,
                surface: Some(surface),
                surface_location: location,
                outer_wire: topods::ShapeRef::NULL,
                inner_wires: Vec::new(),
                sample_point: None,
                uv_domain: None,
                internal_vertices: Vec::new(),
                tolerance: geom_tol,
                natural_restriction,
            },
        ))));
        self.face_shape_idx.push(self.shapes.len() - 1);
        fi
    }

    /// Phase 4: push a wire, populating DSWire and pushing TShape to shapes[].
    /// Returns the wire index.
    pub fn push_wire(&mut self, dw: DSWire, tshape: Option<std::sync::Arc<topods::TShape>>) -> usize {
        let wi = self.wires.len();
        let edges = dw.edges.clone();
        self.wires.push(dw);
        self.shapes.push(tshape.unwrap_or_else(|| std::sync::Arc::new(topods::TShape::Wire(
            topods::TWireData {
                my_shapes: edges.iter().map(|&ei| topods::ShapeRef::synthetic(ei)).collect(),
                flags: topods::tshape_flags::DEFAULT,
                edges: edges.iter().map(|&ei| topods::ShapeRef::synthetic(ei)).collect(),
            },
        ))));
        self.wire_shape_idx.push(self.shapes.len() - 1);
        wi
    }

    /// Phase 4: push a shell, populating DSShell and pushing TShape to shapes[].
    /// Returns the shell index.
    pub fn push_shell(&mut self, dsh: DSShell, tshape: Option<std::sync::Arc<topods::TShape>>) -> usize {
        let shi = self.shells.len();
        let faces = dsh.faces.clone();
        self.shells.push(dsh);
        self.shapes.push(tshape.unwrap_or_else(|| std::sync::Arc::new(topods::TShape::Shell(
            topods::TShellData {
                my_shapes: faces.iter().map(|&fi| topods::ShapeRef::synthetic(fi)).collect(),
                flags: topods::tshape_flags::DEFAULT,
                faces: faces.iter().map(|&fi| topods::ShapeRef::synthetic(fi)).collect(),
            },
        ))));
        self.shell_shape_idx.push(self.shapes.len() - 1);
        shi
    }

    /// Phase 4: push a solid, populating DSSolid and pushing TShape to shapes[].
    /// Returns the solid index.
    pub fn push_solid(&mut self, dso: DSSolid, tshape: Option<std::sync::Arc<topods::TShape>>) -> usize {
        let soi = self.solids.len();
        let shells = dso.shells.clone();
        self.solids.push(dso);
        self.shapes.push(tshape.unwrap_or_else(|| std::sync::Arc::new(topods::TShape::Solid(
            topods::TSolidData {
                my_shapes: shells.iter().map(|&shi| topods::ShapeRef::synthetic(shi)).collect(),
                flags: topods::tshape_flags::DEFAULT,
                shells: shells.iter().map(|&shi| topods::ShapeRef::synthetic(shi)).collect(),
                internal_vertices: Vec::new(),
                internal_edges: Vec::new(),
            },
        ))));
        self.solid_shape_idx.push(self.shapes.len() - 1);
        soi
    }

    /// Phase 4: push a compsolid, populating DSCompSolid and pushing TShape to shapes[].
    /// Returns the compsolid index.
    pub fn push_compsolid(&mut self, dcs: DSCompSolid, tshape: Option<std::sync::Arc<topods::TShape>>) -> usize {
        let csi = self.comp_solids.len();
        let solids = dcs.solids.clone();
        self.comp_solids.push(dcs);
        self.shapes.push(tshape.unwrap_or_else(|| std::sync::Arc::new(topods::TShape::CompSolid(
            solids.iter().map(|&si| topods::ShapeRef::synthetic(si)).collect(),
        ))));
        self.compsolid_shape_idx.push(self.shapes.len() - 1);
        csi
    }

    ///  ?BOPDS_ShapeInfo::HasFlag / Flag.
 /// Returns the flag value for an edge index, or 0 if not set.
 /// Flag is stored in shape_info[nv + edge_idx].flag matching OCCT's
 /// per-shape integer flag (BOPDS_ShapeInfo::myFlag).
 pub fn edge_flag(&self, edge_idx: usize) -> i64 {
 let si_idx = if edge_idx < self.edge_shape_idx.len() { self.edge_shape_idx[edge_idx] } else { self.vertices.len() + edge_idx };
 if si_idx < self.shape_info.len() {
 let f = self.shape_info[si_idx].flag;
 if f >= 0 { f } else { 0 }
 } else { 0 }
 }

 ///  ?BOPDS_ShapeInfo::HasFlag(int&)  ?true if flag is set (>= 0).
 pub fn edge_has_flag(&self, edge_idx: usize) -> bool {
 let si_idx = if edge_idx < self.edge_shape_idx.len() { self.edge_shape_idx[edge_idx] } else { self.vertices.len() + edge_idx };
 si_idx < self.shape_info.len() && self.shape_info[si_idx].flag >= 0
 }

 ///  ?BOPDS_ShapeInfo::SetFlag.
 pub fn set_edge_flag(&mut self, edge_idx: usize, flag: i64) {
 let si_idx = if edge_idx < self.edge_shape_idx.len() { self.edge_shape_idx[edge_idx] } else { self.vertices.len() + edge_idx };
 if si_idx < self.shape_info.len() {
 self.shape_info[si_idx].flag = flag;
 }
 }

 ///  ?BRep_Tool::Degenerated(edge) equivalent.
 /// Checks the shape_info flag: a flagged edge is degenerated.
 /// Also falls back to start==end check for edges loaded before
 /// shape_info initialization.
 pub fn is_edge_degenerated(&self, edge_idx: usize) -> bool {
 if edge_idx < self.edges.len() {
  let e = &self.edges[edge_idx];
  if e.start_vertex == e.end_vertex {
   // Closed edge - check curve type. Circle/Ellipse with start==end
   // is a full geometric curve (not degenerate, matching OCCT BRep_Tool::Degenerated).
   return match &e.curve {
    rcad_kernel::geom::Curve3::Circle(_) | rcad_kernel::geom::Curve3::Ellipse(_) => false,
    _ => true,
   };
  }
  return false;
 }
 self.edge_has_flag(edge_idx)
 }

 /// OCCT: BOPDS_DS::IsValidShrunkData (BOPDS_DS.cxx L1547-1578).
 /// Returns true if the PB's shrunk range endpoints are within tolerance
 /// of their corresponding vertex positions.
 pub fn is_valid_shrunk_data(&self, pb: &crate::bopds::pave::PaveBlock) -> bool {
 use crate::tolerance::CONFUSION;
 if !pb.has_shrunk_data() { return false; }
 let (ts1, ts2, _splittable) = pb.shrunk_data();
 let (v1i, v2i) = pb.indices();
 let an_epsilon = self.edge_tolerance(pb.original_edge) * 0.01;
 for &(vi, ts) in &[(v1i, ts1), (v2i, ts2)] {
  let vp = self.vertex_point(vi);
  let a_tol = self.vertex_tolerance(vi) + CONFUSION;
  if let Some(curve) = self.edges.get(pb.original_edge) {
   let pp = curve.curve.point_at(ts);
   if a_tol - vp.distance(pp) > an_epsilon { return false; }
  }
 }
 true
 }

 // ----- OCCT BOPDS_DS data layer methods -----

 ///  ?BOPDS_DS::IsNewShape (L228-233).
 /// Returns true if the shape index was appended during intersection
 /// (not part of the original source shapes).  In rcad, vertices
 /// with `origin: None` are intersection-created (new shapes).
 /// Edges with `origin: None` carry the same semantics.
 pub fn is_new_vertex(&self, vi: usize) -> bool {
 let si = if vi < self.vertex_shape_idx.len() { self.vertex_shape_idx[vi] } else { vi };
 if si < self.shape_info.len() {
 self.shape_info[si].is_new
 } else {
 self.vertices.get(vi).map_or(true, |v| v.origin.is_none())
 }
 }

 ///  ?BOPDS_DS::Rank (L214-226).
 /// Returns the rank (operand index 0=A, 1=B) of a shape.
 /// 0 for shapes from operand A, 1 for operand B.
 pub fn rank(&self, vi: usize) -> usize {
 if vi < self.a_vertex_count { 0 } else { 1 }
 }

 ///  ?BOPDS_DS::Range (L207-212).
 /// Returns the index range [start, end) for shapes of given type.
 /// rcad: returns [0, a_vertex_count) for A, [a_vertex_count, end) for B.
 pub fn range(&self, is_a: bool) -> (usize, usize) {
 if is_a { (0, self.a_vertex_count) }
 else { (self.a_vertex_count, self.vertices.len()) }
 }

 // ----- PaveBlock pool accessors (BOPDS_DS.hxx L156-177) -----

///  BOPDS_DS::HasPaveBlocks (hxx:162-164, cxx L708: ShapeInfo(theIndex).HasReference()).
/// Returns true if the edge with the given index has PaveBlocks reference set.
pub fn has_pave_blocks(&self, edge_idx: usize) -> bool {
    let si = if edge_idx < self.edge_shape_idx.len() { self.edge_shape_idx[edge_idx] } else { edge_idx };
    si < self.shape_info.len() && self.shape_info[si].has_reference()
}

 ///  ?BOPDS_DS::ChangePaveBlocks (hxx:172-174).
 /// Returns a mutable reference to the PaveBlocks list for an edge.
 pub fn change_pave_blocks(&mut self, edge_idx: usize) -> &mut Vec<SharedPB> {
 &mut self.edges[edge_idx].pave_blocks
 }

 ///  ?BOPDS_DS::InitPaveBlocks (cxx L437-501).
 /// Creates the initial PaveBlock for a source edge, covering the full
 /// parametric range [t_range[0], t_range[1]].  For closed edges (seam
 /// edges where start == end), both endpoint paves are added to ext_paves
 /// via AppendExtPave (fence-protected), matching OCCT's closed-edge
 /// InitPaveBlocks for an edge (OCCT BOPDS_DS.cxx L437-500).
 /// Transfers all existing paves from edge.paves to the PB as ext_paves,
 /// then calls update() to sort and create sub-PBs — one per vertex-to-vertex segment.
 pub fn init_pave_blocks_for_edge(&mut self, edge_idx: usize) {
 if edge_idx >= self.edges.len() { return; }
 // OCCT ChangePaveBlocks internal: if pool is empty, initialize
 if !self.edges[edge_idx].pave_blocks.is_empty() { return; }
 let (sv, ev, tr0, tr1) = {
 let e = &self.edges[edge_idx];
 (e.start_vertex, e.end_vertex, e.t_range[0], e.t_range[1])
 };
 // OCCT L457-467: create PaveBlock and set pave1/pave2/original_edge
 let pv1 = Pave { vertex_idx: sv, param: tr0 };
 let pv2 = Pave { vertex_idx: ev, param: tr1 };
 let mut pb = PaveBlock::new(edge_idx, pv1, pv2);
 // OCCT: BOPDS_PaveBlock::SetShrunkData sets isSplittable=true for initial edge PBs.
 pb.is_splittable = true;
 // OCCT: shrunk data computed later by FillShrunkData (PaveFiller_9.cxx).
 // Don't set it here — let fill_shrunk_data compute the correct range.
 // OCCT: initial edge PaveBlock stores start/end paves in edge's internal list.
 self.edges[edge_idx].paves.push(pv1);
 self.edges[edge_idx].paves.push(pv2);
 if edge_idx < self.edge_paves.len() {
   self.edge_paves[edge_idx].push(pv1);
   self.edge_paves[edge_idx].push(pv2);
 }
 // OCCT L469+: add ALL existing paves (from edge.paves) as ext_paves
 // excluding endpoint paves (they are pave1/pave2 already).
 let existing_paves: Vec<Pave> = self.edges[edge_idx].paves.clone();
 for p in &existing_paves {
  if p.vertex_idx != sv && p.vertex_idx != ev {
   pb.append_ext_pave(*p);
  }
 }
 // OCCT L477-483: closed edges — add both endpoint paves as ext_paves
 if sv == ev {
  pb.append_ext_pave(Pave { vertex_idx: sv, param: tr0 });
  pb.append_ext_pave(Pave { vertex_idx: sv, param: tr1 });
 }
 // OCCT L499: aPaveBlock->Update(...) — sort ext_paves and create sub-PBs
 // For an initial PB with no ext_paves (empty ext_paves, non-closed edge),
 // keep the single PB so FillShrunkData can compute proper shrunk data.
 let sub_pbs = if pb.ext_paves.is_empty() && sv != ev {
  Vec::new()
 } else {
  pb.update(true)
 };
 // Replace the single PB with sub-PBs. If update returned empty (degenerate),
 // keep the original single PB.
 if sub_pbs.is_empty() {
  self.edges[edge_idx].pave_blocks = vec![SharedPB::new(pb)];
 } else {
  self.edges[edge_idx].pave_blocks = sub_pbs.into_iter().map(SharedPB::new).collect();
 }
 // OCCT L500: anEdgeInfo.SetReference(...)
 let si_idx = if edge_idx < self.edge_shape_idx.len() { self.edge_shape_idx[edge_idx] } else { self.vertices.len() + edge_idx };
 if si_idx < self.shape_info.len() {
  self.shape_info[si_idx].reference = edge_idx as i64;
 }
 }

 ///  ?BOPDS_DS::PaveBlocks (hxx:167-169).
 /// Returns a reference to the PaveBlocks list for an edge.
 pub fn pave_blocks(&self, edge_idx: usize) -> &[SharedPB] {
 &self.edges[edge_idx].pave_blocks
 }

 // ----- OCCT HasInterf / HasSubShape equivalents -----

 ///  ?HasSubShape(nV, nE) =check if vertex is a sub-shape of edge.
 /// Returns true when vertex nV is an endpoint of edge nE.
 pub fn edge_has_vertex(&self, nV: usize, nE: usize) -> bool {
 self.edges.get(nE).map_or(false, |e| e.start_vertex == nV || e.end_vertex == nV)
 }

 ///  ?myDS->HasInterf(nV, nE) =checks VE interference exists.
 pub fn has_interf_ve(&self, vi: usize, ei: usize) -> bool {
 self.interf_ve.iter().any(|inf| inf.vertex == vi && inf.edge == ei)
 }

 /// myDS->HasInterf(n1, n2) =checks VV interference exists.
 pub fn has_interf_vv(&self, v1: usize, v2: usize) -> bool {
 self.interf_vv.iter().any(|inf| (inf.v1 == v1 && inf.v2 == v2) || (inf.v1 == v2 && inf.v2 == v1))
 }

 /// AddShapeSD =register dynamic SD mapping between two vertices.
 pub fn add_shape_sd(&mut self, from: usize, to: usize) {
 self.shape_sd.add_sd_vertex(from, to);
 }

  /// HasShapeSD(n, nSD) =find the SD root vertex.
  pub fn has_shape_sd(&self, v: usize) -> Option<usize> {
    self.shape_sd.find_sd_partner(v)
  }

  /// OCCT BOPDS_DS.cxx L1487-1499: InitPaveBlocksForVertex
  ///
  /// Ensures the PaveBlocks pool is initialized for all edges incident to
  /// the given vertex.  OCCT uses myMapVE (vertex-to-edge map) for fast
  /// lookup.  rcad: uses map_ve (built by build_map_ve).
  /// The per-edge pave_blocks Vec is always allocated at creation time in
  /// rcad (OCCT pool is lazy-initialized via ChangePaveBlocks), so the
  /// body primarily records the  call structure.
  pub fn init_pave_blocks_for_vertex(&mut self, vertex_idx: usize) {
    // OCCT L1489: anEdgeIndices = myMapVE.Seek(theVertexIndex)
    let incident_edges = self.map_ve.get(&vertex_idx).cloned().unwrap_or_default();
    // OCCT L1490-1492: if !anEdgeIndices return
    if incident_edges.is_empty() {
      return;
    }
    // OCCT L1495-1498: for each edge, ChangePaveBlocks(anEdgeIndex)
    //   ensures edge has initial PaveBlock (pool entry).
    for &ei in &incident_edges {
      self.init_pave_blocks_for_edge(ei);
    }
  }

  /// 閿?build vertex-to-edge map (BOPDS_DS::myMapVE).
  ///   Populates map_ve from edges array.  Must be called after init_shape_topo
  ///   loads all source shapes.
  pub fn build_map_ve(&mut self) {
  self.map_ve.clear();
  for (ei, e) in self.edges.iter().enumerate() {
  self.map_ve.entry(e.start_vertex).or_default().push(ei);
  if e.start_vertex != e.end_vertex {
  self.map_ve.entry(e.end_vertex).or_default().push(ei);
  }
  }
  }

 ///  ?myDS->HasInterf(nE1, nE2) =checks EE interference exists.
 pub fn has_interf_ee(&self, e1: usize, e2: usize) -> bool {
 self.interf_ee.iter().any(|inf| (inf.e1 == e1 && inf.e2 == e2) || (inf.e1 == e2 && inf.e2 == e1))
 }

 ///  ?myDS->HasInterf(nV, nF) =checks VF interference exists.
 pub fn has_interf_vf(&self, vi: usize, fi: usize) -> bool {
 self.interf_vf.iter().any(|inf| inf.vertex == vi && inf.face == fi)
 }

 ///  ?myDS->HasInterf(nE, nF) =checks EF interference exists.
 pub fn has_interf_ef(&self, ei: usize, fi: usize) -> bool {
 self.interf_ef.iter().any(|inf| inf.edge == ei && inf.face == fi)
 }

 ///  ?myDS->HasInterf(nF1, nF2) =checks FF interference exists.
  pub fn has_interf_ff(&self, f1: usize, f2: usize) -> bool {
  let (a, b) = if f1 < f2 { (f1, f2) } else { (f2, f1) };
  self.interf_ff.iter().any(|ff| { let (fa, fb) = if ff.f1 < ff.f2 { (ff.f1, ff.f2) } else { (ff.f2, ff.f1) }; fa == a && fb == b })
  }

  /// 閿?BOPDS_DS::AddInterf (DS.cxx L410-420).
  ///   Global interference pair fence.  Checks (i1, i2) with i1 < i2 against
  ///   myInterfTB.  Returns true if the pair is NEW (first insertion), false
  ///   if it already exists.  Call before adding to any typed interference vec.
  ///   The fence prevents duplicate shape pairs across all interference types.
  pub fn try_add_interf(&mut self, i1: usize, i2: usize) -> bool {
  let key = if i1 < i2 { (i1, i2) } else { (i2, i1) };
  self.interf_tb.insert(key)
  }

 /// dedup FaceFace interferences by (Fmin,Fmax) pair.
 /// Merges curves/points from duplicate entries in interf_ff.
 pub fn dedup_ff_interferences(&mut self) {
 let mut merged: std::collections::HashMap<(usize, usize), (Vec<usize>, Vec<FFPoint>)> = std::collections::HashMap::new();
 for inf in &self.interf_ff {
 let key = if inf.f1 < inf.f2 { (inf.f1, inf.f2) } else { (inf.f2, inf.f1) };
 let entry = merged.entry(key).or_insert((Vec::new(), Vec::new()));
 for &c in &inf.curves { if !entry.0.contains(&c) { entry.0.push(c); } }
 // OCCT: points are BOPDS_Point stored inline; dedup by vertex_index when assigned.
 for p in &inf.points {
   if !entry.1.iter().any(|ep| ep.vertex_index != usize::MAX && ep.vertex_index == p.vertex_index) {
     entry.1.push(p.clone());
   }
 }
 }
 self.interf_ff.clear();
  for ((f1, f2), (curves, points)) in &merged {
  self.interf_ff.push(InterferenceFF { f1: *f1, f2: *f2, curves: curves.clone(), points: points.clone(), tangent_faces: false });
  self.try_add_interf(*f1, *f2);
  }
 }


 /// myDS->HasInterfShapeSubShapes(nV, nE)  ?checks if
 /// vertex already has interference with any sub-shape (face) of the edge.
 pub fn has_interf_ve_via_faces(&self, vi: usize, ei: usize) -> bool {
 // Check if the vertex has VF interference with any face that references this edge
 if self.interf_vf.iter().any(|inf| inf.vertex == vi) {
 if self.faces.iter().any(|f| {
 f.boundary_edges.contains(&ei)
 || f.inner_boundary_edges.iter().any(|iw| iw.iter().any(|&(e, _)| e == ei))
 }) {
 return true;
 }
 }
 // Check if this edge has EF interference with a face that also has VF with the vertex
 if self.interf_ef.iter().any(|inf| inf.edge == ei) {
 if self.interf_vf.iter().any(|inf| inf.vertex == vi) {
 return true;
 }
 }
 false
 }

 /// Build DS from two topods::BRep shapes, bypassing the deprecated `rcad_kernel::BRep`.
///
/// Delegates to [`topods_builder::new_from_topods`].
pub fn new_from_topods(a: &topods::BRep, b: &topods::BRep, fuzzy_tol: f64) -> Self {
 topods_builder::new_from_topods(a, b, fuzzy_tol)
}

/// Compute the characteristic scale of the model from all vertices.
 /// Returns the diagonal of the bounding box, or 1.0 if empty.
 pub fn model_scale(&self) -> f64 {
 use glam::DVec3;
 let mut min_pt = DVec3::splat(f64::INFINITY);
 let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
 let mut has_vertices = false;

 for v in &self.vertices {
 min_pt = min_pt.min(v.point);
 max_pt = max_pt.max(v.point);
 has_vertices = true;
 }

 if !has_vertices {
 return 1.0;
 }

 let diagonal = (max_pt - min_pt).length();
 diagonal.max(TOLERANCE_MODEL_SCALE_MIN)
 }

 /// Compute UV boundary for all curved faces by projecting 3D boundary
 /// points onto the face surface's parameter domain.
 ///
 /// For each boundary edge, we sample `N_SAMPLES` evenly-spaced points along
 /// the edge curve so that the resulting UV polygon is well-defined even when
 /// the wire has very few vertices (e.g. a sphere with only 2 poles).
 pub fn compute_uv_boundaries(&mut self) {
 use std::f64::consts::PI;
 const N_SAMPLES: usize = 8;

 for fi in 0..self.faces.len() {
 let surface = self.faces[fi].surface.clone();

 // For sphere and cylinder, the UV boundary is the full parameter
 // domain rectangle. The topological boundary (seam edge) maps to a
 // degenerate line in UV space and cannot be used as a polygon.
 match &surface {
 Surface3::Sphere(_) => {
 // Sphere param from projection: u = longitude [- ? (atan2 range),
 // v = colatitude [0,  . Use the full domain as UV boundary.
 let uv = vec![
 DVec2::new(-PI, 0.0),
 DVec2::new(PI, 0.0),
 DVec2::new(PI, PI),
 DVec2::new(-PI, PI),
 ];
 self.faces[fi].uv_boundary = Some(uv);
 continue;
 }
 Surface3::Cylinder(cyl) => {
 // Cylinder param: u = azimuth [0, 2  (matches CylindricalSurface::point_at),
 // v = height along axis.  Estimate height range from boundary edge samples.
 let boundary_edges = self.faces[fi].boundary_edges.clone();
 let mut h_min = f64::INFINITY;
 let mut h_max = f64::NEG_INFINITY;
 let axis = cyl.axis.normalize();
 let origin = cyl.origin;
 // Seam edges on a cylinder are lines parallel to the axis.  When a face has
 // multiple seam edges (e.g. a cylinder with explicit front/back seams), the
 // u-range is bounded by them rather than the full [0, 2 .
 let u_ax = any_perpendicular(axis);
 let v_ax = axis.cross(u_ax).normalize();
 let mut seam_u_vals: Vec<f64> = Vec::new();
 for ei in &boundary_edges {
 let edge = &self.edges[*ei];
 let [t0, t1] = edge.t_range;
 // Detect seam edges: lines whose direction is parallel to the cylinder axis
 if let Curve3::Line(line) = &edge.curve {
 let dir = line.direction.normalize();
 if dir.dot(axis).abs() > 1.0 - TOLERANCE_ABS {
 let mid = edge.curve.point_at(0.5 * (t0 + t1));
 let v_comp = (mid - origin).dot(axis);
 let radial = mid - origin - v_comp * axis;
 let u = radial.dot(v_ax).atan2(radial.dot(u_ax));
 seam_u_vals.push(u);
 }
 }
 // v-range sampling (same as before)
 for k in 0..=N_SAMPLES {
 let t = t0 + (t1 - t0) * k as f64 / N_SAMPLES as f64;
 let p = edge.curve.point_at(t);
 let h = (p - origin).dot(axis);
 h_min = h_min.min(h);
 h_max = h_max.max(h);
 }
 }
 if !h_min.is_finite() || !h_max.is_finite() {
 h_min = -1.0;
 h_max = 1.0;
 }
 // Deduplicate seam u-values (same edge may appear fwd+rev in the wire)
 seam_u_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
 seam_u_vals.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_MESH_LEGACY);
 let (u_lo, u_hi) = if seam_u_vals.len() >= 2 {
 (seam_u_vals[0], seam_u_vals[seam_u_vals.len() - 1])
 } else {
 (0.0, 2.0 * PI)
 };
 // Add a small margin to the v-range so that intersection polyline
 // endpoints near v=0 or v=h are not clipped by
 // extend_trim_to_uv_boundary.  A 1% margin matches the cone case
 // below and has worked well in practice.
 let v_range = h_max - h_min;
 let margin = if v_range > TOLERANCE_COORD_SUB {
 v_range * 0.01 + TOLERANCE_COORD_SUB
 } else {
 TOLERANCE_COORD_SUB
 };
 let uv = vec![
 DVec2::new(u_lo, h_min - margin),
 DVec2::new(u_hi, h_min - margin),
 DVec2::new(u_hi, h_max + margin),
 DVec2::new(u_lo, h_max + margin),
 ];
 self.faces[fi].uv_boundary = Some(uv);
 continue;
 }
 Surface3::Cone(cone) => {
 // Cone param: u = azimuth [0, 2 , v = slant distance from the
 // reference circle centered at `cone.apex`.
 // Estimate the full slant range from boundary edge samples so
 // reference-circle cones keep the correct UV window.
 let boundary_edges = self.faces[fi].boundary_edges.clone();
 let mut v_min = f64::INFINITY;
 let mut v_max = f64::NEG_INFINITY;
 let ref_point = cone.apex;
 let axis = cone.axis_dir();
 for ei in &boundary_edges {
 let edge = &self.edges[*ei];
 let [t0, t1] = edge.t_range;
 for k in 0..=N_SAMPLES {
 let t = t0 + (t1 - t0) * k as f64 / N_SAMPLES as f64;
 let p = edge.curve.point_at(t);
 let local = p - ref_point;
 let along = local.dot(axis);
 let slant = cone.slant_from_axial(along);
 v_min = v_min.min(slant);
 v_max = v_max.max(slant);
 }
 }
 if !v_min.is_finite() || !v_max.is_finite() {
 v_min = 0.0;
 v_max = 1.0;
 }
 if (v_max - v_min).abs() < TOLERANCE_COORD_SUB {
 v_min -= 0.5;
 v_max += 0.5;
 }
 let margin = (v_max - v_min) * 0.01 + TOLERANCE_COORD_SUB;
 let uv = vec![
 DVec2::new(0.0, v_min - margin),
 DVec2::new(2.0 * PI, v_min - margin),
 DVec2::new(2.0 * PI, v_max + margin),
 DVec2::new(0.0, v_max + margin),
 ];
 self.faces[fi].uv_boundary = Some(uv);
 continue;
 }
 Surface3::Torus(_) => {
 // Torus param: u = major angle [0, 2 , v = minor angle [0, 2 .
 // Full parameter domain is always the UV boundary.
 let uv = vec![
 DVec2::new(0.0, 0.0),
 DVec2::new(2.0 * PI, 0.0),
 DVec2::new(2.0 * PI, 2.0 * PI),
 DVec2::new(0.0, 2.0 * PI),
 ];
 self.faces[fi].uv_boundary = Some(uv);
 continue;
 }
 _ => {}
 }

 let boundary_edges = self.faces[fi].boundary_edges.clone();

 if boundary_edges.is_empty() {
 continue;
 }

 let mut pts_3d: Vec<DVec3> = Vec::new();
 for ei in &boundary_edges {
 let edge = &self.edges[*ei];
 let [t0, t1] = edge.t_range;
 for k in 0..N_SAMPLES {
 let t = t0 + (t1 - t0) * (k as f64) / (N_SAMPLES as f64);
 pts_3d.push(edge.curve.point_at(t));
 }
 }

 if pts_3d.is_empty() {
 continue;
 }

 let uv_pts: Vec<DVec2> = pts_3d
 .iter()
 .map(|&p| {
 let proj = rcad_kernel::projection::closest_point_on_surface(&surface, p, 16);
 DVec2::new(proj.params.0, proj.params.1)
 })
 .collect();

 // For planar surfaces (Plane and planar BSpline), decimate colinear
 // UV points. OCCT's BOPAlgo_BuilderFace uses the face's topological
 // edges directly (each edge has exactly 2 vertices), not sampled UV
 // polylines.  The 8-sample-per-edge UV boundary creates edge fragments.
 // Tolerance TOLERANCE_MESH_LEGACY * edge_length accounts for projection noise from
 // closest_point_on_surface on BSpline surfaces (Newton iteration).
 let is_planar = matches!(&surface, Surface3::Plane(_))
 || (if let Surface3::BSpline(ref bsp) = surface { rcad_kernel::geom::bspline_is_planar(bsp, TOLERANCE_ABS) } else { false });
 let decimated = if is_planar {
 let n = uv_pts.len();
 if n > 2 {
 let mut kept: Vec<DVec2> = Vec::with_capacity(n);
 for i in 0..n {
 let prev = uv_pts[(i + n - 1) % n];
 let curr = uv_pts[i];
 let next = uv_pts[(i + 1) % n];
 let d1 = curr - prev;
 let d2 = next - curr;
 let cross = (d1.x * d2.y - d1.y * d2.x).abs();
 let len1 = d1.length_squared();
 let len2 = d2.length_squared();
 if cross < TOLERANCE_MESH_LEGACY * len1.max(len2).max(f64::MIN_POSITIVE)
 && len1 > f64::MIN_POSITIVE && len2 > f64::MIN_POSITIVE
 {
 continue;
 }
 kept.push(curr);
 }
 if kept.len() >= 3 && kept.len() < uv_pts.len() {
 kept
 } else { uv_pts }
 } else { uv_pts }
 } else {
 uv_pts
 };

 self.faces[fi].uv_boundary = Some(decimated);
 }
 }

 ///  ?find existing vertex within tolerance (PutPaveOnCurve equivalent).
 /// OCCT's IsVertexOnLine checks if a boundary vertex lies on the intersection
 /// curve, then places the EXISTING vertex index on the curve's pave block,
 /// ensuring the section edge reuses the same TopoDS_Vertex.  This tolerance-
 /// based scan achieves the same sharing for rcad's flat vertex array.
 ///  ?access edge's pcurve representation on a specific face.
 /// Returns None when no representation exists for this (edge, face) pair.
 pub fn edge_on_face(&self, edge_idx: usize, face_idx: usize) -> Option<&DSCurveRepOnFace> {
 self.edges.get(edge_idx)?.face_reps.iter().find(|r| r.face_idx == face_idx)
 }

 ///  ?compute pcurve for a boundary edge on its face surface.
 /// Mirrors BRep_Tool::CurveOnSurface for boundary edges.
 /// Returns (pcurve, pcurve_span_length) where pcurve has normalized direction.
 pub(crate) fn compute_edge_pcurve(curve: &Curve3, surface: &Surface3, pl_basis: Option<(DVec3, DVec3, DVec3)>) -> Option<(Curve2d, f64)> {
 match (surface, curve) {
 (Surface3::Plane(p), Curve3::Line(l)) => {
 let (u_axis, v_axis) = if let Some((u, v, _)) = pl_basis { (u, v) }
  else { let u = DVec3::X - p.normal * p.normal.dot(DVec3::X); let u = if u.length_squared() < 1e-24 { DVec3::Z - p.normal * p.normal.dot(DVec3::Z) } else { u }; (u.normalize(), p.normal.cross(u).normalize()) };
 let diff = l.origin - p.origin;
 let origin = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));
 let dir = DVec2::new(l.direction.dot(u_axis), l.direction.dot(v_axis));
 let len = dir.length();
 if len > TOLERANCE_CLAMP_MIN {
 Some((Curve2d::Line(Line2d { origin, direction: dir / len }), len))
 } else { None }
 }
 (Surface3::Plane(p), Curve3::Circle(c)) => {
 let (u_axis, v_axis) = if let Some((u, v, _)) = pl_basis { (u, v) }
  else { let u = DVec3::X - p.normal * p.normal.dot(DVec3::X); let u = if u.length_squared() < 1e-24 { DVec3::Z - p.normal * p.normal.dot(DVec3::Z) } else { u }; (u.normalize(), p.normal.cross(u).normalize()) };
 let diff = c.center - p.origin;
 let center_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));
 let normal_dot = c.normal.dot(p.normal).abs();
 if (normal_dot - 1.0).abs() < TOLERANCE_MESH_LEGACY {
 let perim = std::f64::consts::TAU * c.radius;
 Some((Curve2d::Circle(rcad_kernel::geom::Circle2d { center: center_2d, x_dir: DVec2::X, y_dir: DVec2::Y, radius: c.radius  }), perim))
 } else { None }
 }
 // compute pcurve for curved surfaces by projecting
 // the edge's 3D curve start/end points onto UV space.
 // OCCT BRep_Tool::CurveOnSurface returns a parametric curve on
 // the face surface for every boundary edge (BRep_CurveRepresentation).
 // For edges that are not seam/deg on periodic surfaces, the simple
 // Line2d approximation from endpoint UV projection is equivalent to
 // OCCT's stored pcurve (the edge is short enough that the pcurve is
 // well-approximated by a line segment in UV space).
 (Surface3::Sphere(s), _) => {
 let mut uv_start = s.world_to_uv(curve.point_at(0.0));
 let mut uv_end = s.world_to_uv(curve.point_at(1.0));
 // At sphere poles (V=0 or V= ?, U is undefined (atan2(0,0) ambiguity).
 // Use the midpoint's U which is reliable (midpoint is at or near equator).
 let at_pole = |v: f64| v.abs() < 1e-10 || (v - std::f64::consts::PI).abs() < 1e-10;
 if at_pole(uv_start.y) && at_pole(uv_end.y) {
 let uv_mid = s.world_to_uv(curve.point_at(0.5));
 uv_start.x = uv_mid.x; uv_end.x = uv_mid.x;
 } else {
 if at_pole(uv_start.y) { uv_start.x = uv_end.x; }
 if at_pole(uv_end.y) { uv_end.x = uv_start.x; }
 }
 let delta = uv_end - uv_start;
 let span = delta.length();
 if span < TOLERANCE_CLAMP_MIN || !span.is_finite() { return None; }
 Some((
 Curve2d::Line(Line2d {
 origin: uv_start,
 direction: delta / span,
 }),
 span,
 ))
 }
 (Surface3::Cylinder(c), _) => {
 let axis = c.axis.normalize_or_zero();
 if axis.length_squared() < 0.5 { return None; }
 let uv_start = {
 let local = curve.point_at(0.0) - c.origin;
 let v = local.dot(axis);
 let radial = local - axis * v;
 let u = radial.y.atan2(radial.x);
 DVec2::new(if u < 0.0 { u + std::f64::consts::TAU } else { u }, v)
 };
 let uv_end = {
 let local = curve.point_at(1.0) - c.origin;
 let v = local.dot(axis);
 let radial = local - axis * v;
 let u = radial.y.atan2(radial.x);
 DVec2::new(if u < 0.0 { u + std::f64::consts::TAU } else { u }, v)
 };
 let delta = uv_end - uv_start;
 let span = delta.length();
 if span < TOLERANCE_CLAMP_MIN || !span.is_finite() { return None; }
 Some((
 Curve2d::Line(Line2d {
 origin: uv_start,
 direction: delta / span,
 }),
 span,
 ))
 }
 (Surface3::Cone(c), _) => {
 let axis = c.axis_dir();
 let uv_start = {
 let local = curve.point_at(0.0) - c.apex;
 let along = local.dot(axis);
 let radial = local - axis * along;
 let r = radial.length();
 let u = if r < TOLERANCE_CLAMP_MIN { 0.0 } else { radial.y.atan2(radial.x) };
 DVec2::new(if u < 0.0 { u + std::f64::consts::TAU } else { u }, along)
 };
 let uv_end = {
 let local = curve.point_at(1.0) - c.apex;
 let along = local.dot(axis);
 let radial = local - axis * along;
 let r = radial.length();
 let u = if r < TOLERANCE_CLAMP_MIN { 0.0 } else { radial.y.atan2(radial.x) };
 DVec2::new(if u < 0.0 { u + std::f64::consts::TAU } else { u }, along)
 };
 let delta = uv_end - uv_start;
 let span = delta.length();
 if span < TOLERANCE_CLAMP_MIN || !span.is_finite() { return None; }
 Some((
 Curve2d::Line(Line2d {
 origin: uv_start,
 direction: delta / span,
 }),
 span,
 ))
 }
 //  ?Torus boundary edge pcurve via world_to_uv projection.
 // Same pattern as Sphere/Cylinder/Cone: project endpoints to UV space,
 // construct Line2d approximation.  Non-seam edges are short enough
 // that the chord in UV-space is a valid pcurve approximation.
 (Surface3::Torus(t), _) => {
 let uv_start = t.world_to_uv(curve.point_at(0.0));
 let uv_end = t.world_to_uv(curve.point_at(1.0));
 let delta = uv_end - uv_start;
 let span = delta.length();
 if span < TOLERANCE_CLAMP_MIN || !span.is_finite() { return None; }
 Some((
 Curve2d::Line(Line2d {
 origin: uv_start,
 direction: delta / span,
 }),
 span,
 ))
 }
 // fallback  ?use full sampling+interpolation via
 // make_pcurve_on_surface (IntTools_Curve::MakePCurveOnSurface equivalent).
 // Handles BSpline and other surface types not covered by analytic cases.
 _ => {
 let t_range = [0.0, 1.0];
 let pc = rcad_kernel::projection::make_pcurve_on_surface(curve, t_range, surface, 16)?;
 // Sample start/end UV points for span estimate (BSpline parameter
 // range may not be [0,1]  ?use direct arc evaluation).
 let uv_pts = {
 let mut pts = Vec::new();
 let n = 16usize;
 for i in 0..n {
 let t = i as f64 / (n - 1) as f64;
 pts.push(pc.point_at(t));
 }
 pts
 };
 let span = uv_pts.windows(2)
 .map(|w| (w[1] - w[0]).length())
 .sum::<f64>();
 Some((pc, span))
 }
 }
 }

 ///  ?InitShapeInfo =build flat ShapeInfo array from existing Vecs.
 /// OCCT: BOPDS_DS::InitShapeInfo (BOPDS_DS.cxx L264-309).  Populates myLines
 /// with one BOPDS_ShapeInfo per shape, setting type, sub-shapes, has_brep.
 /// rcad: builds shape_info from vertices/edges/wires/faces/shells arrays.
 /// Flat index: [0..nV) = VERTEX, [nV..nV+nE) = EDGE,
 /// [nV+nE..nV+nE+nW) = WIRE, [nV+nE+nW..nV+nE+nW+nF) = FACE,
 /// [nV+nE+nW+nF..) = SHELL+SOLID+COMPSOLID.
 ///  ?ShapeInfo(index) =access shape info by shapes[] index.
 /// OCCT: BOPDS_DS::ShapeInfo (BOPDS_DS.cxx L255-258).
 pub fn shape_info_at(&self, idx: usize) -> &ShapeInfo {
 &self.shape_info[idx]
 }

 ///  ?NbSourceShapes() =original source shape count.
 /// OCCT: BOPDS_DS::NbSourceShapes (BOPDS_DS.cxx L193-195).
 pub fn nb_source_shapes(&self) -> usize {
 self.nb_source_shapes
 }

 ///  ?ShapeType(index) =type from flat index.
 /// OCCT: ShapeInfo(index).ShapeType().
 pub fn shape_type_of(&self, idx: usize) -> rcad_kernel::topods::ShapeType {
 self.shape_info[idx].shape_type
 }

 ///  ?build per-face pcurve representations for all boundary edges.
 /// Called after edges and faces are loaded (end of DS construction).
 pub fn build_face_reps(&mut self) {
 // For each face, iterate its boundary edges and create DSCurveRepOnFace entries.
 for fi in 0..self.faces.len() {
 let surface = self.faces[fi].surface.clone();
 // precompute a single UV basis per face for planar surfaces
 let pl_basis: Option<(DVec3, DVec3, DVec3)> = match &surface {
  Surface3::Plane(p) => { let u = any_perpendicular(p.normal).normalize(); Some((u, p.normal.cross(u).normalize(), p.origin)) }
  _ => None,
 };
 // Collect all edge indices from outer and inner wires.
 let mut all_ei: Vec<usize> = self.faces[fi].boundary_edges.clone();
 for w in &self.faces[fi].inner_boundary_edges {
 all_ei.extend(w.iter().map(|&(ei, _)| ei));
 }
 for &ei in &all_ei {
 if self.edge_on_face(ei, fi).is_some() { continue; }
 let Some(edge) = self.edges.get_mut(ei) else { continue; };
 let t_range = edge.t_range;
 // Shared UV basis for planar line edges
 if let Some((ref u_axis, ref v_axis, ref pl_origin)) = pl_basis {
  if let Curve3::Line(l) = &edge.curve {
   let diff = l.origin - *pl_origin;
   let origin = DVec2::new(diff.dot(*u_axis), diff.dot(*v_axis));
   let dir = DVec2::new(l.direction.dot(*u_axis), l.direction.dot(*v_axis));
   let len = dir.length();
   if len > TOLERANCE_CLAMP_MIN && (t_range[1] - t_range[0]).abs() > TOLERANCE_CLAMP_MIN {
    let scale = (t_range[1] - t_range[0]) / len;
    edge.face_reps.push(DSCurveRepOnFace {
     face_idx: fi,
     pcurve: Curve2d::Line(Line2d { origin, direction: dir / len * scale }),
     pcurve2: None,
     pcurve_range: [t_range[0], t_range[1]],
     start_param: t_range[0],
     end_param: t_range[1],
    });
    continue;
   }
  }
 }
 // Generic pcurve computation for non-planar or non-line edges
 if let Some((mut pcurve, mut span)) = Self::compute_edge_pcurve(&edge.curve, &surface, pl_basis) {
 let t_span = t_range[1] - t_range[0];
 if t_span.abs() > TOLERANCE_CLAMP_MIN && (span - t_span.abs()).abs() > 1e-12 {
 if let Curve2d::Line(ref mut l) = pcurve { l.direction *= span / t_span; }
 span = t_span;
 }
 edge.face_reps.push(DSCurveRepOnFace {
 face_idx: fi, pcurve, pcurve2: None,
 pcurve_range: [t_range[0], t_range[0] + span],
 start_param: t_range[0],
 end_param: t_range[0] + span,
 });
 }
 }
 }
 }

 pub fn find_vertex_near(&self, point: DVec3, tol: f64) -> Option<usize> {
 let tol2 = tol * tol;
 self.vertices.iter().position(|v| (v.point - point).length_squared() <= tol2)
 }

 /// Add a vertex, deduplicating against existing vertices.
 ///
 /// Coincidence uses `max(fuzzy_tol, TOLERANCE_ABS, each vertex's geom_tol)` so
 /// imported vertex tolerances and pave fuzzy both widen merging.
 pub fn add_vertex(&mut self, point: DVec3) -> usize {
 let new_base = self.fuzzy_tol.max(TOLERANCE_ABS);
 for (i, v) in self.vertices.iter().enumerate() {
 let merge_tol = new_base.max(v.geom_tol);
 if (v.point - point).length() <= merge_tol {
 return i;
 }
 }
 self.add_vertex_no_dedup(point)
 }

 /// Create a new DS vertex at the given position, always creating a
 /// distinct entry even if one exists at the same position.
 /// OCCT creates separate TopoDS_Vertex objects for each intersection
 /// (EF vs VV vs EE), even at the same geometric point.  rcad's
 /// `add_vertex` deduplicates by position, which collapses distinct
 /// intersection vertices into one.  Use `no_dedup` when the caller
 /// needs a distinct vertex entity (e.g. EF intersection at a VV SD
 /// vertex position).
 pub fn add_vertex_no_dedup(&mut self, point: DVec3) -> usize {
 let new_base = self.fuzzy_tol.max(TOLERANCE_ABS);
 let idx = self.vertices.len();
 self.vertices.push(DSVertex {
 point,
 origin: None,
 geom_tol: new_base,
 is_internal: false,
 location: 0,
 });
 // add ShapeInfo for the new vertex (is_new = true)
 // so PutPavesOnCurve's IsNewShape check passes.
 let mut si = crate::bopds::ds::types::ShapeInfo::new(rcad_kernel::topods::ShapeType::Vertex);
 si.is_new = true;
 si.rank = 0;
 si.box_min = Some(point);
 si.box_max = Some(point);
 si.box_gap = new_base + self.fuzzy_tol * 0.5;
 self.shape_info.push(si);
 // Ensure shapes/vertex_shape_idx are consistent (like push_vertex)
 if idx >= self.vertex_shape_idx.len() {
  use topods::tshape_flags;
  self.shapes.push(std::sync::Arc::new(topods::TShape::Vertex(topods::TVertexData {
   my_shapes: Vec::new(), flags: tshape_flags::DEFAULT, point, tolerance: new_base, points: Vec::new(),
  })));
  self.vertex_shape_idx.push(self.shapes.len() - 1);
 }
 idx
 }

 /// add a TopLoc_Location to the DS location pool, return its index.
 /// Returns 0 if the transform is identity (index 0 = identity, implicit).
 pub fn add_location(&mut self, loc: glam::DAffine3) -> u32 {
 if loc == glam::DAffine3::IDENTITY { return 0; }
 // Check if already stored
 for (i, existing) in self.locations.iter().enumerate() {
 if existing.abs_diff_eq(loc, 1e-12) {
 return (i + 1) as u32;
 }
 }
 self.locations.push(loc);
 self.locations.len() as u32
 }

   pub fn allocate_pave_block(&mut self, pb: PaveBlock) -> usize {
   let idx = self.pave_blocks.len();
   self.pave_blocks.push(SharedPB::new(pb));
   idx
   }

   // ----- CommonBlock accessors (BOPDS_DS.hxx L186-193) -----

  pub fn is_common_block(&self, pb: &SharedPB) -> bool {
  pb.0.read().unwrap().common_block_idx.is_some()
  }

  pub fn common_block(&self, pb: &SharedPB) -> Option<&CommonBlock> {
  pb.0.read().unwrap().common_block_idx.and_then(|idx| self.common_blocks.get(idx))
  }

  pub fn common_block_mut(&mut self, pb: &SharedPB) -> Option<&mut CommonBlock> {
  pb.0.read().unwrap().common_block_idx.and_then(|idx| self.common_blocks.get_mut(idx))
  }

  pub fn real_pave_block_edge(&self, edge_idx: usize, pb: &SharedPB) -> Option<usize> {
  let cb = self.common_block(pb)?;
  let first_pb_idx = cb.pave_blocks().first()?.0;
  let pbr = self.edges.get(edge_idx)?.pave_blocks.get(first_pb_idx)?;
  pbr.0.read().unwrap().new_edge
  }

 /// Collect 3D boundary points for a face.
 ///
 /// When the topological wire produces a degenerate polygon (< 3 unique points),
 /// falls back to the UV boundary if available and the face is planar. This ensures
 /// e.g. cylinder caps (2 boundary vertices on a diameter line) produce a proper
 /// circular polygon for face-face intersection clipping.
 pub fn face_boundary_points(&self, face_idx: usize) -> Vec<DVec3> {
 let Some(face) = self.faces.get(face_idx) else { return vec![] };
 let pts: Vec<DVec3> = face
 .boundary_verts
 .iter()
 .map(|&vi| self.vertices[vi].point)
 .collect();
 if pts.len() < 3 {
 if let Some(uv_bnd) = &face.uv_boundary {
 if uv_bnd.len() >= 3 && matches!(face.surface, Surface3::Plane(_)) {
 return uv_bnd
 .iter()
 .map(|uv| face.surface.point_at(uv.x, uv.y))
 .collect();
 }
 }
 }
 pts
 }

 /// Detect shared topology between ShapeA and ShapeB.
 ///
 /// This method populates `self.shared_topology` with information about
 /// coincident vertices, edges, and faces. It should be called after
 /// the DS is fully constructed but before interference detection.
 ///
 /// # Arguments
 /// * `tolerance` - Base distance for glue-style coincidence; combined per pair with
 /// relevant `geom_tol` on vertices, edges, or faces (`max` of all).
 ///
 /// # Returns
 /// A reference to the populated `SharedTopologyInfo`.
 pub fn detect_shared_topology(&mut self, tolerance: f64) -> &SharedTopologyInfo {
 let tol = tolerance.max(TOLERANCE_ABS);

 // Clear any previous data
 self.shared_topology = SharedTopologyInfo::default();

 // Detect shared vertices
 for vi_a in 0..self.a_vertex_count {
 for vi_b in self.a_vertex_count..self.vertices.len() {
 let p_a = self.vertices[vi_a].point;
 let p_b = self.vertices[vi_b].point;
 let pair_tol = tol
 .max(self.vertices[vi_a].geom_tol)
 .max(self.vertices[vi_b].geom_tol);
 let pair_tol_sq = pair_tol * pair_tol;
 if (p_a - p_b).length_squared() <= pair_tol_sq {
 self.shared_topology.shared_vertices.push((vi_a, vi_b));
 }
 }
 }

 // Detect shared edges
 for ei_a in 0..self.a_edge_count {
 for ei_b in self.a_edge_count..self.edges.len() {
 let edge_tol = tol
 .max(self.edges[ei_a].geom_tol)
 .max(self.edges[ei_b].geom_tol);
 if self.edges_geometry_compatible(ei_a, ei_b, edge_tol) {
 self.shared_topology.shared_edges.push((ei_a, ei_b));
 }
 }
 }

 // Detect shared faces (full and partial)
 for fi_a in 0..self.a_face_count {
 for fi_b in self.a_face_count..self.faces.len() {
 if self.faces[fi_a].origin == self.faces[fi_b].origin {
 continue; // Same shape, skip
 }

 let face_tol = tol
 .max(self.faces[fi_a].geom_tol)
 .max(self.faces[fi_b].geom_tol);
 let full_overlap = self.faces_boundary_fully_overlap(fi_a, fi_b, face_tol);
 let partial_overlap =
 !full_overlap && self.faces_share_edges(fi_a, fi_b, face_tol);

 if full_overlap {
 self.shared_topology.fully_glued_faces.push((fi_a, fi_b));
 self.shared_topology.shared_faces.push((fi_a, fi_b));
 } else if partial_overlap {
 self.shared_topology.partially_glued_faces.push((fi_a, fi_b));
 self.shared_topology.shared_faces.push((fi_a, fi_b));
 }
 }
 }

 &self.shared_topology
 }

 /// Check if two edges have compatible geometry.
 fn edges_geometry_compatible(&self, e1: usize, e2: usize, tol: f64) -> bool {
 let edge1 = &self.edges[e1];
 let edge2 = &self.edges[e2];

 // Check curve compatibility
 match (&edge1.curve, &edge2.curve) {
 (Curve3::Line(l1), Curve3::Line(l2)) => {
 // Lines must be collinear
 let d1 = l1.direction.normalize_or_zero();
 let d2 = l2.direction.normalize_or_zero();
 if d1.dot(d2).abs() < 0.999 {
 return false;
 }
 // Origins must be on the same line
 let v = l2.origin - l1.origin;
 let perp = v - d1 * v.dot(d1);
 if perp.length() > tol {
 return false;
 }
 // Check parameter ranges overlap
 let _p1_start = l1.origin + d1 * edge1.t_range[0];
 let _p1_end = l1.origin + d1 * edge1.t_range[1];
 let p2_start = l2.origin + d2 * edge2.t_range[0];
 let p2_end = l2.origin + d2 * edge2.t_range[1];

 // Project edge2 endpoints onto edge1 line
 let t2_start = (p2_start - l1.origin).dot(d1);
 let t2_end = (p2_end - l1.origin).dot(d1);
 let (t2_min, t2_max) = if t2_start < t2_end {
 (t2_start, t2_end)
 } else {
 (t2_end, t2_start)
 };

 // Check for overlap
 t2_min <= edge1.t_range[1] + tol && t2_max >= edge1.t_range[0] - tol
 }
 (Curve3::Circle(c1), Curve3::Circle(c2)) => {
 // Circles must be same radius and coplanar
 (c1.center - c2.center).length() <= tol
 && c1.normal.dot(c2.normal).abs() >= 0.999
 && (c1.radius - c2.radius).abs() <= tol
 }
 (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
 // Ellipses must be same
 (e1.center - e2.center).length() <= tol
 && e1.normal.dot(e2.normal).abs() >= 0.999
 && (e1.major_radius - e2.major_radius).abs() <= tol
 && (e1.minor_radius - e2.minor_radius).abs() <= tol
 }
 _ => false,
 }
 }

 /// Check if two faces have fully overlapping boundaries.
 fn faces_boundary_fully_overlap(&self, f1: usize, f2: usize, tol: f64) -> bool {
 let pts1 = self.face_boundary_points(f1);
 let pts2 = self.face_boundary_points(f2);

 if pts1.len() < 3 || pts2.len() < 3 {
 return false;
 }

 // Each point in pts1 must have a matching point in pts2
 let tol_sq = tol * tol;
 let mut used = vec![false; pts2.len()];

 for p1 in &pts1 {
 let mut found = false;
 for (j, p2) in pts2.iter().enumerate() {
 if used[j] {
 continue;
 }
 if (*p1 - *p2).length_squared() <= tol_sq {
 used[j] = true;
 found = true;
 break;
 }
 }
 if !found {
 return false;
 }
 }

 true
 }

 /// Check if two faces share any edges.
 fn faces_share_edges(&self, f1: usize, f2: usize, tol: f64) -> bool {
 let edges1: std::collections::HashSet<usize> =
 self.faces[f1].boundary_edges.iter().copied().collect();
 let edges2: std::collections::HashSet<usize> =
 self.faces[f2].boundary_edges.iter().copied().collect();

 // Check for geometry-compatible edges
 for &e1 in &edges1 {
 for &e2 in &edges2 {
 if self.edges_geometry_compatible(e1, e2, tol) {
 return true;
 }
 }
 }

 false
 }

 /// Check if a face pair is fully glued (can skip intersection entirely).
 pub fn is_fully_glued_face_pair(&self, f1: usize, f2: usize) -> bool {
 self.shared_topology
 .fully_glued_faces
 .iter()
 .any(|&(a, b)| (a == f1 && b == f2) || (a == f2 && b == f1))
 }

 /// Check if a face pair is partially glued (has shared edges).
 pub fn is_partially_glued_face_pair(&self, f1: usize, f2: usize) -> bool {
 self.shared_topology
 .partially_glued_faces
 .iter()
 .any(|&(a, b)| (a == f1 && b == f2) || (a == f2 && b == f1))
 }

 /// Get shared vertices for a face pair.
 pub fn get_shared_vertices_for_faces(&self, f1: usize, f2: usize) -> Vec<(usize, usize)> {
 let boundary1: std::collections::HashSet<usize> =
 self.faces[f1].boundary_verts.iter().copied().collect();
 let boundary2: std::collections::HashSet<usize> =
 self.faces[f2].boundary_verts.iter().copied().collect();

 self.shared_topology
 .shared_vertices
 .iter()
 .filter(|(v1, v2)| {
 (boundary1.contains(v1) && boundary2.contains(v2))
 || (boundary1.contains(v2) && boundary2.contains(v1))
 })
 .copied()
 .collect()
 }

 ///  ?Build edge images from pave blocks (BOPAlgo_Builder::FillImagesEdges).
 ///
 /// Reads `pb.new_edge` from each source edge's PaveBlocks to populate
 /// `my_images` / `my_origins` mappings.  Sub-edges are already created by
 /// `build_split_edges` (PaveFiller::MakeSplitEdges) =this function only
 /// constructs the mapping table, it does NOT create new edges.
 ///
 /// This must be called after `build_split_edges()` (end of `make_blocks`).
 pub fn build_edge_images(&mut self) {
 let n_edges = self.edges.len();
 self.my_images = vec![Vec::new(); n_edges];
 self.my_origins = Vec::new();

  for ei in 0..n_edges {
  let edge = &self.edges[ei];
  for spb in &edge.pave_blocks {
  let pbb = spb.0.read().unwrap();
  let sub_ei = pbb.new_edge.unwrap_or(ei);
 if sub_ei < self.edges.len() {
 self.my_images[ei].push(sub_ei);
 self.my_origins.push(ei);
 }
 }
 }
 }

 ///  ?FillImagesContainers (BOPAlgo_Builder_1.cxx L172-276).
 /// For each original wire whose edges were split by the PaveFiller,
 /// build a new edge list from the split sub-edges.
 ///
 /// Uses DS internal data only 閿?no external BRep needed.
 pub fn build_container_images(&mut self) {
 // Count total wires across all solids/shells in the DS
 let n_wires: usize = self.solids.iter()
 .flat_map(|s| &s.shells)
 .flat_map(|sh| &self.shells[*sh].faces)
 .map(|&fi| 1 + self.faces[fi].inner_boundary_edges.len())
 .sum();
 self.wire_images = vec![None; n_wires];

 // Shell images: flag shells whose faces have any split edges
 let n_shells: usize = self.shells.len();
 self.shell_images = vec![false; n_shells];
 self.solid_images = vec![false; self.solids.len()];

 for (si, solid) in self.solids.iter().enumerate() {
 for &shi in &solid.shells {
 let shell = &self.shells[shi];
 let shell_has_split = shell.faces.iter().any(|&fi| {
 let face = &self.faces[fi];
 Self::wire_has_split_edges_ds(&face.boundary_edges, &self.my_images)
 || face.inner_boundary_edges.iter().any(|iw| {
 let iw_edges: Vec<usize> = iw.iter().map(|&(ei, _)| ei).collect();
 Self::wire_has_split_edges_ds(&iw_edges, &self.my_images)
 })
 });
 self.shell_images[shi] = shell_has_split;
 if shell_has_split {
 self.solid_images[si] = true;
 }
 }
 }

 let mut wi = 0usize;
 for solid in &self.solids {
 for &shi in &solid.shells {
 let shell = &self.shells[shi];
 for &fi in &shell.faces {
 let face = &self.faces[fi];
 // Outer wire
 let new_outer = Self::rebuild_wire_edges_ds(&face.boundary_edges, &face.boundary_edge_forwards, &self.my_images);
 if new_outer.is_some() {
 self.wire_images[wi] = new_outer;
 }
 wi += 1;

 // Inner wires
 for iw in &face.inner_boundary_edges {
 let iw_edges: Vec<usize> = iw.iter().map(|&(ei, _)| ei).collect();
 let iw_fwd: Vec<bool> = iw.iter().map(|&(_, fwd)| fwd).collect();
 let new_inner = Self::rebuild_wire_edges_ds(&iw_edges, &iw_fwd, &self.my_images);
 if new_inner.is_some() {
 self.wire_images[wi] = new_inner;
 }
 wi += 1;
 }
 }
 }
 }
 }

 /// Rebuild a wire's edge list, replacing split edges with their sub-edges.
 /// Uses DS edge indices directly (no external BRep needed).
 /// Returns None if no edge was split (wire unchanged).
 fn rebuild_wire_edges_ds(
 edges: &[usize],
 forwards: &[bool],
 my_images: &[Vec<usize>],
 ) -> Option<Vec<(usize, bool)>> {
 let mut new_edges = Vec::new();
 let mut changed = false;
 for (&ei, &fwd) in edges.iter().zip(forwards.iter()) {
 if ei < my_images.len() && !my_images[ei].is_empty() {
 changed = true;
 for &sub_ei in &my_images[ei] {
 new_edges.push((sub_ei, fwd));
 }
 } else {
 new_edges.push((ei, fwd));
 }
 }
 if changed { Some(new_edges) } else { None }
 }

 /// Check if any edge in a wire has been split by the PaveFiller.
 /// Uses DS edge indices directly.
 fn wire_has_split_edges_ds(
 edges: &[usize],
 my_images: &[Vec<usize>],
 ) -> bool {
 edges.iter().any(|&ei| ei < my_images.len() && !my_images[ei].is_empty())
 }

 /// Get the Plane surface for a face (panics if face is not a plane).
 /// Used by the PaveFiller to compute pcurves for coplanar overlap ICs.
 pub fn face_plane(&self, fi: usize) -> Plane {
 match &self.faces[fi].surface {
 Surface3::Plane(p) => *p,
 _ => panic!("DS::face_plane: face {} is not a Plane surface", fi),
 }
 }

 ///  ?BOPDS_DS::RefineFaceInfoOn.
 ///
 /// Removes PaveBlocks from the On set that are degenerate
 /// (pave1.vertex_idx == pave2.vertex_idx =start and end vertices are the
 /// same, so the PaveBlock has zero length and does not contribute to face
 /// splitting).
 pub fn refine_face_info_on(&mut self, fi: usize) {
 let pave_blocks = &self.pave_blocks;
 let info = &mut self.faces[fi].face_info;
 info.pave_blocks_on.retain(|&pb_idx| {
 pave_blocks.get(pb_idx).map_or(false, |pb| {
 pb.0.read().unwrap().pave1.vertex_idx != pb.0.read().unwrap().pave2.vertex_idx
 })
 });
 }

 /// OCCT BOPDS_DS::ReleasePaveBlocks (BOPDS_DS.cxx L1503-1550).
 /// Clears PaveBlocks for edges with exactly 1 PB and no CommonBlock,
 /// so untouched edges do not get images. Does NOT clear face_info PBs.
 pub fn release_pave_blocks(&mut self) {
   for ei in 0..self.edges.len() {
     let pb_count = self.edge_pave_blocks(ei).len();
     if pb_count != 1 { continue; }
     let has_cb = self.edge_pave_blocks(ei).iter().any(|spb| {
       spb.0.read().unwrap().common_block_idx.is_some()
     });
     if has_cb { continue; }
     self.edge_pave_blocks_mut(ei).clear();
   }
 }

 ///  ?BOPDS_DS::RefineFaceInfoIn.
 ///
 /// Removes PaveBlocks from the In set that ALSO appear in the On set.
 /// A PaveBlock is considered "the same" if it has the same original edge
 /// index and the same start/end vertices (matching OCCT's IsPaveBlockOn
 /// check of OriginalEdge + Pave1.IsEqual + Pave2.IsEqual).
 ///
 /// The On classification takes priority =a PaveBlock classified as On
 /// does not need to also be classified as In.
 pub fn refine_face_info_in(&mut self, fi: usize) {
 let pave_blocks = &self.pave_blocks;
 let on_set = self.faces[fi].face_info.pave_blocks_on.clone();
 let info = &mut self.faces[fi].face_info;
 info.pave_blocks_in.retain(|&pb_idx| {
 let pb = match pave_blocks.get(pb_idx) {
 Some(pb) => pb,
 None => return false,
 };
 // Keep only if NOT in On (same edge index + pave bounds)
 !on_set.iter().any(|&on_idx| {
 pave_blocks.get(on_idx).map_or(false, |on_pb| {
  on_pb.0.read().unwrap().original_edge == pb.0.read().unwrap().original_edge
  && on_pb.0.read().unwrap().pave1.vertex_idx == pb.0.read().unwrap().pave1.vertex_idx
  && on_pb.0.read().unwrap().pave2.vertex_idx == pb.0.read().unwrap().pave2.vertex_idx
 })
 })
 });
 }

 /// BOPDS_DS::GetSameDomainIndex (BOPDS_DS.cxx L1244-1253).
 /// Resolves the vertex index through the SD chain to its canonical root.
 /// In rcad, SD targets always have smaller indices than sources
 /// (make_sd_vertices_vv picks the minimum index as merge target),
 /// so only follow when partner < result to avoid bidirectional bounce.
 pub fn get_same_domain_index(&self, vi: usize) -> usize {
  let mut result = vi;
  loop {
   match self.shape_sd.find_sd_partner(result) {
    Some(next) if next < result => result = next,
    _ => break,
   }
  }
  result
 }

 /// BOPDS_DS::FaceInfoIn (BOPDS_DS.cxx L837-889).
 /// Clears and repopulates face_info.vertices_in and pave_blocks_in from:
 ///   1. Face boundary vertices with GetSameDomainIndex (OCCT L843-852)
 ///   2. VF interference vertices (OCCT L854-864)
 ///   3. EF interference vertices (OCCT L866-877)
 ///   4. FF IC endpoint vertices (for SubShapesOnIn common-vertex detection)
 /// Must be called before MakeBlocks so SubShapesOnIn projects onto curves.
 /// 閴?Clear+recompute, GetSameDomainIndex on all vertices.
 pub fn update_face_info_in(&mut self, fi: usize) {
  // Resolve all SD indices upfront to avoid borrow conflicts with mutable self access
  let boundary_copy: Vec<usize> = self.faces[fi].boundary_verts.iter()
   .map(|&vi| self.get_same_domain_index(vi))
   .collect();

  let vf_vertices: Vec<usize> = self.interf_vf.iter()
   .filter(|inf| inf.face == fi)
   .map(|inf| self.get_same_domain_index(inf.vertex))
   .collect();

  let ef_new_vertices: Vec<usize> = self.interf_ef.iter()
   .filter(|inf| inf.face == fi
           && inf.new_vertex != usize::MAX
           && inf.new_vertex < self.vertices.len())
   .map(|inf| self.get_same_domain_index(inf.new_vertex))
   .collect();

  let ff_curve_endpoints: Vec<usize> = self.interf_ff.iter()
   .filter(|inf| inf.f1 == fi || inf.f2 == fi)
   .flat_map(|inf| inf.curves.iter().copied())
   .filter(|&ci| ci < self.intersection_curves.len())
   .flat_map(|ci| {
    let mut v = Vec::with_capacity(2);
    if self.intersection_curves[ci].start_vertex < self.vertices.len() {
     v.push(self.intersection_curves[ci].start_vertex);
    }
    if self.intersection_curves[ci].end_vertex < self.vertices.len() {
     v.push(self.intersection_curves[ci].end_vertex);
    }
    v
   })
   .collect();
  let ff_curve_endpoints: Vec<usize> = ff_curve_endpoints.into_iter()
   .map(|vi| self.get_same_domain_index(vi))
   .collect();

  // OCCT L784-787: Clear then refill
  let info = &mut self.faces[fi].face_info;
  info.vertices_in.clear();

  for &vi in &boundary_copy { info.vertices_in.insert(vi); }
  for &vi in &vf_vertices { info.vertices_in.insert(vi); }
  for &n_v in &ef_new_vertices { info.vertices_in.insert(n_v); }
  for &n_v in &ff_curve_endpoints { info.vertices_in.insert(n_v); }
 }

 /// UpdateFaceInfoOn (BOPDS_DS.cxx L792-807 + L811-833).
 pub fn update_face_info_on(&mut self, fi: usize) {
  let boundary_edges: Vec<usize> = {
   let mut edges: Vec<usize> = self.faces[fi].boundary_edges.clone();
   for w in &self.faces[fi].inner_boundary_edges {
    edges.extend(w.iter().map(|&(ei, _)| ei));
   }
   edges
  };
  let boundary_verts: Vec<usize> = self.faces[fi].boundary_verts.iter()
   .map(|&vi| self.get_same_domain_index(vi))
   .collect();
  let mut edge_pb_vertices: Vec<usize> = Vec::new();
  for &ei in &boundary_edges {
   if ei >= self.edges.len() { continue; }
   for spb in &self.edges[ei].pave_blocks {
    let pb = spb.0.read().unwrap();
    edge_pb_vertices.push(self.get_same_domain_index(pb.pave1.vertex_idx));
    edge_pb_vertices.push(self.get_same_domain_index(pb.pave2.vertex_idx));
   }
  }
  let info = &mut self.faces[fi].face_info;
  info.pave_blocks_on.clear();
  info.vertices_on.clear();
  for &vi in &edge_pb_vertices { info.vertices_on.insert(vi); }
  for &vi in &boundary_verts { info.vertices_on.insert(vi); }
 }



 ///  ?batch refine for all faces.
 ///
 /// Calls `refine_face_info_on` and `refine_face_info_in` for every face
 /// in the DS.  This should be called after all interferences have been
 /// computed and before face splitting.
 /// UpdatePaveBlocksWithSDVertices (BOPDS_DS.cxx L200-280).
 ///  ?BOPDS_DS::UpdatePaveBlocksWithSDVertices.
 /// Replace PaveBlock endpoint vertex indices with their SD (same-domain)
 /// canonical equivalents.  SD vertex pairs indicate geometrically coincident
 /// vertices between operands A and B; using the canonical index ensures
 /// PaveBlocks from SD edges share vertex indices for correct connectivity.
 ///
 /// OCCT BOPAlgo_PaveFiller_10.cxx L166-246: iterates PaveBlocks pool,
 /// for each endpoint checks ShapeSD and replaces with the first (lower)
 /// index in the SD pair.  rcad: iterates all pave_blocks on all edges
 /// plus the global pave_blocks pool.
 pub fn update_pave_blocks_with_sd_vertices(&mut self) {
 // OCCT BOPAlgo_PaveFiller_10.cxx L166-246:
 // ShapesSD() is a one-directional DataMap<source, sd_vertex>.
 // For each entry (source → sd_vertex), replace source with sd_vertex
 // in all PaveBlocks.  There is never a reverse entry (sd_vertex → source).
 //
 // rcad's ShapeSD stores bidirectionally for lookup convenience.
 // To match OCCT semantics, we must infer the direction:
 //   - If exactly one vertex is is_new, the new one is the SD target.
 //   - If both are old (shared topology), the lower index is the target.
 // Then only insert source → target (never target → source).
 let sd_pairs: Vec<(usize, usize)> = self.shape_sd.sd_vertices_iter().copied().collect();
 if sd_pairs.is_empty() { return; }
 let mut replace: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 for &(a, b) in &sd_pairs {
 // Determine which vertex is the SD target (the one that survives)
 let target = if self.is_new_vertex(a) && !self.is_new_vertex(b) {
 a
 } else if self.is_new_vertex(b) && !self.is_new_vertex(a) {
 b
 } else {
 a.min(b) // both old or both new: lower index per OCCT convention
 };
 // Only insert source → target (the replacement direction)
 let source = if a == target { b } else { a };
 replace.entry(source).or_insert(target);
 }

  // Apply SD replacement through SharedPB.
  for edge in &mut self.edges {
  for spb in &mut edge.pave_blocks {
  let mut pb = spb.0.write().unwrap();
  if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
  pb.pave1.vertex_idx = rep;
  }
  if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
  pb.pave2.vertex_idx = rep;
  }
  }
  }
  // Apply to global pool (curve and orphan PBs).
  for spb in &mut self.pave_blocks {
  let mut pb = spb.0.write().unwrap();
  if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
  pb.pave1.vertex_idx = rep;
  }
  if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
  pb.pave2.vertex_idx = rep;
  }
  }
 }

 ///  ?BOPAlgo_PaveFiller::UpdateCommonBlocksWithSDVertices.
 /// Update CommonBlocks after PaveBlock SD vertex replacement.
 /// OCCT iterates the PaveBlocks pool and updates each CommonBlock's
 /// referenced PaveBlocks' vertex indices.  rcad: for non-destructive
 /// mode (which is always false in rcad), OCCT calls
 /// UpdatePaveBlocksWithSDVertices and returns =rcad does the same.
 /// (CommonBlocks are rare in rcad and their pave_block indices are
 /// edge-local, making full iterative update non-trivial.)
 pub fn update_common_blocks_with_sd_vertices(&mut self) {
 // OCCT L175-178: if !myNonDestructive =UpdatePaveBlocksWithSDVertices + return
 // rcad: NonDestructive is always false, so re-use the PB update.
 self.update_pave_blocks_with_sd_vertices();
 }

 /// Apply vertex replacement to all PaveBlocks (both edge-local and global).
  /// batch refine for all faces.
 pub fn refine_all_face_info(&mut self) {
 for fi in 0..self.faces.len() {
 self.refine_face_info_on(fi);
 self.refine_face_info_in(fi);
 }
 }

 /// BOPDS_DS::SubShapesOnIn (BOPDS_DS.cxx L1066-1143).
 /// Collects PBs and vertices of two faces that are ON or IN,
 /// and identifies common PBs/vertices between them.
 pub fn sub_shapes_on_in(
 &self,
 the_face_index1: usize,
 the_face_index2: usize,
 the_mv_on_in: &mut std::collections::HashSet<usize>,
 the_mv_common: &mut std::collections::HashSet<usize>,
 the_pb_on_in: &mut std::collections::HashSet<usize>,
 the_common_pave_blocks: &mut std::collections::HashSet<usize>,
 ) {
 let a_face_info1 = &self.faces[the_face_index1].face_info;
 let a_face_info2 = &self.faces[the_face_index2].face_info;

 // Helper: process PBs from a face's pave_blocks_on/in set
 let process_map = |pb_set: &indexmap::IndexSet<usize>,
 pb_on_in: &mut std::collections::HashSet<usize>,
 mv_on_in: &mut std::collections::HashSet<usize>| {
 for &pb_idx in pb_set {
 pb_on_in.insert(pb_idx);
 if pb_idx < self.pave_blocks.len() {
 let (v1, v2) = self.pave_blocks[pb_idx].0.read().unwrap().indices();
 mv_on_in.insert(v1);
 mv_on_in.insert(v2);
 }
 }
 };

 // Process all four pave-block maps
 process_map(&a_face_info1.pave_blocks_on, the_pb_on_in, the_mv_on_in);
 process_map(&a_face_info1.pave_blocks_in, the_pb_on_in, the_mv_on_in);
 process_map(&a_face_info2.pave_blocks_on, the_pb_on_in, the_mv_on_in);
 process_map(&a_face_info2.pave_blocks_in, the_pb_on_in, the_mv_on_in);

 // Find common PBs (Face1 PBs that are also in Face2)
 for &pb_idx in &a_face_info1.pave_blocks_on {
 if a_face_info2.pave_blocks_on.contains(&pb_idx) || a_face_info2.pave_blocks_in.contains(&pb_idx) {
 the_common_pave_blocks.insert(pb_idx);
 if pb_idx < self.pave_blocks.len() {
 let (v1, v2) = self.pave_blocks[pb_idx].0.read().unwrap().indices();
 the_mv_common.insert(v1);
 the_mv_common.insert(v2);
 }
 }
 }
 for &pb_idx in &a_face_info1.pave_blocks_in {
 if a_face_info2.pave_blocks_on.contains(&pb_idx) || a_face_info2.pave_blocks_in.contains(&pb_idx) {
 the_common_pave_blocks.insert(pb_idx);
 if pb_idx < self.pave_blocks.len() {
 let (v1, v2) = self.pave_blocks[pb_idx].0.read().unwrap().indices();
 the_mv_common.insert(v1);
 the_mv_common.insert(v2);
 }
 }
 }

 // OCCT L1124-1142: vertices from Face1 that are also in Face2
 for &vi in &a_face_info1.vertices_on {
 if a_face_info2.vertices_on.contains(&vi) || a_face_info2.vertices_in.contains(&vi) {
 the_mv_on_in.insert(vi);
 the_mv_common.insert(vi);
 }
 }
 for &vi in &a_face_info1.vertices_in {
 if a_face_info2.vertices_on.contains(&vi) || a_face_info2.vertices_in.contains(&vi) {
 the_mv_on_in.insert(vi);
 the_mv_common.insert(vi);
 }
 }
 }

 /// BOPDS_DS::SharedEdges (BOPDS_DS.cxx L1147-1208).
 /// Collects edges that are shared between two faces.
 pub fn shared_edges(
 &self,
 the_face_index1: usize,
 the_face_index2: usize,
 the_edge_list: &mut Vec<usize>,
 ) {
 let mut a_first_face_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

 // Collect edges of the first face.
 for &ei in &self.faces[the_face_index1].boundary_edges {
 let pbs = self.pave_blocks(ei);
 if pbs.is_empty() {
 a_first_face_edges.insert(ei);
 } else {
 for pb in pbs {
 let re = self.real_pave_block_edge(ei, pb)
 .or(pb.0.read().unwrap().new_edge)
 .unwrap_or(ei);
 a_first_face_edges.insert(re);
 }
 }
 }

 // Add edges of the second face if they are in the first one.
 for &ei in &self.faces[the_face_index2].boundary_edges {
 let pbs = self.pave_blocks(ei);
 if pbs.is_empty() {
 if a_first_face_edges.contains(&ei) {
 the_edge_list.push(ei);
 }
 } else {
 for pb in pbs {
 let re = self.real_pave_block_edge(ei, pb)
 .or(pb.0.read().unwrap().new_edge)
 .unwrap_or(ei);
 if a_first_face_edges.contains(&re) {
 the_edge_list.push(re);
 }
 }
 }
 }
 }
}



