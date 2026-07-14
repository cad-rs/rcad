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
            face_normals: Vec::new(), face_origins: Vec::new(), face_info_vec: Vec::new(),
            source_face_idxs: Vec::new(), face_locations: Vec::new(), face_uv_boundary: Vec::new(),
            source_shell_idxs: Vec::new(), source_solid_idxs: Vec::new(), source_compsolid_idxs: Vec::new(),
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
        match std::sync::Arc::make_mut(&mut self.shapes[fi]) {
            topods::TShape::Face(fd) => fd,
            _ => panic!("face_data_mut: {} is not a Face", fi),
        }
    }

    /// Vertex point from TShape (ds.shape(vi) -> TVertexData.point).
    /// NOTE: after shapes reorder, shapes[0..nv) = Vertex entries, so
    /// vertex index vi maps correctly to shapes[vi].
    pub fn vertex_point(&self, vi: usize) -> glam::DVec3 {
        match self.shape(vi) {
            topods::TShape::Vertex(d) => d.point,
            _ => glam::DVec3::ZERO,
        }
    }

    /// Vertex tolerance from TShape.
    pub fn vertex_tolerance(&self, vi: usize) -> f64 {
        match self.shape(vi) {
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

    /// Edge pave_blocks from parallel array.
    pub fn edge_pave_blocks(&self, ei: usize) -> &[crate::bopds::pave::SharedPB] {
        if ei < self.edge_pave_blocks.len() { &self.edge_pave_blocks[ei] }
        else if let Some(e) = self.edges.get(ei) { &e.pave_blocks }
        else { &[] }
    }

    /// Mutable edge pave_blocks from parallel array.
    pub fn edge_pave_blocks_mut(&mut self, ei: usize) -> &mut Vec<crate::bopds::pave::SharedPB> {
        if ei < self.edge_pave_blocks.len() {
            &mut self.edge_pave_blocks[ei]
        } else if ei < self.edges.len() {
            &mut self.edges[ei].pave_blocks
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

    /// Face surface from TShape.
    pub fn face_surface(&self, fi: usize) -> Option<&rcad_kernel::geom::Surface3> {
        match self.shape(fi) {
            topods::TShape::Face(d) => d.surface.as_ref(),
            _ => None,
        }
    }

    /// Face tolerance from TShape.
    pub fn face_tolerance(&self, fi: usize) -> f64 {
        match self.shape(fi) {
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

    /// FaceInfo reference from parallel array.
    pub fn face_info(&self, fi: usize) -> &FaceInfo {
        if fi < self.face_info_vec.len() { &self.face_info_vec[fi] }
        else { &self.faces.get(fi).expect("face_info out of range").face_info }
    }

    /// Mutable FaceInfo from parallel array.
    pub fn face_info_mut(&mut self, fi: usize) -> &mut FaceInfo {
        if fi < self.face_info_vec.len() {
            &mut self.face_info_vec[fi]
        } else if fi < self.faces.len() {
            &mut self.faces[fi].face_info
        } else {
            panic!("face_info_mut: index {} out of bounds", fi);
        }
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

    /// OCCT-aligned: myDS->Shape(n). Returns &TShape at the given flat index.
    pub fn shape(&self, idx: usize) -> &topods::TShape {
        &self.shapes[idx]
    }

    /// OCCT-aligned: mutable shape access for PaveFiller tolerance updates.
    pub fn shape_mut(&mut self, idx: usize) -> &mut topods::TShape {
        std::sync::Arc::make_mut(&mut self.shapes[idx])
    }

    /// OCCT-aligned: myDS->Append — add a new TShape and return its index.
    pub fn append_shape(&mut self, ts: topods::TShape) -> usize {
        let idx = self.shapes.len();
        self.shapes.push(std::sync::Arc::new(ts));
        idx
    }

    /// OCCT-aligned: BOPDS_DS::NbShapes — total count of all shapes.
    pub fn nb_shapes(&self) -> usize {
        self.shapes.len()
    }

    /// OCCT-aligned: myDS->ShapeType(n) == TopAbs_VERTEX.
    pub fn is_vertex(&self, idx: usize) -> bool {
        self.shapes.get(idx).map_or(false, |s| matches!(&**s, topods::TShape::Vertex(_)))
    }

    /// OCCT-aligned: myDS->ShapeType(n) == TopAbs_EDGE.
    pub fn is_edge(&self, idx: usize) -> bool {
        self.shapes.get(idx).map_or(false, |s| matches!(&**s, topods::TShape::Edge(_)))
    }

    /// OCCT-aligned: myDS->ShapeType(n) == TopAbs_FACE.
    pub fn is_face(&self, idx: usize) -> bool {
        self.shapes.get(idx).map_or(false, |s| matches!(&**s, topods::TShape::Face(_)))
    }

    /// OCCT-aligned: push a vertex (myDS->Append) + track in flat array.
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
        vi
    }

    /// OCCT-aligned: push an edge (myDS->Append) + track in flat array.
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
        let face_info = df.face_info.clone();
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
        self.face_info_vec.push(face_info);
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
        fi
    }

    ///  ?OCCT-aligned: BOPDS_ShapeInfo::HasFlag / Flag.
 /// Returns the flag value for an edge index, or 0 if not set.
 /// Flag is stored in shape_info[nv + edge_idx].flag matching OCCT's
 /// per-shape integer flag (BOPDS_ShapeInfo::myFlag).
 pub fn edge_flag(&self, edge_idx: usize) -> i64 {
 let nv = self.vertices.len();
 let si_idx = nv + edge_idx;
 if si_idx < self.shape_info.len() {
 let f = self.shape_info[si_idx].flag;
 if f >= 0 { f } else { 0 }
 } else { 0 }
 }

 ///  ?OCCT-aligned: BOPDS_ShapeInfo::HasFlag(int&)  ?true if flag is set (>= 0).
 pub fn edge_has_flag(&self, edge_idx: usize) -> bool {
 let nv = self.vertices.len();
 let si_idx = nv + edge_idx;
 si_idx < self.shape_info.len() && self.shape_info[si_idx].flag >= 0
 }

 ///  ?OCCT-aligned: BOPDS_ShapeInfo::SetFlag.
 pub fn set_edge_flag(&mut self, edge_idx: usize, flag: i64) {
 let nv = self.vertices.len();
 let si_idx = nv + edge_idx;
 if si_idx < self.shape_info.len() {
 self.shape_info[si_idx].flag = flag;
 }
 }

 ///  ?OCCT-aligned: BRep_Tool::Degenerated(edge) equivalent.
 /// Checks the shape_info flag: a flagged edge is degenerated.
 /// Also falls back to start==end check for edges loaded before
 /// shape_info initialization.
 pub fn is_edge_degenerated(&self, edge_idx: usize) -> bool {
 if edge_idx < self.edges.len() && self.edges[edge_idx].start_vertex == self.edges[edge_idx].end_vertex {
 return true;
 }
 self.edge_has_flag(edge_idx)
 }

 // ----- OCCT BOPDS_DS data layer methods -----

 ///  ?OCCT-aligned: BOPDS_DS::IsNewShape (L228-233).
 /// Returns true if the shape index was appended during intersection
 /// (not part of the original source shapes).  In rcad, vertices
 /// with `origin: None` are intersection-created (new shapes).
 /// Edges with `origin: None` carry the same semantics.
 pub fn is_new_vertex(&self, vi: usize) -> bool {
 if vi < self.shape_info.len() {
 self.shape_info[vi].is_new
 } else {
 self.vertices.get(vi).map_or(true, |v| v.origin.is_none())
 }
 }

 ///  ?OCCT-aligned: BOPDS_DS::Rank (L214-226).
 /// Returns the rank (operand index 0=A, 1=B) of a shape.
 /// 0 for shapes from operand A, 1 for operand B.
 pub fn rank(&self, vi: usize) -> usize {
 if vi < self.a_vertex_count { 0 } else { 1 }
 }

 ///  ?OCCT-aligned: BOPDS_DS::Range (L207-212).
 /// Returns the index range [start, end) for shapes of given type.
 /// rcad: returns [0, a_vertex_count) for A, [a_vertex_count, end) for B.
 pub fn range(&self, is_a: bool) -> (usize, usize) {
 if is_a { (0, self.a_vertex_count) }
 else { (self.a_vertex_count, self.vertices.len()) }
 }

 // ----- OCCT-aligned: PaveBlock pool accessors (BOPDS_DS.hxx L156-177) -----

 ///  ?OCCT-aligned: BOPDS_DS::HasPaveBlocks (hxx:162-164).
 /// Returns true if the edge with the given index has PaveBlocks.
 pub fn has_pave_blocks(&self, edge_idx: usize) -> bool {
 self.edges.get(edge_idx).map_or(false, |e| !e.pave_blocks.is_empty())
 }

 ///  ?OCCT-aligned: BOPDS_DS::ChangePaveBlocks (hxx:172-174).
 /// Returns a mutable reference to the PaveBlocks list for an edge.
 pub fn change_pave_blocks(&mut self, edge_idx: usize) -> &mut Vec<SharedPB> {
 &mut self.edges[edge_idx].pave_blocks
 }

 ///  ?OCCT-aligned: BOPDS_DS::InitPaveBlocks (cxx L437-501).
 /// Creates the initial PaveBlock for a source edge, covering the full
 /// parametric range [t_range[0], t_range[1]].  For closed edges (seam
 /// edges where start == end), both endpoint paves are added to ext_paves
 /// via AppendExtPave (fence-protected), matching OCCT's closed-edge
 /// split preparation (L477-483).
 pub fn init_pave_blocks_for_edge(&mut self, edge_idx: usize) {
 if edge_idx >= self.edges.len() { return; }
 let (sv, ev, tr0, tr1) = {
 let e = &self.edges[edge_idx];
 (e.start_vertex, e.end_vertex, e.t_range[0], e.t_range[1])
 };
 // OCCT L437-445: shape type check  ?rcad edge is always an edge
 // OCCT L447: curve null check  ?rcad edge always has a curve
 // OCCT L449: BRep_Tool::Range  ?rcad: stored in t_range
 // OCCT L451: BRepAdaptor_Curve  ?rcad: not needed
 // OCCT L453-455: TopExp::FirstVertex / LastVertex  ?rcad: start_vertex/end_vertex
 // OCCT L457-467: create PaveBlock and set pave1/pave2/original_edge
 let pv1 = Pave { vertex_idx: sv, param: tr0 };
 let pv2 = Pave { vertex_idx: ev, param: tr1 };
  let pb = SharedPB::new(PaveBlock::new(edge_idx, pv1, pv2));
  // OCCT L469-471: ChangePaveBlocksPool  ?store in DS
  self.edges[edge_idx].pave_blocks = vec![pb];
  // Seed edge.paves with endpoint vertices so split_pave_blocks can
  // create sub-PBs between endpoint and intersection paves.
  self.edges[edge_idx].paves.push(Pave { vertex_idx: sv, param: tr0 });
  self.edges[edge_idx].paves.push(Pave { vertex_idx: ev, param: tr1 });
 // OCCT L473-475: loaded edge check  ?rcad: always new construction
 // OCCT L477-483: closed edges  ?add BOTH endpoint paves to ext_paves
 if sv == ev {
 // OCCT L479: aPB->AppendExtPave(aP1)  ?first endpoint
 self.edges[edge_idx].pave_blocks[0]
 .0.write().unwrap().append_ext_pave(Pave { vertex_idx: sv, param: tr0 });
 // OCCT L481: aPB->AppendExtPave(aP2)  ?second endpoint
 // (fence dedups by vertex_idx, so second push is accepted
 //  because vertex_idx differs from the first push? No, same
 //  vertex  ?the OCCT fence check also uses vertex_idx.
 //  The second AppendExtPave is rejected by fence in both
 //  implementations; OCCT still writes it for form clarity.)
 self.edges[edge_idx].pave_blocks[0]
 .0.write().unwrap().append_ext_pave(Pave { vertex_idx: sv, param: tr1 });
 }
 // OCCT L499: aPaveBlock->Update(myPaveBlocksPool.Appended(), false);
 // OCCT L500: anEdgeInfo.SetReference(myPaveBlocksPool.Length() - 1);
 // rcad: flat index for this edge = nV + edge_idx
 let si_idx = self.vertices.len() + edge_idx;
 if si_idx < self.shape_info.len() {
  self.shape_info[si_idx].reference = edge_idx as i64;
 }
 }

 ///  ?OCCT-aligned: BOPDS_DS::PaveBlocks (hxx:167-169).
 /// Returns a reference to the PaveBlocks list for an edge.
 pub fn pave_blocks(&self, edge_idx: usize) -> &[SharedPB] {
 &self.edges[edge_idx].pave_blocks
 }

 // ----- OCCT HasInterf / HasSubShape equivalents -----

 ///  ?OCCT-aligned: HasSubShape(nV, nE) =check if vertex is a sub-shape of edge.
 /// Returns true when vertex nV is an endpoint of edge nE.
 pub fn edge_has_vertex(&self, nV: usize, nE: usize) -> bool {
 self.edges.get(nE).map_or(false, |e| e.start_vertex == nV || e.end_vertex == nV)
 }

 ///  ?OCCT-aligned: myDS->HasInterf(nV, nE) =checks VE interference exists.
 pub fn has_interf_ve(&self, vi: usize, ei: usize) -> bool {
 self.interf_ve.iter().any(|inf| inf.vertex == vi && inf.edge == ei)
 }

 /// OCCT-aligned: myDS->HasInterf(n1, n2) =checks VV interference exists.
 pub fn has_interf_vv(&self, v1: usize, v2: usize) -> bool {
 self.interf_vv.iter().any(|inf| (inf.v1 == v1 && inf.v2 == v2) || (inf.v1 == v2 && inf.v2 == v1))
 }

 /// OCCT-aligned: AddShapeSD =register dynamic SD mapping between two vertices.
 pub fn add_shape_sd(&mut self, from: usize, to: usize) {
 self.shape_sd.add_sd_vertex(from, to);
 }

  /// OCCT-aligned: HasShapeSD(n, nSD) =find the SD root vertex.
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
  /// body primarily records the form-aligned call structure.
  pub fn init_pave_blocks_for_vertex(&mut self, vertex_idx: usize) {
    // OCCT L1489: anEdgeIndices = myMapVE.Seek(theVertexIndex)
    let incident_edges = self.map_ve.get(&vertex_idx).cloned().unwrap_or_default();
    // OCCT L1490-1492: if !anEdgeIndices return
    if incident_edges.is_empty() {
      return;
    }
    // OCCT L1495-1498: for each edge, ChangePaveBlocks(anEdgeIndex)
    // rcad: per-edge pave_blocks Vec is always allocated at creation time.
    // No additional initialization needed.  The SD vertex will be
    // incorporated into these PaveBlocks by the subsequent
    // update_pave_blocks_with_sd_vertices() call in the pipeline.
    for _ei in &incident_edges {
      // OCCT: ChangePaveBlocks(anEdgeIndex) ensures pool entry exists.
    }
  }

  /// �?OCCT-aligned: build vertex-to-edge map (BOPDS_DS::myMapVE).
  ///   Populates map_ve from edges array.  Called by init_shape_info()
  ///   after all edges are loaded.
  pub fn build_map_ve(&mut self) {
  self.map_ve.clear();
  for (ei, e) in self.edges.iter().enumerate() {
  self.map_ve.entry(e.start_vertex).or_default().push(ei);
  if e.start_vertex != e.end_vertex {
  self.map_ve.entry(e.end_vertex).or_default().push(ei);
  }
  }
  }

 ///  ?OCCT-aligned: myDS->HasInterf(nE1, nE2) =checks EE interference exists.
 pub fn has_interf_ee(&self, e1: usize, e2: usize) -> bool {
 self.interf_ee.iter().any(|inf| (inf.e1 == e1 && inf.e2 == e2) || (inf.e1 == e2 && inf.e2 == e1))
 }

 ///  ?OCCT-aligned: myDS->HasInterf(nV, nF) =checks VF interference exists.
 pub fn has_interf_vf(&self, vi: usize, fi: usize) -> bool {
 self.interf_vf.iter().any(|inf| inf.vertex == vi && inf.face == fi)
 }

 ///  ?OCCT-aligned: myDS->HasInterf(nE, nF) =checks EF interference exists.
 pub fn has_interf_ef(&self, ei: usize, fi: usize) -> bool {
 self.interf_ef.iter().any(|inf| inf.edge == ei && inf.face == fi)
 }

 ///  ?OCCT-aligned: myDS->HasInterf(nF1, nF2) =checks FF interference exists.
  pub fn has_interf_ff(&self, f1: usize, f2: usize) -> bool {
  let (a, b) = if f1 < f2 { (f1, f2) } else { (f2, f1) };
  self.interf_ff.iter().any(|ff| { let (fa, fb) = if ff.f1 < ff.f2 { (ff.f1, ff.f2) } else { (ff.f2, ff.f1) }; fa == a && fb == b })
  }

  /// �?OCCT-aligned: BOPDS_DS::AddInterf (DS.cxx L410-420).
  ///   Global interference pair fence.  Checks (i1, i2) with i1 < i2 against
  ///   myInterfTB.  Returns true if the pair is NEW (first insertion), false
  ///   if it already exists.  Call before adding to any typed interference vec.
  ///   The fence prevents duplicate shape pairs across all interference types.
  pub fn try_add_interf(&mut self, i1: usize, i2: usize) -> bool {
  let key = if i1 < i2 { (i1, i2) } else { (i2, i1) };
  self.interf_tb.insert(key)
  }

 /// OCCT-aligned: dedup FaceFace interferences by (Fmin,Fmax) pair.
 /// Merges curves/points from duplicate entries in interf_ff.
 pub fn dedup_ff_interferences(&mut self) {
 let mut merged: std::collections::HashMap<(usize, usize), (Vec<usize>, Vec<usize>)> = std::collections::HashMap::new();
 for inf in &self.interf_ff {
 let key = if inf.f1 < inf.f2 { (inf.f1, inf.f2) } else { (inf.f2, inf.f1) };
 let entry = merged.entry(key).or_insert((Vec::new(), Vec::new()));
 for &c in &inf.curves { if !entry.0.contains(&c) { entry.0.push(c); } }
 for &p in &inf.points { if !entry.1.contains(&p) { entry.1.push(p); } }
 }
 self.interf_ff.clear();
  for ((f1, f2), (curves, points)) in &merged {
  self.interf_ff.push(InterferenceFF { f1: *f1, f2: *f2, curves: curves.clone(), points: points.clone(), tangent_faces: false });
  self.try_add_interf(*f1, *f2);
  }
 }


 /// OCCT-aligned: myDS->HasInterfShapeSubShapes(nV, nE)  ?checks if
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

 ///  ?OCCT-aligned: find existing vertex within tolerance (PutPaveOnCurve equivalent).
 /// OCCT's IsVertexOnLine checks if a boundary vertex lies on the intersection
 /// curve, then places the EXISTING vertex index on the curve's pave block,
 /// ensuring the section edge reuses the same TopoDS_Vertex.  This tolerance-
 /// based scan achieves the same sharing for rcad's flat vertex array.
 ///  ?OCCT-aligned: access edge's pcurve representation on a specific face.
 /// Returns None when no representation exists for this (edge, face) pair.
 pub fn edge_on_face(&self, edge_idx: usize, face_idx: usize) -> Option<&DSCurveRepOnFace> {
 self.edges.get(edge_idx)?.face_reps.iter().find(|r| r.face_idx == face_idx)
 }

 ///  ?OCCT-aligned: compute pcurve for a boundary edge on its face surface.
 /// Mirrors BRep_Tool::CurveOnSurface for boundary edges.
 /// Returns (pcurve, pcurve_span_length) where pcurve has normalized direction.
 pub(crate) fn compute_edge_pcurve(curve: &Curve3, surface: &Surface3) -> Option<(Curve2d, f64)> {
 match (surface, curve) {
 (Surface3::Plane(p), Curve3::Line(l)) => {
 let u_axis = any_perpendicular(p.normal).normalize();
 let v_axis = p.normal.cross(u_axis).normalize();
 let diff = l.origin - p.origin;
 let origin = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));
 let dir = DVec2::new(l.direction.dot(u_axis), l.direction.dot(v_axis));
 let len = dir.length();
 if len > TOLERANCE_CLAMP_MIN {
 Some((Curve2d::Line(Line2d { origin, direction: dir / len }), len))
 } else { None }
 }
 (Surface3::Plane(p), Curve3::Circle(c)) => {
 let u_axis = any_perpendicular(p.normal).normalize();
 let v_axis = p.normal.cross(u_axis).normalize();
 let diff = c.center - p.origin;
 let center_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));
 let normal_dot = c.normal.dot(p.normal).abs();
 if (normal_dot - 1.0).abs() < TOLERANCE_MESH_LEGACY {
 let perim = std::f64::consts::TAU * c.radius;
 Some((Curve2d::Circle(rcad_kernel::geom::Circle2d { center: center_2d, x_dir: DVec2::X, y_dir: DVec2::Y, radius: c.radius  }), perim))
 } else { None }
 }
 //  ?OCCT-aligned: compute pcurve for curved surfaces by projecting
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
 //  ?OCCT-aligned: Torus boundary edge pcurve via world_to_uv projection.
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
 // OCCT-aligned: fallback  ?use full sampling+interpolation via
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

 ///  ?OCCT-aligned: InitShapeInfo =build flat ShapeInfo array from existing Vecs.
 /// OCCT: BOPDS_DS::InitShapeInfo (BOPDS_DS.cxx L264-309).  Populates myLines
 /// with one BOPDS_ShapeInfo per shape, setting type, sub-shapes, has_brep.
 /// rcad: builds shape_info from vertices/edges/wires/faces/shells arrays.
 /// Flat index: [0..nV) = VERTEX, [nV..nV+nE) = EDGE,
 /// [nV+nE..nV+nE+nW) = WIRE, [nV+nE+nW..nV+nE+nW+nF) = FACE,
 /// [nV+nE+nW+nF..) = SHELL+SOLID+COMPSOLID.
 pub fn init_shape_info(&mut self) {
 self.shape_info.clear();
 let nv = self.vertices.len();
 let ne = self.edges.len();
 let nw = self.wires.len();
 let nf = self.faces.len();
 let nsh = self.shells.len();
 let a_v = self.a_vertex_count;
 let a_e = self.a_edge_count;
 let a_f = self.a_face_count;

  // VERTEX entries
  for vi in 0..nv {
  let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::Vertex);
  si.rank = if vi < a_v { 0 } else { 1 };
  si.source_idx = vi;
  si.is_new = self.vertices[vi].origin.is_none();
  si.box_min = Some(self.vertices[vi].point);
  si.box_max = Some(self.vertices[vi].point);
  // �?OCCT-aligned: Bnd_Box::SetGap(Tol(V) + theAdditionalTolerance)
  //   theAdditionalTolerance = fuzzy_tol * 0.5 (matching OCCT PaveFiller UpdateTolerance).
  si.box_gap = self.vertices[vi].geom_tol + self.fuzzy_tol * 0.5;
  self.shape_info.push(si);
  }
 // EDGE entries
  for ei in 0..ne {
  let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::Edge);
  si.rank = if ei < a_e { 0 } else { 1 };
  si.source_idx = ei;
  si.is_new = ei >= a_e;
  let sv = self.edges[ei].start_vertex;
  let ev = self.edges[ei].end_vertex;
  if sv < self.vertices.len() && ev < self.vertices.len() {
  let p1 = self.vertices[sv].point;
  let p2 = self.vertices[ev].point;
  si.box_min = Some(p1.min(p2));
  si.box_max = Some(p1.max(p2));
  }
  // �?OCCT-aligned: BOPDS_ShapeInfo for an Edge stores its vertex sub-shapes.
  // OCCT stores sub-shape indices as flat shape indices into myLines.
  si.sub_shapes.push(nv + sv);
  if sv != ev {
  si.sub_shapes.push(nv + ev);
  }
  // �?OCCT-aligned: Bnd_Box::SetGap with additional tolerance for edge.
  si.box_gap = self.edges[ei].geom_tol + self.fuzzy_tol * 0.5;
  // �?OCCT-aligned: BOPDS_ShapeInfo::SetFlag for degenerated edges.
  //   OCCT BOPAlgo_Builder_1.cxx prepareFaces: SetFlag(faceIndex) for deg edges.
  //   rcad: set flag >= 0 when start==end vertex (degenerated geometry).
  if sv == ev {
  si.flag = 0;
  }
  self.shape_info.push(si);
  }
 // WIRE entries -- sub-shapes = edge indices (flat).  OCCT order: after EDGE, before FACE.
 // Flat index: [nv+ne .. nv+ne+nw)
 for wi in 0..nw {
 let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::Wire);
 si.has_brep = false;
 si.rank = 0;
 si.source_idx = wi;
 si.is_new = false;
 for &ei in &self.wires[wi].edges {
 si.sub_shapes.push(nv + ei);
 }
 self.shape_info.push(si);
 }
 // FACE entries -- sub-shapes = edge indices (flat)
 for fi in 0..nf {
 let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::Face);
 si.rank = if fi < a_f { 0 } else { 1 };
 si.source_idx = fi;
 si.is_new = false;
 for &ei in &self.faces[fi].boundary_edges {
 si.sub_shapes.push(nv + ei);
 }
 let verts = &self.faces[fi].boundary_verts;
 if !verts.is_empty() {
 let mut mn = glam::DVec3::splat(f64::INFINITY);
 let mut mx = glam::DVec3::splat(f64::NEG_INFINITY);
 for &vi in verts {
 if vi < self.vertices.len() {
 let p = self.vertices[vi].point;
 mn = mn.min(p);
 mx = mx.max(p);
 }
 }
  if mn.is_finite() {
  si.box_min = Some(mn);
  si.box_max = Some(mx);
  }
  }
  // �?OCCT-aligned: Bnd_Box::SetGap with additional tolerance for face.
  si.box_gap = self.faces[fi].geom_tol + self.fuzzy_tol * 0.5;
  self.shape_info.push(si);
  }
 // SHELL entries -- sub-shapes = face indices (flat)
 for shi in 0..nsh {
 let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::Shell);
 si.has_brep = false;
 si.rank = 0;
 si.source_idx = shi;
 si.is_new = false;
 for &fi in &self.shells[shi].faces {
 si.sub_shapes.push(nv + ne + nw + fi);
 }
 self.shape_info.push(si);
 }
 // SOLID entries -- sub-shapes = shell indices (flat)
 let nso = self.solids.len();
 for soi in 0..nso {
 let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::Solid);
 si.has_brep = false;
 si.rank = 0;
 si.source_idx = soi;
 si.is_new = false;
 for &shi in &self.solids[soi].shells {
 si.sub_shapes.push(nv + ne + nw + nf + shi);
 }
 self.shape_info.push(si);
 }
 // COMPSOLID entries -- sub-shapes = solid indices (flat)
 let ncs = self.comp_solids.len();
 for csi in 0..ncs {
 let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::CompSolid);
 si.has_brep = false;
 si.rank = 0;
 si.source_idx = csi;
 si.is_new = false;
 for &soi in &self.comp_solids[csi].solids {
 si.sub_shapes.push(nv + ne + nw + nf + nsh + soi);
 }
  self.shape_info.push(si);
  }
  self.nb_source_shapes = self.shape_info.len();
  // �?OCCT-aligned: build vertex-to-edge map (myMapVE) after all shapes registered.
  self.build_map_ve();
  }

 ///  ?OCCT-aligned: ShapeInfo(index) =access shape info by flat index.
 /// OCCT: BOPDS_DS::ShapeInfo (BOPDS_DS.cxx L255-258).
 pub fn shape_info_at(&self, idx: usize) -> &ShapeInfo {
 &self.shape_info[idx]
 }

 ///  ?OCCT-aligned: NbSourceShapes() =original source shape count.
 /// OCCT: BOPDS_DS::NbSourceShapes (BOPDS_DS.cxx L193-195).
 pub fn nb_source_shapes(&self) -> usize {
 self.nb_source_shapes
 }

 ///  ?OCCT-aligned: ShapeType(index) =type from flat index.
 /// OCCT: ShapeInfo(index).ShapeType().
 pub fn shape_type_of(&self, idx: usize) -> rcad_kernel::topods::ShapeType {
 self.shape_info[idx].shape_type
 }

 ///  ?OCCT-aligned: build per-face pcurve representations for all boundary edges.
 /// Called after edges and faces are loaded (end of DS construction).
 pub fn build_face_reps(&mut self) {
 // For each face, iterate its boundary edges and create DSCurveRepOnFace entries.
 for fi in 0..self.faces.len() {
 let surface = self.faces[fi].surface.clone();
 // Collect all edge indices from outer and inner wires.
 let mut all_ei: Vec<usize> = self.faces[fi].boundary_edges.clone();
 for w in &self.faces[fi].inner_boundary_edges {
 all_ei.extend(w.iter().map(|&(ei, _)| ei));
 }
 for &ei in &all_ei {
 if self.edge_on_face(ei, fi).is_some() { continue; } // already computed
 let Some(edge) = self.edges.get_mut(ei) else { continue; };
 if let Some((mut pcurve, mut span)) = Self::compute_edge_pcurve(&edge.curve, &surface) {
 // Scale pcurve direction to match edge's 3D parameter range.
 // OCCT Geom2d_Curve::Value(t) uses the 3D curve parameter t,
 // not normalized UV-space.  pcurve point_at(raw_t) must give
 // the correct UV for any raw 3D parameter t on the edge.
 let t_span = edge.t_range[1] - edge.t_range[0];
 if t_span.abs() > TOLERANCE_CLAMP_MIN && (span - t_span.abs()).abs() > 1e-12 {
 let scale = span / t_span;
 if let Curve2d::Line(ref mut l) = pcurve {
 l.direction *= scale;
 }
 span = t_span;
 }
 edge.face_reps.push(DSCurveRepOnFace {
 face_idx: fi,
 pcurve,
 pcurve2: None,
 pcurve_range: [0.0, span],
 start_param: 0.0,
 end_param: span,
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
 let idx = self.vertices.len();
 self.vertices.push(DSVertex {
 point,
 origin: None,
 geom_tol: new_base,
 is_internal: false,
 location: 0,
 });
 // OCCT-aligned: add ShapeInfo for the new vertex (is_new = true)
 // so PutPavesOnCurve's IsNewShape check passes.
 let mut si = crate::bopds::ds::types::ShapeInfo::new(rcad_kernel::topods::ShapeType::Vertex);
 si.is_new = true;
 si.rank = 0;
 si.box_min = Some(point);
 si.box_max = Some(point);
 si.box_gap = new_base + self.fuzzy_tol * 0.5;
 self.shape_info.push(si);
 idx
 }

 /// OCCT-aligned: add a TopLoc_Location to the DS location pool, return its index.
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

   // ----- OCCT-aligned: CommonBlock accessors (BOPDS_DS.hxx L186-193) -----

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

 ///  ?OCCT-aligned: Build edge images from pave blocks (BOPAlgo_Builder::FillImagesEdges).
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

 ///  ?OCCT-aligned: FillImagesContainers (BOPAlgo_Builder_1.cxx L172-276).
 /// For each original wire whose edges were split by the PaveFiller,
 /// build a new edge list from the split sub-edges.
 ///
 /// Uses DS internal data only �?no external BRep needed.
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

 ///  ?OCCT-aligned: BOPDS_DS::RefineFaceInfoOn.
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

 ///  ?OCCT-aligned: BOPDS_DS::RefineFaceInfoIn.
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

 /// OCCT-aligned: BOPDS_DS::FaceInfoIn (BOPDS_DS.cxx L837-889).
 /// Populates face_info.vertices_in with:
 ///   1. Face's own boundary vertices (OCCT L843-852)
 ///   2. VF interference vertices (OCCT L854-864)
 ///   3. EF interference vertices (OCCT L866-889)
 /// Must be called before MakeBlocks so SubShapesOnIn projects onto curves.
 pub fn update_face_info_in(&mut self, fi: usize) {
 // Pre-collect data to avoid borrow conflicts
 let boundary_copy: Vec<usize> = self.faces[fi].boundary_verts.clone();
 let vf_vertices: Vec<usize> = self.interf_vf.iter()
 .filter(|inf| inf.face == fi)
 .map(|inf| inf.vertex)
 .collect();
 let ef_new_vertices: Vec<usize> = self.interf_ef.iter()
 .filter(|inf| inf.face == fi && inf.new_vertex != usize::MAX)
 .map(|inf| inf.new_vertex)
 .collect();
 // Collect IC endpoint vertices from FF interferences involving this face.
 // OCCT InitFaceInfoIn's TopoDS_Iterator dynamically captures all VERTEX
 // sub-shapes of a face, including those added during intersection processing.
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
 // Now modify face_info without borrowing self elsewhere
 let info = &mut self.faces[fi].face_info;
 for &vi in &boundary_copy {
 let n_vsd = vi; // OCCT: GetSameDomainIndex �?skip SD lookup for boundary verts
 info.vertices_in.insert(n_vsd);
 }
 for &vi in &vf_vertices {
 let n_vsd = vi;
 info.vertices_in.insert(n_vsd);
 }
 for &n_v in &ef_new_vertices {
 let n_vsd = n_v;
 info.vertices_in.insert(n_vsd);
 }
 for &n_v in &ff_curve_endpoints {
 info.vertices_in.insert(n_v);
 }
 }

 ///  ?OCCT-aligned: batch refine for all faces.
 ///
 /// Calls `refine_face_info_on` and `refine_face_info_in` for every face
 /// in the DS.  This should be called after all interferences have been
 /// computed and before face splitting.
 /// OCCT-aligned: UpdatePaveBlocksWithSDVertices (BOPDS_DS.cxx L200-280).
 ///  ?OCCT-aligned: BOPDS_DS::UpdatePaveBlocksWithSDVertices.
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
 // Build SD vertex replacement map: for each pair (a,b), pick the
 // canonical vertex (the one from ShapeA, i.e., the lower index).
 let sd_pairs: Vec<(usize, usize)> = self.shape_sd.sd_vertices_iter().copied().collect();
 if sd_pairs.is_empty() { return; }
 let mut replace: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 for &(a, b) in &sd_pairs {
 let canon = a.min(b);
 replace.entry(a).or_insert(canon);
 replace.entry(b).or_insert(canon);
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

 ///  ?OCCT-aligned: BOPAlgo_PaveFiller::UpdateCommonBlocksWithSDVertices.
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
  /// OCCT-aligned: batch refine for all faces.
 pub fn refine_all_face_info(&mut self) {
 for fi in 0..self.faces.len() {
 self.refine_face_info_on(fi);
 self.refine_face_info_in(fi);
 }
 }

 /// OCCT-aligned: BOPDS_DS::SubShapesOnIn (BOPDS_DS.cxx L1066-1143).
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

 /// OCCT-aligned: BOPDS_DS::SharedEdges (BOPDS_DS.cxx L1147-1208).
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



