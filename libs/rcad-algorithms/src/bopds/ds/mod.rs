pub mod types;
pub use types::*;

use super::pave::{Pave, PaveBlock, NO_EDGE};
use super::common_block::CommonBlock;
use super::face_info::FaceInfo;
use crate::tolerance::*;
use std::collections::HashMap;
use glam::{DVec2, DVec3};
use rcad_kernel::{BRep, CurveEval, SurfaceEval, WireEdge};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Line2d, Line3, Plane, Surface3, any_perpendicular};

impl DS {
    /// 鉁?OCCT-aligned: BOPDS_ShapeInfo::HasFlag / Flag.
    ///   Returns the flag value for an edge, or 0 if no flag set.
    pub fn edge_flag(&self, edge_idx: usize) -> usize {
        self.edge_flags.get(&edge_idx).copied().unwrap_or(0)
    }

    /// 鉁?OCCT-aligned: BOPDS_ShapeInfo::HasFlag(int&) 鈥?returns true with flag value.
    pub fn edge_has_flag(&self, edge_idx: usize) -> bool {
        self.edge_flags.contains_key(&edge_idx)
    }

    /// 鉁?OCCT-aligned: BOPDS_ShapeInfo::SetFlag.
    pub fn set_edge_flag(&mut self, edge_idx: usize, flag: usize) {
        self.edge_flags.insert(edge_idx, flag);
    }

    /// 鉁?OCCT-aligned: BRep_Tool::Degenerated(edge) equivalent.
    pub fn is_edge_degenerated(&self, edge_idx: usize) -> bool {
        self.edge_has_flag(edge_idx)
            || self.edges[edge_idx].start_vertex == self.edges[edge_idx].end_vertex
    }

    // ----- OCCT BOPDS_DS data layer methods -----

    /// 鉁?OCCT-aligned: BOPDS_DS::IsNewShape (L228-233).
    ///   Returns true if the shape index was appended during intersection
    ///   (not part of the original source shapes).  In rcad, vertices
    ///   with `origin: None` are intersection-created (new shapes).
    ///   Edges with `origin: None` carry the same semantics.
    pub fn is_new_vertex(&self, vi: usize) -> bool {
        if vi < self.shape_info.len() {
            self.shape_info[vi].is_new
        } else {
            self.vertices.get(vi).map_or(true, |v| v.origin.is_none())
        }
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::Rank (L214-226).
    ///   Returns the rank (operand index 0=A, 1=B) of a shape.
    ///   0 for shapes from operand A, 1 for operand B.
    pub fn rank(&self, vi: usize) -> usize {
        if vi < self.a_vertex_count { 0 } else { 1 }
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::Range (L207-212).
    ///   Returns the index range [start, end) for shapes of given type.
    ///   rcad: returns [0, a_vertex_count) for A, [a_vertex_count, end) for B.
    pub fn range(&self, is_a: bool) -> (usize, usize) {
        if is_a { (0, self.a_vertex_count) }
        else { (self.a_vertex_count, self.vertices.len()) }
    }

    // ----- OCCT-aligned: PaveBlock pool accessors (BOPDS_DS.hxx L156-177) -----

    /// 鉁?OCCT-aligned: BOPDS_DS::HasPaveBlocks (hxx:162-164).
    ///   Returns true if the edge with the given index has PaveBlocks.
    pub fn has_pave_blocks(&self, edge_idx: usize) -> bool {
        self.edges.get(edge_idx).map_or(false, |e| !e.pave_blocks.is_empty())
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::ChangePaveBlocks (hxx:172-174).
    ///   Returns a mutable reference to the PaveBlocks list for an edge.
    pub fn change_pave_blocks(&mut self, edge_idx: usize) -> &mut Vec<PaveBlock> {
        &mut self.edges[edge_idx].pave_blocks
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::InitPaveBlocks (cxx L437-501).
    ///   Creates the initial PaveBlock for a source edge, covering the full
    ///   parametric range [t_range[0], t_range[1]]. For closed edges (seam
    ///   edges where start == end), the PB is initialized with two paves at
    ///   different parameters using the same vertex.
    pub fn init_pave_blocks_for_edge(&mut self, edge_idx: usize) {
        if edge_idx >= self.edges.len() { return; }
        let (sv, ev, tr0, tr1) = {
            let e = &self.edges[edge_idx];
            (e.start_vertex, e.end_vertex, e.t_range[0], e.t_range[1])
        };
        if sv >= self.vertices.len() { return; }
        if ev >= self.vertices.len() { return; }
        let pv1 = Pave { vertex_idx: sv, param: tr0 };
        let pv2 = Pave { vertex_idx: ev, param: tr1 };
        let mut pb = PaveBlock::new(edge_idx, pv1, pv2);
        // OCCT L479-483: closed edges 鈥?add the second endpoint with reversed direction
        if sv == ev {
            pb.ext_paves.push(Pave { vertex_idx: sv, param: tr1 });
        }
        self.edges[edge_idx].pave_blocks = vec![pb];
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::PaveBlocks (hxx:167-169).
    ///   Returns a reference to the PaveBlocks list for an edge.
    pub fn pave_blocks(&self, edge_idx: usize) -> &[PaveBlock] {
        &self.edges[edge_idx].pave_blocks
    }

    // ----- OCCT HasInterf / HasSubShape equivalents -----

    /// 鉁?OCCT-aligned: HasSubShape(nV, nE) 鈥?check if vertex is a sub-shape of edge.
    ///   Returns true when vertex nV is an endpoint of edge nE.
    pub fn edge_has_vertex(&self, nV: usize, nE: usize) -> bool {
        self.edges.get(nE).map_or(false, |e| e.start_vertex == nV || e.end_vertex == nV)
    }

    /// 鉁?OCCT-aligned: myDS->HasInterf(nV, nE) 鈥?checks VE interference exists.
    pub fn has_interf_ve(&self, vi: usize, ei: usize) -> bool {
        self.interferences.iter().any(|interf| {
            matches!(interf, Interference::VertexEdge { vertex, edge, .. }
                if *vertex == vi && *edge == ei)
        })
    }

    /// OCCT-aligned: myDS->HasInterf(n1, n2) 鈥?checks VV interference exists.
    pub fn has_interf_vv(&self, v1: usize, v2: usize) -> bool {
        self.interferences.iter().any(|interf| {
            matches!(interf, Interference::VertexVertex { v1: a, v2: b, .. }
                if (*a == v1 && *b == v2) || (*a == v2 && *b == v1))
        })
    }

    /// OCCT-aligned: AddShapeSD 鈥?register dynamic SD mapping between two vertices.
    pub fn add_shape_sd(&mut self, from: usize, to: usize) {
        self.shape_sd.add_sd_vertex(from, to);
    }

    /// OCCT-aligned: HasShapeSD(n, nSD) 鈥?find the SD root vertex.
    pub fn has_shape_sd(&self, v: usize) -> Option<usize> {
        self.shape_sd.find_sd_partner(v)
    }

    /// 鉁?OCCT-aligned: myDS->HasInterf(nE1, nE2) 鈥?checks EE interference exists.
    pub fn has_interf_ee(&self, e1: usize, e2: usize) -> bool {
        self.interferences.iter().any(|interf| {
            matches!(interf, Interference::EdgeEdge { e1: a, e2: b, .. }
                if (*a == e1 && *b == e2) || (*a == e2 && *b == e1))
        })
    }

    /// 鉁?OCCT-aligned: myDS->HasInterf(nV, nF) 鈥?checks VF interference exists.
    pub fn has_interf_vf(&self, vi: usize, fi: usize) -> bool {
        self.interferences.iter().any(|interf| {
            matches!(interf, Interference::VertexFace { vertex, face, .. }
                if *vertex == vi && *face == fi)
        })
    }

    /// 鉁?OCCT-aligned: myDS->HasInterf(nE, nF) 鈥?checks EF interference exists.
    pub fn has_interf_ef(&self, ei: usize, fi: usize) -> bool {
        self.interferences.iter().any(|interf| {
            matches!(interf, Interference::EdgeFace { edge, face, .. }
                if *edge == ei && *face == fi)
        })
    }

    /// 鉁?OCCT-aligned: myDS->HasInterf(nF1, nF2) 鈥?checks FF interference exists.
    pub fn has_interf_ff(&self, f1: usize, f2: usize) -> bool {
        self.interferences.iter().any(|interf| {
            matches!(interf, Interference::FaceFace { f1: a, f2: b, .. }
                if (*a == f1 && *b == f2) || (*a == f2 && *b == f1))
        })
    }

    /// OCCT-aligned: dedup FaceFace interferences by (Fmin,Fmax) pair.
    pub fn dedup_ff_interferences(&mut self) {
        let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        let mut i = 0;
        while i < self.interferences.len() {
            let do_remove = match &self.interferences[i] {
                Interference::FaceFace { f1, f2, .. } => {
                    let key = if *f1 < *f2 { (*f1, *f2) } else { (*f2, *f1) };
                    if !seen.insert(key) {
                        if let Interference::FaceFace { curves, points, .. } = &self.interferences[i] {
                            let c_add = curves.clone();
                            let p_add = points.clone();
                            for e in &mut self.interferences {
                                if let Interference::FaceFace { f1: fa, f2: fb, curves: ec, points: ep } = e {
                                    let ek = if *fa < *fb { (*fa, *fb) } else { (*fb, *fa) };
                                    if ek == key {
                                        for &c in &c_add { if !ec.contains(&c) { ec.push(c); } }
                                        for &p in &p_add { if !ep.contains(&p) { ep.push(p); } }
                                        break;
                                    }
                                }
                            }
                        }
                        true
                    } else { false }
                }
                _ => false,
            };
            if do_remove { self.interferences.swap_remove(i); }
            else { i += 1; }
        }
    }


    /// 鉁?OCCT-aligned: myDS->HasInterfShapeSubShapes(nV, nE) 鈥?checks if
    ///   vertex already has interference with any sub-shape (face) of the edge.
    ///   In OCCT, edges belong to faces; rcad doesn't track edge鈫抐ace ancestry
    ///   directly 鈥?this is a best-effort check using available face data.
    pub fn has_interf_ve_via_faces(&self, vi: usize, ei: usize) -> bool {
        // Check if the vertex has interference with any face that references this edge
        self.interferences.iter().any(|interf| {
            match interf {
                Interference::VertexFace { vertex, face } if *vertex == vi => {
                    self.faces.get(*face).map_or(false, |f| {
                        f.boundary_edges.contains(&ei)
                            || f.inner_boundary_edges.iter().any(|iw| iw.iter().any(|&(e, _)| e == ei))
                    })
                }
                Interference::EdgeFace { edge, face, .. } if *edge == ei => {
                    // Already have EF with this face; check if same vertex also has VF
                    self.interferences.iter().any(|interf2| {
                        matches!(interf2, Interference::VertexFace { vertex, face: f }
                            if *vertex == vi && *f == *face)
                    })
                }
                _ => false,
            }
        })
    }

    /// Build DS from two BReps using the default absolute tolerance.
    pub fn new(a: &BRep, b: &BRep) -> Self {
        Self::new_with_fuzzy(a, b, crate::tolerance::TOLERANCE_ABS)
    }

    /// Build DS with a caller-supplied fuzzy tolerance.
    ///
    /// `fuzzy_tol` must be 鈮?`TOLERANCE_ABS`; smaller values are clamped up.
    pub fn new_with_fuzzy(a: &BRep, b: &BRep, fuzzy_tol: f64) -> Self {
        let tol = fuzzy_tol.max(crate::tolerance::TOLERANCE_ABS);
        let mut ds = DS {
            vertices: Vec::new(),
            edges: Vec::new(),
            wires: Vec::new(),
            shells: Vec::new(),
            faces: Vec::new(),
            // internal V/E tracking: is_internal flag on DSVertex/DSEdge
            //   used instead of separate arrays (removed in favor of flags).
            interferences: Vec::new(),
            intersection_curves: Vec::new(),
            section_edge_refs: Vec::new(),
            fuzzy_tol: tol,
            a_vertex_count: 0,
            a_edge_count: 0,
            a_face_count: 0,
            shared_topology: SharedTopologyInfo::default(),
            shape_sd: ShapeSD::new(0, &SharedTopologyInfo::default()),
            same_domain_overlaps: Vec::new(),
            common_blocks: Vec::new(),
            my_images: Vec::new(),
            my_origins: Vec::new(),
            wire_images: Vec::new(),
            shell_images: Vec::new(),
            solid_images: Vec::new(),
            pave_blocks: Vec::new(),
            edge_flags: EdgeFlagMap::new(),
            increased_ss: std::collections::HashSet::new(),
            shape_info: Vec::new(),
            nb_source_shapes: 0,
        };

        ds.load_brep(a, ShapeOrigin::ShapeA);
        ds.a_vertex_count = ds.vertices.len();
        ds.a_edge_count = ds.edges.len();
        ds.a_face_count = ds.faces.len();
        ds.load_brep(b, ShapeOrigin::ShapeB);
        ds.compute_uv_boundaries();
        ds.build_face_reps();

        // OCCT-aligned: initialize ShapeInfo flat array after all shapes loaded.
        ds.init_shape_info();

        ds
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
                    // Sphere param from projection: u = longitude [-蟺, 蟺] (atan2 range),
                    // v = colatitude [0, 蟺]. Use the full domain as UV boundary.
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
                    // Cylinder param: u = azimuth [0, 2蟺] (matches CylindricalSurface::point_at),
                    // v = height along axis.  Estimate height range from boundary edge samples.
                    let boundary_edges = self.faces[fi].boundary_edges.clone();
                    let mut h_min = f64::INFINITY;
                    let mut h_max = f64::NEG_INFINITY;
                    let axis = cyl.axis.normalize();
                    let origin = cyl.origin;
                    // Seam edges on a cylinder are lines parallel to the axis.  When a face has
                    // multiple seam edges (e.g. a cylinder with explicit front/back seams), the
                    // u-range is bounded by them rather than the full [0, 2蟺].
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
                    seam_u_vals.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
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
                    // Cone param: u = azimuth [0, 2蟺], v = slant distance from the
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
                    // Torus param: u = major angle [0, 2蟺], v = minor angle [0, 2蟺].
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
            // Tolerance 1e-6 * edge_length accounts for projection noise from
            // closest_point_on_surface on BSpline surfaces (Newton iteration).
            let is_planar = matches!(&surface, Surface3::Plane(_))
                || (if let Surface3::BSpline(ref bsp) = surface { rcad_kernel::geom::bspline_is_planar(bsp, 1e-7) } else { false });
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
                        if cross < 1e-6 * len1.max(len2).max(f64::MIN_POSITIVE)
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

    fn load_brep(&mut self, brep: &BRep, origin: ShapeOrigin) {
        let edge_offset = self.edges.len();
        let face_offset = self.faces.len();

        // OCCT-aligned: share vertices at same 3D position between operands.
        let mut local_to_ds: Vec<usize> = Vec::with_capacity(brep.vertices.len());
        for (local_i, v) in brep.vertices.iter().enumerate() {
            let tol = rcad_kernel::vertex_tolerance(brep, local_i).max(TOLERANCE_ABS);
            if let Some(existing) = self.find_vertex_near(v.point, tol) {
                local_to_ds.push(existing);
            } else {
                let vi = self.vertices.len();
                self.vertices.push(DSVertex {
                    point: v.point,
                    origin: Some(origin),
                    geom_tol: rcad_kernel::vertex_tolerance(brep, local_i),
                    is_internal: false,
                });
                local_to_ds.push(vi);
            }
        }

        // Edges
        for (i, edge) in brep.edges.iter().enumerate() {
            let start = local_to_ds[edge.start];
            let end = local_to_ds[edge.end];

            let curve = brep
                .geom
                .edge_curve
                .get(i)
                .and_then(|c| *c)
                .map(|ci| brep.geom.curves[ci].clone())
                .unwrap_or_else(|| {
                    // Fallback: synthesize line from vertices
                    let p0 = brep.vertices[edge.start].point;
                    let p1 = brep.vertices[edge.end].point;
                    let dir = (p1 - p0).normalize();
                    Curve3::Line(Line3 {
                        origin: p0,
                        direction: dir,
                    })
                });

            // Compute parametric range
            let t_range = match &curve {
                Curve3::Line(line) => {
                    let p0 = brep.vertices[edge.start].point;
                    let p1 = brep.vertices[edge.end].point;
                    let t0 = (p0 - line.origin).dot(line.direction);
                    let t1 = (p1 - line.origin).dot(line.direction);
                    [t0, t1]
                }
                _ => brep
                    .geom
                    .edge_curve_range
                    .get(i)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| curve.default_domain()),
            };

            self.edges.push(DSEdge {
                start_vertex: start,
                end_vertex: end,
                curve,
                t_range,
                origin,
                geom_tol: rcad_kernel::edge_tolerance(brep, i),
                paves: Vec::new(),
                pave_blocks: Vec::new(),
                face_reps: Vec::new(),
                is_internal: false,
                vertex_params: {
                    let mut vp = std::collections::HashMap::new();
                    vp.insert(start, t_range[0]);
                    vp.insert(end, t_range[1]);
                    vp
                },
            });
            self.init_pave_blocks_for_edge(self.edges.len() - 1);
        }

        // Faces.  OCCT BOPDS_ShapeInfo tracks source shell/solid/compsolid
        // hierarchy (TopAbs_COMPSOLID 鈫?TopAbs_SOLID 鈫?TopAbs_SHELL 鈫?TopAbs_FACE).
        // rcad: shell_counter, solid_counter assign sequential indices.
        // source_compsolid_idx: use face-by-shell-count matching against
        // brep.compsolid (OCCT preserves TopoDS identity; rcad matches by
        // structural identity: same solid is found by same shell count).
        let mut face_idx = 0usize;
        let mut solid_counter = 0usize;
        let mut shell_counter = 0usize;
        // OCCT-aligned: match flat-iterated solids to compsolid members by
        // sequential index.  OCCT preserves TopoDS identity; rcad's BRep
        // stores compsolid solids value-copied in brep.solids.  When the
        // counts match, solid i 鈫?compsolid solid i.  When no match, None.
        let cs_count = brep.compsolid.as_ref().map_or(0, |cs| cs.solids.len());
        for solid in &brep.solids {
            let compsolid_idx = if cs_count > 0 && solid_counter < cs_count {
                Some(solid_counter)
            } else {
                None
            };
            for shell in &solid.shells {
                let prev_face_count = self.faces.len(); // track DSShell face range
                for face in &shell.faces {
                    let surface = brep
                        .geom
                        .face_surface
                        .get(face_idx)
                        .and_then(|s| *s)
                        .map(|si| brep.geom.surfaces[si].clone())
                        .unwrap_or_else(|| {
                            // Fallback: synthesize plane from face normal
                            // Use first vertex from outer wire, or origin if no wire
                            let origin = if !face.triangles.is_empty() {
                                brep.vertices[face.triangles[0][0]].point
                            } else if !face.outer_wire.edges.is_empty() {
                                let first_edge = &brep.edges[face.outer_wire.edges[0].idx];
                                brep.vertices[first_edge.start].point
                            } else {
                                DVec3::ZERO
                            };
                            Surface3::Plane(Plane {
                                origin,
                                normal: face.normal,
                            })
                        });

                    // OCCT-aligned: wire traversal order (TopExp_Explorer).
                    let boundary_edges_ordered: Vec<(usize, bool)> = Self::reorder_to_wire_order(
                        &face.outer_wire.edges, &brep.edges, edge_offset);
                    let boundary_edges: Vec<usize> = boundary_edges_ordered.iter().map(|&(ei, _)| ei).collect();
                    let boundary_edge_forwards: Vec<bool> = boundary_edges_ordered.iter().map(|&(_, fwd)| fwd).collect();

                    // Trace the wire edges to get ordered boundary vertices.
                    // Wire edges are not necessarily in traversal order;
                    // we must find shared vertices between consecutive edges.
                    let boundary_verts: Vec<usize> = {
                        let edges_in_wire = &face.outer_wire.edges;
                        if edges_in_wire.is_empty() {
                            Vec::new()
                        } else if edges_in_wire.len() == 1 {
                            let e = &brep.edges[edges_in_wire[0].idx];
                            vec![local_to_ds[e.start], local_to_ds[e.end]]
                        } else {
                            // For each consecutive pair of wire edges, find the
                            // shared vertex 鈫?the other vertex of the first edge
                            // is the boundary vertex contributed by that edge.
                            let mut verts = Vec::with_capacity(edges_in_wire.len());
                            for i in 0..edges_in_wire.len() {
                                let next_i = (i + 1) % edges_in_wire.len();
                                let e = &brep.edges[edges_in_wire[i].idx];
                                let en = &brep.edges[edges_in_wire[next_i].idx];

                                // The shared vertex between e and en
                                let shared = if e.start == en.start || e.start == en.end {
                                    e.start
                                } else {
                                    e.end
                                };

                                // The non-shared vertex of e is the boundary vertex
                                let non_shared = if shared == e.start { e.end } else { e.start };
                                verts.push(local_to_ds[non_shared]);
                            }
                            verts
                        }
                    };

                    // 鉁?OCCT-aligned: create DSWire for outer wire (first-class TopAbs_WIRE).
                    let outer_wire_idx = Some(self.wires.len());
                    self.wires.push(DSWire { edges: boundary_edges.clone() });

                    // 鉁?OCCT-aligned: inner wire edges (TopExp_Explorer iterates outer first, then inner).
                    let inner_boundary_edges: Vec<Vec<(usize, bool)>> = face
                        .inner_wires
                        .iter()
                        .map(|wire| {
                            wire.edges
                                .iter()
                                .map(|we| (we.idx + edge_offset, we.forward))
                                .collect()
                        })
                        .collect();

                    // 鉁?OCCT-aligned: create DSWire for each inner wire.
                    let inner_wire_idxs: Vec<usize> = (0..face.inner_wires.len())
                        .map(|_| {
                            let wi = self.wires.len();
                            self.wires.push(DSWire { edges: Vec::new() });
                            wi
                        })
                        .collect();
                    for (ii, wire) in face.inner_wires.iter().enumerate() {
                        self.wires[inner_wire_idxs[ii]].edges = wire
                            .edges
                            .iter()
                            .map(|we| we.idx + edge_offset)
                            .collect();
                    }

                    self.faces.push(DSFace {
                        surface,
                        boundary_verts,
                        boundary_edges,
                        boundary_edge_forwards,
                        inner_boundary_edges,
                        outer_wire_idx,
                        inner_wire_idxs,
                        normal: face.normal,
                        origin,
                        face_info: FaceInfo::default(),
                        source_face_idx: face_idx,
                        geom_tol: rcad_kernel::face_tolerance(brep, face_idx),
                        uv_boundary: None,
                        natural_restriction: true,
                        source_shell_idx: Some(shell_counter),
                        source_solid_idx: Some(solid_counter),
                        source_compsolid_idx: compsolid_idx,
                    });

                    face_idx += 1;
                }
                // 鉁?OCCT-aligned: create DSShell tracking which DS faces belong to each shell.
                let shell_face_idxs: Vec<usize> = (prev_face_count..self.faces.len()).collect();
                if !shell_face_idxs.is_empty() {
                    self.shells.push(DSShell { faces: shell_face_idxs });
                }
                shell_counter += 1;
            }
            solid_counter += 1;
        }

        // OCCT L622-887 (FillInternalShapes Phase 2): internal V/E from source solids.
        //   TopAbs_INTERNAL sub-shapes inside the solid volume.  Currently no source
        //   BRep provides internal shapes (Solid has no internal_* fields); the DS
        //   is_internal flag is reserved for future use when the BRep data model
        //   supports internal sub-shape storage.
        // rcad: no internal shapes to load at this time.

        // 鉁?Transfer pcurves from BRep's edge_pcurves into DS edge face_reps.
        // This preserves the BRep's stored pcurves (proper surface curves) instead
        // of recomputing them via endpoint projection + Line2d approximation.
        // build_face_reps (called after load_brep) skips edges that already have
        // a DSRepOnFace via edge_on_face check.
        let n_brep_faces = brep.solids.iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum::<usize>();
        for ei in 0..brep.edges.len() {
            let Some(pcurves) = brep.geom.edge_pcurves.get(ei) else { continue; };
            if pcurves.is_empty() { continue; }
            let ds_ei = edge_offset + ei;
            let Some(ds_edge) = self.edges.get_mut(ds_ei) else { continue; };
            for pc in pcurves {
                    let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };
                // Find DS face indices whose BRep surface_idx matches this PCurve's
                for bi in 0..n_brep_faces {
                    let Some(Some(si)) = brep.geom.face_surface.get(bi).map(|&s| s) else { continue; };
                    if si != pc.surface_idx { continue; }
                    let ds_fi = face_offset + bi;
                    if ds_edge.face_reps.iter().any(|r| r.face_idx == ds_fi) { continue; }
                    // Compute span as the 2D chord length at the 3D curve's t_range
                    let t_range = ds_edge.t_range;
                    let uv_start = curve2d.point_at(t_range[0]);
                    let uv_end = curve2d.point_at(t_range[1]);
                    let span = (uv_end - uv_start).length();
                    if span < 1e-15 || !span.is_finite() { continue; }
                    ds_edge.face_reps.push(DSRepOnFace {
                        face_idx: ds_fi,
                        pcurve: curve2d.clone(),
                        pcurve2: None,
                        pcurve_range: [0.0, span],
                        start_param: 0.0,
                        end_param: span,
                    });
                }
            }
        }
    }

    /// 鉁?OCCT-aligned: reorder wire edges by traversal order (TopExp_Explorer).
    /// TopExp_Explorer iterates a wire's edges in the order they form the
    /// closed loop: each edge's end_vertex matches the next edge's start_vertex.
    /// The BRep's wire.edges may not be in this order; we rebuild it by
    /// following end -> start vertex adjacency through the edge graph.
    /// Returns (edge_idx + edge_offset, forward_in_wire) pairs in wire traversal order.
    fn reorder_to_wire_order(
        wire_edges: &[rcad_kernel::topology::WireEdge],
        brep_edges: &[rcad_kernel::topology::Edge],
        edge_offset: usize,
    ) -> Vec<(usize, bool)> {
        if wire_edges.len() <= 1 {
            return wire_edges.iter().map(|we| (we.idx + edge_offset, we.forward)).collect();
        }
        // Build vertex -> wire-edge-index adjacency
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, we) in wire_edges.iter().enumerate() {
            let e = &brep_edges[we.idx];
            adj.entry(e.start).or_default().push(i);
            adj.entry(e.end).or_default().push(i);
        }
        // Walk wire by following end -> start adjacency (start from first edge's start vertex)
        let first = wire_edges[0].idx;
        let mut cur = brep_edges[first].start;
        let mut used = vec![false; wire_edges.len()];
        let mut ordered = Vec::with_capacity(wire_edges.len());
        for _ in 0..wire_edges.len() {
            let next_i = adj.entry(cur).or_default().iter().copied()
                .find(|&i| !used[i])
                .expect("wire is not closed -- broken topology");
            used[next_i] = true;
            let we = &wire_edges[next_i];
            ordered.push((we.idx + edge_offset, we.forward));
            let e = &brep_edges[we.idx];
            cur = if e.start == cur { e.end } else { e.start };
        }
        ordered
    }

    /// 鉁?OCCT-aligned: find existing vertex within tolerance (PutPaveOnCurve equivalent).
    /// OCCT's IsVertexOnLine checks if a boundary vertex lies on the intersection
    /// curve, then places the EXISTING vertex index on the curve's pave block,
    /// ensuring the section edge reuses the same TopoDS_Vertex.  This tolerance-
    /// based scan achieves the same sharing for rcad's flat vertex array.
    /// 鉁?OCCT-aligned: access edge's pcurve representation on a specific face.
    ///   Returns None when no representation exists for this (edge, face) pair.
    pub fn edge_on_face(&self, edge_idx: usize, face_idx: usize) -> Option<&DSRepOnFace> {
        self.edges.get(edge_idx)?.face_reps.iter().find(|r| r.face_idx == face_idx)
    }

    /// 鉁?OCCT-aligned: compute pcurve for a boundary edge on its face surface.
    ///   Mirrors BRep_Tool::CurveOnSurface for boundary edges.
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
                if len > 1e-15 {
                    Some((Curve2d::Line(Line2d { origin, direction: dir / len }), len))
                } else { None }
            }
            (Surface3::Plane(p), Curve3::Circle(c)) => {
                let u_axis = any_perpendicular(p.normal).normalize();
                let v_axis = p.normal.cross(u_axis).normalize();
                let diff = c.center - p.origin;
                let center_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));
                let normal_dot = c.normal.dot(p.normal).abs();
                if (normal_dot - 1.0).abs() < 1e-6 {
                    let perim = std::f64::consts::TAU * c.radius;
                    Some((Curve2d::Circle(rcad_kernel::geom::Circle2d { center: center_2d, x_dir: DVec2::X, y_dir: DVec2::Y, radius: c.radius  }), perim))
                } else { None }
            }
            // 鉁?OCCT-aligned: compute pcurve for curved surfaces by projecting
            //   the edge's 3D curve start/end points onto UV space.
            //   OCCT BRep_Tool::CurveOnSurface returns a parametric curve on
            //   the face surface for every boundary edge (BRep_CurveRepresentation).
            //   For edges that are not seam/deg on periodic surfaces, the simple
            //   Line2d approximation from endpoint UV projection is equivalent to
            //   OCCT's stored pcurve (the edge is short enough that the pcurve is
            //   well-approximated by a line segment in UV space).
            (Surface3::Sphere(s), _) => {
                let mut uv_start = s.world_to_uv(curve.point_at(0.0));
                let mut uv_end = s.world_to_uv(curve.point_at(1.0));
                // At sphere poles (V=0 or V=蟺), U is undefined (atan2(0,0) ambiguity).
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
                if span < 1e-15 || !span.is_finite() { return None; }
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
                if span < 1e-15 || !span.is_finite() { return None; }
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
                    let u = if r < 1e-15 { 0.0 } else { radial.y.atan2(radial.x) };
                    DVec2::new(if u < 0.0 { u + std::f64::consts::TAU } else { u }, along)
                };
                let uv_end = {
                    let local = curve.point_at(1.0) - c.apex;
                    let along = local.dot(axis);
                    let radial = local - axis * along;
                    let r = radial.length();
                    let u = if r < 1e-15 { 0.0 } else { radial.y.atan2(radial.x) };
                    DVec2::new(if u < 0.0 { u + std::f64::consts::TAU } else { u }, along)
                };
                let delta = uv_end - uv_start;
                let span = delta.length();
                if span < 1e-15 || !span.is_finite() { return None; }
                Some((
                    Curve2d::Line(Line2d {
                        origin: uv_start,
                        direction: delta / span,
                    }),
                    span,
                ))
            }
            _ => None,
        }
    }

    /// 鉁?OCCT-aligned: InitShapeInfo 鈥?build flat ShapeInfo array from existing Vecs.
    ///   OCCT: BOPDS_DS::InitShapeInfo (BOPDS_DS.cxx L264-309).  Populates myLines
    ///   with one BOPDS_ShapeInfo per shape, setting type, sub-shapes, has_brep.
    ///   rcad: builds shape_info from vertices/edges/faces/shells arrays.
    ///   Flat index: [0..nV) = VERTEX, [nV..nV+nE) = EDGE,
    ///   [nV+nE..nV+nE+nF) = FACE, [nV+nE+nF..) = SHELL+SOLID.
    pub fn init_shape_info(&mut self) {
        self.shape_info.clear();
        let nv = self.vertices.len();
        let ne = self.edges.len();
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
                si.sub_shapes.push(nv + ne + fi);
            }
            self.shape_info.push(si);
        }
        // SOLID entries -- sub-shapes = shell indices (flat)
        for shi in 0..nsh {
            let mut si = ShapeInfo::new(rcad_kernel::topods::ShapeType::Solid);
            si.has_brep = false;
            si.rank = 0;
            si.source_idx = shi;
            si.is_new = false;
            si.sub_shapes.push(nv + ne + nf + shi);
            self.shape_info.push(si);
        }
        self.nb_source_shapes = self.shape_info.len();
    }

    /// 鉁?OCCT-aligned: ShapeInfo(index) 鈥?access shape info by flat index.
    ///   OCCT: BOPDS_DS::ShapeInfo (BOPDS_DS.cxx L255-258).
    pub fn shape_info_at(&self, idx: usize) -> &ShapeInfo {
        &self.shape_info[idx]
    }

    /// 鉁?OCCT-aligned: NbSourceShapes() 鈥?original source shape count.
    ///   OCCT: BOPDS_DS::NbSourceShapes (BOPDS_DS.cxx L193-195).
    pub fn nb_source_shapes(&self) -> usize {
        self.nb_source_shapes
    }

    /// 鉁?OCCT-aligned: ShapeType(index) 鈥?type from flat index.
    ///   OCCT: ShapeInfo(index).ShapeType().
    pub fn shape_type_of(&self, idx: usize) -> rcad_kernel::topods::ShapeType {
        self.shape_info[idx].shape_type
    }

    /// 鉁?OCCT-aligned: build per-face pcurve representations for all boundary edges.
    ///   Called after edges and faces are loaded (end of DS construction).
    pub fn build_face_reps(&mut self) {
        // For each face, iterate its boundary edges and create DSRepOnFace entries.
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
                    if t_span.abs() > 1e-15 && (span - t_span.abs()).abs() > 1e-12 {
                        let scale = span / t_span;
                        if let Curve2d::Line(ref mut l) = pcurve {
                            l.direction *= scale;
                        }
                        span = t_span;
                    }
                    edge.face_reps.push(DSRepOnFace {
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
        });
        idx
    }

    /// 鉁?OCCT-aligned: Append a PaveBlock to the global pool (BOPDS_DS::ChangePaveBlocksPool).
    ///   Returns the index in the global pool.
    pub fn allocate_pave_block(&mut self, pb: PaveBlock) -> usize {
        let idx = self.pave_blocks.len();
        self.pave_blocks.push(pb);
        idx
    }

    // ----- OCCT-aligned: CommonBlock accessors (BOPDS_DS.hxx L186-193) -----

    /// 鉁?OCCT-aligned: BOPDS_DS::IsCommonBlock (hxx:188).
    ///   Returns true if the PaveBlock belongs to a CommonBlock.
    pub fn is_common_block(&self, pb: &PaveBlock) -> bool {
        pb.common_block_idx.is_some()
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::CommonBlock (hxx:192-193).
    ///   Returns a reference to the CommonBlock for a PaveBlock.
    pub fn common_block(&self, pb: &PaveBlock) -> Option<&CommonBlock> {
        pb.common_block_idx.and_then(|idx| self.common_blocks.get(idx))
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::CommonBlock (hxx:192-193) 鈥?mutable.
    pub fn common_block_mut(&mut self, pb: &PaveBlock) -> Option<&mut CommonBlock> {
        pb.common_block_idx.and_then(|idx| self.common_blocks.get_mut(idx))
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::RealPaveBlock (BOPDS_DS.cxx L658-663).
    ///   If the PaveBlock belongs to a CommonBlock, returns the edge index of
    ///   the first PaveBlock in that block (the "real" edge). Otherwise returns
    ///   the given PaveBlock's new_edge.
    pub fn real_pave_block_edge(&self, edge_idx: usize, pb: &PaveBlock) -> Option<usize> {
        let cb = self.common_block(pb)?;
        let first_pb_idx = cb.pave_blocks().first()?.0;
        self.edges.get(edge_idx)
            .and_then(|e| e.pave_blocks.get(first_pb_idx))
            .and_then(|pbr| pbr.new_edge)
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
    ///   relevant `geom_tol` on vertices, edges, or faces (`max` of all).
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

    /// 鉁?OCCT-aligned: Build edge images from pave blocks (BOPAlgo_Builder::FillImagesEdges).
    ///
    /// Reads `pb.new_edge` from each source edge's PaveBlocks to populate
    /// `my_images` / `my_origins` mappings.  Sub-edges are already created by
    /// `build_split_edges` (PaveFiller::MakeSplitEdges) 鈥?this function only
    /// constructs the mapping table, it does NOT create new edges.
    ///
    /// This must be called after `build_split_edges()` (end of `make_blocks`).
    pub fn build_edge_images(&mut self) {
        let n_edges = self.edges.len();
        self.my_images = vec![Vec::new(); n_edges];
        self.my_origins = Vec::new();

        for ei in 0..n_edges {
            let edge = &self.edges[ei];
            for pb in &edge.pave_blocks {
                let sub_ei = pb.new_edge.unwrap_or(ei);
                if sub_ei < self.edges.len() {
                    self.my_images[ei].push(sub_ei);
                    self.my_origins.push(ei);
                }
            }
        }
    }

    /// 鉁?OCCT-aligned: FillImagesContainers (BOPAlgo_Builder_1.cxx L172-276).
    /// For each original wire whose edges were split by the PaveFiller,
    /// build a new edge list from the split sub-edges.
    pub fn build_container_images(&mut self, brep: &BRep) {
        // Count total wires across all solids/shells
        let n_wires: usize = brep.solids.iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .map(|f| 1 + f.inner_wires.len())
            .sum();
        self.wire_images = vec![None; n_wires];

        // Shell images: flag shells whose faces have any split edges
        let n_shells: usize = brep.solids.iter().map(|s| s.shells.len()).sum();
        self.shell_images = vec![false; n_shells];
        self.solid_images = vec![false; brep.solids.len()];

        let mut shi = 0usize;
        for (si, solid) in brep.solids.iter().enumerate() {
            for shell in &solid.shells {
                let shell_has_split = shell.faces.iter().any(|face| {
                    Self::wire_has_split_edges(&face.outer_wire.edges, &self.my_images)
                        || face.inner_wires.iter().any(|iw| {
                            Self::wire_has_split_edges(&iw.edges, &self.my_images)
                        })
                });
                self.shell_images[shi] = shell_has_split;
                if shell_has_split {
                    self.solid_images[si] = true;
                }
                shi += 1;
            }
        }

        let mut wi = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    // Outer wire
                    let new_outer = Self::rebuild_wire_edges(&face.outer_wire.edges, &self.my_images);
                    if new_outer.is_some() {
                        self.wire_images[wi] = new_outer;
                    }
                    wi += 1;

                    // Inner wires
                    for iw in &face.inner_wires {
                        let new_inner = Self::rebuild_wire_edges(&iw.edges, &self.my_images);
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
    /// Returns None if no edge was split (wire unchanged).
    fn rebuild_wire_edges(
        edges: &[WireEdge],
        my_images: &[Vec<usize>],
    ) -> Option<Vec<(usize, bool)>> {
        let mut new_edges = Vec::new();
        let mut changed = false;
        for we in edges {
            if we.idx < my_images.len() && !my_images[we.idx].is_empty() {
                changed = true;
                for &sub_ei in &my_images[we.idx] {
                    new_edges.push((sub_ei, we.forward));
                }
            } else {
                new_edges.push((we.idx, we.forward));
            }
        }
        if changed { Some(new_edges) } else { None }
    }

    /// Check if any edge in a wire has been split by the PaveFiller.
    /// Used by `build_container_images` for shell/solid image computation.
    fn wire_has_split_edges(
        edges: &[WireEdge],
        my_images: &[Vec<usize>],
    ) -> bool {
        edges.iter().any(|we| we.idx < my_images.len() && !my_images[we.idx].is_empty())
    }

    /// Get the Plane surface for a face (panics if face is not a plane).
    /// Used by the PaveFiller to compute pcurves for coplanar overlap ICs.
    pub fn face_plane(&self, fi: usize) -> Plane {
        match &self.faces[fi].surface {
            Surface3::Plane(p) => *p,
            _ => panic!("DS::face_plane: face {} is not a Plane surface", fi),
        }
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::RefineFaceInfoOn.
    ///
    /// Removes PaveBlocks from the On set that are degenerate
    /// (pave1.vertex_idx == pave2.vertex_idx 鈥?start and end vertices are the
    /// same, so the PaveBlock has zero length and does not contribute to face
    /// splitting).
    pub fn refine_face_info_on(&mut self, fi: usize) {
        let pave_blocks = &self.pave_blocks;
        let info = &mut self.faces[fi].face_info;
        info.pave_blocks_on.retain(|&pb_idx| {
            pave_blocks.get(pb_idx).map_or(false, |pb| {
                pb.pave1.vertex_idx != pb.pave2.vertex_idx
            })
        });
    }

    /// 鉁?OCCT-aligned: BOPDS_DS::RefineFaceInfoIn.
    ///
    /// Removes PaveBlocks from the In set that ALSO appear in the On set.
    /// A PaveBlock is considered "the same" if it has the same original edge
    /// index and the same start/end vertices (matching OCCT's IsPaveBlockOn
    /// check of OriginalEdge + Pave1.IsEqual + Pave2.IsEqual).
    ///
    /// The On classification takes priority 鈥?a PaveBlock classified as On
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
                    on_pb.original_edge == pb.original_edge
                        && on_pb.pave1.vertex_idx == pb.pave1.vertex_idx
                        && on_pb.pave2.vertex_idx == pb.pave2.vertex_idx
                })
            })
        });
    }

    /// 鉁?OCCT-aligned: batch refine for all faces.
    ///
    /// Calls `refine_face_info_on` and `refine_face_info_in` for every face
    /// in the DS.  This should be called after all interferences have been
    /// computed and before face splitting.
        /// OCCT-aligned: UpdatePaveBlocksWithSDVertices (BOPDS_DS.cxx L200-280).
    /// 鉁?OCCT-aligned: BOPDS_DS::UpdatePaveBlocksWithSDVertices.
    ///   Replace PaveBlock endpoint vertex indices with their SD (same-domain)
    ///   canonical equivalents.  SD vertex pairs indicate geometrically coincident
    ///   vertices between operands A and B; using the canonical index ensures
    ///   PaveBlocks from SD edges share vertex indices for correct connectivity.
    ///
    ///   OCCT BOPAlgo_PaveFiller_10.cxx L166-246: iterates PaveBlocks pool,
    ///   for each endpoint checks ShapeSD and replaces with the first (lower)
    ///   index in the SD pair.  rcad: iterates all pave_blocks on all edges
    ///   plus the global pave_blocks pool.
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

        // Apply replacement to all PaveBlocks on edges.
        for edge in &mut self.edges {
            for pb in &mut edge.pave_blocks {
                if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
                    pb.pave1.vertex_idx = rep;
                }
                if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
                    pb.pave2.vertex_idx = rep;
                }
            }
        }
        // Apply replacement to global PaveBlocks pool.
        for pb in &mut self.pave_blocks {
            if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
                pb.pave1.vertex_idx = rep;
            }
            if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
                pb.pave2.vertex_idx = rep;
            }
        }
    }

    /// 鉁?OCCT-aligned: BOPAlgo_PaveFiller::UpdateCommonBlocksWithSDVertices.
    ///   Update CommonBlocks after PaveBlock SD vertex replacement.
    ///   OCCT iterates the PaveBlocks pool and updates each CommonBlock's
    ///   referenced PaveBlocks' vertex indices.  rcad: for non-destructive
    ///   mode (which is always false in rcad), OCCT calls
    ///   UpdatePaveBlocksWithSDVertices and returns 鈥?rcad does the same.
    ///   (CommonBlocks are rare in rcad and their pave_block indices are
    ///    edge-local, making full iterative update non-trivial.)
    pub fn update_common_blocks_with_sd_vertices(&mut self) {
        // OCCT L175-178: if !myNonDestructive 鈫?UpdatePaveBlocksWithSDVertices + return
        //   rcad: NonDestructive is always false, so re-use the PB update.
        self.update_pave_blocks_with_sd_vertices();
    }

    /// Apply vertex replacement to all PaveBlocks (both edge-local and global).
    fn apply_sd_vertex_replacement(&mut self, replace: &std::collections::HashMap<usize, usize>) {
        for edge in &mut self.edges {
            for pb in &mut edge.pave_blocks {
                if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
                    pb.pave1.vertex_idx = rep;
                }
                if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
                    pb.pave2.vertex_idx = rep;
                }
            }
        }
        for pb in &mut self.pave_blocks {
            if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
                pb.pave1.vertex_idx = rep;
            }
            if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
                pb.pave2.vertex_idx = rep;
            }
        }
    }

    /// OCCT-aligned: batch refine for all faces.
    pub fn refine_all_face_info(&mut self) {
        for fi in 0..self.faces.len() {
            self.refine_face_info_on(fi);
            self.refine_face_info_in(fi);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom_populate::populate_box_geom;
    use rcad_kernel::tolerance::{set_edge_tolerance, set_face_tolerance, set_vertex_tolerance};
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn ds_load_brep_copies_geom_tolerances_into_pool() {
        let mut a = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        populate_box_geom(&mut a);
        set_vertex_tolerance(&mut a, 0, 2.0 * TOLERANCE_ADAPTIVE_MAX);
        set_edge_tolerance(&mut a, 0, 3e-3);
        set_face_tolerance(&mut a, 0, 4e-3);

        let b = BRep::default();
        let ds = DS::new(&a, &b);

        assert!((ds.vertices[0].geom_tol - 2.0 * TOLERANCE_ADAPTIVE_MAX).abs() < TOLERANCE_LEN_MIN);
        assert!((ds.edges[0].geom_tol - 3e-3).abs() < TOLERANCE_LEN_MIN);
        assert!((ds.faces[0].geom_tol - 4e-3).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn ds_add_vertex_dedup_respects_fuzzy_tol() {
        let empty = BRep::default();
        let mut ds = DS::new_with_fuzzy(&empty, &empty, TOLERANCE_RETRY_LADDER_COARSE);
        let a = ds.add_vertex(DVec3::ZERO);
        let b = ds.add_vertex(DVec3::new(5e-5, 0.0, 0.0));
        assert_eq!(a, b);
        assert!((ds.vertices[a].geom_tol - TOLERANCE_RETRY_LADDER_COARSE).abs() < TOLERANCE_FLOAT_ULTRA);
    }

    #[test]
    fn ds_from_two_boxes() {
        let mut a = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut b = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        populate_box_geom(&mut a);
        populate_box_geom(&mut b);

        let ds = DS::new(&a, &b);
        assert_eq!(ds.vertices.len(), 16); // 8 + 8
        assert_eq!(ds.edges.len(), 24); // 12 + 12
        assert_eq!(ds.faces.len(), 12); // 6 + 6

        // Check origin tags
        assert!(ds.vertices[0].origin == Some(ShapeOrigin::ShapeA));
        assert!(ds.vertices[8].origin == Some(ShapeOrigin::ShapeB));
        assert!(ds.edges[0].origin == ShapeOrigin::ShapeA);
        assert!(ds.edges[12].origin == ShapeOrigin::ShapeB);

        // OCCT-aligned: BRep_Tool::Parameter 鈥?vertex_params must be populated
        for (ei, e) in ds.edges.iter().enumerate() {
            assert!(e.vertex_params.contains_key(&e.start_vertex),
                "edge {} missing start_vertex {} param", ei, e.start_vertex);
            assert!(e.vertex_params.contains_key(&e.end_vertex),
                "edge {} missing end_vertex {} param", ei, e.end_vertex);
            let sv_p = e.vertex_params.get(&e.start_vertex).copied();
            let ev_p = e.vertex_params.get(&e.end_vertex).copied();
            assert!((sv_p.unwrap() - e.t_range[0]).abs() < 1e-15,
                "edge {} start_vertex param mismatch", ei);
            assert!((ev_p.unwrap() - e.t_range[1]).abs() < 1e-15,
                "edge {} end_vertex param mismatch", ei);
        }
        // vertex_param() convenience method
        assert!(ds.edges[0].vertex_param(ds.edges[0].start_vertex).is_some());
        assert!(ds.edges[0].vertex_param(ds.edges[0].end_vertex).is_some());
        assert!(ds.edges[0].vertex_param(999).is_none());
    }

    #[test]
    fn ds_vertex_params_populated() {
        let mut a = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        populate_box_geom(&mut a);
        let b = a.clone();
        let ds = DS::new(&a, &b);
        for (ei, e) in ds.edges.iter().enumerate() {
            assert!(e.vertex_params.contains_key(&e.start_vertex),
                "edge {} missing start_vertex {} param", ei, e.start_vertex);
            assert!(e.vertex_params.contains_key(&e.end_vertex),
                "edge {} missing end_vertex {} param", ei, e.end_vertex);
            let sv_p = e.vertex_params.get(&e.start_vertex).copied();
            let ev_p = e.vertex_params.get(&e.end_vertex).copied();
            assert!((sv_p.unwrap() - e.t_range[0]).abs() < 1e-15,
                "edge {} start_vertex param mismatch", ei);
            assert!((ev_p.unwrap() - e.t_range[1]).abs() < 1e-15,
                "edge {} end_vertex param mismatch", ei);
        }
        assert!(ds.edges[0].vertex_param(ds.edges[0].start_vertex).is_some());
        assert!(ds.edges[0].vertex_param(999).is_none());
    }

    #[test]
    fn ds_sphere_has_uv_boundary() {
        use rcad_modeling::make_sphere_brep;

        let a = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 0.0, 0.0), 1.0).unwrap();
        let ds = DS::new(&a, &b);

        // Sphere faces should have uv_boundary computed
        let sphere_faces: Vec<_> = ds
            .faces
            .iter()
            .filter(|f| matches!(f.surface, Surface3::Sphere(_)))
            .collect();
        assert!(!sphere_faces.is_empty(), "should have sphere faces");
        for f in &sphere_faces {
            assert!(
                f.uv_boundary.is_some(),
                "sphere face should have uv_boundary"
            );
            let uv = f.uv_boundary.as_ref().unwrap();
            assert!(uv.len() >= 3, "uv boundary should have at least 3 points");
        }
    }

    #[test]
    fn ds_cone_uv_boundary_uses_reference_circle_slant_range() {
        use rcad_modeling::make_cone_brep;

        let a = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0).unwrap();
        let b = BRep::default();
        let ds = DS::new(&a, &b);

        let cone_face = ds
            .faces
            .iter()
            .find(|face| matches!(face.surface, Surface3::Cone(_)))
            .expect("should have a cone face");
        let uv = cone_face
            .uv_boundary
            .as_ref()
            .expect("cone face should have uv_boundary");

        let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        assert!(v_min < 0.0, "expected apex-side slant range below the reference circle, got {v_min}");
        assert!(v_max > 0.0, "expected base-side slant range above the reference circle, got {v_max}");
    }
}
