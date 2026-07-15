use std::collections::{HashSet, BTreeMap};

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Surface3, any_perpendicular};
use rcad_kernel::CurveEval;

use crate::bopalgo::{fill_map, make_blocks};
use crate::bopds::ds::{DS, ShapeOrigin, InterferenceVV, InterferenceVE, InterferenceVF, InterferenceEE, InterferenceEF};
use crate::bopds::pave::Pave;
use crate::inttools;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::pave_filler::helpers::*;
use crate::tolerance::*;

/// OCCT IntTools_CommonPrt::Type() 閳?VERTEX or EDGE.
#[derive(Clone, Copy)]
enum EfHit {
  Vertex { point: DVec3, param: f64 },
  Edge { t1: f64, t2: f64 },
}

/// OCCT L55-93: BOPAlgo_EdgeFace (architecture diff: rcad equivalent).
struct EfTask {
  nE: usize, nF: usize,
  nV1: usize, nV2: usize,
  aT1: f64, aT2: f64,
  aTS1: f64, aTS2: f64,
  bExpressCompute: bool,
  hits: Vec<EfHit>,
}

/// OCCT IntTools_EdgeFace (PaveFiller_5.cxx L340-480): compute edge-face intersection hits.
/// Returns VERTEX-type (point) and EDGE-type (range) common parts.
fn compute_ef_hits(
  ds: &DS, edge_idx: usize, face_idx: usize, ef_range: &[f64; 2],
) -> Vec<EfHit> {
  let edge_curve = &ds.edges[edge_idx].curve;
  let face_surface = &ds.faces[face_idx].surface;
  let etf = ds.edge_tolerance(edge_idx).max(ds.face_tolerance(face_idx)).max(CONFUSION);
  // Phase 1: try point-based (VERTEX-type) intersection
  let mut hits: Vec<EfHit> = match (edge_curve, face_surface) {
    (Curve3::Line(line), Surface3::Plane(plane)) =>
      crate::inttools::edge_face::intersect_line_plane_with_tol(line, *ef_range, plane, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.edge_param }).collect(),
    (Curve3::Line(line), Surface3::Cylinder(cyl)) =>
      crate::inttools::curve_surface::intersect_line_cylinder_with_tol(line, *ef_range, cyl, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.curve_param }).collect(),
    (Curve3::Line(line), Surface3::Sphere(sph)) =>
      crate::inttools::curve_surface::intersect_line_sphere_with_tol(line, *ef_range, sph, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.curve_param }).collect(),
    (Curve3::Line(line), Surface3::Cone(cone)) =>
      crate::inttools::curve_surface::intersect_line_cone_with_tol(line, *ef_range, cone, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.curve_param }).collect(),
    (Curve3::Circle(circle), Surface3::Plane(plane)) => {
      let sv = ds.edge_start_vertex_ds(edge_idx);
      let ref_dir = (ds.vertex_point(sv) - circle.center).normalize();
      crate::inttools::curve_surface::intersect_circle_plane_with_ref(
        circle, *ef_range, plane, etf, Some(ref_dir),
      ).into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.curve_param }).collect()
    }
    (Curve3::Circle(circle), Surface3::Cylinder(cyl)) =>
      crate::inttools::curve_surface::intersect_circle_cylinder_with_tol(circle, *ef_range, cyl, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.curve_param }).collect(),
    (Curve3::Circle(circle), Surface3::Sphere(sph)) =>
      crate::inttools::curve_surface::intersect_circle_sphere_with_tol(circle, *ef_range, sph, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.curve_param }).collect(),
    (Curve3::Circle(circle), Surface3::Cone(cone)) =>
      crate::inttools::curve_surface::intersect_circle_cone_with_tol(circle, *ef_range, cone, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.curve_param }).collect(),
    (Curve3::Ellipse(ellipse), Surface3::Plane(plane)) =>
      crate::inttools::ellipse_intersection::intersect_ellipse_plane_with_tol(ellipse, *ef_range, plane, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.ellipse_param }).collect(),
    (Curve3::Ellipse(ellipse), Surface3::Cylinder(cyl)) =>
      crate::inttools::ellipse_intersection::intersect_ellipse_cylinder_with_tol(ellipse, *ef_range, cyl, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.ellipse_param }).collect(),
    (Curve3::Ellipse(ellipse), Surface3::Sphere(sph)) =>
      crate::inttools::ellipse_intersection::intersect_ellipse_sphere_with_tol(ellipse, *ef_range, sph, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.ellipse_param }).collect(),
    (Curve3::Ellipse(ellipse), Surface3::Cone(cone)) =>
      crate::inttools::ellipse_intersection::intersect_ellipse_cone_with_tol(ellipse, *ef_range, cone, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.ellipse_param }).collect(),
    (Curve3::Parabola(parabola), Surface3::Plane(plane)) =>
      crate::inttools::parabola_intersection::intersect_parabola_plane_with_tol(parabola, *ef_range, plane, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.parabola_param }).collect(),
    (Curve3::Parabola(parabola), Surface3::Cylinder(cyl)) =>
      crate::inttools::parabola_intersection::intersect_parabola_cylinder_with_tol(parabola, *ef_range, cyl, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.parabola_param }).collect(),
    (Curve3::Parabola(parabola), Surface3::Sphere(sph)) =>
      crate::inttools::parabola_intersection::intersect_parabola_sphere_with_tol(parabola, *ef_range, sph, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.parabola_param }).collect(),
    (Curve3::Parabola(parabola), Surface3::Cone(cone)) =>
      crate::inttools::parabola_intersection::intersect_parabola_cone_with_tol(parabola, *ef_range, cone, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.parabola_param }).collect(),
    (Curve3::Hyperbola(hyperbola), Surface3::Plane(plane)) =>
      crate::inttools::hyperbola_intersection::intersect_hyperbola_plane_with_tol(hyperbola, *ef_range, plane, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.hyperbola_param }).collect(),
    (Curve3::Hyperbola(hyperbola), Surface3::Cylinder(cyl)) =>
      crate::inttools::hyperbola_intersection::intersect_hyperbola_cylinder_with_tol(hyperbola, *ef_range, cyl, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.hyperbola_param }).collect(),
    (Curve3::Hyperbola(hyperbola), Surface3::Sphere(sph)) =>
      crate::inttools::hyperbola_intersection::intersect_hyperbola_sphere_with_tol(hyperbola, *ef_range, sph, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.hyperbola_param }).collect(),
    (Curve3::Hyperbola(hyperbola), Surface3::Cone(cone)) =>
      crate::inttools::hyperbola_intersection::intersect_hyperbola_cone_with_tol(hyperbola, *ef_range, cone, etf)
      .into_iter().map(|h| EfHit::Vertex { point: h.point, param: h.hyperbola_param }).collect(),
    _ => {
      let pt_hits = intersect_edge_face_numeric(edge_curve, face_surface, *ef_range, etf);
      pt_hits.into_iter().map(|(p, t)| EfHit::Vertex { point: p, param: t }).collect()
    }
  };
  // Phase 2: if no VERTEX hits, check for EDGE-type (edge coincident with face)
  if hits.is_empty() {
    let on_face = match (edge_curve, face_surface) {
      (Curve3::Line(l), Surface3::Plane(p)) => {
        (l.direction.dot(p.normal)).abs() <= etf
          && crate::inttools::vertex_ops::vertex_on_plane_with_tol(l.origin, p, etf)
      }
      (Curve3::Circle(c), Surface3::Plane(p)) => {
        (c.normal.dot(p.normal)).abs() >= 1.0 - etf
          && crate::inttools::vertex_ops::vertex_on_plane_with_tol(c.center, p, etf)
      }
      _ => {
        let mut all_on = true;
        for i in 0..5 {
          let t = ef_range[0] + (ef_range[1] - ef_range[0]) * (i as f64 / 4.0);
          let pt = edge_curve.point_at(t);
          let proj = rcad_kernel::projection::closest_point_on_surface(face_surface, pt, 8);
          if proj.distance > etf * 10.0 { all_on = false; break; }
        }
        all_on
      }
    };
    if on_face {
      hits.push(EfHit::Edge { t1: ef_range[0], t2: ef_range[1] });
    }
  }
  hits
}

/// OCCT L55-93: BOPAlgo_VertexFace (architecture diff: rcad equivalent).
struct VfTask {
  nV: usize, nF: usize,
  is_on: bool,
  is_on_boundary: bool,
  proj_u: f64, proj_v: f64,
  proj_dist: f64,
}

/// OCCT BOPAlgo_VertexFace: compute vertex projection on face.
fn compute_vf_on_face(ds: &DS, vf_tol: f64, fuzzy_tol: f64, nV: usize, nF: usize) -> VfTask {
  let point = ds.vertex_point(nV);
  let face = &ds.faces[nF];
  let mut r = VfTask { nV, nF, is_on: false, is_on_boundary: false, proj_u: 0.0, proj_v: 0.0, proj_dist: f64::MAX };
  match &face.surface {
    Surface3::Plane(plane) => {
      if crate::inttools::vertex_ops::vertex_on_plane_with_tol(point, plane, vf_tol) {
        let fv = ds.face_boundary_points(nF);
        let on_face = crate::inttools::edge_face::point_in_planar_face_with_tol(point, plane, &fv, vf_tol);
        let u_ax = plane.u_dir; let v_ax = plane.v_dir; let diff = point - plane.origin;
        let u = diff.dot(u_ax); let v = diff.dot(v_ax);
        if on_face { r.is_on = true; r.proj_u = u; r.proj_v = v; }
        else {
          let on_bdy = fv.iter().any(|&p| (point - p).length() <= vf_tol);
          if on_bdy { r.is_on = true; r.is_on_boundary = true; r.proj_u = u; r.proj_v = v; }
        }
      }
    }
    surface => {
      let proj = rcad_kernel::projection::closest_point_on_surface(surface, point, 16);
      let a_tol_v = ds.vertex_tolerance(nV);
      let a_tol_f = face.geom_tol;
      let a_tol_sum = a_tol_v + a_tol_f + fuzzy_tol.max(vf_tol);
      if proj.distance <= a_tol_sum {
        use crate::inttools::fclass2d::State;
        let uv = DVec2::new(proj.params.0, proj.params.1);
        let fclass = crate::inttools::fclass2d::FClass2d::new(ds, nF, vf_tol);
        let st = fclass.perform(uv, false);
        if st == State::In { r.is_on = true; r.proj_u = proj.params.0; r.proj_v = proj.params.1; r.proj_dist = proj.distance; }
        else if st == State::On { r.is_on = true; r.is_on_boundary = true; r.proj_u = proj.params.0; r.proj_v = proj.params.1; r.proj_dist = proj.distance; }
      }
    }
  }
  r
}

impl<'a> super::PaveFiller<'a> {
 /// OCCT PaveFiller L145: BOPDS_Iterator  ?BVH pair enumeration
 pub(crate) fn build_ds_bvh(&self, is_a: bool, is_edge: bool) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let (ds_start, end) = if is_edge {
 if is_a { (0, self.ds.a_edge_count) }
 else { (self.ds.a_edge_count, self.ds.edges.len()) }
 } else {
 if is_a { (0, self.ds.a_vertex_count) }
 else { (self.ds.a_vertex_count, self.ds.vertices.len()) }
 };
 let n = end - ds_start;
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);

 for local_i in 0..n {
 let ds_i = ds_start + local_i;
 indices.push(ds_i);
 let aabb = if is_edge {
 let e = &self.ds.edges[ds_i];
 let pts = [self.ds.vertices[e.start_vertex].point,
 self.ds.vertices[e.end_vertex].point];
 let mut a = Aabb::empty();
 for &p in &pts { a.expand_point(p); }
 // Expand for edge tolerance
 let tol = e.geom_tol.max(CONFUSION);
 a.min -= DVec3::splat(tol);
 a.max += DVec3::splat(tol);
 a
 } else {
 let pt = self.ds.vertex_point(ds_i);
 let tol = self.ds.vertex_tolerance(ds_i).max(CONFUSION);
 Aabb { min: pt - DVec3::splat(tol), max: pt + DVec3::splat(tol) }
 };
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 ///  BOPDS_Iterator  ?build a single BVH for all elements
 /// of the given shape type (both operands A and B combined), used for
 /// single-pass cross-operand pair traversal.
 pub(crate) fn build_ds_bvh_combined(&self, is_edge: bool) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let n = if is_edge { self.ds.edges.len() } else { self.ds.vertices.len() };
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);
 for ds_i in 0..n {
 indices.push(ds_i);
 let aabb = if is_edge {
 let e = &self.ds.edges[ds_i];
 let pts = [self.ds.vertices[e.start_vertex].point,
 self.ds.vertices[e.end_vertex].point];
 let mut a = Aabb::empty();
 for &p in &pts { a.expand_point(p); }
 let tol = e.geom_tol.max(CONFUSION);
 a.min -= DVec3::splat(tol);
 a.max += DVec3::splat(tol);
 a
 } else {
 let pt = self.ds.vertex_point(ds_i);
 let tol = self.ds.vertex_tolerance(ds_i).max(CONFUSION);
 Aabb { min: pt - DVec3::splat(tol), max: pt + DVec3::splat(tol) }
 };
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 /// OCCT PaveFiller_2.cxx L141-206: PerformVE
 /// 閴?BOPDS_Iterator::Initialize(VERTEX, EDGE) 閳?single pass.
 /// Cross-operand filtering is done by BOPDS_Iterator; this function
 /// applies remaining type-specific filters and calls intersect_ve.
 pub(crate) fn perform_ve_bvh(&mut self, pairs: &[(usize, usize)]) {
 use rayon::prelude::*;
 self.fill_shrunk_data();
 // pairs come pre-computed from BOPDS_Iterator
 let ds = &self.ds;
 let a_vc = ds.a_vertex_count;
 let a_ec = ds.a_edge_count;
 if pairs.is_empty() { return; }
 let filtered: Vec<(usize, usize)> = pairs.par_iter()
 .filter(|&(vi, ei)| {
 if (*vi < a_vc) == (*ei < a_ec) { return false; }
 !ds.edge_has_vertex(*vi, *ei) && !ds.edge_has_flag(*ei)
 && !ds.has_interf_ve(*vi, *ei) && !ds.has_interf_ve_via_faces(*vi, *ei)
 && !ds.is_edge_degenerated(*ei)
 && !ds.edges[*ei].pave_blocks.is_empty()
 && ds.edges[*ei].pave_blocks[0].0.read().unwrap().is_splittable
 })
 .copied()
 .collect();
 self.intersect_ve(&filtered);
 }

 ///  IntersectVE (PaveFiller_2.cxx L212-394).
 /// Processes vertex-edge pairs with SD vertex resolution, PB endpoint
 /// dedup (aMVPB), and aDMVSD fence map, matching OCCT's structure.
 fn intersect_ve(&mut self, pairs: &[(usize, usize)]) {
 if pairs.is_empty() { return; }
 // rcad: interferences Vec (no pre-allocation needed)
 // Group vertices by edge, then resolve SD and dedup per (nVSD, nE)
 let mut edge_verts: std::collections::HashMap<usize, Vec<usize>> =
 std::collections::HashMap::new();
 for &(vi, ei) in pairs {
 edge_verts.entry(ei).or_default().push(vi);
 }

 // aDMVSD map: (nVSD, nE)  ?list of original vertices deduped to same SD root
 let mut a_dmv_sd: std::collections::HashMap<(usize, usize), Vec<usize>> =
 std::collections::HashMap::new();
 let mut a_m_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for (&ei, verts) in &edge_verts {
 let a_mv_pb: std::collections::HashSet<usize> = self.ds.edge_paves(ei).iter()
 .map(|p| p.vertex_idx)
 .collect();

 for &vi in verts {
 let n_vsd = self.ds.has_shape_sd(vi).unwrap_or(vi);
 if a_mv_pb.contains(&n_vsd) { continue; }
 let key = (n_vsd, ei);
 a_dmv_sd.entry(key).or_default().push(vi);
 }
 }
 for (&(n_vsd, ei), original_verts) in &a_dmv_sd {
 let point = self.ds.vertex_point(n_vsd);
 let edge = &self.ds.edges[ei];
 let te = self.ve_tol(n_vsd, ei);

 let t_opt = crate::pave_filler::helpers::project_vertex_to_curve(
 point, &edge.curve, te);
 let t = match t_opt {
 Some(t) if t >= edge.t_range[0] && t <= edge.t_range[1] => t,
 _ => continue,
 };
 let dist_3d = edge.curve.point_at(t).distance(point);
 if dist_3d > self.ds.vertex_tolerance(n_vsd) {
 self.ds.vertex_data_mut(n_vsd).tolerance = dist_3d;
 self.ds.increased_ss.insert(n_vsd);
 }
 // OCCT adds pave via aPave.SetIndex(nVx) using the UpdateVertex result.
 // rcad: push Pave directly to edge's pave list.
 let edge_had_paves = !self.ds.edge_paves(ei).is_empty();
 for &vi in original_verts {
 let has_vertex_at_t = self.ds.edge_paves(ei).iter()
 .any(|p| (p.param - t).abs() < TOLERANCE_ABS && p.vertex_idx == vi);
 if !has_vertex_at_t {
 self.ds.edge_paves[ei].push(Pave { vertex_idx: vi, param: t });
 }
 }
 if !edge_had_paves || self.ds.edge_paves(ei).len() > 1 {
 a_m_edges.insert(ei);
 }
 for &vi in original_verts {
 let already_interfered = self.ds.interf_ve.iter().any(|inf| {
 inf.vertex == vi && inf.edge == ei
 });
 if !already_interfered {
 self.ds.interf_ve.push(InterferenceVE{
 vertex: vi,
 edge: ei,
 param: t,
 });
 }
 }
 }
 if !a_m_edges.is_empty() {
 self.split_pave_blocks(&a_m_edges, true);
 }
 }
 /// OCCT PaveFiller_3.cxx L145-244: PerformEE
 /// 閴?BOPDS_Iterator::Initialize(EDGE, EDGE) 閳?single pass.
 /// Cross-operand filtering via a_edge_count.
 /// PerformEE (PaveFiller_3.cxx L145-590).
 /// Returns the set of edges that were modified (got new paves).
 /// The caller should call split_pave_blocks for remaining edges
 /// after treat_new_vertices has processed the new vertex edges.
 pub(crate) fn perform_ee_bvh(&mut self, pairs: &[(usize, usize)])
 -> std::collections::HashSet<usize>
 {
 use rayon::prelude::*;
 self.fill_shrunk_data();
 // pairs come pre-computed from BOPDS_Iterator
 let ds = &self.ds;
 let a_ec = ds.a_edge_count;
 let blocks: Vec<(usize, usize, [f64; 2], [f64; 2])> = pairs.par_iter()
 .filter(|&(ae, be)| {
 if (*ae < a_ec) == (*be < a_ec) { return false; }
 !ds.edge_has_flag(*ae) && !ds.edge_has_flag(*be)
 && !ds.has_interf_ee(*ae, *be)
 && !ds.is_edge_degenerated(*ae) && !ds.is_edge_degenerated(*be)
 })
 .flat_map(|&(ae, be)| {
 let ra = Self::get_pb_boxes(ds, ae, ds.edge_range(ae));
 let rb = Self::get_pb_boxes(ds, be, ds.edge_range(be));
 let mut v = Vec::new();
 for &r1 in &ra { for &r2 in &rb { v.push((ae, be, r1, r2)); } }
 v
 })
 .collect();
 let mut modified: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for &(ae, be, r1, r2) in &blocks {
 self.intersect_ee(ae, be, r1, r2, &mut modified);
 }
 // SplitPaveBlocks) is handled in the caller (perform() in mod.rs):
 // - treat_new_vertices()  ?PerformNewVertices
 // - split_pave_blocks() for remaining modified edges
 modified
 }
 /// OCCT: PaveBlock range extraction (GetPBBox equivalent)
 pub(crate) fn get_pb_boxes(ds: &DS, edge_idx: usize, edge_t_range: [f64; 2]) -> Vec<[f64; 2]> {
 let paves = &ds.edge_paves(edge_idx);
 if paves.is_empty() { return vec![edge_t_range]; }
 let mut params: Vec<f64> = paves.iter().map(|p| p.param).filter(|p| p.is_finite()).collect();
 params.sort_by(|a, b| a.partial_cmp(b).unwrap());
 params.dedup();
 let tol = ds.edge_tolerance(edge_idx).max(crate::tolerance::TOLERANCE_ABS);
 let mut ranges = Vec::new();
 let mut prev = edge_t_range[0];
 for &p in &params {
 if (p - prev).abs() > tol { ranges.push([prev, p]); }
 prev = p;
 }
 if (edge_t_range[1] - prev).abs() > tol { ranges.push([prev, edge_t_range[1]]); }
 ranges
 }
 /// OCCT PaveFiller_4.cxx L139-301: PerformVF
 pub(crate) fn perform_vf_bvh(&mut self, pairs: &[(usize, usize)]) {
 // OCCT L141: myIterator->Initialize(TopAbs_VERTEX, TopAbs_FACE)
 let a_vc = self.ds.a_vertex_count;
 let a_fc = self.ds.a_face_count;
 // OCCT L142-146: iSize check
 let i_size = pairs.len();
 if i_size == 0 { return; }
 //
 // ------------------------------------------------------------------
 // Phase 1: Collect VF tasks (OCCT L181-232)
 // ------------------------------------------------------------------
 // OCCT L174: BOPAlgo_VectorOfVertexFace aVVF
 let mut a_vv_f: Vec<VfTask> = Vec::new();
 // OCCT L180: aMVFPairs dedup map
 let mut a_mvf_pairs: std::collections::HashMap<(usize, usize), Vec<usize>> =
 std::collections::HashMap::new();
 //
 for &(nV, nF) in pairs {
 // OCCT L187: myIterator->Value(nV, nF)
 let same_range = (nV < a_vc) == (nF < a_fc);
 if same_range { continue; }
 // OCCT L194-197: if (myDS->HasInterf(nV, nF)) continue;
 if self.ds.has_interf_vf(nV, nF) { continue; }
 if self.ds.has_interf_ve_via_faces(nV, nF) { continue; }
 // OCCT L205-209: SD resolution
 let nVx = self.ds.has_shape_sd(nV).unwrap_or(nV);
 // OCCT L211-220: aMVFPairs dedup (key = nVx, nF)
 let key = (nVx, nF);
 let entry = a_mvf_pairs.entry(key).or_default();
 entry.push(nV);
 if entry.len() > 1 { continue; }
 // OCCT L222-230: Create BOPAlgo_VertexFace task
 a_vv_f.push(VfTask {
 nV: nVx, nF,
 is_on: false, is_on_boundary: false,
 proj_u: 0.0, proj_v: 0.0, proj_dist: f64::MAX,
 });
 } // for (; myIterator->More(); myIterator->Next()) {
 //
 // ------------------------------------------------------------------
 // Phase 2: BOPAlgo_VertexFace computation (OCCT L234-243)
 // ------------------------------------------------------------------
 // OCCT L242: BOPTools_Parallel::Perform(myRunParallel, aVVF, myContext);
 for task in &mut a_vv_f {
 let tf = self.vf_tol(task.nV, task.nF);
 *task = compute_vf_on_face(&self.ds, tf, self.fuzzy_tolerance, task.nV, task.nF);
 }
 //
 // ------------------------------------------------------------------
 // Phase 3: Process results (OCCT L249-298)
 // ------------------------------------------------------------------
 for task in &a_vv_f {
 if !task.is_on { continue; }
 let nVx = task.nV;
 let nF = task.nF;
 let a_tol_vnew = task.proj_dist;
 // OCCT L272-273: get all original vertices for this SD-face task
 let key = (nVx, nF);
 let orig_verts: Vec<usize> = match a_mvf_pairs.get(&key) {
 Some(v) => v.clone(),
 None => vec![nVx],
 };
 // OCCT L275-297: for each original vertex
 for &nV_orig in &orig_verts {
 // OCCT L277-283: Create InterfVF
 if !self.ds.interf_vf.iter().any(|inf| inf.vertex == nV_orig && inf.face == nF) {
 self.ds.interf_vf.push(InterferenceVF {
 vertex: nV_orig, face: nF, u: task.proj_u, v: task.proj_v,
 });
 }
 // OCCT L285-292: UpdateVertex (simplified)
 if a_tol_vnew.is_finite() && a_tol_vnew > self.ds.vertex_tolerance(nV_orig) {
 if nV_orig < self.ds.vertices.len() {
 self.ds.vertex_data_mut(nV_orig).tolerance = a_tol_vnew;
 self.ds.increased_ss.insert(nV_orig);
 }
 }
 }
 // OCCT L295-297: Add nVx to face's VerticesIn/VerticesOn
 if task.is_on_boundary {
 self.ds.faces[nF].face_info.vertices_on.insert(nVx);
 } else {
 self.ds.face_info_mut(nF).vertices_in.insert(nVx);
 }
 }
 // OCCT L300: TreatVerticesEE() 閳?not yet ported
 }
 /// OCCT PaveFiller_5.cxx L165-592: PerformEF
 pub(crate) fn perform_ef(&mut self, pairs: &[(usize, usize)]) {
 self.fill_shrunk_data(); // OCCT L167
 let a_edge_count = self.ds.a_edge_count;
 let a_face_count = self.ds.a_face_count;
 // OCCT L171-175: iSize check
 if pairs.is_empty() { return; }
 //
 // ------------------------------------------------------------------
 // Phase 1: Collect EF tasks (OCCT L219-307)
 // ------------------------------------------------------------------
 let mut a_v_edge_face: Vec<EfTask> = Vec::new();
 // OCCT L208: aMIEFC
 let mut a_mi_efc: std::collections::HashSet<usize> = std::collections::HashSet::new();
 //
 for &(nE, nF) in pairs {
 let same_range = (nE < a_edge_count && nF < a_face_count)
 || (nE >= a_edge_count && nF >= a_face_count);
 if same_range { continue; }
 // OCCT L227-231: HasFlag / degenerated
 if self.ds.edge_has_flag(nE) || self.ds.is_edge_degenerated(nE) { continue; }
 if self.ds.has_interf_ef(nE, nF) { continue; }
 // OCCT L235: face AABB
 let face_min: DVec3;
 let face_max: DVec3;
 {
 let f = &self.ds.faces[nF];
 let mut mn = DVec3::splat(f64::INFINITY);
 let mut mx = DVec3::splat(f64::NEG_INFINITY);
 for &vi in &f.boundary_verts {
 if vi < self.ds.vertices.len() {
 let p = self.ds.vertex_point(vi); mn = mn.min(p); mx = mx.max(p);
 }
 }
 if let Surface3::Sphere(s) = &f.surface {
 let r = s.radius.abs();
 mn = mn.min(s.center - DVec3::splat(r)); mx = mx.max(s.center + DVec3::splat(r));
 }
 let tol = f.geom_tol.max(CONFUSION);
 face_min = mn - DVec3::splat(tol); face_max = mx + DVec3::splat(tol);
 }
 // OCCT L237-241: FaceInfo
 let face_pbon: Vec<usize> = self.ds.face_info(nF).pave_blocks_on.iter().copied().collect();
 let face_vin: Vec<usize> = self.ds.face_info(nF).vertices_in.iter().copied().collect();
 let face_von: Vec<usize> = self.ds.faces[nF].face_info.vertices_on.iter().copied().collect();
 // OCCT L246-248: iterate edge PBs
 let n_pbs = self.ds.edge_pave_blocks(nE).len();
 for pb_idx in 0..n_pbs {
 // OCCT L256-259: aPBR = RealPaveBlock; if aMPBF.Contains(aPBR) continue;
 let pb = &self.ds.edge_pave_blocks(nE)[pb_idx];
 let pb_ref = pb.0.read().unwrap();
 if face_pbon.contains(&pb_ref.original_edge) { continue; }
 // OCCT L262-266: GetPBBox
 let (aT1, aT2) = pb_ref.range();
 let (aTS1, aTS2, bb_min, bb_max) = match pb_ref.shrunk_range {
 Some(sr) => {
 let (bb_mn, bb_mx) = match pb_ref.my_shrunk_box {
 Some((mn, mx)) => (mn, mx),
 None => {
 let p1 = self.ds.edges[nE].curve.point_at(sr[0]);
 let p2 = self.ds.edges[nE].curve.point_at(sr[1]);
 let tol = self.ds.edge_tolerance(nE).max(CONFUSION);
 (p1.min(p2) - DVec3::splat(tol), p1.max(p2) + DVec3::splat(tol))
 }
 };
 (sr[0], sr[1], bb_mn, bb_mx)
 }
 None => { continue; }
 };
 // OCCT L268-271: AABB overlap check
 if bb_max.x < face_min.x || bb_min.x > face_max.x
 || bb_max.y < face_min.y || bb_min.y > face_max.y
 || bb_max.z < face_min.z || bb_min.z > face_max.z
 { continue; }
 // OCCT L273-276: bExpressCompute
 let (nV1, nV2) = pb_ref.indices();
 let bV1 = face_vin.contains(&nV1) || face_von.contains(&nV1);
 let bV2 = face_vin.contains(&nV2) || face_von.contains(&nV2);
 let b_express_compute = bV1 && bV2;
 // OCCT L278-297: Create task with range correction
 // OCCT L289-292: CorrectRange for shrunk range (aSR 閳?anewSR)
 // OCCT L294-297: CorrectRange for PB range
 //   NOT PORTED: BOPTools_AlgoTools::CorrectRange depends on BRep geometry.
 //   rcad: correct_range_for_face is a simplified approximation.
 drop(pb_ref);
 let tol_ef = self.ef_tol(nE, nF);
 let pb_range = [aT1.min(aT2), aT1.max(aT2)];
 let ef_range = Self::correct_range_for_face(
 &self.ds.edges[nE].curve, tol_ef, pb_range);
 if ef_range[1] - ef_range[0] <= tol_ef { continue; }
 // OCCT L299-305: Save to myFPBDone
 self.fpbdone.entry(nF).or_default().insert(nE);
 //
 a_v_edge_face.push(EfTask {
 nE, nF, nV1, nV2, aT1, aT2, aTS1, aTS2,
 bExpressCompute: b_express_compute,
 hits: Vec::new(),
 });
 } // for (; aIt.More(); aIt.Next())
 } // for (; myIterator->More(); myIterator->Next())
 //
 // ------------------------------------------------------------------
 // Phase 2: EdgeFace computation (OCCT L309-317)
 // ------------------------------------------------------------------
 for task in &mut a_v_edge_face {
 let pb_range = [task.aT1.min(task.aT2), task.aT1.max(task.aT2)];
 let etr = self.ds.edges[task.nE].t_range;
 let ef_range = [pb_range[0].max(etr[0]), pb_range[1].min(etr[1])];
 let tol_ef = self.ef_tol(task.nE, task.nF);
 if ef_range[1] - ef_range[0] <= tol_ef { continue; }
 let ef_corr = Self::correct_range_for_face(
 &self.ds.edges[task.nE].curve, tol_ef, ef_range);
 if ef_corr[1] - ef_corr[0] <= tol_ef { continue; }
 task.hits = compute_ef_hits(&self.ds, task.nE, task.nF, &ef_corr);
 }
 //
 // ------------------------------------------------------------------
 // Phase 3: Process results 閳?VERTEX/EDGE dispatch (OCCT L324-571)
 // ------------------------------------------------------------------
 for k in 0..a_v_edge_face.len() {
 let task = &a_v_edge_face[k];
 let nE = task.nE; let nF = task.nF;
 // OCCT L367-371: PB range, indices, splittable, shrunk range
 // OCCT L382: aR1(aT1, aTS1), aR2(aTS2, aT2)
 let r1 = [task.aT1.min(task.aTS1), task.aT1.max(task.aTS1)];
 let r2 = [task.aTS2.min(task.aT2), task.aTS2.max(task.aT2)];
 // OCCT L384-386: FaceInfo On/In
 let face_vin: std::collections::HashSet<usize> =
 self.ds.face_info(nF).vertices_in.iter().copied().collect();
 let face_von: std::collections::HashSet<usize> =
 self.ds.faces[nF].face_info.vertices_on.iter().copied().collect();
 // OCCT L388-394: bLinePlane
 let b_line_plane = matches!(
 (&self.ds.edges[nE].curve, &self.ds.faces[nF].surface),
 (Curve3::Line(_), Surface3::Plane(_))
 );
 //
 if task.hits.is_empty() { continue; }
 let tol_ef = self.ef_tol(nE, nF);
 //
 // OCCT L373-380: ReduceIntersectionRange 閳?if PB endpoints come from EE
 // intersection with face edges, clip the shrunk range aTS1/aTS2.
 // This shrinks aR1/aR2, making more near-endpoint intersections be
 // treated as valid EF (not erroneously skipped by IsInRange).
 let mut a_ts1_adj = r1[1]; // aTS1 from the r1 range
 let mut a_ts2_adj = r2[0]; // aTS2 from the r2 range
 {
 let nV1 = task.nV1;
 let nV2 = task.nV2;
 // OCCT L692-695: check if either vertex is a new shape
 let is_v1_new = nV1 < self.ds.vertices.len() && self.ds.vertex_origin(nV1).is_none();
 let is_v2_new = nV2 < self.ds.vertices.len() && self.ds.vertex_origin(nV2).is_none();
 if (is_v1_new || is_v2_new) && !self.ds.interf_ee.is_empty() {
 // OCCT L712-723: collect face's boundary edges
 let face_edges: std::collections::HashSet<usize> =
 self.ds.face_boundary_edges(nF).iter().copied().collect();
 // OCCT L725-767: iterate EE interferences
 for inf in &self.ds.interf_ee {
 // OCCT L728-731: check if EE has a new vertex
 let new_v = inf.new_vertex;
 if !is_v1_new && new_v != nV1 { continue; }
 if !is_v2_new && new_v != nV2 { continue; }
 if new_v != nV1 && new_v != nV2 { continue; }
 // OCCT L742-746: check EE involves our edge nE AND a face edge
 let involves_our_edge = inf.e1 == nE || inf.e2 == nE;
 let involves_face_edge = face_edges.contains(&inf.e1) || face_edges.contains(&inf.e2);
 if !involves_our_edge || !involves_face_edge { continue; }
 // OCCT L749-766: clip shrunk range by EE intersection parameter
 let ee_param = if inf.e1 == nE { inf.param1 } else { inf.param2 };
 if new_v == nV1 && a_ts1_adj < ee_param { a_ts1_adj = ee_param; }
 if new_v == nV2 && a_ts2_adj > ee_param { a_ts2_adj = ee_param; }
 }
 }
 }
 // Rebuild r1/r2 with adjusted shrunk range
 let r1_adj = [task.aT1.min(a_ts1_adj), task.aT1.max(a_ts1_adj)];
 let r2_adj = [a_ts2_adj.min(task.aT2), a_ts2_adj.max(task.aT2)];
 // OCCT L396-570: for each common part (hit)
 for &hit in &task.hits {
 match hit {
 EfHit::Vertex { point, param: a_t } => {
 // OCCT L406-543: case TopAbs_VERTEX
 // OCCT L412: VertexParameter 閳?a_t is the edge parameter
 // OCCT L415-419: IsInRange
 let a_tol_to_decide = 5e-8;
 let b_is_on_pave0 = (a_t - r1_adj[0]).abs() <= a_tol_to_decide
 || (a_t - r1_adj[1]).abs() <= a_tol_to_decide;
 let b_is_on_pave1 = (a_t - r2_adj[0]).abs() <= a_tol_to_decide
 || (a_t - r2_adj[1]).abs() <= a_tol_to_decide;
 // OCCT L421-439: if near both 閳?EDGE type
 if (b_is_on_pave0 && b_is_on_pave1) || (b_line_plane && (b_is_on_pave0 || b_is_on_pave1)) {
 let (nV1, nV2) = {
 let pb = &self.ds.edge_pave_blocks(nE)[0];
 pb.0.read().unwrap().indices()
 };
 let bV0 = face_von.contains(&nV1) || face_vin.contains(&nV1);
 let bV1 = face_von.contains(&nV2) || face_vin.contains(&nV2);
 if bV0 && bV1 {
 // OCCT L427-437: Create EF + aMIEFC
 self.ds.interf_ef.push(InterferenceEF{
 edge: nE, face: nF, point, edge_param: a_t, new_vertex: 0,
 });
 a_mi_efc.insert(nF); continue;
 }
 }
 // OCCT L442-444: splittable check
 let is_splittable = {
 let pb = &self.ds.edge_pave_blocks(nE)[0];
 pb.0.read().unwrap().is_splittable
 };
 if !is_splittable { continue; }
 // OCCT L447-457: ForceInterfVF
 let mut b_is_on_pave = [b_is_on_pave0, b_is_on_pave1];
 for j in 0..2 {
 if b_is_on_pave[j] {
 let nVx = if j == 0 {
 self.ds.edge_pave_blocks(nE)[0].0.read().unwrap().pave1.vertex_idx
 } else {
 self.ds.edge_pave_blocks(nE)[0].0.read().unwrap().pave2.vertex_idx
 };
 let bV_on_face = face_von.contains(&nVx) || face_vin.contains(&nVx);
 if !bV_on_face {
 // simplified ForceInterfVF
 if !self.ds.has_interf_vf(nVx, nF) {
 self.ds.interf_vf.push(InterferenceVF{
 vertex: nVx, face: nF, u: 0.0, v: 0.0,
 });
 }
 b_is_on_pave[j] = true;
 }
 }
 }
 // OCCT L459-502: if (bIsOnPave[0] || bIsOnPave[1])
 if b_is_on_pave[0] || b_is_on_pave[1] {
 // OCCT L467-472: hasRealIntersection check 閳?project point onto face
 // rcad: simplified (no projection), trust the hit
 // OCCT L482-501: UpdateVertex for near-endpoint vertex
 for j in 0..2 {
 if b_is_on_pave[j] {
 let nVx = if j == 0 {
 self.ds.edge_pave_blocks(nE)[0].0.read().unwrap().pave1.vertex_idx
 } else {
 self.ds.edge_pave_blocks(nE)[0].0.read().unwrap().pave2.vertex_idx
 };
 // OCCT L486-489: get existing vertex, compute distance
 let dist_pp = (self.ds.vertex_point(nVx) - point).length();
 let v_tol = self.ds.vertex_tolerance(nVx);
 // OCCT L490-494: aMaxDist = 1e4 * aTol; capped at 0.1 if tol < 0.01
 let max_dist = (1e4 * v_tol).min(if v_tol < 0.01 { 0.1 } else { f64::MAX });
 if dist_pp < max_dist {
 // OCCT L495-497: if (aDistPP < aMaxDist) UpdateVertex(nV[j], aDistPP)
 if dist_pp > self.ds.vertex_tolerance(nVx) {
 self.ds.vertex_data_mut(nVx).tolerance = dist_pp;
 self.ds.increased_ss.insert(nVx);
 }
 // OCCT L498: myVertsToAvoidExtension.Add(nV[j])
 // rcad: not tracked
 }
 }
 }
 continue;
 }
 // OCCT L505-508: CheckFacePaves(aVnew, aMIFOn)
 let near_face_vx = face_von.iter().any(|&vi| {
 vi < self.ds.vertices.len()
 && (self.ds.vertex_point(vi) - point).length() <= tol_ef
 }) || face_vin.iter().any(|&vi| {
 vi < self.ds.vertices.len()
 && (self.ds.vertex_point(vi) - point).length() <= tol_ef
 });
 if near_face_vx { continue; }
 // OCCT L510-519: aTolVnew = max(aTolVnew, aTolE, aTolF)
 //   bLinePlane 閳?increase tolerance by (aCR.Last() - aCR.First()) / 2.
 let mut a_tol_vnew = tol_ef.max(self.ds.edge_tolerance(nE)).max(self.ds.face_tolerance(nF));
 if b_line_plane {
 // OCCT L517-518: aTolVnew = max(aTolVnew, (aCR.Last() - aCR.First()) / 2.)
 a_tol_vnew = a_tol_vnew.max(tol_ef * 10.0);
 }
 // OCCT L523-526: if (!myContext->IsPointInFace(aPnew, aF, aTolVnew)) continue;
 if !self.is_point_in_face(point, nF, a_tol_vnew) { continue; }
 // OCCT L528-542: Create EF interference
 let new_v = self.ds.add_vertex(point);
 self.ds.interf_ef.push(InterferenceEF{
 edge: nE, face: nF, point, edge_param: a_t, new_vertex: new_v,
 });
 self.ds.faces[nF].face_info.vertices_on.insert(new_v);
 self.ds.edge_paves[nE].push(Pave { vertex_idx: new_v, param: a_t });
 a_mi_efc.insert(nF);
 // OCCT L537-542: aMVCPB coupling 閳?couples new vertex with PB for PerformNewVertices
 //   rcad: implicit via InterfEF::new_vertex, consumed by treat_new_vertices()
 } // EfHit::Vertex
 EfHit::Edge { t1, t2 } => {
 // ================================================================
 // OCCT L545-565: case TopAbs_EDGE
 // ================================================================
 // OCCT L546: aMIEFC.Add(nF)
 a_mi_efc.insert(nF);
 // OCCT L549-551: BOPDS_InterfEF& aEF = aEFs.Appended(); SetIndices(nE, nF)
 // rcad: use midpoint param for the edge-on-face common part
 let mid_t = (t1 + t2) / 2.0;
 let mid_pt = self.ds.edges[nE].curve.point_at(mid_t);
 // OCCT L553-555: bV[0]=CheckFacePaves(nV[0],aMIFOn,aMIFIn); if !bV[0]||!bV[1]
 let (nV1, nV2) = (task.nV1, task.nV2);
 let bV0 = face_von.contains(&nV1) || face_vin.contains(&nV1);
 let bV1 = face_von.contains(&nV2) || face_vin.contains(&nV2);
 if !bV0 || !bV1 {
 // OCCT L557: myDS->AddInterf(nE, nF) 閳?just mark as interfering
 if !self.ds.has_interf_ef(nE, nF) {
 self.ds.interf_ef.push(InterferenceEF{
 edge: nE, face: nF, point: mid_pt, edge_param: mid_t, new_vertex: 0,
 });
 }
 } else {
 // OCCT L560-564: aEF.SetCommonPart + AddInterf + FillMap(aPB, nF, aMPBLI)
 if !self.ds.has_interf_ef(nE, nF) {
 self.ds.interf_ef.push(InterferenceEF{
 edge: nE, face: nF, point: mid_pt, edge_param: mid_t, new_vertex: 0,
 });
 }
 // rcad: FillMap (PB閳姁ace list for common blocks) 閳?not yet aligned
 }
 } // EfHit::Edge
 } // match hit
 } // for each hit
 } // for each task
 //
 // ------------------------------------------------------------------
 // Phase 4: Post-treatment (OCCT L576-592)
 // ------------------------------------------------------------------
 // OCCT L576: BOPAlgo_Tools::PerformCommonBlocks(aMPBLI, ...)
 //   rcad: calls perform_common_blocks which scans all edges for coincident PBs
 // OCCT L577: UpdateVerticesOfCB() 閳?handled inside perform_common_blocks
 // OCCT L578: PerformNewVertices(aMVCPB) 閳?handled by treat_new_vertices() in mod.rs
 // OCCT L585: myDS->UpdateFaceInfoIn(aMIEFC)
 if !a_v_edge_face.is_empty() {
 crate::bopds::tools::perform_common_blocks(&mut self.ds);
 }
 for &fi in &a_mi_efc {
 self.ds.update_face_info_in(fi);
 }
 for &fi in &a_mi_efc {
 self.ds.update_face_info_in(fi);
 }
 }
 /// OCCT BOPDS_Iterator: face BVH construction
 pub(crate) fn build_ds_bvh_face(&self, is_a: bool) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let (start, end) = if is_a {
 (0, self.ds.a_face_count)
 } else {
 (self.ds.a_face_count, self.ds.faces.len())
 };
 let n = end - start;
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);
 for local_i in 0..n {
 let fi = start + local_i;
 indices.push(fi);
 let f = &self.ds.faces[fi];
 let mut aabb = Aabb::empty();
 // Boundary vertices
 for &vi in &f.boundary_verts {
 if vi < self.ds.vertices.len() {
 aabb.expand_point(self.ds.vertex_point(vi));
 }
 }
 // OCCT BndLib_AddSurface: expand AABB for curved surfaces.
 // Sphere: full sphere AABB = center  ?radius (face boundary
 // vertices only cover a patch, not the whole sphere volume).
 // Cylinder/Cone: boundary vertices already span the full
 // parametric extent  ?no extra expansion needed.
 if let Surface3::Sphere(s) = &f.surface {
 let r = s.radius.abs();
 aabb.expand_point(s.center + DVec3::splat(r));
 aabb.expand_point(s.center - DVec3::splat(r));
 }
 let tol = f.geom_tol.max(CONFUSION);
 aabb.min -= DVec3::splat(tol);
 aabb.max += DVec3::splat(tol);
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 ///  BOPDS_Iterator  ?combined face BVH (both operands).
 pub(crate) fn build_ds_bvh_face_all(&self) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let n = self.ds.faces.len();
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);
 for fi in 0..n {
 indices.push(fi);
 let f = &self.ds.faces[fi];
 let mut aabb = Aabb::empty();
 for &vi in &f.boundary_verts {
 if vi < self.ds.vertices.len() {
 aabb.expand_point(self.ds.vertex_point(vi));
 }
 }
 if let Surface3::Sphere(s) = &f.surface {
 let r = s.radius.abs();
 aabb.expand_point(s.center + DVec3::splat(r));
 aabb.expand_point(s.center - DVec3::splat(r));
 }
 let tol = f.geom_tol.max(CONFUSION);
 aabb.min -= DVec3::splat(tol);
 aabb.max += DVec3::splat(tol);
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ve_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all vertices are shared, skip V-E pass
 let shared_verts = &self.ds.shared_topology.shared_vertices;
 if shared_verts.is_empty() {
 return false;
 }

 // Check if all vertices from shape A have matches in shape B
 let a_verts: std::collections::HashSet<usize> = self.ds.vertices
 .iter()
 .enumerate()
 .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeA))
 .map(|(i, _)| i)
 .collect();

 let matched_a: std::collections::HashSet<usize> = shared_verts
 .iter()
 .map(|(a, _)| *a)
 .collect();

 a_verts == matched_a && !a_verts.is_empty()
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ee_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 let shared_edges = &self.ds.shared_topology.shared_edges;
 if shared_edges.is_empty() {
 return false;
 }

 // Check if all edges from shape A have matches in shape B
 let a_edges: std::collections::HashSet<usize> = self.ds.edges
 .iter()
 .enumerate()
 .filter(|(_, e)| e.origin == ShapeOrigin::ShapeA)
 .map(|(i, _)| i)
 .collect();

 let matched_a: std::collections::HashSet<usize> = shared_edges
 .iter()
 .map(|(a, _)| *a)
 .collect();

 a_edges == matched_a && !a_edges.is_empty()
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_vf_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all faces are fully glued, skip V-F pass
 !self.ds.shared_topology.fully_glued_faces.is_empty()
 && self.ds.shared_topology.fully_glued_faces.len()
 == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ef_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all faces are fully glued, skip E-F pass
 !self.ds.shared_topology.fully_glued_faces.is_empty()
 && self.ds.shared_topology.fully_glued_faces.len()
 == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ff_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all faces are fully glued, skip F-F pass
 let total_face_pairs = self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count);
 self.ds.shared_topology.fully_glued_faces.len() == total_face_pairs && total_face_pairs > 0
 }
 /// 閴?PerformVV (PaveFiller_1.cxx L45-132).
 /// Builds vertex-vertex connection map (FillMap), groups connected
 /// vertices (MakeBlocks), then creates SD vertices for each group.
 /// Pairs come pre-computed from BOPDS_Iterator (cross-operand, AABB-filtered).
  pub(crate) fn perform_vv(&mut self, pairs: &[(usize, usize)]) {
  // OCCT L47-56: n1, n2, iFlag, aSize; iterator init + early return
  let a_vc = self.ds.a_vertex_count;
  let a_size = a_vc * (self.ds.vertices.len() - a_vc);
  if a_size == 0 { return; }

  // 閴?BOPDS_Iterator(VERTEX, VERTEX) 閳?BVH-based pair enumeration.
  // OCCT L68-76: myIterator->Initialize(VERTEX, VERTEX) returns overlapping AABB pairs.
  let mut a_mili: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

  // OCCT L68-98: 1. Map V/LV 閳?build connection map of close vertex pairs.
  for &(n1, n2) in pairs {
  // Skip same-operand pairs (OCCT: cross-operand only)
  if (n1 < a_vc) == (n2 < a_vc) { continue; }

  // OCCT L77-81: if HasInterf 閳?FillMap + continue
  if self.ds.has_interf_vv(n1, n2) {
  fill_map(&mut a_mili, n1, n2);
  continue;
  }

  // OCCT L84-91: Resolve SD vertices (HasShapeSD) + ComputeVV
  let n1sd = self.ds.has_shape_sd(n1).unwrap_or(n1);
  let n2sd = self.ds.has_shape_sd(n2).unwrap_or(n2);

  // OCCT L93: ComputeVV(aV1, aV2, myFuzzyValue) 閳?tolerance-based distance check
  let tol = self.vv_pair_tol(n1, n2);
  let dist = (self.ds.vertex_point(n1sd) - self.ds.vertex_point(n2sd)).length();
  let i_flag = if dist <= tol { 0 } else { 1 };

  // OCCT L94-97: if !iFlag (vertices interfere) 閳?FillMap
  if i_flag == 0 {
  fill_map(&mut a_mili, n1, n2);
  }
  }
 let a_m_blocks = make_blocks(&a_mili);
 for block in &a_m_blocks {
 if block.len() < 2 { continue; }
 self.make_sd_vertices_vv(block);
 }
  // OCCT ShapesSD is a DataMap<source, target> (one direction).
  // rcad stores (source, target) bidirectionally; dedup via HashSet.
  let a_dmii: std::collections::HashSet<usize> =
    self.ds.shape_sd.sd_vertices_iter().map(|&(k, _)| k).collect();
  for &n1 in &a_dmii {
    self.ds.init_pave_blocks_for_vertex(n1);
  }
  }

 ///  MakeSDVertices (PaveFiller_1.cxx L136-233).
 /// Merges a connected group of vertices into a single SD vertex.
 /// The first vertex in the block becomes the merge target; all other
 /// vertices in the block are remapped to it via AddShapeSD and
 /// Interference::VertexVertex entries.
 pub(super) fn make_sd_vertices_vv(&mut self, block: &[usize]) {
 if block.len() < 2 { return; }
 // nSD tracks the existing SD vertex index; others go into aLV.
 let mut n_sd: Option<usize> = None;
 let mut a_lv: Vec<usize> = Vec::with_capacity(block.len());
 for &n_x in block {
 if let Some(n_sd1) = self.ds.has_shape_sd(n_x) {
 if n_sd.is_none() {
 n_sd = Some(n_sd1);
 }
 }
 a_lv.push(n_x);
 }
 // rcad: use the minimum-index vertex as the merged target.
 let n_v = n_sd.unwrap_or_else(|| *a_lv.iter().min().unwrap());
 // rcad: geom_tol of the merged vertex is max of all members' tolerances.
 if let Some(&target) = block.iter().max_by(|&&a, &&b| {
 self.ds.vertex_tolerance(a).partial_cmp(&self.ds.vertex_tolerance(b)).unwrap()
 }) {
 if n_v < self.ds.vertices.len() {
 self.ds.vertex_data_mut(n_v).tolerance = self.ds.vertex_tolerance(n_v)
 .max(self.ds.vertex_tolerance(target));
 }
 }
 for i in 0..block.len() {
 let n1 = block[i];
 self.ds.add_shape_sd(n1, n_v);
 // rcad: ShapeOrigin check.
 // (OCCT L208-218: self-interfering shape warning  ?skipped for brevity)
 for j in (i + 1)..block.len() {
 let n2 = block[j];
 // rcad: push VertexVertex interference
 self.ds.interf_vv.push(InterferenceVV{
 v1: n1,
 v2: n2,
 merged_vertex: n_v,
 });
 }
 }
 }
 /// OCCT PaveFiller_2.cxx L141-206: PerformVE
 pub(crate) fn perform_ve(&mut self) {
 // with HasSubShape / HasFlag / HasInterf / HasInterfShapeSubShapes skips.
 self.fill_shrunk_data();
 let a_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeA);
 let b_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeB);
 // shrink data is computed on-the-fly in compute_ve via ve_tol().
 //
 // rcad: manual O(n ? loop (see PairIterator in perform_ee for BVH pattern).

 for &vi in &a_verts {
 for &ei in &b_edges {
 if self.ds.edge_has_vertex(vi, ei) { continue; }
 if self.ds.edge_has_flag(ei) { continue; }
 if self.ds.has_interf_ve(vi, ei) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, ei) { continue; }
 if self.ds.is_edge_degenerated(ei) { continue; }
 if self.ds.edge_pave_blocks(ei).is_empty() { continue; }
 if !self.ds.edge_pave_blocks(ei)[0].0.read().unwrap().is_splittable { continue; }
 self.compute_ve(vi, ei);
 }
 }

 let b_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeB);
 let a_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeA);

 for &vi in &b_verts {
 for &ei in &a_edges {
 if self.ds.edge_has_vertex(vi, ei) { continue; }
 if self.ds.edge_has_flag(ei) { continue; }
 if self.ds.has_interf_ve(vi, ei) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, ei) { continue; }
 if self.ds.is_edge_degenerated(ei) { continue; }
 if self.ds.edge_pave_blocks(ei).is_empty() { continue; }
 if !self.ds.edge_pave_blocks(ei)[0].0.read().unwrap().is_splittable { continue; }
 self.compute_ve(vi, ei);
 }
 }
 }
 /// OCCT PaveFiller_2.cxx L104-121: ComputeVE
 pub(crate) fn compute_ve(&mut self, vi: usize, ei: usize) {
 let fuzz = self.fuzzy_tolerance;
 if let Ok(res) = self.context.compute_ve(self.ds, vi, ei, fuzz) {
  self.ds.interf_ve.push(InterferenceVE{vertex: vi, edge: ei, param: res.param});
  self.ds.edge_paves[ei].push(Pave {vertex_idx: vi, param: res.param});
 }
 }
 /// OCCT PaveFiller_3.cxx L145-244: PerformEE
 pub(crate) fn perform_ee(&mut self) {
 // with HasFlag / PaveBlock emptiness / GetPBBox skip conditions.
 self.fill_shrunk_data();
 let a_count = self.ds.a_edge_count;

 // Build a set of shared edge pairs for fast lookup when glue is enabled
 let shared_edge_set: std::collections::HashSet<(usize, usize)> = if self.use_glue() {
 self.ds
 .shared_topology
 .shared_edges
 .iter()
 .map(|(e1, e2)| (*e1, *e2))
 .collect()
 } else {
 std::collections::HashSet::new()
 };
 // rcad: PairIterator for cross-group pairs (A-edges  ?B-edges).
 // For PaveBlock-level precision (OCCT L200-232), iterate sub-ranges
 // of each edge defined by existing paves (from VE or prior intersections).
 // Each sub-range = one logical PaveBlock.
 let mut it = crate::bopds::ds::PairIterator::prepare_ab(a_count, self.ds.edges.len());
 while it.more() {
 let pk = it.value();
 let ae = pk.i1; let be = pk.i2;
 if self.ds.edge_has_flag(ae) || self.ds.edge_has_flag(be) {
 it.next(); continue;
 }
 if self.ds.has_interf_ee(ae, be) {
 it.next(); continue;
 }

 if self.ds.is_edge_degenerated(ae) || self.ds.is_edge_degenerated(be) {
 it.next(); continue;
 }
 // PaveBlocks of each edge).  rcad: build sub-ranges from existing
 // paves to limit intersection to relevant sub-segments.
 let ranges_a = self.collect_paveblock_ranges(ae, self.ds.edge_range(ae));
 let ranges_b = self.collect_paveblock_ranges(be, self.ds.edge_range(be));

 if ranges_a.is_empty() || ranges_b.is_empty() {
 it.next(); continue;
 }

 if self.use_glue() && shared_edge_set.contains(&(ae, be)) {
 // Glue: use first pave point as shared vertex
 let pv = self.ds.edge_start_vertex_ds(ae);
 if !self.ds.has_interf_ee(ae, be) {
 self.ds.interf_ee.push(InterferenceEE{
 e1: ae, e2: be,
 point: self.ds.vertex_point(pv),
 param1: self.ds.edge_range(ae)[0],
 param2: self.ds.edge_range(be)[0],
 new_vertex: pv,
 });
 }
 } else {
 let mut _ee_modified: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for ra in &ranges_a {
 for rb in &ranges_b {
 self.intersect_ee(ae, be, *ra, *rb, &mut _ee_modified);
 }
 }
 }
 it.next();
 }
 }
 /// OCCT PaveFiller_3.cxx L580-640: CheckEdgeEdge
 pub(crate) fn intersect_ee(&mut self, e1: usize, e2: usize,
 range1: [f64; 2], range2: [f64; 2],
 modified: &mut std::collections::HashSet<usize>) {
 let edge1 = &self.ds.edges[e1];
 let edge2 = &self.ds.edges[e2];
 let tol = self.ee_tol(e1, e2);
 // Capture all edge data before mutable borrow
 let e1_curve = edge1.curve.clone();
 let e2_curve = edge2.curve.clone();
 drop(edge1);
 drop(edge2);
 // 閴?FillShrunkData computes shrunk ranges for each pave block.
 // If shrunk_range fails (edge too short), skip this pair entirely
 // (=OCCT BOPAlgo_PaveFiller_3: !aPB->IsSplittable() 閳?continue).
 let sr1 = match crate::inttools::curve_range::shrunk_range(
  &e1_curve, range1, tol, tol, tol) {
  Some(sr) => sr,
  None => return,
 };
 let sr2 = match crate::inttools::curve_range::shrunk_range(
  &e2_curve, range2, tol, tol, tol) {
  Some(sr) => sr,
  None => return,
 };

 // Compute intersections restricted to shrunk sub-ranges.
 let hits: Vec<(f64, f64, DVec3)> = match (&e1_curve, &e2_curve) {
 (Curve3::Line(l1), Curve3::Line(l2)) => {
 intersect_line_line(l1, sr1, l2, sr2, tol)
 .into_iter().map(|(t1, t2, p)| (t1, t2, p)).collect()
 }
 (Curve3::Line(l), Curve3::Circle(c)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, sr1, tol) && in_range(*t_circle, sr2, tol)
 })
 .map(|(t_line, t_circle, p)| (t_line, t_circle, p))
 .collect(),
 (Curve3::Circle(c), Curve3::Line(l)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, sr2, tol) && in_range(*t_circle, sr1, tol)
 })
 .map(|(t_line, t_circle, p)| (t_circle, t_line, p))
 .collect(),
  (Curve3::Circle(c1), Curve3::Circle(c2)) => intersect_circle_circle(c1, c2, tol)
  .into_iter()
  .filter_map(|p| {
  let t1 = circle_param(p, c1);
  let t2 = circle_param(p, c2);
  if in_range(t1, sr1, tol) && in_range(t2, sr2, tol) {
  Some((t1, t2, p))
  } else { None }
  })
  .collect(),
  _ => {
  // Fallback: use EdgeEdgeIntersector for non-analytic curve pairs
  let mut ee = crate::inttools::edge_edge::EdgeEdgeIntersector::new();
  ee.set_edges(e1, sr1, e2, sr2, self.ds);
  ee.set_fuzzy_value(tol);
  ee.perform();
  ee.common_parts().iter().map(|cp| {
  let t1 = cp.vertex_param1;
  let t2 = cp.vertex_param2;
  Some((t1, t2, cp.bounding_point1))
  }).flatten().collect()
  }
 };

 //  Process each intersection result (PaveFiller_3.cxx L682-750).
 // For each valid intersection, create a new vertex and record EE interference.
 // OCCT's UpdateVertex handles proximity via tolerance merging; rcad creates
 // vertices directly (architecture diff: rcad DSVertex has no UpdateVertex).
 for (t1, t2, point) in hits {
 // 閴?restrict to shrunk range.  IntTools_EdgeEdge computes
 // within the shrunk range; results at/outside the boundary are
 // endpoint-coincident (handled by VV/VE/VF) or coincide with an existing
 // pave vertex 閳?neither should create a new EE interference.
 if t1 < sr1[0] || t1 > sr1[1] || t2 < sr2[0] || t2 > sr2[1] { continue; }
 // 閴?skip tangent/colinear edge pairs.  OCCT
 // (PaveFiller_3.cxx) checks aEECP.TangentEdges() after EE computation
 // and defers the entire range to ForceInterfEE (CommonBlocks).
 let is_parallel = match (&e1_curve, &e2_curve) {
  (Curve3::Line(l1), Curve3::Line(l2)) =>
   l1.direction.cross(l2.direction).length() <= tol,
  _ => false,
 };
 if is_parallel { continue; }
 let new_v = self.ds.add_vertex(point);
 let new_v = self.ds.add_vertex(point);
 self.ds.interf_ee.push(InterferenceEE{
 e1, e2, point, param1: t1, param2: t2, new_vertex: new_v,
 });
 self.ds.edge_paves[e1].push(Pave { vertex_idx: new_v, param: t1 });
 self.ds.edge_paves[e2].push(Pave { vertex_idx: new_v, param: t2 });
 modified.insert(e1);
 modified.insert(e2);
 }
 }
 /// OCCT PaveFiller_3.cxx L580-640: CheckEdgeEdge
 pub(crate) fn check_edge_edge(&mut self, e1: usize, e2: usize) {
 let edge1 = &self.ds.edges[e1];
 let edge2 = &self.ds.edges[e2];
 let tol = self.ee_tol(e1, e2);

 let hits: Vec<(f64, f64, DVec3)> = match (&edge1.curve, &edge2.curve) {
 (Curve3::Line(l1), Curve3::Line(l2)) => {
 intersect_line_line(l1, edge1.t_range, l2, edge2.t_range, tol)
 .into_iter()
 .map(|(t1, t2, p)| (t1, t2, p))
 .collect()
 }
 (Curve3::Line(l), Curve3::Circle(c)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, edge1.t_range, tol)
 && in_range(*t_circle, edge2.t_range, tol)
 })
 .map(|(t_line, t_circle, p)| (t_line, t_circle, p))
 .collect(),
 (Curve3::Circle(c), Curve3::Line(l)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, edge2.t_range, tol)
 && in_range(*t_circle, edge1.t_range, tol)
 })
 .map(|(t_line, t_circle, p)| (t_circle, t_line, p))
 .collect(),
 (Curve3::Circle(c1), Curve3::Circle(c2)) => intersect_circle_circle(c1, c2, tol)
 .into_iter()
 .filter_map(|p| {
 let t1 = circle_param(p, c1);
 let t2 = circle_param(p, c2);
 if in_range(t1, edge1.t_range, tol) && in_range(t2, edge2.t_range, tol) {
 Some((t1, t2, p))
 } else {
 None
 }
 })
 .collect(),
 _ => vec![],
 };

 for (t1, t2, point) in hits {
 let new_v = self.ds.add_vertex(point);
 self.ds.interf_ee.push(InterferenceEE{
 e1,
 e2,
 point,
 param1: t1,
 param2: t2,
 new_vertex: new_v,
 });
 self.ds.edge_paves[e1].push(Pave {
 vertex_idx: new_v,
 param: t1,
 });
 self.ds.edge_paves[e2].push(Pave {
 vertex_idx: new_v,
 param: t2,
 });
 }
 }
 /// OCCT PaveFiller L575-590: TreatNewVertices (merge EE/EF new vertices)
 pub(crate) fn treat_new_vertices(&mut self) -> Vec<usize> {
 // = =  Phase 1: Collect new vertices (OCCT L696-702) = = = = = = = = = = = = = = = 
 #[derive(Clone, Copy)]
 struct NewVertInfo { idx: usize, pos: DVec3, tol: f64 }
 let mut new_verts: Vec<NewVertInfo> = Vec::new();
 let mut seen = std::collections::BTreeSet::new();
 for inf in &self.ds.interf_ee {
 if seen.insert(inf.new_vertex) {
 let v_tol = self.ds.vertices[inf.new_vertex].geom_tol.max(self.ds.fuzzy_tol);
 new_verts.push(NewVertInfo { idx: inf.new_vertex, pos: inf.point, tol: v_tol });
 }
 }
 for inf in &self.ds.interf_ef {
 if seen.insert(inf.new_vertex) {
 let v_tol = self.ds.vertices[inf.new_vertex].geom_tol.max(self.ds.fuzzy_tol);
 new_verts.push(NewVertInfo { idx: inf.new_vertex, pos: inf.point, tol: v_tol });
 }
 }
 if new_verts.len() < 2 { return vec![]; }

 // = =  Phase 2: IntersectVertices (BOPAlgo_Tools.hxx L1119-1205) = = =
 let gap = self.ds.fuzzy_tol / 2.0;

 use crate::bvh::{Aabb, DsBvh};
 let nv = new_verts.len();
 let mut bvh_indices: Vec<usize> = Vec::with_capacity(nv);
 let mut bvh_aabbs: Vec<Aabb> = Vec::with_capacity(nv);
 for i in 0..nv {
 let v = &new_verts[i];
 bvh_indices.push(i);
 let half = v.tol + gap;
 bvh_aabbs.push(Aabb {
 min: v.pos - DVec3::splat(half),
 max: v.pos + DVec3::splat(half),
 });
 }
 let bvh = DsBvh::build(bvh_indices, bvh_aabbs);

 let pairs = DsBvh::candidate_pairs(&bvh, &bvh);
 let mut a_mili: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
 for &(ia, ib) in &pairs {
   let geom_vi = self.ds.vertices[new_verts[ia].idx].geom_tol;
   let geom_vj = self.ds.vertices[new_verts[ib].idx].geom_tol;
   let merge_tol = geom_vi + geom_vj + self.ds.fuzzy_tol;
   if (new_verts[ia].pos - new_verts[ib].pos).length() <= merge_tol {
     fill_map(&mut a_mili, ia, ib);
   }
 }
 let a_blocks = make_blocks(&a_mili);
 let mut groups: Vec<Vec<usize>> = a_blocks.iter()
   .map(|block| block.iter().map(|&i| new_verts[i].idx).collect())
   .collect();
 {
   let mut taken: HashSet<usize> = HashSet::new();
   for group in &groups {
     for &vi in group {
       taken.insert(vi);
     }
   }
   for i in 0..nv {
     if !taken.contains(&new_verts[i].idx) {
       groups.push(vec![new_verts[i].idx]);
     }
   }
 }

 // = =  Phase 3: MakeVertex for each chain (BOPTools_AlgoTools::MakeVertex) = = =
 // Single element -> reuse the vertex.
 // Multiple elements -> BRepLib::BoundingVertex computes center + tolerance.
 let mut survivors: Vec<usize> = Vec::new();
 for members in &groups {
 if members.len() < 2 {
 survivors.push(members[0]);
 continue;
 }

 // OCCT: BRepLib::BoundingVertex computes center + tolerance.
 // rcad: compute centroid of all member vertex positions.
 let centroid = members.iter()
 .map(|&vi| self.ds.vertex_point(vi))
 .sum::<DVec3>() / members.len() as f64;
 let max_tol = members.iter()
 .map(|&vi| self.ds.vertex_tolerance(vi))
 .max_by(|a, b| a.partial_cmp(b).unwrap())
 .unwrap_or(self.ds.fuzzy_tol);

 // OCCT BRep_Builder::MakeVertex: create new vertex at centroid.
 let new_vi = self.ds.add_vertex(centroid);
 self.ds.vertex_data_mut(new_vi).tolerance = max_tol;
 // OCCT: myIncreasedSS.Add(nV)  ?mark tolerance as increased.
 self.ds.increased_ss.insert(new_vi);

 // Update interferences and paves to point to the new vertex
 for &old_vi in members {
 if old_vi == new_vi { continue; }
 for edge in &mut self.ds.edges {
 for pave in &mut edge.paves {
 if pave.vertex_idx == old_vi { pave.vertex_idx = new_vi; }
 }
 }
 for inf in &mut self.ds.interf_ee {
 if inf.new_vertex == old_vi { inf.new_vertex = new_vi; }
 }
 for inf in &mut self.ds.interf_ef {
 if inf.new_vertex == old_vi { inf.new_vertex = new_vi; }
 }
 for face in &mut self.ds.faces {
 if face.face_info.vertices_on.remove(&old_vi) {
 face.face_info.vertices_on.insert(new_vi);
 }
 if face.face_info.vertices_in.remove(&old_vi) {
 face.face_info.vertices_in.insert(new_vi);
 }
 }
 }
 survivors.push(new_vi);
 }
 survivors
 }
 /// OCCT PaveFiller_5.cxx L359-420: RepeatIntersection
 pub(crate) fn repeat_intersection(&mut self) {
 if self.ds.increased_ss.is_empty() { return; }
 let candidates: Vec<usize> = self.ds.increased_ss.iter().copied().collect();

 // Build set of existing interferences for dedup
 //  OCCT L398-413: PerformVV  ?PerformVE  ?PerformVF
 use std::collections::BTreeSet;
 let mut ve_done: BTreeSet<(usize, usize)> = BTreeSet::new();
 let mut vf_done: BTreeSet<(usize, usize)> = BTreeSet::new();
 for inf in &self.ds.interf_ve {
 ve_done.insert((inf.vertex, inf.edge));
 }
 for inf in &self.ds.interf_vf {
 vf_done.insert((inf.vertex, inf.face));
 }

 // = =  VV: check survivors against vertices on the other side = = = = = = 
 //  OCCT L398: PerformVV(aPS.Next())
 // VV safe: if pair already in interferences, add_vertex will dedup
 for &vi in &candidates {
 let vi_origin = self.ds.vertex_origin(vi);
 let other_verts: Vec<usize> = self.ds.vertices.iter().enumerate()
 .filter(|(j, v)| {
 if *j == vi { return false; }
 match (vi_origin, v.origin) {
 (Some(ShapeOrigin::ShapeA), Some(ShapeOrigin::ShapeB)) => true,
 (Some(ShapeOrigin::ShapeB), Some(ShapeOrigin::ShapeA)) => true,
 _ => false,
 }
 })
 .map(|(j, _)| j)
 .collect();
 for &vj in &other_verts {
 let tol = self.vv_pair_tol(vi, vj);
 let dist = (self.ds.vertex_point(vi) - self.ds.vertex_point(vj)).length();
 if dist <= tol {
 self.ds.interf_vv.push(InterferenceVV{
 v1: vi, v2: vj, merged_vertex: vi,
 });
 }
 }
 }

 // = =  VE: check survivors against edges on the other side = = = = = = 
 //  OCCT L403: PerformVE(aPS.Next())
 for &vi in &candidates {
 let vi_origin = self.ds.vertex_origin(vi);
 let other_edges: Vec<usize> = match vi_origin {
 Some(ShapeOrigin::ShapeA) => self.edges_of(ShapeOrigin::ShapeB),
 Some(ShapeOrigin::ShapeB) => self.edges_of(ShapeOrigin::ShapeA),
 _ => continue,
 };
 for &ei in &other_edges {
 if ve_done.contains(&(vi, ei)) { continue; }
 self.compute_ve(vi, ei);
 }
 }

 // = =  VF: check survivors against faces on the other side = = = = = = 
 //  OCCT L408: PerformVF(aPS.Next())
 for &vi in &candidates {
 let vi_origin = self.ds.vertex_origin(vi);
 let other_faces: Vec<usize> = match vi_origin {
 Some(ShapeOrigin::ShapeA) => self.faces_of(ShapeOrigin::ShapeB),
 Some(ShapeOrigin::ShapeB) => self.faces_of(ShapeOrigin::ShapeA),
 _ => continue,
 };
 for &fi in &other_faces {
 if vf_done.contains(&(vi, fi)) { continue; }
 self.check_vertex_face(vi, fi);
 }
 }
 }
 /// OCCT PaveFiller_4.cxx: PerformVF
 pub(crate) fn perform_vf(&mut self) {
 // OCCT PaveFiller_4.cxx: FillShrunkData + BVH pair iteration
 // with HasInterf skip condition.
 self.fill_shrunk_data(); // OCCT: FillShrunkData(VERTEX, FACE)
 let a_verts = self.verts_of(ShapeOrigin::ShapeA);
 let b_faces = self.faces_of(ShapeOrigin::ShapeB);
 for &vi in &a_verts {
 for &fi in &b_faces {
 if self.ds.has_interf_vf(vi, fi) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, fi) { continue; }
 self.check_vertex_face(vi, fi);
 }
 }
 let b_verts = self.verts_of(ShapeOrigin::ShapeB);
 let a_faces = self.faces_of(ShapeOrigin::ShapeA);
 for &vi in &b_verts {
 for &fi in &a_faces {
 if self.ds.has_interf_vf(vi, fi) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, fi) { continue; }
 self.check_vertex_face(vi, fi);
 }
 }
 }

 ///  CheckVertexFace (PaveFiller_4.cxx L249-298).
 /// Vertex/Face proximity check with SD vertex resolution.
 /// OCCT: BOPAlgo_VertexFace parallel solver + result processing;
 /// rcad: sequential equivalent with same projection logic.
 pub(crate) fn check_vertex_face(&mut self, vi: usize, fi: usize) {
 let n_vsd = self.ds.has_shape_sd(vi).unwrap_or(vi);
 let point = self.ds.vertex_point(n_vsd);
 let face = &self.ds.faces[fi];
 let tf = self.vf_tol(n_vsd, fi);
    let (is_on, proj_dist, proj_u, proj_v): (bool, f64, f64, f64) = match &face.surface {
        Surface3::Plane(plane) => {
            if inttools::vertex_ops::vertex_on_plane_with_tol(point, plane, tf) {
                let face_verts = self.ds.face_boundary_points(fi);
                let on_face = inttools::edge_face::point_in_planar_face_with_tol(
                    point, plane, &face_verts, tf);
                // UV: project point onto plane UV coordinate system
                let u_axis = plane.u_dir;
                let v_axis = plane.v_dir;
                let diff = point - plane.origin;
                let u = diff.dot(u_axis);
                let v = diff.dot(v_axis);
                (on_face, 0.0, u, v)
            } else {
                (false, f64::MAX, 0.0, 0.0)
            }
        }
        surface => {
            let proj = rcad_kernel::projection::closest_point_on_surface(
                surface, point, 16);
            let a_tol_v = self.ds.vertex_tolerance(n_vsd);
            let a_tol_f = face.geom_tol;
            let a_tol_sum = a_tol_v + a_tol_f + self.ds.fuzzy_tol.max(tf);
            if proj.distance <= a_tol_sum {
                let uv = DVec2::new(proj.params.0, proj.params.1);
                let fclass = FClass2d::new(self.ds, fi, tf);
                let inside = fclass.perform(uv, false) == State::In;
                (inside, proj.distance, proj.params.0, proj.params.1)
            } else {
                (false, f64::MAX, 0.0, 0.0)
            }
        }
    };

    if is_on {
        self.ds.interf_vf.push(InterferenceVF{
            vertex: n_vsd,
            face: fi,
            u: proj_u,
            v: proj_v,
        });
 if proj_dist > 0.0 && proj_dist < f64::MAX
 && proj_dist > self.ds.vertex_tolerance(n_vsd)
 {
 self.ds.vertex_data_mut(n_vsd).tolerance = proj_dist;
 self.ds.increased_ss.insert(n_vsd);
 }

 //  ALL VF vertices go to VerticesIn (OCCT L297: aMVIn.Add)
 self.ds.face_info_mut(fi).vertices_in.insert(n_vsd);
 }
 }
 /// PerformEF (PaveFiller_5.cxx L165-300) 閳?LEGACY non-BVH path.
 /// Only kept for reference; use perform_ef() (the aligned BVH path) instead.
 /// OCCT PaveFiller_3.cxx L222-228: GetPBBox (PaveBlock range)
 pub(crate) fn collect_paveblock_ranges(&self, edge_idx: usize, edge_t_range: [f64; 2]) -> Vec<[f64; 2]> {
 let paves = &self.ds.edge_paves(edge_idx);
 if paves.is_empty() {
 return vec![edge_t_range];
 }
 let mut params: Vec<f64> = paves.iter().map(|p| p.param).filter(|p| p.is_finite()).collect();
 params.sort_by(|a, b| a.partial_cmp(b).unwrap());
 let edge_tol = self.ds.edge_tolerance(edge_idx).max(self.tol());
 // Deduplicate
 params.dedup_by(|a, b| (*a - *b).abs() < edge_tol);
 // Include endpoints
 let mut bounds = vec![edge_t_range[0]];
 bounds.extend(params);
 bounds.push(edge_t_range[1]);
 // Build ranges
 let mut ranges = Vec::new();
 for w in bounds.windows(2) {
 if w[1] - w[0] > edge_tol {
 ranges.push([w[0], w[1]]);
 }
 }
 ranges
 }
 /// OCCT BOPTools_AlgoTools::CorrectRange (AlgoTools_2.cxx L364-434).
 /// Shrinks the range by the face tolerance converted to parametric space.
 pub(crate) fn correct_range_for_face(edge_curve: &Curve3, etf: f64, range: [f64; 2]) -> [f64; 2] {
 const DT: f64 = 1e-12;
 let a_tf = range[0];
 let a_tl = range[1];
 let mut a_new_first = a_tf;
 let mut a_new_last = a_tl;
 // OCCT L387-433: for (i = 0; i < 2; ++i)
 for i in 0..2 {
 let t = if i == 0 { a_tf } else { a_tl };
 // OCCT L389: aRes = aTolF; then convert to parametric space
 let a_res = match edge_curve {
 // OCCT L416-417: analytic 閳?aBC.Resolution(aRes)
 Curve3::Line(l) => {
 let dir_len = l.direction.length();
 if dir_len > 1e-12 { etf / dir_len } else { etf }
 }
 Curve3::Circle(c) => etf / c.radius.max(TOLERANCE_ABS),
 Curve3::Ellipse(e) => etf / e.major_radius.max(TOLERANCE_ABS),
 // OCCT L391-413: BSpline/Bezier 閳?aRes / |derivative|
 _ => {
 let dt = 1e-7;
 let p1 = edge_curve.point_at(t - dt);
 let p2 = edge_curve.point_at(t + dt);
 let dm = (p2 - p1).length() / (2.0 * dt);
 if dm > 1e-12 { etf / dm } else { etf }
 }
 };
 // OCCT L420-427: shrink endpoint
 if i == 0 { a_new_first = a_tf + a_res; }
 else { a_new_last = a_tl - a_res; }
 // OCCT L429-432: if too small, restore original
 if (a_new_last - a_new_first) < DT { return range; }
 }
 [a_new_first, a_new_last]
 }
 /// OCCT IntTools_FClass2d: point-in-face check
 pub(crate) fn is_point_in_face(&self, point: DVec3, face_idx: usize, tol: f64) -> bool {
 let face = &self.ds.faces[face_idx];
 match &face.surface {
 Surface3::Plane(plane) => {
 let verts = self.ds.face_boundary_points(face_idx);
 inttools::edge_face::point_in_planar_face_with_tol(point, plane, &verts, tol)
 }
 Surface3::Sphere(sphere) => {
 let uv = sphere.world_to_uv(point);
 uv.x >= -tol && uv.x <= std::f64::consts::TAU + tol
 && uv.y >= -std::f64::consts::FRAC_PI_2 - tol
 && uv.y <= std::f64::consts::FRAC_PI_2 + tol
 }
 Surface3::Cylinder(cyl) => {
 let axis = cyl.axis.normalize_or_zero();
 if axis.length_squared() < 0.5 { return true; }
 let local = point - cyl.origin;
 let v = local.dot(axis);
 let radial = local - axis * v;
 let u = radial.y.atan2(radial.x);
 let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
 u >= -tol && u <= std::f64::consts::TAU + tol
 }
 Surface3::Cone(cone) => {
 let axis = cone.axis_dir();
 let ap = point - cone.apex;
 let v = ap.dot(axis);
 let radial = ap - axis * v;
 let u = radial.y.atan2(radial.x);
 let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
 u >= -tol && u <= std::f64::consts::TAU + tol
 }
 Surface3::Torus(torus) => {
 let axis = torus.axis.normalize_or_zero();
 if axis.length_squared() < 0.5 { return true; }
 let local = point - torus.center;
 let v = local.dot(axis);
 let radial = local - axis * v;
 let r = radial.length();
 let u = radial.y.atan2(radial.x);
 let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
 let v_angle = ((r - torus.major_radius) / torus.minor_radius.max(tol)).acos();
 let v = if (r - torus.major_radius).abs() <= torus.minor_radius + tol { v_angle } else { 0.0 };
 u >= -tol && u <= std::f64::consts::TAU + tol
 && v >= -tol && v <= std::f64::consts::TAU + tol
 }
 _ => true,
 }
 }
 /// OCCT PaveFiller_5.cxx L340-480: IntersectEdgeFace
 pub(crate) fn intersect_ef(&mut self, edge_idx: usize, face_idx: usize, pb_range: &[f64; 2]) {
 let edge_curve = self.ds.edges[edge_idx].curve.clone();
 let edge_t_range = self.ds.edge_range(edge_idx);

 // Use PaveBlock range to constrain intersection interval (OCCT L262: SetRange(aPBRange))
 let ef_range = [
 pb_range[0].max(edge_t_range[0]),
 pb_range[1].min(edge_t_range[1]),
 ];
 let etf = self.ef_tol(edge_idx, face_idx);
 if ef_range[1] - ef_range[0] <= etf {
 return;
 }
 let ef_range = Self::correct_range_for_face(&edge_curve, etf, ef_range);
 if ef_range[1] - ef_range[0] <= etf {
 return;
 }
 let face_surface = self.ds.faces[face_idx].surface.clone();

 // Dispatch based on curve type  ?surface type
 let hits: Vec<(DVec3, f64)> = match (&edge_curve, &face_surface) {
 (Curve3::Line(line), Surface3::Plane(plane)) => {
 inttools::edge_face::intersect_line_plane_with_tol(
 line,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.edge_param))
 .collect()
 }
 (Curve3::Line(line), Surface3::Cylinder(cyl)) => {
 inttools::curve_surface::intersect_line_cylinder_with_tol(
 line,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Line(line), Surface3::Sphere(sph)) => {
 inttools::curve_surface::intersect_line_sphere_with_tol(
 line,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Line(line), Surface3::Cone(cone)) => {
 inttools::curve_surface::intersect_line_cone_with_tol(
 line,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Circle(circle), Surface3::Plane(plane)) => {
 // Use edge start vertex as reference direction for  ?0
 let sv = self.ds.edge_start_vertex_ds(edge_idx);
 let ref_dir = (self.ds.vertex_point(sv) - circle.center).normalize();
 inttools::curve_surface::intersect_circle_plane_with_ref(
 circle, ef_range, plane, etf, Some(ref_dir),
 )
 .into_iter().map(|h| (h.point, h.curve_param)).collect()
 }
 (Curve3::Circle(circle), Surface3::Cylinder(cyl)) => {
 inttools::curve_surface::intersect_circle_cylinder_with_tol(
 circle,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Circle(circle), Surface3::Sphere(sph)) => {
 inttools::curve_surface::intersect_circle_sphere_with_tol(
 circle,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Circle(circle), Surface3::Cone(cone)) => {
 inttools::curve_surface::intersect_circle_cone_with_tol(
 circle,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Plane(plane)) => {
 //  IntAna_IntConicQuad Ellipse  ?Plane
 inttools::ellipse_intersection::intersect_ellipse_plane_with_tol(
 ellipse,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Cylinder(cyl)) => {
 //  ?Partially aligned: numeric fallback, same as OCCT for rare cases
 inttools::ellipse_intersection::intersect_ellipse_cylinder_with_tol(
 ellipse,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Sphere(sph)) => {
 inttools::ellipse_intersection::intersect_ellipse_sphere_with_tol(
 ellipse,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Cone(cone)) => {
 inttools::ellipse_intersection::intersect_ellipse_cone_with_tol(
 ellipse,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Plane(plane)) => {
 //  IntAna_IntConicQuad Parabola  ?Plane
 inttools::parabola_intersection::intersect_parabola_plane_with_tol(
 parabola,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Cylinder(cyl)) => {
 //  ?Partially aligned: numeric fallback
 inttools::parabola_intersection::intersect_parabola_cylinder_with_tol(
 parabola,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Sphere(sph)) => {
 inttools::parabola_intersection::intersect_parabola_sphere_with_tol(
 parabola,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Cone(cone)) => {
 inttools::parabola_intersection::intersect_parabola_cone_with_tol(
 parabola,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Plane(plane)) => {
 //  IntAna_IntConicQuad Hyperbola  ?Plane
 inttools::hyperbola_intersection::intersect_hyperbola_plane_with_tol(
 hyperbola,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Cylinder(cyl)) => {
 //  ?Partially aligned: numeric fallback
 inttools::hyperbola_intersection::intersect_hyperbola_cylinder_with_tol(
 hyperbola,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Sphere(sph)) => {
 inttools::hyperbola_intersection::intersect_hyperbola_sphere_with_tol(
 hyperbola,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Cone(cone)) => {
 inttools::hyperbola_intersection::intersect_hyperbola_cone_with_tol(
 hyperbola,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 _ => {
 // Numeric fallback: sample the curve, find sign changes of the
 // surface implicit function. Works for any Curve3  ?Surface3 pair.
 intersect_edge_face_numeric(&edge_curve, &face_surface, ef_range, etf)
 }
 };

 for (point, edge_param) in hits {
 //  IsPointInFace check for ALL surface types (PaveFiller_5.cxx L523)
 let in_face = self.is_point_in_face(point, face_idx, etf);
 if !in_face {
 let near_face_vert = match &face_surface {
 Surface3::Plane(_) => {
 self.ds.face_boundary_points(face_idx).iter().any(|&vp| {
 (vp - point).length() <= etf
 })
 }
 _ => false,
 };
 if !near_face_vert { continue; }
 }

 //  Always create EF interference for intersection hits.
 // OCCT IntTools_EdgeFace creates a new vertex for each hit, even when
 // the hit coincides with an existing edge endpoint.  SD vertex merging
 // handles near-coincident vertices later (MakeSDVerticesFF in PostTreat).
 // rcad: do NOT skip endpoint-coincident hits  ?they are needed for
 // PutPaveOnCurve to split intersection curve pave blocks.
 let new_v = self.ds.add_vertex(point);
 // Register vertices_on for the new vertex if it's near the edge boundary
 let sv = self.ds.edge_start_vertex_ds(edge_idx);
 let ev = self.ds.edge_end_vertex_ds(edge_idx);
 let tol = etf
 .max(self.ds.vertex_tolerance(sv))
 .max(self.ds.vertex_tolerance(ev));
 if (point - self.ds.vertex_point(sv)).length() <= tol
 || (point - self.ds.vertex_point(ev)).length() <= tol
 {
 self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
 }
 //  Create EF interference for EVERY hit, even at edge endpoints.
 // OCCT IntTools_EdgeFace creates a new vertex for each hit (no dedup).
 // rcad: remove the vertices_on skip check  ?always push interference.
 self.ds.interf_ef.push(InterferenceEF{
 edge: edge_idx,
 face: face_idx,
 point,
 edge_param,
 new_vertex: new_v,
 });
 // Only mark vertices_on if actually inserted (avoid duplicate insert msg)
 if !self.ds.faces[face_idx].face_info.vertices_on.contains(&new_v) {
 self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
 }
 self.ds.edge_paves[edge_idx].push(Pave {
 vertex_idx: new_v,
 param: edge_param,
 });
 }
 }
 /// OCCT HasInterf: skip already-processed pairs
 pub fn skip_redundant_interferences(&self) -> std::collections::HashSet<(usize, usize, u8)> {
 let mut skip_set = std::collections::HashSet::new();

 if !self.use_glue() {
 return skip_set;
 }

 // Skip V-V for shared vertices
 for &(va, vb) in &self.ds.shared_topology.shared_vertices {
 skip_set.insert((va, vb, 0)); // 0 = V-V
 }

 // Skip E-E for shared edges
 for &(ea, eb) in &self.ds.shared_topology.shared_edges {
 skip_set.insert((ea, eb, 2)); // 2 = E-E
 }

 // Skip F-F for fully glued faces
 for &(fa, fb) in &self.ds.shared_topology.fully_glued_faces {
 skip_set.insert((fa, fb, 5)); // 5 = F-F
 }

 skip_set
 }

 /// OCCT myDS->HasInterf: check existing EE interference
 pub(crate) fn has_ee_interf(&self, e1: usize, e2: usize) -> bool {
 self.ds.interf_ee.iter().any(|inf| {
 (inf.e1 == e1 && inf.e2 == e2) || (inf.e1 == e2 && inf.e2 == e1)
 })
 }
}

