use std::collections::{HashSet, BTreeMap};

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Surface3, any_perpendicular};
use rcad_kernel::CurveEval;

use crate::bopalgo::{GlueEnum, fill_map, make_blocks};
use crate::bopds::ds::{DS, DSVertex, ShapeOrigin, InterferenceVV, InterferenceVE, InterferenceVF, InterferenceEE, InterferenceEF};
use crate::bopds::pave::Pave;
use crate::inttools;
use crate::inttools::edge_edge::compute_curve_aabb;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::pave_filler::helpers::*;
use crate::tolerance::*;

/// OCCT IntTools_CommonPrt::Type() ??VERTEX or EDGE.
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
  bIsPBSplittable: bool,
  pb_local_idx: usize,
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
 // OCCT PaveFiller_2.cxx L141-208: PerformVE
 // Groups cross-operand (nV, nE) pairs by edge, then calls IntersectVE.
 // rcad: receives pre-computed cross-operand pairs from BOPDS_Iterator.
 pub(crate) fn perform_ve_bvh(&mut self, pairs: &[(usize, usize)]) {
  // OCCT L143: FillShrunkData(VERTEX, EDGE)
  self.fill_shrunk_data();

  // OCCT L148-152: iSize = myIterator->ExpectedLength()
  let i_size = pairs.len();
  if i_size == 0 {
   return;
  }

  // OCCT L155: NCollection_IndexedDataMap<handle<PaveBlock>, NCollection_List<int>> aMVEPairs
  // rcad: HashMap<edge_idx, Vec<vertex_idx> (edge identified by its first PB)
  let a_vc = self.ds.a_vertex_count;
  let a_ec = self.ds.a_edge_count;
  let mut a_mve_pairs: std::collections::HashMap<usize, Vec<usize>> =
   std::collections::HashMap::new();

  // OCCT L156-205: for (; myIterator->More(); myIterator->Next())
  for &(n_v, n_e) in pairs {
   // rcad: cross-operand filter (OCCT: BOPDS_Iterator enforces this)
   if (n_v < a_vc) == (n_e < a_ec) {
    continue;
   }

   // OCCT L165-168: aSIE.HasSubShape(nV)
   if self.ds.edge_has_vertex(n_e, n_v) {
    continue;
   }

   // OCCT L171-174: aSIE.HasFlag()
   if self.ds.edge_has_flag(n_e) {
    continue;
   }

   // OCCT L176-179: myDS->HasInterf(nV, nE)
   if self.ds.has_interf_ve(n_v, n_e) {
    continue;
   }

   // OCCT L181-184: myDS->HasInterfShapeSubShapes(nV, nE)
   if self.ds.has_interf_ve_via_faces(n_v, n_e) {
    continue;
   }

   // OCCT L186-190: aLPB empty
   if self.ds.edge_pave_blocks(n_e).is_empty() {
    continue;
   }

   // OCCT L192-197: first PB not splittable
   if !self.ds.edge_pave_blocks(n_e)[0].0.read().unwrap().is_splittable {
    continue;
   }

   // OCCT L199-204: group vertices by edge (keyed by first PB)
   a_mve_pairs.entry(n_e).or_default().push(n_v);
  }

  // OCCT L207: IntersectVE(aMVEPairs, ...)
  self.intersect_ve(&a_mve_pairs, true);
 }

 // OCCT PaveFiller_2.cxx L212-395: IntersectVE
 fn intersect_ve(
  &mut self,
  the_ve_pairs: &std::collections::HashMap<usize, Vec<usize>>,
  the_add_interfs: bool,
 ) {
  // OCCT L217-221: aNbVE = theVEPairs.Extent()
  let a_nb_ve = the_ve_pairs.len();
  if a_nb_ve == 0 {
   return;
  }

  // OCCT L223-227: aVEs.SetIncrement(aNbVE)
  if the_add_interfs {
   self.ds.interf_ve.reserve(a_nb_ve);
  }

  // OCCT L230: BOPAlgo_VectorOfVertexEdge aVVE
  // rcad: Vec storing (nV, nE) task data
  struct VeTask {
   n_v: usize,
   n_e: usize,
  }
  let mut a_vve: Vec<VeTask> = Vec::new();

  // OCCT L235: aDMVSD ??map (nVSD, nE) -> list of original vertices
  let mut a_dmv_sd: std::collections::HashMap<(usize, usize), Vec<usize>> =
   std::collections::HashMap::new();

  // OCCT L238-291: for (i = 1; i <= aNbVE; ++i)
  for (&n_e, verts) in the_ve_pairs {
   // OCCT L244: nE = aPB->OriginalEdge() ??rcad: n_e is the edge index directly

   // OCCT L247-254: build aMVPB from all PBs of this edge
   let mut a_mv_pb: std::collections::HashSet<usize> = std::collections::HashSet::new();
   for spb in self.ds.edge_pave_blocks(n_e) {
    let pb = spb.0.read().unwrap();
    a_mv_pb.insert(pb.pave1.vertex_idx);
    a_mv_pb.insert(pb.pave2.vertex_idx);
   }

   // OCCT L256-291: iterate vertex list for this PB
   for &n_v in verts {
    // OCCT L262-263: resolve SD vertex
    let n_vsd = self.ds.has_shape_sd(n_v).unwrap_or(n_v);

    // OCCT L265-268: skip if nVSD is a PB endpoint
    if a_mv_pb.contains(&n_vsd) {
     continue;
    }

    // OCCT L270-277: check if (nVSD, nE) already in aDMVSD
    let a_pair = (n_vsd, n_e);
    if let Some(p_li) = a_dmv_sd.get_mut(&a_pair) {
     // Already added ??just append the original vertex
     p_li.push(n_v);
     continue;
    }

    // OCCT L279-291: new pair ??create solver task
    a_dmv_sd.insert(a_pair, vec![n_v]);
    a_vve.push(VeTask { n_v: n_vsd, n_e });
   }
  }

  // OCCT L294: aNbVE = aVVE.Length()
  let a_nb_ve = a_vve.len();

  // OCCT L302-304: BOPTools_Parallel::Perform(myRunParallel, aVVE, myContext)
  // rcad: sequential execution

  // OCCT L312: NCollection_Map<int> aMEdges
  let mut a_m_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

  // OCCT L315-387: for (i = 0; i < aNbVE; ++i)
  for task in &a_vve {
   // OCCT L321-329: if flag != 0 ??skip / warn
   let res = match self.context.compute_ve(
    self.ds, task.n_v, task.n_e, self.fuzzy_tolerance,
   ) {
    Ok(res) => res,
    Err(_) => {
     // OCCT L324-328: HasErrors ??AddIntersectionFailedWarning
     self.add_intersection_failed_warning(task.n_v, task.n_e);
     continue;
    }
   };

   // OCCT L332-338: extract result
   let a_t = res.param;
   let a_tol_v_new = res.tolerance;
   // OCCT L338: nVx = UpdateVertex(nV, aTolVNew)
   let n_vx = self.update_vertex(task.n_v, a_tol_v_new);

   // OCCT L341-354: Find PB on edge containing aT
   let a_lpb = self.ds.edge_pave_blocks(task.n_e);
   let pb_idx = a_lpb.iter().position(|spb| {
    let pb = spb.0.read().unwrap();
    let (a_t1, a_t2) = pb.range();
    a_t > a_t1 && a_t < a_t2
   });
   let pb_idx = match pb_idx {
    Some(i) => i,
    None => continue,
   };

   // OCCT L360-363: AppendExtPave
   let a_pave = Pave { vertex_idx: n_vx, param: a_t };
   a_lpb[pb_idx].0.write().unwrap().append_ext_pave(a_pave);
   a_m_edges.insert(task.n_e);

   // OCCT L366-387: create interferences
   if the_add_interfs {
    // OCCT L369: BOPDS_Pair aPair(nV, nE)
    let a_pair = (task.n_v, task.n_e);
    // OCCT L370: aDMVSD.Find(aPair)
    if let Some(a_li) = a_dmv_sd.get(&a_pair) {
     // OCCT L371-386: for each original vertex
     for &n_v_old in a_li {
      // OCCT L376-378: create VE interference
      let b_new = self.ds.is_new_vertex(n_vx);
      self.ds.interf_ve.push(InterferenceVE {
       vertex: n_v_old,
       edge: task.n_e,
       param: a_t,
      });
      // OCCT L380: myDS->AddInterf(nVOld, nE)
      self.ds.try_add_interf(n_v_old, task.n_e);
      // OCCT L382-385: if new shape, SetIndexNew
      // rcad: no index_new field on InterferenceVE ??skip SetIndexNew
      _ = b_new; // marker for potential future alignment
     }
    }
   }
  }

  // OCCT L394: SplitPaveBlocks(aMEdges, theAddInterfs)
  if !a_m_edges.is_empty() {
   self.split_pave_blocks(&a_m_edges, the_add_interfs);
  }
 }
 // OCCT BOPAlgo_PaveFiller_3.cxx L145-590: PerformEE
 //
 // OCCT structure:
 //   L147: FillShrunkData(EDGE, EDGE)
 //   L149-150: Iterator init, iSize check
 //   L157-175: variable declarations (aEEs, aMEdges, allocators)
 //   L181-267: Phase 1 -- collect BOPAlgo_EdgeEdge tasks
 //   L269-278: Phase 2 -- parallel execution (BOPTools_Parallel)
 //   L285-556: Phase 3 -- process CommonPrt (VERTEX/EDGE types)
 //   L558-585: Phase 4 -- PerformCommonBlocks + PerformNewVertices + SplitPaveBlocks
 //
 // rcad architecture differences:
 //   - intersect_ee combines computation + InterfEE creation (no CommonPrt)
 //   - treat_new_vertices() called separately from perform() (not inside)
 //   - No aMVCPB / aMPBLPB coupling (common blocks handled elsewhere)
 //   - Sequential execution (no BOPTools_Parallel)
 pub(crate) fn perform_ee_bvh(&mut self, pairs: &[(usize, usize)]) {
  // OCCT L147: FillShrunkData(EDGE, EDGE)
  self.fill_shrunk_data();

  // OCCT L149-150: myIterator->Initialize(EDGE, EDGE)
  // iSize = myIterator->ExpectedLength()
  let i_size = pairs.len();

  // OCCT L152-155: if (!iSize) return
  if i_size == 0 {
   return;
  }

  // rcad EeTask replaces BOPAlgo_EdgeEdge (no TopoDS_Shape, no handle types)
  struct EeTask {
   nE1: usize,
   nE2: usize,
   aT11: f64,
   aT12: f64,
   aTS11: f64,
   aTS12: f64,
   aT21: f64,
   aT22: f64,
   aTS21: f64,
   aTS22: f64,
   nV11: usize,
   nV12: usize,
   nV21: usize,
   nV22: usize,
   b_express_compute: bool,
   b_is_pb_splittable1: bool,
   b_is_pb_splittable2: bool,
  }

  // OCCT L167: NCollection_Map<int> aMEdges
  let mut a_m_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

  // OCCT L178-179: aEEs.SetIncrement(iSize)
  self.ds.interf_ee.reserve(i_size);

  // OCCT L181: for (; myIterator->More(); myIterator->Next())
  let a_ec = self.ds.a_edge_count;
  let mut a_vee: Vec<EeTask> = Vec::new();

  for &(nE1, nE2) in pairs {
   // rcad: cross-operand filter (OCCT: done by BOPDS_Iterator)
   if (nE1 < a_ec) == (nE2 < a_ec) {
    continue;
   }

   // OCCT L189-192: myDS->ShapeInfo(nE1).HasFlag()
   if self.ds.edge_has_flag(nE1) || self.ds.edge_has_flag(nE2) {
    continue;
   }

   // OCCT L200-204: myDS->ChangePaveBlocks(nE1).IsEmpty()
   let a_lpb1 = self.ds.edge_pave_blocks(nE1);
   if a_lpb1.is_empty() {
    continue;
   }

   // OCCT L206-210: myDS->ChangePaveBlocks(nE2).IsEmpty()
   let a_lpb2 = self.ds.edge_pave_blocks(nE2);
   if a_lpb2.is_empty() {
    continue;
   }

   // rcad: additional skip conditions (OCCT applies these earlier)
   if self.ds.has_interf_ee(nE1, nE2) {
    continue;
   }
   if self.ds.is_edge_degenerated(nE1) || self.ds.is_edge_degenerated(nE2) {
    continue;
   }

   // OCCT L215-266: PB pair iteration
   for pb1 in a_lpb1.iter() {
    let pb1_r = pb1.0.read().unwrap();

    // OCCT L222-229: GetPBBox
    let (aT11, aT12) = pb1_r.range();
    let (aTS11, aTS12, b_is_pb_splittable1) = if pb1_r.has_shrunk_data() {
     let (ts1, ts2, spl) = pb1_r.shrunk_data();
     (ts1, ts2, spl)
    } else {
     (aT11, aT12, false)
    };

    // OCCT L231: aPB1->Indices(nV11, nV12)
    let (nV11, nV12) = pb1_r.indices();
    drop(pb1_r);

    // OCCT L233-265: aIt2.Initialize(aLPB2); for (; aIt2.More(); aIt2.Next())
    for pb2 in a_lpb2.iter() {
     let pb2_r = pb2.0.read().unwrap();

     // OCCT L238-243: GetPBBox
     let (aT21, aT22) = pb2_r.range();
     let (aTS21, aTS22, b_is_pb_splittable2) = if pb2_r.has_shrunk_data() {
      let (ts1, ts2, spl) = pb2_r.shrunk_data();
      (ts1, ts2, spl)
     } else {
      (aT21, aT22, false)
     };

     // OCCT L245-248: if (aBB1.IsOut(aBB2)) continue
     // rcad: compute AABB of each PB's curve segment (equivalent to GetPBBox + bbox cache).
     let e1_curve = &self.ds.edges[nE1].curve;
     let e2_curve = &self.ds.edges[nE2].curve;
     let bbox_tol = self.ds.edge_tolerance(nE1).max(self.ds.edge_tolerance(nE2))
         + crate::tolerance::TOLERANCE_ABS;
     let bbox1 = compute_curve_aabb(e1_curve, aTS11.min(aTS12), aTS11.max(aTS12), bbox_tol);
     let bbox2 = compute_curve_aabb(e2_curve, aTS21.min(aTS22), aTS21.max(aTS22), bbox_tol);
     if bbox1.0.x > bbox2.1.x || bbox1.1.x < bbox2.0.x
      || bbox1.0.y > bbox2.1.y || bbox1.1.y < bbox2.0.y
      || bbox1.0.z > bbox2.1.z || bbox1.1.z < bbox2.0.z
     {
      drop(pb2_r);
      continue;
     }

     // OCCT L250: aPB2->Indices(nV21, nV22)
     let (nV21, nV22) = pb2_r.indices();

     // OCCT L252: bExpressCompute = same vertex bounds
     let b_express_compute =
      (nV11 == nV21 && nV12 == nV22) || (nV12 == nV21 && nV11 == nV22);

     drop(pb2_r);

     a_vee.push(EeTask {
      nE1,
      nE2,
      aT11,
      aT12,
      aTS11,
      aTS12,
      aT21,
      aT22,
      aTS21,
      aTS22,
      nV11,
      nV12,
      nV21,
      nV22,
      b_express_compute,
      b_is_pb_splittable1,
      b_is_pb_splittable2,
     });
    }
   }
  }

  // OCCT L269: aNbEdgeEdge = aVEdgeEdge.Length()
  let a_nb_edge_edge = a_vee.len();

  // OCCT L285-556: Process results
  for k in 0..a_nb_edge_edge {
   let task = &a_vee[k];
   let nE1 = task.nE1;
   let nE2 = task.nE2;

   // Compute shrunk ranges for this task
   let sr1 = [task.aTS11.min(task.aTS12), task.aTS11.max(task.aTS12)];
   let sr2 = [task.aTS21.min(task.aTS22), task.aTS21.max(task.aTS22)];

   let mut modified: std::collections::HashSet<usize> = std::collections::HashSet::new();
   self.intersect_ee(nE1, nE2, sr1, sr2, &mut modified);

   if !modified.is_empty() {
    for &e in &modified {
     a_m_edges.insert(e);
    }
   }
  }

  // OCCT L558-560: PerformCommonBlocks + UpdateVerticesOfCB
  if !a_vee.is_empty() {
   crate::bopds::tools::perform_common_blocks(&mut self.ds);
  }
  self.update_vertices_of_cb();
  // OCCT L565: PerformNewVertices (TreatNewVertices + IntersectVE)
  self.treat_new_vertices();

  // OCCT L571-585: if (aMEdges.Extent()) { SplitPaveBlocks(aMEdges, false); }
  if !a_m_edges.is_empty() {
   self.split_pave_blocks(&a_m_edges, false);
  }
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
/// OCCT PaveFiller_4.cxx L139-301: PerformVF
pub(crate) fn perform_vf_bvh(&mut self, pairs: &[(usize, usize)]) {
 // OCCT L141: myIterator->Initialize(TopAbs_VERTEX, TopAbs_FACE)
 let a_vc = self.ds.a_vertex_count;
 let a_fc = self.ds.a_face_count;
 // OCCT L142: iSize = myIterator->ExpectedLength()
 let i_size = pairs.len();
 // OCCT L147-160: myGlue == GlueFull handled in mod.rs
 //
 // OCCT L162: InterfVF()  -- aVFs
 // OCCT L163-170: if (!iSize)
 if i_size == 0 {
   // OCCT L165: iSize = 10
   // OCCT L166: aVFs.SetIncrement(iSize)
   self.ds.interf_vf.reserve(10);
   // OCCT L168: TreatVerticesEE()
   self.treat_vertices_ee();
   return;
 }
 // OCCT L172-174: variable declarations
 // OCCT L174: BOPAlgo_VectorOfVertexFace aVVF
 let mut a_vv_f: Vec<VfTask> = Vec::new();
 //
 // OCCT L176: aVFs.SetIncrement(iSize)
 self.ds.interf_vf.reserve(i_size);
 //
 // OCCT L178-180: NCollection_DataMap<BOPDS_Pair, NCollection_Map<int>> aMVFPairs
 let mut a_mvf_pairs: std::collections::HashMap<(usize, usize), Vec<usize>> =
   std::collections::HashMap::new();
 //
 // OCCT L181: for (; myIterator->More(); myIterator->Next())
 for &(nV, nF) in pairs {
   // OCCT L183-186: UserBreak check (not ported)
   //
   // OCCT L187: myIterator->Value(nV, nF)
   //
   // OCCT L189-192: IsSubShape
   let same_range = (nV < a_vc) == (nF < a_fc);
   if same_range { continue; }
   //
   // OCCT L194-197: if (myDS->HasInterf(nV, nF)) continue;
   if self.ds.has_interf_vf(nV, nF) { continue; }
   //
   // OCCT L199: myDS->ChangeFaceInfo(nF)
   self.ds.face_info_mut(nF);
   //
   // OCCT L200-203: if (myDS->HasInterfShapeSubShapes(nV, nF)) continue;
   //   Checks if nV has interference with any sub-shape (edge) of nF.
   {
     let mut has_interf = false;
     let face_edges = self.ds.face_boundary_edges(nF).to_vec();
     for &ei in &face_edges {
       if self.ds.has_interf_ve(nV, ei) {
         has_interf = true;
         break;
       }
     }
     if !has_interf {
       let inner = self.ds.face_inner_boundary(nF);
       for iw in inner {
         for &(ei, _) in iw {
           if self.ds.has_interf_ve(nV, ei) {
             has_interf = true;
             break;
           }
         }
         if has_interf { break; }
       }
     }
     if has_interf { continue; }
   }
   //
   // OCCT L205-209: SD resolution
   let nVx = self.ds.has_shape_sd(nV).unwrap_or(nV);
   //
   // OCCT L211-220: aMVFPairs dedup (key = nVx, nF)
   let key = (nVx, nF);
   let entry = a_mvf_pairs.entry(key).or_default();
   entry.push(nV);
   if entry.len() > 1 {
     // OCCT L216: continue - already have a task for this SD-face pair
     continue;
   }
   //
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
 // OCCT L234-240: SetProgressRange (not ported)
 // OCCT L242: BOPTools_Parallel::Perform(myRunParallel, aVVF, myContext);
 for task in &mut a_vv_f {
   let tf = self.vf_tol(task.nV, task.nF);
   *task = compute_vf_on_face(&self.ds, tf, self.fuzzy_tolerance, task.nV, task.nF);
 }
 // OCCT L244-247: UserBreak check (not ported)
 //
 // ------------------------------------------------------------------
 // Phase 3: Process results (OCCT L249-298)
 // ------------------------------------------------------------------
 // OCCT L249: for (k = 0; k < aNbVF; ++k)
 for task in &a_vv_f {
   // OCCT L251-254: UserBreak check (not ported)
   //
   // OCCT L257-265: iFlag = aVertexFace.Flag();
   if !task.is_on { continue; }
   //
   // OCCT L268: aVertexFace.Indices(nVx, nF)
   let mut nVx = task.nV;
   let nF = task.nF;
   // OCCT L269: aVertexFace.Parameters(aT1, aT2)
   let a_t1 = task.proj_u;
   let a_t2 = task.proj_v;
   // OCCT L270: double aTolVNew = aVertexFace.VertexNewTolerance()
   let a_tol_vnew = task.proj_dist;
   //
   // OCCT L272-273: aMVFPairs.Find(aVFPair) - get all original vertices
   let key = (nVx, nF);
   let orig_verts: Vec<usize> = match a_mvf_pairs.get(&key) {
     Some(v) => v.clone(),
     None => vec![nVx],
   };
   // OCCT L275: for (; itMV.More(); itMV.Next())
   for &nV in &orig_verts {
     // OCCT L279-281: BOPDS_InterfVF& aVF = aVFs.Appended();
     //   aVF.SetIndices(nV, nF);
     //   aVF.SetUV(aT1, aT2);
     // OCCT L283: myDS->AddInterf(nV, nF);
     self.ds.try_add_interf(nV, nF);
     // OCCT L286: nVx = UpdateVertex(nV, aTolVNew);  [shadows outer nVx]
     nVx = self.update_vertex(nV, a_tol_vnew);
     // OCCT L289-292: if (myDS->IsNewShape(nVx)) { aVF.SetIndexNew(nVx); }
     let idx_new = if self.ds.is_new_vertex(nVx) { Some(nVx) } else { None };
     // OCCT L279-281 (InterfVF pushed after UpdateVertex due to Rust borrow checker)
     self.ds.interf_vf.push(InterferenceVF {
       vertex: nV, face: nF,
       u: a_t1, v: a_t2,
       index_new: idx_new,
     });
   }
   // OCCT L295-297: FaceInfo VerticesIn (nVx = last shadowed value from UpdateVertex)
   let a_fi = self.ds.face_info_mut(nF);
   a_fi.vertices_in.insert(nVx);
 } // for (k=0; k < aNbVF; ++k) {
 //
 // OCCT L300: TreatVerticesEE()
 self.treat_vertices_ee();
}
 /// OCCT PaveFiller_5.cxx L165-592: PerformEF
 pub(crate) fn perform_ef(&mut self, pairs: &[(usize, usize)]) {
 self.fill_shrunk_data(); // OCCT L167
 let i_size = pairs.len();
 if i_size == 0 {
 return;
 }
 // OCCT L179-192: GlueFull early return
 if self.glue == GlueEnum::GlueFull {
 for &(nE, nF) in pairs {
 if !self.ds.edge_has_flag(nE) {
 self.ds.face_info_mut(nF);
 }
 }
 return;
 }
 // OCCT L194-217: variable declarations
 let mut a_mi_efc: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_v_edge_face: Vec<EfTask> = Vec::new();
 self.ds.interf_ef.reserve(i_size);
 let a_vc = self.ds.a_vertex_count;
 let a_fc = self.ds.a_face_count;
 //
 // ==================================================================
 // Phase 1: Collect EF tasks (OCCT L219-307)
 // ==================================================================
 for &(nE, nF) in pairs {
 // OCCT L227-231: HasFlag / degenerated
 if self.ds.edge_has_flag(nE) || self.ds.is_edge_degenerated(nE) {
 continue;
 }
 // OCCT L235: aBBF = myDS->ShapeInfo(nF).Box() -- face bounding box from DS shape info.
 let tol = self.ds.edge_tolerance(nE).max(CONFUSION);
 let (a_bbf_min, a_bbf_max): (DVec3, DVec3) = self.ds.shape_info.get(nF)
 .and_then(|si| {
 let mn = si.box_min?;
 let mx = si.box_max?;
 Some((mn - DVec3::splat(tol), mx + DVec3::splat(tol)))
 })
 .unwrap_or_else(|| {
 // Fallback: compute from boundary vertices
 let f = &self.ds.faces[nF];
 let mut mn = DVec3::splat(f64::INFINITY);
 let mut mx = DVec3::splat(f64::NEG_INFINITY);
 for &vi in &f.boundary_verts {
 if vi < self.ds.vertices.len() {
 let p = self.ds.vertex_point(vi);
 mn = mn.min(p); mx = mx.max(p);
 }
 }
 if let Surface3::Sphere(s) = &f.surface {
 let r = s.radius.abs();
 mn = mn.min(s.center - DVec3::splat(r));
 mx = mx.max(s.center + DVec3::splat(r));
 }
 (mn - DVec3::splat(tol), mx + DVec3::splat(tol))
 });
 // OCCT L237-241: FaceInfo -- ChangeFaceInfo + On/In sets
 let a_mpbf: Vec<usize> = {
 let fi = self.ds.face_info(nF);
 fi.pave_blocks_on.iter().copied().collect()
 };
 let a_mv_in: Vec<usize> = {
 let fi = self.ds.face_info(nF);
 fi.vertices_in.iter().copied().collect()
 };
 let a_mv_on: Vec<usize> = {
 let fi = self.ds.face_info(nF);
 fi.vertices_on.iter().copied().collect()
 };
 // OCCT L243-244: aTolE, aTolF
 let a_tol_e = self.ds.edge_tolerance(nE);
 let a_tol_f = self.ds.face_tolerance(nF);
 // OCCT L246-248: ChangePaveBlocks + PB iterator
 let n_pbs = self.ds.edge_pave_blocks(nE).len();
 for pb_local_idx in 0..n_pbs {
 let pb = &self.ds.edge_pave_blocks(nE)[pb_local_idx];
 let pb_ref = pb.0.read().unwrap();
 // OCCT L256-259: aPBR = RealPaveBlock(aPB)
 let pb_key = pb_ref.new_edge.unwrap_or(pb_ref.original_edge);
 if a_mpbf.contains(&pb_key) {
 continue;
 }
 // OCCT L262-266: GetPBBox
 let (aT1, aT2) = pb_ref.range();
 let (aTS1, aTS2, bb_min, bb_max, has_box) = match pb_ref.shrunk_range {
 Some(sr) => {
 let (b_mn, b_mx) = match pb_ref.my_shrunk_box {
 Some((mn, mx)) => (mn, mx),
 None => {
 let tol = a_tol_e.max(CONFUSION);
 let p1 = self.ds.edges[nE].curve.point_at(sr[0]);
 let p2 = self.ds.edges[nE].curve.point_at(sr[1]);
 (p1.min(p2) - DVec3::splat(tol), p1.max(p2) + DVec3::splat(tol))
 }
 };
 (sr[0], sr[1], b_mn, b_mx, true)
 }
 None => (0.0, 0.0, DVec3::ZERO, DVec3::ZERO, false)
 };
 if !has_box {
 continue;
 }
 // OCCT L268-271: AABB overlap check
 if bb_max.x < a_bbf_min.x || bb_min.x > a_bbf_max.x
 || bb_max.y < a_bbf_min.y || bb_min.y > a_bbf_max.y
 || bb_max.z < a_bbf_min.z || bb_min.z > a_bbf_max.z
 {
 continue;
 }
 // OCCT L273-276: bExpressCompute
 let (nV1, nV2) = pb_ref.indices();
 let bV1 = a_mv_in.contains(&nV1) || a_mv_on.contains(&nV1);
 let bV2 = a_mv_in.contains(&nV2) || a_mv_on.contains(&nV2);
 let b_express_compute = bV1 && bV2;
 let b_is_pb_splittable = pb_ref.is_splittable;
 let a_pb_range = [aT1.min(aT2), aT1.max(aT2)];
 drop(pb_ref);
 // OCCT L289-292: CorrectRange for shrunk range
 let tol_ef = self.ef_tol(nE, nF);
 let _a_sr_corrected = Self::correct_range_for_face(
 &self.ds.edges[nE].curve, tol_ef, [aTS1.min(aTS2), aTS1.max(aTS2)]);
 // OCCT L294-297: CorrectRange for PB range
 let a_pb_corrected = Self::correct_range_for_face(
 &self.ds.edges[nE].curve, tol_ef, a_pb_range);
 if a_pb_corrected[1] - a_pb_corrected[0] <= tol_ef {
 continue;
 }
 // OCCT L299-305: Save to myFPBDone
 self.fpbdone.entry(nF).or_default().insert(pb_key);
 // OCCT L278: aEdgeFace = aVEdgeFace.Appended()
 a_v_edge_face.push(EfTask {
 nE, nF, nV1, nV2,
 aT1, aT2,
 aTS1, aTS2,
 bExpressCompute: b_express_compute,
 bIsPBSplittable: b_is_pb_splittable,
 pb_local_idx,
 hits: Vec::new(),
 });
 } // for PB
 } // for pairs
 //
 // ==================================================================
 // Phase 2: EF computation (OCCT L309-317)
 // ==================================================================
 // OCCT L317: BOPTools_Parallel::Perform(myRunParallel, aVEdgeFace, myContext);
 let a_nb_edge_face = a_v_edge_face.len();
 for index in 0..a_nb_edge_face {
 let task = &a_v_edge_face[index];
 let pb_range = [task.aT1.min(task.aT2), task.aT1.max(task.aT2)];
 let etr = self.ds.edges[task.nE].t_range;
 let ef_range = [pb_range[0].max(etr[0]), pb_range[1].min(etr[1])];
 let tol_ef = self.ef_tol(task.nE, task.nF);
 if ef_range[1] - ef_range[0] <= tol_ef { continue; }
 let ef_corr = Self::correct_range_for_face(
 &self.ds.edges[task.nE].curve, tol_ef, ef_range);
 if ef_corr[1] - ef_corr[0] <= tol_ef { continue; }
 let hits = compute_ef_hits(&self.ds, task.nE, task.nF, &ef_corr);
 drop(task);
 a_v_edge_face[index].hits = hits;
 }
 //
 // ==================================================================
 // Phase 3: Process results (OCCT L324-571)
 // ==================================================================
 for k in 0..a_nb_edge_face {
 // OCCT L330-336: aEdgeFace.IsDone() / HasErrors() check
 let nE = a_v_edge_face[k].nE;
 let nF = a_v_edge_face[k].nF;
 // OCCT L340-344: aE, aF, aTolE, aTolF
 let a_tol_e = self.ds.edge_tolerance(nE);
 let a_tol_f = self.ds.face_tolerance(nF);
 // OCCT L346-362: aCPrts, aNbCPrts
 let a_nb_cprts = a_v_edge_face[k].hits.len();
 if a_nb_cprts == 0 {
 // OCCT L350-361: MinimalDistance handling
 {
 let a_t1_ef = a_v_edge_face[k].aT1;
 let a_t2_ef = a_v_edge_face[k].aT2;
 let span = (a_t2_ef - a_t1_ef).abs();
 if span > crate::tolerance::TOLERANCE_ABS {
 let mut min_dist = f64::MAX;
 let edge_curve = &self.ds.edges[nE].curve;
 for s in 0..5 {
 let t = a_t1_ef + span * (s as f64 / 4.0);
 let pt = edge_curve.point_at(t);
 if let Some((_, _, dist)) = self.context.proj_ps(self.ds, nF, pt) {
 if dist < min_dist { min_dist = dist; }
 }
 }
 if min_dist < f64::MAX && min_dist > a_tol_e + a_tol_f {
 let entry = self.distances.entry((nE, nF)).or_default();
 entry.push(crate::pave_filler::EdgeRangeDistance {
 first: a_t1_ef,
 last: a_t2_ef,
 distance: min_dist,
 });
 }
 }
 }
 continue;
 }
 let tol_ef = self.ef_tol(nE, nF);
 // OCCT L364-371: anewSR, aPB, Range, Indices, IsSplittable
 let pb_local_idx = a_v_edge_face[k].pb_local_idx;
 let nV = [a_v_edge_face[k].nV1, a_v_edge_face[k].nV2];
 let b_is_pb_splittable = a_v_edge_face[k].bIsPBSplittable;
 let a_t1 = a_v_edge_face[k].aT1;
 let a_t2 = a_v_edge_face[k].aT2;
 let mut a_ts1 = a_v_edge_face[k].aTS1;
 let mut a_ts2 = a_v_edge_face[k].aTS2;
 // OCCT L373-380: ReduceIntersectionRange for VERTEX type
 if a_nb_cprts > 0 {
 let first_hit = &a_v_edge_face[k].hits[0];
 if matches!(first_hit, EfHit::Vertex { .. }) {
 self.reduce_ef_intersection_range(nV[0], nV[1], nE, nF, &mut a_ts1, &mut a_ts2);
 }
 }
 // OCCT L382: IntTools_Range aR1(aT1, aTS1), aR2(aTS2, aT2)
 let a_r1 = [a_t1.min(a_ts1), a_t1.max(a_ts1)];
 let a_r2 = [a_ts2.min(a_t2), a_ts2.max(a_t2)];
 // OCCT L384-386: FaceInfo On/In
 let a_fi = self.ds.face_info_mut(nF);
 let a_mif_on: std::collections::HashSet<usize> = a_fi.vertices_on.iter().copied().collect();
 let a_mif_in: std::collections::HashSet<usize> = a_fi.vertices_in.iter().copied().collect();
 drop(a_fi);
 // OCCT L388-394: bLinePlane
 let b_line_plane = matches!(
 (&self.ds.edges[nE].curve, &self.ds.faces[nF].surface),
 (Curve3::Line(_), Surface3::Plane(_))
 );
 // OCCT L396-570: for each CommonPrt
 for i in 0..a_nb_cprts {
 let a_cpart = &a_v_edge_face[k].hits[i];
 match a_cpart {
 EfHit::Vertex { point, param: a_t } => {
 // ============================================================
 // OCCT L406-543: case TopAbs_VERTEX
 // ============================================================
 // OCCT L415-419: IsInRange(aR1, aR, aTolToDecide)
 let a_tol_to_decide = 5e-8;
 let mut b_is_on_pave = [false, false];
 b_is_on_pave[0] = (a_t - a_r1[0]).abs() <= a_tol_to_decide
 || (a_t - a_r1[1]).abs() <= a_tol_to_decide;
 b_is_on_pave[1] = (a_t - a_r2[0]).abs() <= a_tol_to_decide
 || (a_t - a_r2[1]).abs() <= a_tol_to_decide;
 // OCCT L421-439: if both on pave or (bLinePlane && one on pave)
 if (b_is_on_pave[0] && b_is_on_pave[1])
 || (b_line_plane && (b_is_on_pave[0] || b_is_on_pave[1]))
 {
 // OCCT L423-425: CheckFacePaves(nV[0], aMIFOn, aMIFIn)
 let bv0 = a_mif_on.contains(&nV[0]) || a_mif_in.contains(&nV[0]);
 let bv1 = a_mif_on.contains(&nV[1]) || a_mif_in.contains(&nV[1]);
 if bv0 && bv1 {
 // OCCT L427-437: EDGE-type treatment — edge lies on face
 self.ds.interf_ef.push(InterferenceEF {
 edge: nE, face: nF,
 point: *point,
 edge_param: *a_t,
 new_vertex: 0,
 });
 self.ds.try_add_interf(nE, nF);
 a_mi_efc.insert(nF);
 continue;
 }
 // OCCT L448-455: one vertex NOT on face → mark as processed, no EF
 self.ds.try_add_interf(nE, nF);
 continue;
 }
 // OCCT L442-444: splittable check
 if !b_is_pb_splittable {
 continue;
 }
 // OCCT L447-457: ForceInterfVF for on-pave vertices
 for j in 0..2 {
 if b_is_on_pave[j] {
 let bv = a_mif_on.contains(&nV[j]) || a_mif_in.contains(&nV[j]);
 if !bv {
 b_is_on_pave[j] = self.force_interf_vf_pair(nV[j], nF);
 }
 }
 }
 // OCCT L459-502: if on-pave -> real intersection check + update vertex
 if b_is_on_pave[0] || b_is_on_pave[1] {
 // OCCT L482-501: UpdateVertex for near-endpoint vertices
 for j in 0..2 {
 if b_is_on_pave[j] {
 let dist_pp = (self.ds.vertex_point(nV[j]) - *point).length();
 let a_tol = self.ds.vertex_tolerance(nV[j]);
 let mut a_max_dist = 1e4 * a_tol;
 if a_tol < 0.01 {
 a_max_dist = a_max_dist.min(0.1);
 }
 if dist_pp < a_max_dist {
 self.update_vertex(nV[j], dist_pp);
 self.verts_to_avoid_extension.insert(nV[j]);
 }
 }
 }
 continue;
 }
 // OCCT L505-508: CheckFacePaves(aVnew, aMIFOn)
 {
 let near_face_vx = a_mif_on.iter().chain(a_mif_in.iter()).any(|&vi| {
 vi < self.ds.vertices.len()
 && (self.ds.vertex_point(vi) - *point).length() <= a_tol_e.max(a_tol_f)
 });
 if near_face_vx {
 continue;
 }
 }
 // OCCT L510-519: aTolVnew computation
 let mut a_tol_vnew = a_tol_e.max(a_tol_f);
 if b_line_plane {
 a_tol_vnew = a_tol_vnew.max(tol_ef * 10.0);
 }
 // OCCT L523-526: IsPointInFace check
 if !self.is_point_in_face(*point, nF, a_tol_vnew) {
 continue;
 }
 // OCCT L528-542: Create EF interference
 a_mi_efc.insert(nF);
 let new_v = self.ds.add_vertex(*point);
 self.ds.interf_ef.push(InterferenceEF {
 edge: nE, face: nF,
 point: *point,
 edge_param: *a_t,
 new_vertex: new_v,
 });
 self.ds.try_add_interf(nE, nF);
 // rcad: update face info and edge paves
 self.ds.faces[nF].face_info.vertices_on.insert(new_v);
 if nE < self.ds.edge_paves.len() {
 self.ds.edge_paves[nE].push(Pave {
 vertex_idx: new_v,
 param: *a_t,
 });
 }
 } // EfHit::Vertex
 EfHit::Edge { t1, t2 } => {
 // ============================================================
 // OCCT L545-565: case TopAbs_EDGE
 // ============================================================
 a_mi_efc.insert(nF);
 let mid_t = (t1 + t2) * 0.5;
 let mid_pt = self.ds.edges[nE].curve.point_at(mid_t);
 // OCCT L553-554: CheckFacePaves
 let bv0 = a_mif_on.contains(&nV[0]) || a_mif_in.contains(&nV[0]);
 let bv1 = a_mif_on.contains(&nV[1]) || a_mif_in.contains(&nV[1]);
 // OCCT always appends the interference
 self.ds.interf_ef.push(InterferenceEF {
 edge: nE, face: nF,
 point: mid_pt,
 edge_param: mid_t,
 new_vertex: 0,
 });
 // OCCT L555-558: if (!bV[0] || !bV[1]) { myDS->AddInterf; break; }
 if !bv0 || !bv1 {
 self.ds.try_add_interf(nE, nF);
 } else {
 // OCCT L560-564: SetCommonPart + AddInterf + FillMap
 self.ds.try_add_interf(nE, nF);
 }
 } // EfHit::Edge
 } // match
 } // for each hit
 } // for each task
 //
 // ==================================================================
 // Phase 4: Post-treatment (OCCT L576-592)
 // ==================================================================
 // OCCT L576: BOPAlgo_Tools::PerformCommonBlocks(aMPBLI, ...)
 // OCCT L577: UpdateVerticesOfCB()
 // OCCT L578: PerformNewVertices(aMVCPB, ...)
 // OCCT L585: myDS->UpdateFaceInfoIn(aMIEFC)
 if !a_v_edge_face.is_empty() {
 crate::bopds::tools::perform_common_blocks(&mut self.ds);
 }
 self.update_vertices_of_cb();
 self.treat_new_vertices();
 for &fi in &a_mi_efc {
 self.ds.update_face_info_in(fi);
 }
 }
 /// OCCT PaveFiller_5.cxx L685-768: ReduceIntersectionRange
 fn reduce_ef_intersection_range(
 &self,
 the_v1: usize, the_v2: usize,
 the_e: usize, the_f: usize,
 the_ts1: &mut f64, the_ts2: &mut f64,
 ) {
 if !self.ds.is_new_vertex(the_v1) && !self.ds.is_new_vertex(the_v2) {
 return;
 }
 let has_interf_shape_sub_shapes = {
 let face_edges: std::collections::HashSet<usize> =
 self.ds.face_boundary_edges(the_f).iter().copied().collect();
 self.ds.interf_ee.iter().any(|inf| {
 let involves_our_edge = inf.e1 == the_e || inf.e2 == the_e;
 let involves_face_edge = face_edges.contains(&inf.e1) || face_edges.contains(&inf.e2);
 involves_our_edge && involves_face_edge
 })
 };
 if !has_interf_shape_sub_shapes {
 return;
 }
 let a_nb_ees = self.ds.interf_ee.len();
 if a_nb_ees == 0 {
 return;
 }
 let face_edges: std::collections::HashSet<usize> =
 self.ds.face_boundary_edges(the_f).iter().copied().collect();
 for inf in &self.ds.interf_ee {
 let nv = inf.new_vertex;
 if nv != the_v1 && nv != the_v2 {
 continue;
 }
 let involves_our_edge = inf.e1 == the_e || inf.e2 == the_e;
 let involves_face_edge = face_edges.contains(&inf.e1) || face_edges.contains(&inf.e2);
 if !involves_our_edge || !involves_face_edge {
 continue;
 }
 let ee_param = if inf.e1 == the_e { inf.param1 } else { inf.param2 };
 if nv == the_v1 {
 if *the_ts1 < ee_param {
 *the_ts1 = ee_param;
 }
 } else {
 if *the_ts2 > ee_param {
 *the_ts2 = ee_param;
 }
 }
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
 // OCCT BOPAlgo_PaveFiller_1.cxx L45-132: PerformVV
 pub(crate) fn perform_vv(&mut self, pairs: &[(usize, usize)]) {
   // L47-51: n1, n2, iFlag, aSize; myIterator->Initialize(VERTEX, VERTEX)
   // L50-51: myIterator->Initialize(TopAbs_VERTEX, TopAbs_VERTEX);
   //         aSize = myIterator->ExpectedLength();
   let a_size = pairs.len();
   // L53-56: if (!aSize) return
   if a_size == 0 {
     return;
   }
   // L58-59: InterfVV().SetIncrement(aSize)
   self.ds.interf_vv.reserve(a_size);
   // L62-64: aAllocator, aMILI, aMBlocks
   let mut a_mili: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
   // L66-98: 1. Map V/LV
   for &(n1, n2) in pairs {
     // L71-74: UserBreak check (not ported)
     // L75: myIterator->Value(n1, n2)
     //
     // L77-81: if HasInterf -> FillMap + continue
     // OCCT: myDS->HasInterf checks global fence myInterfTB
     let key = if n1 < n2 { (n1, n2) } else { (n2, n1) };
     if self.ds.interf_tb.contains(&key) {
       fill_map(&mut a_mili, n1, n2);
       continue;
     }
     // L84-88: Resolve SD vertices (HasShapeSD)
     let n1sd: usize = self.ds.has_shape_sd(n1).unwrap_or(n1);
     let n2sd: usize = self.ds.has_shape_sd(n2).unwrap_or(n2);
     // L90-93: ComputeVV(aV1, aV2, myFuzzyValue)
     //   OCCT: BRep_Tool::Tolerance(aV1) + BRep_Tool::Tolerance(aV2) + myFuzzyValue
     let a_tol = self.ds.vertex_tolerance(n1sd) + self.ds.vertex_tolerance(n2sd) + self.tol();
     let a_sq_dist = (self.ds.vertex_point(n1sd) - self.ds.vertex_point(n2sd)).length_squared();
     let i_flag = if a_sq_dist <= a_tol * a_tol { 0 } else { 1 };
     // L94-97: if !iFlag -> FillMap
     if i_flag == 0 {
       fill_map(&mut a_mili, n1, n2);
     }
   }
   // L100-101: 2. Make blocks
   let a_m_blocks: Vec<Vec<usize>> = make_blocks(&a_mili);
   // L103-113: 3. Make SD vertices
   for block in &a_m_blocks {
     // L107-110: UserBreak check (not ported)
     // L111-112: MakeSDVertices(aLI)
     self.make_sd_vertices_vv(block);
   }
   // L115-127: 4. InitPaveBlocksForVertex for each SD vertex source
   // L117: ShapesSD()
   let a_dmii: std::collections::HashSet<usize> =
     self.ds.shape_sd.sd_vertices_iter().map(|&(k, _)| k).collect();
   for &n1_key in &a_dmii {
     // L121-124: UserBreak check (not ported)
     // L125-126: InitPaveBlocksForVertex(n1)
     self.ds.init_pave_blocks_for_vertex(n1_key);
   }
   // L129-131: aMBlocks.Clear(); aMILI.Clear() -- handled by Rust Drop
 }

 /// OCCT BOPAlgo_PaveFiller::MakeSDVertices (PaveFiller_1.cxx L136-233).
 /// Merges a connected group of vertices into a single SD vertex.
 /// If any member already has an SD partner (nSD), that SD vertex is
 /// updated in-place.  Otherwise a new vertex is appended to the DS.
 /// Every pair in the block gets AddShapeSD + a VV interference record
 /// pointing to the merged vertex.
 pub(super) fn make_sd_vertices_vv(&mut self, block: &[usize]) {
 // L136-138: return early if fewer than 2 vertices
 if block.len() < 2 {
   return;
 }
 // L141-161: 1. Collect vertices + track existing SD partner
 let mut n_sd: Option<usize> = None;
 let mut a_lv: Vec<usize> = Vec::with_capacity(block.len());
 for &n_x in block {
   // L145-158: check if vertex already has an SD partner
   if let Some(n_sd1) = self.ds.has_shape_sd(n_x) {
     if n_sd.is_none() {
       // L148-153: keep the first SD vertex as the merge target
       n_sd = Some(n_sd1);
     }
   }
   // L159-160: add vertex to aLV list
   a_lv.push(n_x);
 }
 // L162: MakeVertex(aLV, aVn) — compute centroid + bounding tolerance.
 // OCCT calls BRepLib::BoundingVertex to compute centroid and tolerance
 // large enough to enclose all input vertices.
 let centroid: DVec3 = a_lv.iter()
   .map(|&vi| self.ds.vertex_point(vi))
   .fold(DVec3::ZERO, |acc, p| acc + p) / a_lv.len() as f64;
 let bounding_tol: f64 = a_lv.iter()
   .map(|&vi| (self.ds.vertex_point(vi) - centroid).length() + self.ds.vertex_tolerance(vi))
   .fold(TOLERANCE_ABS, |acc, d| acc.max(d));
 // L163-179: 2. Determine nV — either update existing SD or append new
 let n_v: usize;
 if let Some(n_sd_idx) = n_sd {
   // L166-171: update existing SD vertex in-place (position + tolerance)
   self.ds.vertex_data_mut(n_sd_idx).point = centroid;
   self.ds.vertex_data_mut(n_sd_idx).tolerance = bounding_tol;
   n_v = n_sd_idx;
 } else {
   // L176-179: append new vertex to DS
   n_v = self.ds.vertices.len();
   self.ds.push_vertex(DSVertex {
     point: centroid,
     origin: None,
     geom_tol: bounding_tol,
     is_internal: true,
     location: 0,
   }, None);
   // push_vertex does not maintain shape_info; push a matching entry
   if n_v < self.ds.vertex_shape_idx.len() {
     let si = self.ds.vertex_shape_idx[n_v];
     if si >= self.ds.shape_info.len() {
       let mut new_si = crate::bopds::ds::types::ShapeInfo::new(
         rcad_kernel::topods::ShapeType::Vertex);
       new_si.is_new = true;
       new_si.rank = 0;
       new_si.box_min = Some(centroid);
       new_si.box_max = Some(centroid);
       new_si.box_gap = bounding_tol + TOLERANCE_ABS;
       self.ds.shape_info.push(new_si);
     }
   }
 }
 // L181-184: update bounding box for the SD vertex (both nSD and new)
 if let Some(si_mut) = self.ds.shape_info.get_mut(n_v) {
   si_mut.box_min = Some(centroid);
   si_mut.box_max = Some(centroid);
   si_mut.box_gap = bounding_tol + TOLERANCE_ABS;
 }
 // L186-231: 3. Record SD mappings + VV interferences for every pair
 for i in 0..block.len() {
   let n1 = block[i];
   // L197: AddShapeSD(n1, nV)
   self.ds.add_shape_sd(n1, n_v);
   // L199-218: self-interfering shape warning
   let i_r1 = self.ds.rank(n1);
   // L221-228: VV interference for each pair (n1, n2)
   for j in (i + 1)..block.len() {
     let n2 = block[j];
     // OCCT L199-218: if same rank, add self-interfering shape warning
     if i_r1 == self.ds.rank(n2) {
       self.my_report.add_alert(crate::bopalgo::Alert::SelfInterferingShape(n1, n2));
     }
     // OCCT L223: if (myDS->AddInterf(n1, n2)) — fence check
     let key = if n1 < n2 { (n1, n2) } else { (n2, n1) };
     if self.ds.interf_tb.insert(key) {
       // L225-227: aVV.SetIndices(n1, n2); aVV.SetIndexNew(nV)
       self.ds.interf_vv.push(InterferenceVV {
         v1: n1,
         v2: n2,
         merged_vertex: n_v,
       });
     }
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
 // OCCT L180-184: HasInterfShapeSubShapes(nV, nE) — VV with edge endpoints
 if let Some(edge) = self.ds.edges.get(ei) {
   if self.ds.has_interf_vv(vi, edge.start_vertex)
   || self.ds.has_interf_vv(vi, edge.end_vertex) { continue; }
 }
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
 // OCCT L180-184: HasInterfShapeSubShapes(nV, nE)
 if let Some(edge) = self.ds.edges.get(ei) {
   if self.ds.has_interf_vv(vi, edge.start_vertex)
   || self.ds.has_interf_vv(vi, edge.end_vertex) { continue; }
 }
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
 // ??FillShrunkData computes shrunk ranges for each pave block.
 // If shrunk_range fails (edge too short), skip this pair entirely
 // (=OCCT BOPAlgo_PaveFiller_3: !aPB->IsSplittable() ??continue).
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

 // OCCT L529-551: EDGE-type common parts — coincident edges create EE interferences
 if hits.is_empty() {
  let b_coincident = match (&e1_curve, &e2_curve) {
  (Curve3::Line(l1), Curve3::Line(l2)) => {
   let cross = l1.direction.cross(l2.direction);
   if cross.length() > tol { false }
   else { (l2.origin - l1.origin).cross(l1.direction).length() <= tol }
  }
  _ => false,
  };
  if b_coincident {
   let mid_t = (range1[0] + range1[1]) * 0.5;
   let mid_pt = e1_curve.point_at(mid_t);
   let new_v = self.ds.add_vertex(mid_pt);
   self.ds.interf_ee.push(InterferenceEE{
    e1, e2, point: mid_pt,
    param1: range1[0], param2: range2[0],
    new_vertex: new_v,
   });
   self.ds.edge_paves[e1].push(Pave { vertex_idx: new_v, param: range1[0] });
   self.ds.edge_paves[e2].push(Pave { vertex_idx: new_v, param: range2[0] });
   modified.insert(e1);
   modified.insert(e2);
   return;
  }
  return;
 }

 //  Process each intersection result (PaveFiller_3.cxx L682-750).
 // For each valid intersection, create a new vertex and record EE interference.
 for (t1, t2, point) in hits {
 if t1 < sr1[0] || t1 > sr1[1] || t2 < sr2[0] || t2 > sr2[1] { continue; }
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
 self.my_increased_ss.insert(new_vi);

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
 // OCCT BOPAlgo_PaveFiller.cxx L376-441: RepeatIntersection
 pub(crate) fn repeat_intersection(&mut self) {
 let mut a_extra_map: HashSet<usize> = HashSet::new();
 // OCCT L382-407: iterate source vertices (0..NbSourceShapes, VERTEX only)
 // whose tolerance was increased, or whose SD root tolerance was increased.
 let a_nb_s = self.ds.nb_source_shapes;
 for vi in 0..self.ds.vertices.len() {
   let si = self.ds.vertex_shape_idx.get(vi).copied().unwrap_or(usize::MAX);
   if si >= a_nb_s {
     continue; // skip non-source vertices (OCCT L385: ShapeType check)
   }
   // L390-393: vertex directly in myIncreasedSS
   if self.my_increased_ss.contains(&vi) {
     a_extra_map.insert(vi);
     continue;
   }
   // L396-406: SD root whose tolerance was increased
   if let Some(n_vsd) = self.ds.has_shape_sd(vi) {
     if self.my_increased_ss.contains(&n_vsd) {
       a_extra_map.insert(vi);
     }
   }
 }
 // Build VV pairs: cross-operand involving extra vertices
 // OCCT L414: myIterator->IntersectExt(anExtraInterfMap) uses BVH to find
 // candidate pairs involving the extra vertices only.
 let a_vc = self.ds.a_vertex_count;
 let vv_bvh = self.build_ds_bvh_combined(false);
 let all_candidates = crate::bvh::DsBvh::candidate_pairs(&vv_bvh, &vv_bvh);
 let mut vv_pairs: Vec<(usize, usize)> = Vec::new();
 for &(vi, vj) in &all_candidates {
  let vi_x = a_extra_map.contains(&vi);
  let vj_x = a_extra_map.contains(&vj);
  // Only cross-operand pairs where exactly one vertex is "extra"
  if vi_x != vj_x && ((vi < a_vc) != (vj < a_vc)) {
   vv_pairs.push((vi, vj));
  }
 }
 self.perform_vv(&vv_pairs);
 self.ds.update_pave_blocks_with_sd_vertices();
 // Build VE pairs: cross-operand where vertex is in the extra set
 let a_ec = self.ds.a_edge_count;
 let n_edges = self.ds.edges.len();
 let mut ve_pairs: Vec<(usize, usize)> = Vec::new();
 for &vi in &a_extra_map {
 if vi < a_vc {
 for ei in a_ec..n_edges {
 ve_pairs.push((vi, ei));
 }
 } else {
 for ei in 0..a_ec {
 ve_pairs.push((vi, ei));
 }
 }
 }
 self.perform_ve_bvh(&ve_pairs);
 self.ds.update_pave_blocks_with_sd_vertices();
 // Build VF pairs: cross-operand where vertex is in the extra set
 let a_fc = self.ds.a_face_count;
 let n_faces = self.ds.faces.len();
 let mut vf_pairs: Vec<(usize, usize)> = Vec::new();
 for &vi in &a_extra_map {
 if vi < a_vc {
 for fi in a_fc..n_faces {
 vf_pairs.push((vi, fi));
 }
 } else {
 for fi in 0..a_fc {
 vf_pairs.push((vi, fi));
 }
 }
 }
 self.perform_vf_bvh(&vf_pairs);
 self.ds.update_pave_blocks_with_sd_vertices();
 self.update_interfs_with_sd_vertices();
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
            index_new: None,
        });
 if proj_dist > 0.0 && proj_dist < f64::MAX
 && proj_dist > self.ds.vertex_tolerance(n_vsd)
 {
 self.ds.vertex_data_mut(n_vsd).tolerance = proj_dist;
 self.my_increased_ss.insert(n_vsd);
 }

 //  ALL VF vertices go to VerticesIn (OCCT L297: aMVIn.Add)
 self.ds.face_info_mut(fi).vertices_in.insert(n_vsd);
 }
 }
 /// PerformEF (PaveFiller_5.cxx L165-300) ??LEGACY non-BVH path.
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
 // OCCT L416-417: analytic ??aBC.Resolution(aRes)
 Curve3::Line(l) => {
 let dir_len = l.direction.length();
 if dir_len > 1e-12 { etf / dir_len } else { etf }
 }
 Curve3::Circle(c) => etf / c.radius.max(TOLERANCE_ABS),
 Curve3::Ellipse(e) => etf / e.major_radius.max(TOLERANCE_ABS),
 // OCCT L391-413: BSpline/Bezier ??aRes / |derivative|
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


