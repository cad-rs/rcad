//! Boolean **Union** (fuse) aligned with Open CASCADEâ€™s high-level phases.
//!
//! This module is the dedicated entry for [`crate::BooleanOpType::Union`]. The steps mirror
//! OCCTâ€™s `BOPAlgo_Builder` / `BRepAlgoAPI_Fuse` control flow:
//!
//! 1. **Prepare arguments** â€?build the interference descriptor from both operands
//! ([`bopds::ds::DS::new`]), analogous to loading shapes into the BOP data structure.
//! 2. **Intersection and paving** â€?[`crate::pave_filler::PaveFiller::perform`], analogous to
//! `BOPAlgo_PaveFiller` (edge/face interferences, splits, pave sets).
//! 3. **Build fuse result** â€?[`crate::builder::BooleanBuilder`] with
//! [`crate::BooleanOpType::Union`], analogous to building the fused solid from classified
//! pieces.
//! 4. **Post-process** (serial `fuse` only) â€?[`crate::geom_populate::recompute_plane_surfaces`]
//! then an iterated [`unify_same_domain_faces`](crate::unify_same_domain_faces) /
//! [`crate::unify_same_domain_faces`] pass to merge coplanar box fragments (OCCT
//! `UnifySameDomain` + same-domain analog; wire order can still leave analytic area ~5â€?0% low
//! on some merged planes until the kernelâ€™s planar integrator is tightened).
//!
//! History and parallel-history APIs intentionally skip step 4 to match the existing behavior
//! of [`crate::boolean_op_with_history`] and [`crate::boolean_op_par`].
//!
//! ## In-pipeline validation
//!
//! Each phase runs consistency checks before the next: operand [`BRep`] pool indices and
//! finite coordinates, full [`bopds::ds::DS`] internal consistency, then the same index/finite
//! checks on the result (and again after plane recompute on the serial path). These mirror
//! invariants the implementation is expected to maintain; they are intentionally weaker than
//! a full [`crate::brep_check::check`] pass. Failures surface as [`BooleanError::InvalidResult`],
//! [`BooleanError::EmptyInput`], or [`BooleanError::NumericalFailure`] with message prefix `union:`.

use glam::{DVec2, DVec3};
use crate::bopds;
use crate::bopds::ds::{DS};
use crate::bopds::pave::NO_EDGE;
use crate::builder;
use crate::bvh;
use crate::geom_populate;
use crate::history::{BooleanHistory, FaceOrigin};
use crate::pave_filler;
use crate::tolerance::*;
use crate::total_surface_area;
use crate::BooleanError;
use crate::BooleanOpType;
use rcad_kernel::geom::Surface3;
use rcad_kernel::topods;
use crate::bvh::{Bvh, Aabb};

/// DSU (Union-Find) for building connected SameDomain groups â€?equivalent to OCCT FillMap + MakeBlocks.
struct DSU {
 parent: Vec<usize>,
 rank: Vec<u32>,
}
impl DSU {
 fn new(n: usize) -> Self {
 Self { parent: (0..n).collect(), rank: vec![0; n] }
 }
 fn find(&mut self, x: usize) -> usize {
 if self.parent[x] != x { self.parent[x] = self.find(self.parent[x]); }
 self.parent[x]
 }
 fn union(&mut self, a: usize, b: usize) {
 let ra = self.find(a);
 let rb = self.find(b);
 if ra != rb {
 if self.rank[ra] < self.rank[rb] { self.parent[ra] = rb; }
 else if self.rank[ra] > self.rank[rb] { self.parent[rb] = ra; }
 else { self.parent[rb] = ra; self.rank[ra] += 1; }
 }
 }
}

/// Operands must be usable as boolean arguments: non-empty face set and **index-consistent**
/// topology (plus finite vertex coordinates). We intentionally do **not** require a watertight
/// shell here â€?downstream feature code may pass intermediate shapes that are not yet
/// `brep_check`â€‘clean.
/// Invariants on the BOP DS after prepare ([`DS::new`]) and after paving ([`PaveFiller::perform`]).
fn validate_ds_invariants(ds: &DS) -> Result<(), BooleanError> {
 let nv = ds.vertices.len();
 let ne = ds.edges.len();
 let nf = ds.faces.len();
 let nic = ds.intersection_curves.len();

 if !ds.fuzzy_tol.is_finite() || ds.fuzzy_tol < 0.0 {
 return Err(BooleanError::NumericalFailure("union: DS fuzzy_tol invalid"));
 }
 if ds.a_vertex_count > nv || ds.a_edge_count > ne || ds.a_face_count > nf {
 return Err(BooleanError::InvalidResult(
 "union: DS shape A extent vs pool size inconsistent",
 ));
 }

 for v in &ds.vertices {
 let p = v.point;
 if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
 return Err(BooleanError::NumericalFailure(
 "union: DS vertex coordinate non-finite",
 ));
 }
 }

 for e in &ds.edges {
 if e.start_vertex >= nv || e.end_vertex >= nv {
 return Err(BooleanError::InvalidResult(
 "union: DS edge references vertex out of range",
 ));
 }
 let [t0, t1] = e.t_range;
 if !(t0.is_finite() && t1.is_finite()) {
 // Degenerate edges (start==end) on periodic surfaces like sphere
 // poles may have NaN t_range (zero 3D length). Skip t_range check.
 if e.start_vertex == e.end_vertex { continue; }
 return Err(BooleanError::NumericalFailure(
 "union: DS edge t_range non-finite",
 ));
 }
 for p in &e.paves {
 if p.vertex_idx >= nv {
 return Err(BooleanError::InvalidResult(
 "union: DS pave vertex index out of range",
 ));
 }
 if !p.param.is_finite() {
 return Err(BooleanError::NumericalFailure(
 "union: DS pave param non-finite",
 ));
 }
 }
 for pb in &e.pave_blocks {
 if pb.0.read().unwrap().original_edge != NO_EDGE && pb.0.read().unwrap().original_edge >= ne {
 return Err(BooleanError::InvalidResult(
 "union: DS pave_block original_edge out of range",
 ));
 }
 if pb.0.read().unwrap().pave1.vertex_idx >= nv || pb.0.read().unwrap().pave2.vertex_idx >= nv {
 return Err(BooleanError::InvalidResult(
 "union: DS pave_block vertex out of range",
 ));
 }
 if let Some(ni) = pb.0.read().unwrap().new_edge
 && ni >= ne {
 return Err(BooleanError::InvalidResult(
 "union: DS pave_block new_edge out of range",
 ));
 }
 }
 }

 for f in &ds.faces {
 for &vi in &f.boundary_verts {
 if vi >= nv {
 return Err(BooleanError::InvalidResult(
 "union: DS face boundary_verts out of range",
 ));
 }
 }
 for &ei in &f.boundary_edges {
 if ei >= ne {
 return Err(BooleanError::InvalidResult(
 "union: DS face boundary_edges out of range",
 ));
 }
 }
 let n = f.normal;
 // Paving may still refine orientation; only reject NaNs / infinities here.
 if !n.x.is_finite() || !n.y.is_finite() || !n.z.is_finite() {
 return Err(BooleanError::InvalidResult(
 "union: DS face normal non-finite",
 ));
 }
 for &v in &f.face_info.vertices_on {
 if v >= nv {
 return Err(BooleanError::InvalidResult(
 "union: DS face_info.vertices_on out of range",
 ));
 }
 }
 for &v in &f.face_info.vertices_in {
 if v >= nv {
 return Err(BooleanError::InvalidResult(
 "union: DS face_info.vertices_in out of range",
 ));
 }
 }
 for &ci in &f.face_info.curves_sc {
 if ci >= nic {
 return Err(BooleanError::InvalidResult(
 "union: DS face_info.curves_sc out of range",
 ));
 }
 }
 }

 for ic in &ds.intersection_curves {
 if ic.start_vertex >= nv || ic.end_vertex >= nv {
 return Err(BooleanError::InvalidResult(
 "union: DS intersection_curve vertex out of range",
 ));
 }
 for p in &ic.polyline {
 if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
 return Err(BooleanError::NumericalFailure(
 "union: DS intersection polyline non-finite",
 ));
 }
 }
 }

 for inf in &ds.interf_vv {
 if inf.v1 >= nv || inf.v2 >= nv || inf.merged_vertex >= nv {
 return Err(BooleanError::InvalidResult(
 "union: interference VertexVertex index out of range",
 ));
 }
 }
 for inf in &ds.interf_ve {
 if inf.vertex >= nv || inf.edge >= ne || !inf.param.is_finite() {
 return Err(BooleanError::InvalidResult(
 "union: interference VertexEdge index/param invalid",
 ));
 }
 }
 for inf in &ds.interf_ee {
 if inf.e1 >= ne
 || inf.e2 >= ne
 || inf.new_vertex >= nv
 || !inf.point.x.is_finite()
 || !inf.point.y.is_finite()
 || !inf.point.z.is_finite()
 || !inf.param1.is_finite()
 || !inf.param2.is_finite()
 {
 return Err(BooleanError::InvalidResult(
 "union: interference EdgeEdge invalid",
 ));
 }
 }
 for inf in &ds.interf_vf {
 if inf.vertex >= nv || inf.face >= nf {
 return Err(BooleanError::InvalidResult(
 "union: interference VertexFace index out of range",
 ));
 }
 }
 for inf in &ds.interf_ef {
 if inf.edge >= ne
 || inf.face >= nf
 || inf.new_vertex >= nv
 || !inf.point.x.is_finite()
 || !inf.point.y.is_finite()
 || !inf.point.z.is_finite()
 || !inf.edge_param.is_finite()
 {
 return Err(BooleanError::InvalidResult(
 "union: interference EdgeFace invalid",
 ));
 }
 }
 for inf in &ds.interf_ff {
 if inf.f1 >= nf || inf.f2 >= nf {
 return Err(BooleanError::InvalidResult(
 "union: interference FaceFace face index out of range",
 ));
 }
 for &c in &inf.curves {
 if c >= nic {
 return Err(BooleanError::InvalidResult(
 "union: interference FaceFace curve index out of range",
 ));
 }
 }
 for &pv in &inf.points {
 if pv >= nv {
 return Err(BooleanError::InvalidResult(
 "union: interference FaceFace point vertex out of range",
 ));
 }
 }
 }

 Ok(())
}

fn optional_bvhs(a: &topods::BRep, b: &topods::BRep) -> (Option<bvh::Bvh>, Option<bvh::Bvh>) {
 let has_faces = |brep: &topods::BRep| -> bool {
 brep.tshapes.iter().any(|ts| matches!(ts.as_ref(), topods::TShape::Face(_)))
 };
 (
 if has_faces(a) { Some(build_topods_bvh(a)) } else { None },
 if has_faces(b) { Some(build_topods_bvh(b)) } else { None },
 )
}

/// Build a simple BVH from a topods::BRep by computing face AABBs from wire vertices.
fn build_topods_bvh(brep: &topods::BRep) -> bvh::Bvh {
 use crate::tolerance::TOLERANCE_LINEAR_ULTRA_STRICT;
 let mut face_aabbs = Vec::new();

 for ts in &brep.tshapes {
 if let topods::TShape::Face(fd) = ts.as_ref() {
 let mut aabb = Aabb::empty();
 let wire = &brep.wire(fd.outer_wire);
 for edge_sr in &wire.edges {
 if let topods::TShape::Edge(ed) = brep.tshapes[edge_sr.index].as_ref() {
 let vd = brep.vertex(ed.first);
 aabb.expand_point(vd.point);
 let vd2 = brep.vertex(ed.last);
 aabb.expand_point(vd2.point);
 }
 }
 // Expand for surface type
 if let Some(ref surf) = fd.surface {
 match surf {
 Surface3::Sphere(s) => {
 let r = s.radius.abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
 aabb.expand_point(s.center - glam::DVec3::splat(r));
 aabb.expand_point(s.center + glam::DVec3::splat(r));
 }
 _ => {}
 }
 }
 face_aabbs.push(aabb);
 }
 }

 let indices: Vec<usize> = (0..face_aabbs.len()).collect();
 bvh::Bvh::build(indices, face_aabbs)
}
 let has_faces_a = a
 .solids
 .first()
 .and_then(|s| s.shells.first())
 .is_some_and(|sh| !sh.faces.is_empty());
 let has_faces_b = b
 .solids
 .first()
 .and_then(|s| s.shells.first())
 .is_some_and(|sh| !sh.faces.is_empty());
 (
 if has_faces_a {
 Some(bvh::Bvh::build(a))
 } else {
 None
 },
 if has_faces_b {
 Some(bvh::Bvh::build(b))
 } else {
 None
 },
 )
}

/// `use_bvh`: match [`crate::boolean_op`] (`true`) or [`crate::brep_algo_api::BRepAlgoAPI_Fuse`]
/// when BVH acceleration is toggled off (`false` â†?plain [`pave_filler::PaveFiller::new`]).
/// âœ?OCCT-aligned: PaveFiller creation + configuration + Perform.
/// OCCT BOPAlgo_BOP::Perform L395-405: new PaveFiller + config + Perform.
fn pave_fill(ds: &mut bopds::ds::DS, a: &topods::BRep, b: &topods::BRep, use_bvh: bool,
 brep: &mut rcad_kernel::topods::BRep) -> (Vec<rcad_kernel::topods::ShapeRef>, Vec<Option<rcad_kernel::topods::ShapeRef>>)
{
 let (bvh_a, bvh_b) = if use_bvh { optional_bvhs(a, b) } else { (None, None) };
 let fuzzy_tol = ds.fuzzy_tol;
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(ds, ba, bb, brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(ds);
 f.brep = Some(brep);
 f
 }
 };
 filler.set_run_parallel(false);
 filler.configure_fuzzy(fuzzy_tol);
 filler.set_non_destructive(false);
 filler.configure_glue(false, TOLERANCE_ABS);
 filler.set_use_obb(false);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
}

/// Sum of boundary-edge counts from [`crate::brep_check::validate_solid_closure`].
/// Larger means a **less** closed manifold shell.
/// Union: DS â†?PaveFiller â†?BooleanBuilder(Union) â†?recompute plane surfaces.
///
/// Uses BVH when both operands have faces, matching [`crate::boolean_op`].
pub(crate) fn fuse(a: &topods::BRep, b: &topods::BRep) -> Result<topods::BRep, BooleanError> {
 fuse_with_bvh(a, b, true)
}

/// âœ?OCCT-aligned: DS â†?PaveFiller â†?BooleanBuilder(Union) â†?result.
/// OCCT BOPAlgo_BOP::Perform L395-408: PaveFiller config + Perform + PerformInternal1.
pub(crate) fn fuse_with_bvh(a: &topods::BRep, b: &topods::BRep, use_bvh: bool) -> Result<topods::BRep, BooleanError> {
 let mut ds = bopds::ds::DS::new_from_topods(a, b, crate::tolerance::TOLERANCE_ABS);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = pave_fill(&mut ds, a, b, use_bvh, &mut brep);
 let builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Union, brep, face_refs, ic_edge_map);
 let (result, _history) = builder.build_with_history_topods()?;
 Ok(result)
}

/// Same phases as [`fuse`], but returns [`BooleanHistory`] and does not run plane recompute
/// (matches legacy `boolean_op_with_history` for Union).
pub(crate) fn fuse_with_history(a: &topods::BRep, b: &topods::BRep) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 fuse_with_history_bvh(a, b, true)
}

pub(crate) fn fuse_with_history_bvh(
 a: &topods::BRep,
 b: &topods::BRep,
 use_bvh: bool,
) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 let mut ds = bopds::ds::DS::new_from_topods(a, b, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = pave_fill(&mut ds, a, b, use_bvh, &mut brep);
 validate_ds_invariants(&ds)?;
 let builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Union, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}

/// Parallel classification path; same OCCT phase structure as [`fuse_with_history`].
pub(crate) fn fuse_with_history_par(a: &topods::BRep, b: &topods::BRep) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 fuse_with_history_par_bvh(a, b, true)
}

pub(crate) fn fuse_with_history_par_bvh(
 a: &topods::BRep,
 b: &topods::BRep,
 use_bvh: bool,
) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 let mut ds = bopds::ds::DS::new_from_topods(a, b, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = pave_fill(&mut ds, a, b, use_bvh, &mut brep);
 validate_ds_invariants(&ds)?;
 let builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Union, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}

/// âœ?OCCT : edge set (BOPTools_Set)  ã€?
/// OCCT  : BOPAlgo_Builder_2.cxx L571-L832 (FillSameDomainFaces)
///
/// (  OCCT):
/// 1. edge key set ( , )
/// 2. edge key set â†?anESetFaces
/// 3.  â‰? ,  SameDomain:
/// a.  (Plane/planar BSpline) â†?(OCCT L697-701:  )
/// b.  (Cylinder+Cylinder  ) â†?(OCCT L703-708)
/// c. â†?TODO: AreFacesSameDomain (OCCT L703-708)
/// 4. MakeBlocks: Faceâ†’Face (OCCT L741: BOPAlgo_Tools::MakeBlocks)
/// 5.  :  ( DS ) >  ,  flat index
/// (OCCT L758-788: nFMin â†?myDS->Index(aF) >= 0  )
/// 6.  ,  geom slots
///
///  ( ) ã€?
#
