//! Boolean operations (Union / Common / Cut) aligned with OCCT's high-level phases.
//!
//! This module provides the PaveFiller-based boolean pipeline for all three operation types,
//! mirroring OCCT's `BOPAlgo_BOP` / `BRepAlgoAPI_Fuse` / `BRepAlgoAPI_Common` / `BRepAlgoAPI_Cut`
//! control flow:
//!
//! 1. **Prepare arguments** -- build the interference descriptor from both operands
//!    ([`bopds::ds::DS::new`]), analogous to loading shapes into the BOP data structure.
//! 2. **Intersection and paving** -- [`crate::pave_filler::PaveFiller::perform`], analogous to
//!    `BOPAlgo_PaveFiller` (edge/face interferences, splits, pave sets).
//! 3. **Build result** -- [`crate::builder::BooleanBuilder`] with the requested
//!    [`BooleanOpType`], analogous to building the result solid from classified pieces.
//! 4. **Post-process** (serial `fuse` only) -- [`crate::geom_populate::recompute_plane_surfaces`]
//!    then an iterated [`unify_same_domain_faces`](crate::unify_same_domain_faces) /
//!    [`crate::unify_same_domain_faces`] pass to merge coplanar box fragments (OCCT
//!    `UnifySameDomain` + same-domain analog).
//!
//! History and parallel-history APIs intentionally skip step 4 to match the existing behavior
//! of [`boolean_op_with_history_generic`] and [`crate::boolean_op_par`].
//!
//! ## OCCT mapping
//!
//! | Rust                                | OCCT                                     | Align |
//! |-------------------------------------|-----------------------------------------|-------|
//! | `boolean_op_generic`                | `BRepAlgoAPI_Fuse/Common/Cut::Shape()`  |   |
//! | `boolean_op_with_history_generic`   | `BRepAlgoAPI_BuilderOperation::Build()` |   |
//! | `fuse` / `fuse_with_history`        | `BRepAlgoAPI_Fuse` shortcut             |   |
//! | `common_with_history`               | `BRepAlgoAPI_Common` shortcut           |   |
//! | `cut_with_history`                  | `BRepAlgoAPI_Cut` shortcut              |   |
//! | `validate_ds_invariants`            | rcad-specific (no direct OCCT eq)       |   |
//!
//! ## In-pipeline validation
//!
//! Each phase runs consistency checks before the next: operand [`BRep`] pool indices and
//! finite coordinates, full [`bopds::ds::DS`] internal consistency, then the same index/finite
//! checks on the result (and again after plane recompute on the serial path). These mirror
//! invariants the implementation is expected to maintain; they are intentionally weaker than
//! a full [`crate::brep_check::check`] pass. Failures surface as [`BooleanError::InvalidResult`],
//! [`BooleanError::EmptyInput`], or [`BooleanError::NumericalFailure`] with message prefix `bop:`.

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

/// DSU (Union-Find) for building connected SameDomain groups -- equivalent to OCCT FillMap + MakeBlocks.
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
/// shell here -- downstream feature code may pass intermediate shapes that are not yet
/// `brep_check`-clean.
/// Invariants on the BOP DS after prepare ([`DS::new`]) and after paving ([`PaveFiller::perform]).
fn validate_ds_invariants(ds: &DS) -> Result<(), BooleanError> {
 let nv = ds.vertices.len();
 let ne = ds.edges.len();
 let nf = ds.faces.len();
 let nic = ds.intersection_curves.len();

 if !ds.fuzzy_tol.is_finite() || ds.fuzzy_tol < 0.0 {
 return Err(BooleanError::NumericalFailure("bop: DS fuzzy_tol invalid"));
 }
 if ds.a_vertex_count > nv || ds.a_edge_count > ne || ds.a_face_count > nf {
 return Err(BooleanError::InvalidResult(
 "bop: DS shape A extent vs pool size inconsistent",
 ));
 }

 for v in &ds.vertices {
 let p = v.point;
 if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
 return Err(BooleanError::NumericalFailure(
 "bop: DS vertex coordinate non-finite",
 ));
 }
 }

 for e in &ds.edges {
 if e.start_vertex >= nv || e.end_vertex >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: DS edge references vertex out of range",
 ));
 }
 let [t0, t1] = e.t_range;
 if !(t0.is_finite() && t1.is_finite()) {
 // Degenerate edges (start==end) on periodic surfaces like sphere
 // poles may have NaN t_range (zero 3D length). Skip t_range check.
 if e.start_vertex == e.end_vertex { continue; }
 return Err(BooleanError::NumericalFailure(
 "bop: DS edge t_range non-finite",
 ));
 }
 for p in &e.paves {
 if p.vertex_idx >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: DS pave vertex index out of range",
 ));
 }
 if !p.param.is_finite() {
 return Err(BooleanError::NumericalFailure(
 "bop: DS pave param non-finite",
 ));
 }
 }
 for pb in &e.pave_blocks {
 if pb.0.read().unwrap().original_edge != NO_EDGE && pb.0.read().unwrap().original_edge >= ne {
 return Err(BooleanError::InvalidResult(
 "bop: DS pave_block original_edge out of range",
 ));
 }
 if pb.0.read().unwrap().pave1.vertex_idx >= nv || pb.0.read().unwrap().pave2.vertex_idx >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: DS pave_block vertex out of range",
 ));
 }
 if let Some(ni) = pb.0.read().unwrap().new_edge
 && ni >= ne {
 return Err(BooleanError::InvalidResult(
 "bop: DS pave_block new_edge out of range",
 ));
 }
 }
 }

 for f in &ds.faces {
 for &vi in &f.boundary_verts {
 if vi >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: DS face boundary_verts out of range",
 ));
 }
 }
 for &ei in &f.boundary_edges {
 if ei >= ne {
 return Err(BooleanError::InvalidResult(
 "bop: DS face boundary_edges out of range",
 ));
 }
 }
 let n = f.normal;
 // Paving may still refine orientation; only reject NaNs / infinities here.
 if !n.x.is_finite() || !n.y.is_finite() || !n.z.is_finite() {
 return Err(BooleanError::InvalidResult(
 "bop: DS face normal non-finite",
 ));
 }
 for &v in &f.face_info.vertices_on {
 if v >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: DS face_info.vertices_on out of range",
 ));
 }
 }
 for &v in &f.face_info.vertices_in {
 if v >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: DS face_info.vertices_in out of range",
 ));
 }
 }
 for &ci in &f.face_info.curves_sc {
 if ci >= nic {
 return Err(BooleanError::InvalidResult(
 "bop: DS face_info.curves_sc out of range",
 ));
 }
 }
 }

 for ic in &ds.intersection_curves {
 if ic.start_vertex >= nv || ic.end_vertex >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: DS intersection_curve vertex out of range",
 ));
 }
 for p in &ic.polyline {
 if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
 return Err(BooleanError::NumericalFailure(
 "bop: DS intersection polyline non-finite",
 ));
 }
 }
 }

 for inf in &ds.interf_vv {
 if inf.v1 >= nv || inf.v2 >= nv || inf.merged_vertex >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: interference VertexVertex index out of range",
 ));
 }
 }
 for inf in &ds.interf_ve {
 if inf.vertex >= nv || inf.edge >= ne || !inf.param.is_finite() {
 return Err(BooleanError::InvalidResult(
 "bop: interference VertexEdge index/param invalid",
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
 "bop: interference EdgeEdge invalid",
 ));
 }
 }
 for inf in &ds.interf_vf {
 if inf.vertex >= nv || inf.face >= nf {
 return Err(BooleanError::InvalidResult(
 "bop: interference VertexFace index out of range",
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
 "bop: interference EdgeFace invalid",
 ));
 }
 }
 for inf in &ds.interf_ff {
 if inf.f1 >= nf || inf.f2 >= nf {
 return Err(BooleanError::InvalidResult(
 "bop: interference FaceFace face index out of range",
 ));
 }
 for &c in &inf.curves {
 if c >= nic {
 return Err(BooleanError::InvalidResult(
 "bop: interference FaceFace curve index out of range",
 ));
 }
 }
 for pv in &inf.points {
 if pv.vertex_index != usize::MAX && pv.vertex_index >= nv {
 return Err(BooleanError::InvalidResult(
 "bop: interference FaceFace point vertex out of range",
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
 bvh::Bvh::build(brep)
}

/// --PaveFiller creation + configuration + Perform.
/// OCCT BOPAlgo_BOP::Perform L395-405: new PaveFiller + config + Perform.
/// DS is created empty and populated by PaveFiller::init()
/// inside perform(), mirroring BOPAlgo_PaveFiller::Init.
pub fn pave_fill(a: &topods::BRep, b: &topods::BRep, use_bvh: bool,
 brep: &mut rcad_kernel::topods::BRep,
 fuzzy_tol: f64) -> (bopds::ds::DS,
 Vec<rcad_kernel::topods::ShapeRef>, Vec<Option<rcad_kernel::topods::ShapeRef>>)
{
 let (bvh_a, bvh_b) = if use_bvh { optional_bvhs(a, b) } else { (None, None) };
 let mut ds = bopds::ds::DS::new_empty();
 {
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(brep);
 f
 }
 };
 filler.set_run_parallel(false);
 filler.configure_fuzzy(fuzzy_tol);
 filler.set_non_destructive(false);
 filler.configure_glue(false, TOLERANCE_ABS);
 filler.set_use_obb(false);
 filler.perform(a, b);
 // extract outputs before filler drops (ends borrow on ds)
 let face_refs = std::mem::take(&mut filler.face_refs);
 let ic_edge_map = std::mem::take(&mut filler.ic_edge_map);
 (ds, face_refs, ic_edge_map)
 }
}

/// Sum of boundary-edge counts from [`crate::brep_check::validate_solid_closure`].
/// Larger means a **less** closed manifold shell.
/// Uses DS -- PaveFiller -- BooleanBuilder(op) pipeline.
///
/// Same phases as [`boolean_op_generic`], returns history.
///
/// `BRepAlgoAPI_BuilderOperation::Build()`
pub(crate) fn boolean_op_with_history_generic(
 op: BooleanOpType, a: &topods::BRep, b: &topods::BRep,
) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, true, &mut brep, crate::tolerance::TOLERANCE_ABS);
 let mut builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}

/// Generic boolean operation for Union / Common / Cut.
///
/// `BRepAlgoAPI_Fuse/Common/Cut::Shape()`
/// Uses BVH when both operands have faces.
/// Generic boolean operation for Union / Common / Cut.
///
/// `BRepAlgoAPI_Fuse/Common/Cut::Shape()`
/// Uses BVH when both operands have faces.
pub fn boolean_op_generic(
 op: BooleanOpType, a: &topods::BRep, b: &topods::BRep,
) -> Result<topods::BRep, BooleanError> {
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, true, &mut brep, crate::tolerance::TOLERANCE_ABS);
 let mut builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
 let (result, _history) = builder.build_with_history_topods()?;
 Ok(result)
}

// ── Union (Fuse) shortcuts ─────────────────────────────────────────────

/// OCCT shortcut: `BRepAlgoAPI_Fuse(a, b).Shape()`
pub(crate) fn fuse(a: &topods::BRep, b: &topods::BRep) -> Result<topods::BRep, BooleanError> {
 fuse_with_bvh(a, b, true)
}

/// --DS -- PaveFiller -- BooleanBuilder(Union) -- result.
/// OCCT BOPAlgo_BOP::Perform L395-408: PaveFiller config + Perform + PerformInternal1.
pub(crate) fn fuse_with_bvh(a: &topods::BRep, b: &topods::BRep, use_bvh: bool) -> Result<topods::BRep, BooleanError> {
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, use_bvh, &mut brep, crate::tolerance::TOLERANCE_ABS);
 let mut builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Union, brep, face_refs, ic_edge_map);
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
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, use_bvh, &mut brep, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Union, brep, face_refs, ic_edge_map);
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
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, use_bvh, &mut brep, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Union, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}

// ── Common (Intersection) shortcuts ─────────────────────────────────────

/// OCCT shortcut: `BRepAlgoAPI_Common(a, b).Shape()` with history.
pub(crate) fn common_with_history(a: &topods::BRep, b: &topods::BRep) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 common_with_history_bvh(a, b, true)
}

pub(crate) fn common_with_history_bvh(
 a: &topods::BRep,
 b: &topods::BRep,
 use_bvh: bool,
) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, use_bvh, &mut brep, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Intersection, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}

/// Parallel classification path for Common (Intersection).
pub(crate) fn common_with_history_par(a: &topods::BRep, b: &topods::BRep) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 common_with_history_par_bvh(a, b, true)
}

pub(crate) fn common_with_history_par_bvh(
 a: &topods::BRep,
 b: &topods::BRep,
 use_bvh: bool,
) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, use_bvh, &mut brep, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Intersection, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}

// ── Cut (Difference) shortcuts ──────────────────────────────────────────

/// OCCT shortcut: `BRepAlgoAPI_Cut(a, b).Shape()` with history.
pub(crate) fn cut_with_history(a: &topods::BRep, b: &topods::BRep) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 cut_with_history_bvh(a, b, true)
}

pub(crate) fn cut_with_history_bvh(
 a: &topods::BRep,
 b: &topods::BRep,
 use_bvh: bool,
) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, use_bvh, &mut brep, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Difference, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}

/// Parallel classification path for Cut (Difference).
pub(crate) fn cut_with_history_par(a: &topods::BRep, b: &topods::BRep) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 cut_with_history_par_bvh(a, b, true)
}

pub(crate) fn cut_with_history_par_bvh(
 a: &topods::BRep,
 b: &topods::BRep,
 use_bvh: bool,
) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 let mut brep = rcad_kernel::topods::BRep::new();
 let (mut ds, face_refs, ic_edge_map) = pave_fill(a, b, use_bvh, &mut brep, crate::tolerance::TOLERANCE_ABS);
 validate_ds_invariants(&ds)?;
 let mut builder = builder::BooleanBuilder::with_brep(&ds, BooleanOpType::Difference, brep, face_refs, ic_edge_map);
 let (result_brep, hist) = builder.build_with_history()?;
 Ok((result_brep, hist))
}


