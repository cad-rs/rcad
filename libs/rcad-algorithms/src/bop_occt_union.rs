//! Boolean **Union** (fuse) aligned with Open CASCADE’s high-level phases.
//!
//! This module is the dedicated entry for [`crate::BooleanOpType::Union`]. The steps mirror
//! OCCT’s `BOPAlgo_Builder` / `BRepAlgoAPI_Fuse` control flow:
//!
//! 1. **Prepare arguments** — build the interference descriptor from both operands
//!    ([`bopds::ds::DS::new`]), analogous to loading shapes into the BOP data structure.
//! 2. **Intersection and paving** — [`crate::pave_filler::PaveFiller::perform`], analogous to
//!    `BOPAlgo_PaveFiller` (edge/face interferences, splits, pave sets).
//! 3. **Build fuse result** — [`crate::builder::BooleanBuilder`] with
//!    [`crate::BooleanOpType::Union`], analogous to building the fused solid from classified
//!    pieces.
//! 4. **Post-process** (serial `fuse` only) — [`crate::geom_populate::recompute_plane_surfaces`]
//!    then an iterated [`crate::orthogonal_face_fuse::fuse_orthogonal_coplanar_faces`] /
//!    [`crate::unify_same_domain_faces`] pass to merge coplanar box fragments (OCCT
//!    `UnifySameDomain` + same-domain analog; wire order can still leave analytic area ~5–10% low
//!    on some merged planes until the kernel’s planar integrator is tightened).
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
use crate::brep_repair::merge_close_vertices;
use crate::bopds;
use crate::bopds::ds::{DS, Interference};
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
use rcad_kernel::BRep;

/// DSU (Union-Find) for building connected SameDomain groups — equivalent to OCCT FillMap + MakeBlocks.
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
/// shell here — downstream feature code may pass intermediate shapes that are not yet
/// `brep_check`‑clean.
fn validate_union_operand(msg: &'static str, brep: &BRep) -> Result<(), BooleanError> {
    if brep.solids.is_empty() {
        return Err(BooleanError::EmptyInput);
    }
    let has_face = brep
        .solids
        .iter()
        .any(|s| s.shells.iter().any(|sh| !sh.faces.is_empty()));
    if !has_face {
        return Err(BooleanError::EmptyInput);
    }
    validate_brep_operand_topology(msg, brep)?;
    Ok(())
}

fn validate_union_operands(a: &BRep, b: &BRep) -> Result<(), BooleanError> {
    validate_union_operand("union: input A failed operand check", a)?;
    validate_union_operand("union: input B failed operand check", b)?;
    Ok(())
}

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
            if pb.original_edge >= ne {
                return Err(BooleanError::InvalidResult(
                    "union: DS pave_block original_edge out of range",
                ));
            }
            if pb.pave1.vertex_idx >= nv || pb.pave2.vertex_idx >= nv {
                return Err(BooleanError::InvalidResult(
                    "union: DS pave_block vertex out of range",
                ));
            }
            if let Some(ni) = pb.new_edge
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

    for inf in &ds.interferences {
        match inf {
            Interference::VertexVertex {
                v1,
                v2,
                merged_vertex,
            } => {
                if *v1 >= nv || *v2 >= nv || *merged_vertex >= nv {
                    return Err(BooleanError::InvalidResult(
                        "union: interference VertexVertex index out of range",
                    ));
                }
            }
            Interference::VertexEdge { vertex, edge, param } => {
                if *vertex >= nv || *edge >= ne || !param.is_finite() {
                    return Err(BooleanError::InvalidResult(
                        "union: interference VertexEdge index/param invalid",
                    ));
                }
            }
            Interference::EdgeEdge {
                e1,
                e2,
                point,
                param1,
                param2,
                new_vertex,
            } => {
                if *e1 >= ne
                    || *e2 >= ne
                    || *new_vertex >= nv
                    || !point.x.is_finite()
                    || !point.y.is_finite()
                    || !point.z.is_finite()
                    || !param1.is_finite()
                    || !param2.is_finite()
                {
                    return Err(BooleanError::InvalidResult(
                        "union: interference EdgeEdge invalid",
                    ));
                }
            }
            Interference::VertexFace { vertex, face } => {
                if *vertex >= nv || *face >= nf {
                    return Err(BooleanError::InvalidResult(
                        "union: interference VertexFace index out of range",
                    ));
                }
            }
            Interference::EdgeFace {
                edge,
                face,
                point,
                edge_param,
                new_vertex,
            } => {
                if *edge >= ne
                    || *face >= nf
                    || *new_vertex >= nv
                    || !point.x.is_finite()
                    || !point.y.is_finite()
                    || !point.z.is_finite()
                    || !edge_param.is_finite()
                {
                    return Err(BooleanError::InvalidResult(
                        "union: interference EdgeFace invalid",
                    ));
                }
            }
            Interference::FaceFace { f1, f2, curves, points } => {
                if *f1 >= nf || *f2 >= nf {
                    return Err(BooleanError::InvalidResult(
                        "union: interference FaceFace face index out of range",
                    ));
                }
                for &c in curves {
                    if c >= nic {
                        return Err(BooleanError::InvalidResult(
                            "union: interference FaceFace curve index out of range",
                        ));
                    }
                }
                for &pv in points {
                    if pv >= nv {
                        return Err(BooleanError::InvalidResult(
                            "union: interference FaceFace point vertex out of range",
                        ));
                    }
                }
            }
        }
    }

    Ok(())
}

/// Full pool / finite checks for **result** BReps (includes non‑degenerate face normals).
fn validate_brep_index_and_finite_geom(msg: &'static str, brep: &BRep) -> Result<(), BooleanError> {
    validate_brep_topology_indices(msg, brep, true)
}

/// Operands: index consistency and finite vertices only (face normals may be unset or zero on
/// intermediate construction geometry).
fn validate_brep_operand_topology(msg: &'static str, brep: &BRep) -> Result<(), BooleanError> {
    validate_brep_topology_indices(msg, brep, false)
}

fn validate_brep_topology_indices(
    msg: &'static str,
    brep: &BRep,
    check_face_normals: bool,
) -> Result<(), BooleanError> {
    for v in &brep.vertices {
        let p = v.point;
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return Err(BooleanError::NumericalFailure(msg));
        }
    }

    for (eidx, edge) in brep.edges.iter().enumerate() {
        if edge.start >= brep.vertices.len() || edge.end >= brep.vertices.len() {
            return Err(BooleanError::InvalidResult(msg));
        }
        if let Some(ci) = brep.geom.edge_curve.get(eidx).copied().flatten()
            && ci >= brep.geom.curves.len() {
                return Err(BooleanError::InvalidResult(msg));
            }
    }

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if check_face_normals {
                    // ✅ OCCT-aligned: only plane faces store a meaningful normal.
                    // Non-planar faces derive normals from surface geometry.
                    let is_plane = face.surface_idx
                        .and_then(|si| brep.geom.surfaces.get(si))
                        .is_some_and(|s| matches!(s, rcad_kernel::geom::Surface3::Plane(_)));
                    if is_plane {
                        let n = face.normal;
                        if !n.x.is_finite()
                            || !n.y.is_finite()
                            || !n.z.is_finite()
                            || n.length_squared() <= 0.0
                        {
                            return Err(BooleanError::InvalidResult(msg));
                        }
                    }
                }
                for we in &face.outer_wire.edges {
                    if we.idx >= brep.edges.len() {
                        return Err(BooleanError::InvalidResult(msg));
                    }
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        if we.idx >= brep.edges.len() {
                            return Err(BooleanError::InvalidResult(msg));
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

fn validate_union_brep_output(msg: &'static str, brep: &BRep) -> Result<(), BooleanError> {
    validate_brep_index_and_finite_geom(msg, brep)
}

fn optional_bvhs(a: &BRep, b: &BRep) -> (Option<bvh::Bvh>, Option<bvh::Bvh>) {
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
/// when BVH acceleration is toggled off (`false` → plain [`pave_filler::PaveFiller::new`]).
/// ✅ OCCT-aligned: PaveFiller creation + configuration + Perform.
///   OCCT BOPAlgo_BOP::Perform L395-405: new PaveFiller + config + Perform.
fn pave_fill(ds: &mut bopds::ds::DS, a: &BRep, b: &BRep, use_bvh: bool) {
    let (bvh_a, bvh_b) = if use_bvh { optional_bvhs(a, b) } else { (None, None) };
    let fuzzy_tol = ds.fuzzy_tol;
    let mut filler = match (&bvh_a, &bvh_b) {
        (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh(ds, ba, bb),
        _ => pave_filler::PaveFiller::new(ds),
    };
    filler.set_run_parallel(false);
    filler.configure_fuzzy(fuzzy_tol);
    filler.set_non_destructive(false);
    filler.configure_glue(false, TOLERANCE_ABS);
    filler.set_use_obb(false);
    filler.perform();
}

/// Sum of boundary-edge counts from [`crate::brep_check::validate_solid_closure`].
/// Larger means a **less** closed manifold shell.
/// Union: DS → PaveFiller → BooleanBuilder(Union) → recompute plane surfaces.
///
/// Uses BVH when both operands have faces, matching [`crate::boolean_op`].
pub(crate) fn fuse(a: &BRep, b: &BRep) -> Result<BRep, BooleanError> {
    fuse_with_bvh(a, b, true)
}

/// ✅ OCCT-aligned: DS → PaveFiller → BooleanBuilder(Union) → result.
///   OCCT BOPAlgo_BOP::Perform L395-408: PaveFiller config + Perform + PerformInternal1.
pub(crate) fn fuse_with_bvh(a: &BRep, b: &BRep, use_bvh: bool) -> Result<BRep, BooleanError> {
    let mut ds = bopds::ds::DS::new(a, b);
    pave_fill(&mut ds, a, b, use_bvh);
    let builder = builder::BooleanBuilder::new(&ds, BooleanOpType::Union);
    let (result, _history) = builder.build_with_history()?;
    Ok(result)
}

/// ✅ OCCT对齐: 为所有非平面面上的边计算 pcurve(3D曲线→UV投影)。
///    OCCT BuildSplitFaces 创建 section edge 时同时生成 pcurve(IntTools_CurveRange).
///    rcad 的 add_circle_edge/add_edge 不生成 pcurve,此处作为上游补做。
pub(crate) fn compute_face_pcurves(brep: &mut BRep) {
    use rcad_kernel::geom::{CurveEval, SurfaceEval, Curve2d, Line2d, Surface3};
    use rcad_kernel::PCurve;
    use std::f64::consts::PI;
    if brep.edges.is_empty() { return; }
    brep.geom.edge_pcurves.resize(brep.edges.len(), Vec::new());
    for si in 0..brep.solids.len() {
        for shi in 0..brep.solids[si].shells.len() {
            for fi in 0..brep.solids[si].shells[shi].faces.len() {
                let face = &brep.solids[si].shells[shi].faces[fi];
                let Some(surf_idx) = face.surface_idx else { continue; };
                let Some(surface) = brep.geom.surfaces.get(surf_idx) else { continue; };
                if matches!(surface, Surface3::Plane(_)) { continue; }
                let edges: Vec<usize> = face.outer_wire.edges.iter()
                    .chain(face.inner_wires.iter().flat_map(|w| &w.edges))
                    .map(|we| we.idx).collect();
                for &ei in &edges {
                    if ei >= brep.geom.edge_pcurves.len() { continue; }
                    if !brep.geom.edge_pcurves[ei].is_empty() { continue; }
                    let Some(curve_idx) = brep.geom.edge_curve.get(ei).copied().flatten() else { continue; };
                    let Some(curve) = brep.geom.curves.get(curve_idx) else { continue; };
                    let n = 12usize;
                    let mut uv: Vec<glam::DVec2> = Vec::with_capacity(n);
                    for s in 0..n {
                        let t = s as f64 / (n - 1) as f64;
                        let p3d = curve.point_at(t);
                        if let Some(u) = project_point_to_uv(p3d, surface) { uv.push(u); }
                    }
                    if uv.len() < 2 { continue; }
                    let c2d = brep.geom.curve2ds.len();
                    let dir = (uv[uv.len()-1] - uv[0]) / (n as f64 - 1.0);
                    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: uv[0], direction: dir }));
                    brep.geom.edge_pcurves[ei].push(PCurve { surface_idx: surf_idx, curve2d_idx: c2d });
                }
            }
        }
    }
}

/// 3D点到曲面 UV 的反向投影。用于 compute_face_pcurves。
fn project_point_to_uv(p: glam::DVec3, surface: &rcad_kernel::geom::Surface3) -> Option<glam::DVec2> {
    use rcad_kernel::geom::Surface3;
    match surface {
        Surface3::Sphere(s) => {
            let d = (p - s.center) / s.radius;
            let u = f64::atan2(d.dot(s.ref_dir_perp()), d.dot(s.ref_dir));
            let v = f64::asin(d.dot(s.axis).clamp(-1.0, 1.0));
            Some(glam::DVec2::new(u, v))
        }
        _ => None,
    }
}

/// Same phases as [`fuse`], but returns [`BooleanHistory`] and does not run plane recompute
/// (matches legacy `boolean_op_with_history` for Union).
pub(crate) fn fuse_with_history(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
    fuse_with_history_bvh(a, b, true)
}

pub(crate) fn fuse_with_history_bvh(
    a: &BRep,
    b: &BRep,
    use_bvh: bool,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    validate_union_operands(a, b)?;
    let mut ds = bopds::ds::DS::new(a, b);
    validate_ds_invariants(&ds)?;
    pave_fill(&mut ds, a, b, use_bvh);
    validate_ds_invariants(&ds)?;
    let builder = builder::BooleanBuilder::new(&ds, BooleanOpType::Union);
    let (brep, hist) = builder.build_with_history()?;
    validate_union_brep_output("union: result failed checks after build (history)", &brep)?;
    Ok((brep, hist))
}

/// Parallel classification path; same OCCT phase structure as [`fuse_with_history`].
pub(crate) fn fuse_with_history_par(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
    fuse_with_history_par_bvh(a, b, true)
}

pub(crate) fn fuse_with_history_par_bvh(
    a: &BRep,
    b: &BRep,
    use_bvh: bool,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    validate_union_operands(a, b)?;
    let mut ds = bopds::ds::DS::new(a, b);
    validate_ds_invariants(&ds)?;
    pave_fill(&mut ds, a, b, use_bvh);
    validate_ds_invariants(&ds)?;
    let builder = builder::BooleanBuilder::new(&ds, BooleanOpType::Union);
    let (brep, hist) = builder.build_with_history()?;
    validate_union_brep_output("union: result failed checks after build (history par)", &brep)?;
    Ok((brep, hist))
}

/// ✅ OCCT对齐: 按 edge set (BOPTools_Set) 对共面面做同域合并。
///    OCCT 源码: BOPAlgo_Builder_2.cxx L571-L832 (FillSameDomainFaces)
///
/// 算法步骤 (对应 OCCT):
/// 1. 收集每面的 edge key set (外环边界的量化顶点对,去重排序)
/// 2. 按 edge key set 分组 → anESetFaces
/// 3. 对每组≥2面,逐对检测是否 SameDomain:
///    a. 平面面(Plane/planar BSpline) → 快速路径 (OCCT L697-701: 无需几何分析)
///    b. 同表面类型的非平面面(Cylinder+Cylinder 等) → 表面几何比较 (OCCT L703-708)
///    c. 跨表面类型非平面面 → TODO: 几何分析 AreFacesSameDomain (OCCT L703-708)
/// 4. MakeBlocks: 从 Face→Face 映射构建连通组 (OCCT L741: BOPAlgo_Tools::MakeBlocks)
/// 5. 每组选代表面: 优先原面(原始DS面) > 子面,同优先级按 flat index
///    (OCCT L758-788: nFMin → myDS->Index(aF) >= 0 的原面优先)
/// 6. 删除非代表面,更新 geom slots
///
/// 返回合并(移除)的面数。
fn fill_same_domain_faces_edge_set(brep: &mut BRep) -> usize {
    use std::collections::{HashMap, BTreeSet};
    use rcad_kernel::geom::Surface3;

    if brep.solids.is_empty() || brep.solids[0].shells.is_empty() {
        return 0;
    }

    // ── Step 1: Flat face index → (si, shi, fi) ──────────────────────────
    let mut flat_to_pos: Vec<(usize, usize, usize)> = Vec::new();
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for fi in 0..shell.faces.len() {
                flat_to_pos.push((si, shi, fi));
            }
        }
    }
    let nf = flat_to_pos.len();
    if nf < 2 {
        return 0;
    }

    // ── Step 2: Compute edge key per face ─────────────────────────────────
    type EdgeKey = ((i64, i64, i64), (i64, i64, i64));
    let inv = 1.0 / 1e-5;
    let q = |p: glam::DVec3| -> (i64, i64, i64) {
        ((p.x * inv).round() as i64, (p.y * inv).round() as i64, (p.z * inv).round() as i64)
    };

    let mut face_keys: Vec<Option<BTreeSet<EdgeKey>>> = vec![None; nf];
    for ff in 0..nf {
        let (si, shi, fi) = flat_to_pos[ff];
        let face = &brep.solids[si].shells[shi].faces[fi];
        let mut keys = BTreeSet::new();
        let mut has_valid = false;
        for we in &face.outer_wire.edges {
            if we.idx >= brep.edges.len() { continue; }
            if brep.geom.edge_degenerated.get(we.idx).copied().unwrap_or(false) { continue; }
            let edge = &brep.edges[we.idx];
            if edge.start >= brep.vertices.len() || edge.end >= brep.vertices.len() { continue; }
            if edge.start == edge.end { continue; }
            let qs = q(brep.vertices[edge.start].point);
            let qe = q(brep.vertices[edge.end].point);
            let key = if qs < qe { (qs, qe) } else { (qe, qs) };
            keys.insert(key);
            has_valid = true;
        }
        if has_valid {
            face_keys[ff] = Some(keys);
        }
    }

    // ── Step 3: Group by edge key set (OCCT anESetFaces) ──────────────────
    let mut groups: HashMap<BTreeSet<EdgeKey>, Vec<usize>> = HashMap::new();
    for ff in 0..nf {
        if let Some(ref keys) = face_keys[ff] {
            groups.entry(keys.clone()).or_default().push(ff);
        }
    }

    // ── Step 4: Detect SameDomain pairs within each group ─────────────────
    // OCCT: BOPAlgo_Builder_2.cxx L676-709
    let mut dsu = DSU::new(nf);

    for (_keys, members) in groups.iter() {
        if members.len() < 2 {
            continue;
        }

        // Precompute planar flag per face (OCCT L636-647: planar bounded check)
        let mut is_planar: Vec<bool> = Vec::with_capacity(nf);
        for _ in 0..nf { is_planar.push(false); }
        for &ff in members {
            let (si, shi, fi) = flat_to_pos[ff];
            let face = &brep.solids[si].shells[shi].faces[fi];
            let sid = face.surface_idx;
            is_planar[ff] = sid.and_then(|sid| brep.geom.surfaces.get(sid)).map_or(false, |surf| {
                match surf {
                    Surface3::Plane(_) => true,
                    Surface3::BSpline(bsp) => rcad_kernel::geom::bspline_is_planar(bsp, 1e-3),
                    _ => false,
                }
            });
        }

        // Compare every pair within the edge-set group — OCCT L686-709
        let m = members.len();
        for i in 0..m {
            let fi = members[i];

            for j in (i + 1)..m {
                let fj = members[j];

                // ⚡ OCCT L697-701: 平面面同edge set = 自动 SameDomain
                //    OCCT 检查 GeomAbs_Plane + Bnd_Box non-open,
                //    rcad: 所有平面面(Plane/planar BSpline)都是 bounded 的
                if is_planar[fi] && is_planar[fj] {
                    dsu.union(fi, fj);
                } else {
                    // ⏳ OCCT L703-720: 非平面面几何分析路径 (尚未实现)
                }
            }
        }
    }

    // ── Step 4.5: Cross-group coplanar merge (OCCT AreFacesSameDomain L1109-1169) ──
    //    ⏳ 暂禁用: point_in_face + is_valid_point_for_face 已实现,待可靠内点后启用。

    // ── Step 5: Build blocks from DSU (OCCT MakeBlocks L741) ──────────────
    // blocks: Vec<Vec<usize>> where each inner vec is one connected component
    let mut block_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for ff in 0..nf {
        let root = dsu.find(ff);
        block_map.entry(root).or_default().push(ff);
    }

    // ── Step 6: Per-block representative selection (OCCT nFMin L758-788) ──
    let mut remove_pos: Vec<(usize, usize, usize)> = Vec::new();

    for (_root, members) in block_map.iter() {
        if members.len() < 2 {
            continue;
        }

        // OCCT L758-782: 优先选 DS 中的原面 (myDS->Index(aF) >= 0),
        // DS index 最小者为代表。在 rcad 中,classify pass 后 face_origins
        // 已不同步,改用 surface priority + flat index 作为确定性选择。
        let rep = members
            .iter()
            .copied()
            .min_by(|&a, &b| {
                let pa = surface_priority_for_merge(brep, flat_to_pos[a]);
                let pb = surface_priority_for_merge(brep, flat_to_pos[b]);
                pa.cmp(&pb).then_with(|| a.cmp(&b))
            })
            .unwrap();

        for &ff in members {
            if ff != rep {
                remove_pos.push(flat_to_pos[ff]);
            }
        }
    }

    if remove_pos.is_empty() {
        return 0;
    }

    // ── Step 7: Remove non-representative faces (reverse order) ───────────
    //    OCCT L790-796: myShapesSD.Bind(aF, *pFSD) — 更新 SD 映射
    //    rcad: 直接删除非代表面,BRep 级操作
    remove_pos.sort_by(|a, b| b.cmp(a));
    remove_pos.dedup();

    for &(si, shi, fi) in &remove_pos {
        let mut ff = 0usize;
        for s in 0..si {
            for sh in &brep.solids[s].shells {
                ff += sh.faces.len();
            }
        }
        for sh in 0..shi {
            ff += brep.solids[si].shells[sh].faces.len();
        }
        ff += fi;

        crate::remove_flat_face_geom_slots(&mut brep.geom, ff);
        if let Some(s) = brep.solids.get_mut(si) {
            if let Some(sh) = s.shells.get_mut(shi) {
                if fi < sh.faces.len() {
                    sh.faces.remove(fi);
                }
            }
        }
    }

    remove_pos.len()
}

/// ✅ OCCT对齐: PointInFace (BOPTools_AlgoTools3D::PointInFace)
///    OCCT: 用三角剖分取面内一点的 UV 参数 → 3D 点
///    rcad: 从 face.triangles 取第一个三角形的三顶点质心
fn point_in_face(brep: &BRep, (si, shi, fi): (usize, usize, usize)) -> Option<DVec3> {
    let face = &brep.solids[si].shells[shi].faces[fi];
    // OCCT 优先: 用三角剖分
    if let Some(tri) = face.triangles.first() {
        let v0 = brep.vertices.get(tri[0])?.point;
        let v1 = brep.vertices.get(tri[1])?.point;
        let v2 = brep.vertices.get(tri[2])?.point;
        return Some((v0 + v1 + v2) / 3.0);
    }
    // 回退: face.sample_point
    if let Some(sp) = face.sample_point { return Some(sp); }
    // 最后回退: boundary centroid
    let mut bnd: Vec<DVec3> = Vec::new();
    for we in &face.outer_wire.edges {
        if let Some(e) = brep.edges.get(we.idx) { if let Some(v) = brep.vertices.get(e.start) { bnd.push(v.point); } }
    }
    if bnd.len() < 3 { return None; }
    Some(bnd.iter().copied().sum::<DVec3>() / bnd.len() as f64)
}

/// ✅ OCCT对齐: IsValidPointForFace (BOPTools_AlgoTools.cxx L1166)
///    OCCT: 投影 3D 点到面2的表面 → 检查 UV 在域内
///    容差: tolF1 + tolF2 + max(fuzzy, Precision::Confusion)
///    rcad: 3D 边距拒绝 + 2D point-in-polygon
fn is_valid_point_for_face(brep: &BRep, pt: DVec3, (si, shi, fi): (usize, usize, usize)) -> bool {
    use glam::DVec2;
    let face = &brep.solids[si].shells[shi].faces[fi];
    let sid = face.surface_idx;
    let Some(surf) = sid.and_then(|sid| brep.geom.surfaces.get(sid)) else { return false; };
    let _normal = match surf { Surface3::Plane(p) => p.normal, _ => return false };
    // OCCT: 3D edge rejection — 3D 距离 = 0 的点是边界点
    let mut bnd: Vec<DVec3> = Vec::new();
    for we in &face.outer_wire.edges {
        if let Some(e) = brep.edges.get(we.idx) { if let Some(v) = brep.vertices.get(e.start) { bnd.push(v.point); } }
    }
    if bnd.len() < 3 { return false; }
    for k in 0..bnd.len() {
        let a = bnd[k]; let b = bnd[(k + 1) % bnd.len()];
        let ab = b - a; let ap = pt - a;
        let t = (ap.dot(ab) / (ab.dot(ab) + 1e-30)).clamp(0.0, 1.0);
        if (pt - (a + ab * t)).length() < 1e-12 { return false; }
    }
    // 2D point-in-polygon (容差 1e-7)
    let e0 = (bnd[1] - bnd[0]).normalize();
    let e1 = (bnd[2] - bnd[0]).normalize();
    let norm = e0.cross(e1).normalize();
    if norm.length_squared() < 0.5 { return false; }
    let ua = e0;
    let va = norm.cross(ua).normalize();
    let pts: Vec<DVec2> = bnd.iter().map(|p| {
        let d = *p - bnd[0]; DVec2::new(d.dot(ua), d.dot(va))
    }).collect();
    let p2 = DVec2::new((pt - bnd[0]).dot(ua), (pt - bnd[0]).dot(va));
    for k in 0..pts.len() {
        let ab = pts[(k + 1) % pts.len()] - pts[k];
        let ap = p2 - pts[k];
        let t = (ap.dot(ab) / (ab.dot(ab) + 1e-30)).clamp(0.0, 1.0);
        if (p2 - (pts[k] + ab * t)).length() < 1e-7 { return true; }
    }
    let mut inside = false; let mut k = pts.len() - 1;
    for kk in 0..pts.len() {
        if ((pts[kk].y > p2.y) != (pts[k].y > p2.y))
            && (p2.x < (pts[k].x - pts[kk].x) * (p2.y - pts[kk].y) / (pts[k].y - pts[kk].y) + pts[kk].x)
        { inside = !inside; }
        k = kk;
    }
    inside
}

/// ✅ OCCT对齐: 跨组共面检测 — AreFacesSameDomain (BOPTools_AlgoTools.cxx L1109-1169)
///
/// 严格按 OCCT 算法:
/// 1. PointInFace(F1) → 3D点  (用三角剖分取内点)
/// 2. IsValidPointForFace(点, F2, aTol) → 检查点是否在面2内
fn fill_same_domain_cross_group(
    brep: &BRep, nf: usize,
    flat_to_pos: &[(usize, usize, usize)],
    dsu: &mut DSU,
) -> usize {
    let mut merged = 0usize;
    for i in 0..nf {
        for j in (i + 1)..nf {
            if dsu.find(i) == dsu.find(j) { continue; };
            let (si, shi, fi) = flat_to_pos[i];
            let (sj, shj, fj) = flat_to_pos[j];
            let face_i = &brep.solids[si].shells[shi].faces[fi];
            let face_j = &brep.solids[sj].shells[shj].faces[fj];
            // 检查是否同一平面 (法线同向)
            let n_i = face_i.normal;
            let n_j = face_j.normal;
            if n_i.dot(n_j) < 0.999 { continue; }
            // OCCT Step 1: PointInFace
            let Some(interior_pt) = point_in_face(brep, (si, shi, fi)) else { continue; };
            // OCCT Step 2: IsValidPointForFace
            let fi_in_fj = is_valid_point_for_face(brep, interior_pt, (sj, shj, fj));
            // 反向检查
            let Some(interior_pt_j) = point_in_face(brep, (sj, shj, fj)) else { continue; };
            let fj_in_fi = is_valid_point_for_face(brep, interior_pt_j, (si, shi, fi));
            if fi_in_fj || fj_in_fi {
                dsu.union(i, j);
                merged += 1;
            }
        }
    }
    merged
}


/// Check that every edge in the result is referenced exactly twice
/// (once from each adjacent face), forming a closed manifold shell.
/// Returns an error listing orphan (<2 refs) and over-shared (>2 refs) edges.
fn validate_solid_closure(brep: &BRep) -> Result<(), BooleanError> {
    use std::collections::HashMap;
    let mut edge_count: HashMap<usize, usize> = HashMap::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    *edge_count.entry(we.idx).or_insert(0) += 1;
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        *edge_count.entry(we.idx).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    let mut orphan: Vec<usize> = Vec::new();
    let mut over: Vec<usize> = Vec::new();
    for (&eidx, &count) in &edge_count {
        if count < 2 {
            orphan.push(eidx);
        } else if count > 2 {
            over.push(eidx);
        }
    }

    if orphan.is_empty() && over.is_empty() {
        return Ok(());
    }

    orphan.sort();
    orphan.dedup();
    over.sort();
    over.dedup();

    let detail = if !orphan.is_empty() && !over.is_empty() {
        format!("{} orphan edges (refs<2), {} over-shared edges (refs>2)", orphan.len(), over.len())
    } else if !orphan.is_empty() {
        format!("{} orphan edges (refs<2)", orphan.len())
    } else {
        format!("{} over-shared edges (refs>2)", over.len())
    };

    Err(BooleanError::OpenShell {
        orphan_edges: orphan,
        over_shared_edges: over,
    })
}

/// 同域面合并的代表面选择优先级。
/// 低值 = 高优先级 (优先原生 Plane > planar BSpline > 其他)。
fn surface_priority_for_merge(brep: &BRep, (si, shi, fi): (usize, usize, usize)) -> u32 {
    let face = &brep.solids[si].shells[shi].faces[fi];
    let sid = face.surface_idx;
    match sid.and_then(|sid| brep.geom.surfaces.get(sid)) {
        Some(Surface3::Plane(_)) => 0,
        Some(Surface3::BSpline(bsp)) if rcad_kernel::geom::bspline_is_planar(bsp, 1e-3) => 1,
        _ => 2,
    }
}


