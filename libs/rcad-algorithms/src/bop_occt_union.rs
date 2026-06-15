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
        for &ci in &f.face_info.curves_in {
            if ci >= nic {
                return Err(BooleanError::InvalidResult(
                    "union: DS face_info.curves_in out of range",
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
                    let n = face.normal;
                    if !n.x.is_finite()
                        || !n.y.is_finite()
                        || !n.z.is_finite()
                        || n.length_squared() <= 0.0
                    {
                        return Err(BooleanError::InvalidResult(msg));
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
fn pave_fill(ds: &mut bopds::ds::DS, a: &BRep, b: &BRep, use_bvh: bool) {
    if use_bvh {
        let (bvh_a, bvh_b) = optional_bvhs(a, b);
        let mut filler = match (&bvh_a, &bvh_b) {
            (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh(ds, ba, bb),
            _ => pave_filler::PaveFiller::new(ds),
        };
        filler.perform();
    } else {
        let mut filler = pave_filler::PaveFiller::new(ds);
        filler.perform();
    }
}

/// Sum of boundary-edge counts from [`crate::brep_check::validate_solid_closure`].
/// Larger means a **less** closed manifold shell.
fn solid_closure_boundary_penalty(brep: &BRep) -> usize {
    let r = crate::brep_check::validate_solid_closure(brep);
    r.issues
        .iter()
        .filter_map(|i| match i {
            crate::brep_check::CheckIssue::SolidNotClosed {
                boundary_edge_count,
                ..
            } => Some(*boundary_edge_count),
            _ => None,
        })
        .sum()
}

fn shell_face_total(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

/// Union: DS → PaveFiller → BooleanBuilder(Union) → recompute plane surfaces.
///
/// Uses BVH when both operands have faces, matching [`crate::boolean_op`].
pub(crate) fn fuse(a: &BRep, b: &BRep) -> Result<BRep, BooleanError> {
    fuse_with_bvh(a, b, true)
}

pub(crate) fn fuse_with_bvh(a: &BRep, b: &BRep, use_bvh: bool) -> Result<BRep, BooleanError> {
    validate_union_operands(a, b)?;
    let mut ds = bopds::ds::DS::new(a, b);
    validate_ds_invariants(&ds)?;
    pave_fill(&mut ds, a, b, use_bvh);
    validate_ds_invariants(&ds)?;
    let builder = builder::BooleanBuilder::new(&ds, BooleanOpType::Union);
    let (mut result, mut history) = builder.build_with_history()?;
    if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
        let mut fi = 0usize;
        for solid in &result.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let area = rcad_kernel::face_surface_area(&result, face, fi);
                    eprintln!("[FACE_SA] fi={} sa={:.6}", fi, area);
                    fi += 1;
                }
            }
        }
    }
    validate_union_brep_output("union: result failed checks after build", &result)?;
    // ✅ OCCT对齐: 为所有非平面面上的边创建 pcurve,供 merge 函数的 seam edge 检测。
    //    OCCT BuildSplitFaces 创建 section edge 时同时生成 pcurve。
    //    rcad 的 add_circle_edge/add_edge 不生成,此处在上游补做。
    compute_face_pcurves(&mut result);
    // Sew boolean seam vertices so coplanar edge-adjacent patches share endpoints; orthogonal
    // fuse / `unify_same_domain_faces` need coincident topology to merge remaining fragments.
    let (sewn, _) = merge_close_vertices(&result, crate::tolerance::TOLERANCE_ABS * 1000.0);
    result = sewn;
    // Position-based edge dedup: merge edges at the same geometric endpoints.
    // After merge_close_vertices, most vertices are merged but some edge pairs
    // still connect different vertex indices at the same positions (PaveFiller
    // noise > 6.4e-6 tolerance).  Use 1e-5 quantization to catch these.
    // ✅ OCCT-aligned: dedup includes curve identity — edges with the same
    //    endpoint positions but different curves are kept separate
    //    (TopoDS_Edge uniqueness requires both vertex pair AND curve).
    {
        use std::collections::HashMap;
        let inv = 1.0 / 1e-5;
        let q = |p: glam::DVec3| -> (i64, i64, i64) {
            ((p.x * inv).round() as i64, (p.y * inv).round() as i64, (p.z * inv).round() as i64)
        };
        let vpos: Vec<glam::DVec3> = result.vertices.iter().map(|v| v.point).collect();
        let mut canon: Vec<usize> = (0..result.edges.len()).collect();
        // ✅ OCCT-aligned: key = (quantized_start, quantized_end, edge_curve_idx)
        //    to prevent merging edges with different curves (sphere seam vs box edge).
        let edge_curves = &result.geom.edge_curve;
        let mut geom: HashMap<((i64,i64,i64),(i64,i64,i64),Option<usize>), usize> = HashMap::new();
        for ei in 0..result.edges.len() {
            if let Some(e) = result.edges.get(ei) {
                let a = q(vpos[e.start]); let b = q(vpos[e.end]);
                let key = if a < b { (a, b, edge_curves.get(ei).copied().flatten()) } else { (b, a, edge_curves.get(ei).copied().flatten()) };
                let entry = geom.entry(key).or_insert(ei);
                if *entry != ei { canon[ei] = *entry; }
            }
        }
        for s in &mut result.solids {
            for sh in &mut s.shells {
                for face in &mut sh.faces {
                    for we in &mut face.outer_wire.edges { we.idx = canon[we.idx]; }
                    for w in &mut face.inner_wires { for we in &mut w.edges { we.idx = canon[we.idx]; } }
                }
            }
        }
    }
    // Index-based edge dedup (catches any remaining duplicates).
    result = crate::deduplicate_edges(result);
    geom_populate::recompute_plane_surfaces(&mut result);
    validate_union_brep_output("union: result failed checks after vertex merge", &result)?;
    // Deduplicate edges so adjacent sub-faces share edge topology.
    result = crate::deduplicate_edges(result);
    // OCCT ClassifyFaces (BuildSolid): remove interior faces classified as
    // State_IN for the OTHER operand.  Runs BEFORE the second butterfly merge
    // so history.face_origins (from build_with_history's internal merge) are valid.
    {
        use rcad_kernel::geom::Surface3;
        use crate::classify::{classify_point, Classification};
        use crate::history::FaceOrigin;
        use crate::tolerance::TOLERANCE_ABS;
        let a_ids: Vec<usize> = (0..ds.a_face_count).collect();
        let b_ids: Vec<usize> = (ds.a_face_count..ds.faces.len()).collect();
        let mut to_remove: Vec<(usize, usize, usize)> = Vec::new();
        let origins = &history.face_origins;
        for (si, solid) in result.solids.iter().enumerate() {
            for (shi, shell) in solid.shells.iter().enumerate() {
                for (fi, face) in shell.faces.iter().enumerate() {
                    let origin = origins.get(fi).copied();
                    let classify_ids = match origin {
                        Some(FaceOrigin::FromA(_)) => &b_ids,
                        Some(FaceOrigin::FromB(_)) => &a_ids,
                        _ => {
                            let is_plane = result.geom.face_surface.get(fi).copied().flatten()
                                .and_then(|si| result.geom.surfaces.get(si))
                                .is_some_and(|s| matches!(s, Surface3::Plane(_)));
                            if is_plane { &b_ids } else { &a_ids }
                        }
                    };
                    // ✅ OCCT对齐: BOPTools_AlgoTools::ComputeState for Face.
                    // 遍历面的每条边，找到一条不在对面实体边界上的边，用其中点做分类。
                    // OCCT 源码: BOPTools_AlgoTools.cxx L650-L674
                    let mut classify_pt = face.sample_point.unwrap_or(DVec3::ZERO);
                    let mut found_pt = face.sample_point.is_some();
                    if !found_pt {
                        // 优先从边中点做分类（OCCT 主路径）
                        for we in &face.outer_wire.edges {
                            if let Some(edge) = result.edges.get(we.idx) {
                                let p1 = result.vertices.get(edge.start).map(|v| v.point);
                                let p2 = result.vertices.get(edge.end).map(|v| v.point);
                                if let (Some(p1), Some(p2)) = (p1, p2) {
                                    let mid = (p1 + p2) * 0.5;
                                    let cls = classify_point(mid, classify_ids, &ds);
                                    if cls != Classification::On {
                                        // 边的中点不在对方实体边界上→可用的分类点
                                        classify_pt = mid;
                                        found_pt = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if !found_pt {
                        // 所有边都在边界上→回退到顶点质心（OCCT PointInFace fallback）
                        let mut _pts: Vec<DVec3> = Vec::new();
                        for we in &face.outer_wire.edges {
                            if let Some(edge) = result.edges.get(we.idx) {
                                _pts.push(result.vertices[edge.start].point);
                                _pts.push(result.vertices[edge.end].point);
                            }
                        }
                        if _pts.len() < 3 { continue; }
                        classify_pt = _pts.iter().copied().sum::<DVec3>() / _pts.len() as f64;
                    }
                    let cls = classify_point(classify_pt, classify_ids, &ds);
                    if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
                        let surf_name = result.geom.face_surface.get(fi).copied().flatten()
                            .and_then(|si| result.geom.surfaces.get(si))
                            .map(|s| match s {
                                Surface3::Plane(_) => "Plane",
                                Surface3::BSpline(_) => "BSpline",
                                _ => "Other",
                            }).unwrap_or("None");
                        eprintln!("[CLASSIFY_CLS] fi={fi} origin={origin:?} surf={surf_name} pt=({:.6},{:.6},{:.6}) cls={:?}",
                            classify_pt.x, classify_pt.y, classify_pt.z, cls);
                    }
                    let interior = cls == Classification::In;
                    if interior {
                        if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
                            let surf_name = result.geom.face_surface.get(fi).copied().flatten()
                                .and_then(|si| result.geom.surfaces.get(si))
                                .map(|s| match s {
                                    Surface3::Plane(_) => "Plane",
                                    Surface3::BSpline(_) => "BSpline",
                                    _ => "Other",
                                }).unwrap_or("None");
                            eprintln!("[CLASSIFY] REMOVE fi={fi} origin={origin:?} surf={surf_name} pt=({:.6},{:.6},{:.6}) n=({:.4},{:.4},{:.4})",
                                classify_pt.x, classify_pt.y, classify_pt.z,
                                face.normal.x, face.normal.y, face.normal.z);
                        }
                        to_remove.push((si, shi, fi));
                    } else if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
                        let surf_name = result.geom.face_surface.get(fi).copied().flatten()
                            .and_then(|si| result.geom.surfaces.get(si))
                            .map(|s| match s {
                                Surface3::Plane(_) => "Plane",
                                Surface3::BSpline(_) => "BSpline",
                                _ => "Other",
                            }).unwrap_or("None");
                        eprintln!("[CLASSIFY] KEEP   fi={fi} origin={origin:?} surf={surf_name} pt=({:.6},{:.6},{:.6}) n=({:.4},{:.4},{:.4})",
                            classify_pt.x, classify_pt.y, classify_pt.z,
                            face.normal.x, face.normal.y, face.normal.z);
                    }
                }
            }
        }
        if !to_remove.is_empty() {
            eprintln!("[CLASSIFY_FACES] removing {} interior faces (from {})", to_remove.len(),
                result.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count());
            to_remove.sort_by(|a, b| b.cmp(a));
            for &(si, shi, fi) in &to_remove {
                // Compute flat face index: sum faces in preceding solids/shells + fi
                let mut ff = 0usize;
                for s in 0..si {
                    for sh in &result.solids[s].shells {
                        ff += sh.faces.len();
                    }
                }
                for sh in 0..shi {
                    ff += result.solids[si].shells[sh].faces.len();
                }
                ff += fi;
                crate::remove_flat_face_geom_slots(&mut result.geom, ff);
                if let Some(s) = result.solids.get_mut(si) {
                    if let Some(sh) = s.shells.get_mut(shi) {
                        if fi < sh.faces.len() { sh.faces.remove(fi); }
                    }
                }
            }
        }
    }

    // Second OCCT FillSameDomainFaces pass (butterfly merge).
    // ✅ OCCT对齐: 用 edge set (BOPTools_Set) 对共面面做跨类型(Plane+BSpline)分组合并。
    if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
        let nf = result.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        eprintln!("[CLASSIFY] after classify: {} faces", nf);
    }
    let merged_count = fill_same_domain_faces_edge_set(&mut result);
    if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
        let nf = result.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        eprintln!("[CLASSIFY] edge-set merge: removed {} faces, remaining {}", merged_count, nf);
    }
    result = crate::prune_unused_topology(result);
    result = remove_interior_faces(result);
    result = crate::prune_unused_topology(result);
    result = crate::deduplicate_edges(result);
    result.geom.edge_degenerated.resize(result.edges.len(), false);
    for (i, e) in result.edges.iter().enumerate() {
        if e.start == e.end { result.geom.edge_degenerated[i] = true; }
    }
    // OCCT FillSameDomainFaces: edge-set grouping + surface comparison
    // (edge-index based, handles all surface types).
    {
        let (merged, _cnt) = crate::occt_fill_same_domain_faces(&result);
        result = merged;
    }
    // Legacy surface-index based merge (BuildSolid loop/area equivalent).
    let (merged, _cnt) = crate::occt_merge_same_surface_faces(&result);
    result = merged;
    // compact_brep after merge removes edges/vertices that were only referenced by
    // the now-removed faces — OCCT FillSameDomainFaces does not clear the edge list,
    // but compact_brep is the RCAD equivalent of rebuilding the shape after merge.
    result = crate::prune_unused_topology(result);
    result = crate::deduplicate_edges(result);
    // OCCT aligned: mark start==end edges as degenerated
    result.geom.edge_degenerated.resize(result.edges.len(), false);
    for (i, e) in result.edges.iter().enumerate() {
        if e.start == e.end {
            result.geom.edge_degenerated[i] = true;
        }
    }
    // Remove degenerate edges (start==end) that survive compact_brep due to
    // face-replacement index juggling in the merge function.
    {
        use rcad_kernel::topology::WireEdge;
        let nz: Vec<_> = result.edges.iter().map(|e| e.start != e.end).collect();
        if nz.iter().any(|&k| !k) {
            let mut remap: Vec<Option<usize>> = vec![None; result.edges.len()];
            let mut nedges: Vec<_> = Vec::new();
            for (i, e) in result.edges.iter().enumerate() {
                let is_degen = result.geom.edge_degenerated.get(i).copied().unwrap_or(false);
                if e.start != e.end || is_degen {
                    remap[i] = Some(nedges.len());
                    nedges.push(*e);
                }
            }
            for s in &mut result.solids {
                for sh in &mut s.shells {
                    for f in &mut sh.faces {
                        f.outer_wire.edges.retain(|we| remap.get(we.idx).copied().flatten().is_some());
                        for we in &mut f.outer_wire.edges { we.idx = remap[we.idx].unwrap(); }
                        for w in &mut f.inner_wires {
                            w.edges.retain(|we| remap.get(we.idx).copied().flatten().is_some());
                            for we in &mut w.edges { we.idx = remap[we.idx].unwrap(); }
                        }
                    }
                }
            }
            result.edges = nedges;
            result = crate::prune_unused_topology(result);
        }
    }
    validate_union_brep_output("union: result failed checks after same-plane merge", &result)?;
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
    let (brep, hist) = builder.build_with_history_par()?;
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

/// ✅ OCCT对齐: BuildSolid — 排除被覆盖的内部面 (OCCT BOPAlgo_BuilderSolid)
///    OCCT ShellSplitter 从分裂面重建 shell 时,不构成外环的面被排除。
///    rcad: 对同平面同法线的面对,如果较小面的 centroid 在较大面内 → 移除较小面。
///    返回修改后的 BRep。
fn remove_interior_faces(mut brep: BRep) -> BRep {
    for si in 0..brep.solids.len() {
        for shi in 0..brep.solids[si].shells.len() {
            let nf = brep.solids[si].shells[shi].faces.len();
            if nf < 2 { continue; }

            // Collect planar face info
            struct Pf { area: f64, cx: DVec3, bd: Vec<DVec3>, nm: DVec3 }
            let mut planes: Vec<(usize, Pf)> = Vec::new();
            for fi in 0..nf {
                let face = &brep.solids[si].shells[shi].faces[fi];
                let sid = face.surface_idx;
                let Some(surf) = sid.and_then(|sid| brep.geom.surfaces.get(sid)) else { continue; };
                let is_p = match surf { Surface3::Plane(_) => true, Surface3::BSpline(b) => rcad_kernel::geom::bspline_is_planar(b,1e-3), _ => false };
                if !is_p { continue; }
                let nm = match surf { Surface3::Plane(p) => p.normal, _ => face.normal };
                let mut bd: Vec<DVec3> = Vec::new();
                for we in &face.outer_wire.edges { if let Some(e)=brep.edges.get(we.idx) { if let Some(v)=brep.vertices.get(e.start) { bd.push(v.point); } } }
                if bd.len() < 3 { continue; }
                let cx = bd.iter().copied().sum::<DVec3>() / bd.len() as f64;
                // 2D polygon area
                let (ua, va) = {
                    let e0 = (bd[1]-bd[0]).normalize(); let e1 = (bd[2]-bd[0]).normalize();
                    let n = e0.cross(e1).normalize(); (e0, n.cross(e0).normalize())
                };
                let area: f64 = {
                    let pts: Vec<DVec2> = bd.iter().map(|p| { let d=*p-bd[0]; DVec2::new(d.dot(ua), d.dot(va)) }).collect();
                    let m=pts.len(); let mut a=0.0f64; let mut k=m-1;
                    for i in 0..m { a += pts[i].x*pts[k].y - pts[k].x*pts[i].y; k=i; }
                    (a/2.0f64).abs()
                };
                let _ = area;
                planes.push((fi, Pf { area, cx, bd, nm }));
            }

            // Detect interior faces: smaller face inside larger face on same plane
            let mut to_remove: Vec<usize> = Vec::new();
            for i in 0..planes.len() {
                for j in (i+1)..planes.len() {
                    let (fi_a, ref a) = planes[i];
                    let (fi_b, ref b) = planes[j];
                    if a.nm.dot(b.nm) < 0.999 { continue; }
                    // 3D-boundary-check + pip
                    let inside = |cx: DVec3, bd: &[DVec3]| -> bool {
                        for k in 0..bd.len() {
                            let a = bd[k]; let b = bd[(k+1)%bd.len()];
                            let ab = b-a; let ap = cx-a;
                            let t = (ap.dot(ab)/(ab.dot(ab)+1e-30)).clamp(0.0,1.0);
                            if (cx-(a+ab*t)).length() < 1e-7 { return false; }
                        }
                        let e0 = (bd[1]-bd[0]).normalize(); let e1 = (bd[2]-bd[0]).normalize();
                        let n = e0.cross(e1).normalize(); if n.length_squared()<0.5 { return false; }
                        let ua=e0; let va=n.cross(ua).normalize();
                        let pts: Vec<DVec2> = bd.iter().map(|p| { let d=*p-bd[0]; DVec2::new(d.dot(ua),d.dot(va)) }).collect();
                        let p2 = DVec2::new((cx-bd[0]).dot(ua), (cx-bd[0]).dot(va));
                        let mut inside=false; let mut k=pts.len()-1;
                        for kk in 0..pts.len() {
                            if ((pts[kk].y>p2.y)!=(pts[k].y>p2.y))&&(p2.x<(pts[k].x-pts[kk].x)*(p2.y-pts[kk].y)/(pts[k].y-pts[kk].y)+pts[kk].x) { inside=!inside; }
                            k=kk;
                        }
                        inside
                    };
                    let a_in_b = inside(a.cx, &b.bd);
                    let b_in_a = inside(b.cx, &a.bd);
                    if a_in_b && !b_in_a && a.area < b.area { to_remove.push(fi_a); }
                    if b_in_a && !a_in_b && b.area < a.area { to_remove.push(fi_b); }
                }
            }

            // Remove interior faces
            to_remove.sort_unstable_by(|a,b| b.cmp(a));
            to_remove.dedup();
            for &fi in &to_remove {
                if fi >= brep.solids[si].shells[shi].faces.len() { continue; }
                // compute flat index
                let mut ff = 0usize;
                for s in 0..si { for sh in &brep.solids[s].shells { ff += sh.faces.len(); } }
                for sh in 0..shi { ff += brep.solids[si].shells[sh].faces.len(); }
                ff += fi;
                crate::remove_flat_face_geom_slots(&mut brep.geom, ff);
                brep.solids[si].shells[shi].faces.remove(fi);
            }
        }
    }
    brep
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

