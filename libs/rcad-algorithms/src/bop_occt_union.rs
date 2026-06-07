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

use glam::DVec3;
use crate::brep_repair::merge_close_vertices;
use crate::bopds;
use crate::bopds::ds::{DS, Interference};
use crate::builder;
use crate::bvh;
use crate::geom_populate;
use crate::history::BooleanHistory;
use crate::pave_filler;
use crate::tolerance::*;
use crate::total_surface_area;
use crate::BooleanError;
use crate::BooleanOpType;
use rcad_kernel::BRep;

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
    let mut result = builder.build()?;
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
    // Sew boolean seam vertices so coplanar edge-adjacent patches share endpoints; orthogonal
    // fuse / `unify_same_domain_faces` need coincident topology to merge remaining fragments.
    let (sewn, _) = merge_close_vertices(&result, crate::tolerance::TOLERANCE_ABS * 64.0);
    result = sewn;
    geom_populate::recompute_plane_surfaces(&mut result);
    validate_union_brep_output("union: result failed checks after vertex merge", &result)?;
    // Deduplicate edges so adjacent sub-faces share edge topology (OCCT:
    // edges from BOPAlgo_BuilderFace are shared between split faces).
    // Without this, unify_same_domain_faces cannot detect shared edges.
    result = crate::deduplicate_edges(result);
    let checkpoint = result.clone();
    merge_coplanar_orthogonal_unify(&mut result);
    geom_populate::recompute_plane_surfaces(&mut result);

    // OCCT ClassifyFaces (BuildSolid): for each face, classify its center
    // against the OTHER operand using the DS.  Plane faces (from A-side)
    // are tested against B; BSpline faces (from B-side) against A.
    // If State_IN → the face is interior to the union → remove.
    {
        use rcad_kernel::geom::Surface3;
        use crate::classify::{classify_point, Classification};
        let a_ids: Vec<usize> = (0..ds.a_face_count).collect();
        let b_ids: Vec<usize> = (ds.a_face_count..ds.faces.len()).collect();
        let mut to_remove: Vec<(usize, usize, usize)> = Vec::new();
        for (si, solid) in result.solids.iter().enumerate() {
            for (shi, shell) in solid.shells.iter().enumerate() {
                for (fi, face) in shell.faces.iter().enumerate() {
                    let si2 = result.geom.face_surface.get(fi).copied().flatten();
                    let is_plane = si2.and_then(|i| result.geom.surfaces.get(i))
                        .is_some_and(|s| matches!(s, Surface3::Plane(_)));
                    // Compute centroid from boundary vertices
                    let mut pts: Vec<DVec3> = Vec::new();
                    for we in &face.outer_wire.edges {
                        if let Some(edge) = result.edges.get(we.idx) {
                            pts.push(result.vertices[edge.start].point);
                            pts.push(result.vertices[edge.end].point);
                        }
                    }
                    if pts.len() < 3 { continue; }
                    let centroid = pts.iter().copied().sum::<DVec3>() / pts.len() as f64;
                    // Classify against the other operand
                    let interior = if is_plane {
                        classify_point(centroid, &b_ids, &ds) == Classification::In
                    } else {
                        classify_point(centroid, &a_ids, &ds) == Classification::In
                    };
                    if interior {
                        to_remove.push((si, shi, fi));
                    }
                }
            }
        }
        if !to_remove.is_empty() {
            eprintln!("[CLASSIFY_FACES] removing {} interior faces", to_remove.len());
            to_remove.sort_by(|a, b| b.cmp(a));
            for &(si, shi, fi) in &to_remove {
                if let Some(s) = result.solids.get_mut(si) {
                    if let Some(sh) = s.shells.get_mut(shi) {
                        if fi < sh.faces.len() { sh.faces.remove(fi); }
                    }
                }
            }
        }
    }

    let pen_before = solid_closure_boundary_penalty(&checkpoint);
    let pen_after = solid_closure_boundary_penalty(&result);
    let n_before = shell_face_total(&checkpoint);
    let n_after = shell_face_total(&result);
    let sa_before = total_surface_area(&checkpoint);
    let sa_after = total_surface_area(&result);
    // Revert if the merge changes SA significantly (in either direction) with a large
    // face-count drop — indicates that unify_same_domain_faces merged cylinder sub-faces
    // whose UV projection overcounts (or undercounts) area in try_cylinder_trimmed_face_area.
    let large_sa_change = {
        let abs_delta = (sa_after - sa_before).abs();
        abs_delta > 1e-6 && abs_delta > 0.005 * sa_before.max(1.0)
    };
    let suspicious_planar_snarl = (pen_after > pen_before
        && (n_after <= 10 || n_after.saturating_add(12) < n_before))
        || (large_sa_change && n_after < n_before && n_after.saturating_add(12) < n_before)
        || (large_sa_change && (sa_after - sa_before).abs() > 0.10 * sa_before.max(1.0));
    if suspicious_planar_snarl {
        result = checkpoint;
        geom_populate::recompute_plane_surfaces(&mut result);
    }

    validate_union_brep_output("union: result failed checks after same-plane merge", &result)?;
    Ok(result)
}

/// Merge coplanar orthogonal panels left by boolean split (re-run groups after each success).
fn merge_coplanar_orthogonal_unify(brep: &mut BRep) {
    let tol = TOLERANCE_ABS;
    for _ in 0..12 {
        let (next, m) = crate::orthogonal_face_fuse::fuse_orthogonal_coplanar_faces(brep, tol);
        *brep = next;
        geom_populate::recompute_plane_surfaces(brep);
        // OCCT FillSameDomainFaces: butterfly merge only (no adjacency merge)
        crate::unify_same_domain_faces(brep);
        geom_populate::recompute_plane_surfaces(brep);
        if m == 0 { break; }
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
