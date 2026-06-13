//! Fuse coplanar **axis-aligned** rectangular patches into one [`Face`] using a 2D
//! axis-aligned rectangle union on a grid, producing one outer boundary and optional
//! inner wires (holes). Complements [`unify_same_domain_faces`](crate::unify_same_domain_faces),
//! which only merges along shared **edges** and leaves corner-only adjacency split.

use glam::DVec3;
use std::collections::{HashMap, HashSet};
use rcad_kernel::geom::{Plane, Surface3};
use rcad_kernel::topology::{Edge, Face, Vertex, Wire, WireEdge};
use rcad_kernel::{face_surface_area, surface_area, volume, BRep};

use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::{
    TOLERANCE_ABS, TOLERANCE_ADAPTIVE_MAX, TOLERANCE_AREA_REL, TOLERANCE_COORD_SUB,
    TOLERANCE_FLOAT_DEDUP, TOLERANCE_FLOAT_ULTRA, TOLERANCE_LEN_MIN, TOLERANCE_MESH_LEGACY,
    TOLERANCE_RETRY_LADDER_COARSE, TOLERANCE_RETRY_LADDER_MID, TOLERANCE_TOL_SCALE_MICRO,
    TOLERANCE_VEC_SQ_MIN,
};

type Pt = (i64, i64);

/// Minimum [`face_surface_area`] / axis-aligned UV bbox area before treating strict bbox containment
/// as 鈥渋nner duplicate rectangle vs outer trimmed cap鈥?(`bcommon_simple/C8`). Trimmed intersection
/// patches (`bcommon_simple/G5`) often occupy much less than their axis UV bbox when the boundary is
/// skewed.
const MIN_REDUNDANT_AXIS_UV_FILL: f64 = 0.97;

/// When one axis-aligned bbox is strictly inside another on the same plane, require the inner bbox
/// **area** to be at most this fraction of the outer (`area(inner)/area(outer)`). Otherwise the pair
/// is not treated as an untrimmed duplicate inside a trimmed patch (`bcommon_simple/G5`).
const MAX_STRICT_INNER_BBOX_AREA_FRAC: f64 = 0.55;

fn qpt(x: f64, y: f64, scale: f64) -> Pt {
    ((x * scale).round() as i64, (y * scale).round() as i64)
}

/// Remove redundant axis-aligned faces on the same infinite plane when one face's **world-UV
/// axis-aligned bounding box** is strictly contained in another's (OCCT `bcommon_simple/C8`:
/// an untrimmed `1脳1` top patch is kept alongside the true trimmed cap; they do not pass
/// [`rects_2d_bbox_positive_area_overlap`] for orthogonal fuse, but the smaller bbox lies inside
/// the larger).
///
/// Only considers faces with **no inner wires** and normals snapped to 卤X/卤Y/卤Z. Returns the
/// cleaned BRep and how many faces were removed.
pub fn remove_axis_coplanar_redundant_child_faces(brep: &BRep, tol: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let mut total_removed = 0usize;
    let t = tol.max(TOLERANCE_ABS);

    for si in 0..out.solids.len() {
        for shi in 0..out.solids[si].shells.len() {
            loop {
                let n = out.solids[si].shells[shi].faces.len();
                if n < 2 {
                    break;
                }
                let mut remove_flat: Option<usize> = None;
                'outer: for fi in 0..n {
                    for fj in (fi + 1)..n {
                        if let Some(rm_flat) =
                            try_pick_redundant_axis_coplanar_face(&out, si, shi, fi, fj, t)
                        {
                            remove_flat = Some(rm_flat);
                            break 'outer;
                        }
                    }
                }
                let Some(flat_rm) = remove_flat else {
                    break;
                };
                let Some((rsi, rshi, local_rm)) = flat_index_to_local_shell_face(&out, flat_rm) else {
                    break;
                };
                debug_assert_eq!((rsi, rshi), (si, shi));
                out.solids[rsi].shells[rshi].faces.remove(local_rm);
                crate::remove_flat_face_geom_slots(&mut out.geom, flat_rm);
                total_removed += 1;
            }
        }
    }

    (out, total_removed)
}

fn flat_index_to_local_shell_face(brep: &BRep, flat: usize) -> Option<(usize, usize, usize)> {
    let mut k = 0usize;
    for si in 0..brep.solids.len() {
        for shi in 0..brep.solids[si].shells.len() {
            for fi in 0..brep.solids[si].shells[shi].faces.len() {
                if k == flat {
                    return Some((si, shi, fi));
                }
                k += 1;
            }
        }
    }
    None
}

fn try_pick_redundant_axis_coplanar_face(
    brep: &BRep,
    si: usize,
    shi: usize,
    fi: usize,
    fj: usize,
    tol: f64,
) -> Option<usize> {
    let fi_flat = flat_face_index(brep, si, shi, fi);
    let fj_flat = flat_face_index(brep, si, shi, fj);
    let shell = &brep.solids[si].shells[shi];
    let face_i = &shell.faces[fi];
    let face_j = &shell.faces[fj];
    if !face_i.inner_wires.is_empty() || !face_j.inner_wires.is_empty() {
        return None;
    }
    // Both faces must be on planar surfaces 鈥?a cylindrical wall whose normal happens
    // to snap to an axis (e.g. (-0,-1,0) 鈫?-Y) would be falsely matched with a planar
    // face on the same axis, causing incorrect removal of valid geometry (H8 box鈭ヽyl).
    if !face_is_plane(brep, fi_flat) || !face_is_plane(brep, fj_flat) {
        return None;
    }
    let n_i = snap_almost_axis(face_i.normal.normalize_or_zero());
    let n_j = snap_almost_axis(face_j.normal.normalize_or_zero());
    if axis_aligned_world_plane_uv_axes(n_i).is_none()
        || axis_aligned_world_plane_uv_axes(n_j).is_none()
    {
        return None;
    }
    let p_i = face_first_point(brep, face_i)?;
    let p_j = face_first_point(brep, face_j)?;
    let d_i = n_i.dot(p_i);
    let d_j = n_j.dot(p_j);
    let (n_i_c, d_i_c) = canonicalize_plane_n_d(n_i, d_i);
    let (n_j_c, d_j_c) = canonicalize_plane_n_d(n_j, d_j);
    let key_i = plane_key(n_i_c, d_i_c, tol);
    let key_j = plane_key(n_j_c, d_j_c, tol);
    if key_i != key_j {
        return None;
    }

    let bi = face_axis_world_bbox(brep, face_i, n_i)?;
    let bj = face_axis_world_bbox(brep, face_j, n_j)?;
    let scale = (bi.1 - bi.0)
        .abs()
        .max(bi.3 - bi.2)
        .abs()
        .max(bj.1 - bj.0)
        .abs()
        .max(bj.3 - bj.2)
        .abs()
        .max(1.0);
    let eps = (TOLERANCE_COORD_SUB * scale).max(tol * TOLERANCE_TOL_SCALE_MICRO);

    let subset = |a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)| -> bool {
        let (au0, au1, av0, av1) = a;
        let (bu0, bu1, bv0, bv1) = b;
        au0 >= bu0 - eps && au1 <= bu1 + eps && av0 >= bv0 - eps && av1 <= bv1 + eps
    };

    let bbox_area = |b: (f64, f64, f64, f64)| -> f64 {
        let (u0, u1, v0, v1) = b;
        (u1 - u0).abs().max(0.0) * (v1 - v0).abs().max(0.0)
    };

    let ni = face_i.outer_wire.edges.len();
    let nj = face_j.outer_wire.edges.len();

    // Strictly smaller bbox inside larger: drop the inner (untrimmed cap vs trimmed cap).
    if subset(bi, bj) && !subset(bj, bi) {
        let aj = bbox_area(bj).max(eps * eps);
        let ar = bbox_area(bi) / aj;
        if ar > MAX_STRICT_INNER_BBOX_AREA_FRAC {
            return None;
        }
        let r = redundant_axis_uv_bbox_fill_ratio(brep, face_i, n_i, fi_flat, scale).unwrap_or(0.0);
        if r < MIN_REDUNDANT_AXIS_UV_FILL {
            return None;
        }
        return Some(fi_flat);
    }
    if subset(bj, bi) && !subset(bi, bj) {
        let ai = bbox_area(bi).max(eps * eps);
        let ar = bbox_area(bj) / ai;
        if ar > MAX_STRICT_INNER_BBOX_AREA_FRAC {
            return None;
        }
        let r = redundant_axis_uv_bbox_fill_ratio(brep, face_j, n_j, fj_flat, scale).unwrap_or(0.0);
        if r < MIN_REDUNDANT_AXIS_UV_FILL {
            return None;
        }
        return Some(fj_flat);
    }
    // Equal bboxes (duplicate patches): keep richer topology.
    if subset(bi, bj) && subset(bj, bi) {
        let ri = redundant_axis_uv_bbox_fill_ratio(brep, face_i, n_i, fi_flat, scale).unwrap_or(0.0);
        let rj = redundant_axis_uv_bbox_fill_ratio(brep, face_j, n_j, fj_flat, scale).unwrap_or(0.0);
        if ri < MIN_REDUNDANT_AXIS_UV_FILL || rj < MIN_REDUNDANT_AXIS_UV_FILL {
            return None;
        }
        if ni < nj {
            return Some(fi_flat);
        }
        if nj < ni {
            return Some(fj_flat);
        }
        if bbox_area(bi) + eps < bbox_area(bj) {
            return Some(fi_flat);
        }
        if bbox_area(bj) + eps < bbox_area(bi) {
            return Some(fj_flat);
        }
        // Identical patches: drop the higher shell-local index.
        return Some(fj_flat);
    }
    None
}

/// If removing a single face leaves [`volume`] unchanged (within `vol_abs_tol`) but lowers
/// [`surface_area`] by a noticeable amount, treat that face as a duplicate boundary sheet and drop
/// it (OCCT `bcommon_simple/C8`: one axis patch is redundant while the diagonal patch carries the
/// same material boundary).
///
/// Only considers faces whose normals snap to 卤X/卤Y/卤Z and whose [`face_surface_area`] fills their
/// axis-aligned UV bbox enough ([`MIN_REDUNDANT_AXIS_UV_FILL`]); other faces are skipped so trimmed
/// intersections (`bcommon_simple/G5`) are not peeled incorrectly.
///
/// Scans faces in **flat index order** and returns on the **first** match so behaviour stays
/// deterministic. Intended only for post-processing **plane鈥損lane intersections** (callers gate).
pub fn remove_spurious_intersection_face_preserving_volume(
    brep: &BRep,
    vol_abs_tol: f64,
) -> (BRep, usize) {
    let v0 = volume(brep);
    let a0 = surface_area(brep);
    let n = brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum::<usize>();
    let vtol = vol_abs_tol.max(TOLERANCE_FLOAT_DEDUP * v0.abs().max(1.0));
    for flat_rm in 0..n {
        let Some((si, shi, local_rm)) = flat_index_to_local_shell_face(brep, flat_rm) else {
            continue;
        };
        let face = &brep.solids[si].shells[shi].faces[local_rm];
        if !face_is_plane(brep, flat_rm) {
            continue;
        }
        let n_axis = snap_almost_axis(face.normal.normalize_or_zero());
        if axis_aligned_world_plane_uv_axes(n_axis).is_none() {
            continue;
        }
        let Some(bb) = face_axis_world_bbox(brep, face, n_axis) else {
            continue;
        };
        let scale = (bb.1 - bb.0).abs().max((bb.3 - bb.2).abs()).max(1.0);
        let r =
            redundant_axis_uv_bbox_fill_ratio(brep, face, n_axis, flat_rm, scale).unwrap_or(0.0);
        if r < MIN_REDUNDANT_AXIS_UV_FILL {
            continue;
        }

        let fp = face_first_point(brep, face).unwrap_or(DVec3::ZERO);
        let d = n_axis.dot(fp);
        let (n_i_c, d_i_c) = canonicalize_plane_n_d(n_axis, d);
        let pk = plane_key(n_i_c, d_i_c, vol_abs_tol.max(TOLERANCE_COORD_SUB));
        // Count faces on the same infinite plane in this shell.  If this is the only one,
        // removing it would leave a hole 鈥?it is structurally necessary even if the
        // divergence-theorem volume contribution happens to be zero (e.g. a face at x=0
        // with normal (-1,0,0) where r路n = 0 for all points).
        // Only planar faces are counted 鈥?curved faces (e.g. cylinder wall) whose normals
        // happen to snap to an axis are excluded to prevent false plane matches (H8 box鈭ヽyl).
        let same_plane_count = (0..brep.solids[si].shells[shi].faces.len())
            .filter(|&local_fi| {
                let of_flat = flat_face_index(brep, si, shi, local_fi);
                if !face_is_plane(brep, of_flat) {
                    return false;
                }
                let of = &brep.solids[si].shells[shi].faces[local_fi];
                let on = snap_almost_axis(of.normal.normalize_or_zero());
                let _axes = match axis_aligned_world_plane_uv_axes(on) {
                    Some(a) => a,
                    None => return false,
                };
                let op = face_first_point(brep, of).unwrap_or(DVec3::ZERO);
                let od = on.dot(op);
                let (on_c, od_c) = canonicalize_plane_n_d(on, od);
                plane_key(on_c, od_c, vol_abs_tol.max(TOLERANCE_COORD_SUB)) == pk
            })
            .count();
        if same_plane_count <= 1 {
            continue;
        }

        let mut br = brep.clone();
        let Some((si, shi, local_rm)) = flat_index_to_local_shell_face(&br, flat_rm) else {
            continue;
        };
        br.solids[si].shells[shi].faces.remove(local_rm);
        crate::remove_flat_face_geom_slots(&mut br.geom, flat_rm);
        let v1 = volume(&br);
        let a1 = surface_area(&br);
        if (v1 - v0).abs() <= vtol && a1 < a0 - 0.25 {
            return (br, 1);
        }
    }
    (brep.clone(), 0)
}

/// Merge groups of coplanar orthogonal faces in each shell into single faces (with holes).
pub fn fuse_orthogonal_coplanar_faces(brep: &BRep, tol: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let mut total = 0usize;
    let t = tol.max(TOLERANCE_ABS);

    for si in 0..out.solids.len() {
        for shi in 0..out.solids[si].shells.len() {
            fuse_orthogonal_in_shell(&mut out, si, shi, t, &mut total);
        }
    }

    (out, total)
}

fn fuse_orthogonal_in_shell(brep: &mut BRep, si: usize, shi: usize, tol: f64, total: &mut usize) {
    // After one merge, face indices and shell length change. Rebuild coplanar groups each pass
    // instead of iterating stale `fis` from the previous snapshot.
    loop {
        let n = brep.solids[si].shells[shi].faces.len();
        if n < 2 {
            return;
        }

        let mut groups: std::collections::HashMap<(i64, i64, i64, i64), Vec<usize>> =
            std::collections::HashMap::new();
        for fi in 0..n {
            let face = &brep.solids[si].shells[shi].faces[fi];
            if !face.inner_wires.is_empty() {
                continue;
            }
            // Skip curved surfaces (cylinder, cone, sphere, torus, BSpline, etc.)
            // — the orthogonal coplanar fuse projects faces onto a plane, which
            // flattens curved patches (e.g. tessellated cylinder from a full-wrap
            // TwoEllipses intersection). Curved surfaces are merged by the
            // downstream `unify_same_domain_faces` pass which preserves type.
            if !face_is_plane(brep, flat_face_index(brep, si, shi, fi)) {
                continue;
            }
            let Some(p0) = face_first_point(brep, face) else {
                continue;
            };
            let nrm = face.normal.normalize_or_zero();
            if nrm.length_squared() < TOLERANCE_VEC_SQ_MIN {
                continue;
            }
            let nrm = snap_almost_axis(nrm);
            let d = nrm.dot(p0);
            let (n_key, d_key) = canonicalize_plane_n_d(nrm, d);
            let key = plane_key(n_key, d_key, tol);
            groups.entry(key).or_default().push(fi);
        }

        // Same infinite plane key can include disjoint UV islands. Only merge 2D bbox
        // components that overlap with *positive* area. A second 鈥渆dge only鈥?pass was tried and
        // breaks `boolean_op_healed` on partial 0.5-overlap unions; `unify_same_domain_faces` can
        // still coalesce some edge-coincident fragments afterward.
        let mut group_list: Vec<Vec<usize>> = Vec::new();
        for mut fis in groups.into_values() {
            if fis.len() < 2 {
                continue;
            }
            fis.sort_unstable();
            for sub in split_fis_by_plane_uv_connectivity(brep, si, shi, &fis, tol) {
                if sub.len() >= 2 {
                    group_list.push(sub);
                }
            }
        }
        group_list.sort_by_key(|fis| (fis[0], fis.len()));

        let mut merged = false;
        for fis in group_list {
            if try_fuse_orthogonal_group(brep, si, shi, &fis, tol, false) {
                *total += 1;
                merged = true;
                break;
            }
        }
        if !merged {
            if try_fuse_one_axis_aligned_edge_adjacent_pair(brep, si, shi, tol) {
                *total += 1;
                continue;
            }
            return;
        }
    }
}

fn rects_2d_bbox_positive_area_overlap(
    a: (f64, f64, f64, f64),
    b: (f64, f64, f64, f64),
    gap: f64,
) -> bool {
    let (au0, au1, av0, av1) = a;
    let (bu0, bu1, bv0, bv1) = b;
    let wu = au1.min(bu1) - au0.max(bu0);
    let wv = av1.min(bv1) - av0.max(bv0);
    wu > gap && wv > gap
}

/// Two axis-aligned UV rectangles share a full edge: one overlap dimension is ~0, the other > `tt`.
/// Corner-only (`wu`鈮? and `wv`鈮?) and separated rectangles are excluded.
fn rects_2d_bbox_share_full_edge(a: (f64, f64, f64, f64), b: (f64, f64, f64, f64), tt: f64) -> bool {
    let wu = a.1.min(b.1) - a.0.max(b.0);
    let wv = a.3.min(b.3) - a.2.max(b.2);
    let touch = tt.max(TOLERANCE_COORD_SUB);
    if wu < -touch || wv < -touch {
        return false;
    }
    let near = |x: f64| x.abs() <= touch;
    if near(wu) && near(wv) {
        return false;
    }
    (wu > touch && near(wv)) || (wv > touch && near(wu))
}

fn face_outer_vertex_set(brep: &BRep, face: &Face) -> std::collections::HashSet<usize> {
    let mut s = std::collections::HashSet::new();
    for we in &face.outer_wire.edges {
        let Some(e) = brep.edges.get(we.idx) else {
            continue;
        };
        if we.forward {
            s.insert(e.start);
            s.insert(e.end);
        } else {
            s.insert(e.end);
            s.insert(e.start);
        }
    }
    s
}

fn face_outer_edge_segments(brep: &BRep, face: &Face) -> Vec<(DVec3, DVec3)> {
    let mut out = Vec::new();
    for we in &face.outer_wire.edges {
        let Some(e) = brep.edges.get(we.idx) else {
            continue;
        };
        let Some(va) = brep.vertices.get(if we.forward { e.start } else { e.end }) else {
            continue;
        };
        let Some(vb) = brep.vertices.get(if we.forward { e.end } else { e.start }) else {
            continue;
        };
        out.push((va.point, vb.point));
    }
    out
}

fn segment_coincident(a0: DVec3, a1: DVec3, b0: DVec3, b1: DVec3, tol: f64) -> bool {
    let same = |u: DVec3, v: DVec3| (u - v).length() <= tol;
    (same(a0, b0) && same(a1, b1)) || (same(a0, b1) && same(a1, b0))
}

/// True when the two faces share a full edge (鈮? coincident vertex indices, or the same segment
/// within `geom_tol` 鈥?booleans often duplicate vertex indices on seams).
fn faces_share_full_edge_geom(
    brep: &BRep,
    face_i: &Face,
    face_j: &Face,
    geom_tol: f64,
) -> bool {
    let vi = face_outer_vertex_set(brep, face_i);
    let vj = face_outer_vertex_set(brep, face_j);
    if vi.intersection(&vj).count() >= 2 {
        return true;
    }
    let segsi = face_outer_edge_segments(brep, face_i);
    let segsj = face_outer_edge_segments(brep, face_j);
    for &(p0, p1) in &segsi {
        for &(q0, q1) in &segsj {
            if segment_coincident(p0, p1, q0, q1, geom_tol) {
                return true;
            }
        }
    }
    false
}

fn faces_share_full_edge_geom_by_index(
    brep: &BRep,
    si: usize,
    shi: usize,
    fi: usize,
    fj: usize,
    geom_tol: f64,
) -> bool {
    let fi_face = &brep.solids[si].shells[shi].faces[fi];
    let fj_face = &brep.solids[si].shells[shi].faces[fj];
    faces_share_full_edge_geom(brep, fi_face, fj_face, geom_tol)
}

/// Merge one coplanar axis-aligned pair that shares a geometric full edge in UV but has no 2D area
/// overlap (not in the same `split_fis_by_plane_uv_connectivity` component). Skips corner-only
/// contacts via [`rects_2d_bbox_share_full_edge`] and requires [`faces_share_full_edge_geom`].
fn try_fuse_one_axis_aligned_edge_adjacent_pair(brep: &mut BRep, si: usize, shi: usize, tol: f64) -> bool {
    let n = brep.solids[si].shells[shi].faces.len();
    if n < 2 {
        return false;
    }
    let t = tol.max(TOLERANCE_ABS);
    let gap = (t * 1e2).max(TOLERANCE_MESH_LEGACY);
    let touch = (t * 1e3).max(TOLERANCE_RETRY_LADDER_MID);
    let geom_edge_tol = (t * 1e4).max(TOLERANCE_RETRY_LADDER_COARSE);

    let mut meta: Vec<Option<((i64, i64, i64, i64), (f64, f64, f64, f64))>> = vec![None; n];
    for fi in 0..n {
        let face = &brep.solids[si].shells[shi].faces[fi];
        if !face.inner_wires.is_empty() {
            continue;
        }
        let Some(p0) = face_first_point(brep, face) else {
            continue;
        };
        let nrm = snap_almost_axis(face.normal.normalize_or_zero());
        if axis_aligned_world_plane_uv_axes(nrm).is_none() {
            continue;
        }
        let d = nrm.dot(p0);
        let (nk, dk) = canonicalize_plane_n_d(nrm, d);
        let key = plane_key(nk, dk, t);
        let Some(bb) = face_axis_world_bbox(brep, face, nrm) else {
            continue;
        };
        meta[fi] = Some((key, bb));
    }

    for fi in 0..n {
        let Some((ki, bi)) = meta[fi] else {
            continue;
        };
        for fj in (fi + 1)..n {
            let Some((kj, bj)) = meta[fj] else {
                continue;
            };
            if ki != kj {
                continue;
            }
            // Skip non-Plane faces (BSpline, etc.) — only merge axis-aligned
            // Plane faces.  The main fuse_orthogonal_in_shell already filters
            // by face_is_plane; this edge-adjacent fallback must do the same
            // (bfuse_simple B4/B5: planar BSpline faces were being merged with
            // adjacent Plane faces, losing surface type).
            if !face_is_plane(brep, flat_face_index(brep, si, shi, fi))
                || !face_is_plane(brep, flat_face_index(brep, si, shi, fj))
            {
                continue;
            }
            if rects_2d_bbox_positive_area_overlap(bi, bj, gap) {
                continue;
            }
            if !rects_2d_bbox_share_full_edge(bi, bj, touch) {
                continue;
            }
            if !faces_share_full_edge_geom_by_index(brep, si, shi, fi, fj, geom_edge_tol) {
                continue;
            }
            // XOR / difference lumps can meet along a planar interface with geometrically coincident
            // edges but **opposing** outward normals on the two sheets; do not orthogonal-fuse those.
            let n1 = brep.solids[si].shells[shi].faces[fi]
                .normal
                .normalize_or_zero();
            let n2 = brep.solids[si].shells[shi].faces[fj]
                .normal
                .normalize_or_zero();
            const MIN_SAME_SHELL_COS: f64 = 1e-6;
            if n1.dot(n2) <= MIN_SAME_SHELL_COS {
                continue;
            }
            let cur_n = brep.solids[si].shells[shi].faces.len();
            if fi >= cur_n || fj >= cur_n {
                continue;
            }
            if try_fuse_orthogonal_group(brep, si, shi, &[fi, fj], tol, true) {
                return true;
            }
        }
    }
    false
}

fn uf_find(p: &mut [usize], mut i: usize) -> usize {
    while p[i] != i {
        p[i] = p[p[i]];
        i = p[i];
    }
    i
}

fn uf_unite(p: &mut [usize], i: usize, j: usize) {
    let a = uf_find(p, i);
    let b = uf_find(p, j);
    if a != b {
        p[a] = b;
    }
}

fn face_axis_world_bbox(brep: &BRep, face: &Face, n: DVec3) -> Option<(f64, f64, f64, f64)> {
    let [i, j] = axis_aligned_world_plane_uv_axes(n)?;
    let poly = face_outer_points(brep, face);
    if poly.is_empty() {
        return None;
    }
    let uv: Vec<(f64, f64)> = poly.iter().map(|p| (p[i], p[j])).collect();
    Some(bbox2d(&uv))
}

fn axis_uv_bbox_rect_area(b: (f64, f64, f64, f64)) -> f64 {
    let (u0, u1, v0, v1) = b;
    (u1 - u0).abs() * (v1 - v0).abs()
}

/// Physical face area divided by axis-aligned UV bbox area for 卤axis planes (same projection as
/// [`face_axis_world_bbox`]). Near 1 鈬?patch fills its bbox rectangle (typical redundant caps).
fn redundant_axis_uv_bbox_fill_ratio(
    brep: &BRep,
    face: &Face,
    n: DVec3,
    flat_idx: usize,
    scale: f64,
) -> Option<f64> {
    let b = face_axis_world_bbox(brep, face, n)?;
    let ab = axis_uv_bbox_rect_area(b);
    let s = scale.max(1.0);
    let eps_area = (TOLERANCE_FLOAT_DEDUP * s * s).max(TOLERANCE_FLOAT_ULTRA);
    if !ab.is_finite() || ab <= eps_area {
        return None;
    }
    let a = face_surface_area(brep, face, flat_idx);
    Some(a / ab)
}

/// Split coplanar face indices into groups that are each connected in 2D world UV.
/// Non鈥揳xis-aligned planes keep a single bucket (previous behavior).
fn split_fis_by_plane_uv_connectivity(
    brep: &BRep,
    si: usize,
    shi: usize,
    fis: &[usize],
    tol: f64,
) -> Vec<Vec<usize>> {
    if fis.len() < 2 {
        return Vec::new();
    }
    let shell = &brep.solids[si].shells[shi];
    let f0 = &shell.faces[fis[0]];
    let n = snap_almost_axis(f0.normal.normalize_or_zero());
    if axis_aligned_world_plane_uv_axes(n).is_none() {
        return vec![fis.to_vec()];
    }
    let bboxes: Vec<Option<_>> = fis
        .iter()
        .map(|&fi| face_axis_world_bbox(brep, &shell.faces[fi], n))
        .collect();
    if bboxes.iter().any(|b| b.is_none()) {
        return vec![fis.to_vec()];
    }
    let bboxes: Vec<_> = bboxes.into_iter().map(|b| b.unwrap()).collect();
    let gap = (tol * 1e2).max(TOLERANCE_MESH_LEGACY);
    let mut parent: Vec<usize> = (0..fis.len()).collect();
    for i in 0..fis.len() {
        for j in (i + 1)..fis.len() {
            if rects_2d_bbox_positive_area_overlap(bboxes[i], bboxes[j], gap) {
                uf_unite(&mut parent, i, j);
            }
        }
    }
    let mut buckets: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..fis.len() {
        let r = uf_find(&mut parent, i);
        buckets.entry(r).or_default().push(fis[i]);
    }
    let mut out: Vec<Vec<usize>> = buckets
        .into_values()
        .filter(|g| g.len() >= 2)
        .collect();
    for g in &mut out {
        g.sort_unstable();
    }
    out
}

fn plane_key(n: DVec3, d: f64, tol: f64) -> (i64, i64, i64, i64) {
    let s = 1.0 / tol.max(TOLERANCE_COORD_SUB);
    (
        (n.x * s).round() as i64,
        (n.y * s).round() as i64,
        (n.z * s).round() as i64,
        (d * s).round() as i64,
    )
}

/// If `n` is within `tol_dir` of an axis, snap to exact 卤X/卤Y/卤Z.
/// When `n` is 卤X/卤Y/卤Z, return the two world axis indices that span the plane (e.g. Z虃 鈫?(x,y)).
fn axis_aligned_world_plane_uv_axes(n: DVec3) -> Option<[usize; 2]> {
    let a = n.abs();
    if a.x > 1.0 - 2.0 * TOLERANCE_ADAPTIVE_MAX {
        Some([1, 2])
    } else if a.y > 1.0 - 2.0 * TOLERANCE_ADAPTIVE_MAX {
        Some([0, 2])
    } else if a.z > 1.0 - 2.0 * TOLERANCE_ADAPTIVE_MAX {
        Some([0, 1])
    } else {
        None
    }
}

/// Inverse of [`axis_aligned_world_plane_uv_axes`]: build 3D from two free world coordinates.
fn point_from_axis_plane_world_uv(n: DVec3, o: DVec3, u: f64, v2: f64) -> DVec3 {
    let a = n.abs();
    if a.x > 1.0 - 2.0 * TOLERANCE_ADAPTIVE_MAX {
        DVec3::new(o.x, u, v2)
    } else if a.y > 1.0 - 2.0 * TOLERANCE_ADAPTIVE_MAX {
        DVec3::new(u, o.y, v2)
    } else {
        DVec3::new(u, v2, o.z)
    }
}

fn snap_almost_axis(n: DVec3) -> DVec3 {
    let t = 2.0 * TOLERANCE_ADAPTIVE_MAX;
    for i in 0..3 {
        if n[i].abs() > 1.0 - t {
            let mut o = DVec3::ZERO;
            o[i] = n[i].signum();
            return o;
        }
    }
    n
}

/// Map `n路x = d` to a canonical `n` so that `(n, d)` and `(-n, -d)` (same
/// infinite plane) share the same key when bucketing.
fn canonicalize_plane_n_d(n: DVec3, d: f64) -> (DVec3, f64) {
    const E: f64 = TOLERANCE_LEN_MIN;
    let mut n = n;
    let mut d = d;
    if n.x < -E
        || (n.x.abs() <= E && n.y < -E)
        || (n.x.abs() <= E && n.y.abs() <= E && n.z < -E)
    {
        n = -n;
        d = -d;
    }
    (n, d)
}

fn face_first_point(brep: &BRep, face: &Face) -> Option<DVec3> {
    let we = face.outer_wire.edges.first()?;
    let e = brep.edges.get(we.idx)?;
    let vi = if we.forward { e.start } else { e.end };
    Some(brep.vertices.get(vi)?.point)
}

/// `allow_two_patch_full_edge_touch_only`: must be true only for pairs verified by
/// `faces_share_full_edge_geom_by_index` in [`try_fuse_one_axis_aligned_edge_adjacent_pair`].
/// Batch groups from UV-overlap connectivity keep the overlap-only two-patch gate so XOR / glued
/// shells do not wrongly coalesce disjoint surface sheets.
fn try_fuse_orthogonal_group(
    brep: &mut BRep,
    si: usize,
    shi: usize,
    fis: &[usize],
    tol: f64,
    allow_two_patch_full_edge_touch_only: bool,
) -> bool {
    let n_faces_shell = brep.solids[si].shells[shi].faces.len();
    if fis.is_empty() || fis[0] >= n_faces_shell {
        return false;
    }
    let plane = {
        let f0 = &brep.solids[si].shells[shi].faces[fis[0]];
        let n = f0.normal.normalize_or_zero();
        let Some(p0) = face_first_point(brep, f0) else {
            return false;
        };
        Plane {
            origin: p0,
            normal: n,
        }
    };
    let n_unit = snap_almost_axis(plane.normal.normalize_or_zero());
    let (u_axis, v_axis) = plane_local_basis(&plane);
    let world_axes = axis_aligned_world_plane_uv_axes(n_unit);
    let project_uv = |p: DVec3| -> (f64, f64) {
        if let Some([i, j]) = world_axes {
            (p[i], p[j])
        } else {
            let d = p - plane.origin;
            (d.dot(u_axis), d.dot(v_axis))
        }
    };

    let mut rects: Vec<(f64, f64, f64, f64)> = Vec::with_capacity(fis.len());
    for &fi in fis {
        if fi >= n_faces_shell {
            return false;
        }
        let face = &brep.solids[si].shells[shi].faces[fi];
        if !face.inner_wires.is_empty() {
            return false;
        }
        let poly = face_outer_points(brep, face);
        if poly.len() < 3 || poly.len() > 64 {
            return false;
        }
        let uv: Vec<(f64, f64)> = poly.iter().map(|&p| project_uv(p)).collect();
        let t = tol.max(TOLERANCE_ABS);
        let orth = polygon_is_orthogonal(&poly, &project_uv, tol);
        let ok_rect = orth
            || (poly.len() == 4 && four_uv_points_on_rectangle_corners(&uv, t))
            || (world_axes.is_some() && poly.len() == 4 && orth);
        if !ok_rect {
            return false;
        }
        let (umin, umax, vmin, vmax) = bbox2d(&uv);
        rects.push((umin, umax, vmin, vmax));
    }

    let Some(geo_rings) = union_rects_to_rings_grid(&rects, tol) else {
        return false;
    };
    if geo_rings.is_empty() {
        return false;
    }
    if geo_rings.len() > 1 {
        return false;
    }
    // Grid union rings have collinear samples along edges. After removal: 4 corners for a merged
    // rectangle, 6+ for an L (e.g. two overlapping boxes with different y/z extents).
    let t = tol.max(TOLERANCE_ABS);
    let simplified_outer = simplify_ring_collinear_uv_closed(&geo_rings[0], t);
    if simplified_outer.len() < 4 {
        return false;
    }
    if !ring_is_axis_aligned_orthogonal_uv(&simplified_outer, t) {
        return false;
    }
    let rings_for_mesh = vec![simplified_outer];
    if fis.len() == 2 {
        let g = (t * 1e2).max(TOLERANCE_MESH_LEGACY);
        let overlap = rects_2d_bbox_positive_area_overlap(rects[0], rects[1], g);
        if !overlap {
            if !allow_two_patch_full_edge_touch_only {
                return false;
            }
            let touch = (t * 1e3).max(TOLERANCE_RETRY_LADDER_MID);
            if !rects_2d_bbox_share_full_edge(rects[0], rects[1], touch) {
                return false;
            }
        }
    }

    let normal = brep.solids[si].shells[shi].faces[fis[0]].normal;

    let ring_vertices = if world_axes.is_some() {
        add_vertices_for_rings_with_eval(brep, &rings_for_mesh, |u, v| {
            point_from_axis_plane_world_uv(n_unit, plane.origin, u, v)
        }, tol)
    } else {
        add_vertices_for_rings(brep, &rings_for_mesh, &plane, u_axis, v_axis, tol)
    };
    if ring_vertices.is_empty() || ring_vertices[0].len() < 3 {
        return false;
    }

    let mut edge_pairs: Vec<(usize, usize)> = Vec::new();
    for ring_v in &ring_vertices {
        let n = ring_v.len();
        if n < 3 {
            return false;
        }
        for i in 0..n {
            edge_pairs.push((ring_v[i], ring_v[(i + 1) % n]));
        }
    }

    let base_ei = brep.edges.len();
    push_new_edges(brep, edge_pairs);

    let mut ei = base_ei;
    let n0 = ring_vertices[0].len();
    let outer_wire = Wire {
        edges: (0..n0).map(|i| {
            
            WireEdge::fwd(ei + i)
        })
        .collect(),
    };
    ei += n0;

    let mut inner_wires: Vec<Wire> = Vec::new();
    for ring_v in ring_vertices.iter().skip(1) {
        let n = ring_v.len();
        if n < 3 {
            return false;
        }
        inner_wires.push(Wire {
            edges: (0..n)
                .map(|i| {
                    
                    WireEdge::fwd(ei + i)
                })
                .collect(),
        });
        ei += n;
    }

    let merged_face = Face {
        outer_wire,
        inner_wires,
        normal,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };

    let surf_idx = {
        let p0 = face_first_point(brep, &merged_face).unwrap_or(plane.origin);
        let idx = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: p0,
            normal,
        }));
        idx
    };

    let flat_indices: Vec<usize> = fis
        .iter()
        .map(|&fi| flat_face_index(brep, si, shi, fi))
        .collect();

    replace_shell_faces_and_geom(brep, si, shi, fis, merged_face, surf_idx, &flat_indices);
    true
}

/// Union of axis-aligned rectangles 鈫?outer ring first, then hole rings (UV coords).
fn union_rects_to_rings_grid(rects: &[(f64, f64, f64, f64)], tol: f64) -> Option<Vec<Vec<(f64, f64)>>> {
    let t = tol.max(TOLERANCE_ABS);
    let (occ_ext, xs, ys) = build_padded_occ_grid(rects, t)?;
    let nx = occ_ext.len();
    let ny = occ_ext.first()?.len();
    if nx < 3 || ny < 3 {
        return None;
    }

    let mut outside = vec![vec![false; ny]; nx];
    let mut stack = vec![(0usize, 0usize)];
    outside[0][0] = true;
    while let Some((i, j)) = stack.pop() {
        for (di, dj) in [(0, 1isize), (0, -1), (1, 0), (-1, 0)] {
            let ni = i as isize + di;
            let nj = j as isize + dj;
            if ni < 0 || nj < 0 || ni >= nx as isize || nj >= ny as isize {
                continue;
            }
            let ni = ni as usize;
            let nj = nj as usize;
            if occ_ext[ni][nj] || outside[ni][nj] {
                continue;
            }
            outside[ni][nj] = true;
            stack.push((ni, nj));
        }
    }

    let mut outer_segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    let mut hole_segs: Vec<((f64, f64), (f64, f64))> = Vec::new();

    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            let l = occ_ext[i][j];
            let r = occ_ext[i + 1][j];
            if l == r {
                continue;
            }
            let x = xs[i + 1];
            let y0 = ys[j];
            let y1 = ys[j + 1];
            let seg = ((x, y0), (x, y1));
            if l && !r {
                if outside[i + 1][j] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            } else if !l && r {
                if outside[i][j] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            }
        }
    }

    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            let b = occ_ext[i][j];
            let t = occ_ext[i][j + 1];
            if b == t {
                continue;
            }
            let y = ys[j + 1];
            let x0 = xs[i];
            let x1 = xs[i + 1];
            let seg = ((x0, y), (x1, y));
            if b && !t {
                if outside[i][j + 1] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            } else if !b && t {
                if outside[i][j] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            }
        }
    }

    let scale = 1.0 / t.max(TOLERANCE_LEN_MIN);
    let mut rings = Vec::new();
    if let Some(r0) = segments_to_ring(&outer_segs, scale) {
        rings.push(r0);
    } else {
        return None;
    }
    let mut hs = segments_to_rings(&hole_segs, scale)?;
    rings.append(&mut hs);
    Some(rings)
}

fn build_padded_occ_grid(
    rects: &[(f64, f64, f64, f64)],
    tol: f64,
) -> Option<(Vec<Vec<bool>>, Vec<f64>, Vec<f64>)> {
    let mut xs: Vec<f64> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();
    for &(u0, u1, v0, v1) in rects {
        xs.push(u0);
        xs.push(u1);
        ys.push(v0);
        ys.push(v1);
    }
    xs.sort_by(|a, b| a.total_cmp(b));
    xs.dedup_by(|a, b| (*a - *b).abs() <= tol * 10.0);
    ys.sort_by(|a, b| a.total_cmp(b));
    ys.dedup_by(|a, b| (*a - *b).abs() <= tol * 10.0);
    if xs.len() < 2 || ys.len() < 2 {
        return None;
    }
    let dx = ((xs[xs.len() - 1] - xs[0]).abs()).max(1.0) * TOLERANCE_AREA_REL;
    let dy = ((ys[ys.len() - 1] - ys[0]).abs()).max(1.0) * TOLERANCE_AREA_REL;
    let mut xs_e = vec![xs[0] - dx];
    xs_e.extend(xs.iter().cloned());
    xs_e.push(xs[xs.len() - 1] + dx);
    let mut ys_e = vec![ys[0] - dy];
    ys_e.extend(ys.iter().cloned());
    ys_e.push(ys[ys.len() - 1] + dy);

    let nx = xs_e.len() - 1;
    let ny = ys_e.len() - 1;
    let mut occ = vec![vec![false; ny]; nx];
    for i in 0..nx {
        for j in 0..ny {
            let cu = (xs_e[i] + xs_e[i + 1]) * 0.5;
            let cv = (ys_e[j] + ys_e[j + 1]) * 0.5;
            occ[i][j] = rects.iter().any(|&(u0, u1, v0, v1)| {
                cu + tol >= u0 && cu - tol <= u1 && cv + tol >= v0 && cv - tol <= v1
            });
        }
    }

    let nx2 = nx + 2;
    let ny2 = ny + 2;
    let mut occ_ext = vec![vec![false; ny2]; nx2];
    for i in 0..nx {
        for j in 0..ny {
            occ_ext[i + 1][j + 1] = occ[i][j];
        }
    }

    let mut xs2 = vec![xs_e[0] - dx];
    xs2.extend(xs_e.iter().cloned());
    xs2.push(xs_e[xs_e.len() - 1] + dx);
    let mut ys2 = vec![ys_e[0] - dy];
    ys2.extend(ys_e.iter().cloned());
    ys2.push(ys_e[ys_e.len() - 1] + dy);

    Some((occ_ext, xs2, ys2))
}

fn segments_to_ring(segs: &[((f64, f64), (f64, f64))], scale: f64) -> Option<Vec<(f64, f64)>> {
    let mut rings = segments_to_rings(segs, scale)?;
    if rings.is_empty() {
        return None;
    }
    rings.sort_by(|a, b| {
        ring_area_uv(b)
            .partial_cmp(&ring_area_uv(a))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    rings.into_iter().next()
}

fn segments_to_rings(segs: &[((f64, f64), (f64, f64))], scale: f64) -> Option<Vec<Vec<(f64, f64)>>> {
    if segs.is_empty() {
        return Some(vec![]);
    }
    let mut adj: HashMap<Pt, Vec<Pt>> = HashMap::new();
    let mut undir: HashSet<(Pt, Pt)> = HashSet::new();
    for &((x0, y0), (x1, y1)) in segs {
        let a = qpt(x0, y0, scale);
        let b = qpt(x1, y1, scale);
        let k0 = if a <= b { (a, b) } else { (b, a) };
        if !undir.insert(k0) {
            continue;
        }
        adj.entry(a).or_default().push(b);
        adj.entry(b).or_default().push(a);
    }

    let mut out: Vec<Vec<(f64, f64)>> = Vec::new();
    loop {
        let Some(&v0) = adj.keys().next() else { break };
        if adj.get(&v0).is_none_or(|n| n.is_empty()) {
            adj.remove(&v0);
            continue;
        }
        let Some(ring) = trace_one_ring(&mut adj, v0, scale) else {
            return None;
        };
        if ring.len() >= 3 {
            out.push(ring);
        }
    }
    Some(out)
}

fn remove_adj(adj: &mut HashMap<Pt, Vec<Pt>>, a: Pt, b: Pt) {
    if let Some(v) = adj.get_mut(&a) {
        if let Some(p) = v.iter().position(|&x| x == b) {
            v.remove(p);
        }
        if v.is_empty() {
            adj.remove(&a);
        }
    }
    if let Some(v) = adj.get_mut(&b) {
        if let Some(p) = v.iter().position(|&x| x == a) {
            v.remove(p);
        }
        if v.is_empty() {
            adj.remove(&b);
        }
    }
}

fn trace_one_ring(adj: &mut HashMap<Pt, Vec<Pt>>, v0: Pt, scale: f64) -> Option<Vec<(f64, f64)>> {
    let v1 = adj.get_mut(&v0)?.pop()?;
    remove_adj(adj, v0, v1);
    let to_f = |p: Pt| (p.0 as f64 / scale, p.1 as f64 / scale);
    let mut pts = vec![to_f(v0), to_f(v1)];
    let mut prev = v0;
    let mut cur = v1;
    loop {
        let nbs = adj.get_mut(&cur)?;
        if nbs.is_empty() {
            return None;
        }
        let pos = nbs.iter().position(|&x| x != prev)?;
        let nxt = nbs.remove(pos);
        if nbs.is_empty() {
            adj.remove(&cur);
        }
        remove_adj(adj, cur, nxt);
        if nxt == v0 {
            return Some(pts);
        }
        pts.push(to_f(nxt));
        prev = cur;
        cur = nxt;
        if pts.len() > 100_000 {
            return None;
        }
    }
}

fn ring_area_uv(ring: &[(f64, f64)]) -> f64 {
    if ring.len() < 3 {
        return 0.0;
    }
    let mut a = 0.0;
    for i in 0..ring.len() {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % ring.len()];
        a += x0 * y1 - x1 * y0;
    }
    0.5 * a.abs()
}

fn face_is_plane(brep: &BRep, flat_idx: usize) -> bool {
    brep.geom.face_surface.get(flat_idx)
        .and_then(|&o| o)
        .and_then(|sidx| brep.geom.surfaces.get(sidx))
        .is_some_and(|s| matches!(s, Surface3::Plane(_)))
}

fn flat_face_index(brep: &BRep, si: usize, shi: usize, fi: usize) -> usize {
    let mut idx = 0usize;
    for s in 0..si {
        for sh in &brep.solids[s].shells {
            idx += sh.faces.len();
        }
    }
    for sh in 0..shi {
        idx += brep.solids[si].shells[sh].faces.len();
    }
    idx + fi
}

fn replace_shell_faces_and_geom(
    brep: &mut BRep,
    si: usize,
    shi: usize,
    remove_fis: &[usize],
    merged_face: Face,
    surf_idx: usize,
    flat_indices: &[usize],
) {
    let remove: std::collections::HashSet<usize> = remove_fis.iter().copied().collect();
    let min_keep = *remove_fis.iter().min().unwrap();
    let shell = &mut brep.solids[si].shells[shi];
    let mut new_faces: Vec<Face> = Vec::with_capacity(shell.faces.len() - remove_fis.len() + 1);
    for fi in 0..shell.faces.len() {
        if remove.contains(&fi) {
            if fi == min_keep {
                new_faces.push(merged_face.clone());
            }
        } else {
            new_faces.push(shell.faces[fi].clone());
        }
    }
    shell.faces = new_faces;

    let mut flats: Vec<usize> = flat_indices.to_vec();
    flats.sort_unstable();
    let insert_at = flats[0];
    for &f in flats.iter().rev() {
        if f != insert_at {
            crate::remove_flat_face_geom_slots(&mut brep.geom, f);
        }
    }
    crate::remove_flat_face_geom_slots(&mut brep.geom, insert_at);
    brep.geom.face_surface.insert(insert_at, Some(surf_idx));
    // Insert at insert_at (not resize) so that entries after the removed index
    // are pushed right rather than left-shifted and then padded at the end 鈥?
    // the old resize path caused face_surface_range entries to drift by one
    // position (cylinder wall UV domain 鈫?planar face, H8 box 鈭?cylinder).
    brep.geom.face_surface_range.insert(insert_at, None);
    // NOTE: `insert_at` is a shell-local face index, while
    // remove_flat_face_geom_slots expects geom-store indices.  They can diverge
    // when a shell has fewer faces than the geom store (e.g. multi-shell results
    // after boolean difference).  Re-sync if needed to avoid a debug panic.
    if brep.geom.face_surface_range.len() != brep.geom.face_surface.len() {
        let max_len = brep.geom.face_surface.len().max(brep.geom.face_surface_range.len());
        brep.geom.face_surface.resize(max_len, None);
        brep.geom.face_surface_range.resize(max_len, None);
    }
    // Use resize for both under- and over-shoot to avoid the shell-local
    // insert_at issue (same as above).
    if brep.geom.face_tolerance.len() != brep.geom.face_surface.len() {
        brep.geom
            .face_tolerance
            .resize(brep.geom.face_surface.len(), 0.0);
    } else {
        brep.geom.face_tolerance[insert_at] = 0.0;
    }
}

fn push_new_edges(brep: &mut BRep, edges: Vec<(usize, usize)>) {
    for (a, b) in edges {
        brep.edges.push(Edge { start: a, end: b });
        let ei = brep.edges.len() - 1;
        if brep.geom.edge_curve.len() <= ei {
            brep.geom.edge_curve.resize(ei + 1, None);
        }
        if brep.geom.edge_curve_range.len() <= ei {
            brep.geom.edge_curve_range.resize(ei + 1, None);
        }
        if brep.geom.edge_degenerated.len() <= ei {
            brep.geom.edge_degenerated.resize(ei + 1, false);
        }
        let p0 = brep.vertices[a].point;
        let p1 = brep.vertices[b].point;
        let delta = p1 - p0;
        let len = delta.length();
        let dir = if len > TOLERANCE_LEN_MIN { delta / len } else { DVec3::X };
        let curve_idx = brep.geom.curves.len();
        brep.geom.curves.push(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
            origin: p0,
            direction: dir,
        }));
        brep.geom.edge_curve[ei] = Some(curve_idx);
        brep.geom.edge_curve_range[ei] = Some([0.0, (p1 - p0).dot(dir)]);
        brep.geom.edge_degenerated[ei] = len <= TOLERANCE_LEN_MIN;
    }
}

fn add_vertices_for_rings_with_eval(
    brep: &mut BRep,
    rings: &[Vec<(f64, f64)>],
    eval: impl Fn(f64, f64) -> DVec3,
    tol: f64,
) -> Vec<Vec<usize>> {
    let match_tol = tol.max(TOLERANCE_ABS) * 50.0;
    let mut out: Vec<Vec<usize>> = Vec::new();
    for ring in rings {
        let mut vi: Vec<usize> = Vec::new();
        for &(u, v) in ring {
            let p = eval(u, v);
            let mut found = None;
            for (i, vtx) in brep.vertices.iter().enumerate() {
                if (vtx.point - p).length() <= match_tol {
                    found = Some(i);
                    break;
                }
            }
            let idx = found.unwrap_or_else(|| {
                let i = brep.vertices.len();
                brep.vertices.push(Vertex { point: p });
                i
            });
            vi.push(idx);
        }
        out.push(vi);
    }
    out
}

fn add_vertices_for_rings(
    brep: &mut BRep,
    rings: &[Vec<(f64, f64)>],
    plane: &Plane,
    u_axis: DVec3,
    v_axis: DVec3,
    tol: f64,
) -> Vec<Vec<usize>> {
    let match_tol = tol.max(TOLERANCE_ABS) * 50.0;
    let mut out: Vec<Vec<usize>> = Vec::new();
    for ring in rings {
        let mut vi: Vec<usize> = Vec::new();
        for &(u, v) in ring {
            let p = plane.origin + u_axis * u + v_axis * v;
            let mut found = None;
            for (i, vtx) in brep.vertices.iter().enumerate() {
                if (vtx.point - p).length() <= match_tol {
                    found = Some(i);
                    break;
                }
            }
            let idx = found.unwrap_or_else(|| {
                let i = brep.vertices.len();
                brep.vertices.push(Vertex { point: p });
                i
            });
            vi.push(idx);
        }
        out.push(vi);
    }
    out
}

// Fix try_fuse: use add_vertices_for_rings first, then build wires from vertex indices only

fn face_outer_points(brep: &BRep, face: &Face) -> Vec<DVec3> {
    let mut pts = Vec::new();
    for we in &face.outer_wire.edges {
        let Some((u, _)) = crate::oriented_edge_vertices(brep, *we) else {
            continue;
        };
        if let Some(v) = brep.vertices.get(u) {
            pts.push(v.point);
        }
    }
    pts
}

/// Four UV samples are the four corners of their axis-aligned bounding box (rect), allowing
/// unordered / diagonal wire traversal on plane projections.
fn four_uv_points_on_rectangle_corners(uv: &[(f64, f64)], t: f64) -> bool {
    if uv.len() != 4 {
        return false;
    }
    let (umin, umax, vmin, vmax) = bbox2d(uv);
    if (umax - umin) <= t * 10.0 || (vmax - vmin) <= t * 10.0 {
        return false;
    }
    let corners = [
        (umin, vmin),
        (umax, vmin),
        (umax, vmax),
        (umin, vmax),
    ];
    let mut used = [false; 4];
    for p in uv {
        let mut matched = false;
        for i in 0..4 {
            if used[i] {
                continue;
            }
            let c = corners[i];
            if (p.0 - c.0).abs() <= t * 100.0 && (p.1 - c.1).abs() <= t * 100.0 {
                used[i] = true;
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }
    used.iter().all(|&b| b)
}

fn polygon_is_orthogonal(poly: &[DVec3], to_uv: &dyn Fn(DVec3) -> (f64, f64), tol: f64) -> bool {
    let t = tol.max(TOLERANCE_ABS);
    let uv: Vec<(f64, f64)> = poly.iter().map(|p| to_uv(*p)).collect();
    // Consecutive duplicate vertices in UV (split collinear edges) would make both
    // `du` and `dv` ~0 and incorrectly fail; collapse them first.
    let mut compact: Vec<(f64, f64)> = Vec::with_capacity(uv.len());
    for q in uv {
        if let Some(&last) = compact.last() {
            if (q.0 - last.0).abs() <= t * 10.0 && (q.1 - last.1).abs() <= t * 10.0 {
                continue;
            }
        }
        compact.push(q);
    }
    if compact.len() >= 2 {
        let f = compact[0];
        let l = *compact.last().unwrap();
        if (f.0 - l.0).abs() <= t * 10.0 && (f.1 - l.1).abs() <= t * 10.0 {
            compact.pop();
        }
    }
    let n = compact.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        let (u0, v0) = compact[i];
        let (u1, v1) = compact[(i + 1) % n];
        let du = (u1 - u0).abs();
        let dv = (v1 - v0).abs();
        if du > t && dv > t {
            return false;
        }
        if du <= t && dv <= t {
            return false;
        }
    }
    true
}

fn bbox2d(uv: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for &(u, v) in uv {
        umin = umin.min(u);
        umax = umax.max(u);
        vmin = vmin.min(v);
        vmax = vmax.max(v);
    }
    (umin, umax, vmin, vmax)
}

/// Collapse 180掳 vertices on a closed UV ring from [`union_rects_to_rings_grid`].
fn simplify_ring_collinear_uv_closed(ring: &[(f64, f64)], tol: f64) -> Vec<(f64, f64)> {
    if ring.len() < 3 {
        return ring.to_vec();
    }
    let collinear_abs = (tol * 1e2).max(TOLERANCE_COORD_SUB);
    let mut pts: Vec<(f64, f64)> = ring.to_vec();
    if pts.len() > 1 {
        let a = pts[0];
        let b = *pts.last().unwrap();
        if (a.0 - b.0).abs() + (a.1 - b.1).abs() < (tol * 10.0).max(TOLERANCE_COORD_SUB) {
            pts.pop();
        }
    }
    if pts.len() < 3 {
        return pts;
    }
    let max_rounds = pts.len() * 8 + 8;
    for _ in 0..max_rounds {
        let n = pts.len();
        if n < 3 {
            break;
        }
        let mut remove: Option<usize> = None;
        for i in 0..n {
            let p0 = pts[(i + n - 1) % n];
            let p1 = pts[i];
            let p2 = pts[(i + 1) % n];
            let c = (p1.0 - p0.0) * (p2.1 - p0.1) - (p1.1 - p0.1) * (p2.0 - p0.0);
            if c.abs() <= collinear_abs {
                remove = Some(i);
                break;
            }
        }
        match remove {
            Some(i) => {
                pts.remove(i);
            }
            None => break,
        }
    }
    pts
}

/// Closed ring with only axis-aligned edges (each step moves in u **or** v, not both).
fn ring_is_axis_aligned_orthogonal_uv(compact: &[(f64, f64)], tol: f64) -> bool {
    let t = tol.max(TOLERANCE_ABS);
    let n = compact.len();
    if n < 4 {
        return false;
    }
    for i in 0..n {
        let (u0, v0) = compact[i];
        let (u1, v1) = compact[(i + 1) % n];
        let du = (u1 - u0).abs();
        let dv = (v1 - v0).abs();
        if du > t && dv > t {
            return false;
        }
        if du <= t && dv <= t {
            return false;
        }
    }
    true
}

/// Vertices of a closed ring from [`union_rects_to_rings_grid`] include many collinear samples.
/// This removes 180掳 vertices until only corners remain. A true merged rectangle has 4; an L
/// (two edge-adjacent quads) keeps 6+.
fn ring_corner_count_after_collinear_removal(ring: &[(f64, f64)], tol: f64) -> usize {
    simplify_ring_collinear_uv_closed(ring, tol).len()
}

// 鈹€鈹€ Coplanar overlap clipping for Intersection results 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Ensure a 2D polygon has counter鈥慶lockwise winding (SH expects CCW for the clip polygon).
/// Sutherland鈥慔odgman's `is_inside` uses left鈥憃f鈥慹dge as "inside", which requires the clip
/// polygon to be CCW.  Faces with negative normals (e.g. 鈭抁) produce CW UV projections, so
/// we reverse them here.
fn ensure_ccw(poly: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if poly.len() < 3 {
        return poly.to_vec();
    }
    let mut area2 = 0.0;
    for i in 0..poly.len() {
        let j = (i + 1) % poly.len();
        area2 += poly[i][0] * poly[j][1] - poly[j][0] * poly[i][1];
    }
    if area2 < 0.0 {
        poly.iter().copied().rev().collect()
    } else {
        poly.to_vec()
    }
}

/// For Intersection results, replace partially-overlapping coplanar axis-aligned face
/// pairs with a single face covering their exact 2D overlap polygon.
///
/// `handle_coplanar_faces` in `PaveFiller` records a FaceFace interference but does **not**
/// create intersection curves, so both coplanar faces remain un-split and both get emitted
/// in `build_with_history`.  When the two faces overlap (not strict bbox subset, which
/// [`remove_axis_coplanar_redundant_child_faces`] already handles), both are kept in the
/// shell and the surface area is inflated by the duplicate region.
///
/// This pass finds such pairs, computes the 2D polygon intersection via Sutherland鈥揌odgman,
/// and replaces both faces with a single face covering only the overlap.
pub fn clip_coplanar_overlap_for_intersection(brep: &BRep, a: &BRep, b: &BRep, tol: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let mut total = 0usize;
    let t = tol.max(TOLERANCE_ABS);

    // Phase 1: pair-based clipping 鈥?find coplanar face pairs in the result and
    // replace each pair with a single face covering their 2D polygon overlap.
    for si in 0..out.solids.len() {
        for shi in 0..out.solids[si].shells.len() {
            loop {
                let n = out.solids[si].shells[shi].faces.len();
                if n < 2 {
                    break;
                }
                let mut found = false;
                'pair: for fi in 0..n {
                    for fj in (fi + 1)..n {
                        if clip_one_coplanar_pair(&mut out, si, shi, fi, fj, t).is_some() {
                            total += 1;
                            found = true;
                            break 'pair;
                        }
                    }
                }
                if !found {
                    break;
                }
            }
        }
    }

    // Phase 2: clip remaining faces against input solids.
    // When the boolean classifier removes one of a coplanar pair (e.g. classifying it
    // as "Out"), that pair is invisible to Phase 1 鈥?the surviving face is too large.
    // We build a map of axis-aligned faces from the input solids by plane, then clip
    // each surviving result face against every input face on the same plane.
    let mut input_map: HashMap<(i64, i64, i64, i64), Vec<Vec<[f64; 2]>>> = HashMap::new();
    for input in [a, b] {
        for s in &input.solids {
            for sh in &s.shells {
                for f in &sh.faces {
                    if !f.inner_wires.is_empty() {
                        continue;
                    }
                    let n = snap_almost_axis(f.normal.normalize_or_zero());
                    let Some(axes) = axis_aligned_world_plane_uv_axes(n) else {
                        continue;
                    };
                    let Some(p) = face_first_point(input, f) else {
                        continue;
                    };
                    let d = n.dot(p);
                    let (n_c, d_c) = canonicalize_plane_n_d(n, d);
                    let pk = plane_key(n_c, d_c, t);

                    let [i, j] = axes;
                    let uv: Vec<[f64; 2]> = face_outer_points(input, f)
                        .iter()
                        .map(|p| [p[i], p[j]])
                        .collect();
                    if uv.len() < 3 {
                        continue;
                    }

                    input_map.entry(pk).or_default().push(ensure_ccw(&uv));
                }
            }
        }
    }

    // Pre-scan is no longer needed 鈥?the area check (sa_before / max_input_area < 0.97)
    // correctly distinguishes already-trimmed faces from full coplanar faces.

    for si in 0..out.solids.len() {
        for shi in 0..out.solids[si].shells.len() {
            let mut fi = 0;
            while fi < out.solids[si].shells[shi].faces.len() {
                let mut replaced = false;

                // Read face data for the current result face
                let (face, normal, _axes, point, pk, poly_uv, sa_before) = {
                    let shell = &out.solids[si].shells[shi];
                    if fi >= shell.faces.len() {
                        break;
                    }
                    let f = &shell.faces[fi];
                    if !f.inner_wires.is_empty() {
                        fi += 1;
                        continue;
                    }

                    let n = snap_almost_axis(f.normal.normalize_or_zero());
                    let Some(axes) = axis_aligned_world_plane_uv_axes(n) else {
                        fi += 1;
                        continue;
                    };
                    let Some(p) = face_first_point(&out, f) else {
                        fi += 1;
                        continue;
                    };
                    let d = n.dot(p);
                    let (n_c, d_c) = canonicalize_plane_n_d(n, d);
                    let pk = plane_key(n_c, d_c, t);

                    // Skip curved (non-plane) faces 鈥?the axis-aligned normal check
                    // is not sufficient for cylinder walls that get a planar
                    // classification after boolean face dedup.
                    let flat = flat_face_index(&out, si, shi, fi);
                    let sidx = out.geom.face_surface.get(flat).copied().flatten();
                    if let Some(sidx) = sidx {
                        if let Some(surf) = out.geom.surfaces.get(sidx) {
                            if !matches!(surf, Surface3::Plane(_)) {
                                fi += 1;
                                continue;
                            }
                        }
                    }

                    let [i, j] = axes;
                    let uv: Vec<[f64; 2]> = face_outer_points(&out, f)
                        .iter()
                        .map(|p| [p[i], p[j]])
                        .collect();
                    if uv.len() < 3 {
                        fi += 1;
                        continue;
                    }

                    let sa = face_surface_area(&out, f, 0);

                    (f.clone(), n, axes, p, pk, uv, sa)
                };

                if let Some(input_polys) = input_map.get(&pk) {
                    // Skip faces that were already trimmed by the boolean (area significantly
                    // smaller than all input faces on this plane).  Further clipping would
                    // over-trim (e.g. G5 faces with SA=0.75 on a 1.0 plane).
                    let max_input_area = input_polys
                            .iter()
                            .map(|poly| {
                                let mut a2 = 0.0_f64;
                                for k in 0..poly.len() {
                                    let j = (k + 1) % poly.len();
                                    a2 += poly[k][0] * poly[j][1] - poly[j][0] * poly[k][1];
                                }
                                0.5 * a2.abs()
                            })
                            .reduce(f64::max)
                            .unwrap_or(0.0);
                        if max_input_area > 0.0 && sa_before / max_input_area < 0.97 {
                            fi += 1;
                            continue;
                        }

                    // Clip the result face polygon against each input polygon
                    let mut clipped = poly_uv.clone();
                    for (_pi, input_poly) in input_polys.iter().enumerate() {
                        let sh = crate::inttools::coplanar::sutherland_hodgman_clip(&clipped, input_poly);
                        clipped = sh;
                        if clipped.len() < 3 {
                            break;
                        }
                    }

                    if clipped.len() >= 3 {
                        // Check whether the clipping actually changed the polygon
                        let changed = clipped.len() != poly_uv.len()
                            || clipped.iter().zip(poly_uv.iter()).any(|(a, b)| {
                                (a[0] - b[0]).abs() > t || (a[1] - b[1]).abs() > t
                            });

                        if changed {
                            // Minimum area guard
                            let mut area2 = 0.0_f64;
                            for k in 0..clipped.len() {
                                let (x0, y0) = (clipped[k][0], clipped[k][1]);
                                let (x1, y1) = (clipped[(k + 1) % clipped.len()][0], clipped[(k + 1) % clipped.len()][1]);
                                area2 += x0 * y1 - x1 * y0;
                            }
                            let area = 0.5 * area2.abs();
                            let min_area = (t * t).max(TOLERANCE_FLOAT_ULTRA);

                            if area >= min_area {

                                let rings = vec![clipped
                                    .iter()
                                    .map(|&c| (c[0], c[1]))
                                    .collect::<Vec<_>>()];

                                let ring_vertices =
                                    add_vertices_for_rings_with_eval(&mut out, &rings, |u, v| {
                                        point_from_axis_plane_world_uv(normal, point, u, v)
                                    }, t);

                                if !ring_vertices.is_empty() && ring_vertices[0].len() >= 3 {
                                    let mut edge_pairs: Vec<(usize, usize)> = Vec::new();
                                    for rv in &ring_vertices {
                                        let nv = rv.len();
                                        for k in 0..nv {
                                            edge_pairs.push((rv[k], rv[(k + 1) % nv]));
                                        }
                                    }

                                    let base_ei = out.edges.len();
                                    push_new_edges(&mut out, edge_pairs);

                                    let outer_wire = Wire {
                                        edges: (0..ring_vertices[0].len())
                                            .map(|k| WireEdge::fwd(base_ei + k))
                                            .collect(),
                                    };

                                    let new_face = Face {
                                        outer_wire,
                                        inner_wires: vec![],
                                        normal: face.normal,
                                        triangles: vec![],
                                        sample_point: None,
                                        mesh_dirty: true,
                surface_idx: None,
                                    };

                                    let surf_idx = {
                                        let p0 =
                                            face_first_point(&out, &new_face).unwrap_or(point);
                                        let idx = out.geom.surfaces.len();
                                        out.geom
                                            .surfaces
                                            .push(Surface3::Plane(Plane { origin: p0, normal: face.normal }));
                                        idx
                                    };

                                    let flat_idx = flat_face_index(&out, si, shi, fi);
                                    replace_shell_faces_and_geom(
                                        &mut out,
                                        si,
                                        shi,
                                        &[fi],
                                        new_face,
                                        surf_idx,
                                        &[flat_idx],
                                    );
                                    total += 1;
                                    replaced = true;
                                }
                            }
                        }
                    }
                }

                if !replaced {
                    fi += 1;
                }
            }
        }
    }

    (out, total)
}

/// Try to clip one overlapping coplanar pair.  Returns `Some(())` when a replacement was made.
///
/// # Read-then-write approach
///
/// 1. Read everything from `brep` into local variables while holding only shared references.
/// 2. Compute the SH intersection.
/// 3. Mutate `brep` 鈥?no borrow conflicts because the "read" borrows are dropped.
fn clip_one_coplanar_pair(
    brep: &mut BRep,
    si: usize,
    shi: usize,
    fi: usize,
    fj: usize,
    tol: f64,
) -> Option<()> {
    // Skip curved (non-plane) faces 鈥?a cylinder wall classified with an
    // axis-aligned face normal (e.g. after boolean On鈥慺ace dedup) must not
    // be treated as a planar face.  Clipping its 2D axis projection corrupts
    // face_surface_range and destroys the curved鈥憇urface sub鈥慺ace.
    for &f in &[fi, fj] {
        let flat = flat_face_index(brep, si, shi, f);
        let sidx = brep.geom.face_surface.get(flat).copied().flatten()?;
        if !matches!(brep.geom.surfaces.get(sidx)?, rcad_kernel::geom::Surface3::Plane(_)) {
            return None;
        }
    }

    // 鈹€鈹€ Read phase (shared borrows only) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let (face_i, face_j, n, axes, poly_i_uv, poly_j_uv, plane_origin, normal);

    {
        let shell = &brep.solids[si].shells[shi];
        if fi >= shell.faces.len() || fj >= shell.faces.len() {
            return None;
        }
        face_i = shell.faces[fi].clone();
        face_j = shell.faces[fj].clone();

        // No holes
        if !face_i.inner_wires.is_empty() || !face_j.inner_wires.is_empty() {
            return None;
        }

        // Both must snap to 卤axis
        let n_i = snap_almost_axis(face_i.normal.normalize_or_zero());
        let n_j = snap_almost_axis(face_j.normal.normalize_or_zero());
        let axes_i = axis_aligned_world_plane_uv_axes(n_i)?;
        let axes_j = axis_aligned_world_plane_uv_axes(n_j)?;
        // Must be on the same infinite plane
        if axes_i != axes_j {
            // Different axis families cannot be the same oriented plane.
            // (e.g. one Z-plane and one X-plane are different no matter what.)
            return None;
        }

        let p_i = face_first_point(brep, &face_i)?;
        let p_j = face_first_point(brep, &face_j)?;
        let d_i = n_i.dot(p_i);
        let d_j = n_j.dot(p_j);
        let (n_i_c, d_i_c) = canonicalize_plane_n_d(n_i, d_i);
        let (n_j_c, d_j_c) = canonicalize_plane_n_d(n_j, d_j);
        if plane_key(n_i_c, d_i_c, tol) != plane_key(n_j_c, d_j_c, tol) {
            return None;
        }

        // UV bbox overlap with positive area
        let bi = face_axis_world_bbox(brep, &face_i, n_i)?;
        let bj = face_axis_world_bbox(brep, &face_j, n_j)?;
        let gap = (tol * 1e2).max(TOLERANCE_MESH_LEGACY);
        if !rects_2d_bbox_positive_area_overlap(bi, bj, gap) {
            return None;
        }

        // Skip strict bbox subset 鈥?already handled by remove_axis_coplanar_redundant_child_faces.
        let scale = (bi.1 - bi.0)
            .abs()
            .max(bi.3 - bi.2)
            .abs()
            .max(bj.1 - bj.0)
            .abs()
            .max(bj.3 - bj.2)
            .abs()
            .max(1.0);
        let eps = (TOLERANCE_COORD_SUB * scale).max(tol * TOLERANCE_TOL_SCALE_MICRO);
        let subset = |a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)| -> bool {
            a.0 >= b.0 - eps && a.1 <= b.1 + eps && a.2 >= b.2 - eps && a.3 <= b.3 + eps
        };
        let s_ij = subset(bi, bj) && !subset(bj, bi);
        let s_ji = subset(bj, bi) && !subset(bi, bj);
        if s_ij || s_ji {
            return None;
        }
        if subset(bi, bj) && subset(bj, bi) {
            // Equal bboxes 鈥?also handled by the subset pass.
            return None;
        }

        // Project both face boundaries to world-axis UV
        let [i_axis, j_axis] = axes_i;
        poly_i_uv = face_outer_points(brep, &face_i)
            .iter()
            .map(|p| [p[i_axis], p[j_axis]])
            .collect::<Vec<_>>();
        poly_j_uv = face_outer_points(brep, &face_j)
            .iter()
            .map(|p| [p[i_axis], p[j_axis]])
            .collect::<Vec<_>>();

        if poly_i_uv.len() < 3 || poly_j_uv.len() < 3 {
            return None;
        }

        n = n_i;
        axes = axes_i;
        plane_origin = p_i;
        normal = face_i.normal;
    }

    // 鈹€鈹€ Compute 2D polygon intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // SH expects the clip polygon to be CCW (left-of-edge = inside).  Faces with negative
    // normals (e.g. 鈭抁) project to CW in world-axis UV, so we normalise the clip polygon.
    let poly_j_uv = ensure_ccw(&poly_j_uv);
    let overlap =
        crate::inttools::coplanar::sutherland_hodgman_clip(&poly_i_uv, &poly_j_uv);
    if overlap.len() < 3 {
        return None;
    }

    // Minimum area guard
    {
        let mut area2 = 0.0_f64;
        for k in 0..overlap.len() {
            let (x0, y0) = (overlap[k][0], overlap[k][1]);
            let (x1, y1) = (overlap[(k + 1) % overlap.len()][0], overlap[(k + 1) % overlap.len()][1]);
            area2 += x0 * y1 - x1 * y0;
        }
        let area = 0.5 * area2.abs();
        let min_area = (tol * tol).max(TOLERANCE_FLOAT_ULTRA);
        if area < min_area {
            return None;
        }
    }

    // 鈹€鈹€ Write phase: create face from overlap polygon 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let [_i_ax, _j_ax] = axes;

    // Convert to (f64, f64) rings for add_vertices_for_rings_with_eval
    let overlap_uv: Vec<(f64, f64)> = overlap.iter().map(|&c| (c[0], c[1])).collect();
    let rings = vec![overlap_uv];

    let ring_vertices = add_vertices_for_rings_with_eval(brep, &rings, |u, v| {
        point_from_axis_plane_world_uv(n, plane_origin, u, v)
    }, tol);

    if ring_vertices.is_empty() || ring_vertices[0].len() < 3 {
        return None;
    }

    // Build edges
    let mut edge_pairs: Vec<(usize, usize)> = Vec::new();
    for rv in &ring_vertices {
        let nv = rv.len();
        for k in 0..nv {
            edge_pairs.push((rv[k], rv[(k + 1) % nv]));
        }
    }

    let base_ei = brep.edges.len();
    push_new_edges(brep, edge_pairs);

    // Outer wire
    let n0 = ring_vertices[0].len();
    let outer_wire = Wire {
        edges: (0..n0).map(|k| WireEdge::fwd(base_ei + k)).collect(),
    };

    let merged_face = Face {
        outer_wire,
        inner_wires: vec![],
        normal,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };

    // Plane surface
    let surf_idx = {
        let p0 = face_first_point(brep, &merged_face).unwrap_or(plane_origin);
        let idx = brep.geom.surfaces.len();
        brep.geom
            .surfaces
            .push(Surface3::Plane(Plane { origin: p0, normal }));
        idx
    };

    // Replace both faces with the overlap face
    let remove_fis = [fi, fj];
    let flat_indices: Vec<usize> = remove_fis
        .iter()
        .map(|&f| flat_face_index(brep, si, shi, f))
        .collect();
    replace_shell_faces_and_geom(brep, si, shi, &remove_fis, merged_face, surf_idx, &flat_indices);

    Some(())
}

#[cfg(test)]
mod orth_union_tests {
    use super::ring_corner_count_after_collinear_removal;
    use super::rects_2d_bbox_positive_area_overlap;
    use super::union_rects_to_rings_grid;
    use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_MESH_LEGACY, TOLERANCE_RETRY_LADDER_COARSE};

    #[test]
    fn union_rects_three_adjacent_strips_forms_outer_ring() {
        let rects = [
            (0.0, 5.0, 0.0, 10.0),
            (5.0, 10.0, 0.0, 10.0),
            (10.0, 15.0, 0.0, 10.0),
        ];
        let rings = union_rects_to_rings_grid(&rects, TOLERANCE_ABS);
        assert!(rings.is_some(), "expected grid union to succeed for three strips");
        let rings = rings.unwrap();
        assert!(!rings.is_empty(), "expected at least one ring");
    }

    /// Same bucket key as disjoint islands on one plane: no 2D area overlap in UV.
    #[test]
    fn bbox_positive_area_overlap_distinguishes_disjoint_corner_edge() {
        let gap = TOLERANCE_RETRY_LADDER_COARSE;
        let a = (0.0, 1.0, 0.0, 1.0);
        let b_corner = (2.0, 3.0, 2.0, 3.0);
        assert!(!rects_2d_bbox_positive_area_overlap(a, b_corner, gap));
        let b_edge = (1.0, 2.0, 0.0, 1.0);
        assert!(!rects_2d_bbox_positive_area_overlap(a, b_edge, gap));
        let c_overlap = (0.5, 1.5, 0.0, 1.0);
        assert!(rects_2d_bbox_positive_area_overlap(a, c_overlap, gap));
    }

    /// L-shaped outline keeps >4 corners; a 3脳1 rectangle of samples collapses to 4 corners.
    #[test]
    fn ring_collinear_simplify_rect_vs_l() {
        let tol = TOLERANCE_MESH_LEGACY;
        let l_ring: Vec<(f64, f64)> = vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (0.0, 2.0),
        ];
        assert!(ring_corner_count_after_collinear_removal(&l_ring, tol) >= 5);

        let rect_dense: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)];
        assert_eq!(ring_corner_count_after_collinear_removal(&rect_dense, tol), 4);
    }
}

/// OCCT `bcommon_simple/G1` intersection: document axis-UV bbox overlap vs strict containment between coplanar 卤axis faces.
///
/// Run with `cargo test -p rcad-algorithms g1_intersection_axis_bbox_relationship_probe -- --nocapture` to print counts.
/// Intended for diagnosing +1 `checkprops -s` gaps when duplicate caps do **not** satisfy strict bbox subset.
#[cfg(test)]
mod bcommon_g1_bbox_probe_tests {
    use glam::{DAffine3, DVec3};
    use rcad_kernel::BRep;
    use rcad_modeling::make_box_brep;

    use crate::boolean_op;
    use crate::tolerance::{
        TOLERANCE_ABS, TOLERANCE_COORD_SUB, TOLERANCE_MESH_LEGACY, TOLERANCE_TOL_SCALE_MICRO,
    };
    use crate::BooleanOpType;

    use super::{
        axis_aligned_world_plane_uv_axes, canonicalize_plane_n_d, face_axis_world_bbox, face_first_point,
        plane_key, rects_2d_bbox_positive_area_overlap, snap_almost_axis,
    };

    fn g1_operands() -> (BRep, BRep) {
        let ba = make_box_brep(
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0).normalize(),
            DVec3::new(0.0, 1.0, 0.0).normalize(),
            1.0,
            1.0,
            1.0,
        )
        .expect("ba");
        let bb = make_box_brep(
            DVec3::new(0.0, 0.7071067811865476, 0.0),
            DVec3::new(0.0, 0.0, 1.0).normalize(),
            DVec3::new(0.0, -1.0, 0.0).normalize(),
            1.0,
            0.7071067811865476,
            1.4142135623730951,
        )
        .expect("bb");
        let bb = {
            let mut shape = bb;
            let pivot = DVec3::new(0.0, 0.0, 0.0);
            let axis = DVec3::new(0.0, 0.0, 1.0).normalize_or(DVec3::Z);
            let rot = DAffine3::from_axis_angle(axis, (45.0_f64).to_radians());
            let xf = DAffine3::from_translation(pivot) * rot * DAffine3::from_translation(-pivot);
            shape.apply_transform(xf);
            shape
        };
        (ba, bb)
    }

    #[test]
    fn g1_intersection_axis_bbox_relationship_probe() {
        let (ba, bb) = g1_operands();
        let brep = boolean_op(BooleanOpType::Intersection, &ba, &bb).expect("g1 intersection");
        let t = TOLERANCE_ABS;
        let mut strict_ij = 0usize;
        let mut strict_ji = 0usize;
        let mut overlap_only = 0usize;
        let gap = (t * 1e2).max(TOLERANCE_MESH_LEGACY);

        let subset_containment = |a: (f64, f64, f64, f64), b: (f64, f64, f64, f64), scale: f64| -> bool {
            let eps = (TOLERANCE_COORD_SUB * scale).max(t * TOLERANCE_TOL_SCALE_MICRO);
            let (au0, au1, av0, av1) = a;
            let (bu0, bu1, bv0, bv1) = b;
            au0 >= bu0 - eps && au1 <= bu1 + eps && av0 >= bv0 - eps && av1 <= bv1 + eps
        };

        for si in 0..brep.solids.len() {
            for shi in 0..brep.solids[si].shells.len() {
                let shell = &brep.solids[si].shells[shi];
                let n = shell.faces.len();
                for fi in 0..n {
                    for fj in (fi + 1)..n {
                        let fa = &shell.faces[fi];
                        let fb = &shell.faces[fj];
                        if !fa.inner_wires.is_empty() || !fb.inner_wires.is_empty() {
                            continue;
                        }
                        let n_i = snap_almost_axis(fa.normal.normalize_or_zero());
                        let n_j = snap_almost_axis(fb.normal.normalize_or_zero());
                        if axis_aligned_world_plane_uv_axes(n_i).is_none()
                            || axis_aligned_world_plane_uv_axes(n_j).is_none()
                        {
                            continue;
                        }
                        let Some(p_i) = face_first_point(&brep, fa) else {
                            continue;
                        };
                        let Some(p_j) = face_first_point(&brep, fb) else {
                            continue;
                        };
                        let d_i = n_i.dot(p_i);
                        let d_j = n_j.dot(p_j);
                        let (n_i_c, d_i_c) = canonicalize_plane_n_d(n_i, d_i);
                        let (n_j_c, d_j_c) = canonicalize_plane_n_d(n_j, d_j);
                        let key_i = plane_key(n_i_c, d_i_c, t);
                        let key_j = plane_key(n_j_c, d_j_c, t);
                        if key_i != key_j {
                            continue;
                        }
                        let Some(bi) = face_axis_world_bbox(&brep, fa, n_i) else {
                            continue;
                        };
                        let Some(bj) = face_axis_world_bbox(&brep, fb, n_j) else {
                            continue;
                        };
                        let scale = (bi.1 - bi.0)
                            .abs()
                            .max((bi.3 - bi.2).abs())
                            .max((bj.1 - bj.0).abs())
                            .max((bj.3 - bj.2).abs())
                            .max(1.0);

                        let s_ij =
                            subset_containment(bi, bj, scale) && !subset_containment(bj, bi, scale);
                        let s_ji =
                            subset_containment(bj, bi, scale) && !subset_containment(bi, bj, scale);
                        if s_ij {
                            strict_ij += 1;
                        }
                        if s_ji {
                            strict_ji += 1;
                        }
                        if rects_2d_bbox_positive_area_overlap(bi, bj, gap)
                            && !(subset_containment(bi, bj, scale) && subset_containment(bj, bi, scale))
                            && !s_ij
                            && !s_ji
                        {
                            overlap_only += 1;
                        }
                    }
                }
            }
        }

        eprintln!(
            "G1 intersection axis-plane face pairs (same shell): strict-subset one-way ij={strict_ij} ji={strict_ji}; overlap-not-mutual-subset={overlap_only}"
        );
    }
}

 
