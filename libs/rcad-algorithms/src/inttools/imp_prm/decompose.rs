//! OCCT-aligned: DecomposeResult + helpers for quadric-surface post-processing.
//!
//! OCCT IntPatch_ImpPrmIntersection.cxx L3146-3730 + helper functions.
//! Splits intersection lines at seam/pole boundaries for sphere/cone/cylinder/torus.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, CurveEval, Surface3, SurfaceEval, Line3};
use crate::inttools::int_surf_quadric::Quadric;
use crate::inttools::geom_abs_surface_type::GeomAbsSurfaceType;
use crate::inttools::int_patch_line::{IntPatchLine, WLinePnt, WLineType};
use crate::inttools::int_patch_type::IntPatchIType;

// ── OCCT IntPatch_SpecPntType ──────────────────────────────────────
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum SpecPntType { None, SeamU, SeamV, SeamUV, PoleSeamU, Pole }

const DELTA_U_MAX: f64 = std::f64::consts::FRAC_PI_2;
const TOL_3D: f64 = 1e-10;
const TOL_2D: f64 = 1e-12;

// ── internal point type ────────────────────────────────────────────
#[derive(Clone, Copy)]
struct Pt {
    p3d: DVec3, u1: f64, v1: f64, u2: f64, v2: f64,
}

fn pt_from_wp(wp: &WLinePnt) -> Pt { Pt { p3d: wp.p3d, u1: wp.u1, v1: wp.v1, u2: wp.u2, v2: wp.v2 } }

// ═══════════════════════════════════════════════════════════════════
// OCCT L2057-2142: GetVertices — collect non-duplicate vertices
// ═══════════════════════════════════════════════════════════════════
fn get_vertices(line: &IntPatchLine) -> Vec<Pt> {
    // For each pair of consecutive wline_pnts, check if 3D distance
    // AND all 4 UV coords are within tolerance. If so, they are the
    // same vertex — skip the duplicate.
    let mut verts: Vec<Pt> = Vec::new();
    for i in 0..line.wline_pnts.len() {
        let pi = pt_from_wp(&line.wline_pnts[i]);
        let mut dup = false;
        for k in i + 1..line.wline_pnts.len() {
            let pk = pt_from_wp(&line.wline_pnts[k]);
            if pi.p3d.distance(pk.p3d) <= TOL_3D
                && (pi.u1 - pk.u1).abs() <= TOL_2D
                && (pi.v1 - pk.v1).abs() <= TOL_2D
                && (pi.u2 - pk.u2).abs() <= TOL_2D
                && (pi.v2 - pk.v2).abs() <= TOL_2D
            {
                dup = true; break;
            }
        }
        if !dup { verts.push(pi); }
    }
    verts
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L2144-2165: SearchVertices — classify line points vs vertices
// ═══════════════════════════════════════════════════════════════════
fn search_vertices(sline: &[Pt], vertices: &[Pt]) -> Vec<i32> {
    // For each point in sline, find which vertex it corresponds to
    // (or 0 if none)
    let mut ptypes = vec![0i32; sline.len()];
    for (ip, p) in sline.iter().enumerate() {
        for (iv, v) in vertices.iter().enumerate() {
            if p.p3d.distance(v.p3d) <= TOL_3D
                && (p.u1 - v.u1).abs() <= TOL_2D
                && (p.v1 - v.v1).abs() <= TOL_2D
                && (p.u2 - v.u2).abs() <= TOL_2D
                && (p.v2 - v.v2).abs() <= TOL_2D
            {
                ptypes[ip] = iv as i32 + 1; // 1-indexed
                break;
            }
        }
    }
    ptypes
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: InsertSeamVertices — add seam vertices where needed
// ═══════════════════════════════════════════════════════════════════
fn insert_seam_vertices(
    sline: &mut Vec<Pt>,
    _is_reversed: bool,
    vertices: &[Pt],
    ptypes: &[i32],
    _tol2d: f64,
) -> bool {
    // Check if any vertex lies between two consecutive sline points
    // but is not already a point in sline. If so, insert it.
    let mut inserted = false;
    let mut i = 0;
    while i < sline.len() - 1 {
        let p1 = sline[i];
        let p2 = sline[i + 1];
        for v in vertices {
            // Check if v is between p1 and p2 in the line
            let along = (v.p3d - p1.p3d).dot(p2.p3d - p1.p3d);
            let dist2 = (p2.p3d - p1.p3d).length_squared();
            if along > 0.0 && along < dist2 {
                let ratio = along / dist2;
                // Check UV is intermediate
                let u_mid = p1.u1 + ratio * (p2.u1 - p1.u1);
                let v_mid = p1.v1 + ratio * (p2.v1 - p1.v1);
                if (v.u1 - u_mid).abs() < _tol2d && (v.v1 - v_mid).abs() < _tol2d {
                    // Check not already in ptypes
                    let already = ptypes.iter().any(|&t| t > 0);
                    if !already {
                        sline.insert(i + 1, Pt {
                            p3d: v.p3d, u1: v.u1, v1: v.v1, u2: v.u2, v2: v.v2,
                        });
                        inserted = true;
                        break;
                    }
                }
            }
        }
        i += 1;
    }
    inserted
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L2226+: AdjustLine — adjust UV parametrization for quadric
// ═══════════════════════════════════════════════════════════════════
fn adjust_line(sline: &mut [Pt], is_reversed: bool, quad: &Quadric) {
    for p in sline.iter_mut() {
        let (uq, vq) = quad.parameters(p.p3d);
        if !is_reversed {
            p.u2 = uq; p.v2 = vq;
        } else {
            p.u1 = uq; p.v1 = vq;
        }
    }
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L85-177: IsSeamOrPole — detect seam/pole at line index
// ═══════════════════════════════════════════════════════════════════
fn is_seam_or_pole(
    quad_type: GeomAbsSurfaceType,
    line: &[Pt],
    ref_idx: usize,
    tol3d: f64,
    _delta_max: f64,
) -> SpecPntType {
    if ref_idx + 1 >= line.len() { return SpecPntType::None; }

    // OCCT L116-135: sphere pole detection
    if quad_type == GeomAbsSurfaceType::Sphere || quad_type == GeomAbsSurfaceType::Cone {
        let a_p3d = line[ref_idx + 1].p3d;
        let d2 = a_p3d.distance_squared(line[ref_idx].p3d);
        // If consecutive points are very close, near singular
        if d2 < tol3d * tol3d { return SpecPntType::PoleSeamU; }
        return SpecPntType::PoleSeamU;
    }

    // OCCT L137-177: seam detection for cylinder/torus
    let a_du = (line[ref_idx].u1 - line[ref_idx + 1].u1).abs();
    let a_dv = (line[ref_idx].v1 - line[ref_idx + 1].v1).abs();
    // For quadric, check U parameter jump (periodic boundary)
    let a_uq_ref = line[ref_idx].u2;
    let a_uq_next = line[ref_idx + 1].u2;
    let du_q = (a_uq_ref - a_uq_next).abs();

    match quad_type {
        GeomAbsSurfaceType::Cylinder => {
            if du_q > std::f64::consts::FRAC_PI_2 * 0.5 { SpecPntType::SeamU }
            else { SpecPntType::None }
        }
        GeomAbsSurfaceType::Torus => {
            if du_q > std::f64::consts::FRAC_PI_2 {
                if a_dv > std::f64::consts::FRAC_PI_2 { SpecPntType::SeamUV }
                else { SpecPntType::SeamU }
            } else { SpecPntType::None }
        }
        GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Cone => {
            if du_q > std::f64::consts::FRAC_PI_2 * 0.5 { SpecPntType::PoleSeamU }
            else { SpecPntType::None }
        }
        _ => {
            let _ = (a_du, a_dv);
            SpecPntType::None
        }
    }
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L1968-2054: AdjustUFirst
// ═══════════════════════════════════════════════════════════════════
fn adjust_u_first(u1: f64, u2: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    if u1 == 0.0 || u1.abs() <= 1e-9 {
        if u2 > 0.0 && u2 < two_pi {
            return if u2 < (two_pi - u2) { 0.0 } else { two_pi };
        }
        let mut uu = u2;
        while uu > two_pi { uu -= two_pi; } while uu < 0.0 { uu += two_pi; }
        return if uu < (two_pi - uu) { 0.0 } else { two_pi };
    }
    if u1 == two_pi || (two_pi - u1.abs()).abs() <= 1e-9 {
        if u2 > 0.0 && u2 < two_pi {
            return if u2 < (two_pi - u2) { 0.0 } else { two_pi };
        }
        let mut uu = u2;
        while uu > two_pi { uu -= two_pi; } while uu < 0.0 { uu += two_pi; }
        return if uu < (two_pi - uu) { 0.0 } else { two_pi };
    }
    let mut u = u1;
    while u < 0.0 { u += two_pi; } while u > two_pi { u -= two_pi; }
    u
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: VerifyVertices
// ═══════════════════════════════════════════════════════════════════
fn verify_vertices(
    sline: &[Pt], vertices: &[Pt],
) -> (Option<Pt>, bool, Option<Pt>, bool) {
    let mut add_vf = false; let mut add_vl = false;
    let mut vf = None; let mut vl = None;

    if sline.is_empty() || vertices.is_empty() { return (None, false, None, false); }

    let fst = sline[0]; let lst = sline[sline.len() - 1];
    // Check if first/last point matches any vertex within tolerance
    for v in vertices {
        if fst.p3d.distance(v.p3d) <= TOL_3D { add_vf = true; vf = Some(*v); }
        if lst.p3d.distance(v.p3d) <= TOL_3D { add_vl = true; vl = Some(*v); }
    }
    (vf, add_vf, vl, add_vl)
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: HasInternals
// ═══════════════════════════════════════════════════════════════════
fn has_internals(sline: &[Pt], vertices: &[Pt]) -> bool {
    for p in sline {
        for v in vertices {
            if p.p3d.distance(v.p3d) <= TOL_3D { return true; }
        }
    }
    false
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: SplitOnSegments — split at sharp angle changes
// ═══════════════════════════════════════════════════════════════════
fn split_on_segments(wp: &[WLinePnt], angle_tol: f64) -> Vec<Vec<WLinePnt>> {
    if wp.len() < 3 { return vec![wp.to_vec()]; }

    let mut segments: Vec<Vec<WLinePnt>> = Vec::new();
    let mut cur = vec![wp[0], wp[1]];

    for i in 2..wp.len() {
        let d1 = (wp[i - 1].p3d - wp[i - 2].p3d).normalize_or_zero();
        let d2 = (wp[i].p3d - wp[i - 1].p3d).normalize_or_zero();
        let angle = d1.dot(d2).clamp(-1.0, 1.0).acos();

        if angle > angle_tol && cur.len() >= 2 {
            cur.push(wp[i - 1]);
            segments.push(cur);
            cur = vec![wp[i - 1]];
        }
        cur.push(wp[i]);
    }

    if cur.len() >= 2 { segments.push(cur); }
    if segments.is_empty() { vec![wp.to_vec()] } else { segments }
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L3146-3730: DecomposeResult
// ═══════════════════════════════════════════════════════════════════
pub fn decompose_result(
    line: &IntPatchLine,
    quad: &Quadric,
    _quad_surf: &Surface3,
    _param_surf: &Surface3,
    tol_arc: f64,
    tol_tang: f64,
) -> Option<Vec<IntPatchLine>> {
    // OCCT L3156-3172
    if line.line_type == IntPatchIType::Restriction {
        if !line.is_wline() { return None; }
    }

    // OCCT L3178-3183
    if !line.is_wline() || line.wline_pnts.len() <= 2 { return None; }

    // OCCT L3186: get vertices
    let a_v_line = get_vertices(line);

    // OCCT L3188: copy the line
    let mut a_ss_line: Vec<Pt> = line.wline_pnts.iter().map(pt_from_wp).collect();
    if a_ss_line.len() <= 1 { return None; }

    // OCCT L3195: adjust
    adjust_line(&mut a_ss_line, false, quad);

    // OCCT L3197-3207: insert seam vertices
    if line.line_type == IntPatchIType::Walking {
        let mut is_inserted = true;
        while is_inserted {
            let ptypes = search_vertices(&a_ss_line, &a_v_line);
            is_inserted = insert_seam_vertices(&mut a_ss_line, false, &a_v_line, &ptypes, TOL_2D);
        }
    }

    let a_lindex = a_ss_line.len();
    let mut a_findex = 0usize;
    let mut a_bindex = 0usize;
    let mut fl_next_line = true;
    let mut has_been_decomposed = false;
    let mut a_pre_point_exist = SpecPntType::None;
    let mut pre_point = a_ss_line[0];
    let mut out_lines: Vec<IntPatchLine> = Vec::new();

    while fl_next_line {
        fl_next_line = false;
        let mut is_decomposited = false;
        let mut sline: Vec<Pt> = Vec::new();

        if a_lindex <= a_findex && a_pre_point_exist == SpecPntType::None { break; }

        // OCCT L3233-3269: handle pre-point
        if a_pre_point_exist != SpecPntType::None {
            sline.push(pre_point);
            let mut dup = true;
            while dup && a_findex < a_lindex {
                if pre_point.p3d.distance(a_ss_line[a_findex].p3d) >= tol_tang {
                    dup = false; break;
                }
                a_findex += 1;
            }
            a_pre_point_exist = SpecPntType::None;
        }

        // OCCT L3274-3525: analyze points
        for k in a_findex..a_lindex {
            if k == a_findex { pre_point = a_ss_line[k]; sline.push(pre_point); continue; }

            // OCCT L3292: detect seam/pole
            a_pre_point_exist = is_seam_or_pole(quad.surface_type(), &a_ss_line, k - 1, tol_tang, DELTA_U_MAX);

            if a_pre_point_exist != SpecPntType::None {
                a_bindex = k; is_decomposited = true;
                if a_pre_point_exist == SpecPntType::SeamU || a_pre_point_exist == SpecPntType::PoleSeamU {
                    let mut new_pt = a_ss_line[k - 1];
                    new_pt.u2 = adjust_u_first(new_pt.u2, a_ss_line[k].u2);
                    new_pt.v2 = adjust_u_first(new_pt.v2, a_ss_line[k].v2);
                    sline.push(new_pt);
                    a_pre_point_exist = SpecPntType::SeamU;
                    pre_point = new_pt;
                }
                break;
            }

            pre_point = a_ss_line[k];
            sline.push(pre_point);
        }

        // OCCT L3527-3544
        if sline.len() == 1 {
            fl_next_line = true;
            if a_findex < a_bindex { a_findex = a_bindex; }
            continue;
        }

        // OCCT L3559: check internal vertices
        let _has_int = has_internals(&sline, &a_v_line);

        // OCCT L3583: create WLine
        let wps: Vec<WLinePnt> = sline.iter().map(|p| WLinePnt {
            p3d: p.p3d, u1: p.u1, v1: p.v1, u2: p.u2, v2: p.v2,
        }).collect();

        // OCCT L3614: split on segments
        let segs = split_on_segments(&wps, DELTA_U_MAX);
        for seg in segs {
            let mut wl = IntPatchLine::walking(seg, WLineType::ImpPrm);
            wl.tolerance = tol_arc; wl.tang_tolerance = tol_tang;
            out_lines.push(wl);
        }
        has_been_decomposed = true;

        // OCCT L3722-3726
        if is_decomposited {
            a_findex = a_bindex;
            fl_next_line = true;
        }
    }

    if has_been_decomposed { Some(out_lines) } else { None }
}

// ═══════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════

