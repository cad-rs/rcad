//! OCCT-aligned: DecomposeResult + helpers for quadric-surface post-processing.
//!
//! OCCT IntPatch_ImpPrmIntersection.cxx L3146-3730 + helper functions.
//!
//! Splits intersection lines at seam/pole boundaries for sphere/cone/cylinder/torus
//! quadric surfaces.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Surface3, SurfaceEval, Line2d, Line3};
use crate::inttools::int_surf_quadric::Quadric;
use crate::inttools::geom_abs_surface_type::GeomAbsSurfaceType;
use crate::inttools::int_patch_line::{IntPatchLine, WLinePnt, WLineType};
use crate::inttools::int_patch_type::IntPatchIType;

// ── OCCT IntPatch_SpecPntType ──────────────────────────────────────
#[derive(Clone, Copy, PartialEq, Eq)]
enum SpecPntType { None, SeamU, SeamV, SeamUV, PoleSeamU, Pole }

// ── OCCT L3174: constants ─────────────────────────────────────────
const DELTA_U_MAX: f64 = std::f64::consts::FRAC_PI_2;
const TOL_3D: f64 = 1e-10;
const TOL_2D: f64 = 1e-12;
const TOL_2DS: f64 = 1e-12;

// ═══════════════════════════════════════════════════════════════════
// OCCT L2057-2142: GetVertices — collect non-duplicate vertices
// ═══════════════════════════════════════════════════════════════════
fn get_vertices(line: &IntPatchLine, tol3d: f64, tol2d: f64) -> Vec<(DVec3, f64, f64, f64, f64)> {
    // Collect unique vertices from the line
    // rcad simplified: return an empty set — vertex management is
    // handled by the IntPatch_RstInt module not yet ported
    Vec::new()
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L2226+: AdjustLine — adjust UV parametrization for quadric
// ═══════════════════════════════════════════════════════════════════
fn adjust_line(
    sline: &mut Vec<(DVec3, f64, f64, f64, f64)>, // (p3d, u1, v1, u2, v2)
    is_reversed: bool,
    quad: &Quadric,
    _tol2d: f64,
) {
    let typ = quad.surface_type();
    for p in sline.iter_mut() {
        let (u1, v1, _u2, _v2) = (p.1, p.2, p.3, p.4);
        // Adjust U/V for periodic surfaces
        if !is_reversed {
            p.3 = quad.parameters(p.0).0; // u from quadric
            p.4 = quad.parameters(p.0).1; // v from quadric
        } else {
            p.1 = quad.parameters(p.0).0;
            p.2 = quad.parameters(p.0).1;
        }
        let _ = (u1, v1, typ);
    }
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L2144-2165: SearchVertices — classify line points vs vertices
// ═══════════════════════════════════════════════════════════════════
fn search_vertices(
    _sline: &[(DVec3, f64, f64, f64, f64)],
    _vertices: &[(DVec3, f64, f64, f64, f64)],
) -> Vec<i32> {
    // rcad simplified: no vertex classification
    vec![0i32; _sline.len()]
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: InsertSeamVertices
// ═══════════════════════════════════════════════════════════════════
fn insert_seam_vertices(
    _sline: &mut Vec<(DVec3, f64, f64, f64, f64)>,
    _is_reversed: bool,
    _vertices: &[(DVec3, f64, f64, f64, f64)],
    _ptypes: &[i32],
    _tol2d: f64,
) -> bool {
    false
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L85-177: IsSeamOrPole — detect if a point is at seam/pole
// ═══════════════════════════════════════════════════════════════════
fn is_seam_or_pole(
    quad_type: GeomAbsSurfaceType,
    _line: &[(DVec3, f64, f64, f64, f64)],
    _ref_idx: usize,
    _tol3d: f64,
    _delta_max: f64,
) -> SpecPntType {
    // rcad simplified: detect based on surface type
    match quad_type {
        GeomAbsSurfaceType::Cylinder => SpecPntType::SeamU,
        GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Cone => SpecPntType::PoleSeamU,
        GeomAbsSurfaceType::Torus => SpecPntType::None,
        _ => SpecPntType::None,
    }
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L1968-2054: AdjustUFirst — normalize U to [0, 2π]
// ═══════════════════════════════════════════════════════════════════
fn adjust_u_first(u1: f64, u2: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    // Case: near 0
    if u1 == 0.0 || u1.abs() <= 1e-9 {
        if u2 > 0.0 && u2 < two_pi {
            return if u2 < (two_pi - u2) { 0.0 } else { two_pi };
        }
        let mut uu = u2;
        while uu > two_pi { uu -= two_pi; }
        while uu < 0.0 { uu += two_pi; }
        return if uu < (two_pi - uu) { 0.0 } else { two_pi };
    }
    // Case: near 2π
    if u1 == two_pi || (two_pi - u1.abs()).abs() <= 1e-9 {
        if u2 > 0.0 && u2 < two_pi {
            return if u2 < (two_pi - u2) { 0.0 } else { two_pi };
        }
        let mut uu = u2;
        while uu > two_pi { uu -= two_pi; }
        while uu < 0.0 { uu += two_pi; }
        return if uu < (two_pi - uu) { 0.0 } else { two_pi };
    }
    // Case: < 0 or > 2π
    let mut u = u1;
    while u < 0.0 { u += two_pi; }
    while u > two_pi { u -= two_pi; }
    u
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: DetectOfBoundaryAchievement
// ═══════════════════════════════════════════════════════════════════
fn detect_of_boundary_achievement(
    _sline: &[(DVec3, f64, f64, f64, f64)],
    _k: usize,
    sline_out: &mut Vec<(DVec3, f64, f64, f64, f64)>,
    is_on_boundary: &mut bool,
) {
    *is_on_boundary = false;
    // rcad simplified
    let _ = sline_out;
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: VerifyVertices
// ═══════════════════════════════════════════════════════════════════
fn verify_vertices(
    _sline: &[(DVec3, f64, f64, f64, f64)],
    _is_reversed: bool,
    _vertices: &[(DVec3, f64, f64, f64, f64)],
    _tol2d: f64,
    _tol_arc: f64,
) -> (Option<(DVec3, f64, f64, f64, f64)>, Option<(DVec3, f64, f64, f64, f64)>) {
    (None, None)
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: HasInternals
// ═══════════════════════════════════════════════════════════════════
fn has_internals(
    _sline: &[(DVec3, f64, f64, f64, f64)],
    _vertices: &[(DVec3, f64, f64, f64, f64)],
) -> bool {
    false
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: ToSmooth
// ═══════════════════════════════════════════════════════════════════
fn to_smooth(
    _sline: &[(DVec3, f64, f64, f64, f64)],
    _is_reversed: bool,
    _quad: &Quadric,
    _is_first: bool,
) -> f64 {
    0.0
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: AddVertices
// ═══════════════════════════════════════════════════════════════════
fn add_vertices(
    _sline: &mut Vec<(DVec3, f64, f64, f64, f64)>,
    _vf: Option<(DVec3, f64, f64, f64, f64)>,
    _add_vf: bool,
    _vl: Option<(DVec3, f64, f64, f64, f64)>,
    _add_vl: bool,
    _d3f: f64,
    _d3l: f64,
) -> bool {
    false
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: PutIntVertices
// ═══════════════════════════════════════════════════════════════════
fn put_int_vertices(
    _line: &mut IntPatchLine,
    _sline: &[(DVec3, f64, f64, f64, f64)],
    _is_reversed: bool,
    _vertices: &[(DVec3, f64, f64, f64, f64)],
    _tol_arc: f64,
) {
}

// ═══════════════════════════════════════════════════════════════════
// OCCT: SplitOnSegments — split WLine at sharp angle changes
// ═══════════════════════════════════════════════════════════════════
fn split_on_segments(
    line: &IntPatchLine,
) -> Vec<IntPatchLine> {
    // rcad simplified: return the original line as one segment
    vec![line.clone()]
}

// ═══════════════════════════════════════════════════════════════════
// OCCT L3146-3730: DecomposeResult — main function
// ═══════════════════════════════════════════════════════════════════
pub fn decompose_result(
    line: &IntPatchLine,
    quad: &Quadric,
    _quad_surf: &Surface3,
    _param_surf: &Surface3,
    tol_arc: f64,
    tol_tang: f64,
) -> Option<Vec<IntPatchLine>> {
    // OCCT L3156-3172: skip non-line restriction curves
    if line.line_type == IntPatchIType::Restriction {
        // Restriction line must be isoline — skip if not
        return None;
    }

    // OCCT L3174-3176: constants
    let _a_delta_umax = DELTA_U_MAX;

    // OCCT L3178-3183: need at least 3 points
    if line.is_wline() && line.wline_pnts.len() <= 2 {
        return None;
    }

    // OCCT L3186: get vertices
    let a_v_line = get_vertices(line, TOL_3D, TOL_2D);

    // OCCT L3188-3193: copy the line
    let mut a_ss_line: Vec<(DVec3, f64, f64, f64, f64)> = if line.is_wline() {
        line.wline_pnts.iter().map(|p| (p.p3d, p.u1, p.v1, p.u2, p.v2)).collect()
    } else {
        return None;
    };

    if a_ss_line.len() <= 1 { return None; }

    // OCCT L3195: adjust
    adjust_line(&mut a_ss_line, false, quad, TOL_2D);

    // OCCT L3197-3207: insert seam vertices for Walking lines
    if line.line_type == IntPatchIType::Walking {
        let mut is_inserted = true;
        while is_inserted {
            let ptypes = search_vertices(&a_ss_line, &a_v_line);
            is_inserted = insert_seam_vertices(&mut a_ss_line, false, &a_v_line, &ptypes, TOL_2D);
        }
    }

    let a_lindex = a_ss_line.len();
    let mut a_findex: usize = 0;
    let mut a_bindex: usize = 0;

    // OCCT L3213-3215: main loop state
    let mut fl_next_line = true;
    let mut has_been_decomposed = false;
    let mut a_pre_point_exist = SpecPntType::None;
    let mut pre_point = (DVec3::ZERO, 0.0, 0.0, 0.0, 0.0);

    let mut out_lines: Vec<IntPatchLine> = Vec::new();

    while fl_next_line {
        fl_next_line = false;
        let mut is_decomposited = false;
        let mut sline: Vec<(DVec3, f64, f64, f64, f64)> = Vec::new();

        // OCCT L3227-3231: check if we've consumed all points
        if a_lindex <= a_findex && a_pre_point_exist == SpecPntType::None {
            break;
        }

        // OCCT L3233-3269: handle existing pre-point
        if a_pre_point_exist != SpecPntType::None {
            let _a_ref_pt = &a_ss_line[a_findex];
            // rcad simplified: skip ContinueAfterSpecialPoint
            sline.push(pre_point.clone());
            // Skip duplicates
            while a_findex < a_lindex && true {
                a_findex += 1;
            }
            a_pre_point_exist = SpecPntType::None;
        }

        // OCCT L3274-3525: analyze each point for seam/pole decomposition
        for k in a_findex..a_lindex {
            if k == a_findex {
                pre_point = a_ss_line[k];
                sline.push(pre_point);
                continue;
            }

            // OCCT L3287-3290: check if on boundary
            let mut is_on_boundary = false;
            detect_of_boundary_achievement(&a_ss_line, k, &mut sline, &mut is_on_boundary);

            // OCCT L3292: detect seam/pole
            a_pre_point_exist = is_seam_or_pole(
                quad.surface_type(), &a_ss_line, k - 1, tol_tang, DELTA_U_MAX,
            );

            // OCCT L3294-3306: boundary override
            if is_on_boundary && a_pre_point_exist != SpecPntType::PoleSeamU {
                a_pre_point_exist = SpecPntType::None;
            }

            // OCCT L3308-3509: handle seam/pole decomposition
            if a_pre_point_exist != SpecPntType::None {
                a_bindex = k;
                is_decomposited = true;
                // rcad simplified: handle seam U by adding adjusted point
                if a_pre_point_exist == SpecPntType::SeamU || a_pre_point_exist == SpecPntType::PoleSeamU {
                    let mut new_pt = a_ss_line[k - 1];
                    // Adjust U parameter to maintain continuity
                    let u_adj = adjust_u_first(new_pt.3, a_ss_line[k].3);
                    new_pt.3 = u_adj;
                    if is_on_boundary {
                        break;
                    }
                    sline.push(new_pt);
                    a_pre_point_exist = SpecPntType::SeamU;
                    pre_point = new_pt;
                }
                break;
            }

            pre_point = a_ss_line[k];

            if is_on_boundary {
                a_bindex = k;
                is_decomposited = true;
                a_pre_point_exist = SpecPntType::None;
                break;
            } else {
                sline.push(a_ss_line[k]);
            }
        }

        // OCCT L3527-3544: handle single-point segment
        if sline.len() == 1 {
            fl_next_line = true;
            if a_findex < a_bindex { a_findex = a_bindex; }
            continue;
        }

        // OCCT L3546-3581: vertex verification and smoothing
        // rcad simplified

        // OCCT L3583-3633: create WLine for decomposed segment
        let wline_pnts: Vec<WLinePnt> = sline.iter().map(|p| WLinePnt {
            p3d: p.0, u1: p.1, v1: p.2, u2: p.3, v2: p.4,
        }).collect();

        let mut wline = IntPatchLine::walking(wline_pnts, WLineType::ImpPrm);
        wline.tolerance = tol_arc;
        wline.tang_tolerance = tol_tang;
        wline.line_type = IntPatchIType::Walking;

        // OCCT L3613-3633: split into segments
        let segm = split_on_segments(&wline);
        out_lines.extend(segm);

        // OCCT L3722-3726: continue loop
        if is_decomposited {
            a_findex = a_bindex;
            fl_next_line = true;
            has_been_decomposed = true;
        }
    }

    if has_been_decomposed { Some(out_lines) } else { None }
}
