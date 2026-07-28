use super::angle_2d::{angle_2d, clock_wise_angle, dir_to_angle};
use super::types::{WireEdgeSource, WireFace, WireOrientation, WireSegment};
use super::wire_splitter::{
    EdgeInfo, are_verts_coincident, find_angle_at, mark_edge_passed, mark_edge_passed_both_dirs,
    select_best_outgoing, world_to_uv,
};
use crate::bop::algo::builder::point_in_polygon_2d;
use crate::bop::ds::*;
use crate::tolerance::*;
use glam::DVec2;
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::PCurve;
use rcad_kernel::geom::*;
use std::collections::{HashMap, VecDeque};
// use crate::classify::{Classification, classify_point};
use super::intres2d::IntRes2dDomain;
use super::types::BooleanOpType;
use crate::bop::int_tools::context::Context;
use crate::bop::int_tools::fclass2d::{CSLibClass2d, CSLibResult, curve2d_nb_samples};

pub(crate) fn refine_angles(
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    segments: &[WireSegment],
    ds: &DS,
    face_idx: usize,
) {
    let vertices: Vec<usize> = smart_map.keys().copied().collect();
    let face_surface = ds.face_surface(face_idx).unwrap_or_else(|| {
    panic!("wire_path: face {} has no surface", face_idx)
});
    for &v in &vertices {
        let Some(infos) = smart_map.get(&v).cloned() else {
            continue;
        };

        let mut cnt_bnd = 0;
        let mut cnt_int = 0;
        let mut a1_bnd = 0.0; // outgoing boundary angle
        let mut a2_bnd = 0.0; // incoming boundary angle

        for ei in &infos {
            if !ei.is_inside {
                cnt_bnd += 1;
                if !ei.in_flag {
                    a1_bnd = ei.angle;
                }
                // outgoing (in_flag=false)
                else {
                    a2_bnd = ei.angle;
                } // incoming (in_flag=true)
            } else {
                cnt_int += 1;
            }
        }

        // OCCT L965-968: only vertices with exactly 2 boundary edges
        if cnt_bnd != 2 {
            continue;
        }

        let a_delta = clock_wise_angle(a2_bnd, a1_bnd);

        // OCCT L970-1000: refine IC outgoing angles
        // Maps edge index  ?refined angle (OCCT aDMSR)
        let mut refined_map: std::collections::HashMap<usize, f64> =
            std::collections::HashMap::new();
        for ei in &infos {
            if ei.is_inside && !ei.in_flag {
                let a_ic = ei.angle;
                let a_da = clock_wise_angle(a2_bnd, a_ic);
                if a_da < a_delta {
                    continue; // OCCT L986-989: already inside boundary sweep
                }

                // OCCT L991: try pcurve-based refinement first
                let b_refined = refine_angle_2d(
                    v,
                    &segments[ei.seg_idx],
                    segments,
                    ds,
                    face_surface,
                    a1_bnd,
                    a2_bnd,
                    a_delta,
                    a_ic,
                );
                if let Some(refined_angle) = b_refined {
                    refined_map.insert(ei.seg_idx, refined_angle);
                } else if cnt_int == 2 {
                    // OCCT L1012: aA = (aA <= aA1) ? (aA1 + Precision::Angular()) : (aA2 - Precision::Angular());
                    let eps = ANGULAR;
                    let new_angle = if a_ic <= a1_bnd {
                        (a1_bnd + eps) % std::f64::consts::TAU
                    } else {
                        (a2_bnd - eps + std::f64::consts::TAU) % std::f64::consts::TAU
                    };
                    refined_map.insert(ei.seg_idx, new_angle);
                }
            }
        }

        if refined_map.is_empty() {
            continue;
        }

        // OCCT L1008-1028: update angles in SmartMap
        if let Some(infos_mut) = smart_map.get_mut(&v) {
            for ei in infos_mut.iter_mut() {
                if let Some(&new_angle) = refined_map.get(&ei.seg_idx) {
                    ei.angle = new_angle;
                    // OCCT L1022-1024: for incoming edges, adjust by PI
                    if ei.in_flag {
                        ei.angle = (new_angle + std::f64::consts::PI) % std::f64::consts::TAU;
                    }
                }
            }
        }
    }
}

/// Path  (BOPAlgo_WireSplitter_1.cxx L359-618).
///    walk path with ClockWiseAngle steering
///
///    per-EdgeInfo passed per vertex,

/// Get the parameter range [t_min, t_max] of a Curve2d.
/// For Trimmed: uses its stored t_min/t_max.
/// For Line: returns [0.0, 1.0] (segment from origin to origin+direction).
/// For Circle: returns [0.0, 2閿滅.
/// For other types: returns [0.0, 1.0].
pub(crate) fn pc_parameter_range(curve: &Curve2d) -> (f64, f64) {
    match curve {
        Curve2d::Trimmed(tc) => (tc.t_min, tc.t_max),
        Curve2d::Circle(_) => (0.0, std::f64::consts::TAU),
        _ => (0.0, 1.0),
    }
}

/// Geom2dInt_GInter  ?intersect a ray with a 2D curve.
/// Returns (param_on_curve, param_on_ray) for all intersections within [t_min, t_max].
///
/// OCCT dispatch (IntCurve_IntCurveCurveGen.gxx L247-815):
///   Line  ?Line/Circle/Ellipse  ?IntConicConic (analytic)
///   Line  ?BSpline/Bezier/other  ?TheIntConicCurveOfGInter (projection+Newton)
pub(crate) fn intersect_ray_curve_2d(
    ray_origin: DVec2,
    ray_dir: DVec2,
    curve: &Curve2d,
    t_min: f64,
    t_max: f64,
) -> Vec<(f64, f64)> {
    let (base, tr_shift) = match curve {
        Curve2d::Trimmed(tc) => (&*tc.curve, tc.t_min),
        _ => (curve, 0.0),
    };
    let t_min_a = t_min - tr_shift;
    let t_max_a = t_max - tr_shift;
    match base {
        Curve2d::Line(line) => {
            let a = ray_dir.x;
            let b = -line.direction.x;
            let c = ray_dir.y;
            let d = -line.direction.y;
            let det = a * d - b * c;
            if det.abs() < TOLERANCE_CLAMP_MIN {
                return vec![];
            }
            let rhs_x = line.origin.x - ray_origin.x;
            let rhs_y = line.origin.y - ray_origin.y;
            let s = (d * rhs_x - b * rhs_y) / det;
            let t = (a * rhs_y - c * rhs_x) / det;
            if s >= 0.0 && t >= t_min_a && t <= t_max_a {
                vec![(t + tr_shift, s)]
            } else {
                vec![]
            }
        }
        Curve2d::Circle(circle) => {
            let oc = ray_origin - circle.center;
            let a_c = ray_dir.dot(ray_dir);
            let b_c = 2.0 * ray_dir.dot(oc);
            let c_c = oc.dot(oc) - circle.radius * circle.radius;
            let disc = b_c * b_c - 4.0 * a_c * c_c;
            if disc < 0.0 {
                return vec![];
            }
            let sqrt_disc = disc.sqrt();
            let mut result = Vec::new();
            for &s in &[
                (-b_c - sqrt_disc) / (2.0 * a_c),
                (-b_c + sqrt_disc) / (2.0 * a_c),
            ] {
                if s >= 0.0 {
                    let p = ray_origin + s * ray_dir;
                    let mut t = (p.y - circle.center.y).atan2(p.x - circle.center.x);
                    if t < 0.0 {
                        t += std::f64::consts::TAU;
                    }
                    if t >= t_min_a && t <= t_max_a {
                        result.push((t + tr_shift, s));
                    }
                }
            }
            // Dedup: if tangent (double root), the last two entries are identical
            if result.len() >= 2 {
                let (t1, _) = result[result.len() - 1];
                let (t2, _) = result[result.len() - 2];
                if (t1 - t2).abs() < 1e-12_f64 {
                    result.pop();
                }
            }
            result
        }
        Curve2d::Ellipse(ellipse) => {
            // Ray: P = O + s*D, s >= 0
            // Ellipse: ((P-C)鐠虹椃/a)^2 + ((P-C)鐠虹椄/b)^2 = 1
            let u = ellipse.major_dir;
            let v = DVec2::new(-u.y, u.x);
            let a_e = ellipse.major_radius;
            let b_e = ellipse.minor_radius;
            let oc = ray_origin - ellipse.center;
            let du = ray_dir.dot(u) / a_e;
            let dv = ray_dir.dot(v) / b_e;
            let ou = oc.dot(u) / a_e;
            let ov = oc.dot(v) / b_e;
            let a_c = du * du + dv * dv;
            let b_c = 2.0 * (du * ou + dv * ov);
            let c_c = ou * ou + ov * ov - 1.0;
            let disc = b_c * b_c - 4.0 * a_c * c_c;
            if disc < 0.0 {
                return vec![];
            }
            let sqrt_disc = disc.sqrt();
            let mut result = Vec::new();
            for &s in &[
                (-b_c - sqrt_disc) / (2.0 * a_c),
                (-b_c + sqrt_disc) / (2.0 * a_c),
            ] {
                if s >= 0.0 {
                    let p = ray_origin + s * ray_dir;
                    let dp = p - ellipse.center;
                    let mut t = dp.y.atan2(dp.x);
                    if t < 0.0 {
                        t += std::f64::consts::TAU;
                    }
                    if t >= t_min_a && t <= t_max_a {
                        result.push((t + tr_shift, s));
                    }
                }
            }
            // Dedup: if tangent (double root)
            if result.len() >= 2 {
                let (t1, _) = result[result.len() - 1];
                let (t2, _) = result[result.len() - 2];
                if (t1 - t2).abs() < 1e-12_f64 {
                    result.pop();
                }
            }
            result
        }
        // OCCT TheIntConicCurveOfGInter / TheIntPCurvePCurveOfGInter:
        //   For non-conic curves, sample curve  ?find nearest point to ray  ?Newton refine.
        _ => {
            const N_SEG: usize = 256;
            let ray_len2 = ray_dir.length_squared();
            if ray_len2 < TOLERANCE_LEN_SQ_DIV_SAFE {
                return vec![];
            }
            let ray_d = ray_dir / ray_len2.sqrt();
            let mut candidates: Vec<(f64, f64)> = Vec::new();
            for i in 0..=N_SEG {
                let t = t_min_a + (t_max_a - t_min_a) * (i as f64) / (N_SEG as f64);
                let p = curve.point_at(t + tr_shift);
                let delta = p - ray_origin;
                let s = delta.dot(ray_d);
                if s < 0.0 {
                    continue;
                }
                let perp = (delta - ray_d * s).length_squared();
                if perp < 1e-10 {
                    let is_dup = candidates.last().map_or(false, |&(lt, _)| {
                        (t - lt).abs() < 1e-9 * (t_max_a - t_min_a + 1.0)
                    });
                    if !is_dup {
                        candidates.push((t, s));
                    }
                }
            }
            if !candidates.is_empty() {
                return candidates
                    .into_iter()
                    .map(|(t, s)| (t + tr_shift, s))
                    .collect();
            }
            let mut best: Option<(f64, f64, f64)> = None;
            for i in 0..=N_SEG {
                let t = t_min_a + (t_max_a - t_min_a) * (i as f64) / (N_SEG as f64);
                let p = curve.point_at(t + tr_shift);
                let delta = p - ray_origin;
                let s = delta.dot(ray_d);
                if s < 0.0 {
                    continue;
                }
                let perp = (delta - ray_d * s).length();
                if best.map_or(true, |(bp, _, _)| perp < bp) {
                    best = Some((perp, t, s));
                }
            }
            if let Some((perp, t0, s0)) = best {
                if perp > 1e-4 {
                    return vec![];
                }
                let eps_der = 1e-7;
                let mut t = t0;
                for _ in 0..20 {
                    let p = curve.point_at(t + tr_shift);
                    let p_hi = curve.point_at((t + eps_der).clamp(t_min_a, t_max_a) + tr_shift);
                    let p_lo = curve.point_at((t - eps_der).clamp(t_min_a, t_max_a) + tr_shift);
                    let der = (p_hi - p_lo) / (2.0 * eps_der);
                    let dp = p - ray_origin;
                    let f = dp.x * ray_d.y - dp.y * ray_d.x;
                    let df = der.x * ray_d.y - der.y * ray_d.x;
                    if df.abs() < TOLERANCE_CLAMP_MIN {
                        break;
                    }
                    let dt = -f / df;
                    t = (t + dt).clamp(t_min_a, t_max_a);
                    if dt.abs() < 1e-12 {
                        let s = (curve.point_at(t + tr_shift) - ray_origin).dot(ray_d);
                        if s >= 0.0 {
                            return vec![(t + tr_shift, s)];
                        }
                        break;
                    }
                }
            }
            vec![]
        }
    }
}

///  ?project a UV point onto a curve to find the nearest parameter.
/// OCCT ref: BRep_Tool::Parameter (returns the parameter of a vertex on an edge's curve).
pub(crate) fn project_uv_to_curve(
    uv: DVec2,
    curve: &Curve2d,
    t_min: f64,
    t_max: f64,
) -> Option<f64> {
    let (base, tr_shift) = match curve {
        Curve2d::Trimmed(tc) => (&*tc.curve, tc.t_min),
        _ => (curve, 0.0),
    };
    match base {
        Curve2d::Line(line) => {
            // Project UV onto line: t = dot(UV - L0, Ld) / |Ld|^2
            let dir = line.direction;
            let denom = dir.dot(dir);
            if denom < TOLERANCE_LEN_SQ_DIV_SAFE {
                return None;
            }
            let t = (uv - line.origin).dot(dir) / denom;
            let t_clamped = t.clamp(t_min - tr_shift, t_max - tr_shift);
            Some(t_clamped + tr_shift)
        }
        Curve2d::Circle(circle) => {
            let mut t = (uv.y - circle.center.y).atan2(uv.x - circle.center.x);
            if t < 0.0 {
                t += std::f64::consts::TAU;
            }
            // Normalize to [t_min, t_min + period) by wrapping
            let period = std::f64::consts::TAU;
            let t_norm = if t < t_min {
                t + period * ((t_min - t) / period).ceil()
            } else if t > t_max {
                t - period * ((t - t_max) / period).floor()
            } else {
                t
            };
            let t_clamped = t_norm.clamp(t_min, t_max);
            Some(t_clamped + tr_shift)
        }
        _ => {
            // Fallback: discrete search for nearest parameter
            const N_SEG: usize = 256;
            let mut best_t = t_min;
            let mut best_d2 = (curve.point_at(t_min) - uv).length_squared();
            for i in 1..=N_SEG {
                let t = t_min + (t_max - t_min) * (i as f64) / (N_SEG as f64);
                let d2 = (curve.point_at(t) - uv).length_squared();
                if d2 < best_d2 {
                    best_d2 = d2;
                    best_t = t;
                }
            }
            Some(best_t)
        }
    }
}

///  ?RefineAngle2D (BOPAlgo_WireSplitter_1.cxx L1032-1124).
///
/// For an IC outgoing edge outside the boundary sweep, compute a refined
/// angle by intersecting the edge's UV pcurve with rays along the boundary
/// directions (aA1 = outgoing, aA2+PI = incoming opposite).  The nearest
/// intersection point inside the sweep gives the corrected angle.
///
/// OCCT algorithm:
///   1. Get edge pcurve and vertex parameter (L1057-1061)
///   2. Determine "other end" parameter direction (L1063)
///   3. For each boundary direction aA1, aA2+M_PI (L1070):
///      a. Create ray from vertex UV
///      b. Intersect ray with edge pcurve (L1080)
///      c. Find furthest intersection within MaxDT of vertex param (L1095)
///      d. Sample curve slightly before intersection (L1110)
///      e. Compute angle and check CWA < aDelta (L1115-1121)
pub(crate) fn refine_angle_2d(
    vertex_idx: usize,
    seg: &WireSegment,
    _segments: &[WireSegment],
    ds: &DS,
    face_surface: &Surface3,
    a1_bnd: f64,
    a2_bnd: f64,
    _a_delta: f64,
    _current_angle: f64,
) -> Option<f64> {
    // OCCT L1057-1061: use vertex parameter on edge's pcurve (BRep_Tool::Parameter).
    //   rcad does NOT store per-vertex-on-edge parameters, but for IC arcs we can
    //   use the pcurve endpoint UV directly instead of re-projecting 3D闁愁偅澹廣 via
    //   world_to_uv (which gives wrong UV at periodic surface singularities).
    // OCCT L1062-1068: get pcurve and range
    let (curve2d, t_min, t_max): (Curve2d, f64, f64) = match &seg.source {
        WireEdgeSource::IntersectionCurve(ci) => {
            let ic = &ds.intersection_curves[*ci];
            if let Some(ref pc) = ic.pcurve_on_a {
                let (t_a, t_b) = pc_parameter_range(pc);
                (pc.clone(), t_a, t_b)
            } else if let Some(ref pc) = ic.pcurve_on_b {
                let (t_a, t_b) = pc_parameter_range(pc);
                (pc.clone(), t_a, t_b)
            } else {
                // Fallback: construct Line from vertex UVs
                let uv_s = world_to_uv(face_surface, ds.vertex_point(seg.start_vertex))?;
                let uv_e = world_to_uv(face_surface, ds.vertex_point(seg.end_vertex))?;
                let dir = uv_e - uv_s;
                if dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                    return None;
                }
                (
                    Curve2d::Line(Line2d {
                        origin: uv_s,
                        direction: dir,
                    }),
                    0.0,
                    1.0,
                )
            }
        }
        WireEdgeSource::DsEdge(_ei) => {
            //  ?L1057: use actual pcurve (BRep_Tool::CurveOnSurface)
            //   from WireSegment when available.  Seam/deg edges on periodic
            //   surfaces store their DoSplitSEAMOnFace pcurves in first_pcurve
            //   (native U side) and second_pcurve (shifted U side).  The
            //   forward flag selects the correct pcurve per orientation.
            let pc = if seg.orientation == WireOrientation::Forward {
                seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref())
            } else {
                seg.second_pcurve.as_ref().or(seg.first_pcurve.as_ref())
            };
            if let Some(pc) = pc {
                (pc.clone(), 0.0, 1.0)
            } else {
                let uv_s = world_to_uv(face_surface, ds.vertex_point(seg.start_vertex))?;
                let uv_e = world_to_uv(face_surface, ds.vertex_point(seg.end_vertex))?;
                let dir = uv_e - uv_s;
                if dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                    return None;
                }
                (
                    Curve2d::Line(Line2d {
                        origin: uv_s,
                        direction: dir,
                    }),
                    0.0,
                    1.0,
                )
            }
        }
        WireEdgeSource::SeamEdge => {
            let pc = if seg.orientation == WireOrientation::Forward {
                seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref())
            } else {
                seg.second_pcurve.as_ref().or(seg.first_pcurve.as_ref())
            };
            if let Some(pc) = pc {
                (pc.clone(), 0.0, 1.0)
            } else {
                let uv_s = world_to_uv(face_surface, ds.vertex_point(seg.start_vertex))?;
                let uv_e = world_to_uv(face_surface, ds.vertex_point(seg.end_vertex))?;
                let dir = uv_e - uv_s;
                if dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                    return None;
                }
                (
                    Curve2d::Line(Line2d {
                        origin: uv_s,
                        direction: dir,
                    }),
                    0.0,
                    1.0,
                )
            }
        }
    };

    //  ?L1060-1061: BRep_Tool::Parameter(aV, aE, myFace).
    //   For DSEdge/SeamEdge: use DSEdge.vertex_params directly.
    //   For IntersectionCurve: vertex param = curve endpoint (t_min or t_max).
    let t_v = match &seg.source {
        WireEdgeSource::DsEdge(ei) => ds.edge_vertex_params(*ei).get(&vertex_idx).copied().unwrap_or_else(|| {
            if vertex_idx == seg.start_vertex {
                t_min
            } else {
                t_max
            }
        }),
        WireEdgeSource::SeamEdge | WireEdgeSource::IntersectionCurve(_) => {
            if vertex_idx == seg.start_vertex {
                t_min
            } else {
                t_max
            }
        }
    };

    //  ?OCCT L1060: vertex UV for ray origin (aGAC1.D0(aTV, aPV)).
    let v_uv = match &seg.source {
        WireEdgeSource::DsEdge(_) | WireEdgeSource::SeamEdge => {
            world_to_uv(face_surface, ds.vertex_point(vertex_idx))
        }
        WireEdgeSource::IntersectionCurve(_) => Some(curve2d.point_at(t_v)),
    }
    .unwrap_or(DVec2::ZERO);

    //  ?OCCT L1063-1065: determine "other end" direction and MaxDT
    let t_op = if (t_v - t_min).abs() < (t_v - t_max).abs() {
        t_max
    } else {
        t_min
    };
    let max_dt = 0.3 * (t_max - t_min);
    let a_tol_int = 1e-10;
    let a_cf = 0.01;

    //  ?OCCT L1080-1082: create IntRes2dDomain for the curve (aDomain1).
    let p1 = curve2d.point_at(t_min);
    let p2 = curve2d.point_at(t_max);
    let mut domain_curve = IntRes2dDomain::new();
    domain_curve.set_values_bounded(p1, t_min, a_tol_int, p2, t_max, a_tol_int);

    //  ?OCCT L1070: try both boundary directions (aA1, aA2+M_PI)
    let a_delta = clock_wise_angle(a2_bnd, a1_bnd);
    for i in 0..2 {
        let a_ai = if i == 0 {
            a1_bnd
        } else {
            a2_bnd + std::f64::consts::PI
        };
        let ray_dir = DVec2::new(a_ai.cos(), a_ai.sin());
        if ray_dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
            continue;
        }

        //  ?OCCT L1084-1094: create ray line + call Geom2dInt_GInter.
        let ray_line = Curve2d::Line(Line2d {
            origin: v_uv,
            direction: ray_dir,
        });
        let mut domain_ray = IntRes2dDomain::new(); // infinite domain (no bounds)
        // OCCT uses Geom2dInt_GInter::Perform with two domains.
        let hits = crate::bop::algo::builder::intersection::intersect_curves_2d_ginter(
            &curve2d,
            &domain_curve,
            &ray_line,
            &domain_ray,
            a_tol_int,
            a_tol_int,
        );
        // hits: (param_on_curve, param_on_ray)  ?swap to (t_on_curve, t_on_ray)
        let hits: Vec<(f64, f64)> = hits.into_iter().map(|(tc, tr)| (tr, tc)).collect();

        if hits.is_empty() {
            continue;
        }

        //  ?OCCT L1100-1114: find best intersection (max param_on_ray, within MaxDT)
        let mut best: Option<(f64, f64)> = None;
        for &(t_c, t_r) in &hits {
            let is_better = match best {
                Some((_, best_r)) => t_r > best_r,
                None => true,
            };
            if is_better && (t_c - t_v).abs() < max_dt {
                best = Some((t_c, t_r));
            }
        }

        if let Some((t_1max, _t_2max)) = best {
            // OCCT L1104-1108: skip if intersection is at far end
            let dt = t_op - t_1max;
            if dt.abs() < a_tol_int {
                continue;
            }

            // OCCT L1110-1113: sample curve slightly before intersection
            let t_sample = t_1max + a_cf * dt;
            let p_sample = curve2d.point_at(t_sample);
            let dir = p_sample - v_uv;
            if dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                continue;
            }

            // OCCT L1115-1121: compute angle and check if inside boundary wedge
            let a_angle = dir.y.atan2(dir.x);
            let a_angle = if a_angle < 0.0 {
                a_angle + std::f64::consts::TAU
            } else {
                a_angle
            };
            let a_da = clock_wise_angle(a2_bnd, a_angle);
            if a_da < a_delta {
                return Some(a_angle);
            }
        }
    }
    None
}
///  ?Walk a path extracting closed wires (BOPAlgo_WireSplitter_1.cxx L359-618).
///
/// Key differences from the previous implementation:
/// 1. Tracks UV coordinates of each visited vertex (aCoordVa).
/// 2. Loop detection uses 2D UV distance for closed/degenerate vertices,
///    preventing false loops at seam/IC junctions on periodic surfaces.
/// 3. Sequence truncation matches OCCT L488-521.
pub(crate) fn walk_path_extract_wires(
    start_si: usize,
    segments: &[WireSegment],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    wires: &mut Vec<Vec<usize>>,
    ds: &DS,
    face_idx: usize,
) {
    let start_seg = &segments[start_si];
    // If this segment has no EdgeInfo, it cannot be walked.
    // Mark a dummy EdgeInfo as passed so the outer loop skips it.
    let has_info = smart_map
        .values()
        .any(|v| v.iter().any(|ei| ei.seg_idx == start_si));
    if std::env::var("RCAD_DEBUG_IC").is_ok()
        && matches!(ds.face_surface(face_idx).unwrap_or_else(|| {
            panic!("wire_path: face {} has no surface", face_idx)
        }), Surface3::Sphere(_))
    {
        let seg_src = match &segments[start_si].source {
            WireEdgeSource::DsEdge(ei) => format!("Ds({})", ei),
            WireEdgeSource::IntersectionCurve(ci) => format!("IC({})", ci),
            WireEdgeSource::SeamEdge => "Seam".to_string(),
        };
        eprintln!(
            "[WALK_START] face={} si={} seg={} start_v={} end_v={} fwd={:?} seam={} has_info={}",
            face_idx,
            start_si,
            seg_src,
            segments[start_si].start_vertex,
            segments[start_si].end_vertex,
            segments[start_si].orientation,
            segments[start_si].is_closed_on_face,
            has_info
        );
    }
    if !has_info {
        smart_map
            .entry(start_seg.start_vertex)
            .or_default()
            .push(EdgeInfo {
                seg_idx: start_si,
                passed: true,
                in_flag: false,
                is_inside: false,
                is_circle_arc: false,
                angle: 0.0,
            });
        smart_map
            .entry(start_seg.end_vertex)
            .or_default()
            .push(EdgeInfo {
                seg_idx: start_si,
                passed: true,
                in_flag: true,
                is_inside: false,
                is_circle_arc: false,
                angle: 0.0,
            });
        return;
    }

    let face_surface = ds.face_surface(face_idx).unwrap_or_else(|| {
    panic!("wire_path: face {} has no surface", face_idx)
});
    let two_pi = std::f64::consts::TAU;

    // OCCT: aLS (edge sequence), aVertVa (vertex sequence), aCoordVa (UV coordinates)
    // OCCT L389: anInfoSeq (EdgeInfo pointer sequence for isBoundary tracking)
    let mut edge_seq: Vec<usize> = Vec::new();
    let mut vert_seq: Vec<usize> = Vec::new();
    let mut uv_seq: Vec<DVec2> = Vec::new();
    let mut info_seq: Vec<usize> = Vec::new();

    let mut ci = start_si;
    let mut arrived_vertex = start_seg.end_vertex;
    let mut current_vertex = start_seg.start_vertex;
    let max_iter = segments.len() * 4 + 200; // increased safety limit

    // Build a per-vertex map: does this vertex belong to a closed/degenerate edge?
    // OCCT L424: bIsClosed = aVertMap.Find(aVb)
    let is_vert_closed = |smart_map: &IndexMap<usize, Vec<EdgeInfo>>, v: usize| -> bool {
        smart_map.get(&v).map_or(false, |infos| {
            infos.iter().any(|ei| {
                let seg = &segments[ei.seg_idx];
                seg.start_vertex == seg.end_vertex || seg.is_closed_on_face
            })
        })
    };

    // Coord2d (BOPAlgo_WireSplitter_1.cxx L663-674).
    // Gets UV of a vertex on a specific edge by evaluating the edge's pcurve
    // at the vertex parameter.  Different edges at the same 3D vertex can
    // return DIFFERENT UVs if their pcurves are on different sides of the
    // parametric seam (e.g. U=0 vs U=2 ?on a sphere).
    let vertex_uv = |vi: usize, segment: &WireSegment, at_start: bool| -> Option<DVec2> {
        // Use pcurve-based UV when available (OCCT Coord2d path)
        let pc_uv = match &segment.source {
            WireEdgeSource::IntersectionCurve(ci) => {
                let ic = &ds.intersection_curves[*ci];
                let pc = ic.pcurve_on_a.as_ref().or(ic.pcurve_on_b.as_ref())?;
                // OCCT BRep_Tool::Parameter(aV, aE, aF): vertex parameter on
                // edge's pcurve.  vi == ic.start_vertex  ?t_range[0];
                // vi == ic.end_vertex  ?t_range[1].
                //  ?compare by 3D position, not index.  rcad's DS
                //   assigns different vertex indices to the same 3D point (remap_ic_v),
                //   so vi == ic.start_vertex fails silently for remapped vertices.
                //   Use geometric distance at remap_ic_v's tolerance.
                let vi_at_pole = ds.vertex_point(vi);
                let t = if ds.vertex_point(ic.start_vertex)
                    .distance_squared(vi_at_pole)
                    <= TOLERANCE_ABS_SQ * 1_000_000.0
                {
                    ic.t_range[0]
                } else {
                    ic.t_range[1]
                };
                Some(pc.point_at(t))
            }
            WireEdgeSource::DsEdge(_) if segment.is_closed_on_face => {
                //  ?Coord2d (WireSplitter_1.cxx L663-674) uses the
                //   edge's own pcurve, selected by orientation per CurveOnSurface
                //   (BRep_Tool.cxx L354-361): FORWARD  ?PCurve (native U side),
                //   REVERSED  ?PCurve2 (shifted U side).  rcad models a closed
                //   seam edge as a FWD/REV WireSegment pair; the REVERSED segment
                //   carries the shifted pcurve in `second_pcurve`.
                //
                //   A degenerate pole edge (start==end) is a self-loop that bridges
                //   the parametric seam at the pole.  Its UV goes from (0, Vpole) at
                //   the "out" end to (2 ? Vpole) at the "in" end, spanning the full
                //   U circle at Vpole  ?exactly matching OCCT's pcurve for a sphere
                //   degenerated edge.
                //  ?CurveOnSurface returns PCurve for FORWARD (L354-361),
                //   PCurve2 for REVERSED.  vertex_uv uses first_pcurve (PCurve) for
                //   FORWARD segments, second_pcurve (PCurve2) for REVERSED, matching
                //   Coord2d per-edge pcurve evaluation (WireSplitter_1.cxx L663-674).
                //   Self-loop deg edges store a full-span line in second_pcurve.
                if segment.start_vertex == segment.end_vertex {
                    match &segment.second_pcurve {
                        Some(Curve2d::Line(l)) => {
                            let t = if at_start {
                                segment.t_range[0]
                            } else {
                                segment.t_range[1]
                            };
                            Some(l.point_at(t))
                        }
                        _ => {
                            // OCCT: Coord2d always expects a pcurve  ?fall back to
                            // world_to_uv when unavailable (e.g. degenerated edge).
                            world_to_uv(face_surface, ds.vertex_point(vi))
                        }
                    }
                } else if segment.orientation == WireOrientation::Forward {
                    match (&segment.first_pcurve, &segment.second_pcurve) {
                        (Some(Curve2d::Line(l)), _) => {
                            let t = if at_start {
                                segment.t_range[0]
                            } else {
                                segment.t_range[1]
                            };
                            Some(l.point_at(t))
                        }
                        _ => world_to_uv(face_surface, ds.vertex_point(vi)),
                    }
                } else {
                    // for REVERSED seam traversal, use second_pcurve
                    //   (shifted pcurve).  Fall back to world_to_uv when unavailable
                    //   (e.g. degenerated seam edge at sphere pole).
                    match &segment.second_pcurve {
                        Some(Curve2d::Line(l)) => {
                            let t = if at_start {
                                segment.t_range[0]
                            } else {
                                segment.t_range[1]
                            };
                            Some(l.point_at(t))
                        }
                        _ => world_to_uv(face_surface, ds.vertex_point(vi)),
                    }
                }
            }
            _ => None,
        };
        if let Some(uv) = pc_uv {
            return Some(uv);
        }
        //  ?non-seam DsEdge vertex_uv from first_pcurve (OCCT:
        //   BRep_Tool::CurveOnSurface  ?C2D->D0(BRep_Tool::Parameter(aV,aE,aF), aP2D1)).
        if let WireEdgeSource::DsEdge(_) = &segment.source {
            if !segment.is_closed_on_face {
                if let Some(pc) = &segment.first_pcurve {
                    let t = if at_start {
                        segment.t_range[0]
                    } else {
                        segment.t_range[1]
                    };
                    return Some(pc.point_at(t));
                }
            }
        }

        // OCCT: Coord2d always expects a valid pcurve  ?this fallback should never
        // be reached in OCCT (the edge would not be in the wire).  Release builds
        // use world_to_uv as a best-effort approximation.
        let v_pt = ds.vertex_point(vi);
        match face_surface {
            Surface3::Sphere(s) => Some(s.world_to_uv(v_pt)),
            Surface3::Cylinder(c) => {
                let ax = c.axis.normalize_or_zero();
                let v = (v_pt - c.origin).dot(ax);
                let to_axis = v_pt - (c.origin + ax * v);
                let u = to_axis
                    .dot(c.ref_dir)
                    .atan2(to_axis.dot(c.ref_dir.cross(ax)));
                Some(DVec2::new(u, v))
            }
            Surface3::Plane(p) => {
                let x_axis = p.u_dir;
                let y_axis = p.v_dir;
                let local = v_pt - p.origin;
                Some(DVec2::new(local.dot(x_axis), local.dot(y_axis)))
            }
            _ => None,
        }
    };

    // OCCT Tolerance2D/UTolerance2D/VTolerance2D (BOPAlgo_WireSplitter_1.cxx L859-901).
    let vtol = |vi: usize| -> f64 { ds.vertex_tolerance(vi).max(TOLERANCE_ABS) };
    let u_resolution = |vt: f64| -> f64 {
        match face_surface {
            Surface3::Sphere(s) => vt / s.radius.max(TOLERANCE_CLAMP_MIN),
            Surface3::Cylinder(c) => vt / c.radius.max(TOLERANCE_CLAMP_MIN),
            Surface3::Cone(_) => vt * 1e-3,
            Surface3::Torus(t) => vt / t.major_radius.max(TOLERANCE_CLAMP_MIN),
            _ => vt,
        }
    };
    let v_resolution = |vt: f64| -> f64 {
        match face_surface {
            Surface3::Sphere(s) => vt / s.radius.max(TOLERANCE_CLAMP_MIN),
            Surface3::Cylinder(_) => vt,
            Surface3::Cone(_) => vt,
            Surface3::Torus(t) => vt / t.minor_radius.max(TOLERANCE_CLAMP_MIN),
            _ => vt,
        }
    };
    // OCCT L859-881: Tolerance2D  ?max(UResolution, VResolution, tolV3D)
    let tolerance_2d = |vi: usize| -> f64 {
        let vt = vtol(vi);
        let mut t2d = u_resolution(vt).max(v_resolution(vt)).max(vt);
        if matches!(face_surface, Surface3::BSpline(_) | Surface3::Bezier(_)) {
            t2d *= 1.1;
        }
        t2d
    };
    // OCCT L885-891: UTolerance2D = UResolution(aTolV3D)
    let u_tolerance_2d = |vi: usize| -> f64 { u_resolution(vtol(vi)) };
    // OCCT L895-901: VTolerance2D = VResolution(aTolV3D)
    let v_tolerance_2d = |vi: usize| -> f64 { v_resolution(vtol(vi)) };
    let uv_tolerance = |vi: usize| -> f64 { 2.0 * tolerance_2d(vi) };

    for _iter in 0..max_iter {
        // Do not escape through edge from which you enter
        if edge_seq.len() == 1 {
            let same_edge = match (&segments[edge_seq[0]].source, &segments[ci].source) {
                (WireEdgeSource::DsEdge(ea), WireEdgeSource::DsEdge(eb)) => ea == eb,
                (WireEdgeSource::IntersectionCurve(ca), WireEdgeSource::IntersectionCurve(cb)) => {
                    ca == cb
                }
                (WireEdgeSource::SeamEdge, WireEdgeSource::SeamEdge) => true,
                _ => false,
            };
            if ci == edge_seq[0] || same_edge {
                return;
            }
        }

        // Mark edge as passed
        let seg = &segments[ci];
        mark_edge_passed(smart_map, ci, seg.start_vertex, false);

        edge_seq.push(ci);
        vert_seq.push(current_vertex);
        let cur_uv = vertex_uv(current_vertex, seg, true);
        uv_seq.push(cur_uv.unwrap_or(DVec2::ZERO));
        info_seq.push(ci);

        // Loop Detection */
        let b_is_closed = is_vert_closed(smart_map, arrived_vertex);
        let a_tol_2d = uv_tolerance(arrived_vertex);
        let a_tol_2d_sq = a_tol_2d * a_tol_2d;
        let a_pb = vertex_uv(arrived_vertex, &segments[ci], false).unwrap_or(DVec2::ZERO);

        let mut b_has_edge = false;
        let a_nb = edge_seq.len();
        for i in (0..a_nb).rev() {
            let prev_v = vert_seq[i];
            let prev_uv = uv_seq[i];
            let prev_si = edge_seq[i];

            // OCCT L449-458: bHasEdge  ?skip degenerate-only wires
            if !b_has_edge {
                b_has_edge = match &segments[prev_si].source {
                    WireEdgeSource::DsEdge(ei) => !ds.is_edge_degenerated(*ei),
                    _ => true,
                };
                if !b_has_edge {
                    continue;
                }
            }

            let is_same_v = prev_v == arrived_vertex;
            let mut is_same_v_2d = is_same_v;

            if is_same_v {
                if b_is_closed {
                    // 2D distance check for closed vertices
                    let a_d2 = prev_uv.distance_squared(a_pb);
                    is_same_v_2d = a_d2 < a_tol_2d_sq;
                    if is_same_v_2d {
                        let u_dist = (prev_uv.x - a_pb.x).abs();
                        let v_dist = (prev_uv.y - a_pb.y).abs();
                        let a_tol_u = 2.0 * u_tolerance_2d(arrived_vertex);
                        let a_tol_v = 2.0 * v_tolerance_2d(arrived_vertex);
                        if u_dist > a_tol_u || v_dist > a_tol_v {
                            is_same_v_2d = false;
                        }
                    }
                }
            }

            if is_same_v && is_same_v_2d {
                // Extract wire from edge_seq[i..]
                let wire: Vec<usize> = edge_seq[i..].to_vec();

                // Skip 2-edge wires where both edges are the same
                let mut is_valid = true;
                if wire.len() == 2 {
                    let a = &segments[wire[0]];
                    let b = &segments[wire[1]];
                    let same_edge = match (&a.source, &b.source) {
                        (WireEdgeSource::DsEdge(ea), WireEdgeSource::DsEdge(eb)) => ea == eb,
                        (
                            WireEdgeSource::IntersectionCurve(ca),
                            WireEdgeSource::IntersectionCurve(cb),
                        ) => ca == cb,
                        (WireEdgeSource::SeamEdge, WireEdgeSource::SeamEdge) => true,
                        _ => false,
                    };
                    if same_edge {
                        is_valid = false;
                    }
                }
                if is_valid {
                    wires.push(wire);
                }

                let a_nbj = i;
                if a_nbj == 0 {
                    edge_seq.clear();
                    vert_seq.clear();
                    uv_seq.clear();
                    return;
                }

                // Keep first a_nbj entries, truncate the rest
                let continue_vertex = vert_seq[i];
                edge_seq.truncate(a_nbj);
                vert_seq.truncate(a_nbj);
                uv_seq.truncate(a_nbj);
                info_seq.truncate(a_nbj);

                ci = *info_seq.last().unwrap();
                arrived_vertex = continue_vertex;
                break;
            }
        }

        // Outgoing Edge Selection
        let angle_in = match find_angle_at(smart_map, ci, arrived_vertex, true) {
            Some(a) => a,
            None => return,
        };

        let raw_candidates: Vec<&EdgeInfo> = if let Some(infos) = smart_map.get(&arrived_vertex) {
            infos
                .iter()
                .filter(|ei| !ei.passed && !ei.in_flag)
                .collect()
        } else {
            return;
        };

        // OCCT L571-582: 2D distance check (Coord2dVf vs aPb) applies to ALL
        //   candidates.  Compute a_pb (UV of arrived vertex on current edge) and
        //   b_is_closed before candidate filtering/selection.
        let b_is_closed = is_vert_closed(smart_map, arrived_vertex);
        let a_pb = vertex_uv(arrived_vertex, &segments[ci], false).unwrap_or(DVec2::ZERO);
        let a_tol_2d_sq = {
            let tol = uv_tolerance(arrived_vertex);
            tol * tol
        };

        let raw_candidates: Vec<&EdgeInfo> = if let Some(infos) = smart_map.get(&arrived_vertex) {
            infos
                .iter()
                .filter(|ei| !ei.passed && !ei.in_flag)
                .collect()
        } else {
            return;
        };

        let i_cnt = raw_candidates.len();
        if i_cnt == 0 {
            return;
        }

        // Single candidate  ?take it before 2D UV check (which may reject it at seam)
        if i_cnt == 1 {
            let best = raw_candidates[0];
            current_vertex = arrived_vertex;
            ci = best.seg_idx;
            arrived_vertex = segments[ci].end_vertex;
            continue;
        }

        // For closed vertices, filter multi-candidates by 2D UV distance
        let candidates: Vec<&EdgeInfo> = if b_is_closed {
            let is_periodic = matches!(
                ds.face_surface(face_idx).unwrap_or_else(|| {
            panic!("wire_path: face {} has no surface", face_idx)
        }),
                Surface3::Sphere(_) | Surface3::Cylinder(_) | Surface3::Torus(_)
            );
            let pb_uv = if is_periodic {
                world_to_uv(&ds.face_surface(face_idx).unwrap_or_else(|| {
            panic!("wire_path: face {} has no surface", face_idx)
        }), ds.vertex_point(arrived_vertex))
                    .unwrap_or(DVec2::ZERO)
            } else {
                a_pb
            };
            raw_candidates
                .into_iter()
                .filter(|ei| {
                    let cand_uv = if is_periodic {
                        world_to_uv(&ds.face_surface(face_idx).unwrap_or_else(|| {
            panic!("wire_path: face {} has no surface", face_idx)
        }), ds.vertex_point(arrived_vertex))
                            .unwrap_or(DVec2::ZERO)
                    } else {
                        vertex_uv(arrived_vertex, &segments[ei.seg_idx], true)
                            .unwrap_or(DVec2::ZERO)
                    };
                    let a_d2 = cand_uv.distance_squared(pb_uv);
                    a_d2 < a_tol_2d_sq
                })
                .collect()
        } else {
            raw_candidates
        };

        if candidates.is_empty() {
            return;
        }

        // Use is_inside flag from EdgeInfo for boundary check (OCCT: !anEdgeInfo->IsInside())
        let incoming_info = smart_map
            .get(&arrived_vertex)
            .and_then(|infos| infos.iter().find(|ei| ei.seg_idx == ci && ei.in_flag));
        let incoming_is_boundary = incoming_info.map_or(false, |ei| !ei.is_inside);
        let best = match select_best_outgoing(&candidates, angle_in, incoming_is_boundary, ci) {
            Some(e) => e,
            None => return,
        };

        current_vertex = arrived_vertex;
        ci = best.seg_idx;
        arrived_vertex = segments[ci].end_vertex;
    }
    // Mark visited segments as passed when max_iter is exhausted
    for &si in &edge_seq {
        mark_all_edge_infos_passed(smart_map, si);
    }
}

/// Mark ALL EdgeInfo entries for a segment as passed (both in_flag values).
pub(crate) fn mark_all_edge_infos_passed(
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    seg_idx: usize,
) {
    for infos in smart_map.values_mut() {
        for ei in infos.iter_mut() {
            if ei.seg_idx == seg_idx {
                ei.passed = true;
            }
        }
    }
}

///  ?wire 3D boundary polygon
///     DS  3D

///  ?classify wires into growth/outer and holes
/// (BOPAlgo_BuilderFace::PerformAreas L387-606).
///
/// OCCT creates a TopoDS_Face from each wire via BRepBuilderAPI_MakeFace,
/// then uses IntTools_FClass2d to test if a sample point is IsHole().
/// Growth wires (sample point is NOT in a hole) form the outer boundary.
///
/// rcad equivalent: map 3D wire boundary to UV space, build a UV polygon,
/// then use ray-casting point-in-polygon.  Full-wrap wires (<3 unique
/// vertices, spanning the full periodic domain) use the surface's full
/// UV rectangle as their polygon.
///  ?merge sphere wires by interleaving seam+IC segments.
///    OCCT's DoSplitSEAMOnFace produces a single wire alternating between
///    seam sub-segments and IC arcs.  rcad produces 2 wires (one IC-loop,
///    one seam-loop) on the same vertices but opposite directions.
///    This function interleaves them: seam闁愁偅澧矯闁愁偅濮眅am闁愁偅澧矯闁愁偅濮眅am闁愁偅澧矯.

/// Compute the signed area of a UV polygon using the shoelace formula.
/// Used for sorting wires by size  ?the largest wire is the outer boundary.
pub(crate) fn uv_polygon_area(poly: &[DVec2]) -> f64 {
    if poly.len() < 3 {
        return 0.0;
    }
    let n = poly.len();
    (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            poly[i].x * poly[j].y - poly[j].x * poly[i].y
        })
        .sum::<f64>()
        .abs()
        * 0.5
}

/// Test whether a UV point is inside a UV polygon using the ray casting method.
/// Handles periodic U wrapping for values in [0, 2pi).
pub(crate) fn point_in_uv_polygon(pt: DVec2, poly: &[DVec2]) -> bool {
    if poly.len() < 3 {
        return false;
    }
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let j = (n + i - 1) % n;
        let vi = poly[i];
        let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x
        {
            inside = !inside;
        }
    }
    inside
}

/// Compute the projected area of a 3D polygon onto the given coordinate plane.
pub(crate) fn projected_area_on(b: &[DVec3], u_idx: usize, v_idx: usize) -> f64 {
    let pick = |p: DVec3, i: usize| -> f64 {
        match i {
            0 => p.x,
            1 => p.y,
            _ => p.z,
        }
    };
    (0..b.len())
        .map(|i| {
            let j = (i + 1) % b.len();
            pick(b[i], u_idx) * pick(b[j], v_idx) - pick(b[j], u_idx) * pick(b[i], v_idx)
        })
        .sum::<f64>()
        .abs()
        * 0.5
}

/// Compute the maximum projected area across XY, YZ, and XZ planes.
pub(crate) fn projected_area_max(b: &[DVec3]) -> f64 {
    let xy = projected_area_on(b, 0, 1);
    let yz = projected_area_on(b, 1, 2);
    let xz = projected_area_on(b, 0, 2);
    xy.max(yz).max(xz)
}

/// Test whether a point projects inside a polygon on the XY plane.
/// Falls back to YZ or XZ if the polygon is degenerate in XY.
pub(crate) fn point_in_polygon_best(pt: DVec3, poly: &[DVec3]) -> bool {
    let xy_area = projected_area_on(poly, 0, 1);
    if xy_area > TOLERANCE_CLAMP_MIN {
        return point_in_polygon_xy_impl(pt, poly, 0, 1);
    }
    let yz_area = projected_area_on(poly, 1, 2);
    if yz_area > TOLERANCE_CLAMP_MIN {
        return point_in_polygon_xy_impl(pt, poly, 1, 2);
    }
    point_in_polygon_xy_impl(pt, poly, 0, 2) // XZ fallback
}

/// Ray casting point-in-polygon test in the given 2D projection (u,v).
pub(crate) fn point_in_polygon_xy_impl(
    pt: DVec3,
    poly: &[DVec3],
    u_idx: usize,
    v_idx: usize,
) -> bool {
    let pu = |p: DVec3| -> f64 {
        match u_idx {
            0 => p.x,
            1 => p.y,
            _ => p.z,
        }
    };
    let pv = |p: DVec3| -> f64 {
        match v_idx {
            0 => p.x,
            1 => p.y,
            _ => p.z,
        }
    };
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let j = (n + i - 1) % n;
        let vi = poly[i];
        let vj = poly[j];
        if ((pv(vi) > pv(pt)) != (pv(vj) > pv(pt)))
            && pu(pt) < (pu(vj) - pu(vi)) * (pv(pt) - pv(vi)) / (pv(vj) - pv(vi)) + pu(vi)
        {
            inside = !inside;
        }
    }
    inside
}

/// Legacy: XY-only projection (replaced by projected_area_max / point_in_polygon_best).
/// Kept for callers that explicitly need XY projection.
pub(crate) fn projected_area_xy(b: &[DVec3]) -> f64 {
    projected_area_on(b, 0, 1)
}
