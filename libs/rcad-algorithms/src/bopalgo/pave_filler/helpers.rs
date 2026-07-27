use crate::bnd_lib;
use crate::bopds::ds::{
    DS, DSCurveRepOnFace, DSEdge, Interference, IntersectionCurve, NearTangentType, ShapeOrigin,
};
use crate::bopds::pave::*;
use crate::inttools;
use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::closest_point_on_curve;
use rcad_kernel::geom::*;
use std::collections::HashSet;

/// IntPatch_Intersection surface category (L1264-1294).
///   ParamParam = ts1==ts2==0 -> PrmPrmIntersection (parametric-parametric)

// ---- Phase 2a helpers: vertex -> curve parameter projection ----

/// ✅ OCCT-aligned: delegates to rcad_kernel::closest_point_on_curve (128 samples + analytic dispatch).
pub(crate) fn project_vertex_to_curve(pt: DVec3, curve: &Curve3, tol: f64) -> Option<f64> {
    let result = closest_point_on_curve(curve, pt, 128);
    if result.distance <= tol {
        Some(result.param)
    } else {
        None
    }
}

// ---- FindValidRange / ShrunkData helper functions ----
// OCCT references:
//   IntTools_ShrunkRange::Perform()  IntTools_ShrunkRange.cxx L107-191
//   BRepLib::FindValidRange          BRepLib_1.cxx L173-258
//   findNearestValidPoint            BRepLib_1.cxx L31-148

/// Curve parameter step: parameter increment needed to move tol distance along curve.
/// OCCT: Adaptor3d_Curve::Resolution(theTol) (BRepLib_1.cxx L61, IntTools_ShrunkRange.cxx L162)
/// Note: rcad uses `tol` directly in the formula (tol / speed), while OCCT also
/// applies `* 1.01` in findNearestValidPoint (L61).
pub(crate) fn curve_resolution(curve: &Curve3, t: f64, tol: f64) -> f64 {
    use rcad_kernel::geom::CurveEval;
    let speed = curve.tangent_at(t).length();
    if speed < TOLERANCE_CLAMP_MIN {
        tol
    } else {
        tol / speed
    }
}

/// (core logic): findNearestValidPoint (BRepLib_1.cxx L31-148)
/// Step along the curve from one end until outside the vertex tolerance sphere,
/// then binary-search to refine the exit parameter.
///
/// OCCT differences:
/// 1. OCCT uses `theCurve.Resolution(theTol) * 1.01` (L61) -- rcad now adds * 1.01.
/// 2. OCCT has BSpline/Bezier specific handling (aD1Mag threshold, L70-81) to
///    accelerate through near-singular derivative regions -- rcad does not implement this.
/// 3. OCCT checks `aP.SquareDistance(theVertPnt) > aSqTol` as the exit condition -- rcad matches.
/// 4. OCCT mid-point refinement exits when `aDelta <= theEps` -- rcad matches.
pub(crate) fn find_nearest_valid_point(
    curve: &Curve3,
    first: f64,
    last: f64,
    is_first: bool,
    vert_pt: DVec3,
    vert_tol: f64,
    eps: f64,
) -> Option<f64> {
    use rcad_kernel::geom::CurveEval;
    let (start_u, end_u) = if is_first {
        (first, last)
    } else {
        (last, first)
    };
    let tol_sq = vert_tol * vert_tol;

    // 1. Check if endpoint is inside tolerance sphere
    if curve.point_at(start_u).distance_squared(vert_pt) > tol_sq {
        return None;
    }

    // 2. Step until outside tolerance sphere (OCCT L55-63: `theCurve.Resolution(theTol) * 1.01`)
    let step = (curve_resolution(curve, start_u, vert_tol) * 1.01).max(eps);
    let step = if is_first { step } else { -step };

    // OCCT L70-82: aD1Mag threshold for BSpline/Bezier singularity detection
    //   aD1Mag = (1. / theCurve.Resolution(1.)) * 0.01;  aD1Mag *= aD1Mag;
    //   = (speed_at_start * 0.01)^2  since curve_resolution(c,u,1) = 1/speed
    let a_d1_mag = if matches!(curve, Curve3::BSpline(_) | Curve3::Bezier(_)) {
        let speed = curve.tangent_at(start_u).length();
        let d1mag = speed * 0.01;
        d1mag * d1mag
    } else {
        0.0
    };

    let (mut u_in, mut u_out) = (start_u, start_u);
    let mut is_out = false;
    while !is_out {
        u_in = u_out;
        u_out += step;
        if (is_first && u_out > end_u) || (!is_first && u_out < end_u) {
            if curve.point_at(end_u).distance_squared(vert_pt) <= tol_sq {
                return None;
            }
            u_out = end_u;
            is_out = true;
            break;
        }
        if a_d1_mag > 0.0 {
            // OCCT L108-137: singularity handling — double step through near-zero derivative regions
            let mut step_local = step.abs();
            loop {
                let a_p = curve.point_at(u_out);
                let is_out_local = a_p.distance_squared(vert_pt) > tol_sq;
                if is_out_local {
                    is_out = true;
                    break;
                }
                let a_d1 = curve.derivative_at(u_out);
                if a_d1.length_squared() < a_d1_mag {
                    step_local *= 2.0;
                    let new_out = if is_first {
                        u_out + step_local
                    } else {
                        u_out - step_local
                    };
                    if (is_first && new_out > end_u) || (!is_first && new_out < end_u) {
                        // went out of range — check if end point is outside
                        if curve.point_at(end_u).distance_squared(vert_pt) <= tol_sq {
                            return None;
                        }
                        u_out = end_u;
                        is_out = true;
                        break;
                    }
                    u_out = new_out;
                    continue;
                }
                if is_out_local {
                    is_out = true;
                }
                break;
            }
        } else {
            is_out = curve.point_at(u_out).distance_squared(vert_pt) > tol_sq;
        }
    }

    // 3. Bisection refinement
    while (u_out - u_in).abs() > eps {
        let mid = (u_in + u_out) * 0.5;
        if curve.point_at(mid).distance_squared(vert_pt) > tol_sq {
            u_out = mid;
        } else {
            u_in = mid;
        }
    }
    Some(if is_first { u_out } else { u_in })
}

/// BRepLib::FindValidRange (BRepLib_1.cxx L173-258)
/// Compute the valid (shrunk) range of curve segment [t0, t1] excluding endpoint tolerance spheres.
/// `theTolE`  ?edge tolerance used in Resolution (OCCT L201: curve.Resolution(theTolE * 0.1)).
/// Returns (first, last); returns None if fully covered by tolerance spheres (micro edge).
pub(crate) fn find_valid_range(
    curve: &Curve3,
    t0: f64,
    t1: f64,
    theTolE: f64,
    sv_pt: DVec3,
    sv_tol: f64,
    ev_pt: DVec3,
    ev_tol: f64,
) -> Option<(f64, f64)> {
    use rcad_kernel::geom::CurveEval;
    if (t1 - t0).abs() < rcad_kernel::tolerance::CONFUSION {
        return None;
    }
    let abs_max = t0.abs().max(t1.abs()).max(1.0);
    let eps = curve_resolution(curve, (t0 + t1) * 0.5, theTolE * 0.1)
        .max(abs_max * f64::EPSILON)
        .max(rcad_kernel::tolerance::CONFUSION);
    let first = if t0.is_infinite() {
        t0
    } else {
        match find_nearest_valid_point(curve, t0, t1, true, sv_pt, sv_tol, eps) {
            Some(f) => {
                if t1 - f < eps {
                    return None;
                }
                f
            }
            None => {
                return None;
            }
        }
    };
    let last = if t1.is_infinite() {
        t1
    } else {
        match find_nearest_valid_point(curve, t0, t1, false, ev_pt, ev_tol, eps) {
            Some(l) => {
                if l - t0 < eps {
                    return None;
                }
                l
            }
            None => {
                return None;
            }
        }
    };
    if first > last {
        None
    } else {
        Some((first, last))
    }
}

// ---- Seam Edge Shift Struct ----

/// Result of checking whether a seam edge shift is needed between two faces.
/// BOPAlgo_PaveFiller_6.cxx L393-479
pub(crate) struct SeamEdgeShift {
    /// Translation vector to apply to one face's surface.
    pub(crate) shift_vector: DVec3,
    /// Distance of the shift (used for tolerance contribution).
    pub(crate) shift_value: f64,
    /// Which face is shifted: 1 = f1, 2 = f2.
    pub(crate) shifted_face: u8,
}

// ---- Free Helper Functions ----

/// Apply a translation to a surface's position.
/// The shift modifies the surface's origin (or center) so that the surface
/// appears to move in 3D space. Surface normals and parameterization are
/// preserved.
///
/// gp_Trsf.SetTranslation -- moving the face before intersection
pub(crate) fn apply_shift_to_surface(surface: &Surface3, shift: DVec3) -> Surface3 {
    match *surface {
        Surface3::Plane(p) => Surface3::Plane(Plane {
            origin: p.origin + shift,
            ..p
        }),
        Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
            origin: c.origin + shift,
            ..c
        }),
        Surface3::Sphere(s) => Surface3::Sphere(SphericalSurface {
            center: s.center + shift,
            ..s
        }),
        Surface3::Torus(t) => Surface3::Torus(ToroidalSurface {
            center: t.center + shift,
            ..t
        }),
        Surface3::Cone(c) => Surface3::Cone(ConicalSurface {
            apex: c.apex + shift,
            ..c
        }),
        Surface3::BSpline(ref bs) => {
            let mut bs = bs.clone();
            for row in &mut bs.control_points {
                for cp in row {
                    *cp += shift;
                }
            }
            Surface3::BSpline(bs)
        }
        Surface3::Bezier(ref bz) => {
            let mut bz = bz.clone();
            for row in &mut bz.control_points {
                for cp in row {
                    *cp += shift;
                }
            }
            Surface3::Bezier(bz)
        }
        Surface3::LinearExtrusion(ref le) => {
            let mut le = le.clone();
            le.direction = le.direction; // direction unchanged
            // The profile curve's origin is not directly accessible as a field;
            // clone without position modification for now
            Surface3::LinearExtrusion(le)
        }
        ref other => other.clone(),
    }
}

/// Translate a 3D curve by a displacement vector.
/// All control points and origin/center positions are shifted.
///
/// aFaceFace.ApplyTrsf() -- reversing the shift after intersection
pub(crate) fn translate_curve3(curve: &Curve3, shift: DVec3) -> Curve3 {
    match *curve {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin + shift,
            ..l
        }),
        Curve3::Circle(c) => Curve3::Circle(Circle3 {
            center: c.center + shift,
            ..c
        }),
        Curve3::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            center: e.center + shift,
            ..e
        }),
        Curve3::BSpline(ref bs) => {
            let mut bs = bs.clone();
            for cp in &mut bs.control_points {
                *cp += shift;
            }
            Curve3::BSpline(bs)
        }
        Curve3::Bezier(ref bz) => {
            let mut bz = bz.clone();
            for cp in &mut bz.control_points {
                *cp += shift;
            }
            Curve3::Bezier(bz)
        }
        Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
            center: h.center + shift,
            ..h
        }),
        Curve3::Parabola(p) => Curve3::Parabola(Parabola3 {
            vertex: p.vertex + shift,
            ..p
        }),
        Curve3::Offset(ref o) => {
            let mut o = o.clone();
            o.basis = Box::new(translate_curve3(&o.basis, shift));
            Curve3::Offset(o)
        }
        Curve3::CircularHelix(ref h) => {
            let mut h = h.clone();
            h.origin += shift;
            Curve3::CircularHelix(h)
        }
        Curve3::SineWave(ref sw) => {
            let mut sw = sw.clone();
            sw.origin += shift;
            Curve3::SineWave(sw)
        }
        Curve3::Trimmed(ref tc) => {
            let (inner, first, last) = (tc.curve.clone(), tc.first, tc.last);
            let translated = translate_curve3(&inner, shift);
            Curve3::Trimmed(Box::new(TrimmedCurve3::new(translated, first, last)))
        }
    }
}

// ---- Phase 2a: MakeBlocks candidate injection helpers ----

/// Find up-to-2 face indices that reference a given intersection curve.
/// checks curves_sc (PaveBlocksSc).
pub(crate) fn find_face_idxs_for_curve(ds: &DS, ci: usize) -> [usize; 2] {
    let mut result = [usize::MAX; 2];
    let mut idx = 0;
    for fi in 0..ds.face_count() {
        if ds.face_info(fi).curves_sc.contains(&ci) {
            if idx < 2 {
                result[idx] = fi;
                idx += 1;
            }
        }
    }
    result
}

/// PutPaveOnCurve (BOPAlgo_PaveFiller_6.cxx L833-900)
///    OCCT: EF vertices first (theMVEF), then ON/IN vertices (theMVOnIn) with
///    BBox filtering (aBoxC.IsOut(aBoxV), L2409) and IsNewShape check (L2413-2415).
///    This prevents projecting too many vertices onto each IC, which would cause
///    excessive edge splitting (bfuse_simple B3: 539 edges -> 28 ref).
pub(crate) fn put_pave_on_curve_full(
    ds: &DS,
    curve_idx: usize,
    face_idxs: &[usize; 2],
) -> Vec<(f64, usize)> {
    let ic = &ds.intersection_curves[curve_idx];
    let [t0, t1] = ic.t_range;
    let a_tol_r3d = ic.geom_tol;
    let mut paves: Vec<(f64, usize)> = Vec::new();

    // compute curve bounding box for vertex filtering (L2409: aBoxC.IsOut(aBoxV)).
    let curve_bbox = curve_bounding_box_simple(&ic.curve, a_tol_r3d);

    // GetStickVertices (PaveFiller_6.cxx L2847-2905) collects EF vertex set
    //   per FF pair.  Only EF vertices belonging to this specific pair are added to aMVEF.
    //   rcad: filter EF vertices by checking if the interference's face is in this pair.
    let ef_vertices: std::collections::HashSet<usize> = ds
        .interf_ef
        .iter()
        .filter_map(|inf| {
            //   Both sub-shapes belong to the two faces -- the EF vertex involves this pair.
            if inf.face == face_idxs[0] || inf.face == face_idxs[1] {
                Some(inf.new_vertex)
            } else {
                None
            }
        })
        .collect();

    for &fi in face_idxs.iter().filter(|&&fi| fi != usize::MAX) {
        for &vi in &ds.face_info(fi).vertices_on {
            if !ef_vertices.contains(&vi) {
                continue;
            } // OCCT GetStickVertices: skip non-pair EF
            if vi == ic.start_vertex || vi == ic.end_vertex {
                continue;
            }
            if paves.iter().any(|&(_, v)| v == vi) {
                continue;
            }
            let pt = ds.vertex_point(vi);
            if let Some(t) = project_vertex_to_curve(pt, &ic.curve, a_tol_r3d) {
                if t >= t0 - a_tol_r3d && t <= t1 + a_tol_r3d {
                    paves.push((t, vi));
                }
            }
        }
        for &vi in &ds.face_info(fi).vertices_in {
            if vi == ic.start_vertex || vi == ic.end_vertex {
                continue;
            }
            if ef_vertices.contains(&vi) {
                continue;
            }
            if paves.iter().any(|&(_, v)| v == vi) {
                continue;
            }
            if let Some([c_min, c_max]) = curve_bbox {
                let v_pt = ds.vertex_point(vi);
                let v_tol = ds.vertex_tolerance(vi).max(a_tol_r3d);
                let v_min = v_pt - DVec3::splat(v_tol);
                let v_max = v_pt + DVec3::splat(v_tol);
                if v_max.x < c_min.x
                    || v_min.x > c_max.x
                    || v_max.y < c_min.y
                    || v_min.y > c_max.y
                    || v_max.z < c_min.z
                    || v_min.z > c_max.z
                {
                    continue;
                }
            }
            if !ds.is_new_vertex(vi) {
                continue;
            }

            let pt = ds.vertex_point(vi);
            if let Some(t) = project_vertex_to_curve(pt, &ic.curve, a_tol_r3d) {
                if t >= t0 - a_tol_r3d && t <= t1 + a_tol_r3d {
                    paves.push((t, vi));
                }
            }
        }
    }

    // Sort by parameter
    paves.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    // Deduplicate by parameter or vertex idx
    paves.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-12 || a.1 == b.1);

    // PutClosingPaveOnCurve for closed curves
    put_closing_pave_on_curve(&mut paves, matches!(&ic.curve, Curve3::Circle(_)));

    paves
}

/// ✅ OCCT-aligned: delegates to bnd_lib::curve_bounds (OCCT GeomBndLib per-type dispatch).
pub(crate) fn curve_bounding_box_simple(curve: &Curve3, tol: f64) -> Option<[DVec3; 2]> {
    let bbox = bnd_lib::curve_bounds(curve, tol);
    if bbox.is_valid() {
        Some([bbox.min, bbox.max])
    } else {
        None
    }
}

/// FilterPavesOnCurves (PaveFiller_6.cxx L2437-2538).
/// OCCT uses a multi-candidate distance comparison + sin-angle check.
/// rcad simplified: single-threshold filter against curve tolerance + fuzzy.
/// OCCT L2449: aTolR3D = max(curve.Tolerance(), curve.TangentialTolerance())
pub(crate) fn filter_paves_on_curves(
    ds: &DS,
    curve_idx: usize,
    paves: &[(f64, usize)],
) -> Vec<(f64, usize)> {
    let ic = &ds.intersection_curves[curve_idx];
    // curve tolerance + fuzzy (SUM matching PutPaveOnCurve L2976 aTolR3D + myFuzzyValue)
    let tol = ic.geom_tol + TOLERANCE_ABS;
    let tol_sq = tol * tol;
    paves
        .iter()
        .filter(|&&(_, vi)| {
            let pt = ds.vertex_point(vi);
            let dist_sq = match &ic.curve {
                Curve3::Line(line) => {
                    let to_pt = pt - line.origin;
                    let proj = line.origin + line.direction * to_pt.dot(line.direction);
                    proj.distance_squared(pt)
                }
                Curve3::Circle(circ) => {
                    let center_dist = pt.distance(circ.center);
                    (center_dist - circ.radius).powi(2)
                }
                _ => 0.0,
            };
            dist_sq < tol_sq
        })
        .copied()
        .collect()
}

/// PutClosingPaveOnCurve (L828-833)
///    Only replace the last vertex when the curve spans a full closed period (parameter diff ~ 2*pi or full curve range).
///    Arc segments (parameter diff < pi) are not replaced, to avoid incorrectly changing arc endpoints to start points.
pub(crate) fn put_closing_pave_on_curve(paves: &mut Vec<(f64, usize)>, is_closed: bool) {
    if paves.len() < 2 {
        return;
    }
    if is_closed {
        let first_t = paves[0].0;
        let last_t = paves[paves.len() - 1].0;
        let span = last_t - first_t;
        // Only replace if the curve spans at least one full period (~2*pi for circles)
        if (span - std::f64::consts::TAU).abs() < 0.1 {
            let first_vi = paves[0].1;
            let last_idx = paves.len() - 1;
            paves[last_idx].1 = first_vi;
        }
    }
}

/// Intersect two bounded line segments in 3D. Returns (t1, t2, point) if they
/// cross within tolerance.
pub(crate) fn intersect_line_line(
    l1: &Line3,
    r1: [f64; 2],
    l2: &Line3,
    r2: [f64; 2],
    coincidence_tol: f64,
) -> Option<(f64, f64, DVec3)> {
    let tol = coincidence_tol.max(TOLERANCE_ABS);
    let tol_sq = tol * tol;
    let d1 = l1.direction;
    let d2 = l2.direction;
    let w0 = l1.origin - l2.origin;

    let a = d1.dot(d1);
    let b = d1.dot(d2);
    let c = d2.dot(d2);
    let d = d1.dot(w0);
    let e = d2.dot(w0);

    let denom = a * c - b * b;
    if denom.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // Parallel lines. Check if they are colinear (on the same line).
        // If colinear, compute the overlap of their ranges and return the midpoint.
        let cross_sq = d1.cross(w0).length_squared();
        let d1_sq = d1.length_squared();
        if cross_sq > tol_sq * d1_sq.max(1.0) {
            return None; // parallel but not colinear -- no intersection
        }
        // Colinear: map l2's range into l1's parameter space.
        // l1: P(t) = l1.origin + t * d1
        // l2: P(s) = l2.origin + s * d2
        // For colinear lines, d2 = +/- d1 (parallel). Map s-parameter to t:
        // t = (l2.origin - l1.origin).dot(d1) / d1_sq + s * (d2.dot(d1) / d1_sq)
        let sign = if d1.dot(d2) > 0.0 { 1.0 } else { -1.0 };
        let origin_offset = (l2.origin - l1.origin).dot(d1) / d1_sq;
        let t2_lo = origin_offset + r2[0] * sign;
        let t2_hi = origin_offset + r2[1] * sign;
        let overlap_lo = r1[0].max(t2_lo.min(t2_hi));
        let overlap_hi = r1[1].min(t2_lo.max(t2_hi));
        if overlap_hi <= overlap_lo + tol {
            return None; // no overlap
        }
        let t_mid = (overlap_lo + overlap_hi) * 0.5;
        let s_mid = (t_mid - origin_offset) * sign;
        let p = l1.origin + d1 * t_mid;
        return Some((t_mid, s_mid, p));
    }

    let t1 = (b * e - c * d) / denom;
    let t2 = (a * e - b * d) / denom;

    // Check within ranges
    if t1 < r1[0] - tol || t1 > r1[1] + tol || t2 < r2[0] - tol || t2 > r2[1] + tol {
        return None;
    }

    let p1 = l1.origin + d1 * t1;
    let p2 = l2.origin + d2 * t2;

    if (p1 - p2).length_squared() > tol_sq {
        return None; // skew, don't actually intersect
    }

    Some((t1, t2, (p1 + p2) * 0.5))
}

// ---- Sampling helpers for marching seed-point generation ----

/// Sample a flat plane (infinite) over a 2D square of side `half_extent*2`
/// centred at `plane.origin`.
pub(crate) fn sample_plane(plane: &Plane, half_extent: f64, n: usize) -> Vec<DVec3> {
    let u = plane.u_dir;
    let v = plane.v_dir;
    let mut pts = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let su = -half_extent + 2.0 * half_extent * i as f64 / (n - 1).max(1) as f64;
            let sv = -half_extent + 2.0 * half_extent * j as f64 / (n - 1).max(1) as f64;
            pts.push(plane.origin + u * su + v * sv);
        }
    }
    pts
}

/// Sample a cone surface between heights `h_min` and `h_max` along its axis.
pub(crate) fn sample_cone(
    cone: &ConicalSurface,
    h_min: f64,
    h_max: f64,
    n_theta: usize,
    n_h: usize,
) -> Vec<DVec3> {
    let u = rcad_kernel::any_perpendicular(cone.axis);
    let v = cone.axis.cross(u);
    let tan_h = cone.half_angle_rad.tan();
    let mut pts = Vec::with_capacity(n_theta * n_h);
    for ih in 0..n_h {
        let h = h_min + (h_max - h_min) * ih as f64 / (n_h - 1).max(1) as f64;
        let r = h * tan_h;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = cone.apex + cone.axis * h + (u * theta.cos() + v * theta.sin()) * r;
            pts.push(p);
        }
    }
    pts
}

/// Sample `n` points on a circular arc from `t_start` to `t_end`.
pub(crate) fn sample_circle_arc(
    circle: &Circle3,
    t_start: f64,
    t_end: f64,
    n: usize,
) -> Vec<DVec3> {
    use rcad_kernel::CurveEval;
    use rcad_kernel::geom::Curve3;
    let curve = Curve3::Circle(*circle);
    (0..n)
        .map(|i| {
            let t = t_start + (t_end - t_start) * i as f64 / (n - 1).max(1) as f64;
            curve.point_at(t)
        })
        .collect()
}

/// Compute the angular parameter of `point` on `circle` in [0, 2*pi).
/// uses deterministic reference direction (normal x DVec3::X).
pub(crate) fn circle_param(point: DVec3, circle: &Circle3) -> f64 {
    let nm = circle.normal.normalize();
    let x_ax = if nm.x.abs() < 0.9 {
        nm.cross(DVec3::X).normalize()
    } else {
        nm.cross(DVec3::Y).normalize()
    };
    let y_ax = nm.cross(x_ax);
    let d = point - circle.center;
    let mut theta = d.dot(y_ax).atan2(d.dot(x_ax));
    if theta < 0.0 {
        theta += std::f64::consts::TAU;
    }
    theta
}

/// Intersect a 3D line with a 3D circle.
/// Returns `(t_on_line, t_on_circle, point)` for each intersection found.
pub(crate) fn intersect_line_circle(
    line: &Line3,
    circle: &Circle3,
    tol: f64,
) -> Vec<(f64, f64, DVec3)> {
    let mut results = Vec::new();
    let d = line.direction;
    let o = line.origin;
    let c = circle.center;
    let n = circle.normal;
    let r = circle.radius;
    let r_sq = r * r;

    // Planarity constraint: every point on the circle satisfies (P - c).n = 0.
    let dn = d.dot(n);
    let w = o - c;
    let wn = w.dot(n);

    if dn.abs() > tol {
        // Line pierces the circle plane at one point.
        let t = -wn / dn;
        let p = o + d * t;
        // check distance to circle circumference, not inside-circle.
        // (p - c).length_squared <= r_sq allows points at the circle CENTER (false positive).
        let dist = (p - c).length();
        if (dist - r).abs() <= tol {
            results.push((t, circle_param(p, circle), p));
        }
    } else if wn.abs() <= tol {
        // Line lies in the circle plane -- solve 2D line-circle.
        let t_closest = -w.dot(d);
        let perp_dist_sq = ((o + d * t_closest) - c).length_squared();

        if perp_dist_sq <= r_sq + tol * tol {
            let along = (r_sq - perp_dist_sq).max(0.0).sqrt();
            let t1 = t_closest - along;
            let p1 = o + d * t1;
            results.push((t1, circle_param(p1, circle), p1));

            let t2 = t_closest + along;
            if (t2 - t1).abs() > tol {
                let p2 = o + d * t2;
                results.push((t2, circle_param(p2, circle), p2));
            }
        }
    }

    results
}

/// Intersect two coplanar 3D circles (their planes are parallel/coincident).
pub(crate) fn intersect_coplanar_circles(c1: &Circle3, c2: &Circle3, tol: f64) -> Vec<DVec3> {
    let d_vec = c2.center - c1.center;
    let d = d_vec.length();
    let r1 = c1.radius;
    let r2 = c2.radius;

    // Disjoint or concentric -- no isolated intersection points
    if d > r1 + r2 + tol || d < (r1 - r2).abs() - tol || d < tol {
        return vec![];
    }

    // 2D circle-circle intersection
    // x = projection of intersection point onto the line of centers
    let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let y_sq = r1 * r1 - x * x;

    if y_sq < -tol * tol {
        return vec![];
    }
    let y = y_sq.max(0.0).sqrt();

    let dir = d_vec / d;
    let perp = c1.normal.cross(dir).try_normalize().unwrap_or(DVec3::ZERO);

    let mid = c1.center + dir * x;
    if y < tol || perp == DVec3::ZERO {
        vec![mid]
    } else {
        vec![mid + perp * y, mid - perp * y]
    }
}

/// Intersect two 3D circles that may lie in different planes.
/// Returns up to 2 intersection points.
pub(crate) fn intersect_circle_circle(c1: &Circle3, c2: &Circle3, tol: f64) -> Vec<DVec3> {
    let n1 = c1.normal;
    let n2 = c2.normal;
    let cross = n1.cross(n2);
    let cross_len_sq = cross.length_squared();

    // Parallel/coincident planes -- coplanar circle-circle case
    if cross_len_sq < TOLERANCE_ANG * TOLERANCE_ANG {
        let offset = (c2.center - c1.center).dot(n1).abs();
        if offset > tol {
            return vec![];
        }
        return intersect_coplanar_circles(c1, c2, tol);
    }

    // Planes intersect in a line L along the cross-product direction.
    let line_dir = cross / cross_len_sq.sqrt();
    let b = n1.dot(n2);
    let denom = 1.0 - b * b; // sin^2(theta) > 0 (not parallel)
    let h1 = c1.center.dot(n1);
    let h2 = c2.center.dot(n2);
    let alpha = (h1 - h2 * b) / denom;
    let beta = (h2 - h1 * b) / denom;
    let base = n1 * alpha + n2 * beta; // a point on line L

    // Intersect sphere of circle1 (center=c1.center, radius=r1) with line L.
    let w = base - c1.center;
    let a = line_dir.dot(line_dir); // = 1 for unit direction
    let b2 = 2.0 * w.dot(line_dir);
    let c = w.dot(w) - c1.radius * c1.radius;
    let disc = b2 * b2 - 4.0 * a * c;

    if disc < -tol * tol {
        return vec![];
    }
    if disc < tol * tol {
        let t = -b2 / (2.0 * a);
        let p = base + line_dir * t;
        return if (p - c2.center).length_squared() <= (c2.radius + tol) * (c2.radius + tol) {
            vec![p]
        } else {
            vec![]
        };
    }

    let sqrt_disc = disc.sqrt();
    let t1 = (-b2 - sqrt_disc) / (2.0 * a);
    let t2 = (-b2 + sqrt_disc) / (2.0 * a);
    let p1 = base + line_dir * t1;
    let p2 = base + line_dir * t2;

    let r2_tol_sq = (c2.radius + tol) * (c2.radius + tol);
    let mut results = Vec::with_capacity(2);
    if (p1 - c2.center).length_squared() <= r2_tol_sq {
        results.push(p1);
    }
    if (p2 - p1).length_squared() > tol * tol && (p2 - c2.center).length_squared() <= r2_tol_sq {
        results.push(p2);
    }
    results
}

/// Check if a parameter `t` falls within `range` (inclusive, with tolerance).
pub(crate) fn in_range(t: f64, range: [f64; 2], tol: f64) -> bool {
    let lo = range[0].min(range[1]) - tol;
    let hi = range[0].max(range[1]) + tol;
    t >= lo && t <= hi
}

pub(crate) fn point_in_sphere_face(pt: DVec3, boundary_verts: &[DVec3], _ds: &DS) -> bool {
    if boundary_verts.is_empty() {
        return false;
    }
    // OCCT-style single-seam sphere: only two pole vertices. An axis-aligned hull of those
    // poles rejects almost every real point on the sphere (e.g. equator vs poles on +/-Y),
    // so plane-sphere tangent handling never records `FaceFace` points and downstream
    // trimming misses imprint geometry (see OCCT `bcommon_simple/A4`).
    if boundary_verts.len() == 2 {
        let a = boundary_verts[0];
        let b = boundary_verts[1];
        let diam = (a - b).length();
        let r = diam * 0.5;
        if r < TOLERANCE_LEN_MIN {
            return false;
        }
        let c = (a + b) * 0.5;
        let radial_err = ((pt - c).length() - r).abs();
        return radial_err < (TOLERANCE_ABS * 500.0).max(TOLERANCE_COORD_SUB * r);
    }
    // Convex hull approximation for faces with a full boundary polygon.
    let cx = boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::max) + TOLERANCE_COORD_SUB);
    let cy = boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::max) + TOLERANCE_COORD_SUB);
    let cz = boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::max) + TOLERANCE_COORD_SUB);
    cx.contains(&pt.x) && cy.contains(&pt.y) && cz.contains(&pt.z)
}

/// Generic UV-grid sampling for any surface type via `SurfaceEval::default_domain()`.
/// Works for BSpline, Bezier, Offset, Revolution, Trimmed, LinearExtrusion.
pub(crate) fn sample_surface_generic(surface: &Surface3, n_u: usize, n_v: usize) -> Vec<DVec3> {
    use rcad_kernel::geom::SurfaceEval;
    let [u0, u1, v0, v1] = surface.default_domain();
    let mut pts = Vec::with_capacity(n_u * n_v);
    for iu in 0..n_u {
        for iv in 0..n_v {
            let u = u0 + (u1 - u0) * iu as f64 / (n_u - 1).max(1) as f64;
            let v = v0 + (v1 - v0) * iv as f64 / (n_v - 1).max(1) as f64;
            let p = surface.point_at(u, v);
            if p.is_finite() {
                pts.push(p);
            }
        }
    }
    pts
}

/// Numeric edge-face intersection: sample the curve, find sign changes of the
/// surface implicit function, then refine via bisection.
///
/// Used as fallback for unsupported curve x surface combinations (Ellipse,
/// Hyperbola, Parabola, BSpline, Bezier, OffsetCurve -- any surface).
pub(crate) fn intersect_edge_face_numeric(
    curve: &Curve3,
    surface: &Surface3,
    t_range: [f64; 2],
    geom_tol: f64,
) -> Vec<(DVec3, f64)> {
    use rcad_kernel::CurveEval;
    use rcad_kernel::geom::SurfaceEval;
    use rcad_kernel::projection::closest_point_on_surface;
    const N_SAMPLES: usize = 64;
    const MAX_BISECT: usize = 30;

    let eps = geom_tol.max(TOLERANCE_ABS);
    let zero_tol = (eps * TOLERANCE_AREA_REL).max(TOLERANCE_LEN_MIN);

    let [t0, t1] = t_range;
    let mut values = Vec::with_capacity(N_SAMPLES + 1);
    let mut points = Vec::with_capacity(N_SAMPLES + 1);

    for i in 0..=N_SAMPLES {
        let t = t0 + (t1 - t0) * i as f64 / N_SAMPLES as f64;
        let p = curve.point_at(t);
        if !p.is_finite() {
            values.push(f64::NAN);
            points.push(p);
            continue;
        }
        values.push(inttools::marching::surface_implicit(surface, p));
        points.push(p);
    }

    let mut hits = Vec::new();
    for i in 0..N_SAMPLES {
        let va = values[i];
        let vb = values[i + 1];
        if va.is_nan() || vb.is_nan() {
            continue;
        }
        if va * vb > 0.0 {
            continue;
        }
        // Bisection refinement (Stage 1 -- coarse detection)
        let mut ta = t0 + (t1 - t0) * i as f64 / N_SAMPLES as f64;
        let mut tb = t0 + (t1 - t0) * (i + 1) as f64 / N_SAMPLES as f64;
        let mut fa = va;
        let mut fb = vb;
        for _ in 0..MAX_BISECT {
            let tm = (ta + tb) * 0.5;
            let pm = curve.point_at(tm);
            if !pm.is_finite() {
                break;
            }
            let fm = inttools::marching::surface_implicit(surface, pm);
            if fm.abs() < zero_tol {
                hits.push((pm, tm));
                break;
            }
            if (tb - ta).abs() < zero_tol {
                hits.push((pm, tm));
                break;
            }
            if fa * fm < 0.0 {
                tb = tm;
                fb = fm;
            } else {
                ta = tm;
                fa = fm;
            }
        }
        // If bisection didn't converge well, use midpoint
        let tm = (ta + tb) * 0.5;
        let pm = curve.point_at(tm);
        let dedup_dt = (TOLERANCE_MESH_LEGACY).max(eps * 10.0);
        if pm.is_finite() && !hits.iter().any(|(_, t)| (t - tm).abs() < dedup_dt) {
            hits.push((pm, tm));
        }
    }

    // Stage 2 -- Newton refinement: polish each bisection result
    // IntCurveSurface_TheExactHInter two-stage approach
    //   coarse sign-change detection -> Newton-Raphson refinement.
    for (point, t) in hits.iter_mut() {
        let initial_t = *t;
        let initial_point = *point;

        // Get initial UV guess via closest-point projection
        let proj = closest_point_on_surface(surface, initial_point, 8);
        let initial_uv = DVec2::new(proj.params.0, proj.params.1);

        if let Some((refined_t, refined_uv)) = inttools::curve_surface::newton_refine_curve_surface(
            curve, initial_t, surface, initial_uv, 20, eps,
        ) {
            // validation (Stage 3):
            //   1. t within the curve's parametric range
            if refined_t < t_range[0] - eps || refined_t > t_range[1] + eps {
                continue; // Keep bisection result
            }

            //   2. uv within the surface's natural UV domain (if bounded)
            let [u0, u1, v0, v1] = surface.default_domain();
            let u_ok = u0.is_infinite()
                || u1.is_infinite()
                || (refined_uv.x >= u0 - eps && refined_uv.x <= u1 + eps);
            let v_ok = v0.is_infinite()
                || v1.is_infinite()
                || (refined_uv.y >= v0 - eps && refined_uv.y <= v1 + eps);
            if !u_ok || !v_ok {
                continue; // Keep bisection result
            }

            //   3. Distance |C(t) - S(uv)| within tolerance
            let refined_point = curve.point_at(refined_t);
            let surface_point = surface.point_at(refined_uv.x, refined_uv.y);
            if (refined_point - surface_point).length() > eps * 10.0 {
                continue; // Keep bisection result
            }

            // Newton refinement passed all checks -- replace the hit
            *t = refined_t;
            *point = refined_point;
        }
    }

    hits
}

/// Result of partial face overlap analysis.
#[derive(Debug, Clone)]
pub(crate) struct PartialOverlapInfo {
    /// Face index in shape A.
    pub face_a: usize,
    /// Face index in shape B.
    pub face_b: usize,
    /// Estimated overlap ratio (0.0 to 1.0).
    pub overlap_ratio: f64,
    /// Overlap type.
    pub overlap_type: PartialOverlapType,
}

/// Type of partial overlap between faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PartialOverlapType {
    /// Faces are coplanar with partial boundary overlap.
    CoplanarBoundary,
    /// Faces share an edge partially.
    EdgeOverlap,
    /// One face is contained within another.
    Contained,
}

/// Result of edge overlap detection between two edges.
#[derive(Debug, Clone)]
pub(crate) struct EdgeOverlapResult {
    /// Edge index in shape A.
    pub edge_a: usize,
    /// Edge index in shape B.
    pub edge_b: usize,
    /// Type of overlap detected.
    pub overlap_type: EdgeOverlapType,
    /// Overlap ratio for the first edge (0.0 to 1.0).
    pub overlap_ratio_a: f64,
    /// Overlap ratio for the second edge (0.0 to 1.0).
    pub overlap_ratio_b: f64,
    /// Parameter range of overlap on edge A [t_start, t_end].
    pub param_range_a: Option<[f64; 2]>,
    /// Parameter range of overlap on edge B [t_start, t_end].
    pub param_range_b: Option<[f64; 2]>,
    /// Maximum distance between edges in the overlap region.
    pub max_distance: f64,
}

/// Type of overlap between two edges.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EdgeOverlapType {
    /// No overlap - edges are on different curves or don't intersect.
    None,
    /// Partial overlap - edges share part of their parameter range.
    Partial,
    /// Full overlap - one edge completely overlaps the other.
    Full,
    /// Edge A is contained within edge B's parameter range.
    AContainedInB,
    /// Edge B is contained within edge A's parameter range.
    BContainedInA,
}

/// Result of edge containment detection.
#[derive(Debug, Clone)]
pub(crate) struct EdgeContainmentResult {
    /// Edge index that is contained.
    pub contained_edge: usize,
    /// Edge index that contains.
    pub containing_edge: usize,
    /// Containment ratio (how much of the contained edge is inside).
    pub containment_ratio: f64,
    /// Whether the containment is exact within tolerance.
    pub is_exact: bool,
}

/// Parameter overlap result for two parameter ranges.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ParamOverlap {
    /// Overlap type.
    pub overlap_type: ParamOverlapType,
    /// Overlap range [min, max] if any overlap exists.
    pub overlap_range: Option<[f64; 2]>,
    /// Ratio of first range that overlaps.
    pub ratio_a: f64,
    /// Ratio of second range that overlaps.
    pub ratio_b: f64,
}

/// Type of parameter range overlap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ParamOverlapType {
    /// No overlap.
    None,
    /// Partial overlap - ranges partially intersect.
    Partial,
    /// Range A contains range B entirely.
    AContainsB,
    /// Range B contains range A entirely.
    BContainsA,
    /// Exact match - ranges are identical.
    Exact,
}

/// Result of near-tangent face detection.
#[derive(Debug, Clone)]
pub(crate) struct NearTangentFaceInfo {
    /// Face index in shape A.
    pub face_a: usize,
    /// Face index in shape B.
    pub face_b: usize,
    /// Distance between faces at closest point.
    pub distance: f64,
    /// Type of tangency.
    pub tangent_type: NearTangentType,
    /// Whether the faces should be merged.
    pub should_merge: bool,
}

/// Result of near-coincident face detection.
#[derive(Debug, Clone)]
pub(crate) struct NearCoincidentFaceInfo {
    /// Face index in shape A.
    pub face_a: usize,
    /// Face index in shape B.
    pub face_b: usize,
    /// Maximum distance between faces in overlap region.
    pub max_distance: f64,
    /// Area of overlap region (approximate).
    pub overlap_area: f64,
    /// Whether faces should be merged.
    pub should_merge: bool,
}

/// Result of micro-gap detection.
#[derive(Debug, Clone)]
pub(crate) struct MicroGapInfo {
    /// Edge index on shape A.
    pub edge_a: usize,
    /// Edge index on shape B.
    pub edge_b: usize,
    /// Gap distance.
    pub gap_distance: f64,
    /// Whether the gap can be bridged.
    pub can_bridge: bool,
}

/// Result of coincident edge detection.
#[derive(Debug, Clone)]
pub(crate) struct CoincidentEdgeInfo {
    /// Edge index in shape A.
    pub edge_a: usize,
    /// Edge index in shape B.
    pub edge_b: usize,
    /// Maximum distance between edges.
    pub max_distance: f64,
    /// Overlap ratio (0.0 to 1.0).
    pub overlap_ratio: f64,
    /// Whether edges should be merged.
    pub should_merge: bool,
}
