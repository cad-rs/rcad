use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};

/// ✅ OCCT-aligned: IntTools_EdgeEdge::TypeToInteger (cxx L1456-1482).
/// Maps curve type to integer priority for edge swapping (lower = simpler).
pub fn curve_type_to_integer(curve: &Curve3) -> i32 {
    match curve {
        Curve3::Line(_) => 0,
        Curve3::Hyperbola(_) | Curve3::Parabola(_) => 1,
        Curve3::Circle(_) | Curve3::Ellipse(_) => 2,
        Curve3::BSpline(_) | Curve3::Bezier(_) => 3,
        _ => 4,
    }
}

/// ✅ OCCT-aligned: IntTools_EdgeEdge::PointBoxDistance (cxx L1423-1452).
/// Computes min distance from a point to an axis-aligned bounding box.
pub fn point_box_distance(p: DVec3, box_min: DVec3, box_max: DVec3) -> f64 {
    let mut dist = 0.0;
    for i in 0..3 {
        let c = [p.x, p.y, p.z][i];
        let bmin = [box_min.x, box_min.y, box_min.z][i];
        let bmax = [box_max.x, box_max.y, box_max.z][i];
        if c < bmin {
            let d = bmin - c;
            dist += d * d;
        } else if c > bmax {
            let d = c - bmax;
            dist += d * d;
        }
    }
    dist.sqrt()
}

/// ✅ OCCT-aligned: IntTools_EdgeEdge::SplitRangeOnSegments (cxx L1366-1406).
/// Splits range [aT1, aT2] into segments based on resolution. Returns number of segments.
pub fn split_range_on_segments(t1: f64, t2: f64, resolution: f64, nb_seg: i32) -> (i32, Vec<[f64; 2]>) {
    let diff = t2 - t1;
    if diff < resolution || nb_seg == 1 {
        return (1, vec![[t1, t2]]);
    }
    let mut a_nb_segments = nb_seg;
    let mut a_dt = diff / a_nb_segments as f64;
    if a_dt < resolution {
        let seg = (diff / resolution) as i32;
        a_nb_segments = seg + 1;
        a_dt = diff / a_nb_segments as f64;
    }
    let mut segments = Vec::new();
    let mut t1x = t1;
    for _ in 1..a_nb_segments {
        let t2x = t1x + a_dt;
        segments.push([t1x, t2x]);
        t1x = t2x;
    }
    segments.push([t1x, t2]);
    (a_nb_segments, segments)
}

/// ✅ OCCT-aligned: IntTools_EdgeEdge::Resolution (cxx L1561-1607).
/// Computes curve resolution (parameter step for a given 3D tolerance).
/// For lines: returns theR3D directly.
/// For circles: 2*asin(res_coeff * theR3D).
/// For BSpline/Bezier: delegates to curve resolution method.
pub fn curve_resolution_edge(curve: &Curve3, res_coeff: f64, r3d: f64) -> f64 {
    match curve {
        Curve3::Line(_) => r3d,
        Curve3::Circle(c) => {
            let dt = res_coeff * r3d;
            if dt <= 1.0 { 2.0 * dt.asin() } else { std::f64::consts::TAU }
        }
        Curve3::Ellipse(e) => {
            let dt = res_coeff * r3d;
            if dt <= 1.0 { 2.0 * dt.asin() } else { std::f64::consts::TAU }
        }
        Curve3::BSpline(_) | Curve3::Bezier(_) => {
            // OCCT: theCurve->Resolution(theR3D, aRes)
            // rcad: approximate using curve resolution
            crate::inttools::curve_range::curve_resolution(curve, 0.0, r3d)
        }
        _ => res_coeff * r3d,
    }
}

/// ✅ OCCT-aligned: IntTools_EdgeEdge::ResolutionCoeff (cxx L1486-1557).
/// Computes the resolution coefficient for a curve type.
/// For circles: 1/(2*radius). For ellipses: 1/major_radius.
pub fn resolution_coeff(curve: &Curve3, t_range: [f64; 2]) -> f64 {
    match curve {
        Curve3::Circle(c) => 1.0 / (2.0 * c.radius.max(1e-30)),
        Curve3::Ellipse(e) => 1.0 / e.major_radius.max(1e-30),
        _ => {
            // OCCT: sample 30 points, find min dt/dist ratio
            let nb_p = 30;
            let t1 = t_range[0];
            let t2 = t_range[1];
            let dt = (t2 - t1) / nb_p as f64;
            let mut t = t1;
            let mut p1 = curve.point_at(t1);
            let mut k_min = 10.0;
            for _ in 1..=nb_p {
                t += dt;
                let p2 = curve.point_at(t);
                let dist = (p1 - p2).length();
                if dist > 1e-30 {
                    let k = dt / dist;
                    if k < k_min { k_min = k; }
                }
                p1 = p2;
            }
            k_min
        }
    }
}

/// ✅ OCCT-aligned: IntTools_EdgeEdge::CurveDeflection (cxx L1611-1638).
/// Computes total angular deflection of a curve over its range by sampling.
pub fn curve_deflection(curve: &Curve3, t_range: [f64; 2]) -> f64 {
    let nb_p = 10;
    let t1 = t_range[0];
    let t2 = t_range[1];
    let dt = (t2 - t1) / nb_p as f64;
    let mut t = t1;
    let mut v1 = curve.tangent_at(t1);
    let mut defl = 0.0;
    for _ in 1..=nb_p {
        t += dt;
        let v2 = curve.tangent_at(t);
        let len1 = v1.length_squared();
        let len2 = v2.length_squared();
        if len1 > 1e-30 && len2 > 1e-30 {
            let d1 = v1 / len1.sqrt();
            let d2 = v2 / len2.sqrt();
            defl += d1.dot(d2).acos();
        }
        v1 = v2;
    }
    defl
}

/// ✅ OCCT-aligned: IsClosed (IntTools_EdgeEdge.cxx L1642-1659).
/// Checks if the curve segment between aT1 and aT2 is closed.
pub fn is_curve_segment_closed(curve: &Curve3, t1: f64, t2: f64, tol: f64, res: f64) -> bool {
    if (t1 - t2).abs() < res { return false; }
    let p1 = curve.point_at(t1);
    let p2 = curve.point_at(t2);
    (p1 - p2).length() < tol
}

/// ✅ OCCT-aligned: ComputeLineLine common part detection (cxx L902-1056).
/// Determines if two line segments intersect, returning parameters if they do.
pub fn intersect_line_line_3d(
    l1_origin: DVec3, l1_dir: DVec3, t1_range: [f64; 2],
    l2_origin: DVec3, l2_dir: DVec3, t2_range: [f64; 2],
    tol: f64,
) -> Option<([f64; 2], [f64; 2], bool)> {
    let d1 = l1_dir.normalize();
    let d2 = l2_dir.normalize();
    let angle = d1.dot(d2).acos();
    let ang_tol = 1e-12;
    let is_coincide = angle < ang_tol || (std::f64::consts::PI - angle).abs() < ang_tol;

    if is_coincide {
        // OCCT L916-919: check distance between lines
        let dist = (l2_origin - l1_origin).cross(d1).length();
        if dist > tol { return None; }
        // Project both endpoints onto line1
        let t21 = (l2_origin - l1_origin).dot(d1);
        let t22 = t21 + (l2_origin + d2 - l1_origin).dot(d1);
        let (mut t21, mut t22) = if t21 < t22 { (t21, t22) } else { (t22, t21) };
        let [t11, t12] = t1_range;
        if (t21 > t12 && t22 > t12) || (t21 < t11 && t22 < t11) { return None; }
        t21 = t21.max(t11);
        t22 = t22.min(t12);
        let range1 = [t11, t12];
        let range2 = [t21, t22];
        return Some((range1, range2, true));
    }

    // Non-coincident lines: find intersection point
    let cross = d1.cross(d2);
    let cross_len2 = cross.length_squared();
    if cross_len2 < 1e-30 { return None; }
    let o1o2 = l2_origin - l1_origin;
    let dist_ll = o1o2.dot(cross / cross_len2.sqrt()).abs();
    if dist_ll > tol { return None; }

    // Find parameters of closest approach
    let a = d1.dot(d2);
    let b = d1.dot(o1o2);
    let c = d2.dot(o1o2);
    let denom = 1.0 - a * a;
    if denom.abs() < 1e-15 { return None; }
    let t1 = (b - a * c) / denom;
    let t2 = (a * b - c) / denom;

    let [t11, t12] = t1_range;
    let [t21, t22] = t2_range;
    if t1 < t11 || t1 > t12 || t2 < t21 || t2 > t22 { return None; }

    let p1 = l1_origin + d1 * t1;
    let p2 = l2_origin + d2 * t2;
    if (p1 - p2).length_squared() > tol * tol { return None; }

    // OCCT L1047-1055: compute intersection range with ComputeIntRange
    let a_tol_1 = tol;
    let a_tol_2 = tol;
    let a_dt1 = crate::boptools::compute_int_range(a_tol_1, a_tol_2, angle);
    let a_dt2 = crate::boptools::compute_int_range(a_tol_2, a_tol_1, angle);

    let range1 = [t1 - a_dt1, t1 + a_dt1];
    let range2 = [t2 - a_dt2, t2 + a_dt2];
    Some((range1, range2, false))
}
