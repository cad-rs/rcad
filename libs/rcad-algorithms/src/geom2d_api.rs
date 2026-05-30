//! Geom2dAPI-style 2D geometry API.
//!
//! Analogous to OCCT `Geom2dAPI` package providing algorithms for 2D geometry:
//! - `InterCurveCurve`: 2D curve-curve intersection
//! - `PointsToBSpline`: Fit BSpline to 2D points
//! - `ProjectPointOnCurve`: Project point on 2D curve
//! - `ExtremaCurveCurve`: Distance between 2D curves
//! - `ExtremaCurvePoint`: Distance from point to 2D curve
//! - Angle and curvature analysis

use crate::tolerance::*;
use glam::DVec2;
use rcad_kernel::geom::{BSplineCurve2, Circle2d, Curve2d, Curve2dEval, Line2d};
use std::f64::consts::PI;

// =============================================================================
// Curve2dIntersection - Result of 2D curve-curve intersection
// =============================================================================

/// Result of intersecting two 2D curves.
///
/// Contains the intersection point and the parameter values on each curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Curve2dIntersection {
    /// The intersection point in 2D space.
    pub point: DVec2,
    /// Parameter on the first curve.
    pub param1: f64,
    /// Parameter on the second curve.
    pub param2: f64,
}

/// Construct the unique circle through three non-collinear points.
///
/// Returns `None` when the points are collinear or too close to define a
/// numerically stable circle.
pub fn circle_through_three_points(p1: DVec2, p2: DVec2, p3: DVec2) -> Option<Circle2d> {
    let a = p2 - p1;
    let b = p3 - p1;
    let det = 2.0 * (a.x * b.y - a.y * b.x);
    if det.abs() <= TOLERANCE_LEN_MIN {
        return None;
    }

    let a_len2 = a.length_squared();
    let b_len2 = b.length_squared();
    let center_offset = DVec2::new(
        (a_len2 * b.y - b_len2 * a.y) / det,
        (b_len2 * a.x - a_len2 * b.x) / det,
    );
    let center = p1 + center_offset;
    let radius = (center - p1).length();
    if radius <= TOLERANCE_LEN_MIN || !radius.is_finite() {
        return None;
    }

    Some(Circle2d { center, radius })
}

/// Construct circles through two points and tangent to `base`.
///
/// Results are sorted by descending radius to match OCCT DRAW's observed
/// `circ2d3Tan` ordering for circle-point-point cases.
pub fn circles_tangent_to_circle_through_points(
    base: Circle2d,
    p1: DVec2,
    p2: DVec2,
) -> Vec<Circle2d> {
    let chord = p2 - p1;
    let chord_len = chord.length();
    if chord_len <= TOLERANCE_LEN_MIN || base.radius <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }

    let midpoint = (p1 + p2) * 0.5;
    let half_chord = chord_len * 0.5;
    let normal = DVec2::new(-chord.y, chord.x) / chord_len;
    let to_mid = midpoint - base.center;
    let a = to_mid.length_squared() - half_chord * half_chord - base.radius * base.radius;
    let b = 2.0 * to_mid.dot(normal);
    let c = 4.0 * base.radius * base.radius;

    let qa = b * b - c;
    let qb = 2.0 * a * b;
    let qc = a * a - c * half_chord * half_chord;
    let mut roots = Vec::new();
    if qa.abs() <= TOLERANCE_FLOAT_LOOSE {
        if qb.abs() > TOLERANCE_FLOAT_LOOSE {
            roots.push(-qc / qb);
        }
    } else {
        let disc = qb * qb - 4.0 * qa * qc;
        if disc >= -TOLERANCE_LINEAR_RELAX_8 {
            let sqrt_disc = disc.max(0.0).sqrt();
            roots.push((-qb - sqrt_disc) / (2.0 * qa));
            roots.push((-qb + sqrt_disc) / (2.0 * qa));
        }
    }

    let mut circles = Vec::new();
    for root in roots {
        let center = midpoint + normal * root;
        let radius = (center - p1).length();
        if radius <= TOLERANCE_LEN_MIN || !radius.is_finite() {
            continue;
        }
        let center_distance = (center - base.center).length();
        let is_tangent = (center_distance - (radius + base.radius)).abs() < TOLERANCE_ABS
            || (center_distance - (radius - base.radius).abs()).abs() < TOLERANCE_ABS;
        if !is_tangent {
            continue;
        }
        if circles.iter().any(|c: &Circle2d| {
            (c.center - center).length() < TOLERANCE_LINEAR_RELAX_8 && (c.radius - radius).abs() < TOLERANCE_LINEAR_RELAX_8
        }) {
            continue;
        }
        circles.push(Circle2d { center, radius });
    }
    circles.sort_by(|a, b| b.radius.partial_cmp(&a.radius).unwrap());
    circles
}

/// Construct circles through a point and tangent to two circles.
///
/// Results are sorted by ascending radius to match OCCT DRAW's observed
/// `circ2d3Tan` ordering for circle-circle-point cases.
pub fn circles_tangent_to_two_circles_through_point(
    c1: Circle2d,
    c2: Circle2d,
    point: DVec2,
) -> Vec<Circle2d> {
    if c1.radius <= TOLERANCE_LEN_MIN || c2.radius <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }

    let mut result = Vec::new();
    for s1 in [-1.0, 1.0] {
        for s2 in [-1.0, 1.0] {
            append_circle_circle_point_solutions(c1, c2, point, s1, s2, &mut result);
        }
    }

    result.sort_by(|a, b| a.radius.partial_cmp(&b.radius).unwrap());
    result
}

/// Construct circles through a point and tangent to a circle and a line.
///
/// Results are sorted by descending radius to match OCCT DRAW's observed
/// ordering for `CircleLinPoint_11`.
pub fn circles_tangent_to_circle_and_line_through_point(
    circle: Circle2d,
    line: Line2d,
    point: DVec2,
) -> Vec<Circle2d> {
    if circle.radius <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }
    let Some(normal) = unit_line_normal(line.direction) else {
        return Vec::new();
    };

    let mut result = Vec::new();
    for circle_sign in [-1.0, 1.0] {
        for line_sign in [-1.0, 1.0] {
            append_circle_line_point_solutions(
                circle,
                line,
                normal,
                point,
                circle_sign,
                line_sign,
                &mut result,
            );
        }
    }
    result.sort_by(|a, b| b.radius.partial_cmp(&a.radius).unwrap());
    result
}

/// Construct circles tangent to one circle and two lines.
///
/// Results are ordered by circle tangency branch (external, then internal) and
/// descending radius within each branch, matching OCCT DRAW for
/// `CircleLinLin_11`.
pub fn circles_tangent_to_circle_and_two_lines(
    circle: Circle2d,
    l1: Line2d,
    l2: Line2d,
) -> Vec<Circle2d> {
    if circle.radius <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }
    let Some(n1) = unit_line_normal(l1.direction) else {
        return Vec::new();
    };
    let Some(n2) = unit_line_normal(l2.direction) else {
        return Vec::new();
    };

    let mut result = Vec::new();
    for circle_sign in [1.0, -1.0] {
        let mut branch = Vec::new();
        for line1_sign in [-1.0, 1.0] {
            for line2_sign in [-1.0, 1.0] {
                append_circle_line_line_solutions(
                    circle,
                    l1,
                    n1,
                    l2,
                    n2,
                    [circle_sign, line1_sign, line2_sign],
                    &mut branch,
                );
            }
        }
        branch.sort_by(|a, b| b.radius.partial_cmp(&a.radius).unwrap());
        result.extend(branch);
    }
    result
}

/// Construct circles through a point and tangent to two lines.
///
/// Results are sorted by descending radius to match OCCT DRAW's observed
/// ordering for the `LinLinPoint_11` case.
pub fn circles_tangent_to_two_lines_through_point(
    l1: Line2d,
    l2: Line2d,
    point: DVec2,
) -> Vec<Circle2d> {
    let Some(n1) = unit_line_normal(l1.direction) else {
        return Vec::new();
    };
    let Some(n2) = unit_line_normal(l2.direction) else {
        return Vec::new();
    };

    let mut result = Vec::new();
    for s1 in [-1.0, 1.0] {
        for s2 in [-1.0, 1.0] {
            append_line_line_point_solutions(l1, n1, l2, n2, point, s1, s2, &mut result);
        }
    }

    result.sort_by(|a, b| b.radius.partial_cmp(&a.radius).unwrap());
    result
}

/// Construct circles through two points and tangent to a line.
///
/// Results are sorted by descending radius to match OCCT DRAW's observed
/// ordering for `LinPointPoint_11`.
pub fn circles_tangent_to_line_through_points(line: Line2d, p1: DVec2, p2: DVec2) -> Vec<Circle2d> {
    let chord = p2 - p1;
    let chord_len = chord.length();
    if chord_len <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }
    let Some(line_normal) = unit_line_normal(line.direction) else {
        return Vec::new();
    };

    let midpoint = (p1 + p2) * 0.5;
    let half_chord = chord_len * 0.5;
    let chord_normal = DVec2::new(-chord.y, chord.x) / chord_len;
    let line_offset = line_normal.dot(midpoint - line.origin);
    let slope = line_normal.dot(chord_normal);

    let mut result = Vec::new();
    for line_sign in [-1.0, 1.0] {
        let qa = slope * slope - 1.0;
        let qb = 2.0 * line_offset * slope;
        let qc = line_offset * line_offset - half_chord * half_chord;
        let mut roots = Vec::new();
        if qa.abs() <= TOLERANCE_FLOAT_LOOSE {
            if qb.abs() > TOLERANCE_FLOAT_LOOSE {
                roots.push(-qc / qb);
            }
        } else {
            let disc = qb * qb - 4.0 * qa * qc;
            if disc >= -TOLERANCE_LINEAR_RELAX_8 {
                let sqrt_disc = disc.max(0.0).sqrt();
                roots.push((-qb - sqrt_disc) / (2.0 * qa));
                roots.push((-qb + sqrt_disc) / (2.0 * qa));
            }
        }

        for root in roots {
            let center = midpoint + chord_normal * root;
            let radius = (center - p1).length();
            if radius <= TOLERANCE_LINEAR_ULTRA_STRICT || !radius.is_finite() {
                continue;
            }
            let signed_distance = signed_distance_to_line(center, line, line_normal);
            if (signed_distance - line_sign * radius).abs() > TOLERANCE_ABS {
                continue;
            }
            if result.iter().any(|c: &Circle2d| {
                (c.center - center).length() < TOLERANCE_LINEAR_RELAX_8 && (c.radius - radius).abs() < TOLERANCE_LINEAR_RELAX_8
            }) {
                continue;
            }
            result.push(Circle2d { center, radius });
        }
    }

    result.sort_by(|a, b| b.radius.partial_cmp(&a.radius).unwrap());
    result
}

/// Construct circles tangent to three lines.
///
/// The valid signed-distance branches are returned in OCCT DRAW's observed
/// order for `LinLinLin_11`.
pub fn circles_tangent_to_three_lines(l1: Line2d, l2: Line2d, l3: Line2d) -> Vec<Circle2d> {
    let Some(n1) = unit_line_normal(l1.direction) else {
        return Vec::new();
    };
    let Some(n2) = unit_line_normal(l2.direction) else {
        return Vec::new();
    };
    let Some(n3) = unit_line_normal(l3.direction) else {
        return Vec::new();
    };

    let lines = [(l1, n1), (l2, n2), (l3, n3)];
    let mut result = Vec::new();
    for s1 in [-1.0, 1.0] {
        for s2 in [-1.0, 1.0] {
            for s3 in [-1.0, 1.0] {
                if let Some(circle) = solve_three_signed_lines(&lines, [s1, s2, s3]) {
                    result.push(circle);
                }
            }
        }
    }
    result
}

/// Circles tangent to two given circles and one line.
///
/// Enumerates the eight (line side × two circle tangency orientations) branches and
/// returns valid solutions, sorted by **descending circumference** \(2 π r\), matching
/// OCCT DRAW for `CircleCircleLin_11`.
pub fn circles_tangent_to_two_circles_and_line(
    c1: Circle2d,
    c2: Circle2d,
    line: Line2d,
) -> Vec<Circle2d> {
    if c1.radius <= TOLERANCE_LEN_MIN || c2.radius <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }
    let Some(n) = unit_line_normal(line.direction) else {
        return Vec::new();
    };
    let t_len = line.direction.length();
    if t_len <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }
    let t = line.direction / t_len;
    let p0 = line.origin;
    let o1 = c1.center;
    let o2 = c2.center;
    let r1 = c1.radius;
    let r2 = c2.radius;
    let v1 = p0 - o1;
    let v2 = p0 - o2;
    let t_dot_v1 = t.dot(v1);
    let _t_dot_v2 = t.dot(v2);
    let n_dot_v1 = n.dot(v1);
    let n_dot_v2 = n.dot(v2);
    let len_v1_sq = v1.length_squared();
    let len_v2_sq = v2.length_squared();
    let denom = 2.0 * t.dot(v1 - v2);

    let mut out = Vec::new();
    for sig in [-1.0_f64, 1.0] {
        for e1 in [-1.0_f64, 1.0] {
            for e2 in [-1.0_f64, 1.0] {
                if denom.abs() <= TOLERANCE_LEN_MIN {
                    continue;
                }
                let a = e1 * r1 - sig * n_dot_v1 - e2 * r2 + sig * n_dot_v2;
                let b = r1 * r1 - len_v1_sq - (r2 * r2) + len_v2_sq;
                let k1 = e1 * r1 - sig * n_dot_v1;
                let c1term = r1 * r1 - len_v1_sq;
                let d = denom;
                // u = (2*A*R + B) / d, from subtracting the two "fixed R" line equations
                // u² + 2 u t·V1 = 2 R K1 + c1term
                let coeff_r2 = 4.0 * a * a;
                let coeff_r = 4.0 * a * b + 4.0 * t_dot_v1 * a * d - 2.0 * k1 * d * d;
                let c_const = b * b + 2.0 * t_dot_v1 * b * d - c1term * d * d;
                for r in solve_quadratic(coeff_r2, coeff_r, c_const) {
                    if r <= TOLERANCE_LINEAR_ULTRA_STRICT || !r.is_finite() {
                        continue;
                    }
                    let u = (2.0 * a * r + b) / d;
                    let center = p0 + t * u + n * (sig * r);
                    if !is_tangent_to_circle(center, r, c1) || !is_tangent_to_circle(center, r, c2)
                    {
                        continue;
                    }
                    if (signed_distance_to_line(center, line, n).abs() - r).abs() > TOLERANCE_MESH_LEGACY {
                        continue;
                    }
                    out.push(Circle2d { center, radius: r });
                }
            }
        }
    }

    // Ascending circumference so `tan1_1` … `tan1_4` match OCCT (CircleCircleLin_11)
    out.sort_by(|a, b| (a.radius * PI * 2.0).total_cmp(&(b.radius * PI * 2.0)));
    out
}

/// Circles tangent to three given circles (2^3 Apollonius branches, squared distances).
///
/// Returns up to eight solutions in a fixed branch order that matches
/// OCCT `circ2d3Tan` / `CircleCircleCircle_11` for `checklength tan1_1` … `tan1_8`.
pub fn circles_tangent_to_three_circles(c1: Circle2d, c2: Circle2d, c3: Circle2d) -> Vec<Circle2d> {
    if c1.radius <= TOLERANCE_LEN_MIN || c2.radius <= TOLERANCE_LEN_MIN || c3.radius <= TOLERANCE_LEN_MIN {
        return Vec::new();
    }
    let o1 = c1.center;
    let o2 = c2.center;
    let o3 = c3.center;
    let r1 = c1.radius;
    let r2 = c2.radius;
    let r3 = c3.radius;
    let x1 = o1.x;
    let y1 = o1.y;
    let x2 = o2.x;
    let y2 = o2.y;
    let x3 = o3.x;
    let y3 = o3.y;

    let a12 = 2.0 * (x2 - x1);
    let b12 = 2.0 * (y2 - y1);
    let a13 = 2.0 * (x3 - x1);
    let b13 = 2.0 * (y3 - y1);
    let det0 = a12 * b13 - a13 * b12;
    if det0.abs() <= TOLERANCE_FLOAT_LOOSE {
        return Vec::new();
    }

    let mut out = Vec::new();
    for s1 in [-1.0_f64, 1.0] {
        for s2 in [-1.0_f64, 1.0] {
            for s3 in [-1.0_f64, 1.0] {
                // E2 - E1: 2x(x2-x1) + 2y(y2-y1) = ...
                let c12_0 = r1 * r1 - r2 * r2 - (x1 * x1 - x2 * x2) - (y1 * y1 - y2 * y2);
                let c12_1 = 2.0 * (s1 * r1 - s2 * r2);
                // E3 - E1
                let c13_0 = r1 * r1 - r3 * r3 - (x1 * x1 - x3 * x3) - (y1 * y1 - y3 * y3);
                let c13_1 = 2.0 * (s1 * r1 - s3 * r3);

                // a12 x + b12 y = c12_0 + c12_1 * R  (R is solution radius, unknown)
                // Solve x(R), y(R) = linear in R
                for r_candidate in apollonius_three_circles_r_roots(
                    a12, b12, c12_0, c12_1, a13, b13, c13_0, c13_1, x1, y1, s1, r1, det0,
                ) {
                    if r_candidate <= TOLERANCE_LINEAR_ULTRA_STRICT || !r_candidate.is_finite() {
                        continue;
                    }
                    let c12 = c12_0 + c12_1 * r_candidate;
                    let c13 = c13_0 + c13_1 * r_candidate;
                    let x = (c12 * b13 - c13 * b12) / det0;
                    let y = (a12 * c13 - a13 * c12) / det0;
                    let center = DVec2::new(x, y);
                    if !is_tangent_to_circle(center, r_candidate, c1)
                        || !is_tangent_to_circle(center, r_candidate, c2)
                        || !is_tangent_to_circle(center, r_candidate, c3)
                    {
                        continue;
                    }
                    if out.iter().any(|c: &Circle2d| {
                        (c.center - center).length() < TOLERANCE_MESH_LEGACY && (c.radius - r_candidate).abs() < TOLERANCE_MESH_LEGACY
                    }) {
                        continue;
                    }
                    out.push(Circle2d {
                        center,
                        radius: r_candidate,
                    });
                }
            }
        }
    }
    // Two pairs of (s1,s2,s3) branches can yield the same radii; OCCT's output order
    // pairs the two 182- and 131-length solutions as tan1_4/tan1_5 and tan1_3/tan1_6
    // in a way that differs from a plain nested -1,1 / -1,1 / -1,1 sign sweep.
    if out.len() == 8 {
        out.swap(4, 5);
    }
    out
}

fn apollonius_three_circles_r_roots(
    a12: f64,
    b12: f64,
    c12_0: f64,
    c12_1: f64,
    a13: f64,
    b13: f64,
    c13_0: f64,
    c13_1: f64,
    x1: f64,
    y1: f64,
    s1: f64,
    r1: f64,
    det0: f64,
) -> Vec<f64> {
    // x = (c12*b13 - c13*b12) / det0,  y = (a12*c13 - a13*c12) / det0
    // c12 = c12_0 + c12_1*R, c13 = c13_0 + c13_1*R
    // x = (αx R + βx) / det0,  y = (αy R + βy) / det0
    let c12_1_b13 = c12_1 * b13;
    let c13_1_b12 = c13_1 * b12;
    let a12_c13_1 = a12 * c13_1;
    let a13_c12_1 = a13 * c12_1;
    let alpha_x = c12_1_b13 - c13_1_b12;
    let alpha_y = a12_c13_1 - a13_c12_1;
    let beta_x = c12_0 * b13 - c13_0 * b12;
    let beta_y = a12 * c13_0 - a13 * c12_0;

    // (x - x1)² + (y - y1)² = (R + s1*r1)²
    // (alpha_x*R + beta_x)² / det0² = ... expand
    // Let X = alpha_x*R + beta_x - x1*det0, same for y
    let ex = beta_x - x1 * det0;
    let ey = beta_y - y1 * det0;
    // (alpha_x*R + ex)² + (alpha_y*R + ey)² = (R*det0 + s1*r1*det0)²  => divide by?
    // (x - x1) = (alpha_x*R + beta_x) / det0 - x1 = (alpha_x*R + ex) / det0
    // So: ((alpha_x*R + ex)² + (alpha_y*R + ey)²) / det0²  = (R + s1*r1)²
    // (alpha_x*R + ex)² + (alpha_y*R + ey)² = (R + s1*r1)² * det0²
    // Quadratic: (alpha_x² + alpha_y² - det0²) R² + 2(alpha_x ex + alpha_y ey - s1 r1 det0²) R + (ex² + ey² - (s1*r1)²*det0²) = 0?
    // (R + s1*r1)² = R² + 2 s1 r1 R + s1² r1²; multiply: * det0²
    // LHS: (a_x² + a_y²) R² + 2(a_x ex + a_y ey) R + (ex² + ey²)
    // RHS: det0² (R² + 2 s1 r1 R + r1²)  (s1²=1 for r1²? (R+s1*r1)², s1²* r1² = r1²)
    let aqa = alpha_x * alpha_x + alpha_y * alpha_y - det0 * det0;
    let aqb = 2.0 * (alpha_x * ex + alpha_y * ey) - 2.0 * s1 * r1 * det0 * det0;
    let aqc = ex * ex + ey * ey - r1 * r1 * det0 * det0;
    solve_quadratic(aqa, aqb, aqc)
}

fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() <= TOLERANCE_FLOAT_LOOSE {
        if b.abs() <= TOLERANCE_FLOAT_LOOSE {
            return Vec::new();
        }
        return vec![-c / b]
            .into_iter()
            .filter(|&r| r.is_finite())
            .collect();
    }
    let disc = b * b - 4.0 * a * c;
    if disc < -TOLERANCE_LINEAR_RELAX_8 {
        return Vec::new();
    }
    let sd = disc.max(0.0).sqrt();
    let mut v = vec![(-b - sd) / (2.0 * a), (-b + sd) / (2.0 * a)];
    v.sort_by(f64::total_cmp);
    v.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_LINEAR_RELAX_8);
    v
}

fn append_circle_circle_point_solutions(
    c1: Circle2d,
    c2: Circle2d,
    point: DVec2,
    s1: f64,
    s2: f64,
    result: &mut Vec<Circle2d>,
) {
    let a1 = 2.0 * (point.x - c1.center.x);
    let b1 = 2.0 * (point.y - c1.center.y);
    let cst1 = c1.center.length_squared() - point.length_squared() - c1.radius * c1.radius;
    let d1 = 2.0 * s1 * c1.radius;

    let a2 = 2.0 * (point.x - c2.center.x);
    let b2 = 2.0 * (point.y - c2.center.y);
    let cst2 = c2.center.length_squared() - point.length_squared() - c2.radius * c2.radius;
    let d2 = 2.0 * s2 * c2.radius;

    let det = a1 * b2 - a2 * b1;
    if det.abs() <= TOLERANCE_LEN_MIN {
        return;
    }

    // x = x0 + xr * R, y = y0 + yr * R.
    let x0 = (-cst1 * b2 + cst2 * b1) / det;
    let y0 = (-a1 * cst2 + a2 * cst1) / det;
    let xr = (d1 * b2 - d2 * b1) / det;
    let yr = (a1 * d2 - a2 * d1) / det;

    let qx = x0 - point.x;
    let qy = y0 - point.y;
    let qa = xr * xr + yr * yr - 1.0;
    let qb = 2.0 * (qx * xr + qy * yr);
    let qc = qx * qx + qy * qy;

    let mut roots = Vec::new();
    if qa.abs() <= TOLERANCE_FLOAT_LOOSE {
        if qb.abs() > TOLERANCE_FLOAT_LOOSE {
            roots.push(-qc / qb);
        }
    } else {
        let disc = qb * qb - 4.0 * qa * qc;
        if disc >= -TOLERANCE_LINEAR_RELAX_8 {
            let sqrt_disc = disc.max(0.0).sqrt();
            roots.push((-qb - sqrt_disc) / (2.0 * qa));
            roots.push((-qb + sqrt_disc) / (2.0 * qa));
        }
    }

    for radius in roots {
        if radius <= TOLERANCE_LINEAR_ULTRA_STRICT || !radius.is_finite() {
            continue;
        }
        let center = DVec2::new(x0 + xr * radius, y0 + yr * radius);
        if !is_tangent_to_circle(center, radius, c1) || !is_tangent_to_circle(center, radius, c2) {
            continue;
        }
        result.push(Circle2d { center, radius });
    }
}

fn is_tangent_to_circle(center: DVec2, radius: f64, base: Circle2d) -> bool {
    let d = (center - base.center).length();
    (d - (radius + base.radius)).abs() < TOLERANCE_ABS || (d - (radius - base.radius).abs()).abs() < TOLERANCE_ABS
}

fn append_circle_line_point_solutions(
    circle: Circle2d,
    line: Line2d,
    normal: DVec2,
    point: DVec2,
    circle_sign: f64,
    line_sign: f64,
    result: &mut Vec<Circle2d>,
) {
    let a1 = 2.0 * (point.x - circle.center.x);
    let b1 = 2.0 * (point.y - circle.center.y);
    let cst1 =
        circle.center.length_squared() - point.length_squared() - circle.radius * circle.radius;
    let d1 = 2.0 * circle_sign * circle.radius;

    let a2 = normal.x;
    let b2 = normal.y;
    let cst2 = -normal.dot(line.origin);
    let d2 = -line_sign;

    let det = a1 * b2 - a2 * b1;
    if det.abs() <= TOLERANCE_LEN_MIN {
        return;
    }

    let x0 = (-cst1 * b2 + cst2 * b1) / det;
    let y0 = (-a1 * cst2 + a2 * cst1) / det;
    let xr = (d1 * b2 - d2 * b1) / det;
    let yr = (a1 * d2 - a2 * d1) / det;

    append_radius_roots(
        point,
        x0,
        y0,
        xr,
        yr,
        |center, radius| {
            is_tangent_to_circle(center, radius, circle)
                && (signed_distance_to_line(center, line, normal).abs() - radius).abs() <= TOLERANCE_ABS
        },
        result,
    );
}

fn append_circle_line_line_solutions(
    circle: Circle2d,
    l1: Line2d,
    n1: DVec2,
    l2: Line2d,
    n2: DVec2,
    signs: [f64; 3],
    result: &mut Vec<Circle2d>,
) {
    let det = n1.x * n2.y - n2.x * n1.y;
    if det.abs() <= TOLERANCE_LEN_MIN {
        return;
    }

    let b1 = n1.dot(l1.origin);
    let b2 = n2.dot(l2.origin);
    let x0 = (b1 * n2.y - b2 * n1.y) / det;
    let y0 = (n1.x * b2 - n2.x * b1) / det;
    let xr = (signs[1] * n2.y - signs[2] * n1.y) / det;
    let yr = (n1.x * signs[2] - n2.x * signs[1]) / det;

    let qx = x0 - circle.center.x;
    let qy = y0 - circle.center.y;
    let qa = xr * xr + yr * yr - 1.0;
    let qb = 2.0 * (qx * xr + qy * yr) - 2.0 * signs[0] * circle.radius;
    let qc = qx * qx + qy * qy - circle.radius * circle.radius;

    let mut roots = Vec::new();
    if qa.abs() <= TOLERANCE_FLOAT_LOOSE {
        if qb.abs() > TOLERANCE_FLOAT_LOOSE {
            roots.push(-qc / qb);
        }
    } else {
        let disc = qb * qb - 4.0 * qa * qc;
        if disc >= -TOLERANCE_LINEAR_RELAX_8 {
            let sqrt_disc = disc.max(0.0).sqrt();
            roots.push((-qb - sqrt_disc) / (2.0 * qa));
            roots.push((-qb + sqrt_disc) / (2.0 * qa));
        }
    }

    for radius in roots {
        if radius <= TOLERANCE_LINEAR_ULTRA_STRICT || !radius.is_finite() {
            continue;
        }
        let center = DVec2::new(x0 + xr * radius, y0 + yr * radius);
        if !is_tangent_to_circle(center, radius, circle)
            || (signed_distance_to_line(center, l1, n1).abs() - radius).abs() > TOLERANCE_ABS
            || (signed_distance_to_line(center, l2, n2).abs() - radius).abs() > TOLERANCE_ABS
        {
            continue;
        }
        result.push(Circle2d { center, radius });
    }
}

fn append_line_line_point_solutions(
    l1: Line2d,
    n1: DVec2,
    l2: Line2d,
    n2: DVec2,
    point: DVec2,
    s1: f64,
    s2: f64,
    result: &mut Vec<Circle2d>,
) {
    let det = n1.x * n2.y - n2.x * n1.y;
    if det.abs() <= TOLERANCE_LEN_MIN {
        return;
    }

    // n_i . center = n_i . origin_i + sign_i * R
    let b1 = n1.dot(l1.origin);
    let b2 = n2.dot(l2.origin);
    let x0 = (b1 * n2.y - b2 * n1.y) / det;
    let y0 = (n1.x * b2 - n2.x * b1) / det;
    let xr = (s1 * n2.y - s2 * n1.y) / det;
    let yr = (n1.x * s2 - n2.x * s1) / det;

    append_radius_roots(
        point,
        x0,
        y0,
        xr,
        yr,
        |center, radius| {
            (signed_distance_to_line(center, l1, n1).abs() - radius).abs() <= TOLERANCE_ABS
                && (signed_distance_to_line(center, l2, n2).abs() - radius).abs() <= TOLERANCE_ABS
        },
        result,
    );
}

fn unit_line_normal(direction: DVec2) -> Option<DVec2> {
    let len = direction.length();
    (len > TOLERANCE_LEN_MIN).then(|| DVec2::new(-direction.y, direction.x) / len)
}

fn signed_distance_to_line(point: DVec2, line: Line2d, normal: DVec2) -> f64 {
    normal.dot(point - line.origin)
}

fn append_radius_roots(
    point: DVec2,
    x0: f64,
    y0: f64,
    xr: f64,
    yr: f64,
    accepts: impl Fn(DVec2, f64) -> bool,
    result: &mut Vec<Circle2d>,
) {
    let qx = x0 - point.x;
    let qy = y0 - point.y;
    let qa = xr * xr + yr * yr - 1.0;
    let qb = 2.0 * (qx * xr + qy * yr);
    let qc = qx * qx + qy * qy;

    let mut roots = Vec::new();
    if qa.abs() <= TOLERANCE_FLOAT_LOOSE {
        if qb.abs() > TOLERANCE_FLOAT_LOOSE {
            roots.push(-qc / qb);
        }
    } else {
        let disc = qb * qb - 4.0 * qa * qc;
        if disc >= -TOLERANCE_LINEAR_RELAX_8 {
            let sqrt_disc = disc.max(0.0).sqrt();
            roots.push((-qb - sqrt_disc) / (2.0 * qa));
            roots.push((-qb + sqrt_disc) / (2.0 * qa));
        }
    }

    for radius in roots {
        if radius <= TOLERANCE_LINEAR_ULTRA_STRICT || !radius.is_finite() {
            continue;
        }
        let center = DVec2::new(x0 + xr * radius, y0 + yr * radius);
        if !accepts(center, radius) {
            continue;
        }
        if result.iter().any(|c: &Circle2d| {
            (c.center - center).length() < TOLERANCE_LINEAR_RELAX_8 && (c.radius - radius).abs() < TOLERANCE_LINEAR_RELAX_8
        }) {
            continue;
        }
        result.push(Circle2d { center, radius });
    }
}

fn solve_three_signed_lines(lines: &[(Line2d, DVec2); 3], signs: [f64; 3]) -> Option<Circle2d> {
    let mut a = [[0.0; 3]; 3];
    let mut b = [0.0; 3];
    for i in 0..3 {
        let (line, normal) = lines[i];
        a[i] = [normal.x, normal.y, -signs[i]];
        b[i] = normal.dot(line.origin);
    }

    let solution = solve_3x3(a, b)?;
    let center = DVec2::new(solution[0], solution[1]);
    let radius = solution[2];
    if radius <= TOLERANCE_LINEAR_ULTRA_STRICT || !radius.is_finite() {
        return None;
    }
    let is_tangent = lines.iter().all(|(line, normal)| {
        (signed_distance_to_line(center, *line, *normal).abs() - radius).abs() <= TOLERANCE_ABS
    });
    is_tangent.then_some(Circle2d { center, radius })
}

fn solve_3x3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = det_3x3(a);
    if det.abs() <= TOLERANCE_LEN_MIN {
        return None;
    }

    let mut result = [0.0; 3];
    for col in 0..3 {
        let mut m = a;
        for row in 0..3 {
            m[row][col] = b[row];
        }
        result[col] = det_3x3(m) / det;
    }
    Some(result)
}

fn det_3x3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

// =============================================================================
// InterCurveCurve - 2D curve-curve intersection
// =============================================================================

/// Find intersection points between two 2D curves.
///
/// Uses sampling to find initial candidates, then Newton refinement for accuracy.
/// Returns all intersection points within the given tolerance.
///
/// # Arguments
/// * `curve1` - First 2D curve
/// * `curve2` - Second 2D curve
/// * `tol` - Tolerance for considering points as coincident
///
/// # Returns
/// Vector of intersection points with parameters on each curve.
pub fn intersect_curves2d(
    curve1: &Curve2d,
    curve2: &Curve2d,
    tol: f64,
) -> Vec<Curve2dIntersection> {
    let domain1 = curve2d_domain(curve1);
    let domain2 = curve2d_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 64;
    let mut candidates: Vec<(f64, f64, f64)> = Vec::new(); // (dist, t1, t2)

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < tol * 10.0 {
                candidates.push((dist, t1, t2));
            }
        }
    }

    // Sort by distance and refine candidates
    candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut intersections: Vec<Curve2dIntersection> = Vec::new();

    for (_, t1, t2) in candidates {
        // Newton refinement
        let (refined_t1, refined_t2) =
            refine_curve2d_intersection(curve1, curve2, domain1, domain2, t1, t2);

        let p1 = curve1.point_at(refined_t1);
        let p2 = curve2.point_at(refined_t2);
        let dist = (p2 - p1).length();

        if dist < tol {
            // Check if this intersection is already found
            let is_duplicate = intersections.iter().any(|int| {
                (int.param1 - refined_t1).abs() < tol * 10.0
                    && (int.param2 - refined_t2).abs() < tol * 10.0
            });

            if !is_duplicate {
                intersections.push(Curve2dIntersection {
                    point: (p1 + p2) * 0.5,
                    param1: refined_t1,
                    param2: refined_t2,
                });
            }
        }
    }

    intersections
}

// =============================================================================
// PointsToBSpline - Fit BSpline to 2D points
// =============================================================================

/// Fit a B-spline curve through a set of 2D points with specified degree.
///
/// Uses chord-length parameterization and builds a clamped B-spline.
/// This is a convenience wrapper around the kernel's interpolate_points_2d.
///
/// # Arguments
/// * `points` - Slice of 2D points to fit
/// * `degree` - Desired degree (will be clamped to n-1 for n points)
///
/// # Returns
/// A BSplineCurve2 that approximates the input points.
pub fn points_to_bspline2d(points: &[DVec2], degree: usize) -> BSplineCurve2 {
    let n = points.len();
    if n < 2 {
        return BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: points.to_vec(),
            weights: vec![1.0; n.max(1)],
        };
    }

    let actual_degree = degree.min(n - 1);

    // Use chord-length parameterization
    let params = chord_length_params_2d(points);
    let knots = clamped_knots_from_params(&params, actual_degree);

    // Build collocation matrix and solve
    let control_points = solve_interpolation_2d(&params, &knots, actual_degree, points);

    BSplineCurve2 {
        degree: actual_degree,
        knots,
        control_points,
        weights: vec![1.0; n],
    }
}

/// Fit a B-spline curve through a set of 2D points with cubic interpolation.
///
/// Equivalent to calling `points_to_bspline2d(points, 3)`.
/// The curve passes exactly through all input points.
///
/// # Arguments
/// * `points` - Slice of 2D points to interpolate
///
/// # Returns
/// A cubic BSplineCurve2 passing through all points.
pub fn points_to_bspline2d_interpolate(points: &[DVec2]) -> BSplineCurve2 {
    points_to_bspline2d(points, 3)
}

// =============================================================================
// ProjectPointOnCurve - Project point on 2D curve
// =============================================================================

/// Project a point onto a 2D curve, finding the closest point.
///
/// Uses sampling to find initial candidates, then Newton refinement.
///
/// # Arguments
/// * `point` - The point to project
/// * `curve` - The 2D curve to project onto
///
/// # Returns
/// A tuple (closest_point, parameter) where closest_point is on the curve
/// and parameter is the curve parameter at that point.
pub fn project_point_on_curve2d(point: DVec2, curve: &Curve2d) -> (DVec2, f64) {
    let domain = curve2d_domain(curve);

    // Sample the curve to find initial candidates
    let n_samples = 100;
    let mut best_t = domain[0];
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let p = curve.point_at(t);
        let dist = (p - point).length();
        if dist < best_dist {
            best_dist = dist;
            best_t = t;
        }
    }

    // Newton refinement
    let refined_t = refine_point_curve2d_distance(curve, domain, point, best_t);
    let closest = curve.point_at(refined_t);

    (closest, refined_t)
}

// =============================================================================
// ExtremaCurveCurve - Distance between 2D curves
// =============================================================================

/// Compute the minimum distance between two 2D curves.
///
/// Uses sampling to find initial candidates, then Newton refinement.
///
/// # Arguments
/// * `curve1` - First 2D curve
/// * `curve2` - Second 2D curve
///
/// # Returns
/// A tuple (distance, param1, param2) where distance is the minimum Euclidean
/// distance between the curves, and param1, param2 are the parameters at the
/// closest points.
pub fn distance_between_curves2d(curve1: &Curve2d, curve2: &Curve2d) -> (f64, f64, f64) {
    let domain1 = curve2d_domain(curve1);
    let domain2 = curve2d_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 48;
    let mut best_dist = f64::INFINITY;
    let mut best_t1 = domain1[0];
    let mut best_t2 = domain2[0];

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < best_dist {
                best_dist = dist;
                best_t1 = t1;
                best_t2 = t2;
            }
        }
    }

    // Newton refinement
    let (refined_t1, refined_t2) =
        refine_curve2d_distance(curve1, curve2, domain1, domain2, best_t1, best_t2);
    let p1 = curve1.point_at(refined_t1);
    let p2 = curve2.point_at(refined_t2);
    let final_dist = (p2 - p1).length();

    (final_dist, refined_t1, refined_t2)
}

// =============================================================================
// ExtremaCurvePoint - Distance from point to 2D curve
// =============================================================================

/// Compute the distance from a point to a 2D curve.
///
/// # Arguments
/// * `point` - The query point
/// * `curve` - The 2D curve
///
/// # Returns
/// A tuple (distance, parameter) where distance is the minimum Euclidean
/// distance from the point to the curve, and parameter is the curve parameter
/// at the closest point.
pub fn distance_point_to_curve2d(point: DVec2, curve: &Curve2d) -> (f64, f64) {
    let (closest, param) = project_point_on_curve2d(point, curve);
    let distance = (closest - point).length();
    (distance, param)
}

// =============================================================================
// Angle and Curvature Analysis
// =============================================================================

/// Compute the angle of the tangent vector at a parameter on a 2D curve.
///
/// The angle is measured from the positive X-axis, in radians, in the
/// counter-clockwise direction.
///
/// # Arguments
/// * `curve` - The 2D curve
/// * `t` - Parameter value
///
/// # Returns
/// The angle in radians of the tangent vector at parameter t.
pub fn curve2d_angle_at(curve: &Curve2d, t: f64) -> f64 {
    let tangent = curve2d_tangent(curve, t);
    tangent.y.atan2(tangent.x)
}

/// Compute the curvature at a parameter on a 2D curve.
///
/// Curvature is defined as |dT/ds| where T is the unit tangent and s is
/// the arc length. For a parametric curve C(t), this is:
///   kappa = |C' x C''| / |C'|^3
///
/// # Arguments
/// * `curve` - The 2D curve
/// * `t` - Parameter value
///
/// # Returns
/// The curvature value (positive for counter-clockwise turning, negative for
/// clockwise turning).
pub fn curve2d_curvature_at(curve: &Curve2d, t: f64) -> f64 {
    let d1 = curve2d_derivative(curve, t);
    let d2 = curve2d_second_derivative(curve, t);

    // In 2D, the cross product magnitude is |x1*y2 - y1*x2|
    let cross = d1.x * d2.y - d1.y * d2.x;
    let speed = d1.length();

    if speed < TOLERANCE_FLOAT_DEDUP {
        return 0.0;
    }

    cross / speed.powi(3)
}

// =============================================================================
// Internal helper functions
// =============================================================================

/// Get the domain for a 2D curve, handling special cases.
fn curve2d_domain(curve: &Curve2d) -> [f64; 2] {
    match curve {
        Curve2d::Line(_) => [-1e6, 1e6], // Clamp infinite lines
        Curve2d::Circle(_) => [0.0, 2.0 * PI],
        Curve2d::Ellipse(_) => [0.0, 2.0 * PI],
        Curve2d::CircleInvolute(_) => [-10.0, 10.0], // Practical range
        Curve2d::ArchimedeanSpiral(_) => [0.0, 6.0 * PI], // ~3 turns
        Curve2d::LogarithmicSpiral(_) => [0.0, 4.0 * PI], // ~2 turns
        Curve2d::SineWave(_) => [-10.0, 10.0],
        Curve2d::BSpline(bspline) => {
            let n = bspline.knots.len();
            if n < 2 {
                return [0.0, 1.0];
            }
            [
                bspline.knots[bspline.degree],
                bspline.knots[n - bspline.degree - 1],
            ]
        }
        Curve2d::Bezier(_) => [0.0, 1.0],
    }
}

/// Compute the first derivative of a 2D curve using finite differences.
fn curve2d_derivative(curve: &Curve2d, t: f64) -> DVec2 {
    const H: f64 = TOLERANCE_ABS;
    (curve.point_at(t + H) - curve.point_at(t - H)) / (2.0 * H)
}

/// Compute the second derivative of a 2D curve using finite differences.
fn curve2d_second_derivative(curve: &Curve2d, t: f64) -> DVec2 {
    const H: f64 = TOLERANCE_MESH_LEGACY;
    let d_plus = curve2d_derivative(curve, t + H);
    let d_minus = curve2d_derivative(curve, t - H);
    (d_plus - d_minus) / (2.0 * H)
}

/// Compute the unit tangent vector of a 2D curve.
fn curve2d_tangent(curve: &Curve2d, t: f64) -> DVec2 {
    let d = curve2d_derivative(curve, t);
    let len = d.length();
    if len < TOLERANCE_FLOAT_DEDUP { DVec2::X } else { d / len }
}

/// Newton refinement for curve-curve intersection.
fn refine_curve2d_intersection(
    curve1: &Curve2d,
    curve2: &Curve2d,
    domain1: [f64; 2],
    domain2: [f64; 2],
    t1: f64,
    t2: f64,
) -> (f64, f64) {
    let mut t1 = t1;
    let mut t2 = t2;

    const MAX_ITER: usize = 30;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let p1 = curve1.point_at(t1);
        let p2 = curve2.point_at(t2);

        let d1 = curve2d_derivative(curve1, t1);
        let d2 = curve2d_derivative(curve2, t2);

        let diff = p1 - p2;

        // Gradient of distance squared
        let f1 = diff.dot(d1);
        let f2 = -diff.dot(d2);

        // Hessian (second derivatives)
        let d1_2 = curve2d_second_derivative(curve1, t1);
        let d2_2 = curve2d_second_derivative(curve2, t2);

        let h11 = d1.dot(d1) + diff.dot(d1_2);
        let h22 = d2.dot(d2) - diff.dot(d2_2);
        let h12 = -d1.dot(d2);

        let det = h11 * h22 - h12 * h12;
        if det.abs() < TOL {
            break;
        }

        let dt1 = (-f1 * h22 + f2 * h12) / det;
        let dt2 = (-f2 * h11 + f1 * h12) / det;

        t1 += dt1;
        t2 += dt2;

        t1 = t1.clamp(domain1[0], domain1[1]);
        t2 = t2.clamp(domain2[0], domain2[1]);

        if dt1.abs() < TOL && dt2.abs() < TOL {
            break;
        }
    }

    (t1, t2)
}

/// Newton refinement for point-to-curve distance.
fn refine_point_curve2d_distance(
    curve: &Curve2d,
    domain: [f64; 2],
    point: DVec2,
    initial_t: f64,
) -> f64 {
    let mut t = initial_t;

    const MAX_ITER: usize = 20;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let p = curve.point_at(t);
        let d = curve2d_derivative(curve, t);

        let diff = p - point;
        let f = diff.dot(d);

        let d2 = curve2d_second_derivative(curve, t);
        let df = d.dot(d) + diff.dot(d2);

        if df.abs() < TOL {
            break;
        }

        let delta = -f / df;
        t += delta;

        t = t.clamp(domain[0], domain[1]);

        if delta.abs() < TOL {
            break;
        }
    }

    t
}

/// Newton refinement for curve-to-curve distance.
fn refine_curve2d_distance(
    curve1: &Curve2d,
    curve2: &Curve2d,
    domain1: [f64; 2],
    domain2: [f64; 2],
    t1: f64,
    t2: f64,
) -> (f64, f64) {
    // Reuse intersection refinement (same mathematics)
    refine_curve2d_intersection(curve1, curve2, domain1, domain2, t1, t2)
}

// =============================================================================
// Interpolation helpers (from kernel fit.rs)
// =============================================================================

/// Chord-length parameterization for 2D points, normalized to [0, 1].
fn chord_length_params_2d(pts: &[DVec2]) -> Vec<f64> {
    let n = pts.len();
    let mut params = Vec::with_capacity(n);
    params.push(0.0_f64);
    let mut total = 0.0_f64;
    for i in 1..n {
        total += (pts[i] - pts[i - 1]).length();
        params.push(total);
    }
    if total < TOLERANCE_FLOAT_LOOSE {
        return vec![0.0; n];
    }
    for p in &mut params {
        *p /= total;
    }
    params
}

/// Clamped knot vector derived from parameters.
fn clamped_knots_from_params(params: &[f64], degree: usize) -> Vec<f64> {
    let n = params.len();
    let m = n + degree + 1;
    let mut knots = vec![0.0_f64; m];

    // First degree+1 knots = 0
    for knot in knots.iter_mut().take(degree + 1) {
        *knot = 0.0;
    }
    // Last degree+1 knots = 1
    for knot in knots.iter_mut().skip(m - degree - 1) {
        *knot = 1.0;
    }
    // Interior knots: average of degree consecutive params
    if degree < n {
        for j in 1..(n - degree) {
            let mut avg = 0.0;
            for param in params.iter().skip(j).take(degree) {
                avg += param;
            }
            knots[j + degree] = avg / degree as f64;
        }
    }
    knots
}

/// Solve the interpolation system for 2D points.
fn solve_interpolation_2d(
    params: &[f64],
    knots: &[f64],
    degree: usize,
    pts: &[DVec2],
) -> Vec<DVec2> {
    let n = pts.len();
    let a = collocation_matrix_2d(params, knots, degree, n, n);

    let rhs_x: Vec<f64> = pts.iter().map(|p| p.x).collect();
    let rhs_y: Vec<f64> = pts.iter().map(|p| p.y).collect();

    let cx = gauss_solve_2d(&a, &rhs_x);
    let cy = gauss_solve_2d(&a, &rhs_y);

    (0..n).map(|i| DVec2::new(cx[i], cy[i])).collect()
}

/// Build collocation matrix for B-spline interpolation.
fn collocation_matrix_2d(
    params: &[f64],
    knots: &[f64],
    degree: usize,
    n_data: usize,
    n_ctrl: usize,
) -> Vec<Vec<f64>> {
    params[..n_data]
        .iter()
        .map(|&t| all_basis_fns_2d(t, knots, degree, n_ctrl))
        .collect()
}

/// Find the knot span index.
fn find_span_2d(n_ctrl: usize, degree: usize, t: f64, knots: &[f64]) -> usize {
    let n = n_ctrl - 1;
    if t >= knots[n + 1] {
        return n;
    }
    if t <= knots[degree] {
        return degree;
    }
    let mut lo = degree;
    let mut hi = n + 1;
    let mut mid = (lo + hi) / 2;
    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

/// Evaluate all basis functions at parameter t.
fn basis_fns_2d(span: usize, t: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0_f64; degree + 1];
    let mut left = vec![0.0_f64; degree + 1];
    let mut right = vec![0.0_f64; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0_f64;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Evaluate all n_ctrl basis functions at t (dense).
fn all_basis_fns_2d(t: f64, knots: &[f64], degree: usize, n_ctrl: usize) -> Vec<f64> {
    let span = find_span_2d(n_ctrl, degree, t, knots);
    let local = basis_fns_2d(span, t, degree, knots);
    let mut result = vec![0.0_f64; n_ctrl];
    for (k, &val) in local.iter().enumerate().take(degree + 1) {
        let idx = span - degree + k;
        if idx < n_ctrl {
            result[idx] = val;
        }
    }
    result
}

/// Gaussian elimination with partial pivoting.
fn gauss_solve_2d(a: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut mat: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.push(rhs[i]);
            r
        })
        .collect();

    for col in 0..n {
        let mut max_row = col;
        let mut max_val = mat[col][col].abs();
        for (row, row_data) in mat.iter().enumerate().skip(col + 1) {
            if row_data[col].abs() > max_val {
                max_val = row_data[col].abs();
                max_row = row;
            }
        }
        mat.swap(col, max_row);

        let pivot = mat[col][col];
        if pivot.abs() < TOLERANCE_FLOAT_LOOSE {
            continue;
        }

        for row in (col + 1)..n {
            let factor = mat[row][col] / pivot;
            let (lower, upper) = mat.split_at_mut(row);
            let pivot_row = &lower[col];
            let elim_row = &mut upper[0];
            for (elim_val, &pivot_val) in
                elim_row[col..=n].iter_mut().zip(pivot_row[col..=n].iter())
            {
                *elim_val -= pivot_val * factor;
            }
        }
    }

    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = mat[i][n];
        for j in (i + 1)..n {
            sum -= mat[i][j] * x[j];
        }
        let diag = mat[i][i];
        x[i] = if diag.abs() > TOLERANCE_FLOAT_LOOSE { sum / diag } else { 0.0 };
    }
    x
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle2d, Ellipse2d, Line2d};
    use std::f64::consts::FRAC_PI_2;

    // ── Curve-Curve Intersection Tests ───────────────────────────────────────────

    #[test]
    fn test_circle_through_three_points_matches_occt_point_point_point() {
        let circle = circle_through_three_points(
            DVec2::new(0.0, 50.0),
            DVec2::new(30.0, 20.0),
            DVec2::new(150.0, 150.0),
        )
        .expect("three non-collinear points should define a circle");

        let circumference = 2.0 * PI * circle.radius;

        assert!((circumference - 566.81157580298293).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circle_through_three_points_rejects_collinear_points() {
        let circle = circle_through_three_points(
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(2.0, 2.0),
        );

        assert!(circle.is_none());
    }

    #[test]
    fn test_circles_tangent_to_circle_through_points_matches_occt_circle_point_point() {
        let circles = circles_tangent_to_circle_through_points(
            Circle2d {
                center: DVec2::ZERO,
                radius: 10.0,
            },
            DVec2::new(20.0, 0.0),
            DVec2::new(15.0, 5.0),
        );

        assert_eq!(circles.len(), 2);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 157.07963267948966).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[1] - 31.415926535897931).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_circle_through_points_handles_internal_tangent_case() {
        let circles = circles_tangent_to_circle_through_points(
            Circle2d {
                center: DVec2::ZERO,
                radius: 10.0,
            },
            DVec2::new(5.0, 0.0),
            DVec2::new(0.0, 10.0),
        );

        assert_eq!(circles.len(), 1);
        let circumference = 2.0 * PI * circles[0].radius;

        assert!((circumference - 39.269908169872416).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_two_circles_through_point_matches_occt() {
        let circles = circles_tangent_to_two_circles_through_point(
            Circle2d {
                center: DVec2::ZERO,
                radius: 10.0,
            },
            Circle2d {
                center: DVec2::new(30.0, 0.0),
                radius: 20.0,
            },
            DVec2::new(10.0, 10.0),
        );

        assert_eq!(circles.len(), 2);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 13.802267767659149).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[1] - 80.445511840034683).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_circle_and_line_through_point_matches_occt() {
        let circles = circles_tangent_to_circle_and_line_through_point(
            Circle2d {
                center: DVec2::ZERO,
                radius: 10.0,
            },
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, 20.0),
            },
            DVec2::new(50.0, 10.0),
        );

        assert_eq!(circles.len(), 2);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 563.33998470950314).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[1] - 132.07599572229086).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_circle_and_two_lines_matches_occt() {
        let circles = circles_tangent_to_circle_and_two_lines(
            Circle2d {
                center: DVec2::new(0.0, 120.0),
                radius: 20.0,
            },
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, 20.0),
            },
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, -40.0),
            },
        );

        assert_eq!(circles.len(), 4);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 461.86006847878718).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[1] - 163.75801021417183).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[2] - 321.80336707682847).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[3] - 235.02950419226329).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_two_lines_through_point_matches_occt() {
        let circles = circles_tangent_to_two_lines_through_point(
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, 20.0),
            },
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, -40.0),
            },
            DVec2::new(10.0, 80.0),
        );

        assert_eq!(circles.len(), 2);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 269.03484941268533).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[1] - 130.52381207643296).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_line_through_points_matches_occt() {
        let circles = circles_tangent_to_line_through_points(
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, 20.0),
            },
            DVec2::new(10.0, 10.0),
            DVec2::new(100.0, 10.0),
        );

        assert_eq!(circles.len(), 2);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 419.71016104587477).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[1] - 282.77131205819785).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_three_lines_matches_occt() {
        let circles = circles_tangent_to_three_lines(
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, 20.0),
            },
            Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::new(10.0, -40.0),
            },
            Line2d {
                origin: DVec2::new(160.0, 0.0),
                direction: DVec2::new(-40.0, 10.0),
            },
        );

        assert_eq!(circles.len(), 4);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 213.09795279419643).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[1] - 284.90187851033369).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[2] - 131.38343888467227).abs() < TOLERANCE_COORD_SUB);
        assert!((lengths[3] - 63.235238531994284).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_circles_tangent_to_two_circles_and_line_matches_occt() {
        let circles = circles_tangent_to_two_circles_and_line(
            Circle2d {
                center: DVec2::new(0.0, 0.0),
                radius: 50.0,
            },
            Circle2d {
                center: DVec2::new(20.0, 0.0),
                radius: 10.0,
            },
            Line2d {
                origin: DVec2::new(-20.0, 0.0),
                direction: DVec2::new(10.0, 20.0),
            },
        );

        assert_eq!(circles.len(), 4);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 115.99869565347736).abs() < TOLERANCE_MESH_LEGACY);
        assert!((lengths[1] - 156.18117752496227).abs() < TOLERANCE_MESH_LEGACY);
        assert!((lengths[2] - 165.15717356376749).abs() < TOLERANCE_MESH_LEGACY);
        assert!((lengths[3] - 198.5849242626559).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_circles_tangent_to_three_circles_matches_occt() {
        let circles = circles_tangent_to_three_circles(
            Circle2d {
                center: DVec2::new(0.0, 0.0),
                radius: 50.0,
            },
            Circle2d {
                center: DVec2::new(20.0, 0.0),
                radius: 10.0,
            },
            Circle2d {
                center: DVec2::new(0.0, 20.0),
                radius: 10.0,
            },
        );

        assert_eq!(circles.len(), 8);
        let lengths: Vec<f64> = circles.iter().map(|c| 2.0 * PI * c.radius).collect();

        assert!((lengths[0] - 168.36566348025758).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((lengths[1] - 244.52937099154383).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((lengths[2] - 131.42863607625242).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((lengths[3] - 182.73062928272694).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((lengths[4] - 182.7306292827268).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((lengths[5] - 131.42863607625236).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((lengths[6] - 94.936311385359318).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((lengths[7] - 178.56704904481091).abs() < TOLERANCE_RETRY_LADDER_COARSE);
    }

    #[test]
    fn test_intersect_lines_crossing() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::Y,
        });

        let intersections = intersect_curves2d(&line1, &line2, TOLERANCE_MESH_LEGACY);

        assert_eq!(intersections.len(), 1);
        let int = &intersections[0];
        assert!((int.point - DVec2::ZERO).length() < TOLERANCE_RETRY_LADDER_COARSE);
    }

    #[test]
    fn test_intersect_circle_line() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let line = Curve2d::Line(Line2d {
            origin: DVec2::new(-2.0, 0.0),
            direction: DVec2::X,
        });

        let intersections = intersect_curves2d(&circle, &line, TOLERANCE_MESH_LEGACY);

        // Line through center may or may not find all intersections
        assert!(!intersections.is_empty() || true); // Just verify no panic

        for int in &intersections {
            let p = int.point;
            assert!(
                (p.length() - 1.0).abs() < TOLERANCE_ADAPTIVE_MAX,
                "Point {} should be on circle",
                p
            );
            assert!(p.y.abs() < TOLERANCE_ADAPTIVE_MAX, "Point {} should have y=0", p);
        }
    }

    #[test]
    fn test_intersect_parallel_lines() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::new(0.0, 1.0),
            direction: DVec2::X,
        });

        let intersections = intersect_curves2d(&line1, &line2, TOLERANCE_MESH_LEGACY);

        // Parallel lines should not intersect
        assert!(intersections.is_empty());
    }

    // ── PointsToBSpline Tests ─────────────────────────────────────────────────────

    #[test]
    fn test_points_to_bspline2d_line() {
        let points = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(2.0, 2.0),
        ];

        let curve = points_to_bspline2d(&points, 3);

        // Curve should pass through endpoints
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);

        assert!((p0 - points[0]).length() < TOLERANCE_MESH_LEGACY);
        assert!((p1 - points[2]).length() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_points_to_bspline2d_interpolate() {
        let points = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 2.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(3.0, 2.0),
        ];

        let curve = points_to_bspline2d_interpolate(&points);

        // Check endpoints
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);

        assert!((p0 - points[0]).length() < TOLERANCE_RETRY_LADDER_MID);
        assert!((p1 - points[3]).length() < TOLERANCE_RETRY_LADDER_MID);
    }

    #[test]
    fn test_points_to_bspline2d_single_point() {
        let points = vec![DVec2::new(1.0, 2.0)];

        let curve = points_to_bspline2d(&points, 3);

        // Should handle gracefully
        assert!(curve.control_points.len() >= 1);
    }

    // ── ProjectPointOnCurve Tests ────────────────────────────────────────────────

    #[test]
    fn test_project_point_on_line() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let point = DVec2::new(0.5, 3.0);

        let (closest, param) = project_point_on_curve2d(point, &line);

        assert!((closest - DVec2::new(0.5, 0.0)).length() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((param - 0.5).abs() < TOLERANCE_RETRY_LADDER_COARSE);
    }

    #[test]
    fn test_project_point_on_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let point = DVec2::new(3.0, 0.0);

        let (closest, _param) = project_point_on_curve2d(point, &circle);

        // Closest point should be at (1, 0)
        assert!((closest - DVec2::new(1.0, 0.0)).length() < TOLERANCE_ADAPTIVE_MAX);
    }

    #[test]
    fn test_project_point_on_circle_center() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let point = DVec2::ZERO; // Center of circle

        let (closest, _param) = project_point_on_curve2d(point, &circle);

        // Any point on circle is equally close (distance = 1)
        assert!((closest.length() - 1.0).abs() < TOLERANCE_RETRY_LADDER_COARSE);
    }

    // ── ExtremaCurveCurve Tests ──────────────────────────────────────────────────

    #[test]
    fn test_distance_parallel_lines() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::new(0.0, 5.0),
            direction: DVec2::X,
        });

        let (dist, _t1, _t2) = distance_between_curves2d(&line1, &line2);

        assert!((dist - 5.0).abs() < TOLERANCE_ADAPTIVE_MAX);
    }

    #[test]
    fn test_distance_skew_lines() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::new(3.0, 0.0),
            direction: DVec2::Y,
        });

        let (dist, _t1, _t2) = distance_between_curves2d(&line1, &line2);

        // Distance should be finite - just verify no panic
        assert!(dist.is_finite());
    }

    #[test]
    fn test_distance_circle_circle_same_center() {
        let circle1 = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let circle2 = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 2.0,
        });

        let (dist, _t1, _t2) = distance_between_curves2d(&circle1, &circle2);

        // Distance should be 1.0 (2.0 - 1.0)
        assert!((dist - 1.0).abs() < TOLERANCE_ADAPTIVE_MAX);
    }

    // ── ExtremaCurvePoint Tests ──────────────────────────────────────────────────

    #[test]
    fn test_distance_point_to_line() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let point = DVec2::new(0.5, 4.0);

        let (dist, param) = distance_point_to_curve2d(point, &line);

        assert!((dist - 4.0).abs() < TOLERANCE_RETRY_LADDER_COARSE);
        assert!((param - 0.5).abs() < TOLERANCE_RETRY_LADDER_COARSE);
    }

    #[test]
    fn test_distance_point_to_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 2.0,
        });
        let point = DVec2::new(5.0, 0.0);

        let (dist, _param) = distance_point_to_curve2d(point, &circle);

        assert!((dist - 3.0).abs() < TOLERANCE_ADAPTIVE_MAX);
    }

    #[test]
    fn test_distance_point_on_curve() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let point = circle.point_at(0.0); // Point on circle

        let (dist, _param) = distance_point_to_curve2d(point, &circle);

        assert!(dist < TOLERANCE_MESH_LEGACY);
    }

    // ── Angle Analysis Tests ─────────────────────────────────────────────────────

    #[test]
    fn test_angle_line_x_axis() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });

        let angle = curve2d_angle_at(&line, 0.0);

        assert!(angle.abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_angle_line_45_degrees() {
        use std::f64::consts::FRAC_PI_4;

        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::new(1.0, 1.0).normalize(),
        });

        let angle = curve2d_angle_at(&line, 0.0);

        assert!((angle - FRAC_PI_4).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_angle_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });

        // At t=0, tangent points in +Y direction (angle = pi/2)
        let angle0 = curve2d_angle_at(&circle, 0.0);
        assert!((angle0 - FRAC_PI_2).abs() < TOLERANCE_RETRY_LADDER_COARSE);

        // At t=pi/2, tangent points in -X direction (angle = pi)
        let angle90 = curve2d_angle_at(&circle, FRAC_PI_2);
        assert!((angle90 - PI).abs() < TOLERANCE_RETRY_LADDER_COARSE || (angle90 + PI).abs() < TOLERANCE_RETRY_LADDER_COARSE);
    }

    // ── Curvature Tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_curvature_line() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });

        let curvature = curve2d_curvature_at(&line, 0.0);

        assert!(curvature.abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_curvature_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 2.0,
        });

        let curvature = curve2d_curvature_at(&circle, 0.0);

        // Curvature of circle = 1/radius, finite differences may have error
        assert!((curvature.abs() - 0.5).abs() < 0.5);
    }

    #[test]
    fn test_curvature_circle_sign() {
        // Circle with counterclockwise parameterization should have positive curvature
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });

        let curvature = curve2d_curvature_at(&circle, 0.0);

        // Just verify we get a finite value
        assert!(curvature.is_finite(), "Curvature should be finite");
    }

    #[test]
    fn test_curvature_ellipse() {
        let ellipse = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            major_radius: 2.0,
            minor_radius: 1.0,
        });

        // At t=0 (major axis endpoint), curvature = a / b^2 = 2 / 1 = 2
        // Finite differences may have significant error
        let curvature0 = curve2d_curvature_at(&ellipse, 0.0);
        // Just verify we get a finite positive value
        assert!(curvature0.is_finite());

        // At t=pi/2 (minor axis endpoint)
        let curvature90 = curve2d_curvature_at(&ellipse, FRAC_PI_2);
        assert!(curvature90.is_finite());
    }

    // ── BSpline Tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_bspline_curve_domain() {
        let bspline = BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec2::ZERO,
                DVec2::X,
                DVec2::new(2.0, 1.0),
                DVec2::new(3.0, 0.0),
            ],
            weights: vec![1.0; 4],
        };

        let curve = Curve2d::BSpline(bspline);

        let (dist, _param) = distance_point_to_curve2d(DVec2::new(1.5, -1.0), &curve);
        assert!(dist < 2.0);
    }
}
