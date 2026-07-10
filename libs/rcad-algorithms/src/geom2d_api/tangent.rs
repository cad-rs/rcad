use super::*;

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

    Some(Circle2d { center, radius, x_dir: DVec2::X, y_dir: DVec2::Y })
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
        circles.push(Circle2d { center, radius, x_dir: DVec2::X, y_dir: DVec2::Y });
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
            result.push(Circle2d { center, radius, x_dir: DVec2::X, y_dir: DVec2::Y });
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
                    out.push(Circle2d { center, x_dir: DVec2::X, y_dir: DVec2::Y, radius: r });
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
                        x_dir: DVec2::X,
                        y_dir: DVec2::Y,
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
        result.push(Circle2d { center, radius, x_dir: DVec2::X, y_dir: DVec2::Y });
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
        result.push(Circle2d { center, radius, x_dir: DVec2::X, y_dir: DVec2::Y });
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
        result.push(Circle2d { center, radius, x_dir: DVec2::X, y_dir: DVec2::Y });
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
    is_tangent.then_some(Circle2d { center, radius, x_dir: DVec2::X, y_dir: DVec2::Y })
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

