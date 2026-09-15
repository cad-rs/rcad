//! B-spline / Bezier basis algorithms shared by the evaluation impls: the de
//! Boor recursions (homogeneous and rational, 3D and 2D), the analytic
//! tangent forms and the de Casteljau evaluations.
//!
//! Extracted verbatim from `geom/eval.rs` (project Rule 5 file-size split);
//! every function stays reachable as `crate::geom::eval::<name>` through the
//! `pub(crate) use` re-exports in `eval.rs`.
use glam::{DVec2, DVec3};

/// De Boor's algorithm in homogeneous 4D space.
/// Returns `[wx, wy, wz, w]` (not divided by w yet).
pub(crate) fn de_boor_homo(
    degree: usize,
    knots: &[f64],
    points: &[DVec3],
    weights: &[f64],
    t: f64,
) -> [f64; 4] {
    let n = points.len();
    if n == 0 {
        return [0.0; 4];
    }
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };
    // OCCT BSplCLib: a NULL weight array (non-rational curve) is treated as
    // all weights equal to 1.0. Missing trailing entries are likewise 1.0.
    let mut d: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights.get(idx).copied().unwrap_or(1.0);
            [points[idx].x * w, points[idx].y * w, points[idx].z * w, w]
        })
        .collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }
    d[degree]
}

/// De Boor's algorithm in homogeneous 3D space for 2D rational curves.
/// Returns `[wx, wy, w]` (not divided by w yet).
pub(crate) fn de_boor_homo_2d(
    degree: usize,
    knots: &[f64],
    points: &[DVec2],
    weights: &[f64],
    t: f64,
) -> [f64; 3] {
    let n = points.len();
    if n == 0 {
        return [0.0; 3];
    }
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };
    // OCCT BSplCLib: a NULL weight array (non-rational curve) is treated as
    // all weights equal to 1.0. Missing trailing entries are likewise 1.0.
    let mut d: Vec<[f64; 3]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights.get(idx).copied().unwrap_or(1.0);
            [points[idx].x * w, points[idx].y * w, w]
        })
        .collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }
    d[degree]
}

/// De Boor's algorithm for rational B-spline evaluation.
/// Returns the 3D point at parameter `t`.
pub(crate) fn de_boor(degree: usize, knots: &[f64], points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 {
        return DVec3::ZERO;
    }

    // Find knot span index k such that knots[k] <= t < knots[k+1]
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };

    // Initialize homogeneous control points for the span
    let mut d: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let idx = k - degree + j;
            let idx = idx.min(n - 1);
            let w = weights[idx];
            [points[idx].x * w, points[idx].y * w, points[idx].z * w, w]
        })
        .collect();

    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }

    let w = d[degree][3];
    if w.abs() < 1e-15 {
        DVec3::ZERO
    } else {
        DVec3::new(d[degree][0] / w, d[degree][1] / w, d[degree][2] / w)
    }
}

/// De Boor's algorithm for rational B-spline evaluation in 2D parameter space.
/// Returns the 2D point at parameter `t`. Identical logic to `de_boor` with DVec2.
pub(crate) fn de_boor_2d(degree: usize, knots: &[f64], points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 {
        return DVec2::ZERO;
    }

    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };

    // Homogeneous control points [x*w, y*w, w]
    // OCCT BSplCLib: a NULL weight array (non-rational curve) is treated as
    // all weights equal to 1.0. Missing trailing entries are likewise 1.0.
    let mut d: Vec<[f64; 3]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights.get(idx).copied().unwrap_or(1.0);
            [points[idx].x * w, points[idx].y * w, w]
        })
        .collect();

    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }

    let w = d[degree][2];
    if w.abs() < 1e-15 {
        DVec2::ZERO
    } else {
        DVec2::new(d[degree][0] / w, d[degree][1] / w)
    }
}

/// Analytic tangent for a rational B-Spline curve (NURBS) using the quotient rule.
///
/// The derivative of C(t) = A(t)/W(t) is:
///   C'(t) = (A'(t) - W'(t)*C(t)) / W(t)
///
/// A'(t) and W'(t) are degree-(p-1) B-Splines with control points:
///   A'_i = p * (w_{i+1}*P_{i+1} - w_i*P_i) / (t_{i+p+1} - t_{i+1})
///   W'_i = p * (w_{i+1} - w_i)              / (t_{i+p+1} - t_{i+1})
///
/// Returns the unnormalised derivative vector (caller normalises if needed).
pub(crate) fn bspline_tangent_analytic(
    degree: usize,
    knots: &[f64],
    points: &[DVec3],
    weights: &[f64],
    t: f64,
) -> DVec3 {
    let n = points.len();
    if n < 2 || degree == 0 {
        return DVec3::ZERO;
    }

    // OCCT: a NULL weight array (non-rational curve) is treated as all weights
    // equal to 1.0 (BSplCLib weight accessor, cf. math/bspl.rs wgt).
    let ws: Vec<f64> = if weights.is_empty() { vec![1.0; n] } else { weights.to_vec() };

    let p = degree as f64;
    let m = n - 1; // number of derivative control points

    let mut a_prime: Vec<DVec3> = Vec::with_capacity(m);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(m); // scalar stored in .x
    for i in 0..m {
        let denom = knots[i + degree + 1] - knots[i + 1];
        if denom.abs() < 1e-15 {
            a_prime.push(DVec3::ZERO);
            w_prime.push(DVec3::ZERO);
        } else {
            let s = p / denom;
            a_prime.push(s * (ws[i + 1] * points[i + 1] - ws[i] * points[i]));
            w_prime.push(DVec3::new(s * (ws[i + 1] - ws[i]), 0.0, 0.0));
        }
    }

    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0f64; m];

    // A'(t): non-rational B-Spline of degree p-1
    let a_prime_t = de_boor(degree - 1, deriv_knots, &a_prime, &unit, t);
    // W'(t): scalar B-Spline of degree p-1 (embedded in .x)
    let w_prime_t = de_boor(degree - 1, deriv_knots, &w_prime, &unit, t).x;

    // W(t) and C(t) from the homogeneous evaluation
    let h = crate::math::bspl::de_boor_homo(degree, knots, points, weights, t);
    let w_t = h[3];
    if w_t.abs() < 1e-15 {
        return DVec3::ZERO;
    }
    let c_t = DVec3::new(h[0] / w_t, h[1] / w_t, h[2] / w_t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

pub(crate) fn bspline_tangent_analytic_2d(
    degree: usize,
    knots: &[f64],
    points: &[DVec2],
    weights: &[f64],
    t: f64,
) -> DVec2 {
    let n = points.len();
    if n < 2 || degree == 0 {
        return DVec2::ZERO;
    }

    // OCCT: a NULL weight array (non-rational curve) is treated as all weights
    // equal to 1.0 (BSplCLib weight accessor, cf. math/bspl.rs wgt).
    let ws: Vec<f64> = if weights.is_empty() { vec![1.0; n] } else { weights.to_vec() };

    let p = degree as f64;
    let m = n - 1;

    let mut a_prime = Vec::with_capacity(m);
    let mut w_prime = Vec::with_capacity(m);
    for i in 0..m {
        let denom = knots[i + degree + 1] - knots[i + 1];
        if denom.abs() < 1e-15 {
            a_prime.push(DVec2::ZERO);
            w_prime.push(DVec2::ZERO);
        } else {
            let s = p / denom;
            a_prime.push(s * (ws[i + 1] * points[i + 1] - ws[i] * points[i]));
            w_prime.push(DVec2::new(s * (ws[i + 1] - ws[i]), 0.0));
        }
    }

    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0; m];
    let a_prime_t = de_boor_2d(degree - 1, deriv_knots, &a_prime, &unit, t);
    let w_prime_t = de_boor_2d(degree - 1, deriv_knots, &w_prime, &unit, t).x;

    let h = de_boor_homo_2d(degree, knots, points, weights, t);
    let w_t = h[2];
    if w_t.abs() < 1e-15 {
        return DVec2::ZERO;
    }
    let c_t = DVec2::new(h[0] / w_t, h[1] / w_t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

/// Analytic tangent for a rational Bezier curve using the quotient rule.
///
/// The derivative of a degree-n Bezier is a degree-(n-1) Bezier with:
///   A'_i = n*(w_{i+1}*P_{i+1} - w_i*P_i)
///   W'_i = n*(w_{i+1} - w_i)
pub(crate) fn bezier_tangent_analytic(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n < 2 {
        return DVec3::ZERO;
    }
    let deg = (n - 1) as f64;

    let mut a_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        a_prime.push(deg * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
        w_prime.push(DVec3::new(deg * (weights[i + 1] - weights[i]), 0.0, 0.0));
    }

    let unit = vec![1.0f64; n - 1];
    let a_prime_t = de_casteljau_3d(&a_prime, &unit, t);
    let w_prime_t = de_casteljau_3d(&w_prime, &unit, t).x;

    // W(t): evaluate weights as scalar Bezier (embed in .x with unit weights)
    let w_pts: Vec<DVec3> = weights.iter().map(|&w| DVec3::new(w, 0.0, 0.0)).collect();
    let w_unit = vec![1.0f64; n]; // n elements to match w_pts
    let w_t = de_casteljau_3d(&w_pts, &w_unit, t).x;
    if w_t.abs() < 1e-15 {
        return DVec3::ZERO;
    }

    // C(t) from the standard rational evaluation
    let c_t = de_casteljau_3d(points, weights, t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

/// Analytic derivative for a rational Bezier curve in 2D.
/// Same formula as `bezier_tangent_analytic` but operating on DVec2.
pub(crate) fn bezier_tangent_analytic_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n < 2 {
        return DVec2::ZERO;
    }
    let deg = (n - 1) as f64;
    let mut a_prime: Vec<DVec2> = Vec::with_capacity(n - 1);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        a_prime.push(deg * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
        w_prime.push(DVec3::new(deg * (weights[i + 1] - weights[i]), 0.0, 0.0));
    }
    let unit = vec![1.0f64; n - 1];
    let a_prime_t = de_casteljau_2d(&a_prime, &unit, t);
    // Evaluate w'(t) — use DVec3 to embed w' scalar in .x
    let w_prime_t = de_casteljau_3d(&w_prime, &unit, t).x;
    let w_pts: Vec<DVec3> = weights.iter().map(|&w| DVec3::new(w, 0.0, 0.0)).collect();
    let w_unit = vec![1.0f64; n];
    let w_t = de_casteljau_3d(&w_pts, &w_unit, t).x;
    if w_t.abs() < 1e-15 {
        return DVec2::ZERO;
    }
    let c_t = de_casteljau_2d(points, weights, t);
    a_prime_t - (w_prime_t * c_t) / w_t
}

// --- Curve2dEval implementations ---

pub(crate) fn de_casteljau_3d(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 {
        return DVec3::ZERO;
    }
    // Work in homogeneous coordinates [x*w, y*w, z*w, w]
    let mut d: Vec<[f64; 4]> = points
        .iter()
        .zip(weights)
        .map(|(p, &w)| [p.x * w, p.y * w, p.z * w, w])
        .collect();
    for r in 1..n {
        for j in 0..n - r {
            let next = d[j + 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(next.iter()) {
                *elem = (1.0 - t) * *elem + t * p;
            }
        }
    }
    let w = d[0][3];
    if w.abs() < 1e-15 {
        DVec3::ZERO
    } else {
        DVec3::new(d[0][0] / w, d[0][1] / w, d[0][2] / w)
    }
}

/// De Casteljau algorithm for rational Bezier curve evaluation in 2D.
pub(crate) fn de_casteljau_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 {
        return DVec2::ZERO;
    }
    let mut d: Vec<[f64; 3]> = points
        .iter()
        .zip(weights)
        .map(|(p, &w)| [p.x * w, p.y * w, w])
        .collect();
    for r in 1..n {
        for j in 0..n - r {
            let next = d[j + 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(next.iter()) {
                *elem = (1.0 - t) * *elem + t * p;
            }
        }
    }
    let w = d[0][2];
    if w.abs() < 1e-15 {
        DVec2::ZERO
    } else {
        DVec2::new(d[0][0] / w, d[0][1] / w)
    }
}
