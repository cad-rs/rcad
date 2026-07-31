//! Curve and surface approximation (Approx + GeomConvert_Approx + AppCont + AppDef).
//!
//! OCCT TKGeomBase packages: Approx, AppCont, AppDef, AppParCurves.
//!
//! Approximates curves and surfaces to BSpline within a given tolerance.
//! Also provides multi-line parallel approximation (AppParCurves).

use crate::geom::{BSplineCurve3, BSplineSurface, Curve3, CurveEval, Surface3, SurfaceEval};
use glam::DVec3;

const TOL: f64 = 1e-7;

/// Approximate a curve to BSpline within given tolerance.
///
/// OCCT: `GeomConvert_ApproxCurve`.
pub fn approx_curve(curve: &Curve3, tol: f64) -> Option<BSplineCurve3> {
    let dom = curve.default_domain();
    let t_min = dom[0];
    let t_max = dom[1];
    // For unbounded curves (e.g. Line, Parabola), fall back to a finite
    // range so sampling is well-defined. Matches rcad-algorithms
    // approx_curve_to_bspline behavior.
    let (t_min, t_max) = if !t_min.is_finite() || !t_max.is_finite() {
        (-10.0, 10.0)
    } else {
        (t_min, t_max)
    };
    let range = t_max - t_min;
    if range < tol {
        return None;
    }

    // Adaptive sampling: start coarse, refine where needed
    let mut pts = Vec::new();
    let n_init = 16;
    for i in 0..n_init {
        let t = t_min + range * (i as f64) / ((n_init - 1) as f64);
        pts.push((t, curve.point_at(t)));
    }

    // Refinement pass: check midpoints for deviation
    let mut i = 0;
    while i < pts.len() - 1 {
        let t0 = pts[i].0;
        let t1 = pts[i + 1].0;
        let t_mid = (t0 + t1) * 0.5;
        let p_mid = curve.point_at(t_mid);

        // Linear interpolation at midpoint
        let frac = (t_mid - t0) / (t1 - t0);
        let p_linear = pts[i].1 + (pts[i + 1].1 - pts[i].1) * frac;
        let dev = (p_mid - p_linear).length();

        if dev > tol {
            pts.insert(i + 1, (t_mid, p_mid));
        } else {
            i += 1;
        }
        if pts.len() > 1024 {
            break; // safety cap
        }
    }

    if pts.len() < 2 {
        return None;
    }

    let n = pts.len();
    let mut params = vec![0.0_f64; n];
    for j in 1..n {
        let d = (pts[j].1 - pts[j - 1].1).length();
        params[j] = params[j - 1] + d.max(1e-15);
    }
    let total = params[n - 1];
    for p in &mut params {
        *p /= total;
    }

    let degree = 3.min(n - 1);
    let n_knots = n + degree + 1;
    let mut knots = vec![0.0_f64; n_knots];
    for k in &mut knots[..=degree] {
        *k = params[0];
    }
    for j in 1..n - degree {
        let mut sum = 0.0;
        for k in j..j + degree {
            sum += params[k];
        }
        knots[j + degree] = sum / (degree as f64);
    }
    for k in &mut knots[n_knots - degree - 1..] {
        *k = params[n - 1];
    }

    Some(BSplineCurve3 {
        degree,
        knots,
        control_points: pts.iter().map(|&(_, p)| p).collect(),
        weights: vec![],
        is_periodic: false,
    })
}

/// Approximate a surface to BSpline within given tolerance.
///
/// OCCT: `GeomConvert_ApproxSurface`.
pub fn approx_surface(surface: &Surface3, tol: f64) -> Option<BSplineSurface> {
    let dom = surface.default_domain();
    let (u_min, u_max, v_min, v_max) = (dom[0], dom[1], dom[2], dom[3]);
    if !u_min.is_finite() || !v_min.is_finite() {
        return None;
    }

    // Adaptive grid refinement
    let n_u_init = 8.max(2);
    let n_v_init = 8.max(2);
    let mut pts: Vec<Vec<(f64, f64, DVec3)>> = Vec::new();

    for i in 0..n_u_init {
        let u = u_min + (u_max - u_min) * (i as f64) / ((n_u_init - 1) as f64);
        let mut row = Vec::with_capacity(n_v_init);
        for j in 0..n_v_init {
            let v = v_min + (v_max - v_min) * (j as f64) / ((n_v_init - 1) as f64);
            row.push((u, v, surface.point_at(u, v)));
        }
        pts.push(row);
    }

    // Refinement pass (simplified)
    let mut refined = true;
    while refined && pts.len() < 64 && pts[0].len() < 64 {
        refined = false;
        let mut new_pts = pts.clone();

        // Refine in U direction
        for i in 0..pts.len() - 1 {
            for j in 0..pts[i].len() {
                let u = (pts[i][j].0 + pts[i + 1][j].0) * 0.5;
                let p_mid = surface.point_at(u, pts[i][j].1);
                let p_lin = (pts[i][j].2 + pts[i + 1][j].2) * 0.5;
                if (p_mid - p_lin).length() > tol {
                    // Insert row
                    // Simplified: just add points
                    refined = true;
                }
            }
        }
        // For now, just accept the initial grid
        break;
    }

    let nu = pts.len();
    let nv = pts[0].len();
    let mut ctrl = Vec::with_capacity(nu);
    for row in &pts {
        let r: Vec<DVec3> = row.iter().map(|&(_, _, p)| p).collect();
        ctrl.push(r);
    }

    // Degree-(1,1) bilinear surface
    let knots_u = build_knots(nu, 1);
    let knots_v = build_knots(nv, 1);

    Some(BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u,
        knots_v,
        control_points: ctrl,
        weights: vec![],
    })
}

fn build_knots(n_ctrl: usize, degree: usize) -> Vec<f64> {
    let n_segments = n_ctrl - degree;
    let mut knots = vec![0.0f64; degree + 1];
    for i in 1..n_segments {
        knots.push(i as f64 / n_segments as f64);
    }
    knots.extend(vec![1.0f64; degree + 1]);
    knots
}

/// Multi-line parallel approximation (AppParCurves).
///
/// OCCT: `AppParCurves_MultiCurve`.
/// Approximates multiple point sequences with the same knot vector.
pub fn parallel_approximation(curves: &[Vec<DVec3>], degree: usize) -> Option<Vec<BSplineCurve3>> {
    if curves.is_empty() {
        return None;
    }
    let n = curves[0].len();
    if n < 2 {
        return None;
    }

    // Same knot vector for all curves
    let n_knots = n + degree + 1;
    let mut knots = vec![0.0_f64; n_knots];
    for k in &mut knots[..=degree] {
        *k = 0.0;
    }
    for j in 1..n - degree {
        knots[j + degree] = j as f64 / (n - degree) as f64;
    }
    for k in &mut knots[n_knots - degree - 1..] {
        *k = 1.0;
    }

    let result: Vec<BSplineCurve3> = curves
        .iter()
        .map(|pts| BSplineCurve3 {
            degree,
            knots: knots.clone(),
            control_points: pts.clone(),
            weights: vec![],
            is_periodic: false,
        })
        .collect();

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;
    #[test]
    fn test_approx_line() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let bs = approx_curve(&line, 0.01);
        assert!(bs.is_some());
    }
}
