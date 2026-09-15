//! Per-kind evaluation impls for the bounded curve/surface kinds: OCCT
//! `Geom_BoundedCurve` / `Geom2d_BoundedCurve` (BSpline, Bezier in 3D and 2D),
//! `Geom_BoundedSurface` (BSpline / Bezier / triangular Bezier surfaces) and
//! the AHT / T-Bezier special 2D forms.
//!
//! Extracted verbatim from `geom/eval.rs` (project Rule 5 file-size split);
//! the evaluation traits stay re-exported from `geom/eval.rs`.
use crate::geom::*;
use super::{
    bezier_surface_d0, bezier_surface_d1, bezier_surface_d2, bezier_tangent_analytic,
    bezier_tangent_analytic_2d, de_casteljau_2d, de_casteljau_3d,
};

// --- BoundedCurveEval implementations ---

impl BoundedCurveEval for BSplineCurve3 {
    fn degree(&self) -> usize { self.degree }
}

impl BoundedCurveEval for BezierCurve3 {
    fn degree(&self) -> usize {
        self.control_points.len().saturating_sub(1)
    }
}

// --- BoundedCurve2dEval implementations ---

impl BoundedCurve2dEval for BSplineCurve2 {
    fn degree(&self) -> usize { self.degree }
}

impl BoundedCurve2dEval for BezierCurve2 {
    fn degree(&self) -> usize {
        self.control_points.len().saturating_sub(1)
    }
}

/// OCCT-aligned: `Geom_ElementarySurface` intermediate abstract class.
///
// --- BoundedSurfaceEval implementations ---

impl BoundedSurfaceEval for BSplineSurface {
    fn degree_u(&self) -> usize { self.degree_u }
    fn degree_v(&self) -> usize { self.degree_v }
}

impl BoundedSurfaceEval for BezierSurface {
    fn degree_u(&self) -> usize {
        self.control_points.len().saturating_sub(1)
    }
    fn degree_v(&self) -> usize {
        self.control_points.first().map_or(0, |r| r.len().saturating_sub(1))
    }
}

impl CurveEval for BSplineCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        crate::math::bspl::de_boor(
            self.degree, &self.knots, &self.control_points, &self.weights, t,
        )
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t,
        ).normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t,
        )
    }
    fn derivative2_at(&self, t: f64) -> DVec3 {
        let h = 1e-5;
        (crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t + h,
        ) - crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t - h,
        )) / (2.0 * h)
    }
    fn derivative3_at(&self, t: f64) -> DVec3 {
        let h = 1e-4;
        (self.derivative2_at(t + h) - self.derivative2_at(t - h)) / (2.0 * h)
    }
    fn default_domain(&self) -> [f64; 2] {
        let d = self.degree;
        let n = self.knots.len();
        if n < 2 * d + 2 { return [0.0, 1.0]; }
        [self.knots[d], self.knots[n - d - 1]]
    }
}

impl SurfaceEval for BSplineSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        // Tensor product rational evaluation (NURBS):
        // 1. For each v-column, evaluate the u-direction NURBS in homogeneous coords
        //    -> get (wx, wy, wz, w) for each column index.
        // 2. Collect column weights and weighted positions.
        // 3. Run de Boor in v on the homogeneous results, then divide by weight.
        let n_u = self.control_points.len();
        if n_u == 0 {
            return DVec3::ZERO;
        }
        let n_v = self.control_points[0].len();
        if n_v == 0 {
            return DVec3::ZERO;
        }
        // Step 1: evaluate each v-column in the u direction -> homogeneous 4-vector
        let col_homo: Vec<[f64; 4]> = (0..n_v)
            .map(|j| {
                let pts: Vec<DVec3> = (0..n_u).map(|i| self.control_points[i][j]).collect();
                let wts: Vec<f64> = (0..n_u).map(|i| self.weights[i][j]).collect();
                crate::math::bspl::de_boor_homo(self.degree_u, &self.knots_u, &pts, &wts, u)
            })
            .collect();
        // Step 2: build the v-direction "control points" and "weights" from col_homo
        let v_pts: Vec<DVec3> = col_homo
            .iter()
            .map(|h| {
                let w = h[3];
                if w.abs() < 1e-15 {
                    DVec3::ZERO
                } else {
                    DVec3::new(h[0] / w, h[1] / w, h[2] / w)
                }
            })
            .collect();
        let v_wts: Vec<f64> = col_homo.iter().map(|h| h[3]).collect();
        // Step 3: rational de Boor in v
        crate::math::bspl::de_boor(self.degree_v, &self.knots_v, &v_pts, &v_wts, v)
    }
    /// OCCT `Geom_BSplineSurface::EvalD1` (Geom_BSplineSurface_1.cxx L146-190)
    /// — `BSplSLib::D1` over the flat knot sequences (`Mults == NoMults`).
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let rational = self.is_rational_u() || self.is_rational_v();
        let weights: Option<&[Vec<f64>]> = if rational { Some(&self.weights) } else { None };
        crate::geom::eval_b::bspl_slib_d1(
            u,
            v,
            0,
            0,
            self.degree_u as i32,
            self.degree_v as i32,
            self.is_rational_u(),
            self.is_rational_v(),
            self.is_periodic_u,
            self.is_periodic_v,
            &self.control_points,
            weights,
            &self.knots_u,
            &self.knots_v,
        )
    }
    /// OCCT `Geom_BSplineSurface::EvalD2` (Geom_BSplineSurface_1.cxx L191-238)
    /// — `BSplSLib::D2`, in the rcad `SurfaceEval::derivatives2` order
    /// `(P, dP/du, dP/dv, d2P/du2, d2P/dudv, d2P/dv2)`.
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let rational = self.is_rational_u() || self.is_rational_v();
        let weights: Option<&[Vec<f64>]> = if rational { Some(&self.weights) } else { None };
        crate::geom::eval_b::bspl_slib_d2(
            u,
            v,
            0,
            0,
            self.degree_u as i32,
            self.degree_v as i32,
            self.is_rational_u(),
            self.is_rational_v(),
            self.is_periodic_u,
            self.is_periodic_v,
            &self.control_points,
            weights,
            &self.knots_u,
            &self.knots_v,
        )
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        // The unit normal from the analytic first partials (`dP/du ^ dP/dv`);
        // the previous one-sided finite difference is gone now that the
        // analytic `Geom_BSplineSurface::EvalD1` body lives here.
        let (_p, du, dv) = self.derivatives(u, v);
        let n = du.cross(dv);
        let len = n.length();
        if len < 1e-15 { DVec3::Z } else { n / len }
    }
    fn default_domain(&self) -> [f64; 4] {
        let du = self.degree_u;
        let dv = self.degree_v;
        let nu = self.knots_u.len();
        let nv = self.knots_v.len();
        let u0 = if nu > du { self.knots_u[du] } else { 0.0 };
        let u1 = if nu > du + 1 {
            self.knots_u[nu - du - 1]
        } else {
            1.0
        };
        let v0 = if nv > dv { self.knots_v[dv] } else { 0.0 };
        let v1 = if nv > dv + 1 {
            self.knots_v[nv - dv - 1]
        } else {
            1.0
        };
        [u0, u1, v0, v1]
    }
}

impl CurveEval for BezierCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        de_casteljau_3d(&self.control_points, &self.weights, t)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        bezier_tangent_analytic(&self.control_points, &self.weights, t).normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        bezier_tangent_analytic(&self.control_points, &self.weights, t)
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }
}

impl SurfaceEval for BezierSurface {
    /// OCCT `Geom_BezierSurface::EvalD0` (Geom_BezierSurface.cxx L1416-1466) —
    /// the rational tensor evaluation `P = N(u, v) / W(u, v)` of `BSplSLib::D0`
    /// (via `PrepareEval` + `BSplCLib::Eval`), NOT the per-column rational
    /// u-then-unit-weight-v scheme this method used before (that scheme equals
    /// the tensor evaluation only when every V-column shares one weight
    /// polynomial `W_j(u)`).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        bezier_surface_d0(self, u, v)
    }
    /// OCCT `Geom_BezierSurface::EvalD1` (Geom_BezierSurface.cxx L1470-1524).
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        bezier_surface_d1(self, u, v)
    }
    /// OCCT `Geom_BezierSurface::EvalD2` (Geom_BezierSurface.cxx L1528-1590).
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        bezier_surface_d2(self, u, v)
    }
    /// The unit normal from the analytic first partials (`dP/du ^ dP/dv`);
    /// OCCT has no `Normal()` on `Geom_BezierSurface` (callers of
    /// `Geom_Surface::D1` build the normal the same way).
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let (_p, du, dv) = self.derivatives(u, v);
        let n = du.cross(dv);
        let len = n.length();
        if len < 1e-15 { DVec3::Z } else { n / len }
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, v| acc * v as f64)
}

fn trinomial_coeff(n: usize, i: usize, j: usize, k: usize) -> f64 {
    factorial(n) / (factorial(i) * factorial(j) * factorial(k))
}

impl SurfaceEval for TriBezierSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let degree = self.control_points.len().saturating_sub(1);
        if self.control_points.is_empty() || self.weights.len() != self.control_points.len() {
            return DVec3::ZERO;
        }

        let w = 1.0 - u - v;
        let mut homo = [0.0; 4];
        for (i, row) in self.control_points.iter().enumerate() {
            if row.len() != degree + 1 - i
                || self.weights.get(i).map(|r| r.len()) != Some(row.len())
            {
                return DVec3::ZERO;
            }
            for (j, point) in row.iter().enumerate() {
                let k = degree - i - j;
                let basis = trinomial_coeff(degree, i, j, k)
                    * u.powi(i as i32)
                    * v.powi(j as i32)
                    * w.powi(k as i32);
                let weight = self.weights[i][j];
                homo[0] += basis * weight * point.x;
                homo[1] += basis * weight * point.y;
                homo[2] += basis * weight * point.z;
                homo[3] += basis * weight;
            }
        }

        if homo[3].abs() < 1e-15 {
            DVec3::ZERO
        } else {
            DVec3::new(homo[0] / homo[3], homo[1] / homo[3], homo[2] / homo[3])
        }
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

impl Curve2dEval for BezierCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        de_casteljau_2d(&self.control_points, &self.weights, t)
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        bezier_tangent_analytic_2d(&self.control_points, &self.weights, t)
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }
}

impl Curve2dEval for BSplineCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        if self.is_periodic {
            // OCCT Geom2d_BSplineCurve::Value -> D0 -> EvalD0
            // (Geom2d_BSplineCurve_1.cxx L175-201): the PeriodicNormalization
            // wrap and the periodic PoleIndex/BuildKnots arms of the BSplCLib
            // evaluation — the legacy de Boor stencil has no pole wrap, so a
            // periodic curve routes through the exact DN engine.
            return crate::geom::bspline2d_dn::eval_d0(self, t);
        }
        crate::math::bspl::de_boor_2d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_BSplineCurve::D1 -> EvalD1 (Geom2d_BSplineCurve_1.cxx
        // L199-226) — the exact BSplCLib evaluation (LocateParameter /
        // BuildKnots / Bohm / PLib::RationalDerivative); the rcad engine
        // lives in bspline2d_dn.
        crate::geom::bspline2d_dn::eval_d1(self, t).d1
    }
    fn derivative2_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_BSplineCurve::D2 -> EvalD2 (Geom2d_BSplineCurve_1.cxx
        // L232-264) — the exact BSplCLib evaluation (bspline2d_dn).
        crate::geom::bspline2d_dn::eval_d2(self, t).d2
    }
    fn derivative3_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_BSplineCurve::D3 -> EvalD3 (Geom2d_BSplineCurve_1.cxx
        // L266-299) — the exact BSplCLib evaluation (bspline2d_dn).
        crate::geom::bspline2d_dn::eval_d3(self, t).d3
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        let d = self.degree;
        let n = self.knots.len();
        if self.is_periodic {
            // OCCT FirstParameter/LastParameter of a periodic curve:
            // myFlatKnots.Value(myDeg + 1) / Value(Upper - myDeg) — on the
            // wrapped OCCT knot sequence these are Knots(1) / Knots(NbKnots);
            // the carrier stores the plain expansion whose first/last entries
            // are exactly those knots (the end multiplicities equal the
            // degree).
            return [self.knots[0], self.knots[n - 1]];
        }
        if n > 2 * d {
            [self.knots[d], self.knots[n - d - 1]]
        } else if n >= 2 {
            [self.knots[0], self.knots[n - 1]]
        } else {
            [0.0, 1.0]
        }
    }
    // OCCT Geom2d_BSplineCurve::IsPeriodic() — the carrier flag.
    fn is_periodic(&self) -> bool {
        self.is_periodic
    }
}

fn aht_basis_values(t: f64, alg_deg: usize, alpha: f64, beta: f64) -> Vec<f64> {
    // Basis: {1, t, ..., t^k, sinh(αt), cosh(αt), sin(βt), cos(βt)}
    let mut basis = Vec::new();
    // Polynomial part: 1, t, t^2, ..., t^k
    let mut tp = 1.0;
    for _ in 0..=alg_deg {
        basis.push(tp);
        tp *= t;
    }
    // Hyperbolic part: sinh(αt), cosh(αt)
    if alpha > 0.0 {
        let a = alpha * t;
        basis.push(a.sinh());
        basis.push(a.cosh());
    }
    // Trigonometric part: sin(βt), cos(βt)
    if beta > 0.0 {
        let b = beta * t;
        basis.push(b.sin());
        basis.push(b.cos());
    }
    basis
}

impl Curve2dEval for AHTBezierCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        let basis = aht_basis_values(t, self.alg_degree, self.alpha, self.beta);
        let n = self.control_points.len().min(basis.len());
        if self.weights.is_empty() {
            // Non-rational: straight sum
            let mut pt = DVec2::ZERO;
            for i in 0..n {
                pt += self.control_points[i] * basis[i];
            }
            pt
        } else {
            // Rational: weighted sum / weight sum
            let mut pt = DVec2::ZERO;
            let mut wsum = 0.0;
            for i in 0..n {
                let w = if i < self.weights.len() {
                    self.weights[i]
                } else {
                    1.0
                };
                pt += self.control_points[i] * (w * basis[i]);
                wsum += w * basis[i];
            }
            if wsum.abs() > 1e-15 { pt / wsum } else { pt }
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }
}

impl Curve2dEval for TBezierCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        // Basis: {1, cos(αt), sin(αt), cos(2αt), sin(2αt), ..., cos(n·αt), sin(n·αt)}
        let n = self.order;
        let at = self.alpha * t;
        let mut pt = DVec2::ZERO;
        let mut wsum = 0.0;
        let has_weights = !self.weights.is_empty();
        // Constant basis = 1
        let w0 = if has_weights { self.weights[0] } else { 1.0 };
        pt += self.control_points[0] * w0;
        wsum += w0;
        for i in 1..=n {
            let fi = i as f64;
            let c = (fi * at).cos();
            let s = (fi * at).sin();
            let idx_c = 2 * i - 1;
            let idx_s = 2 * i;
            if idx_c < self.control_points.len() {
                let wc = if has_weights && idx_c < self.weights.len() {
                    self.weights[idx_c]
                } else {
                    1.0
                };
                pt += self.control_points[idx_c] * (wc * c);
                wsum += wc * c;
            }
            if idx_s < self.control_points.len() {
                let ws = if has_weights && idx_s < self.weights.len() {
                    self.weights[idx_s]
                } else {
                    1.0
                };
                pt += self.control_points[idx_s] * (ws * s);
                wsum += ws * s;
            }
        }
        if has_weights && wsum.abs() > 1e-15 {
            pt / wsum
        } else {
            pt
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, std::f64::consts::PI / self.alpha]
    }
}
