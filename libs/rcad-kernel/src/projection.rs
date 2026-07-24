//! Closest-point projection from a 3D point onto a curve or surface.
//!
//! Analogous to OCCT `GeomAPI_ProjectPointOnCurve` and
//! `GeomAPI_ProjectPointOnSurf`.

/// Solve cubic x³ + a·x² + b·x + c = 0, returning real roots.
fn solve_cubic(a: f64, b: f64, c: f64) -> Vec<f64> {
    let a2 = a * a;
    let p = b - a2 / 3.0;
    let q = c + (2.0 * a2 * a - 9.0 * a * b) / 27.0;
    let disc = q * q / 4.0 + p * p * p / 27.0;
    let shift = -a / 3.0;
    if disc >= 0.0 {
        let s = disc.sqrt();
        let u = (-q / 2.0 + s).cbrt();
        let v = (-q / 2.0 - s).cbrt();
        vec![u + v + shift]
    } else {
        let r = (-p * p * p / 27.0).sqrt();
        let phi = (-q / (2.0 * r)).acos();
        let r_rt = 2.0 * r.cbrt();
        vec![
            r_rt * (phi / 3.0).cos() + shift,
            r_rt * ((phi + 2.0 * std::f64::consts::PI) / 3.0).cos() + shift,
            r_rt * ((phi + 4.0 * std::f64::consts::PI) / 3.0).cos() + shift,
        ]
    }
}

/// Solve quartic a·x⁴ + b·x³ + c·x² + d·x + e = 0.
/// OCCT-aligned: math_DirectPolynomialRoots.
fn solve_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Vec<f64> {
    if a.abs() < 1e-30 { return vec![]; }
    let inv_a = 1.0 / a;
    let ba = b * inv_a; let ca = c * inv_a; let da = d * inv_a; let ea = e * inv_a;
    let p = ca - 3.0 * ba * ba / 8.0;
    let q = da - ba * ca / 2.0 + ba * ba * ba / 8.0;
    let r = ea - ba * da / 4.0 + ba * ba * ca / 16.0 - 3.0 * ba * ba * ba * ba / 256.0;
    if q.abs() < 1e-30 {
        let disc = p * p - 4.0 * r;
        if disc < 0.0 { return vec![]; }
        let sd = disc.sqrt();
        let mut roots = Vec::new();
        for &t in &[(-p + sd) / 2.0, (-p - sd) / 2.0] {
            if t >= 0.0 { roots.push(t.sqrt()); }
            if t > 0.0 { roots.push(-t.sqrt()); }
        }
        let shift = -ba / 4.0;
        return roots.into_iter().map(|x| x + shift).collect();
    }
    let rc = solve_cubic(2.0 * p, p * p - 4.0 * r, -q * q);
    let m = rc.into_iter().find(|&m| m > 0.0).unwrap_or(0.0);
    if m <= 0.0 { return vec![]; }
    let sq = (m * 2.0).sqrt();
    let t1 = -p - m; let t2 = q / sq;
    let disc1 = -t1 - 2.0 * t2; let disc2 = -t1 + 2.0 * t2;
    let shift = -ba / 4.0;
    let mut roots = Vec::new();
    if disc1 >= 0.0 { let s = disc1.sqrt(); roots.push((sq + s) / 2.0 + shift); roots.push((-sq + s) / 2.0 + shift); }
    if disc2 >= 0.0 { let s = disc2.sqrt(); roots.push((s - sq) / 2.0 + shift); roots.push((-sq - s) / 2.0 + shift); }
    roots
}

use glam::DVec3;

use crate::geom::{Curve2d, Curve3, CurveEval, Surface3, SurfaceEval};

// ─────────────────────────────────────────────────────────────────────────────
// Result types
// ─────────────────────────────────────────────────────────────────────────────

/// Result of projecting a point onto a curve.
#[derive(Debug, Clone)]
pub struct CurveProjection {
    /// Nearest point on the curve.
    pub point: DVec3,
    /// Curve parameter at the nearest point.
    pub param: f64,
    /// Distance from the query point to the curve.
    pub distance: f64,
}

/// Result of projecting a point onto a surface.
#[derive(Debug, Clone)]
pub struct SurfaceProjection {
    /// Nearest point on the surface.
    pub point: DVec3,
    /// Surface parameter (u, v) at the nearest point.
    pub params: (f64, f64),
    /// Distance from the query point to the surface.
    pub distance: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Curve projection
// ─────────────────────────────────────────────────────────────────────────────

/// Newton refinement helper (OCCT-aligned: interval-perform Newton step).
/// Solves g(t) = (P - C(t))·C'(t) = 0 via newton:
///   t_{n+1} = t_n - (P-C)·C' / (|C'|² + (P-C)·C'')
/// where C'' is approximated by finite-difference.
fn newton_refine(
    curve: &Curve3,
    t: &mut f64,
    best_dist: &mut f64,
    query: DVec3,
    max_iter: usize,
    clamp: impl Fn(f64) -> f64,
) {
    let dt = 1e-7;
    for _ in 0..max_iter {
        let p = curve.point_at(*t);
        let diff = p - query;
        let deriv = curve.derivative_at(*t);
        let deriv_sq = deriv.dot(deriv);
        if deriv_sq < 1e-20 { break; }
        let curv = (curve.point_at(*t + 2.0 * dt) - 2.0 * p
            + curve.point_at(*t - 2.0 * dt)) / (dt * dt);
        let denom = deriv_sq + diff.dot(curv);
        let delta = diff.dot(deriv) / if denom.abs() > 1e-20 { denom } else { deriv_sq };
        let new_t = clamp(*t - delta);
        let new_dist = (curve.point_at(new_t) - query).length();
        if new_dist < *best_dist {
            *best_dist = new_dist;
            *t = new_t;
        }
        if delta.abs() < 1e-10 { break; }
    }
}

/// Project the point `query` onto `curve`, returning the nearest point on the
/// curve, its parameter value, and the Euclidean distance.
///
/// ✅ OCCT-aligned: dispatches per-type matching Extrema_ExtPC:
///   - Line/Circle: analytic via Extrema_ExtPElC equivalent
///   - Ellipse: analytic init + Newton refinement
///   - BSpline: C2 interval splitting (Extrema_GGExtPC L190-388)
///   - Bezier/Other: uniform sampling + Newton
///
/// `n_samples` is the uniform sampling count (used for Bezier and fallback;
/// for BSpline it's overridden by `degree + 1` per C2 interval).
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Circle3};
/// use rcad_kernel::projection::closest_point_on_curve;
///
/// let circle = Curve3::Circle(Circle3 {
///     center: DVec3::ZERO,
///     normal: DVec3::Z,
///     radius: 1.0,
/// });
/// let q = DVec3::new(2.0, 0.0, 0.0);
/// let result = closest_point_on_curve(&circle, q, 64);
/// assert!((result.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
/// ```
pub fn closest_point_on_curve(curve: &Curve3, query: DVec3, n_samples: usize) -> CurveProjection {
    // ── Analytic fast paths ───────────────────────────────────────────────────
    // Analogous to OCCT `ExtremaPC` per-type dispatch.

    match curve {
        Curve3::Line(l) => {
            // Closest point on an infinite line: project query onto line direction.
            let dir_sq = l.direction.dot(l.direction);
            if dir_sq < 1e-20 {
                // Degenerate line — return origin.
                let pt = l.origin;
                return CurveProjection {
                    point: pt,
                    param: 0.0,
                    distance: (pt - query).length(),
                };
            }
            let t = (query - l.origin).dot(l.direction) / dir_sq;
            let [t0, t1] = curve.default_domain();
            let t_clamped = if t0.is_finite() && t1.is_finite() {
                t.clamp(t0, t1)
            } else {
                t
            };
            let pt = l.origin + t_clamped * l.direction;
            return CurveProjection {
                point: pt,
                param: t_clamped,
                distance: (pt - query).length(),
            };
        }

        Curve3::Circle(circ) => {
            // ✅ OCCT-aligned: Extrema_ExtPElC::Perform(Circle) L92-190.
            // 1. Project query onto circle plane (subtract axial component)
            let o = circ.center;
            let axis = circ.normal.normalize_or_zero();
            let axial = (query - o).dot(axis);
            let pp = query - axis * axial;  // Pp = P projected to circle plane
            let opp = pp - o;               // OPp vector from center to projected point
            let opp_mag = opp.length();
            if opp_mag < 1e-15 {
                // P on axis → infinite solutions; return center(0) (OCCT returns IsDone=false)
                let pt = circ.point_at(0.0);
                return CurveProjection { point: pt, param: 0.0, distance: (pt - query).length() };
            }
            // Circle axes (OCCT: C.XAxis().Direction() and C.YAxis().Direction())
            let cx = circ.x_dir.normalize();
            let cy = circ.y_dir.normalize();
            let x = opp.dot(cx);
            let y = opp.dot(cy);
            // OCCT L138: Usol[0] = XAxis.AngleWithRef(OPp, Axis) → atan2 of (X×OPp)·Axis, X·OPp
            let u_min = y.atan2(x);                          // angle of OPp (minimum distance)
            let half = std::f64::consts::PI;
            let u_max = if u_min < 0.0 { u_min + half } else { u_min - half };  // antipodal (maximum)
            let [t0, t1] = curve.default_domain();
            // Adjust each solution into [t0, t0+TAU) then check against [t0, t1]
            let mut best_t = u_min;
            let mut best_d = f64::INFINITY;
            for &u in &[u_min, u_max] {
                let mut u_adj = u;
                // Shift into [t0, t0+TAU) — OCCT ElCLib::AdjustPeriodic
                if t0.is_finite() {
                    let period = std::f64::consts::TAU;
                    let diff = u_adj - t0;
                    u_adj = t0 + diff - period * (diff / period).floor();
                    if u_adj < t0 - 1e-12 { u_adj += period; }
                    if u_adj > t0 + period + 1e-12 { u_adj -= period; }
                }
                if u_adj >= t0 - 1e-12 && u_adj <= t1 + 1e-12 {
                    let pt = circ.point_at(u_adj);
                    let d = (pt - query).length();
                    if d < best_d { best_d = d; best_t = u_adj; }
                }
            }
            // If no solution in range, try endpoints + fallback to u_min clamped
            if best_d.is_infinite() {
                best_t = u_min.clamp(t0, t1);
                best_d = (circ.point_at(best_t) - query).length();
            }
            let pt = circ.point_at(best_t);
            return CurveProjection { point: pt, param: best_t, distance: (pt - query).length() };
        }

        Curve3::Ellipse(ell) => {
            // ✅ OCCT-aligned: Extrema_ExtPElC::Perform(Ellipse) L203-281.
            // Solve g(u) = (P-C)·C' = 0 using sign-change detection + Newton.
            // g(u) = (B²-A²)/2 * sin(2u) + A*X*sin(u) - B*Y*cos(u)
            // where A=major_radius, B=minor_radius, X=OPp·XAxis, Y=OPp·YAxis
            let o = ell.center;
            let axis = ell.normal.normalize_or_zero();
            let axial = (query - o).dot(axis);
            let pp = query - axis * axial;
            let opp = pp - o;
            if opp.length_squared() < 1e-30 && (ell.major_radius - ell.minor_radius).abs() < 1e-15 {
                // Point at center of circular ellipse → infinite solutions
                let pt = ell.point_at(0.0);
                return CurveProjection { point: pt, param: 0.0, distance: (pt - query).length() };
            }
            let ex = ell.major_dir.normalize();
            let ey = axis.cross(ex).normalize();
            let x = opp.dot(ex);
            let y = opp.dot(ey);
            let a = ell.major_radius;
            let b = ell.minor_radius;
            let [t0, t1] = curve.default_domain();
            let n = 33_usize;
            // Sample g(u) to find sign-change intervals
            let g = |u: f64| {
                let (s, c) = u.sin_cos();
                let sin2u = 2.0 * s * c;
                (b * b - a * a) / 2.0 * sin2u + a * x * s - b * y * c
            };
            let mut u_best = t0;
            let mut d_best = f64::INFINITY;
            for i in 0..=n {
                let u = t0 + (t1 - t0) * i as f64 / n as f64;
                let pt = ell.point_at(u);
                let d = (pt - query).length();
                if d < d_best { d_best = d; u_best = u; }
            }
            // Sign-change-based refinement (OCCT: TrigonometricFunctionRoots per-interval)
            let clamp = |u: f64| u.clamp(t0, t1);
            let mut g_prev = g(t0);
            for i in 1..=n {
                let u = t0 + (t1 - t0) * i as f64 / n as f64;
                let g_cur = g(u);
                if g_prev * g_cur <= 0.0 || g_prev.abs() < 1e-12 || g_cur.abs() < 1e-12 {
                    let u_mid = (u + (t0 + (t1 - t0) * (i - 1) as f64 / n as f64)) * 0.5;
                    let mut t = clamp(u_mid);
                    let mut dist = (ell.point_at(t) - query).length();
                    newton_refine(curve, &mut t, &mut dist, query, 20, clamp);
                    if dist < d_best { d_best = dist; u_best = t; }
                }
                g_prev = g_cur;
            }
            // Fallback: endpoints
            for &u in &[t0, t1] {
                let d = (ell.point_at(u) - query).length();
                if d < d_best { d_best = d; u_best = u; }
            }
            let pt = ell.point_at(u_best);
            return CurveProjection { point: pt, param: u_best, distance: (pt - query).length() };
        }

        // OCCT-aligned: Extrema_ExtPElC::Perform(Hyperbola) L294-389.
        //   Project onto hyperbola plane, solve quartic for v = e^u:
        //     C1·v⁴ + C2·v³ + 0·v² + C3·v - C1 = 0
        //   where C1=(R²+r²)/4, C2=-(X·R+Y·r)/2, C3=(X·R-Y·r)/2.
        Curve3::Hyperbola(hyp) => {
            let o = hyp.center;
            let axis = hyp.normal.normalize_or_zero();
            let pp = query - axis * (query - o).dot(axis);
            let opp = pp - o;
            let hx = hyp.major_dir.normalize();
            let hy = axis.cross(hx).normalize();
            let x = opp.dot(hx);
            let y = opp.dot(hy);
            let r = hyp.semi_major;
            let r2 = hyp.semi_minor;
            let c1 = (r * r + r2 * r2) / 4.0;
            let c2 = -(x * r + y * r2) / 2.0;
            let c3 = (x * r - y * r2) / 2.0;
            let roots = solve_quartic(c1, c2, 0.0, c3, -c1);
            let [t0, t1] = curve.default_domain();
            let mut u_best = t0;
            let mut d_best = f64::INFINITY;
            for v in roots {
                if v > 0.0 { let u = v.ln();
                    if u >= t0 && u <= t1 { let pt = hyp.point_at(u); let d = (pt - query).length(); if d < d_best { d_best = d; u_best = u; } }
                }
            }
            for &u in &[t0, t1] { let d = (hyp.point_at(u) - query).length(); if d < d_best { d_best = d; u_best = u; } }
            let pt = hyp.point_at(u_best);
            return CurveProjection { point: pt, param: u_best, distance: (pt - query).length() };
        }

        // OCCT-aligned: Extrema_ExtPElC::Perform(Parabola) L402-480.
        //   Project onto parabola plane, solve cubic for u:
        //     (1/(4·F))·u³ + (2·F-X)·u - 2·F·Y = 0
        Curve3::Parabola(par) => {
            let o = par.vertex;
            let axis = par.normal.normalize_or_zero();
            let pp = query - axis * (query - o).dot(axis);
            let opp = pp - o;
            let px = par.axis_dir.normalize();
            let py = axis.cross(px).normalize();
            let x = opp.dot(px);
            let y = opp.dot(py);
            let f = par.focal_param;
            let coeff = 1.0 / (4.0 * f);
            let roots = solve_cubic(0.0, (2.0 * f - x) / coeff, -2.0 * f * y / coeff);
            let [t0, t1] = curve.default_domain();
            let mut u_best = t0;
            let mut d_best = f64::INFINITY;
            for u in roots { if u >= t0 && u <= t1 { let pt = par.point_at(u); let d = (pt - query).length(); if d < d_best { d_best = d; u_best = u; } } }
            for &u in &[t0, t1] { let d = (par.point_at(u) - query).length(); if d < d_best { d_best = d; u_best = u; } }
            let pt = par.point_at(u_best);
            return CurveProjection { point: pt, param: u_best, distance: (pt - query).length() };
        }

        // OCCT-aligned: BSpline — C2 interval splitting,
        //   per-interval sampling + Newton refinement
        //   (Extrema_GGExtPC L190-388: knot interval subdivision)
        Curve3::BSpline(bs) => {
            let [t0, t1] = curve.default_domain();
            let intervals = bs.c2_intervals();
            let n_per_int = (bs.degree + 1).max(4);
            let mut best_t = t0;
            let mut best_dist = f64::INFINITY;

            let clamp_t = |t: f64| t.clamp(t0, t1);

            for wi in 0..intervals.len().saturating_sub(1) {
                let lo = intervals[wi];
                let hi = intervals[wi + 1];
                let span = hi - lo;
                if span < 1e-15 { continue; }

                // Sample n_per_int points in this interval
                for i in 0..=n_per_int {
                    let t = lo + span * i as f64 / n_per_int as f64;
                    let p = bs.point_at(t);
                    let d = (p - query).length();
                    if d < best_dist {
                        best_dist = d;
                        best_t = t;
                    }
                }
            }

            // Refine best candidate with Newton (finite-difference curvature)
            if best_dist < f64::INFINITY {
                newton_refine(curve, &mut best_t, &mut best_dist, query, 30, clamp_t);
            }

            let best_point = curve.point_at(best_t);
            return CurveProjection {
                point: best_point,
                param: best_t,
                distance: (best_point - query).length(),
            };
        }

        // Bezier: no C2 splitting needed (analytic)
        Curve3::Bezier(_) => {
            let [t0, t1] = curve.default_domain();
            let mut best_t = t0;
            let mut best_dist = f64::INFINITY;
            let n = n_samples.max(4);
            let clamp_t = |t: f64| t.clamp(t0, t1);
            for i in 0..=n {
                let t = t0 + (t1 - t0) * i as f64 / n as f64;
                let p = curve.point_at(t);
                let d = (p - query).length();
                if d < best_dist {
                    best_dist = d;
                    best_t = t;
                }
            }
            newton_refine(curve, &mut best_t, &mut best_dist, query, 30, clamp_t);
            let best_point = curve.point_at(best_t);
            return CurveProjection {
                point: best_point,
                param: best_t,
                distance: (best_point - query).length(),
            };
        }

        _ => {}
    }

    // ── Numerical fallback for all other curve types ───────────────────────────
    let [t0_raw, t1_raw] = curve.default_domain();
    let n = n_samples.max(4);

    // For infinite domains (lines), use a heuristic finite sampling range
    // centered on the closest parameter analytically (dot product for lines).
    let (t0, t1) = if t0_raw.is_infinite() || t1_raw.is_infinite() {
        // Use the analytical projection for the domain center estimate
        let t_center = 0.0; // non-line types; 0 is a safe fallback
        let span = 100.0_f64; // generous range around t_center
        (t_center - span, t_center + span)
    } else {
        (t0_raw, t1_raw)
    };

    // Step 1: coarse sampling
    let (mut best_t, mut best_dist) = {
        let mut bd = f64::INFINITY;
        let mut bt = t0;
        for i in 0..=n {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            let p = curve.point_at(t);
            let d = (p - query).length();
            if d < bd {
                bd = d;
                bt = t;
            }
        }
        (bt, bd)
    };

    // Step 2: Newton refinement
    let clamp_t = |t: f64| {
        if t0_raw.is_infinite() || t1_raw.is_infinite() {
            t
        } else {
            t.clamp(t0, t1)
        }
    };
    newton_refine(curve, &mut best_t, &mut best_dist, query, 30, clamp_t);

    let best_point = curve.point_at(best_t);
    CurveProjection {
        point: best_point,
        param: best_t,
        distance: (best_point - query).length(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface projection
// ─────────────────────────────────────────────────────────────────────────────

/// Project the point `query` onto `surface`, returning the nearest point on the
/// surface, its (u, v) parameters, and the Euclidean distance.
///
/// Analytic surfaces (Plane, Sphere, Cylinder, Cone, Torus) use closed-form
/// formulae.  All other surfaces fall back to numerical sampling + Newton
/// refinement.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::any_perpendicular;
/// use rcad_kernel::geom::{Surface3, SphericalSurface};
/// use rcad_kernel::projection::closest_point_on_surface;
///
/// let sphere = Surface3::Sphere(SphericalSurface {
///     center: DVec3::ZERO,
///     axis: DVec3::Z,
///     radius: 1.0,
///     ref_dir: any_perpendicular(DVec3::Z),
/// });
/// let q = DVec3::new(3.0, 0.0, 0.0);
/// let result = closest_point_on_surface(&sphere, q, 16);
/// assert!((result.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
/// ```
/// OCCT-aligned: closest_point_on_curve with range restriction [t0, t1].
/// Calls closest_point_on_curve then clamps parameter to [t0, t1] and
/// re-evaluates at the endpoint if the original projection fell outside range.
/// Matches OCCT GeomAPI_ProjectPointOnCurve::Init(curve, f, l).
pub fn closest_point_on_curve_range(curve: &Curve3, query: DVec3, t0: f64, t1: f64, n_samples: usize) -> CurveProjection {
    let mut result = closest_point_on_curve(curve, query, n_samples);
    if result.param < t0 || result.param > t1 {
        let (t_min, t_max) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        let t_clamped = if result.param < t_min { t_min } else { t_max };
        let pt = curve.point_at(t_clamped);
        result = CurveProjection {
            point: pt,
            param: t_clamped,
            distance: (pt - query).length(),
        };
    }
    result
}

pub fn closest_point_on_surface(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    use crate::geom::*;

    match surface {
        // ── Analytic closed-form projections ──────────────────────────────────
        Surface3::Plane(plane) => {
            // Project onto the infinite plane, no clamping needed.
            let d = (query - plane.origin).dot(plane.normal);
            let point = query - plane.normal * d;
            let distance = d.abs();
            // (u, v) = coordinates in plane (not strictly needed, but provide them)
            let u_axis = plane.u_dir;
            let v_axis = plane.v_dir;
            let diff = point - plane.origin;
            SurfaceProjection {
                point,
                params: (diff.dot(u_axis), diff.dot(v_axis)),
                distance,
            }
        }

        Surface3::Sphere(sph) => {
            let v = query - sph.center;
            let len = v.length();
            let point = if len < 1e-14 {
                sph.center + sph.radius * DVec3::X // degenerate: pick arbitrary
            } else {
                sph.center + v / len * sph.radius
            };
            // Compute (theta, phi) from point relative to center
            let w = (point - sph.center).normalize_or_zero();
            let u_axis = sph.ref_dir;
            let v_axis = sph.axis.cross(u_axis);
            let theta = w.dot(sph.axis).clamp(-1.0, 1.0).acos();
            let phi = w.dot(v_axis).atan2(w.dot(u_axis));
            SurfaceProjection {
                point,
                params: (phi, theta),
                distance: (point - query).length(),
            }
        }

        Surface3::Cylinder(cyl) => {
            // Project by collapsing along axis, then normalizing radial component.
            let v = query - cyl.origin;
            let along = v.dot(cyl.axis);
            let radial = v - cyl.axis * along;
            let radial_len = radial.length();
            let point = if radial_len < 1e-14 {
                cyl.origin + cyl.axis * along + cyl.radius * cyl.ref_dir
            } else {
                cyl.origin + cyl.axis * along + radial / radial_len * cyl.radius
            };
            let u_axis = cyl.ref_dir.normalize();
            let v_axis = cyl.axis.cross(u_axis);
            let r = (point - cyl.origin - cyl.axis * along).normalize_or_zero();
            let theta = r.dot(v_axis).atan2(r.dot(u_axis));
            // Map [-π, π] → [0, 2π] to match the canonical cylinder UV domain.
            let theta = if theta < 0.0 { theta + std::f64::consts::TAU } else { theta };
            SurfaceProjection {
                point,
                params: (theta, along),
                distance: (point - query).length(),
            }
        }

        Surface3::Cone(cone) => {
            // Project onto the cone's reference-circle parameterization.
            let axis = cone.axis_dir();
            let x_axis = any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let local = query - cone.apex;
            let along = local.dot(axis);
            let radial = local - axis * along;
            let radial_len = radial.length();
            let half = cone.half_angle_rad;
            let tan_h = half.tan();
            let axial = (along + (radial_len - cone.radius) * tan_h) / (1.0 + tan_h * tan_h);
            let r_hat = if radial_len < 1e-14 {
                x_axis
            } else {
                radial / radial_len
            };
            let point = cone.apex + axis * axial + r_hat * cone.radius_at_axial(axial);
            let slant = cone.slant_from_axial(axial);
            let theta = r_hat.dot(y_axis).atan2(r_hat.dot(x_axis));
            SurfaceProjection {
                point,
                params: (theta, slant),
                distance: (point - query).length(),
            }
        }

        Surface3::Torus(torus) => {
            // Step 1: project onto the major-radius circle in the equatorial plane.
            let v = query - torus.center;
            let along = v.dot(torus.axis);
            let radial = v - torus.axis * along;
            let radial_len = radial.length();
            // Compute a reference direction in the equatorial plane (like any_perpendicular)
            let torus_ref_dir = {
                let ref_dir = if torus.axis.x.abs() > 1.0 - 1e-12 { DVec3::Z } else { DVec3::X };
                (ref_dir - torus.axis * ref_dir.dot(torus.axis)).normalize_or_zero()
            };
            let major_dir = if radial_len < 1e-14 {
                torus_ref_dir
            } else {
                radial / radial_len
            };
            let tube_center = torus.center + major_dir * torus.major_radius;
            // Step 2: project onto the tube circle.
            let w = query - tube_center;
            let w_len = w.length();
            let point = if w_len < 1e-14 {
                tube_center + major_dir * torus.minor_radius
            } else {
                tube_center + w / w_len * torus.minor_radius
            };
            let u = major_dir
                .dot(torus_ref_dir)
                .atan2(major_dir.dot(torus.axis.cross(torus_ref_dir)));
            let w_dir = (point - tube_center).normalize_or_zero();
            let v_param = w_dir.dot(torus.axis).atan2(w_dir.dot(major_dir));
            SurfaceProjection {
                point,
                params: (u, v_param),
                distance: (point - query).length(),
            }
        }

        // ── Planar BSpline: use analytic plane projection ──────────────────────
        Surface3::BSpline(bsp) if bsp.degree_u == 1 && bsp.degree_v == 1
            && bsp.control_points.len() >= 2 && bsp.control_points[0].len() >= 2 =>
        {
            // The degree-1 BSpline from plane_to_bspline is geometrically a plane.
            // Use the analytic Plane projection instead of numeric_surface_projection
            // which may converge to the wrong UV when the query point is far from
            // the initial grid sample.
            let p00 = bsp.control_points[0][0];
            let p10 = bsp.control_points[1][0];
            let p01 = bsp.control_points[0][1];
            let du = p10 - p00;
            let dv = p01 - p00;
            let normal = du.cross(dv).normalize_or_zero();
            if normal.length_squared() > 0.5 {
                let n = normal;
                let d = (query - p00).dot(n);
                let point = query - n * d;
                let diff = point - p00;
                // Project onto the bilinear patch's local axes
                let u_len2 = du.length_squared();
                let v_len2 = dv.length_squared();
                let (u, v) = if u_len2 > 1e-30 && v_len2 > 1e-30 {
                    (diff.dot(du) / u_len2, diff.dot(dv) / v_len2)
                } else {
                    (0.0, 0.0)
                };
                SurfaceProjection { point, params: (u, v), distance: d.abs() }
            } else {
                numeric_surface_projection(surface, query, n_samples)
            }
        }

        // ── Numerical fallback for parametric surfaces ─────────────────────────
        _ => numeric_surface_projection(surface, query, n_samples),
    }
}

/// Numerical closest-point on a parametric surface via uniform sampling +
/// Newton refinement of `f(u,v) = |S(u,v) - Q|²`.
/// ✅ OCCT-aligned: TreatSolution (Extrema_ExtPS L97-135).
///   Wraps periodic UV into [uinf, uinf+period), then checks bounds.
fn treat_solution(
    u: f64, v: f64,
    point: DVec3,
    uinf: f64, usup: f64,
    vinf: f64, vsup: f64,
    u_period: Option<f64>,
    v_period: Option<f64>,
    tolu: f64, tolv: f64,
) -> Option<(f64, f64, DVec3)> {
    let mut u = u;
    let mut v = v;
    if let Some(per) = u_period {
        let diff = u - uinf;
        u = uinf + diff - per * (diff / per).floor();
        if u > usup + tolu { u -= per; }
        if u < uinf - tolu { u += per; }
    }
    if let Some(per) = v_period {
        let diff = v - vinf;
        v = vinf + diff - per * (diff / per).floor();
        if v > vsup + tolv { v -= per; }
        if v < vinf - tolv { v += per; }
    }
    if (uinf - u) <= tolu && (u - usup) <= tolu
        && (vinf - v) <= tolv && (v - vsup) <= tolv
    {
        Some((u, v, point))
    } else {
        None
    }
}

/// ✅ OCCT-aligned: Extrema_GenLocateExtPS (local Newton from initial UV guess).
///   Solves (S-Q)·dS/du = 0  and  (S-Q)·dS/dv = 0 via Gauss-Newton.
///   Returns the projected point and UV, or the initial guess if Newton fails.
pub fn closest_point_on_surface_near(
    surface: &Surface3,
    query: DVec3,
    u0: f64,
    v0: f64,
) -> SurfaceProjection {
    use crate::geom::SurfaceEval;
    let [uinf, usup, vinf, vsup] = surface.default_domain();
    let mut u = u0.clamp(uinf, usup);
    let mut v = v0.clamp(vinf, vsup);
    let mut best_dist = (surface.point_at(u, v) - query).length();
    for _ in 0..30 {
        let (p, du, dv) = surface.derivatives(u, v);
        let diff = p - query;
        let gu = diff.dot(du);
        let gv = diff.dot(dv);
        let huu = du.dot(du);
        let hvv = dv.dot(dv);
        let huv = du.dot(dv);
        let det = huu * hvv - huv * huv;
        if det.abs() < 1e-20 { break; }
        let du_step = (hvv * gu - huv * gv) / det;
        let dv_step = (huu * gv - huv * gu) / det;
        let nu = (u - du_step).clamp(uinf, usup);
        let nv = (v - dv_step).clamp(vinf, vsup);
        let nd = (surface.point_at(nu, nv) - query).length();
        if nd < best_dist { best_dist = nd; u = nu; v = nv; }
        if du_step.abs() < 1e-10 && dv_step.abs() < 1e-10 { break; }
    }
    let point = surface.point_at(u, v);
    SurfaceProjection { point, params: (u, v), distance: (point - query).length() }
}
///   Initialize: grid size = 32 (44 for BSpline, 300 when isoparametric curves degenerate).
///   Perform: grid sampling → Newton refinement with analytic derivatives.
///   OCCT L237-259: grid sizing + IsoIsDeg degenerate isoparametric detection.
fn numeric_surface_projection(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    use crate::geom::SurfaceEval;
    let [u0, u1, v0, v1] = surface.default_domain();

    // OCCT L237-259: grid size depends on surface type and isoparametric degeneracy
    let is_bspline = matches!(surface, Surface3::BSpline(_));
    // OCCT: nbU = isB ? 44 : 32, nbV = isB ? 44 : 32
    let mut nu = if is_bspline { 44usize } else { 32usize }.max(n_samples);
    let mut nv = if is_bspline { 44usize } else { 32usize }.max(n_samples);

    // OCCT L244-257: IsoIsDeg — check if isoparametric curves at domain boundaries
    //   have near-zero or near-infinite first derivative → degenerate → 300 samples
    let step_u = (u1 - u0) / 10.0;
    let step_v = (v1 - v0) / 10.0;
    if step_u.is_finite() && u1 > u0 {
        let d_max = (0..=10).filter_map(|i| {
            let u = u0 + step_u * i as f64;
            if u < u0 || u > u1 { return None; }
            let (_p, du, _dv) = surface.derivatives(u, v0);
            Some(du.length_squared())
        }).fold(0.0_f64, f64::max);
        if d_max < 1e-18 || d_max > 1e9 { nu = nu.max(300); }
    }
    if step_v.is_finite() && v1 > v0 {
        let d_max = (0..=10).filter_map(|i| {
            let v = v0 + step_v * i as f64;
            if v < v0 || v > v1 { return None; }
            let (_p, _du, dv) = surface.derivatives(u0, v);
            Some(dv.length_squared())
        }).fold(0.0_f64, f64::max);
        if d_max < 1e-18 || d_max > 1e9 { nv = nv.max(300); }
    }

    // Coarse sampling (OCCT L261: myExtPS.Initialize with nbU, nbV)
    let (mut best_u, mut best_v, mut best_dist) = {
        let mut bd = f64::INFINITY;
        let (mut bu, mut bv) = (u0, v0);
        for i in 0..=nu {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            for j in 0..=nv {
                let v = v0 + (v1 - v0) * j as f64 / nv as f64;
                let p = surface.point_at(u, v);
                let d = (p - query).length_squared();
                if d < bd { bd = d; bu = u; bv = v; }
            }
        }
        (bu, bv, bd.sqrt())
    };

    // OCCT-aligned: Newton refinement using analytic derivatives (SurfaceEval::derivatives).
    // Solves g(u,v) = (S(u,v)-Q) · dS/du = 0 and h(u,v) = (S(u,v)-Q) · dS/dv = 0.
    for _ in 0..40 {
        let (p, du, dv) = surface.derivatives(best_u, best_v);
        let diff = p - query;
        let gu = diff.dot(du);
        let gv = diff.dot(dv);
        let huu = du.dot(du);
        let hvv = dv.dot(dv);
        let huv = du.dot(dv);
        let det = huu * hvv - huv * huv;
        if det.abs() < 1e-20 { break; }
        let delta_u = (hvv * gu - huv * gv) / det;
        let delta_v = (huu * gv - huv * gu) / det;
        let new_u = (best_u - delta_u).clamp(u0, u1);
        let new_v = (best_v - delta_v).clamp(v0, v1);
        let new_dist = (surface.point_at(new_u, new_v) - query).length();
        if new_dist < best_dist {
            best_dist = new_dist;
            best_u = new_u;
            best_v = new_v;
        }
        if delta_u.abs() < 1e-10 && delta_v.abs() < 1e-10 { break; }
    }

    let best_point = surface.point_at(best_u, best_v);
    SurfaceProjection {
        point: best_point,
        params: (best_u, best_v),
        distance: (best_point - query).length(),
    }
}

/// OCCT-aligned: IntTools_Curve::MakePCurveOnSurface.
///
/// Samples the 3D curve at `n_samples` uniformly-spaced parameter values
/// within `t_range`, projects each onto `surface` via [`closest_point_on_surface`],
/// and fits a [`Curve2d::BSpline`] through the UV points.
///
/// Returns `None` when fewer than 2 valid projections are obtained, or
/// when the curve fit fails (coincident UV points, degenerate sample set).
pub fn make_pcurve_on_surface(
    curve: &Curve3,
    t_range: [f64; 2],
    surface: &Surface3,
    n_samples: usize,
) -> Option<Curve2d> {
    use crate::fit::interpolate_points_2d;
    use crate::geom::CurveEval;
    use glam::DVec2;

    let n = n_samples.max(2);
    let [t0, t1] = t_range;
    let dt = if n > 1 { (t1 - t0) / (n - 1) as f64 } else { 0.0 };

    let mut uv_pts: Vec<DVec2> = Vec::with_capacity(n);
    for i in 0..n {
        let t = t0 + dt * i as f64;
        let pt3d = curve.point_at(t);
        let proj = closest_point_on_surface(surface, pt3d, 8);
        if proj.distance > 1e-4 { continue; }
        uv_pts.push(DVec2::new(proj.params.0, proj.params.1));
    }

    if uv_pts.len() < 2 { return None; }
    uv_pts.dedup_by(|a, b| (*a - *b).length_squared() < 1e-20);

    let bspline = interpolate_points_2d(&uv_pts).ok()?;
    Some(Curve2d::BSpline(bspline))
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn project_onto_plane() {
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Y));
        let q = DVec3::new(3.0, 5.0, -2.0);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!(
            (r.point - DVec3::new(3.0, 0.0, -2.0)).length() < 1e-9,
            "point={}",
            r.point
        );
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Y),
        });
        let q = DVec3::new(5.0, 0.0, 0.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        assert!((r.point - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 3.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_cylinder() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
            ref_dir: DVec3::X,
        });
        let q = DVec3::new(3.0, 2.0, 0.0);
        let r = closest_point_on_surface(&cyl, q, 16);
        assert!((r.point - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 2.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_torus() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        // Point far along +X axis → nearest point on outer equator
        let q = DVec3::new(10.0, 0.0, 0.0);
        let r = closest_point_on_surface(&torus, q, 16);
        // Nearest should be at (4, 0, 0)
        assert!((r.point - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn project_onto_cone_returns_theta_and_slant_params() {
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 30.0_f64.to_radians(),
        });
        let expected_slant = 4.0;
        let on_surface = match &cone {
            Surface3::Cone(surface) => surface.point_at(0.0, expected_slant),
            _ => unreachable!(),
        };
        let query_normal = match &cone {
            Surface3::Cone(surface) => surface.normal_at(0.0, expected_slant),
            _ => unreachable!(),
        };
        let q = on_surface + query_normal * 0.25;
        let r = closest_point_on_surface(&cone, q, 16);

        assert!((r.point - on_surface).length() < 5e-3, "projected point={}", r.point);
        assert!((r.params.1 - expected_slant).abs() < 5e-3, "slant={}", r.params.1);
        let lifted = match &cone {
            Surface3::Cone(surface) => surface.point_at(r.params.0, r.params.1),
            _ => unreachable!(),
        };
        assert!((lifted - r.point).length() < 1e-6, "lifted point={lifted} projected={}", r.point);
    }

    #[test]
    fn project_onto_circle_curve() {
        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0,
        ));
        let q = DVec3::new(2.0, 0.0, 0.0);
        let r = closest_point_on_curve(&circle, q, 64);
        assert!(
            (r.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6,
            "expected (1,0,0) got {}",
            r.point
        );
        assert!((r.distance - 1.0).abs() < 1e-6);
    }

    #[test]
    fn project_onto_line_curve() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        // Line has infinite domain; nearest point to (3, 4, 0) is (3, 0, 0)
        let q = DVec3::new(3.0, 4.0, 0.0);
        let r = closest_point_on_curve(&line, q, 32);
        let expected = DVec3::new(3.0, 0.0, 0.0);
        assert!(
            (r.point - expected).length() < 1e-4,
            "expected {:?} got {}",
            expected,
            r.point
        );
        assert!((r.distance - 4.0).abs() < 1e-4, "distance={}", r.distance);
    }

    #[test]
    fn project_onto_ellipse_curve_analytic() {
        // Ellipse centered at origin in XY plane, semi-axes 3 and 1.
        let ellipse = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        // Query point along +X beyond the ellipse → nearest should be (3, 0, 0).
        let q = DVec3::new(5.0, 0.0, 0.0);
        let r = closest_point_on_curve(&ellipse, q, 64);
        assert!(
            (r.point - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-5,
            "expected (3,0,0) got {}",
            r.point
        );
        assert!((r.distance - 2.0).abs() < 1e-5, "distance={}", r.distance);
    }

    #[test]
    fn project_onto_ellipse_curve_off_plane() {
        // Query lifted off the ellipse plane — projection should still land on ellipse.
        let ellipse = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 2.0,
            minor_radius: 1.0,
        });
        // Query at (2, 0, 5) → closest in-plane point is (2, 0, 0).
        let q = DVec3::new(2.0, 0.0, 5.0);
        let r = closest_point_on_curve(&ellipse, q, 64);
        assert!(r.point.z.abs() < 1e-5, "z should be ~0, got {}", r.point.z);
        assert!(
            (r.point - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-5,
            "expected (2,0,0) got {}",
            r.point
        );
    }

    #[test]
    fn project_onto_line_curve_oblique() {
        // Line along (1,1,0)/sqrt(2), query off axis → test 3-D dot product.
        let dir = DVec3::new(1.0, 1.0, 0.0).normalize();
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: dir,
        });
        let q = DVec3::new(0.0, 1.0, 2.0);
        let r = closest_point_on_curve(&line, q, 32);
        // t = q·dir = 0*0.707 + 1*0.707 + 0 = 0.707, point = t*dir
        let t = q.dot(dir);
        let expected = dir * t;
        assert!(
            (r.point - expected).length() < 1e-9,
            "expected {:?} got {}",
            expected,
            r.point
        );
    }

    #[test]
    fn project_onto_partial_circle_arc() {
        // Arc from 0 to π/2 (first quadrant).
        let arc = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0,
        ));
        // For a full circle, query at (-2, 0, 0) should give t = π, point = (-1, 0, 0).
        let q = DVec3::new(-2.0, 0.0, 0.0);
        let r = closest_point_on_curve(&arc, q, 64);
        assert!(
            (r.point - DVec3::new(-1.0, 0.0, 0.0)).length() < 1e-6,
            "expected (-1,0,0) got {}",
            r.point
        );
    }

    #[test]
    fn project_onto_bspline_surface() {
        // Flat BSpline surface at z=0 over [0,1]²
        use crate::geom::BSplineSurface;
        let surf = Surface3::BSpline(BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0; 2]; 2],
        });
        let q = DVec3::new(0.5, 0.5, 5.0);
        let r = closest_point_on_surface(&surf, q, 8);
        assert!(
            (r.point - DVec3::new(0.5, 0.5, 0.0)).length() < 1e-4,
            "got {}",
            r.point
        );
    }

    // ============================================================================
    // OCCT TKGeomBase Alignment Tests - Projection Edge Cases
    // ============================================================================

    #[test]
    fn project_onto_cone_surface() {
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: std::f64::consts::FRAC_PI_6,
        });
        // Query near the cone surface (not at apex - apex is singular)
        let q = DVec3::new(1.0, 0.0, 1.0); // Near the cone surface
        let r = closest_point_on_surface(&cone, q, 16);
        assert!(r.distance < 0.5, "near-surface projection should be close");

        // Query along axis away from apex
        let q2 = DVec3::new(0.0, 0.0, 5.0);
        let r2 = closest_point_on_surface(&cone, q2, 16);
        assert!(r2.distance > 0.0, "axis projection distance should be positive");
    }

    #[test]
    fn project_onto_torus_surface() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        // Query at center (inside hole)
        let q = DVec3::new(0.0, 0.0, 0.0);
        let r = closest_point_on_surface(&torus, q, 16);
        assert!(r.distance > 0.0, "center distance should be positive");

        // Query on outer ring
        let q2 = DVec3::new(4.0, 0.0, 0.0);
        let r2 = closest_point_on_surface(&torus, q2, 16);
        assert!((r2.distance - 0.0).abs() < 0.1, "on-torus distance should be small");
    }

    #[test]
    fn project_onto_cylinder_interior() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: DVec3::X,
        });
        // Query inside cylinder
        let q = DVec3::new(0.0, 0.0, 1.0);
        let r = closest_point_on_surface(&cyl, q, 16);
        assert!((r.distance - 2.0).abs() < 1e-6, "interior distance should be radius");
    }

    #[test]
    fn project_onto_sphere_interior() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        // Query inside sphere
        let q = DVec3::new(1.0, 1.0, 1.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        assert!(r.distance < 3.0, "interior distance should be less than radius");
    }

    #[test]
    fn project_onto_plane_offset() {
        let plane = Surface3::Plane(Plane::new(DVec3::new(0.0, 5.0, 0.0), DVec3::Y));
        let q = DVec3::new(3.0, 0.0, -2.0);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!((r.point.y - 5.0).abs() < 1e-9);
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_line_at_origin() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let q = DVec3::new(0.0, 5.0, 0.0);
        let r = closest_point_on_curve(&line, q, 16);
        assert!((r.point - DVec3::ZERO).length() < 1e-9);
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_circle_at_parameter() {
        let circle = Curve3::Circle(Circle3::new(DVec3::new(2.0, 3.0, 0.0), DVec3::Z, 1.0));
        // Query at parameter 0 (should be at center + radius * X)
        let q = DVec3::new(3.0, 3.0, 0.0);
        let r = closest_point_on_curve(&circle, q, 32);
        assert!((r.point - DVec3::new(3.0, 3.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn project_distant_point_onto_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        // Very distant query
        let q = DVec3::new(1000.0, 1000.0, 1000.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        // Distance should be approximately |q| - radius
        let expected_dist = q.length() - 1.0;
        assert!((r.distance - expected_dist).abs() < 1.0, "distant projection distance");
    }

    #[test]
    fn project_near_surface_boundary() {
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        // Query very close to plane
        let q = DVec3::new(1.0, 2.0, 1e-10);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!(r.distance < 1e-9, "near-surface distance should be tiny");
    }
}
