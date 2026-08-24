//! Extrema algorithms: point-curve, point-surface, and curve-curve extrema.
//!
//! Analogous to OCCT's Extrema package:
//! - Extrema_ExtPElC — point to elementary curve (analytic: Line/Circle/Ellipse/Hyperbola/Parabola)
//! - Extrema_ExtPC — point to curve (with refinement)
//! - Extrema_ExtPS — point to surface (numerical grid + Newton)
//! - Extrema_GenLocateExtPS — local Newton from an initial UV guess
//! - Extrema_ExtCC — curve-curve extrema
//!
//! This is the low-level math module; higher-level GeomAPI wrappers live in
//! base::geom_api.

use crate::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3, SurfaceEval};
use crate::math::direct_polynomial_roots::DirectPolynomialRoots;
use glam::DVec2;
use glam::DVec3;

// =============================================================================
// Result types
// =============================================================================

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

/// A single local-minimum pair returned by [`extrema_curve_curve`].
#[derive(Debug, Clone)]
pub struct ExtremaPair {
    /// Parameter on the first curve at the closest approach.
    pub param1: f64,
    /// Parameter on the second curve at the closest approach.
    pub param2: f64,
    /// Point on the first curve.
    pub point1: DVec3,
    /// Point on the second curve.
    pub point2: DVec3,
    /// Euclidean distance at closest approach.
    pub distance: f64,
}

/// Point on a curve at an extremum distance (OCCT `Extrema_POnCurv`).
#[derive(Debug, Clone)]
pub struct POnCurve {
    /// Parameter on the curve.
    pub param: f64,
    /// Point on the curve.
    pub point: DVec3,
}

/// Point on a surface at an extremum distance (OCCT `Extrema_POnSurf`).
#[derive(Debug, Clone)]
pub struct POnSurface {
    /// Parameter U on the surface.
    pub u: f64,
    /// Parameter V on the surface.
    pub v: f64,
    /// Point on the surface.
    pub point: DVec3,
}

// ============================================================================
// Extrema_LocateExtPC — local Point-to-Curve extremum from a seed parameter
// ============================================================================

/// OCCT-aligned: local extremum search starting from a seed parameter.
///
/// OCCT: `Extrema_LocateExtPC(Point, Adaptor3d_Curve, Seed, Tol)`.
/// Performs Newton refinement starting from `seed` to find a local minimum
/// of the distance function. Returns None if the search fails to converge
/// or the distance increases.
pub fn extrema_locate_ext_pc(
    point: DVec3,
    curve: &Curve3,
    seed: f64,
    uinf: f64,
    usup: f64,
    tol: f64,
) -> Option<POnCurve> {
    let clamp = |t: f64| t.clamp(uinf, usup);
    let dt = 1e-7;
    let max_iter = 50;
    let mut t = clamp(seed);
    let mut best_t = t;
    let mut best_d = (curve.point_at(t) - point).length();
    for _ in 0..max_iter {
        let p = curve.point_at(t);
        let diff = p - point;
        let deriv = curve.derivative_at(t);
        let deriv_sq = deriv.dot(deriv);
        if deriv_sq < 1e-20 {
            break;
        }
        let curv = (curve.point_at(t + 2.0 * dt) - 2.0 * p + curve.point_at(t - 2.0 * dt))
            / (dt * dt);
        let denom = deriv_sq + diff.dot(curv);
        let delta = diff.dot(deriv) / if denom.abs() > 1e-20 { denom } else { deriv_sq };
        let new_t = clamp(t - delta);
        let new_d = (curve.point_at(new_t) - point).length();
        if new_d < best_d {
            best_d = new_d;
            best_t = new_t;
            t = new_t;
        } else {
            break;
        }
        if delta.abs() < tol {
            break;
        }
    }
    let best_p = curve.point_at(best_t);
    Some(POnCurve {
        param: best_t,
        point: best_p,
    })
}

// ============================================================================
// Extrema_ExtPC — Point-to-Curve extremum (OCCT-aligned class)
// ============================================================================

/// OCCT-aligned: computes all extremum distances between a point and a curve.
///
/// Uses coarse grid sampling to find seeds, then Newton refinement for each
/// local minimum.  Mirrors `Extrema_ExtPC` / `Extrema_GGExtPC`.
///
/// OCCT: `Extrema_ExtPC(Point, Curve, TolC, Uinf, Usup)`.
pub struct ExtPC {
    done: bool,
    points: Vec<POnCurve>,
    sq_dists: Vec<f64>,
    tol: f64,
}

impl ExtPC {
    /// Constructor: compute all extrema between `point` and `curve` on [uinf, usup].
    ///
    /// OCCT: `Extrema_ExtPC(gp_Pnt, Adaptor3d_Curve, TolC, Uinf, Usup)`.
    pub fn new(point: DVec3, curve: &Curve3, tol: f64, uinf: f64, usup: f64) -> Self {
        let mut ext = ExtPC {
            done: false,
            points: Vec::new(),
            sq_dists: Vec::new(),
            tol: tol.max(1e-12),
        };
        ext.perform(point, curve, uinf, usup);
        ext
    }

    /// Compute the extrema (can be called after default construction + Initialize).
    ///
    /// OCCT: `Perform(Point)`.
    pub fn perform(&mut self, point: DVec3, curve: &Curve3, uinf: f64, usup: f64) {
        self.points.clear();
        self.sq_dists.clear();

        // 1. Analytic path: handle elementary curve types directly
        let analytic_result = match curve {
            Curve3::Line(l) => self.ext_pelc_line(point, l, uinf, usup),
            Curve3::Circle(c) => self.ext_pelc_circle(point, c, uinf, usup),
            _ => None,
        };

        if let Some((pts, sqs)) = analytic_result {
            self.points = pts;
            self.sq_dists = sqs;
            self.done = true;
            return;
        }

        // 2. General path: coarse grid + Newton refinement
        let (t_min, t_max) = (uinf.max(f64::NEG_INFINITY), usup.min(f64::INFINITY));
        if !t_min.is_finite() || !t_max.is_finite() || (t_max - t_min).abs() < self.tol {
            self.done = true;
            return;
        }

        const N_GRID: usize = 51;
        let mut candidates: Vec<(f64, f64)> = Vec::new(); // (t, dist²)

        for i in 0..=N_GRID {
            let t = t_min + (t_max - t_min) * (i as f64) / (N_GRID as f64);
            let p = curve.point_at(t);
            let d2 = (p - point).length_squared();
            // Keep local minima
            if (i == 0 || d2 <= candidates.last().map(|&(_, ld)| ld).unwrap_or(f64::INFINITY))
                && (i == N_GRID || {
                    let next_t = t_min + (t_max - t_min) * ((i + 1) as f64) / (N_GRID as f64);
                    let next_p = curve.point_at(next_t);
                    let next_d2 = (next_p - point).length_squared();
                    d2 <= next_d2
                })
            {
                candidates.push((t, d2));
            }
        }

        // Deduplicate close candidates
        candidates.dedup_by(|a, b| (a.0 - b.0).abs() < (t_max - t_min) / (N_GRID as f64) * 0.5);

        // 3. Newton refinement for each candidate
        for &(t0, _) in &candidates {
            let t = self.newton_refine_curve(point, curve, t0, t_min, t_max);
            let p = curve.point_at(t);
            let d2 = (p - point).length_squared();

            // Deduplicate against already-found solutions
            let is_dup = self.points.iter().any(|existing| {
                let dt = (existing.param - t).abs();
                let dp = (existing.point - p).length();
                dt < self.tol && dp < self.tol * 10.0
            });

            if !is_dup {
                self.points.push(POnCurve { param: t, point: p });
                self.sq_dists.push(d2);
            }
        }

        // Sort by distance
        let mut indices: Vec<usize> = (0..self.points.len()).collect();
        indices.sort_by(|&a, &b| self.sq_dists[a].partial_cmp(&self.sq_dists[b]).unwrap());
        self.points = indices.iter().map(|&i| self.points[i].clone()).collect();
        self.sq_dists = indices.iter().map(|&i| self.sq_dists[i]).collect();

        self.done = true;
    }

    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT: `NbExt()`.
    pub fn nb_ext(&self) -> usize {
        self.points.len()
    }

    /// OCCT: `SquareDistance(N)` — 1-indexed.
    pub fn square_distance(&self, n: usize) -> f64 {
        assert!(n >= 1 && n <= self.sq_dists.len(), "ExtPC: index out of range");
        self.sq_dists[n - 1]
    }

    /// OCCT: `Point(N)` — returns the point on curve. 1-indexed.
    pub fn point(&self, n: usize) -> &POnCurve {
        assert!(n >= 1 && n <= self.points.len(), "ExtPC: index out of range");
        &self.points[n - 1]
    }

    // --- Analytic sub-solvers ---

    fn ext_pelc_line(&self, point: DVec3, line: &crate::geom::Line3, uinf: f64, usup: f64) -> Option<(Vec<POnCurve>, Vec<f64>)> {
        // Project point onto infinite line: t = (p - o)·d
        let t = (point - line.origin).dot(line.direction);
        let t_clamped = t.clamp(uinf, usup);
        let p = line.origin + t_clamped * line.direction;
        let d2 = (p - point).length_squared();
        Some((
            vec![POnCurve { param: t_clamped, point: p }],
            vec![d2],
        ))
    }

    fn ext_pelc_circle(&self, point: DVec3, circle: &crate::geom::Circle3, uinf: f64, usup: f64) -> Option<(Vec<POnCurve>, Vec<f64>)> {
        // Project point onto circle center plane, then find angle
        let d = point - circle.center;
        let along = d.dot(circle.normal);
        let planar = d - circle.normal * along;
        let r = planar.length();
        if r < 1e-15 {
            // Point at center: return any point on circle
            let t = uinf;
            let p = circle.center + circle.x_dir * circle.radius * t.cos()
                + circle.y_dir * circle.radius * t.sin();
            let d2 = (p - point).length_squared();
            return Some((vec![POnCurve { param: t, point: p }], vec![d2]));
        }
        let angle = planar.dot(circle.y_dir).atan2(planar.dot(circle.x_dir));
        let t = angle.clamp(uinf, usup);
        let p = circle.center + circle.x_dir * circle.radius * t.cos()
            + circle.y_dir * circle.radius * t.sin();
        let d2 = (p - point).length_squared();

        // May also need to check endpoints if domain is restricted
        let mut pts = vec![POnCurve { param: t, point: p }];
        let mut sqs = vec![d2];

        if uinf.is_finite() && usup.is_finite() {
            for &bound in &[uinf, usup] {
                let bp = circle.center + circle.x_dir * circle.radius * bound.cos()
                    + circle.y_dir * circle.radius * bound.sin();
                let bd2 = (bp - point).length_squared();
                if bd2 < d2 - 1e-12 {
                    pts.push(POnCurve { param: bound, point: bp });
                    sqs.push(bd2);
                }
            }
        }

        Some((pts, sqs))
    }

    fn newton_refine_curve(&self, point: DVec3, curve: &Curve3, t0: f64, t_min: f64, t_max: f64) -> f64 {
        let mut t = t0.clamp(t_min, t_max);
        for _ in 0..20 {
            let p = curve.point_at(t);
            let dp = curve.derivative_at(t);
            let d = p - point;
            let f = d.dot(dp);
            let speed_sq = dp.length_squared();
            if speed_sq < 1e-30 || f.abs() < self.tol {
                break;
            }
            let d2 = if speed_sq > 1e-30 {
                let ddp = curve.derivative2_at(t);
                let df = dp.dot(dp) + d.dot(ddp);
                df
            } else {
                speed_sq
            };
            if d2.abs() < 1e-30 {
                break;
            }
            let dt = -f / d2;
            t = (t + dt).clamp(t_min, t_max);
            if dt.abs() < self.tol {
                break;
            }
        }
        t
    }
}

// ============================================================================
// Extrema_ExtPC2d — Point-to-Curve extremum (2D)
// ============================================================================

/// Point on a 2D curve at an extremum distance (OCCT `Extrema_POnCurv2d`).
#[derive(Debug, Clone)]
pub struct POnCurve2d {
    /// Parameter on the curve.
    pub param: f64,
    /// Point on the curve.
    pub point: DVec2,
}

/// OCCT-aligned: computes the extremum distances between a 2D point and a 2D
/// curve. Mirrors `Extrema_ExtPC2d` =
/// `Extrema_GGExtPC<Adaptor2d_Curve2d, Extrema_Curve2dTool, ...>` (used by
/// BRepClass_Intersector::CheckOn).
pub struct ExtPC2d {
    done: bool,
    points: Vec<POnCurve2d>,
    sq_dists: Vec<f64>,
    tol: f64,
}

impl ExtPC2d {
    /// Constructor: compute all extrema between `point` and `curve` on
    /// `[uinf, usup]`.
    ///
    /// OCCT: `Extrema_ExtPC2d(gp_Pnt2d, Adaptor2d_Curve2d)`.
    pub fn new(point: DVec2, curve: &Curve2d, tol: f64, uinf: f64, usup: f64) -> Self {
        let mut ext = ExtPC2d {
            done: false,
            points: Vec::new(),
            sq_dists: Vec::new(),
            tol: tol.max(1e-12),
        };
        ext.perform(point, curve, uinf, usup);
        ext
    }

    /// Compute the extrema.
    pub fn perform(&mut self, point: DVec2, curve: &Curve2d, uinf: f64, usup: f64) {
        self.points.clear();
        self.sq_dists.clear();

        // Analytic path for line (elementary-curve case).
        if let Curve2d::Line(l) = curve {
            let (pts, sqs) = self.ext_pelc_line2d(point, l, uinf, usup);
            self.points = pts;
            self.sq_dists = sqs;
            self.done = true;
            return;
        }

        // General path: coarse grid + Newton refinement.
        let (t_min, t_max) = (uinf.max(f64::NEG_INFINITY), usup.min(f64::INFINITY));
        if !t_min.is_finite() || !t_max.is_finite() || (t_max - t_min).abs() < self.tol {
            self.done = true;
            return;
        }

        const N_GRID: usize = 51;
        let mut candidates: Vec<(f64, f64)> = Vec::new(); // (t, dist²)
        for i in 0..=N_GRID {
            let t = t_min + (t_max - t_min) * (i as f64) / (N_GRID as f64);
            let p = curve.point_at(t);
            let d2 = (p - point).length_squared();
            // Keep local minima.
            if (i == 0 || d2 <= candidates.last().map(|&(_, ld)| ld).unwrap_or(f64::INFINITY))
                && (i == N_GRID || {
                    let next_t = t_min + (t_max - t_min) * ((i + 1) as f64) / (N_GRID as f64);
                    let next_p = curve.point_at(next_t);
                    let next_d2 = (next_p - point).length_squared();
                    d2 <= next_d2
                })
            {
                candidates.push((t, d2));
            }
        }
        candidates.dedup_by(|a, b| (a.0 - b.0).abs() < (t_max - t_min) / (N_GRID as f64) * 0.5);

        for &(t0, _) in &candidates {
            let t = self.newton_refine_curve2d(point, curve, t0, t_min, t_max);
            let p = curve.point_at(t);
            let d2 = (p - point).length_squared();
            let is_dup = self.points.iter().any(|existing| {
                let dt = (existing.param - t).abs();
                let dp = (existing.point - p).length();
                dt < self.tol && dp < self.tol * 10.0
            });
            if !is_dup {
                self.points.push(POnCurve2d { param: t, point: p });
                self.sq_dists.push(d2);
            }
        }

        // Sort by distance.
        let mut indices: Vec<usize> = (0..self.points.len()).collect();
        indices.sort_by(|&a, &b| self.sq_dists[a].partial_cmp(&self.sq_dists[b]).unwrap());
        self.points = indices.iter().map(|&i| self.points[i].clone()).collect();
        self.sq_dists = indices.iter().map(|&i| self.sq_dists[i]).collect();

        self.done = true;
    }

    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT: `NbExt()`.
    pub fn nb_ext(&self) -> usize {
        self.points.len()
    }

    /// OCCT: `SquareDistance(N)` — 1-indexed.
    pub fn square_distance(&self, n: usize) -> f64 {
        assert!(n >= 1 && n <= self.sq_dists.len(), "ExtPC2d: index out of range");
        self.sq_dists[n - 1]
    }

    /// OCCT: `Point(N)` — returns the point on curve. 1-indexed.
    pub fn point(&self, n: usize) -> &POnCurve2d {
        assert!(n >= 1 && n <= self.points.len(), "ExtPC2d: index out of range");
        &self.points[n - 1]
    }

    fn ext_pelc_line2d(
        &self,
        point: DVec2,
        line: &Line2d,
        uinf: f64,
        usup: f64,
    ) -> (Vec<POnCurve2d>, Vec<f64>) {
        // Project point onto the infinite line: t = (p - o)·d.
        let t = (point - line.origin).dot(line.direction);
        let t_clamped = t.clamp(uinf, usup);
        let p = line.origin + t_clamped * line.direction;
        let d2 = (p - point).length_squared();
        (vec![POnCurve2d { param: t_clamped, point: p }], vec![d2])
    }

    fn newton_refine_curve2d(
        &self,
        point: DVec2,
        curve: &Curve2d,
        t0: f64,
        t_min: f64,
        t_max: f64,
    ) -> f64 {
        let mut t = t0.clamp(t_min, t_max);
        for _ in 0..20 {
            let p = curve.point_at(t);
            let dp = curve.derivative_at(t);
            let d = p - point;
            let f = d.dot(dp);
            let speed_sq = dp.length_squared();
            if speed_sq < 1e-30 || f.abs() < self.tol {
                break;
            }
            let d2 = if speed_sq > 1e-30 {
                let ddp = curve.derivative2_at(t);
                dp.dot(dp) + d.dot(ddp)
            } else {
                speed_sq
            };
            if d2.abs() < 1e-30 {
                break;
            }
            let dt = -f / d2;
            t = (t + dt).clamp(t_min, t_max);
            if dt.abs() < self.tol {
                break;
            }
        }
        t
    }
}

// ============================================================================
// Extrema_GenLocateExtPS — Point-to-Surface local extremum (OCCT-aligned)
// ============================================================================

/// OCCT-aligned: computes the closest point on a surface from an initial (u,v) guess.
///
/// Uses Newton iteration on the surface param space.  Mirrors
/// `Extrema_GenLocateExtPS` / the local-refinement path of `Extrema_ExtPS`.
///
/// OCCT: `Extrema_GenLocateExtPS(Point, Surface, U, V, TolU, TolV, Uinf, Usup, Vinf, Vsup)`.
pub struct GenLocateExtPS {
    done: bool,
    u: f64,
    v: f64,
    point: DVec3,
    sq_dist: f64,
}

impl GenLocateExtPS {
    /// Constructor and compute in one call.
    ///
    /// OCCT: `GenLocateExtPS(gp_Pnt, Surface, U0, V0, TolU, TolV, Uinf, Usup, Vinf, Vsup)`.
    pub fn new(
        point: DVec3,
        surface: &Surface3,
        u0: f64,
        v0: f64,
        tol_u: f64,
        tol_v: f64,
        uinf: f64,
        usup: f64,
        vinf: f64,
        vsup: f64,
    ) -> Self {
        let mut ext = GenLocateExtPS {
            done: false,
            u: u0.clamp(uinf, usup),
            v: v0.clamp(vinf, vsup),
            point: DVec3::ZERO,
            sq_dist: f64::INFINITY,
        };
        ext.perform(point, surface, tol_u, tol_v, uinf, usup, vinf, vsup);
        ext
    }

    /// Perform the local search.
    ///
    /// OCCT: `Perform(Point)`.
    pub fn perform(&mut self, point: DVec3, surface: &Surface3, tol_u: f64, tol_v: f64, uinf: f64, usup: f64, vinf: f64, vsup: f64) {
        let tol_u = tol_u.max(1e-12);
        let tol_v = tol_v.max(1e-12);

        for _ in 0..30 {
            let (p, pu, pv) = surface.derivatives(self.u, self.v);
            self.point = p;

            let d = p - point;
            let gu = d.dot(pu);
            let gv = d.dot(pv);
            self.sq_dist = d.length_squared();

            if gu.abs() < tol_u && gv.abs() < tol_v {
                self.done = true;
                return;
            }

            let huu = pu.dot(pu);
            let hvv = pv.dot(pv);
            let huv = pu.dot(pv);
            let det = huu * hvv - huv * huv;
            if det.abs() < 1e-30 {
                break;
            }

            let du = (hvv * gu - huv * gv) / det;
            let dv = (huu * gv - huv * gu) / det;

            self.u = (self.u - du).clamp(uinf, usup);
            self.v = (self.v - dv).clamp(vinf, vsup);

            if du.abs() < tol_u && dv.abs() < tol_v {
                self.done = true;
                return;
            }
        }
        self.done = true;
    }

    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT: `SquareDistance()`.
    pub fn square_distance(&self) -> f64 {
        self.sq_dist
    }

    /// OCCT: `Point()` — returns the point on surface.
    pub fn point_on_surface(&self) -> DVec3 {
        self.point
    }

    /// Parameter (U, V) of the solution.
    /// OCCT: accessed via `Extrema_POnSurf`.
    pub fn parameters(&self) -> (f64, f64) {
        (self.u, self.v)
    }
}

// ============================================================================
// Extrema_ExtPS — Point-to-Surface full extrema (OCCT-aligned class)
// ============================================================================

/// OCCT-aligned: computes all extremum distances between a point and a surface.
///
/// Uses a coarse (Nu × Nv) grid to find seeds, then Newton refinement for each
/// candidate.  Mirrors `Extrema_ExtPS`.
///
/// OCCT: `Extrema_ExtPS(Point, Surface, TolU, TolV)` with optional domain.
pub struct ExtPS {
    done: bool,
    points: Vec<POnSurface>,
    sq_dists: Vec<f64>,
}

impl ExtPS {
    /// Constructor with automatic domain from the surface.
    ///
    /// OCCT: `Extrema_ExtPS(gp_Pnt, Adaptor3d_Surface, TolU, TolV)`.
    pub fn new(point: DVec3, surface: &Surface3, tol_u: f64, tol_v: f64) -> Self {
        let dom = surface.default_domain();
        let (uinf, usup, vinf, vsup) = (dom[0], dom[1], dom[2], dom[3]);
        Self::with_domain(point, surface, uinf, usup, vinf, vsup, tol_u, tol_v)
    }

    /// Constructor with explicit domain.
    ///
    /// OCCT: `Extrema_ExtPS(gp_Pnt, Adaptor3d_Surface, Uinf, Usup, Vinf, Vsup, TolU, TolV)`.
    pub fn with_domain(
        point: DVec3,
        surface: &Surface3,
        uinf: f64,
        usup: f64,
        vinf: f64,
        vsup: f64,
        tol_u: f64,
        tol_v: f64,
    ) -> Self {
        let mut ext = ExtPS {
            done: false,
            points: Vec::new(),
            sq_dists: Vec::new(),
        };
        ext.perform(point, surface, uinf, usup, vinf, vsup, tol_u, tol_v);
        ext
    }

    /// Perform the computation.
    ///
    /// OCCT: `Perform(Point)` (after Initialize).
    pub fn perform(
        &mut self,
        point: DVec3,
        surface: &Surface3,
        uinf: f64,
        usup: f64,
        vinf: f64,
        vsup: f64,
        tol_u: f64,
        tol_v: f64,
    ) {
        self.points.clear();
        self.sq_dists.clear();

        if !uinf.is_finite() || !usup.is_finite() || !vinf.is_finite() || !vsup.is_finite() {
            self.done = true;
            return;
        }

        let range_u = usup - uinf;
        let range_v = vsup - vinf;
        if range_u < tol_u || range_v < tol_v {
            self.done = true;
            return;
        }

        // 1. Coarse grid to find candidate seeds
        let n_u = (range_u / tol_u).ceil() as usize;
        let n_v = (range_v / tol_v).ceil() as usize;
        let n_u = n_u.clamp(5, 30);
        let n_v = n_v.clamp(5, 30);

        let mut grid: Vec<Vec<f64>> = vec![vec![0.0; n_v]; n_u];
        for i in 0..n_u {
            let u = uinf + range_u * (i as f64) / ((n_u - 1).max(1) as f64);
            for j in 0..n_v {
                let v = vinf + range_v * (j as f64) / ((n_v - 1).max(1) as f64);
                let p = surface.point_at(u, v);
                grid[i][j] = (p - point).length_squared();
            }
        }

        // 2. Find local minima in the grid
        let mut seeds: Vec<(f64, f64)> = Vec::new();
        for i in 0..n_u {
            for j in 0..n_v {
                let val = grid[i][j];
                let is_min = {
                    let i_min = if i > 0 { i - 1 } else { i };
                    let i_max = if i + 1 < n_u { i + 1 } else { i };
                    let j_min = if j > 0 { j - 1 } else { j };
                    let j_max = if j + 1 < n_v { j + 1 } else { j };

                    (i_min..=i_max).all(|ii| {
                        (j_min..=j_max).all(|jj| {
                            if ii == i && jj == j { true } else { grid[ii][jj] >= val }
                        })
                    })
                };

                if is_min {
                    let u = uinf + range_u * (i as f64) / ((n_u - 1).max(1) as f64);
                    let v = vinf + range_v * (j as f64) / ((n_v - 1).max(1) as f64);
                    seeds.push((u, v));
                }
            }
        }

        // Deduplicate seeds
        let tol_u_v = (range_u / n_u as f64 * 0.5).max(tol_u);
        let tol_v_v = (range_v / n_v as f64 * 0.5).max(tol_v);
        seeds.dedup_by(|a, b| (a.0 - b.0).abs() < tol_u_v && (a.1 - b.1).abs() < tol_v_v);

        // 3. Newton refinement from each seed
        for &(u0, v0) in &seeds {
            let mut loc = GenLocateExtPS::new(
                point, surface, u0, v0, tol_u, tol_v, uinf, usup, vinf, vsup,
            );
            if loc.is_done() {
                let u = loc.parameters().0;
                let v = loc.parameters().1;
                let p = loc.point_on_surface();
                let d2 = loc.square_distance();

                // Deduplicate
                let is_dup = self.points.iter().any(|existing| {
                    (existing.u - u).abs() < tol_u && (existing.v - v).abs() < tol_v
                });

                if !is_dup {
                    self.points.push(POnSurface { u, v, point: p });
                    self.sq_dists.push(d2);
                }
            }
        }

        // Sort by distance
        let mut indices: Vec<usize> = (0..self.points.len()).collect();
        indices.sort_by(|&a, &b| self.sq_dists[a].partial_cmp(&self.sq_dists[b]).unwrap());
        self.points = indices.iter().map(|&i| self.points[i].clone()).collect();
        self.sq_dists = indices.iter().map(|&i| self.sq_dists[i]).collect();

        self.done = true;
    }

    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT: `NbExt()`.
    pub fn nb_ext(&self) -> usize {
        self.points.len()
    }

    /// OCCT: `SquareDistance(N)` — 1-indexed.
    pub fn square_distance(&self, n: usize) -> f64 {
        assert!(n >= 1 && n <= self.sq_dists.len(), "ExtPS: index out of range");
        self.sq_dists[n - 1]
    }

    /// OCCT: `Point(N)` — 1-indexed, returns the point on surface.
    pub fn point(&self, n: usize) -> &POnSurface {
        assert!(n >= 1 && n <= self.points.len(), "ExtPS: index out of range");
        &self.points[n - 1]
    }
}

/// Collection of all local minima found between two curves.
#[derive(Debug, Clone)]
pub struct CurveCurveExtrema {
    /// All local minima, sorted by distance ascending.
    pub pairs: Vec<ExtremaPair>,
}

impl CurveCurveExtrema {
    /// Convenience: distance of the global (closest) minimum.
    /// Returns `f64::INFINITY` if no pairs were found.
    pub fn min_distance(&self) -> f64 {
        self.pairs
            .first()
            .map(|p| p.distance)
            .unwrap_or(f64::INFINITY)
    }
}

// =============================================================================
// Point-Curve Extrema (Extrema_ExtPElC + Extrema_ExtPC)
// =============================================================================

/// Newton refinement helper for point-curve projection.
/// Solves g(t) = (P - C(t))·C'(t) = 0 via Newton:
///   t_{n+1} = t_n - (P-C)·C' / (|C'|² + (P-C)·C'')
/// where C'' is approximated by finite-difference.
fn newton_refine_pc(
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
        if deriv_sq < 1e-20 {
            break;
        }
        let curv =
            (curve.point_at(*t + 2.0 * dt) - 2.0 * p + curve.point_at(*t - 2.0 * dt)) / (dt * dt);
        let denom = deriv_sq + diff.dot(curv);
        let delta = diff.dot(deriv) / if denom.abs() > 1e-20 { denom } else { deriv_sq };
        let new_t = clamp(*t - delta);
        let new_dist = (curve.point_at(new_t) - query).length();
        if new_dist < *best_dist {
            *best_dist = new_dist;
            *t = new_t;
        }
        if delta.abs() < 1e-10 {
            break;
        }
    }
}

/// Project the point `query` onto `curve`, returning the nearest point on the
/// curve, its parameter value, and the Euclidean distance.
///
/// OCCT-aligned: dispatches per-type matching Extrema_ExtPC:
///   - Line/Circle: analytic via Extrema_ExtPElC equivalent
///   - Ellipse: analytic init + Newton refinement
///   - Hyperbola/Parabola: analytic via Extrema_ExtPElC (solutions filtered
///     by [uinf, usup], OCCT Extrema_ExtPElC.cxx L365/L457)
///   - BSpline: C2 interval splitting (Extrema_GGExtPC)
///   - Bezier/Other: uniform sampling + Newton
///
/// `n_samples` is the uniform sampling count (used for Bezier and fallback;
/// for BSpline it is overridden by `degree + 1` per C2 interval).
pub fn closest_point_on_curve(curve: &Curve3, query: DVec3, n_samples: usize) -> CurveProjection {
    let [uinf, usup] = curve.default_domain();
    closest_point_on_curve_impl(curve, query, n_samples, uinf, usup)
}

/// Range-restricted projection — OCCT `GeomAPI_ProjectPointOnCurve::Init(P, Curve, Umin, Usup)`:
/// the parameter bounds [uinf, usup] are applied inside the analytic
/// elementary-curve solvers (Extrema_ExtPElC Uinf/Usup), not by clamping the
/// projection afterwards.
pub fn closest_point_on_curve_with_range(
    curve: &Curve3,
    query: DVec3,
    n_samples: usize,
    uinf: f64,
    usup: f64,
) -> CurveProjection {
    closest_point_on_curve_impl(curve, query, n_samples, uinf, usup)
}

fn closest_point_on_curve_impl(
    curve: &Curve3,
    query: DVec3,
    n_samples: usize,
    uinf: f64,
    usup: f64,
) -> CurveProjection {
    match curve {
        Curve3::Line(l) => {
            let dir_sq = l.direction.dot(l.direction);
            if dir_sq < 1e-20 {
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
            // OCCT-aligned: Extrema_ExtPElC::Perform(Circle)
            let o = circ.center;
            let axis = circ.normal.normalize_or_zero();
            let axial = (query - o).dot(axis);
            let pp = query - axis * axial;
            let opp = pp - o;
            let opp_mag = opp.length();
            if opp_mag < 1e-15 {
                let pt = circ.point_at(0.0);
                return CurveProjection {
                    point: pt,
                    param: 0.0,
                    distance: (pt - query).length(),
                };
            }
            let cx = circ.x_dir.normalize();
            let cy = circ.y_dir.normalize();
            let x = opp.dot(cx);
            let y = opp.dot(cy);
            let u_min = y.atan2(x);
            let half = std::f64::consts::PI;
            let u_max = if u_min < 0.0 { u_min + half } else { u_min - half };
            let [t0, t1] = curve.default_domain();
            let mut best_t = u_min;
            let mut best_d = f64::INFINITY;
            for &u in &[u_min, u_max] {
                let mut u_adj = u;
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
            if best_d.is_infinite() {
                best_t = u_min.clamp(t0, t1);
                best_d = (circ.point_at(best_t) - query).length();
            }
            let pt = circ.point_at(best_t);
            return CurveProjection { point: pt, param: best_t, distance: (pt - query).length() };
        }

        Curve3::Ellipse(ell) => {
            // OCCT-aligned: Extrema_ExtPElC::Perform(Ellipse)
            let o = ell.center;
            let axis = ell.normal.normalize_or_zero();
            let pp = query - axis * (query - o).dot(axis);
            let opp = pp - o;
            if opp.length_squared() < 1e-30 && (ell.major_radius - ell.minor_radius).abs() < 1e-15 {
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
            let clamp = |u: f64| u.clamp(t0, t1);
            let mut g_prev = g(t0);
            for i in 1..=n {
                let u = t0 + (t1 - t0) * i as f64 / n as f64;
                let g_cur = g(u);
                if g_prev * g_cur <= 0.0 || g_prev.abs() < 1e-12 || g_cur.abs() < 1e-12 {
                    let u_mid = (u + (t0 + (t1 - t0) * (i - 1) as f64 / n as f64)) * 0.5;
                    let mut t = clamp(u_mid);
                    let mut dist = (ell.point_at(t) - query).length();
                    newton_refine_pc(curve, &mut t, &mut dist, query, 20, clamp);
                    if dist < d_best { d_best = dist; u_best = t; }
                }
                g_prev = g_cur;
            }
            for &u in &[t0, t1] {
                let d = (ell.point_at(u) - query).length();
                if d < d_best { d_best = d; u_best = u; }
            }
            let pt = ell.point_at(u_best);
            return CurveProjection { point: pt, param: u_best, distance: (pt - query).length() };
        }

        // OCCT-aligned: Extrema_ExtPElC::Perform(Hyperbola) (Extrema_ExtPElC.cxx L327-389)
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
            // OCCT L347-348: C1 = (R*R + r*r)/4; Sol(C1, -(X*R+Y*r)/2, 0, (X*R-Y*r)/2, -C1)
            let c1 = (r * r + r2 * r2) / 4.0;
            let sol = DirectPolynomialRoots::new_quartic(
                c1, -(x * r + y * r2) / 2.0, 0.0, (x * r - y * r2) / 2.0, -c1);
            // OCCT L365: the solution is kept when Uinf <= Us <= Usup (the
            // bounds passed by GeomAPI_ProjectPointOnCurve::Init(curve, f, l)).
            let [t0, t1] = [uinf, usup];
            let mut u_best = t0;
            let mut d_best = f64::INFINITY;
            if sol.is_done() {
                // OCCT L355-387: for (NoSol=1..NbSol) { Vs=Value(NoSol); if (Vs>0) { Us=log(Vs); if Uinf<=Us<=Usup ... } }
                for no_sol in 1..=sol.nb_solutions() {
                    let v_s = sol.value(no_sol);
                    if v_s > 0.0 {
                        let u_s = v_s.ln();
                        if u_s >= t0 && u_s <= t1 {
                            let pt = hyp.point_at(u_s);
                            let d = (pt - query).length();
                            if d < d_best { d_best = d; u_best = u_s; }
                        }
                    }
                }
            }
            for &u in &[t0, t1] {
                let d = (hyp.point_at(u) - query).length();
                if d < d_best { d_best = d; u_best = u; }
            }
            let pt = hyp.point_at(u_best);
            return CurveProjection { point: pt, param: u_best, distance: (pt - query).length() };
        }

        // OCCT-aligned: Extrema_ExtPElC::Perform(Parabola) (Extrema_ExtPElC.cxx L427-480)
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
            // OCCT L443: Sol(1/(4*F), 0., 2*F - X, -2*F*Y)
            let sol = DirectPolynomialRoots::new_cubic(
                1.0 / (4.0 * f), 0.0, 2.0 * f - x, -2.0 * f * y);
            // OCCT L457: the solution is kept when Uinf <= Us <= Usup.
            let [t0, t1] = [uinf, usup];
            let mut u_best = t0;
            let mut d_best = f64::INFINITY;
            if sol.is_done() {
                // OCCT L454-478: for (NoSol=1..NbSol) { Us=Value(NoSol); if Uinf<=Us<=Usup ... }
                for no_sol in 1..=sol.nb_solutions() {
                    let u_s = sol.value(no_sol);
                    if u_s >= t0 && u_s <= t1 {
                        let pt = par.point_at(u_s);
                        let d = (pt - query).length();
                        if d < d_best { d_best = d; u_best = u_s; }
                    }
                }
            }
            for &u in &[t0, t1] {
                let d = (par.point_at(u) - query).length();
                if d < d_best { d_best = d; u_best = u; }
            }
            let pt = par.point_at(u_best);
            return CurveProjection { point: pt, param: u_best, distance: (pt - query).length() };
        }

        // OCCT-aligned: BSpline — C2 interval splitting
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
                for i in 0..=n_per_int {
                    let t = lo + span * i as f64 / n_per_int as f64;
                    let p = bs.point_at(t);
                    let d = (p - query).length();
                    if d < best_dist { best_dist = d; best_t = t; }
                }
            }
            if best_dist < f64::INFINITY {
                newton_refine_pc(curve, &mut best_t, &mut best_dist, query, 30, clamp_t);
            }
            let best_point = curve.point_at(best_t);
            return CurveProjection { point: best_point, param: best_t, distance: (best_point - query).length() };
        }

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
                if d < best_dist { best_dist = d; best_t = t; }
            }
            newton_refine_pc(curve, &mut best_t, &mut best_dist, query, 30, clamp_t);
            let best_point = curve.point_at(best_t);
            return CurveProjection { point: best_point, param: best_t, distance: (best_point - query).length() };
        }

        _ => {}
    }

    // Numerical fallback for all other curve types
    let [t0_raw, t1_raw] = curve.default_domain();
    let n = n_samples.max(4);
    let (t0, t1) = if t0_raw.is_infinite() || t1_raw.is_infinite() {
        (0.0 - 100.0, 0.0 + 100.0)
    } else {
        (t0_raw, t1_raw)
    };

    let (mut best_t, mut best_dist) = {
        let mut bd = f64::INFINITY;
        let mut bt = t0;
        for i in 0..=n {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            let p = curve.point_at(t);
            let d = (p - query).length();
            if d < bd { bd = d; bt = t; }
        }
        (bt, bd)
    };

    let clamp_t = |t: f64| {
        if t0_raw.is_infinite() || t1_raw.is_infinite() { t }
        else { t.clamp(t0, t1) }
    };
    newton_refine_pc(curve, &mut best_t, &mut best_dist, query, 30, clamp_t);

    let best_point = curve.point_at(best_t);
    CurveProjection { point: best_point, param: best_t, distance: (best_point - query).length() }
}

// =============================================================================
// Point-Surface Extrema (Extrema_ExtPS + Extrema_GenLocateExtPS)
// =============================================================================

/// OCCT-aligned: TreatSolution (Extrema_ExtPS).
/// Wraps periodic UV into [uinf, uinf+period), then checks bounds.
fn treat_solution(
    u: f64, v: f64, point: DVec3,
    uinf: f64, usup: f64,
    vinf: f64, vsup: f64,
    u_period: Option<f64>, v_period: Option<f64>,
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
    if (uinf - u) <= tolu && (u - usup) <= tolu && (vinf - v) <= tolv && (v - vsup) <= tolv {
        Some((u, v, point))
    } else {
        None
    }
}

/// OCCT-aligned: Extrema_GenLocateExtPS (local Newton from initial UV guess).
/// Solves (S-Q)·dS/du = 0 and (S-Q)·dS/dv = 0 via Gauss-Newton.
/// Returns the projected point and UV, or the initial guess if Newton fails.
pub fn closest_point_on_surface_near(
    surface: &Surface3,
    query: DVec3,
    u0: f64,
    v0: f64,
) -> SurfaceProjection {
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
        if nd < best_dist {
            best_dist = nd;
            u = nu;
            v = nv;
        }
        if du_step.abs() < 1e-10 && dv_step.abs() < 1e-10 { break; }
    }
    let point = surface.point_at(u, v);
    SurfaceProjection { point, params: (u, v), distance: (point - query).length() }
}

/// Numerical closest-point on a parametric surface via uniform sampling +
/// Newton refinement of f(u,v) = |S(u,v) - Q|².
/// OCCT-aligned: Extrema_ExtPS numeric mode (grid + Newton).
pub(crate) fn numeric_surface_projection(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    let [u0, u1, v0, v1] = surface.default_domain();
    // For unbounded surfaces (e.g. Plane, whose default domain is
    // [-inf,inf]x[-inf,inf]), fall back to a finite range so grid sampling and
    // Newton clamping are well-defined. Matches curve_domain below.
    let u0 = if u0.is_finite() { u0 } else { -1e6 };
    let u1 = if u1.is_finite() { u1 } else { 1e6 };
    let v0 = if v0.is_finite() { v0 } else { -1e6 };
    let v1 = if v1.is_finite() { v1 } else { 1e6 };

    let is_bspline = matches!(surface, Surface3::BSpline(_));
    let mut nu = if is_bspline { 44usize } else { 32usize }.max(n_samples);
    let mut nv = if is_bspline { 44usize } else { 32usize }.max(n_samples);

    // IsoIsDeg — detect degenerate isoparametric curves → increase samples
    let step_u = (u1 - u0) / 10.0;
    let step_v = (v1 - v0) / 10.0;
    if step_u.is_finite() && u1 > u0 {
        let d_max = (0..=10)
            .filter_map(|i| {
                let u = u0 + step_u * i as f64;
                if u < u0 || u > u1 { return None; }
                let (_p, du, _dv) = surface.derivatives(u, v0);
                Some(du.length_squared())
            })
            .fold(0.0_f64, f64::max);
        if d_max < 1e-18 || d_max > 1e9 { nu = nu.max(300); }
    }
    if step_v.is_finite() && v1 > v0 {
        let d_max = (0..=10)
            .filter_map(|i| {
                let v = v0 + step_v * i as f64;
                if v < v0 || v > v1 { return None; }
                let (_p, _du, dv) = surface.derivatives(u0, v);
                Some(dv.length_squared())
            })
            .fold(0.0_f64, f64::max);
        if d_max < 1e-18 || d_max > 1e9 { nv = nv.max(300); }
    }

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

    // Newton refinement using analytic derivatives
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
    SurfaceProjection { point: best_point, params: (best_u, best_v), distance: (best_point - query).length() }
}

// =============================================================================
// Curve-Curve Extrema (Extrema_ExtCC)
// =============================================================================

fn curve_domain(c: &Curve3) -> [f64; 2] {
    match c {
        Curve3::Line(_) => [-1e6, 1e6],
        other => other.default_domain(),
    }
}

#[inline]
fn sq_dist(c1: &Curve3, c2: &Curve3, s: f64, t: f64) -> f64 {
    (c1.point_at(s) - c2.point_at(t)).length_squared()
}

const H_CC: f64 = 1e-6; // finite-difference step

fn gradient(c1: &Curve3, c2: &Curve3, s: f64, t: f64) -> [f64; 2] {
    let p1 = c1.point_at(s);
    let p2 = c2.point_at(t);
    let diff = p1 - p2;
    let d1 = (c1.point_at(s + H_CC) - c1.point_at(s - H_CC)) / (2.0 * H_CC);
    let d2 = (c2.point_at(t + H_CC) - c2.point_at(t - H_CC)) / (2.0 * H_CC);
    [2.0 * diff.dot(d1), -2.0 * diff.dot(d2)]
}

fn hessian_diag(c1: &Curve3, c2: &Curve3, s: f64, t: f64) -> [f64; 2] {
    let d1 = (c1.point_at(s + H_CC) - c1.point_at(s - H_CC)) / (2.0 * H_CC);
    let d2 = (c2.point_at(t + H_CC) - c2.point_at(t - H_CC)) / (2.0 * H_CC);
    [
        2.0 * d1.length_squared().max(1e-30),
        2.0 * d2.length_squared().max(1e-30),
    ]
}

const MAX_ITER_CC: usize = 50;
const GRAD_TOL_CC: f64 = 1e-12;
const PARAM_TOL_CC: f64 = 1e-9;

fn newton_refine_cc(
    c1: &Curve3,
    c2: &Curve3,
    dom1: [f64; 2],
    dom2: [f64; 2],
    s0: f64,
    t0: f64,
) -> (f64, f64) {
    let mut s = s0;
    let mut t = t0;
    for _ in 0..MAX_ITER_CC {
        let g = gradient(c1, c2, s, t);
        if g[0].abs() < GRAD_TOL_CC && g[1].abs() < GRAD_TOL_CC { break; }
        let h = hessian_diag(c1, c2, s, t);
        let ds = -g[0] / h[0];
        let dt = -g[1] / h[1];
        let f0 = sq_dist(c1, c2, s, t);
        let mut alpha = 1.0;
        for _ in 0..8 {
            let ns = (s + alpha * ds).clamp(dom1[0], dom1[1]);
            let nt = (t + alpha * dt).clamp(dom2[0], dom2[1]);
            if sq_dist(c1, c2, ns, nt) < f0 {
                s = ns;
                t = nt;
                break;
            }
            alpha *= 0.5;
        }
        if (ds * alpha).abs() < PARAM_TOL_CC && (dt * alpha).abs() < PARAM_TOL_CC { break; }
    }
    (s, t)
}

/// Find all local minima of the curve-curve distance function.
///
/// `n_samples` controls the coarse grid density used to seed the Newton
/// refinement. A value of 16–32 is sufficient for most analytic curves.
/// For complex B-splines use 64+.
pub fn extrema_curve_curve(c1: &Curve3, c2: &Curve3, n_samples: usize) -> CurveCurveExtrema {
    let dom1 = curve_domain(c1);
    let dom2 = curve_domain(c2);
    let n = n_samples.max(2);

    let ss: Vec<f64> = (0..n)
        .map(|i| dom1[0] + (dom1[1] - dom1[0]) * i as f64 / (n - 1) as f64)
        .collect();
    let tt: Vec<f64> = (0..n)
        .map(|j| dom2[0] + (dom2[1] - dom2[0]) * j as f64 / (n - 1) as f64)
        .collect();

    let mut grid = vec![vec![0.0f64; n]; n];
    for (i, &s) in ss.iter().enumerate() {
        for (j, &t) in tt.iter().enumerate() {
            grid[i][j] = sq_dist(c1, c2, s, t);
        }
    }

    let mut seeds: Vec<(f64, f64)> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let v = grid[i][j];
            let is_min = (0usize..2)
                .flat_map(|di| (0usize..2).map(move |dj| (di, dj)))
                .all(|(di, dj)| {
                    let ni = i.wrapping_add(di).wrapping_sub(1);
                    let nj = j.wrapping_add(dj).wrapping_sub(1);
                    if ni == i && nj == j { return true; }
                    if ni >= n || nj >= n { return true; }
                    grid[ni][nj] >= v
                });
            if is_min { seeds.push((ss[i], tt[j])); }
        }
    }
    for &s in &[dom1[0], dom1[1]] {
        for &t in &[dom2[0], dom2[1]] {
            seeds.push((s, t));
        }
    }

    let mut pairs: Vec<ExtremaPair> = seeds
        .iter()
        .map(|&(s0, t0)| {
            let (s, t) = newton_refine_cc(c1, c2, dom1, dom2, s0, t0);
            let p1 = c1.point_at(s);
            let p2 = c2.point_at(t);
            ExtremaPair { param1: s, param2: t, point1: p1, point2: p2, distance: (p1 - p2).length() }
        })
        .collect();

    const DEDUP_TOL: f64 = 1e-4;
    pairs.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<ExtremaPair> = Vec::new();
    'outer: for p in pairs {
        for k in &kept {
            if (p.param1 - k.param1).abs() < DEDUP_TOL && (p.param2 - k.param2).abs() < DEDUP_TOL {
                continue 'outer;
            }
        }
        kept.push(p);
    }

    CurveCurveExtrema { pairs: kept }
}

// =============================================================================
// Extrema_ExtElC — line × elementary curve extrema (OCCT ExtElC 1:1)
// =============================================================================

use crate::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};

/// Frame of an elementary curve: x2/y2/z2 orthonormal + location O2.
struct ElCFrame {
    x2: DVec3,
    y2: DVec3,
    z2: DVec3,
    o2: DVec3,
}

impl ElCFrame {
    fn circle(c: &Circle3) -> Self {
        ElCFrame {
            x2: c.x_dir.normalize_or_zero(),
            y2: c.y_dir.normalize_or_zero(),
            z2: c.normal.normalize_or_zero(),
            o2: c.center,
        }
    }
    fn ellipse(e: &Ellipse3) -> Self {
        let x2 = e.major_dir.normalize_or_zero();
        let z2 = e.normal.normalize_or_zero();
        let y2 = z2.cross(x2).normalize_or_zero();
        ElCFrame { x2, y2, z2, o2: e.center }
    }
    fn hyperbola(h: &Hyperbola3) -> Self {
        let x2 = h.major_dir.normalize_or_zero();
        let z2 = h.normal.normalize_or_zero();
        let y2 = z2.cross(x2).normalize_or_zero();
        ElCFrame { x2, y2, z2, o2: h.center }
    }
    fn parabola(p: &Parabola3) -> Self {
        let x2 = p.axis_dir.normalize_or_zero();
        let z2 = p.normal.normalize_or_zero();
        let y2 = z2.cross(x2).normalize_or_zero();
        ElCFrame { x2, y2, z2, o2: p.vertex }
    }
}

/// Coordinates of the line in the conic frame (OCCT ExtElC: D and O2O1).
struct LineInFrame {
    dx: f64,
    dy: f64,
    dz: f64,
    ox: f64,
    oy: f64,
    oz: f64,
}

fn line_in_frame(line: &Line3, f: &ElCFrame) -> LineInFrame {
    let d = line.direction.normalize_or_zero();
    let o2o1 = line.origin - f.o2;
    LineInFrame {
        dx: d.dot(f.x2),
        dy: d.dot(f.y2),
        dz: d.dot(f.z2),
        ox: o2o1.dot(f.x2),
        oy: o2o1.dot(f.y2),
        oz: o2o1.dot(f.z2),
    }
}

/// OCCT `RefineDir` (Extrema_ExtElC.cxx): re-normalize a direction after
/// expressing it in a rotated frame.
fn refine_dir(v: DVec3) -> DVec3 {
    v.normalize_or_zero()
}

/// OCCT `Extrema_ExtElC(gp_Lin, gp_Circ)` (L471-623) + `PlanarLineCircleExtrema`
/// (L361-439). Returns (distance, u1_on_line, u2_on_circle) interior extrema.
pub fn line_circle_extrema(line: &Line3, circle: &Circle3) -> Vec<(f64, f64, f64)> {
    let f = ElCFrame::circle(circle);
    let lf = line_in_frame(line, &f);
    let r = circle.radius;
    // OCCT L365: if |aDirC.Dot(aDirL)| > Angular -> not planar, use 3D equation.
    if lf.dz.abs() <= 1e-12 {
        // PlanarLineCircleExtrema (L361-439): line parallel to circle plane.
        // 2D line: point (ox, oy), direction (dx, dy) in the circle frame.
        let plx = lf.ox;
        let ply = lf.oy;
        let dlx = lf.dx;
        let dly = lf.dy;
        let dl_sq = dlx * dlx + dly * dly;
        if dl_sq < 1e-30 {
            return vec![];
        }
        let dc = (plx * dly - ply * dlx).abs() / dl_sq.sqrt();
        let h = lf.oz; // plane offset (constant since line ∥ plane)
        let mut cands = Vec::new();
        // ExtElC2d line-circle: closest pair is at the foot.
        let t_foot = -(plx * dlx + ply * dly) / dl_sq;
        let foot = (plx + t_foot * dlx, ply + t_foot * dly);
        // IntAna2d_AnaIntersection (line-circle): intersections when dc <= R.
        if dc <= r {
            // 3D min distance = |h| at the intersection params.
            let s = (r * r - dc * dc).max(0.0).sqrt() / dl_sq.sqrt();
            for &sign in &[1.0, -1.0] {
                let t = t_foot + sign * s;
                let p2d = (plx + t * dlx, ply + t * dly);
                let u2 = p2d.1.atan2(p2d.0);
                cands.push((h.abs(), t, u2));
            }
        } else {
            let u2 = foot.1.atan2(foot.0);
            let dist = (h * h + (dc - r) * (dc - r)).sqrt();
            cands.push((dist, t_foot, u2));
        }
        cands
    } else {
        // Non-planar: trigonometric equation (OCCT L490-621).
        let d = refine_dir(DVec3::new(lf.dx, lf.dy, lf.dz));
        let o2o1 = DVec3::new(lf.ox, lf.oy, lf.oz);
        let v = d * o2o1.dot(d) - o2o1;
        // OCCT L556-585: coefficients (divided by R), zeroed at 1e-12.
        let dx = lf.dx;
        let dy = lf.dy;
        let a5 = r * dx * dy;
        let a1 = -2.0 * a5;
        let a2 = 0.5 * r * (dx * dx - dy * dy);
        let a3 = v.y;
        let a4 = -v.x;
        let mut coeff = [a1, a2, a3, a4, a5];
        for c in coeff.iter_mut() {
            if c.abs() <= 1e-12 {
                *c = 0.0;
            }
        }
        // ExtremaExtElC_TrigonometricRoots (L121-252) -> math_TrigonometricFunctionRoots.
        let res = crate::math::root::trig_function_roots(coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], 0.0, std::f64::consts::TAU);
        if !res.done || res.infinite {
            return vec![];
        }
        let d1 = line.direction.normalize_or_zero();
        let mut cands = Vec::new();
        for u2 in res.roots {
            let p2 = circle.point_at(u2);
            let u1 = (p2 - line.origin).dot(d1);
            let p1 = line.point_at(u1);
            let dist = (p1 - p2).length();
            cands.push((dist, u1, u2));
        }
        cands
    }
}

/// OCCT `Extrema_ExtElC(gp_Lin, gp_Elips)` (L627-753).
pub fn line_ellipse_extrema(line: &Line3, ell: &Ellipse3) -> Vec<(f64, f64, f64)> {
    let f = ElCFrame::ellipse(ell);
    let lf = line_in_frame(line, &f);
    let d = refine_dir(DVec3::new(lf.dx, lf.dy, lf.dz));
    let o2o1 = DVec3::new(lf.ox, lf.oy, lf.oz);
    let v = d * o2o1.dot(d) - o2o1;
    let dx = lf.dx;
    let dy = lf.dy;
    let maj_r = ell.major_radius;
    let min_r = ell.minor_radius;
    let r2 = maj_r * maj_r;
    let m2 = min_r * min_r;
    // OCCT L690-719
    let mut a5 = maj_r * min_r * dx * dy;
    let mut a1 = -2.0 * a5;
    let mut a2 = (r2 * dx * dx - m2 * dy * dy - r2 + m2) / 2.0;
    let mut a3 = min_r * v.y;
    let mut a4 = -maj_r * v.x;
    for c in [&mut a1, &mut a2, &mut a3, &mut a4, &mut a5] {
        if c.abs() <= 1e-12 {
            *c = 0.0;
        }
    }
    let res = crate::math::root::trig_function_roots(a1, a2, a3, a4, a5, 0.0, std::f64::consts::TAU);
    if !res.done || res.infinite {
        return vec![];
    }
    let d1 = line.direction.normalize_or_zero();
    let mut cands = Vec::new();
    for u2 in res.roots {
        let p2 = ell.point_at(u2);
        let u1 = (p2 - line.origin).dot(d1);
        let p1 = line.point_at(u1);
        cands.push(((p1 - p2).length(), u1, u2));
    }
    cands
}

/// OCCT `Extrema_ExtElC(gp_Lin, gp_Hypr)` (L757-858): quartic in v, u2 = ln(v).
pub fn line_hyperbola_extrema(line: &Line3, hyp: &Hyperbola3) -> Vec<(f64, f64, f64)> {
    let f = ElCFrame::hyperbola(hyp);
    let lf = line_in_frame(line, &f);
    let d = refine_dir(DVec3::new(lf.dx, lf.dy, lf.dz));
    let o2o1 = DVec3::new(lf.ox, lf.oy, lf.oz);
    let v_xyz = d * o2o1.dot(d) - o2o1;
    let vx = v_xyz.x;
    let vy = v_xyz.y;
    let dx = lf.dx;
    let dy = lf.dy;
    let r_maj = hyp.semi_major;
    let r_min = hyp.semi_minor;
    // OCCT L823-830
    let a = -2.0 * r_maj * r_min * dx * dy;
    let b = -r_maj * r_maj * dx * dx - r_min * r_min * dy * dy + r_maj * r_maj + r_min * r_min;
    let a1 = a + b;
    let a2 = 2.0 * r_maj * vx + 2.0 * r_min * vy;
    let a4 = -2.0 * r_maj * vx + 2.0 * r_min * vy;
    let a5 = a - b;
    // math_DirectPolynomialRoots (A1, A2, 0, A4, A5)
    let roots_v = crate::math::math_poly::solve_quartic(a1, a2, 0.0, a4, a5);
    let d1 = line.direction.normalize_or_zero();
    let mut cands = Vec::new();
    for v in roots_v {
        if v > 0.0 {
            let u2 = v.ln();
            let p2 = hyp.point_at(u2);
            let u1 = (p2 - line.origin).dot(d1);
            let p1 = line.point_at(u1);
            cands.push(((p1 - p2).length(), u1, u2));
        }
    }
    cands
}

/// OCCT `Extrema_ExtElC(gp_Lin, gp_Parab)` (L862-951): cubic in y.
pub fn line_parabola_extrema(line: &Line3, par: &Parabola3) -> Vec<(f64, f64, f64)> {
    let f = ElCFrame::parabola(par);
    let lf = line_in_frame(line, &f);
    let d = refine_dir(DVec3::new(lf.dx, lf.dy, lf.dz));
    let o2o1 = DVec3::new(lf.ox, lf.oy, lf.oz);
    let v_xyz = d * o2o1.dot(d) - o2o1;
    let dx = lf.dx;
    let dy = lf.dy;
    let p = par.focal_param;
    // OCCT L923-927
    let a1 = (1.0 - dx * dx) / (2.0 * p * p);
    let a2 = -3.0 * dx * dy / (2.0 * p);
    let a3 = 1.0 - dy * dy + v_xyz.x / p;
    let a4 = v_xyz.y;
    // math_DirectPolynomialRoots (A1, A2, A3, A4)
    let roots_y = crate::math::math_poly::solve_cubic(a1, a2, a3, a4);
    let d1 = line.direction.normalize_or_zero();
    let mut cands = Vec::new();
    for u2 in roots_y {
        let p2 = par.point_at(u2);
        let u1 = (p2 - line.origin).dot(d1);
        let p1 = line.point_at(u1);
        cands.push(((p1 - p2).length(), u1, u2));
    }
    cands
}

/// OCCT `Extrema_ExtElC(gp_Lin, gp_Lin)` (L268-357): interior closest pair.
pub fn line_line_extrema(l1: &Line3, l2: &Line3) -> Vec<(f64, f64, f64)> {
    let a_d1 = l1.direction.normalize_or_zero();
    let a_d2 = l2.direction.normalize_or_zero();
    let a_cos_a = a_d1.dot(a_d2);
    let a_sq_sin_a = 1.0 - a_cos_a * a_cos_a;
    let mut result = Vec::new();
    if a_sq_sin_a < 1e-30 || a_d1.cross(a_d2).length() < 1e-12 {
        // Parallel (OCCT L327-347): constant distance at any point of C1.
        // mySqDist[0] = C2.SquareDistance(C1.Location()) — one solution.
        let d = (l2.point_at(0.0) - l1.origin).length();
        result.push((d, 0.0, 0.0));
        return result;
    }
    // OCCT L333-336
    let a_l1l2 = l2.origin - l1.origin;
    let a_d1_l = a_d1.dot(a_l1l2);
    let a_d2_l = a_d2.dot(a_l1l2);
    let a_u1 = (a_d1_l - a_cos_a * a_d2_l) / a_sq_sin_a;
    let a_u2 = (a_cos_a * a_d1_l - a_d2_l) / a_sq_sin_a;
    let p1 = l1.point_at(a_u1);
    let p2 = l2.point_at(a_u2);
    result.push(((p1 - p2).length(), a_u1, a_u2));
    result
}

// =============================================================================
// Extrema_ExtCC — curve-curve extrema + range trimming (OCCT ExtCC 1:1)
// =============================================================================

/// OCCT `Extrema_ExtCC` trimmed square distances for the range corners
/// (mydist11/12/21/22, OCCT L214-245).
pub struct CornerDists {
    pub dist11: f64,
    pub dist12: f64,
    pub dist21: f64,
    pub dist22: f64,
}

/// OCCT `Extrema_ExtCC::Perform` (L177-317) + `PrepareResults` (L832-901)
/// for a line and an elementary curve. Returns the interior extrema whose
/// parameters fall inside the ranges, plus the corner distances.
pub struct ExtCCResult {
    /// Interior extrema (distance, u1, u2) already clipped to the ranges.
    pub interior: Vec<(f64, f64, f64)>,
    pub corners: CornerDists,
}

pub fn ext_cc_line_conic(
    line: &Line3,
    t1: f64,
    t2: f64,
    conic: &Curve3,
    u1: f64,
    u2: f64,
) -> ExtCCResult {
    // OCCT L214-245: corner distances from the 4 range endpoints.
    let p1f = line.point_at(t1);
    let p1l = line.point_at(t2);
    let p2f = conic.point_at(u1);
    let p2l = conic.point_at(u2);
    let corners = CornerDists {
        dist11: p1f.distance_squared(p2f),
        dist12: p1f.distance_squared(p2l),
        dist21: p1l.distance_squared(p2f),
        dist22: p1l.distance_squared(p2l),
    };
    // OCCT L247-294: dispatch to ExtElC for line-conic.
    let all = match conic {
        Curve3::Circle(c) => line_circle_extrema(line, c),
        Curve3::Ellipse(e) => line_ellipse_extrema(line, e),
        Curve3::Hyperbola(h) => line_hyperbola_extrema(line, h),
        Curve3::Parabola(p) => line_parabola_extrema(line, p),
        _ => vec![],
    };
    // OCCT PrepareResults (L832-898): keep extrema whose parameters are in range.
    let is_periodic = matches!(conic, Curve3::Circle(_) | Curve3::Ellipse(_));
    let mut interior = Vec::new();
    for (dist, u, u2c) in all {
        // Periodic wrapping of the conic parameter (OCCT L869-876).
        let u2w = if is_periodic {
            let period = std::f64::consts::TAU;
            let diff = u2c - u1;
            u1 + diff - period * (diff / period).floor()
        } else {
            u2c
        };
        // OCCT L878-879: within ranges (RealEpsilon() margin).
        if u >= t1 - f64::EPSILON && u <= t2 + f64::EPSILON
            && u2w >= u1 - f64::EPSILON && u2w <= u2 + f64::EPSILON
        {
            interior.push((dist, u, u2w));
        }
    }
    ExtCCResult { interior, corners }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle3, Line3};

    fn line(origin: DVec3, dir: DVec3) -> Curve3 {
        Curve3::Line(Line3 { origin, direction: dir.normalize() })
    }

    fn circle(center: DVec3, normal: DVec3, radius: f64) -> Curve3 {
        Curve3::Circle(Circle3::new(center, normal, radius))
    }

    // ── Point-curve extrema ─────────────────────────────────────────────────

    #[test]
    fn project_onto_circle_curve() {
        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
        let q = DVec3::new(2.0, 0.0, 0.0);
        let r = closest_point_on_curve(&circle, q, 64);
        assert!((r.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
        assert!((r.distance - 1.0).abs() < 1e-6);
    }

    #[test]
    fn project_onto_line_curve() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let q = DVec3::new(3.0, 4.0, 0.0);
        let r = closest_point_on_curve(&line, q, 32);
        let expected = DVec3::new(3.0, 0.0, 0.0);
        assert!((r.point - expected).length() < 1e-4);
        assert!((r.distance - 4.0).abs() < 1e-4);
    }

    #[test]
    fn project_onto_ellipse_curve_analytic() {
        let ellipse = Curve3::Ellipse(crate::geom::Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 3.0, minor_radius: 1.0,
        });
        let q = DVec3::new(5.0, 0.0, 0.0);
        let r = closest_point_on_curve(&ellipse, q, 64);
        assert!((r.point - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-5);
        assert!((r.distance - 2.0).abs() < 1e-5);
    }

    #[test]
    fn project_onto_line_curve_oblique() {
        let dir = DVec3::new(1.0, 1.0, 0.0).normalize();
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: dir });
        let q = DVec3::new(0.0, 1.0, 2.0);
        let r = closest_point_on_curve(&line, q, 32);
        let t = q.dot(dir);
        let expected = dir * t;
        assert!((r.point - expected).length() < 1e-9);
    }

    #[test]
    fn project_onto_partial_circle_arc() {
        let arc = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
        let q = DVec3::new(-2.0, 0.0, 0.0);
        let r = closest_point_on_curve(&arc, q, 64);
        assert!((r.point - DVec3::new(-1.0, 0.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn project_onto_distant_circle() {
        let circle = Curve3::Circle(Circle3::new(DVec3::new(2.0, 3.0, 0.0), DVec3::Z, 1.0));
        let q = DVec3::new(3.0, 3.0, 0.0);
        let r = closest_point_on_curve(&circle, q, 32);
        assert!((r.point - DVec3::new(3.0, 3.0, 0.0)).length() < 1e-6);
    }

    // ── Point-surface extrema ───────────────────────────────────────────────

    #[test]
    fn project_near_surface_boundary() {
        use crate::geom::Plane;
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let q = DVec3::new(1.0, 2.0, 1e-10);
        let r = numeric_surface_projection(&plane, q, 8);
        assert!(r.distance < 1e-9);
    }

    // ── Curve-curve extrema ─────────────────────────────────────────────────

    #[test]
    fn parallel_lines() {
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = line(DVec3::new(0.0, 3.0, 0.0), DVec3::X);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        assert!(!ex.pairs.is_empty());
        let d = ex.min_distance();
        assert!((d - 3.0).abs() < 0.01, "expected 3.0, got {d}");
    }

    #[test]
    fn skew_lines() {
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = line(DVec3::new(0.0, 0.0, 5.0), DVec3::Y);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        assert!((d - 5.0).abs() < 0.01, "expected 5.0, got {d}");
    }

    #[test]
    fn line_and_circle() {
        let c1 = line(DVec3::new(5.0, 0.0, 0.0), DVec3::Z);
        let c2 = circle(DVec3::ZERO, DVec3::Z, 2.0);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        assert!((d - 3.0).abs() < 0.05, "expected ~3.0, got {d}");
    }

    #[test]
    fn concentric_circles() {
        let c1 = circle(DVec3::ZERO, DVec3::Z, 2.0);
        let c2 = circle(DVec3::ZERO, DVec3::Z, 5.0);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        assert!((d - 3.0).abs() < 0.01, "expected 3.0, got {d}");
    }

    #[test]
    fn intersecting_lines_have_zero_distance() {
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = line(DVec3::ZERO, DVec3::Y);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        assert!(d < 0.01, "crossing lines should have distance ≈ 0, got {d}");
    }

    #[test]
    fn same_circle_has_zero_min_distance() {
        let c = circle(DVec3::ZERO, DVec3::Z, 3.0);
        let ex = extrema_curve_curve(&c, &c, 32);
        let d = ex.min_distance();
        assert!(d < 0.01, "same circle min distance should be 0, got {d}");
    }

    #[test]
    fn extrema_pairs_are_sorted_by_distance() {
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = circle(DVec3::ZERO, DVec3::Z, 3.0);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let distances: Vec<f64> = ex.pairs.iter().map(|p| p.distance).collect();
        for w in distances.windows(2) {
            assert!(w[0] <= w[1] + 1e-10, "pairs should be sorted ascending by distance");
        }
    }

    // ── Extrema_ExtElC (line × elementary curve) ────────────────────────────

    /// Min over the interior extrema (ignoring ranges).
    fn min_interior(cands: &[(f64, f64, f64)]) -> f64 {
        cands.iter().map(|&(d, _, _)| d).fold(f64::INFINITY, f64::min)
    }

    #[test]
    fn line_circle_far_away() {
        // line y=5 (x-axis at z=0), circle radius 1 at origin -> min dist 4.
        let l = Curve3::Line(crate::geom::Line3 { origin: DVec3::new(0.0, 5.0, 0.0), direction: DVec3::X });
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
        let cands = match &c {
            Curve3::Circle(cc) => line_circle_extrema(match &l { Curve3::Line(ll) => ll, _ => unreachable!() }, cc),
            _ => unreachable!(),
        };
        let d = min_interior(&cands);
        assert!((d - 4.0).abs() < 1e-6, "expected 4.0, got {d}");
    }

    #[test]
    fn line_circle_intersecting() {
        // line y=0.5 passes through circle radius 1 at origin -> min dist 0.
        let l = Curve3::Line(crate::geom::Line3 { origin: DVec3::new(0.0, 0.5, 0.0), direction: DVec3::X });
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
        let cands = match &c {
            Curve3::Circle(cc) => line_circle_extrema(match &l { Curve3::Line(ll) => ll, _ => unreachable!() }, cc),
            _ => unreachable!(),
        };
        let d = min_interior(&cands);
        assert!(d < 1e-6, "expected 0, got {d}");
    }

    #[test]
    fn line_circle_off_plane() {
        // line z=3, x-axis; circle radius 1 at origin in xy-plane.
        // Line is parallel to circle plane at offset 3 -> min dist = 3 (dc2d=0 <= R).
        let l = Curve3::Line(crate::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        // shift the circle to z=3 so the line is at offset -3
        let c = Curve3::Circle(Circle3::new(DVec3::new(0.0, 0.0, 3.0), DVec3::Z, 1.0));
        let cands = match &c {
            Curve3::Circle(cc) => line_circle_extrema(match &l { Curve3::Line(ll) => ll, _ => unreachable!() }, cc),
            _ => unreachable!(),
        };
        let d = min_interior(&cands);
        assert!((d - 3.0).abs() < 1e-6, "expected 3.0, got {d}");
    }

    #[test]
    fn line_ellipse_extrema_smoke() {
        // line x=4 (direction +y); ellipse major 3 (x), minor 1 (y) at origin.
        // Nearest ellipse point is (3, 0) -> dist 1.
        let l = Curve3::Line(crate::geom::Line3 { origin: DVec3::new(4.0, 0.0, 0.0), direction: DVec3::Y });
        let e = Curve3::Ellipse(crate::geom::Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 3.0, minor_radius: 1.0,
        });
        let cands = match &e {
            Curve3::Ellipse(ee) => line_ellipse_extrema(match &l { Curve3::Line(ll) => ll, _ => unreachable!() }, ee),
            _ => unreachable!(),
        };
        let d = min_interior(&cands);
        assert!((d - 1.0).abs() < 1e-6, "expected 1.0, got {d}");
    }
}
