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

use crate::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
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
fn solve_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Vec<f64> {
    if a.abs() < 1e-30 {
        return vec![];
    }
    let inv_a = 1.0 / a;
    let ba = b * inv_a;
    let ca = c * inv_a;
    let da = d * inv_a;
    let ea = e * inv_a;
    let p = ca - 3.0 * ba * ba / 8.0;
    let q = da - ba * ca / 2.0 + ba * ba * ba / 8.0;
    let r = ea - ba * da / 4.0 + ba * ba * ca / 16.0 - 3.0 * ba * ba * ba * ba / 256.0;
    if q.abs() < 1e-30 {
        let disc = p * p - 4.0 * r;
        if disc < 0.0 {
            return vec![];
        }
        let sd = disc.sqrt();
        let mut roots = Vec::new();
        for &t in &[(-p + sd) / 2.0, (-p - sd) / 2.0] {
            if t >= 0.0 {
                roots.push(t.sqrt());
            }
            if t > 0.0 {
                roots.push(-t.sqrt());
            }
        }
        let shift = -ba / 4.0;
        return roots.into_iter().map(|x| x + shift).collect();
    }
    let rc = solve_cubic(2.0 * p, p * p - 4.0 * r, -q * q);
    let m = rc.into_iter().find(|&m| m > 0.0).unwrap_or(0.0);
    if m <= 0.0 {
        return vec![];
    }
    let sq = (m * 2.0).sqrt();
    let t1 = -p - m;
    let t2 = q / sq;
    let disc1 = -t1 - 2.0 * t2;
    let disc2 = -t1 + 2.0 * t2;
    let shift = -ba / 4.0;
    let mut roots = Vec::new();
    if disc1 >= 0.0 {
        let s = disc1.sqrt();
        roots.push((sq + s) / 2.0 + shift);
        roots.push((-sq + s) / 2.0 + shift);
    }
    if disc2 >= 0.0 {
        let s = disc2.sqrt();
        roots.push((s - sq) / 2.0 + shift);
        roots.push((-sq - s) / 2.0 + shift);
    }
    roots
}

/// Project the point `query` onto `curve`, returning the nearest point on the
/// curve, its parameter value, and the Euclidean distance.
///
/// OCCT-aligned: dispatches per-type matching Extrema_ExtPC:
///   - Line/Circle: analytic via Extrema_ExtPElC equivalent
///   - Ellipse: analytic init + Newton refinement
///   - BSpline: C2 interval splitting (Extrema_GGExtPC)
///   - Bezier/Other: uniform sampling + Newton
///
/// `n_samples` is the uniform sampling count (used for Bezier and fallback;
/// for BSpline it is overridden by `degree + 1` per C2 interval).
pub fn closest_point_on_curve(curve: &Curve3, query: DVec3, n_samples: usize) -> CurveProjection {
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

        // OCCT-aligned: Extrema_ExtPElC::Perform(Hyperbola)
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
                if v > 0.0 {
                    let u = v.ln();
                    if u >= t0 && u <= t1 {
                        let pt = hyp.point_at(u);
                        let d = (pt - query).length();
                        if d < d_best { d_best = d; u_best = u; }
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

        // OCCT-aligned: Extrema_ExtPElC::Perform(Parabola)
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
            for u in roots {
                if u >= t0 && u <= t1 {
                    let pt = par.point_at(u);
                    let d = (pt - query).length();
                    if d < d_best { d_best = d; u_best = u; }
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
}
