//! Curve-surface intersection (GeomAPI_IntCS).
//!
//! OCCT TKGeomBase GeomAPI package: GeomAPI_IntCS.
//!
//! Computes intersection points between a 3D curve and a surface.
//! Algorithm (OCCT IntCurveSurface_HInter):
//! 1. Sample the curve at N points
//! 2. For each point compute signed distance to surface via ExtPS
//! 3. Detect sign changes → seeds at curve-surface crossings
//! 4. Newton refinement at each crossing: solve f(t,u,v) = P(t) - S(u,v) = 0

#![allow(clippy::manual_clamp)]

use glam::DVec3;

use crate::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use crate::base::extrema::ExtPS;

const TOL: f64 = 1e-7;

/// Result of a curve-surface intersection point.
#[derive(Debug, Clone)]
pub struct IntCSPoint {
    /// The intersection point in 3D.
    pub point: DVec3,
    /// Parameter on the curve.
    pub w: f64,
    /// Parameters on the surface.
    pub u: f64,
    /// Parameters on the surface.
    pub v: f64,
}

/// Curve-surface intersection algorithm.
///
/// OCCT: `GeomAPI_IntCS`.
pub struct IntCS {
    done: bool,
    points: Vec<IntCSPoint>,
    segments: Vec<Curve3>,
    seg_params: Vec<(f64, f64, f64, f64)>, // (U1,V1, U2,V2) for each segment
}

impl IntCS {
    /// Default constructor.
    ///
    /// OCCT: `GeomAPI_IntCS()`.
    pub fn new() -> Self {
        IntCS {
            done: false,
            points: Vec::new(),
            segments: Vec::new(),
            seg_params: Vec::new(),
        }
    }

    /// Constructor with curve and surface.
    ///
    /// OCCT: `GeomAPI_IntCS(Curve, Surface)`.
    pub fn with_curve_surface(curve: &Curve3, surface: &Surface3) -> Self {
        let mut intcs = IntCS::new();
        intcs.perform(curve, surface);
        intcs
    }

    /// Perform the intersection.
    ///
    /// OCCT: `Perform(Curve, Surface)`.
    pub fn perform(&mut self, curve: &Curve3, surface: &Surface3) {
        self.points.clear();
        self.segments.clear();
        self.seg_params.clear();

        let dom = curve.default_domain();
        let t_min = dom[0];
        let t_max = dom[1];

        if !t_min.is_finite() || !t_max.is_finite() {
            self.done = true;
            return;
        }

        let range = t_max - t_min;
        if range < TOL {
            self.done = true;
            return;
        }

        let surf_dom = surface.default_domain();

        // 1. Sample the curve coarse grid
        const N_SAMPLES: usize = 101;
        let mut dists: Vec<(f64, f64)> = Vec::with_capacity(N_SAMPLES + 1); // (t, distance)

        for i in 0..=N_SAMPLES {
            let t = t_min + range * (i as f64) / (N_SAMPLES as f64);
            let p3d = curve.point_at(t);
            let dist = signed_distance_to_surface(p3d, surface, surf_dom[0], surf_dom[1], surf_dom[2], surf_dom[3]);
            dists.push((t, dist));
        }

        // 2. Find sign changes → seeds for Newton refinement
        let mut seeds: Vec<f64> = Vec::new();
        for i in 0..N_SAMPLES {
            let (t1, d1) = dists[i];
            let (t2, d2) = dists[i + 1];
            if d1 * d2 < 0.0 || d1.abs() < TOL {
                // Sign change or already near zero
                let t_seed = if d1.abs() < TOL {
                    t1
                } else if d2.abs() < TOL {
                    t2
                } else {
                    // Linear interpolation for the zero crossing
                    let frac = d1.abs() / (d1.abs() + d2.abs());
                    t1 + (t2 - t1) * frac
                };
                seeds.push(t_seed);
            }
        }

        // Deduplicate close seeds
        seeds.dedup_by(|a, b| (*a - *b).abs() < range / (N_SAMPLES as f64) * 0.5);

        // 3. Newton refinement at each seed
        for &t0 in &seeds {
            if let Some(pt) = refine_intersection(curve, surface, t0, t_min, t_max, surf_dom[0], surf_dom[1], surf_dom[2], surf_dom[3]) {
                // Deduplicate against already-found points
                let is_dup = self.points.iter().any(|existing| {
                    (existing.point - pt.point).length() < TOL * 10.0
                });
                if !is_dup {
                    self.points.push(pt);
                }
            }
        }

        // Sort by curve parameter
        self.points.sort_by(|a, b| a.w.partial_cmp(&b.w).unwrap());

        self.done = true;
    }

    /// Returns true if the intersection was computed.
    ///
    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Returns the number of intersection points.
    ///
    /// OCCT: `NbPoints()`.
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// Returns the intersection point at 1-based index.
    ///
    /// OCCT: `Point(Index)`.
    pub fn point(&self, index: usize) -> &IntCSPoint {
        assert!(index >= 1 && index <= self.points.len(), "IntCS: index out of range");
        &self.points[index - 1]
    }

    /// Returns parameters (U, V, W) for the point at 1-based index.
    ///
    /// OCCT: `Parameters(Index, U, V, W)`.
    pub fn parameters(&self, index: usize) -> (f64, f64, f64) {
        let pt = self.point(index);
        (pt.u, pt.v, pt.w)
    }

    /// Returns the number of intersection segments (tangential case).
    ///
    /// OCCT: `NbSegments()`.
    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }

    /// Returns the intersection segment at 1-based index.
    ///
    /// OCCT: `Segment(Index)`.
    pub fn segment(&self, index: usize) -> &Curve3 {
        assert!(index >= 1 && index <= self.segments.len(), "IntCS: segment index out of range");
        &self.segments[index - 1]
    }
}

impl Default for IntCS {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute signed distance from a point to a surface using ExtPS.
/// Returns negative when the point is "below" the surface (opposite normal direction).
fn signed_distance_to_surface(
    point: DVec3,
    surface: &Surface3,
    ufirst: f64,
    ulast: f64,
    vfirst: f64,
    vlast: f64,
) -> f64 {
    let ext = ExtPS::with_domain(point, surface, ufirst, ulast, vfirst, vlast, TOL, TOL);
    if ext.nb_ext() > 0 {
        let p_on_surf = ext.point(1);
        let d = p_on_surf.point - point;
        let abs_dist = d.length();

        // Compute surface normal at the closest point for sign
        let n = surface.normal_at(p_on_surf.u, p_on_surf.v);
        let sign = n.dot(d).signum();
        abs_dist * sign
    } else {
        // No projection found — use center of domain as fallback
        let u_mid = (ufirst + ulast) * 0.5;
        let v_mid = (vfirst + vlast) * 0.5;
        let surf_pt = surface.point_at(u_mid, v_mid);
        let n = surface.normal_at(u_mid, v_mid);
        let d_vec = surf_pt - point;
        d_vec.length() * n.dot(d_vec).signum()
    }
}

/// Newton refinement: solve P(t) - S(u,v) = 0 for (t, u, v).
fn refine_intersection(
    curve: &Curve3,
    surface: &Surface3,
    t0: f64,
    t_min: f64,
    t_max: f64,
    ufirst: f64,
    ulast: f64,
    vfirst: f64,
    vlast: f64,
) -> Option<IntCSPoint> {
    let mut t = t0.clamp(t_min, t_max);
    let mut u = (ufirst + ulast) * 0.5;
    let mut v = (vfirst + vlast) * 0.5;

    for _ in 0..30 {
        let p_curve = curve.point_at(t);
        let dp_curve = curve.derivative_at(t);
        let (p_surf, ps_u, ps_v) = surface.derivatives(u, v);

        let f = p_curve - p_surf;

        // 3x2 Jacobian: [dp_curve, -ps_u, -ps_v]
        // Solve for (dt, du, dv): J * [dt, du, dv]^T = -f
        // Using two-step: project onto surface tangent plane

        let fu = -ps_u;
        let fv = -ps_v;

        // Solve normal equation: J^T * J * d = -J^T * f
        // This is a 3x3 system
        let jtj_00 = dp_curve.dot(dp_curve);
        let jtj_01 = dp_curve.dot(fu);
        let jtj_02 = dp_curve.dot(fv);
        let jtj_11 = fu.dot(fu);
        let jtj_12 = fu.dot(fv);
        let jtj_22 = fv.dot(fv);

        let rhs_0 = -(f.dot(dp_curve));
        let rhs_1 = -(f.dot(fu));
        let rhs_2 = -(f.dot(fv));

        // Solve 3x3 via Cramer's rule — the matrix is symmetric
        let det = jtj_00 * (jtj_11 * jtj_22 - jtj_12 * jtj_12)
            - jtj_01 * (jtj_01 * jtj_22 - jtj_12 * jtj_02)
            + jtj_02 * (jtj_01 * jtj_12 - jtj_11 * jtj_02);

        if det.abs() < 1e-30 {
            break;
        }

        let inv_det = 1.0 / det;
        let dt = inv_det * (
            rhs_0 * (jtj_11 * jtj_22 - jtj_12 * jtj_12)
            - rhs_1 * (jtj_01 * jtj_22 - jtj_12 * jtj_02)
            + rhs_2 * (jtj_01 * jtj_12 - jtj_11 * jtj_02)
        );
        let du = inv_det * (
            -rhs_0 * (jtj_01 * jtj_22 - jtj_12 * jtj_02)
            + rhs_1 * (jtj_00 * jtj_22 - jtj_02 * jtj_02)
            - rhs_2 * (jtj_00 * jtj_12 - jtj_01 * jtj_02)
        );
        let dv = inv_det * (
            rhs_0 * (jtj_01 * jtj_12 - jtj_11 * jtj_02)
            - rhs_1 * (jtj_00 * jtj_12 - jtj_01 * jtj_02)
            + rhs_2 * (jtj_00 * jtj_11 - jtj_01 * jtj_01)
        );

        let new_t = (t + dt).clamp(t_min, t_max);
        let new_u = (u + du).clamp(ufirst, ulast);
        let new_v = (v + dv).clamp(vfirst, vlast);

        // Check convergence
        let err = (p_curve - p_surf).length();
        if err < TOL {
            t = new_t;
            u = new_u;
            v = new_v;
            break;
        }

        t = new_t;
        u = new_u;
        v = new_v;

        if dt.abs() < TOL && du.abs() < TOL && dv.abs() < TOL {
            break;
        }
    }

    // Verify the final result
    let p_curve = curve.point_at(t);
    let p_surf = surface.point_at(u, v);
    let err = (p_curve - p_surf).length();
    if err < TOL * 100.0 {
        Some(IntCSPoint {
            point: (p_curve + p_surf) * 0.5,
            w: t,
            u,
            v,
        })
    } else {
        None
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn test_line_intersects_plane() {
        // Line along Z from (0,0,-5) to (0,0,5), plane at z=0
        let line = Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, -5.0), DVec3::Z));
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let intcs = IntCS::with_curve_surface(&line, &plane);
        assert!(intcs.is_done());
        assert_eq!(intcs.nb_points(), 1);
        let pt = intcs.point(1);
        assert!((pt.point - DVec3::ZERO).length() < 1e-6);
    }

    #[test]
    fn test_line_misses_plane() {
        // Line along Y at z=5, plane at z=0
        let line = Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, 5.0), DVec3::Y));
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let intcs = IntCS::with_curve_surface(&line, &plane);
        assert!(intcs.is_done());
        assert_eq!(intcs.nb_points(), 0);
    }
}
