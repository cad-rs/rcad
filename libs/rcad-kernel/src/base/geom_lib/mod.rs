//! Geometry utility library (GeomLib).
//!
//! OCCT TKGeomBase GeomLib package: GeomLib_Tool, GeomLib_IsPlanarSurface,
//! GeomLib_CheckCurveOnSurface, GeomLib::Inertia, GeomLib::AxeOfInertia.

#![allow(clippy::manual_clamp)]

use glam::{DVec2, DVec3};

use crate::geom::{Curve2dEval, Curve3, CurveEval, Plane, Surface3, SurfaceEval};
use crate::math::gp::Ax2;
use crate::math::math_jacobi::MathJacobi;
use crate::math::MatD;

// OCCT Precision::Confusion()
const TOL_CONF: f64 = 1e-7;
// OCCT Precision::PConfusion()
const TOL_PCONF: f64 = 1e-12;

// ============================================================================
// GeomLib_Tool
// ============================================================================

/// Static methods for parameter extraction and curve deviation.
///
/// OCCT: `GeomLib_Tool`.
pub struct Tool;

impl Tool {
    /// Extract the parameter of a 3D point on a 3D curve within MaxDist.
    ///
    /// OCCT: `GeomLib_Tool::Parameter(Curve, Point, MaxDist, U)`.
    /// Returns `Some(u)` if the point is within MaxDist of the curve.
    pub fn parameter_curve(curve: &Curve3, point: DVec3, max_dist: f64) -> Option<f64> {
        // OCCT uses Newton iteration on the curve. We sample a coarse grid
        // and refine with golden-section search.
        let domain = curve.default_domain();
        let t_min = domain[0];
        let t_max = domain[1];

        if !t_min.is_finite() || !t_max.is_finite() {
            // Infinite domain (line): sample at a few candidate points
            let candidates = [
                Tool::param_at_dist(curve, point, 0.0),
                Tool::param_at_dist(curve, point, 1.0),
                Tool::param_at_dist(curve, point, -1.0),
            ];
            let mut best = candidates[0];
            let mut best_d = (curve.point_at(best) - point).length();
            for &c in &candidates[1..] {
                let d = (curve.point_at(c) - point).length();
                if d < best_d {
                    best_d = d;
                    best = c;
                }
            }
            return if best_d < max_dist { Some(best) } else { None };
        }

        // Finite domain: coarse grid then local refinement
        const N_SAMPLES: usize = 101;
        let mut best_t = t_min;
        let mut best_d = (curve.point_at(t_min) - point).length();

        for i in 1..=N_SAMPLES {
            let t = t_min + (t_max - t_min) * (i as f64) / (N_SAMPLES as f64);
            let d = (curve.point_at(t) - point).length();
            if d < best_d {
                best_d = d;
                best_t = t;
            }
        }

        if best_d >= max_dist {
            return None;
        }

        // Refine with Newton iteration
        let mut t = best_t;
        for _ in 0..10 {
            let p = curve.point_at(t);
            let dp = curve.derivative_at(t);
            let d = p - point;
            let err = d.dot(dp);
            let speed_sq = dp.length_squared();
            if speed_sq < TOL_PCONF {
                break;
            }
            let dt = -err / speed_sq;
            t = (t + dt).clamp(t_min, t_max);
            if dt.abs() < TOL_PCONF {
                break;
            }
        }

        let final_d = (curve.point_at(t) - point).length();
        if final_d < max_dist { Some(t) } else { None }
    }

    /// Extract the (u, v) parameters of a 3D point on a surface within MaxDist.
    ///
    /// OCCT: `GeomLib_Tool::Parameters(Surface, Point, MaxDist, U, V)`.
    pub fn parameters_surface(surface: &Surface3, point: DVec3, max_dist: f64) -> Option<(f64, f64)> {
        let domain = surface.default_domain();
        let (u_min, u_max, v_min, v_max) = (domain[0], domain[1], domain[2], domain[3]);

        if !u_min.is_finite() || !v_min.is_finite() {
            // Infinite domain: sample near the closest point estimate
            let (u0, v0) = Tool::estimate_uv(surface, point);
            return Tool::refine_uv(surface, point, u0, v0, u_min, u_max, v_min, v_max, max_dist);
        }

        // Coarse grid
        const N_U: usize = 21;
        const N_V: usize = 21;
        let mut best_u = u_min;
        let mut best_v = v_min;
        let mut best_d = (surface.point_at(u_min, v_min) - point).length();

        for i in 0..=N_U {
            let u = u_min + (u_max - u_min) * (i as f64) / (N_U as f64);
            for j in 0..=N_V {
                let v = v_min + (v_max - v_min) * (j as f64) / (N_V as f64);
                let d = (surface.point_at(u, v) - point).length();
                if d < best_d {
                    best_d = d;
                    best_u = u;
                    best_v = v;
                }
            }
        }

        if best_d >= max_dist {
            return None;
        }

        // Refine with Newton
        let (u, v) = Tool::newton_refine(surface, point, best_u, best_v, u_min, u_max, v_min, v_max);
        let final_d = (surface.point_at(u, v) - point).length();
        if final_d < max_dist { Some((u, v)) } else { None }
    }

    /// Estimate UV near a point by projecting onto the surface's closest analytic form.
    fn estimate_uv(surface: &Surface3, point: DVec3) -> (f64, f64) {
        match surface {
            Surface3::Plane(p) => {
                let d = point - p.origin;
                (d.dot(p.u_dir), d.dot(p.v_dir))
            }
            Surface3::Cylinder(c) => {
                let d = point - c.origin;
                let along = d.dot(c.axis);
                let radial = d - c.axis * along;
                let u = radial.dot(c.ref_dir).atan2(
                    c.axis.cross(c.ref_dir).normalize_or_zero().dot(radial),
                );
                (u, along)
            }
            _ => (0.0, 0.0),
        }
    }

    fn refine_uv(
        surface: &Surface3,
        point: DVec3,
        u0: f64, v0: f64,
        _u_min: f64, _u_max: f64, _v_min: f64, _v_max: f64,
        max_dist: f64,
    ) -> Option<(f64, f64)> {
        let (u, v) = Tool::newton_refine(surface, point, u0, v0, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY);
        let d = (surface.point_at(u, v) - point).length();
        if d < max_dist { Some((u, v)) } else { None }
    }

    fn newton_refine(
        surface: &Surface3,
        point: DVec3,
        u0: f64, v0: f64,
        u_min: f64, u_max: f64,
        v_min: f64, v_max: f64,
    ) -> (f64, f64) {
        let mut u = u0;
        let mut v = v0;
        for _ in 0..15 {
            let (p, pu, pv) = surface.derivatives(u, v);
            let d = p - point;
            let err_u = d.dot(pu);
            let err_v = d.dot(pv);
            let puu = pu.dot(pu);
            let puv = pu.dot(pv);
            let pvv = pv.dot(pv);
            let det = puu * pvv - puv * puv;
            if det.abs() < TOL_PCONF {
                break;
            }
            let du = (-err_u * pvv + err_v * puv) / det;
            let dv = (-err_v * puu + err_u * puv) / det;
            u = (u + du).clamp(u_min, u_max);
            v = (v + dv).clamp(v_min, v_max);
            if du.abs() < TOL_PCONF && dv.abs() < TOL_PCONF {
                break;
            }
        }
        (u, v)
    }

    fn param_at_dist(curve: &Curve3, point: DVec3, t0: f64) -> f64 {
        let mut t = t0;
        for _ in 0..5 {
            let dp = curve.derivative_at(t);
            let d = curve.point_at(t) - point;
            let speed_sq = dp.length_squared();
            if speed_sq < TOL_PCONF {
                break;
            }
            t -= d.dot(dp) / speed_sq;
        }
        t
    }

    /// Extract the 2D parameter on a 2D curve within MaxDist.
    ///
    /// OCCT: `GeomLib_Tool::Parameter(Curve2d, Point2d, MaxDist, U)`.
    pub fn parameter_curve2d(curve: &crate::geom::Curve2d, point: DVec2, max_dist: f64) -> Option<f64> {
        let domain = curve.default_domain();
        let t_min = domain[0];
        let t_max = domain[1];

        if !t_min.is_finite() || !t_max.is_finite() {
            let t = Tool::param_at_dist_2d(curve, point, 0.0);
            let d = (curve.point_at(t) - point).length();
            return if d < max_dist { Some(t) } else { None };
        }

        const N_SAMPLES: usize = 101;
        let mut best_t = t_min;
        let mut best_d = (curve.point_at(t_min) - point).length();

        for i in 1..=N_SAMPLES {
            let t = t_min + (t_max - t_min) * (i as f64) / (N_SAMPLES as f64);
            let d = (curve.point_at(t) - point).length();
            if d < best_d {
                best_d = d;
                best_t = t;
            }
        }

        if best_d >= max_dist {
            return None;
        }

        let mut t = best_t;
        for _ in 0..10 {
            let dp = curve.derivative_at(t);
            let d = curve.point_at(t) - point;
            let speed_sq = dp.length_squared();
            if speed_sq < TOL_PCONF {
                break;
            }
            let dt = -d.dot(dp) / speed_sq;
            t = (t + dt).clamp(t_min, t_max);
            if dt.abs() < TOL_PCONF {
                break;
            }
        }

        let final_d = (curve.point_at(t) - point).length();
        if final_d < max_dist { Some(t) } else { None }
    }

    fn param_at_dist_2d(curve: &crate::geom::Curve2d, point: DVec2, t0: f64) -> f64 {
        let mut t = t0;
        for _ in 0..5 {
            let dp = curve.derivative_at(t);
            let d = curve.point_at(t) - point;
            let speed_sq = dp.length_squared();
            if speed_sq < TOL_PCONF {
                break;
            }
            t -= d.dot(dp) / speed_sq;
        }
        t
    }
}

// ============================================================================
// GeomLib_IsPlanarSurface
// ============================================================================

/// Check if a surface is planar within a given tolerance.
///
/// OCCT: `GeomLib_IsPlanarSurface`.
pub struct IsPlanarSurface {
    plane: Plane,
    is_planar: bool,
}

impl IsPlanarSurface {
    /// Construct and perform the check.
    ///
    /// OCCT: `GeomLib_IsPlanarSurface(Surface, Tol)`.
    pub fn new(surface: &Surface3, tol: f64) -> Self {
        let tol = if tol <= 0.0 { TOL_CONF } else { tol };
        let mut result = IsPlanarSurface {
            plane: Plane::new(DVec3::ZERO, DVec3::Z),
            is_planar: false,
        };

        match surface {
            Surface3::Plane(p) => {
                result.plane = *p;
                result.is_planar = true;
            }
            Surface3::BSpline(bsp) => {
                if crate::geom::bspline_is_planar(bsp, tol) {
                    result.plane = crate::geom::bspline_to_plane(bsp);
                    result.is_planar = true;
                }
            }
            Surface3::Bezier(bz) => {
                let pts: Vec<DVec3> = bz.control_points.iter().flat_map(|row| row.iter()).copied().collect();
                if let Some(pl) = fit_plane_to_points(&pts, tol) {
                    result.plane = pl;
                    result.is_planar = true;
                }
            }
            _ => {
                // For other surface types, sample the surface and check planarity
                let domain = surface.default_domain();
                let (u_min, u_max, v_min, v_max) = (domain[0], domain[1], domain[2], domain[3]);
                if !u_min.is_finite() || !v_min.is_finite() {
                    return result;
                }

                let mut pts = Vec::new();
                for i in 0..5 {
                    let u = u_min + (u_max - u_min) * (i as f64) / 4.0;
                    for j in 0..5 {
                        let v = v_min + (v_max - v_min) * (j as f64) / 4.0;
                        pts.push(surface.point_at(u, v));
                    }
                }

                if let Some(pl) = fit_plane_to_points(&pts, tol) {
                    // Verify all 25 points
                    let all_planar = pts.iter().all(|&p| {
                        let d = (p - pl.origin).dot(pl.normal);
                        d.abs() < tol
                    });
                    if all_planar {
                        result.plane = pl;
                        result.is_planar = true;
                    }
                }
            }
        }

        result
    }

    /// Returns true if the surface is planar.
    ///
    /// OCCT: `IsPlanar()`.
    pub fn is_planar(&self) -> bool {
        self.is_planar
    }

    /// Returns the best-fit plane.
    ///
    /// OCCT: `Plan()`.
    pub fn plan(&self) -> &Plane {
        &self.plane
    }
}

/// Fit a plane to a set of points within tolerance.
fn fit_plane_to_points(pts: &[DVec3], tol: f64) -> Option<Plane> {
    if pts.len() < 3 {
        return None;
    }

    let centroid = pts.iter().copied().sum::<DVec3>() / (pts.len() as f64);

    // Find normal via covariance
    let mut cov = DVec3::ZERO;
    let mut xx = 0.0;
    let mut yy = 0.0;
    let mut zz = 0.0;
    let mut xy = 0.0;
    let mut xz = 0.0;
    let mut yz = 0.0;

    for &p in pts {
        let d = p - centroid;
        xx += d.x * d.x;
        yy += d.y * d.y;
        zz += d.z * d.z;
        xy += d.x * d.y;
        xz += d.x * d.z;
        yz += d.y * d.z;
    }

    // Smallest eigenvector of the covariance matrix
    // Using the formula for a 3x3 symmetric matrix
    let normal = {
        // Use cross product of two largest spread directions as a fallback
        let d1 = pts[1] - pts[0];
        let d2 = pts[2] - pts[0];
        let n = d1.cross(d2);
        if n.length_squared() > tol * tol {
            n.normalize()
        } else {
            return None;
        }
    };

    let plane = Plane::new(centroid, normal);

    // Verify all points
    let max_dev = pts.iter().map(|&p| (p - centroid).dot(normal).abs()).fold(0.0, f64::max);
    if max_dev < tol {
        Some(plane)
    } else {
        None
    }
}

// ============================================================================
// GeomLib_CheckCurveOnSurface
// ============================================================================

/// Compute the max distance between a 3D curve and its 2D representation on a surface.
///
/// OCCT: `GeomLib_CheckCurveOnSurface`.
pub struct CheckCurveOnSurface {
    max_distance: f64,
    max_parameter: f64,
    error_status: i32,
    tol_range: f64,
}

impl CheckCurveOnSurface {
    /// Default constructor.
    ///
    /// OCCT: default constructor.
    pub fn new() -> Self {
        CheckCurveOnSurface {
            max_distance: 0.0,
            max_parameter: 0.0,
            error_status: 1,
            tol_range: TOL_PCONF,
        }
    }

    /// Constructor with curve data.
    ///
    /// OCCT: `CheckCurveOnSurface(Curve, TolRange)`.
    pub fn with_curve(_curve: &Curve3, tol_range: f64) -> Self {
        CheckCurveOnSurface {
            max_distance: 0.0,
            max_parameter: 0.0,
            error_status: 0,
            tol_range: if tol_range <= 0.0 { TOL_PCONF } else { tol_range },
        }
    }

    /// Perform the check: compute max distance between `curve_3d` and
    /// the surface evaluation of its pcurve parametrization.
    ///
    /// OCCT: `Perform(CurveOnSurface)`.
    /// `curve_3d` is the 3D curve, `curve_2d` is the pcurve on `surface`.
    /// Samples the curve and finds the maximum 3D deviation.
    pub fn perform(&mut self, curve_3d: &Curve3, curve_2d: &crate::geom::Curve2d, surface: &Surface3) {
        let domain = curve_3d.default_domain();
        // For unbounded curves (e.g. Line), fall back to a finite range so
        // sampling is well-defined.
        let (t_min, t_max) = if !domain[0].is_finite() || !domain[1].is_finite() {
            (-1e6, 1e6)
        } else {
            (domain[0], domain[1])
        };

        let range = t_max - t_min;
        if range < self.tol_range {
            self.error_status = 2;
            return;
        }

        const N_SAMPLES: usize = 257;
        let mut max_d = 0.0;
        let mut max_t = t_min;

        for i in 0..=N_SAMPLES {
            let t = t_min + range * (i as f64) / (N_SAMPLES as f64);
            let p3d = curve_3d.point_at(t);
            let p2d = curve_2d.point_at(t);
            let psurf = surface.point_at(p2d.x, p2d.y);
            let d = (p3d - psurf).length();
            if d > max_d {
                max_d = d;
                max_t = t;
            }
        }

        self.max_distance = max_d;
        self.max_parameter = max_t;
        self.error_status = 0;
    }

    /// Returns true if the max distance has been found.
    ///
    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.error_status == 0
    }

    /// Returns the error status:
    /// 0 = OK, 1 = null curve/surface, 2 = invalid range, 3 = calculation error.
    ///
    /// OCCT: `ErrorStatus()`.
    pub fn error_status(&self) -> i32 {
        self.error_status
    }

    /// Returns the max distance.
    ///
    /// OCCT: `MaxDistance()`.
    pub fn max_distance(&self) -> f64 {
        self.max_distance
    }

    /// Returns the parameter at which max distance occurs.
    ///
    /// OCCT: `MaxParameter()`.
    pub fn max_parameter(&self) -> f64 {
        self.max_parameter
    }
}

impl Default for CheckCurveOnSurface {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// GeomLib::Inertia / GeomLib::AxeOfInertia
// ============================================================================

/// Compute the principal axes of inertia and dispersion values of an array
/// of points.
///
/// OCCT: `GeomLib::Inertia` (GeomLib.cxx L1976-2093).
/// `NCollection_Array1<gp_Pnt>` -> `&[DVec3]`; output parameters -> `&mut`.
pub fn inertia(
    points: &[DVec3],
    bary: &mut DVec3,
    x_dir: &mut DVec3,
    y_dir: &mut DVec3,
    x_gap: &mut f64,
    y_gap: &mut f64,
    z_gap: &mut f64,
) {
    // gp_XYZ GB(0., 0., 0.); int nb = Points.Length();
    let mut gb = DVec3::ZERO;
    let nb = points.len();
    for point in points {
        gb += *point;
    }
    gb /= nb as f64;

    // math_Matrix M(1, 3, 1, 3); M.Init(0.);
    let mut m = MatD::new(3, 3);
    for point in points {
        // Diff.SetLinearForm(-1, Points(i).XYZ(), GB)  ==  GB - Points(i).
        let diff = gb - *point;
        m.set(1, 1, m.get(1, 1) + diff.x * diff.x);
        m.set(2, 2, m.get(2, 2) + diff.y * diff.y);
        m.set(3, 3, m.get(3, 3) + diff.z * diff.z);
        m.set(1, 2, m.get(1, 2) + diff.x * diff.y);
        m.set(1, 3, m.get(1, 3) + diff.x * diff.z);
        m.set(2, 3, m.get(2, 3) + diff.y * diff.z);
    }
    m.set(2, 1, m.get(1, 2));
    m.set(3, 1, m.get(1, 3));
    m.set(3, 2, m.get(2, 3));

    // M /= nb;
    for i in 1..=3usize {
        for j in 1..=3usize {
            m.set(i, j, m.get(i, j) / nb as f64);
        }
    }

    let jacobi = MathJacobi::new(&m);
    // OCCT: if (!J.IsDone()) — debug dump only, no control-flow effect.

    let n1 = jacobi.value(1);
    let n2 = jacobi.value(2);
    let n3 = jacobi.value(3);

    let r1 = n1.min(n2.min(n3));
    let r2: f64;
    let m1: usize;
    let m2: usize;
    let m3: usize;
    if r1 == n1 {
        m1 = 1;
        r2 = n2.min(n3);
        if r2 == n2 {
            m2 = 2;
            m3 = 3;
        } else {
            m2 = 3;
            m3 = 2;
        }
    } else if r1 == n2 {
        m1 = 2;
        r2 = n1.min(n3);
        if r2 == n1 {
            m2 = 1;
            m3 = 3;
        } else {
            m2 = 3;
            m3 = 1;
        }
    } else {
        m1 = 3;
        r2 = n1.min(n2);
        if r2 == n1 {
            m2 = 1;
            m3 = 2;
        } else {
            m2 = 2;
            m3 = 1;
        }
    }

    let v2 = jacobi.vector(m2);
    let v3 = jacobi.vector(m3);

    // gp_Dir::SetCoord normalizes; the eigenvector columns are unit vectors.
    *bary = gb;
    *x_dir = DVec3::new(v3.get(1), v3.get(2), v3.get(3)).normalize_or_zero();
    *y_dir = DVec3::new(v2.get(1), v2.get(2), v2.get(3)).normalize_or_zero();

    *z_gap = jacobi.value(m1).abs().sqrt();
    *y_gap = jacobi.value(m2).abs().sqrt();
    *x_gap = jacobi.value(m3).abs().sqrt();
}

/// Compute the main axis of inertia of an array of points.
///
/// OCCT: `GeomLib::AxeOfInertia` (GeomLib.cxx L2096-2124).  `Axe.XDirection`
/// is the axis of upper inertia; `Axe.Direction` is the normal to the average
/// plane; `is_singular` is true if the points lie on a line.  OCCT default
/// `Tol = 1.0e-7` (GeomLib.hxx L152-156).
pub fn axe_of_inertia(points: &[DVec3], axe: &mut Ax2, is_singular: &mut bool, tol: f64) {
    let mut bary = DVec3::ZERO;
    let mut ox = DVec3::ZERO;
    let mut oy = DVec3::ZERO;
    let mut gx = 0.0f64;
    let mut gy = 0.0f64;
    let mut gz = 0.0f64;

    inertia(points, &mut bary, &mut ox, &mut oy, &mut gx, &mut gy, &mut gz);

    if gy * points.len() as f64 <= tol {
        // OCCT: gp_Ax2 axe(Bary, OX); OY = axe.XDirection().
        let axe2 = Ax2::from_direction(bary, ox);
        oy = axe2.x_direction;
        *is_singular = true;
    } else {
        *is_singular = false;
    }

    // OZ = OX ^ OY; gp_Ax2 TheAxe(Bary, OZ, OX).
    let oz = ox.cross(oy);
    *axe = Ax2::new(bary, oz, ox);
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn test_is_planar_surface_plane() {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let surf = Surface3::Plane(plane);
        let check = IsPlanarSurface::new(&surf, 1e-7);
        assert!(check.is_planar());
    }

    #[test]
    fn test_tool_parameter_curve() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let p = DVec3::new(5.0, 0.0, 0.0);
        let u = Tool::parameter_curve(&line, p, 1e-6);
        assert!(u.is_some());
        assert!((u.unwrap() - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_tool_parameters_surface() {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let surf = Surface3::Plane(plane);
        let p = DVec3::new(3.0, 4.0, 0.0);
        let uv = Tool::parameters_surface(&surf, p, 1e-6);
        assert!(uv.is_some());
        let (u, v) = uv.unwrap();
        assert!((u - 3.0).abs() < 1e-6);
        assert!((v - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_check_curve_on_surface() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let pcurve = Curve2d::Line(Line2d::new(DVec2::ZERO, DVec2::X));
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let mut check = CheckCurveOnSurface::new();
        check.perform(&line, &pcurve, &plane);
        assert!(check.is_done());
        assert!(check.max_distance() < 1e-10);
    }
}
