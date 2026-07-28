//! Projection of a 3D curve onto a surface (GeomProjLib).
//!
//! OCCT TKGeomBase GeomProjLib package.
//!
//! High-level wrapper over ProjLib: accepts any `Curve3` + `Surface3` pair.
//! For analytic curve+surface pairs it delegates to the ProjLib projector;
//! for general pairs it samples the curve, projects each point onto the
//! surface (via ExtPS), and fits a BSpline.

use glam::{DVec2, DVec3};

use crate::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Plane, Surface3, SurfaceEval,
};
use crate::base::proj_lib;

const TOL_DEFAULT: f64 = 1e-7;

/// Project a 3D curve onto a surface, returning the 2D pcurve.
///
/// OCCT: `GeomProjLib::Curve2d(Curve, First, Last, Surface, UFirst, ULast, VFirst, VLast, Tolerance)`.
pub fn curve2d(
    curve: &Curve3,
    first: f64,
    last: f64,
    surface: &Surface3,
    ufirst: f64,
    ulast: f64,
    vfirst: f64,
    vlast: f64,
) -> Option<Curve2d> {
    let tol = TOL_DEFAULT;

    // OCCT: try direct ProjLib projector first
    if let Some(pc) = try_project_direct(curve, surface, first, last) {
        return Some(pc);
    }

    // Fallback: sample the curve, project points, fit BSpline
    let domain = [first, last];
    if !domain[0].is_finite() || !domain[1].is_finite() {
        return None;
    }

    let n = 33.max(2);
    let mut pts_2d = Vec::with_capacity(n);

    for i in 0..n {
        let t = first + (last - first) * (i as f64) / ((n - 1) as f64);
        let p3d = curve.point_at(t);

        // Project onto surface using ExtPS (grid + Newton)
        let uv = project_point_on_surface(p3d, surface, ufirst, ulast, vfirst, vlast, tol)?;
        pts_2d.push(uv);
    }

    if pts_2d.len() < 2 {
        return None;
    }

    Some(Curve2d::BSpline(BSplineCurve2Approx::approximate(&pts_2d)))
}

/// Simplified overload: uses the surface's default domain.
///
/// OCCT: `GeomProjLib::Curve2d(Curve, First, Last, Surface)`.
pub fn curve2d_simple(curve: &Curve3, first: f64, last: f64, surface: &Surface3) -> Option<Curve2d> {
    let dom = surface.default_domain();
    curve2d(curve, first, last, surface, dom[0], dom[1], dom[2], dom[3])
}

/// Simplified overload: uses default domain for both curve and surface.
///
/// OCCT: `GeomProjLib::Curve2d(Curve, Surface)`.
pub fn curve2d_auto(curve: &Curve3, surface: &Surface3) -> Option<Curve2d> {
    let dom_curve = curve.default_domain();
    let dom_surf = surface.default_domain();
    let first = dom_curve[0].max(-1e6);
    let last = dom_curve[1].min(1e6);
    if !first.is_finite() || !last.is_finite() {
        return None;
    }
    curve2d(
        curve,
        first,
        last,
        surface,
        dom_surf[0],
        dom_surf[1],
        dom_surf[2],
        dom_surf[3],
    )
}

/// Try direct ProjLib projection for analytic curve+surface pairs.
fn try_project_direct(curve: &Curve3, surface: &Surface3, _first: f64, _last: f64) -> Option<Curve2d> {
    match (curve, surface) {
        (Curve3::Line(l), Surface3::Plane(p)) => {
            let mut proj = proj_lib::PlaneProjector::with_plane(p);
            proj.project_line(l);
            Some(proj.projector().to_curve2d())
        }
        (Curve3::Circle(c), Surface3::Plane(p)) => {
            let mut proj = proj_lib::PlaneProjector::with_plane(p);
            proj.project_circle(c);
            Some(proj.projector().to_curve2d())
        }
        (Curve3::Ellipse(e), Surface3::Plane(p)) => {
            let mut proj = proj_lib::PlaneProjector::with_plane(p);
            proj.project_ellipse(e);
            Some(proj.projector().to_curve2d())
        }
        (Curve3::Line(l), Surface3::Cylinder(c)) => {
            let mut proj = proj_lib::CylinderProjector::with_cylinder(c);
            proj.project_line(l);
            Some(proj.projector().to_curve2d())
        }
        (Curve3::Circle(c), Surface3::Cylinder(cyl)) => {
            let mut proj = proj_lib::CylinderProjector::with_cylinder(cyl);
            proj.project_circle(c);
            Some(proj.projector().to_curve2d())
        }
        _ => None,
    }
}

/// Project a 3D point onto a surface using grid search + Newton refinement.
fn project_point_on_surface(
    point: DVec3,
    surface: &Surface3,
    ufirst: f64,
    ulast: f64,
    vfirst: f64,
    vlast: f64,
    tol: f64,
) -> Option<DVec2> {
    use crate::base::extrema::ExtPS;

    let ext = ExtPS::with_domain(point, surface, ufirst, ulast, vfirst, vlast, tol, tol);
    if ext.nb_ext() > 0 {
        let p = ext.point(1);
        if (p.point - point).length() < tol * 100.0 {
            return Some(DVec2::new(p.u, p.v));
        }
    }
    None
}

/// Minimal BSplineCurve2 approximation from points (inline helper).
struct BSplineCurve2Approx;

impl BSplineCurve2Approx {
    fn approximate(pts: &[DVec2]) -> crate::geom::BSplineCurve2 {
        let n = pts.len();
        if n < 2 {
            return crate::geom::BSplineCurve2 {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: pts.to_vec(),
                weights: vec![],
            };
        }

        // Chord-length parametrization
        let mut params = vec![0.0_f64; n];
        for i in 1..n {
            let d = (pts[i] - pts[i - 1]).length();
            params[i] = params[i - 1] + d.max(1e-15);
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
            for i in j..j + degree {
                sum += params[i];
            }
            knots[j + degree] = sum / (degree as f64);
        }
        for k in &mut knots[n_knots - degree - 1..] {
            *k = params[n - 1];
        }

        crate::geom::BSplineCurve2 {
            degree,
            knots,
            control_points: pts.to_vec(),
            weights: vec![],
        }
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
    fn test_curve2d_line_plane() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let result = curve2d_auto(&line, &plane);
        assert!(result.is_some());
        let pc = result.unwrap();
        let p0 = pc.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-10);
        let p5 = pc.point_at(5.0);
        assert!((p5 - DVec2::new(5.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn test_curve2d_circle_plane() {
        let circle = Curve3::Circle(Circle3::new(DVec3::new(0.0, 0.0, 5.0), DVec3::Z, 3.0));
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let result = curve2d_simple(&circle, 0.0, 6.283, &plane);
        assert!(result.is_some());
    }
}
