//! Bounding-box computation for geometric entities (BndLib).
//!
//! OCCT TKGeomBase BndLib package: BndLib_Add3dCurve, BndLib_AddSurface.
//!
//! In OCCT, `Bnd_Box` (the data structure) lives in TKMath while `BndLib`
//! (algorithms that add `Geom_*` curves/surfaces to a `Bnd_Box`) lives in
//! TKGeomBase.  The functions below correspond to `BndLib_Add3dCurve` and
//! `BndLib_AddSurface`.

use glam::DVec3;
use crate::geom::{Curve3, Surface3, SurfaceEval};

/// Conservative bounding box for an analytic curve.
/// OCCT: BndLib_Add3dCurve::Add.
pub fn curve_bounding_box(curve: &Curve3) -> Option<[DVec3; 2]> {
    match curve {
        Curve3::Circle(c) => {
            let n = c.normal.normalize();
            let extent = DVec3::new(
                c.radius * (1.0 - n.x * n.x).sqrt(),
                c.radius * (1.0 - n.y * n.y).sqrt(),
                c.radius * (1.0 - n.z * n.z).sqrt(),
            );
            Some([c.center - extent, c.center + extent])
        }
        Curve3::Ellipse(e) => {
            let n = e.normal.normalize();
            let max_r = e.major_radius.max(e.minor_radius);
            let extent = DVec3::new(
                max_r * (1.0 - n.x * n.x).sqrt(),
                max_r * (1.0 - n.y * n.y).sqrt(),
                max_r * (1.0 - n.z * n.z).sqrt(),
            );
            Some([e.center - extent, e.center + extent])
        }
        Curve3::BSpline(b) => {
            let (mut mn, mut mx) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for &p in &b.control_points { mn = mn.min(p); mx = mx.max(p); }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        Curve3::Bezier(b) => {
            let (mut mn, mut mx) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for &p in &b.control_points { mn = mn.min(p); mx = mx.max(p); }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        _ => None,
    }
}

/// Conservative bounding box for an analytic surface.
/// OCCT: BndLib_AddSurface::Add.
pub fn surface_bounding_box(
    surface: &Surface3,
    vertices: &[crate::topo::topology::Vertex],
) -> Option<[DVec3; 2]> {
    match surface {
        Surface3::Cylinder(cyl) => {
            let r = cyl.radius;
            let axis = cyl.axis.normalize_or_zero();
            if axis.length_squared() < 0.5 { return None; }
            let (mut min_axial, mut max_axial) = (f64::INFINITY, f64::NEG_INFINITY);
            for v in vertices {
                let proj = (v.point - cyl.origin).dot(axis);
                min_axial = min_axial.min(proj); max_axial = max_axial.max(proj);
            }
            if !min_axial.is_finite() { return None; }
            let p_lo = cyl.origin + axis * min_axial;
            let p_hi = cyl.origin + axis * max_axial;
            let rv = DVec3::splat(r);
            Some([p_lo.min(p_hi) - rv, p_lo.max(p_hi) + rv])
        }
        Surface3::Sphere(sph) => Some([sph.center - DVec3::splat(sph.radius), sph.center + DVec3::splat(sph.radius)]),
        Surface3::Torus(tor) => {
            let r = DVec3::splat(tor.major_radius + tor.minor_radius);
            Some([tor.center - r, tor.center + r])
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis_dir();
            let apex = cone.apex_point();
            let (mut min_axial, mut max_axial) = (f64::INFINITY, f64::NEG_INFINITY);
            for v in vertices {
                let proj = (v.point - apex).dot(axis);
                min_axial = min_axial.min(proj); max_axial = max_axial.max(proj);
            }
            if !min_axial.is_finite() { return None; }
            let max_r = cone.radius_at_axial(min_axial).max(cone.radius_at_axial(max_axial));
            let rv = DVec3::splat(max_r.max(cone.radius));
            let p_lo = apex + axis * min_axial;
            let p_hi = apex + axis * max_axial;
            Some([p_lo.min(p_hi) - rv, p_lo.max(p_hi) + rv])
        }
        Surface3::Ellipsoid(e) => {
            let max_r = e.radius_x.max(e.radius_y).max(e.radius_z);
            Some([e.center - DVec3::splat(max_r), e.center + DVec3::splat(max_r)])
        }
        Surface3::BSpline(b) => {
            let (mut mn, mut mx) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for row in &b.control_points { for p in row { mn = mn.min(*p); mx = mx.max(*p); } }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        Surface3::Bezier(b) => {
            let (mut mn, mut mx) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for row in &b.control_points { for &p in row { mn = mn.min(p); mx = mx.max(p); } }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        _ => None,
    }
}
