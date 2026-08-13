//! Planar face construction helpers (legacy API used by generated DRAW tests).
//!
//! OCCT `mkplane` equivalent: build a single rectangular planar face.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};

/// Build a rectangular planar face with the given UV extents.
///
/// Mirrors the removed legacy `make_planar_rect_brep` (originally built on the
/// old builder API) using the current `topods::BRep` pool. The result is a
/// single plane face (OCCT `mkplane` semantics).
pub fn make_planar_rect_brep(
    origin: DVec3,
    u_axis: DVec3,
    v_axis: DVec3,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> Result<BRep, crate::BuildError> {
    let normal = u_axis.cross(v_axis).normalize();
    let surface = Surface3::Plane(Plane::new(origin, normal));

    let c0 = origin + u_axis * umin + v_axis * vmin;
    let c1 = origin + u_axis * umax + v_axis * vmin;
    let c2 = origin + u_axis * umax + v_axis * vmax;
    let c3 = origin + u_axis * umin + v_axis * vmax;

    let mut brep = BRep::new();
    let v0 = brep.add_tvertex(c0);
    let v1 = brep.add_tvertex(c1);
    let v2 = brep.add_tvertex(c2);
    let v3 = brep.add_tvertex(c3);

    let line = |p0: DVec3, p1: DVec3| -> Curve3 {
        let d = p1 - p0;
        Curve3::Line(Line3 {
            origin: p0,
            direction: d.normalize(),
        })
    };
    let e0 = brep.add_tedge(Some(line(c0, c1)), v0.clone(), v1.clone(), [0.0, (c1 - c0).length()]);
    let e1 = brep.add_tedge(Some(line(c1, c2)), v1.clone(), v2.clone(), [0.0, (c2 - c1).length()]);
    let e2 = brep.add_tedge(Some(line(c2, c3)), v2.clone(), v3.clone(), [0.0, (c3 - c2).length()]);
    let e3 = brep.add_tedge(Some(line(c3, c0)), v3.clone(), v0.clone(), [0.0, (c0 - c3).length()]);

    let wire = brep.add_twire(vec![e0, e1, e2, e3]);
    brep.add_tface(Some(surface), wire, vec![], None, None, vec![], false);
    Ok(brep)
}
