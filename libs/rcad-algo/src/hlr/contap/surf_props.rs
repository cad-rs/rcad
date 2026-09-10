// OCCT Contap_SurfProps.cxx (TKHLR) — internal tool to compute the normal
// and its derivatives on a surface.
//
// OCCT L25-351 (Normale / DerivAndNorm / NormAndDn).  The `S` argument is
// the `occ::handle<Adaptor3d_Surface>`; rcad carries `&dyn SurfaceAdapter`.
// The canonic branches evaluate through the ElSLib component functions
// (kernel `math/el.rs`), the default branch through the surface D1/D2.

use glam::DVec3;
use rcad_kernel::geom::{Point3, SurfaceEval, Vec3};
use rcad_kernel::precision::COMPUTATIONAL;

use crate::geomalgo::int_patch::GeomAbsSurfaceType;

use super::surface_adaptor::SurfaceAdapter;

/// OCCT gp_Ax3::Direct — X ^ Y . Z > 0.  The rcad analytic frames carry the
/// same invariant through their stored directions.
fn ax3_direct(x: Vec3, y: Vec3, z: Vec3) -> bool {
    x.cross(y).dot(z) > 0.0
}

/// OCCT Contap_SurfProps::Normale (Contap_SurfProps.cxx L25-127) — computes
/// the point <P> and normal vector <Norm> on <S> at parameters U,V.
pub fn normale(s: &dyn SurfaceAdapter, u: f64, v: f64, p: &mut Point3, norm: &mut Vec3) {
    let typ_s = s.get_type();
    match typ_s {
        GeomAbsSurfaceType::Plane => {
            let pl = s.plane();
            *norm = pl.normal;
            *p = rcad_kernel::math::el::elslib_plane_value(u, v, pl.origin, pl.u_dir, pl.v_dir);
            if !ax3_direct(pl.u_dir, pl.v_dir, pl.normal) {
                *norm = -*norm;
            }
        }
        GeomAbsSurfaceType::Sphere => {
            let sp = s.sphere();
            *p = SurfaceEval::point_at(&rcad_kernel::geom::Surface3::Sphere(sp.clone()), u, v);
            *norm = *p - sp.center;
            if ax3_direct(sp.ref_dir, sp.axis.cross(sp.ref_dir), sp.axis) {
                *norm /= sp.radius;
            } else {
                *norm /= -sp.radius;
            }
        }
        GeomAbsSurfaceType::Cylinder => {
            let cy = s.cylinder();
            *p = SurfaceEval::point_at(&rcad_kernel::geom::Surface3::Cylinder(cy.clone()), u, v);
            *norm = u.cos() * cy.ref_dir.normalize_or_zero() + u.sin() * cy.y_axis();
            if !ax3_direct(cy.ref_dir.normalize_or_zero(), cy.y_axis(), cy.axis) {
                *norm = -*norm;
            }
        }
        GeomAbsSurfaceType::Cone => {
            let co = s.cone();
            *p = SurfaceEval::point_at(&rcad_kernel::geom::Surface3::Cone(co.clone()), u, v);
            let angle = co.half_angle_rad;
            let sina = angle.sin();
            let cosa = angle.cos();
            let rad = co.radius;

            let mut vcalc = v;
            if (v * sina + rad).abs() <= 1e-12 {
                // on est a l'apex
                *norm = DVec3::ZERO;
                return;
            }

            if rad + vcalc * sina < 0. {
                *norm = sina * co.axis
                    + cosa * u.cos() * co.ref_dir.normalize_or_zero()
                    + cosa * u.sin() * co.axis.cross(co.ref_dir).normalize_or_zero();
            } else {
                *norm = -sina * co.axis
                    + cosa * u.cos() * co.ref_dir.normalize_or_zero()
                    + cosa * u.sin() * co.axis.cross(co.ref_dir).normalize_or_zero();
            }
            if !ax3_direct(co.ref_dir.normalize_or_zero(), co.axis.cross(co.ref_dir).normalize_or_zero(), co.axis) {
                *norm = -*norm;
            }
        }
        _ => {
            let (pt, d1u, d1v) = s.d1(u, v);
            *p = pt;
            *norm = d1u.cross(d1v);
        }
    }
}

/// OCCT Contap_SurfProps::DerivAndNorm (Contap_SurfProps.cxx L131-232) —
/// computes the point <P>, first derivatives <d1u>, <d1v> and normal
/// vector <Norm> on <S> at parameters U,V.
pub fn deriv_and_norm(
    s: &dyn SurfaceAdapter,
    u: f64,
    v: f64,
    p: &mut Point3,
    d1u: &mut Vec3,
    d1v: &mut Vec3,
    norm: &mut Vec3,
) {
    let typ_s = s.get_type();
    match typ_s {
        GeomAbsSurfaceType::Plane => {
            let pl = s.plane();
            *norm = pl.normal;
            let (pt, du, dv) =
                rcad_kernel::math::el::elslib_plane_d1(u, v, pl.origin, pl.u_dir, pl.v_dir);
            *p = pt;
            *d1u = du;
            *d1v = dv;
            if !ax3_direct(pl.u_dir, pl.v_dir, pl.normal) {
                *norm = -*norm;
            }
        }
        GeomAbsSurfaceType::Sphere => {
            let sp = s.sphere();
            let (pt, du, dv) = SurfaceEval::derivatives(
                &rcad_kernel::geom::Surface3::Sphere(sp.clone()),
                u,
                v,
            );
            *p = pt;
            *d1u = du;
            *d1v = dv;
            *norm = *p - sp.center;
            if ax3_direct(sp.ref_dir, sp.axis.cross(sp.ref_dir), sp.axis) {
                *norm /= sp.radius;
            } else {
                *norm /= -sp.radius;
            }
        }
        GeomAbsSurfaceType::Cylinder => {
            let cy = s.cylinder();
            let (pt, du, dv) = SurfaceEval::derivatives(
                &rcad_kernel::geom::Surface3::Cylinder(cy.clone()),
                u,
                v,
            );
            *p = pt;
            *d1u = du;
            *d1v = dv;
            *norm = u.cos() * cy.ref_dir.normalize_or_zero() + u.sin() * cy.y_axis();
            if !ax3_direct(cy.ref_dir.normalize_or_zero(), cy.y_axis(), cy.axis) {
                *norm = -*norm;
            }
        }
        GeomAbsSurfaceType::Cone => {
            let co = s.cone();
            let (pt, du, dv) =
                SurfaceEval::derivatives(&rcad_kernel::geom::Surface3::Cone(co.clone()), u, v);
            *p = pt;
            *d1u = du;
            *d1v = dv;
            let angle = co.half_angle_rad;
            let sina = angle.sin();
            let cosa = angle.cos();
            let rad = co.radius;

            let mut vcalc = v;
            if (v * sina + rad).abs() <= COMPUTATIONAL {
                // on est a l'apex
                let vfi = s.first_v_parameter();
                if vfi < -rad / sina {
                    // partie valide pour V < Vapex
                    vcalc = v - 1.0;
                } else {
                    vcalc = v + 1.0;
                }
            }

            if rad + vcalc * sina < 0. {
                *norm = sina * co.axis
                    + cosa * u.cos() * co.ref_dir.normalize_or_zero()
                    + cosa * u.sin() * co.axis.cross(co.ref_dir).normalize_or_zero();
            } else {
                *norm = -sina * co.axis
                    + cosa * u.cos() * co.ref_dir.normalize_or_zero()
                    + cosa * u.sin() * co.axis.cross(co.ref_dir).normalize_or_zero();
            }
            if !ax3_direct(co.ref_dir.normalize_or_zero(), co.axis.cross(co.ref_dir).normalize_or_zero(), co.axis) {
                *norm = -*norm;
            }
        }
        _ => {
            let (pt, du, dv) = s.d1(u, v);
            *p = pt;
            *d1u = du;
            *d1v = dv;
            *norm = du.cross(dv);
        }
    }
}

/// OCCT Contap_SurfProps::NormAndDn (Contap_SurfProps.cxx L236-351) —
/// computes the point <P>, normal vector <Norm> and its derivatives
/// <Dnu>, <Dnv> on <S> at parameters U,V.
pub fn norm_and_dn(
    s: &dyn SurfaceAdapter,
    u: f64,
    v: f64,
    p: &mut Point3,
    norm: &mut Vec3,
    dnu: &mut Vec3,
    dnv: &mut Vec3,
) {
    let typ_s = s.get_type();
    match typ_s {
        GeomAbsSurfaceType::Plane => {
            let pl = s.plane();
            *p = rcad_kernel::math::el::elslib_plane_value(u, v, pl.origin, pl.u_dir, pl.v_dir);
            *norm = pl.normal;
            if !ax3_direct(pl.u_dir, pl.v_dir, pl.normal) {
                *norm = -*norm;
            }
            *dnu = DVec3::ZERO;
            *dnv = DVec3::ZERO;
        }
        GeomAbsSurfaceType::Sphere => {
            let sp = s.sphere();
            let (pt, du, dv) = SurfaceEval::derivatives(
                &rcad_kernel::geom::Surface3::Sphere(sp.clone()),
                u,
                v,
            );
            *p = pt;
            *dnu = du;
            *dnv = dv;
            *norm = *p - sp.center;
            let mut rad = sp.radius;
            if !ax3_direct(sp.ref_dir, sp.axis.cross(sp.ref_dir), sp.axis) {
                rad = -rad;
            }
            *norm /= rad;
            *dnu /= rad;
            *dnv /= rad;
        }
        GeomAbsSurfaceType::Cylinder => {
            let cy = s.cylinder();
            *p = SurfaceEval::point_at(&rcad_kernel::geom::Surface3::Cylinder(cy.clone()), u, v);
            *norm = u.cos() * cy.ref_dir.normalize_or_zero() + u.sin() * cy.y_axis();
            *dnu = -u.sin() * cy.ref_dir.normalize_or_zero() + u.cos() * cy.y_axis();
            if !ax3_direct(cy.ref_dir.normalize_or_zero(), cy.y_axis(), cy.axis) {
                *norm = -*norm;
                *dnu = -*dnu;
            }
            *dnv = DVec3::ZERO;
        }
        GeomAbsSurfaceType::Cone => {
            let co = s.cone();
            *p = SurfaceEval::point_at(&rcad_kernel::geom::Surface3::Cone(co.clone()), u, v);
            let angle = co.half_angle_rad;
            let sina = angle.sin();
            let cosa = angle.cos();
            let rad = co.radius;
            let mut vcalc = v;
            if (v * sina + rad).abs() <= COMPUTATIONAL {
                // on est a l'apex
                let vfi = s.first_v_parameter();
                if vfi < -rad / sina {
                    // partie valide pour V < Vapex
                    vcalc = v - 1.0;
                } else {
                    vcalc = v + 1.0;
                }
            }

            if rad + vcalc * sina < 0. {
                *norm = sina * co.axis
                    + cosa * u.cos() * co.ref_dir.normalize_or_zero()
                    + cosa * u.sin() * co.axis.cross(co.ref_dir).normalize_or_zero();
            } else {
                *norm = -sina * co.axis
                    + cosa * u.cos() * co.ref_dir.normalize_or_zero()
                    + cosa * u.sin() * co.axis.cross(co.ref_dir).normalize_or_zero();
            }
            *dnu = -cosa * u.sin() * co.ref_dir.normalize_or_zero() + cosa * u.cos() * co.axis.cross(co.ref_dir).normalize_or_zero();
            if !ax3_direct(co.ref_dir.normalize_or_zero(), co.axis.cross(co.ref_dir).normalize_or_zero(), co.axis) {
                *norm = -*norm;
                *dnu = -*dnu;
            }
            *dnv = DVec3::ZERO;
        }
        _ => {
            let (pt, d1u, d1v, d2u, d2v, d2uv) = s.d2(u, v);
            *p = pt;
            *norm = d1u.cross(d1v);
            *dnu = d2u.cross(d1v) + d1u.cross(d2uv);
            *dnv = d2uv.cross(d1v) + d1u.cross(d2v);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{CylindricalSurface, Plane, Surface3};

    use crate::hlr::contap::surface_adaptor::GeomSurfaceAdapter;

    /// OCCT anchor: on a direct unit cylinder the normal at (0, 0) is the
    /// X direction; at u = pi it is -X; the point follows
    /// P = O + R*(cos u X + sin u Y) + v Z.
    #[test]
    fn surf_props_cylinder_normale() {
        let s = GeomSurfaceAdapter::new(Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
            y_dir: None,
        }));
        let mut p = Point3::ZERO;
        let mut n = Vec3::ZERO;
        normale(&s, 0.0, 5.0, &mut p, &mut n);
        assert!((p.x - 3.0).abs() < 1e-12);
        assert!((p.y - 2.0).abs() < 1e-12);
        assert!((p.z - 8.0).abs() < 1e-12);
        assert!((n.x - 1.0).abs() < 1e-12 && n.y.abs() < 1e-12 && n.z.abs() < 1e-12);

        normale(&s, std::f64::consts::PI, 0.0, &mut p, &mut n);
        assert!((n.x + 1.0).abs() < 1e-12 && n.y.abs() < 1e-12 && n.z.abs() < 1e-12);
    }

    /// OCCT anchor: plane normal is the plane direction (reversed for an
    /// indirect frame); Dnu = Dnv = 0.
    #[test]
    fn surf_props_plane_norm_and_dn() {
        let s = GeomSurfaceAdapter::new(Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(0.0, 0.0, 1.0),
            u_dir: DVec3::new(1.0, 0.0, 0.0),
            v_dir: DVec3::new(0.0, 1.0, 0.0),
        }));
        let mut p = Point3::ZERO;
        let mut n = Vec3::ZERO;
        let mut dnu = Vec3::ZERO;
        let mut dnv = Vec3::ZERO;
        norm_and_dn(&s, 1.0, 2.0, &mut p, &mut n, &mut dnu, &mut dnv);
        assert!((p.x - 1.0).abs() < 1e-12 && (p.y - 2.0).abs() < 1e-12);
        assert!((n.z - 1.0).abs() < 1e-12);
        assert_eq!(dnu, DVec3::ZERO);
        assert_eq!(dnv, DVec3::ZERO);
    }
}
