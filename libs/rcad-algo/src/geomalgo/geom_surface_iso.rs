//! OCCT Geom_Surface::UIso / VIso (ModelingData/TKG3d/Geom) — the single rcad
//! translation of the virtual isoparametric-curve dispatch, shared by the two
//! isoline consumers (`geomalgo::geom_lib_iso_line` = GeomLib::buildC3dOnIsoLine
//! and `geomalgo::approx_curve_on_surface` = Approx_CurveOnSurface).
//!
//! Architecture difference: OCCT carries this as one virtual function per
//! concrete surface class (`Geom_Plane::UIso`, `Geom_CylindricalSurface::UIso`,
//! ... plus the `Geom_RectangularTrimmedSurface::UIso` complementary trim).
//! The rcad `Surface3` enum lives in `rcad-kernel` (geom), which owns no
//! geometry-constructing surface API, so the single translation of that
//! dispatch lives here and both consumers delegate to it (previously two
//! non-identical partial copies existed, one covering Trimmed, the other
//! Revolution — this is one OCCT function, so the copies are merged).
//!
//! Arm set (the union of both former copies): Trimmed + Plane + Cylinder +
//! Cone + Sphere + Torus + Revolution + BSpline. The remaining `Surface3`
//! variants keep the OCCT-faithful failure path (the panics below):
//! - Bezier / Offset / LinearExtrusion have OCCT overrides
//!   (Geom_BezierSurface.cxx L1769/L1821, Geom_OffsetSurface.cxx L601/L657,
//!   Geom_SurfaceOfLinearExtrusion.cxx L275/L285) that are not translated yet;
//! - Ellipsoid / Helicoid / Pipe / Ruled / Coons / TriBezier are rcad-only
//!   `Surface3` variants with no Geom override in TKG3d/Geom at all.

use glam::DVec3;

use rcad_kernel::base::proj_lib::elslib_iso as el;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, TrimmedCurve3, TrimmedSurface, BSplineCurve3};
use rcad_kernel::math::bspl_lib;

/// OCCT Geom_RectangularTrimmedSurface(S, U1, U2, V1, V2, USense, VSense)
/// (Geom_RectangularTrimmedSurface.cxx L67-112) — the nested trimmed basis is
/// killed and the resulting surface carries `isutrimmed = isvtrimmed = true`.
pub(crate) fn surface_rectangular_trimmed(
    surf: &Surface3,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Surface3 {
    // OCCT: kill trimmed basis surfaces.
    let basis = match surf {
        Surface3::Trimmed(t) => (*t.basis).clone(),
        other => other.clone(),
    };
    Surface3::Trimmed(TrimmedSurface::new(basis, u1, u2, v1, v2))
}

/// OCCT Geom_Surface::UIso(U) — the concrete-class overrides:
/// - Geom_RectangularTrimmedSurface.cxx L444-459 (complementary V trim),
/// - Geom_Plane.cxx L261-265, Geom_CylindricalSurface.cxx L294-298,
///   Geom_ConicalSurface.cxx L337-341, Geom_SphericalSurface.cxx L292-297,
///   Geom_ToroidalSurface.cxx L305-310 (all via ElSLib),
/// - Geom_SurfaceOfRevolution.cxx L372-379 (rotated copy of the basis curve),
/// - Geom_BSplineSurface_1.cxx L598-635 (BSplSLib::Iso on the U direction).
pub(crate) fn surface_u_iso(surf: &Surface3, param: f64) -> Curve3 {
    match surf {
        Surface3::Trimmed(t) => {
            // OCCT L447-458: C = basisSurf->UIso(U); if (isvtrimmed) return
            // new Geom_TrimmedCurve(C, vtrim1, vtrim2, true).
            let c = surface_u_iso(&t.basis, param);
            Curve3::Trimmed(TrimmedCurve3::new(c, t.trim[2], t.trim[3]))
        }
        // OCCT Geom_Plane.cxx L263: new Geom_Line(ElSLib::PlaneUIso(pos, U)).
        Surface3::Plane(p) => Curve3::Line(el::elslib_plane_u_iso(
            &el::Ax3View::from_axes(p.origin, p.normal, p.u_dir),
            param,
        )),
        // OCCT Geom_CylindricalSurface.cxx L296:
        // new Geom_Line(ElSLib::CylinderUIso(pos, radius, U)).
        Surface3::Cylinder(c) => Curve3::Line(el::elslib_cylinder_u_iso(
            &el::Ax3View::from_axes(c.origin, c.axis, c.ref_dir),
            c.radius,
            param,
        )),
        // OCCT Geom_ConicalSurface.cxx L339:
        // new Geom_Line(ElSLib::ConeUIso(pos, radius, semiAngle, U)).
        Surface3::Cone(c) => Curve3::Line(el::elslib_cone_u_iso(
            &el::Ax3View::from_axes(c.apex, c.axis, c.ref_dir),
            c.radius,
            c.half_angle_rad,
            param,
        )),
        // OCCT Geom_SphericalSurface.cxx L294-296: the ElSLib circle is wrapped
        // in new Geom_TrimmedCurve(GC, -M_PI/2., M_PI/2) — the UIso of a sphere
        // is only the meridian semicircle.
        Surface3::Sphere(s) => Curve3::Trimmed(TrimmedCurve3::new(
            Curve3::Circle(el::elslib_sphere_u_iso(
                &el::Ax3View::from_axes(s.center, s.axis, s.ref_dir),
                s.radius,
                param,
            )),
            -std::f64::consts::FRAC_PI_2,
            std::f64::consts::FRAC_PI_2,
        )),
        // OCCT Geom_ToroidalSurface.cxx L307-308:
        // new Geom_Circle(ElSLib::TorusUIso(pos, majorRadius, minorRadius, U)).
        Surface3::Torus(t) => Curve3::Circle(el::elslib_torus_u_iso(
            &el::Ax3View::from_axes(t.center, t.axis, t.ref_dir),
            t.major_radius,
            t.minor_radius,
            param,
        )),
        // OCCT Geom_SurfaceOfRevolution.cxx L374-378:
        //   C = basisCurve->Copy(); C->Rotate(Ax1(loc, direction), U); return C.
        Surface3::Revolution(r) => {
            rotate_curve_about_axis(&r.profile, r.axis_origin, r.axis_dir, param)
        }
        // OCCT Geom_BSplineSurface_1.cxx L598-635: BSplSLib::Iso on the U
        // direction; the result curve's knots/degree/multiplicities are the V
        // direction's, and the curve carries myVPeriodic.
        Surface3::BSpline(b) => {
            let weights = if b.is_rational_u() || b.is_rational_v() {
                Some(b.weights.as_slice())
            } else {
                None
            };
            let (cpoles, cweights) = bspl_lib::bspl_slib_iso(
                param,
                true,
                b.degree_u,
                &b.knots_u,
                &b.control_points,
                weights,
                b.is_periodic_u,
            );
            Curve3::BSpline(BSplineCurve3 {
                degree: b.degree_v,
                knots: b.knots_v.clone(),
                control_points: cpoles,
                weights: cweights,
                is_periodic: b.is_periodic_v,
            })
        }
        _ => panic!("GAP: Geom_Surface::UIso not translated for this surface type"),
    }
}

/// OCCT Geom_Surface::VIso(V) — the V-isoparametric counterpart of
/// [`surface_u_iso`]:
/// - Geom_RectangularTrimmedSurface.cxx L463-478 (complementary U trim),
/// - Geom_Plane.cxx L269-273, Geom_CylindricalSurface.cxx L302-306,
///   Geom_ConicalSurface.cxx L345-349, Geom_SphericalSurface.cxx L301-305,
///   Geom_ToroidalSurface.cxx L314-319,
/// - Geom_SurfaceOfRevolution.cxx L383-410 (the parallel circle of the basis
///   point through the axis),
/// - Geom_BSplineSurface_1.cxx L775-812 (BSplSLib::Iso on the V direction).
pub(crate) fn surface_v_iso(surf: &Surface3, param: f64) -> Curve3 {
    match surf {
        Surface3::Trimmed(t) => {
            // OCCT L466-477: C = basisSurf->VIso(V); if (isutrimmed) return
            // new Geom_TrimmedCurve(C, utrim1, utrim2, true).
            let c = surface_v_iso(&t.basis, param);
            Curve3::Trimmed(TrimmedCurve3::new(c, t.trim[0], t.trim[1]))
        }
        // OCCT Geom_Plane.cxx L271: new Geom_Line(ElSLib::PlaneVIso(pos, V)).
        Surface3::Plane(p) => Curve3::Line(el::elslib_plane_v_iso(
            &el::Ax3View::from_axes(p.origin, p.normal, p.u_dir),
            param,
        )),
        // OCCT Geom_CylindricalSurface.cxx L304:
        // new Geom_Circle(ElSLib::CylinderVIso(pos, radius, V)).
        Surface3::Cylinder(c) => Curve3::Circle(el::elslib_cylinder_v_iso(
            &el::Ax3View::from_axes(c.origin, c.axis, c.ref_dir),
            c.radius,
            param,
        )),
        // OCCT Geom_ConicalSurface.cxx L347:
        // new Geom_Circle(ElSLib::ConeVIso(pos, radius, semiAngle, V)).
        Surface3::Cone(c) => Curve3::Circle(el::elslib_cone_v_iso(
            &el::Ax3View::from_axes(c.apex, c.axis, c.ref_dir),
            c.radius,
            c.half_angle_rad,
            param,
        )),
        // OCCT Geom_SphericalSurface.cxx L303:
        // new Geom_Circle(ElSLib::SphereVIso(pos, radius, V)) — untrimmed.
        Surface3::Sphere(s) => Curve3::Circle(el::elslib_sphere_v_iso(
            &el::Ax3View::from_axes(s.center, s.axis, s.ref_dir),
            s.radius,
            param,
        )),
        // OCCT Geom_ToroidalSurface.cxx L316-317:
        // new Geom_Circle(ElSLib::TorusVIso(pos, majorRadius, minorRadius, V)).
        Surface3::Torus(t) => Curve3::Circle(el::elslib_torus_v_iso(
            &el::Ax3View::from_axes(t.center, t.axis, t.ref_dir),
            t.major_radius,
            t.minor_radius,
            param,
        )),
        // OCCT Geom_SurfaceOfRevolution.cxx L383-410: the circle of the basis
        // point at V about the axis.  Rad = distance from the axis; the circle
        // frame is gp_Ax2(C, direction, D) where C is the projection of the
        // basis point onto the axis and D the unit vector from C to it.
        // OCCT tests the inner `P.Modulus() > gp::Resolution()` as well, but
        // P = Pc - C so |P| == Rad, hence the single distance test carries both
        // branches (the degenerate case falls back to the frame's default X).
        Surface3::Revolution(r) => {
            let pc = r.profile.point_at(param);
            let d = pc - r.axis_origin;
            let rad = (d - r.axis_dir * d.dot(r.axis_dir)).length();
            let c = r.axis_origin + r.axis_dir * d.dot(r.axis_dir);
            let radial = pc - c;
            let (normal, x_dir) = if rad > rcad_kernel::precision::CONFUSION {
                (r.axis_dir, radial.normalize_or_zero())
            } else {
                // OCCT: Rep = gp_Ax2(C, direction) — the zero-radius case uses
                // the frame's default X direction.
                let x = rcad_kernel::geom::any_perpendicular(r.axis_dir);
                (r.axis_dir, x)
            };
            Curve3::Circle(rcad_kernel::geom::Circle3 {
                center: c,
                normal,
                x_dir,
                y_dir: normal.cross(x_dir).normalize_or_zero(),
                radius: rad,
            })
        }
        // OCCT Geom_BSplineSurface_1.cxx L775-812: BSplSLib::Iso on the V
        // direction; the result curve carries the U direction's knots, degree
        // and myUPeriodic.
        Surface3::BSpline(b) => {
            let weights = if b.is_rational_u() || b.is_rational_v() {
                Some(b.weights.as_slice())
            } else {
                None
            };
            let (cpoles, cweights) = bspl_lib::bspl_slib_iso(
                param,
                false,
                b.degree_v,
                &b.knots_v,
                &b.control_points,
                weights,
                b.is_periodic_v,
            );
            Curve3::BSpline(BSplineCurve3 {
                degree: b.degree_u,
                knots: b.knots_u.clone(),
                control_points: cpoles,
                weights: cweights,
                is_periodic: b.is_periodic_u,
            })
        }
        _ => panic!("GAP: Geom_Surface::VIso not translated for this surface type"),
    }
}

/// OCCT Geom_Geometry::Rotate (Geom_Geometry.cxx L47-55) applied to a curve:
/// the per-type `Rotate(Ax1, Ang)` overrides.  The rotation is a rigid motion
/// so a TrimmedCurve rotates its basis and keeps its parameter range
/// (Geom_TrimmedCurve::Transform).
fn rotate_curve_about_axis(
    curve: &Curve3,
    axis_loc: DVec3,
    axis_dir: DVec3,
    angle: f64,
) -> Curve3 {
    let axis_dir = axis_dir.normalize_or_zero();
    let trsf = glam::DAffine3::from_translation(axis_loc)
        * glam::DAffine3::from_axis_angle(axis_dir, angle)
        * glam::DAffine3::from_translation(-axis_loc);
    match curve {
        Curve3::Trimmed(t) => Curve3::Trimmed(TrimmedCurve3::new(
            rotate_curve_about_axis(&t.curve, axis_loc, axis_dir, angle),
            t.first,
            t.last,
        )),
        other => rcad_kernel::geom::transform_curve(other, &trsf),
    }
}
