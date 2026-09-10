// OCCT Geom_Surface UIso/VIso — the per-type iso-curve constructions of the
// analytic quadrics plus the RectangularTrimmedSurface delegation, consumed
// by the EnlargeGeometry statics of brep_offset_tool_c.rs (split out for the
// 2000-line module rule; the architecture-difference numbering continues in
// brep_offset_tool.rs).

use rcad_kernel::geom::{Circle3, Curve3, Line3, Surface3, TrimmedCurve3};

/// OCCT Geom_Surface::UIso(U) — the virtual dispatch over the rcad Surface3
/// carriers (architecture bridge: the closed Surface3 enum replaces the
/// OCCT dynamic dispatch).  The analytic quadrics follow the Geom_*::UIso
/// wrappers over ElSLib: PlaneUIso (ElSLib.cxx L1705-1715, Geom_Plane.cxx
/// L261-265), CylinderUIso (ElSLib.cxx L1716-1726 = CylinderD1(U, 0) ->
/// Lin(P, D1V), Geom_CylindricalSurface.cxx L294-298), ConeUIso
/// (ElSLib.cxx L1727-1737 = ConeD1(U, 0) -> Lin(P, D1V),
/// Geom_ConicalSurface.cxx L337-341), SphereUIso (ElSLib.cxx L1738-1750,
/// Geom_SphericalSurface.cxx), TorusUIso (ElSLib.cxx L1751-1766,
/// Geom_ToroidalSurface.cxx).  The RectangularTrimmedSurface delegates to
/// the basis with the v-trim wrap (Geom_RectangularTrimmedSurface.cxx
/// L444-459).
pub(crate) fn surface_uiso(s: &Surface3, u: f64) -> Curve3 {
    match s {
        Surface3::Plane(a_pln) => {
            // OCCT ElSLib::PlaneUIso: L(Pos.Location(), Pos.YDirection())
            // translated by XDirection()*U.
            Curve3::Line(Line3 {
                origin: a_pln.origin + a_pln.u_dir * u,
                direction: a_pln.v_dir,
            })
        }
        Surface3::Cylinder(a_cyl) => {
            // OCCT ElSLib::CylinderUIso: CylinderD1(U, 0) —
            // P = Loc + R cos U X + R sin U Y, D1V = ZDirection.
            let a_y = a_cyl.y_axis();
            let p = a_cyl.origin + a_cyl.ref_dir * (a_cyl.radius * u.cos()) + a_y * (a_cyl.radius * u.sin());
            Curve3::Line(Line3::new(p, a_cyl.axis))
        }
        Surface3::Cone(a_con) => {
            // OCCT ElSLib::ConeUIso: ConeD1(U, 0) —
            // P = Loc + R cos U X + R sin U Y,
            // D1V = ZDirection*cos(SA) + sin(SA)*(cos U X + sin U Y).
            let a_y = a_con.axis.cross(a_con.ref_dir);
            let a_cos_u = u.cos();
            let a_sin_u = u.sin();
            let p = a_con.apex
                + a_con.ref_dir * (a_con.radius * a_cos_u)
                + a_y * (a_con.radius * a_sin_u);
            let a_dv = a_con.axis * a_con.half_angle_rad.cos()
                + (a_con.ref_dir * a_cos_u + a_y * a_sin_u) * a_con.half_angle_rad.sin();
            Curve3::Line(Line3::new(p, a_dv))
        }
        Surface3::Sphere(a_sph) => {
            // OCCT ElSLib::SphereUIso: cx = cos U X + sin U Y;
            // Circ(Ax2(Location, cx.Crossed(dz), cx), Radius).
            let a_cx = a_sph.ref_dir * u.cos() + a_sph.ref_dir_perp() * u.sin();
            let a_n = a_cx.cross(a_sph.axis);
            Curve3::Circle(Circle3 {
                center: a_sph.center,
                normal: a_n,
                x_dir: a_cx,
                y_dir: a_n.cross(a_cx),
                radius: a_sph.radius,
            })
        }
        Surface3::Torus(a_tor) => {
            // OCCT ElSLib::TorusUIso: cx = cos U X + sin U Y;
            // axes(Ax2(Location, cx.Crossed(dz), cx)) translated by
            // cx*MajorRadius; Circ(axes, MinorRadius).
            let a_y = a_tor.axis.cross(a_tor.ref_dir);
            let a_cx = a_tor.ref_dir * u.cos() + a_y * u.sin();
            let a_n = a_cx.cross(a_tor.axis);
            Curve3::Circle(Circle3 {
                center: a_tor.center + a_cx * a_tor.major_radius,
                normal: a_n,
                x_dir: a_cx,
                y_dir: a_n.cross(a_cx),
                radius: a_tor.minor_radius,
            })
        }
        Surface3::Trimmed(a_ts) => {
            // OCCT Geom_RectangularTrimmedSurface::UIso: the basis UIso,
            // wrapped in the v-trim when the surface is v-trimmed
            // (cxx L444-459; trim = [u1, u2, v1, v2]).
            let a_c = surface_uiso(&a_ts.basis, u);
            let a_v1 = a_ts.trim[2];
            let a_v2 = a_ts.trim[3];
            Curve3::Trimmed(TrimmedCurve3::new(a_c, a_v1, a_v2))
        }
        a_other => panic!(
            "GAP: Geom_Surface::UIso for {:?} (iso-curve construction not translated)",
            a_other
        ),
    }
}

/// OCCT Geom_Surface::VIso(V) — the sibling of surface_uiso: PlaneVIso
/// (ElSLib.cxx L1781-1791), CylinderVIso (L1793-1804), ConeVIso
/// (L1806-1827), SphereVIso (L1829-1855), TorusVIso (L1836-1854), and the
/// RectangularTrimmedSurface u-trim wrap (cxx L463-478).
pub(crate) fn surface_viso(s: &Surface3, v: f64) -> Curve3 {
    match s {
        Surface3::Plane(a_pln) => {
            // OCCT ElSLib::PlaneVIso: L(Pos.Location(), Pos.XDirection())
            // translated by YDirection()*V.
            Curve3::Line(Line3 {
                origin: a_pln.origin + a_pln.v_dir * v,
                direction: a_pln.u_dir,
            })
        }
        Surface3::Cylinder(a_cyl) => {
            // OCCT ElSLib::CylinderVIso: Ax2(Pos) translated by
            // Direction()*V; Circ(axes, Radius).
            Curve3::Circle(Circle3 {
                center: a_cyl.origin + a_cyl.axis * v,
                normal: a_cyl.axis,
                x_dir: a_cyl.ref_dir,
                y_dir: a_cyl.y_axis(),
                radius: a_cyl.radius,
            })
        }
        Surface3::Cone(a_con) => {
            // OCCT ElSLib::ConeVIso: Ax3(Pos) translated by
            // Direction()*(V*cos(SA)); R = Radius + V*sin(SA); when R < 0
            // the X/Y directions reverse and R is negated.
            let a_y = a_con.axis.cross(a_con.ref_dir);
            let a_r = a_con.radius + v * a_con.half_angle_rad.sin();
            let (a_x_dir, a_y_dir, a_r) = if a_r < 0. {
                (-a_con.ref_dir, -a_y, -a_r)
            } else {
                (a_con.ref_dir, a_y, a_r)
            };
            Curve3::Circle(Circle3 {
                center: a_con.apex + a_con.axis * (v * a_con.half_angle_rad.cos()),
                normal: a_con.axis,
                x_dir: a_x_dir,
                y_dir: a_y_dir,
                radius: a_r,
            })
        }
        Surface3::Sphere(a_sph) => {
            // OCCT ElSLib::SphereVIso: Ax2(Pos) translated by
            // Direction()*(Radius*sin(V)); radius = Radius*cos(V); when
            // radius < 0 the axis direction flips and the radius is
            // negated (the #23170 analytic-continuation form).
            let a_r = a_sph.radius * v.cos();
            let (a_n, a_r) = if a_r < 0. {
                (-a_sph.axis, -a_r)
            } else {
                (a_sph.axis, a_r)
            };
            Curve3::Circle(Circle3 {
                center: a_sph.center + a_sph.axis * (a_sph.radius * v.sin()),
                normal: a_n,
                x_dir: a_sph.ref_dir,
                y_dir: a_sph.ref_dir_perp(),
                radius: a_r,
            })
        }
        Surface3::Torus(a_tor) => {
            // OCCT ElSLib::TorusVIso: Ax3(Pos) translated by
            // Direction()*(MinorRadius*sin(V)); R = MajorRadius +
            // MinorRadius*cos(V); when R < 0 the X/Y directions reverse
            // and R is negated.
            let a_y = a_tor.axis.cross(a_tor.ref_dir);
            let a_r = a_tor.major_radius + a_tor.minor_radius * v.cos();
            let (a_x_dir, a_y_dir, a_r) = if a_r < 0. {
                (-a_tor.ref_dir, -a_y, -a_r)
            } else {
                (a_tor.ref_dir, a_y, a_r)
            };
            Curve3::Circle(Circle3 {
                center: a_tor.center + a_tor.axis * (a_tor.minor_radius * v.sin()),
                normal: a_tor.axis,
                x_dir: a_x_dir,
                y_dir: a_y_dir,
                radius: a_r,
            })
        }
        Surface3::Trimmed(a_ts) => {
            // OCCT Geom_RectangularTrimmedSurface::VIso: the basis VIso,
            // wrapped in the u-trim when the surface is u-trimmed
            // (cxx L463-478).
            let a_c = surface_viso(&a_ts.basis, v);
            let a_u1 = a_ts.trim[0];
            let a_u2 = a_ts.trim[1];
            Curve3::Trimmed(TrimmedCurve3::new(a_c, a_u1, a_u2))
        }
        a_other => panic!(
            "GAP: Geom_Surface::VIso for {:?} (iso-curve construction not translated)",
            a_other
        ),
    }
}

