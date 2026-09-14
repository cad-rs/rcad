//! Regression guard for the Geom_OffsetSurface::UIso / VIso arms
//! (Geom_OffsetSurface.cxx L601-688): the equivalent-surface fast path
//! (`directRepSurface` -> `Geom_OffsetSurface::Surface()`), the
//! `GeomAbs_SurfaceOfExtrusion` branch and the general AdvApprox arm
//! driven by `Geom_OffsetSurfaceUtils::EvaluateD0/D1`.

use super::*;
use rcad_kernel::geom::{OffsetSurface, Plane, SphericalSurface};
use rcad_kernel::geom::SurfaceEval as _;

/// Pointwise comparison of an iso curve against the offset surface value
/// at the fixed parameter.
fn check_offset_iso(iso: &Curve3, off: &Surface3, fixed: f64, is_u: bool) {
    for jj in 0..=8 {
        let t = jj as f64 / 8.0;
        let p_iso = iso.point_at(t);
        let p_ref = if is_u {
            off.point_at(fixed, t)
        } else {
            off.point_at(t, fixed)
        };
        assert!(
            (p_iso - p_ref).length() < 1e-7,
            "t={t} iso={p_iso:?} ref={p_ref:?}"
        );
    }
}

/// A plane basis: `Surface()` yields the plane translated by
/// `offsetValue * normal`, so both iso arms delegate to
/// `Geom_Plane::UIso` / `VIso` (Geom_OffsetSurface.cxx L900-907 + L652).
#[test]
fn offset_of_plane_delegates_to_the_translated_plane_iso() {
    let base = Surface3::Plane(Plane {
        origin: DVec3::new(1.0, 2.0, 3.0),
        normal: DVec3::Z,
        u_dir: DVec3::X,
        v_dir: DVec3::Y,
    });
    let off = Surface3::Offset(OffsetSurface {
        basis: Box::new(base),
        offset_distance: 3.0,
    });

    let uiso = surface_uiso(&off, 0.4);
    let Curve3::Line(l) = &uiso else {
        panic!("Geom_Plane::UIso is a Geom_Line")
    };
    assert!((l.origin - DVec3::new(1.4, 2.0, 6.0)).length() < 1e-12);
    assert!((l.direction - DVec3::Y).length() < 1e-12);
    check_offset_iso(&uiso, &off, 0.4, true);

    let viso = surface_viso(&off, 0.7);
    let Curve3::Line(lv) = &viso else {
        panic!("Geom_Plane::VIso is a Geom_Line")
    };
    assert!((lv.origin - DVec3::new(1.0, 2.7, 6.0)).length() < 1e-12);
    assert!((lv.direction - DVec3::X).length() < 1e-12);
    check_offset_iso(&viso, &off, 0.7, false);
}

/// A zero offset returns the basis surface itself
/// (Geom_OffsetSurface.cxx L870-873).
#[test]
fn zero_offset_delegates_to_the_basis_surface_iso() {
    let base = Surface3::Sphere(SphericalSurface {
        center: DVec3::new(1.0, 2.0, 3.0),
        axis: DVec3::Z,
        radius: 2.5,
        ref_dir: DVec3::X,
    });
    let off = Surface3::Offset(OffsetSurface {
        basis: Box::new(base.clone()),
        offset_distance: 0.0,
    });
    let uiso = surface_uiso(&off, 0.4);
    let uiso_base = surface_uiso(&base, 0.4);
    for jj in 0..=8 {
        let t = jj as f64 / 8.0;
        assert!((uiso.point_at(t) - uiso_base.point_at(t)).length() < 1e-12);
    }
}

/// A cylinder basis: `Surface()` yields the cylinder of radius
/// `R + aSign * offsetValue` (Geom_OffsetSurface.cxx L908-925), so the
/// VIso is the parallel circle of the offset radius.
#[test]
fn offset_of_cylinder_viso_is_the_offset_parallel_circle() {
    use rcad_kernel::geom::CylindricalSurface;
    let base = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(1.0, 2.0, 3.0),
        axis: DVec3::Z,
        radius: 2.5,
        ref_dir: DVec3::X,
        y_dir: None,
    });
    let off = Surface3::Offset(OffsetSurface {
        basis: Box::new(base),
        offset_distance: 2.0,
    });

    let viso = surface_viso(&off, 0.7);
    let Curve3::Circle(c) = &viso else {
        panic!("Geom_CylindricalSurface::VIso is a Geom_Circle")
    };
    assert!((c.radius - 4.5).abs() < 1e-12);
    check_offset_iso(&viso, &off, 0.7, false);

    let uiso = surface_uiso(&off, 0.4);
    assert!(matches!(uiso, Curve3::Line(_)));
    check_offset_iso(&uiso, &off, 0.4, true);
}

/// An extrusion basis takes the `GeomAbs_SurfaceOfExtrusion` arm
/// (Geom_OffsetSurface.cxx L606-624): the ruling is translated by
/// `offsetValue * normalized(D1U ^ D1V)`.
#[test]
fn offset_of_extrusion_uiso_translates_the_ruling() {
    use rcad_kernel::geom::LinearExtrusionSurface;
    let profile = Curve3::Circle(Circle3 {
        center: DVec3::new(1.0, 2.0, 3.0),
        normal: DVec3::Z,
        x_dir: DVec3::X,
        y_dir: DVec3::Y,
        radius: 2.0,
    });
    let off = Surface3::Offset(OffsetSurface {
        basis: Box::new(Surface3::LinearExtrusion(LinearExtrusionSurface {
            profile: Box::new(profile),
            direction: DVec3::Z,
        })),
        offset_distance: 1.5,
    });

    let uiso = surface_uiso(&off, 0.0);
    let Curve3::Line(l) = &uiso else {
        panic!("the extrusion UIso arm returns the basis ruling line")
    };
    // C(0) = (3, 2, 3) and the radial direction at u = 0 is +X.
    assert!((l.origin - DVec3::new(4.5, 2.0, 3.0)).length() < 1e-12);
    assert!((l.direction - DVec3::Z).length() < 1e-12);
    check_offset_iso(&uiso, &off, 0.0, true);
}

/// A BSpline basis has no canonical equivalent surface, so the iso arms
/// run the general `Geom_OffsetSurface_UIsoEvaluator` / AdvApprox body
/// (Geom_OffsetSurface.cxx L625-654) whose Value/D1 calls route through
/// `Geom_OffsetSurfaceUtils::EvaluateD0/EvaluateD1`.
#[test]
fn offset_of_bspline_basis_runs_the_approximation_arm() {
    let patch = BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
            vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
        ],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        is_periodic_u: false,
        is_periodic_v: false,
    };
    let off = Surface3::Offset(OffsetSurface {
        basis: Box::new(Surface3::BSpline(patch)),
        offset_distance: 0.5,
    });

    let uiso = surface_uiso(&off, 0.3);
    assert!(
        matches!(uiso, Curve3::BSpline(_)),
        "the approximation arm returns a Geom_BSplineCurve"
    );
    check_offset_iso(&uiso, &off, 0.3, true);

    let viso = surface_viso(&off, 0.7);
    assert!(matches!(viso, Curve3::BSpline(_)));
    check_offset_iso(&viso, &off, 0.7, false);
}
