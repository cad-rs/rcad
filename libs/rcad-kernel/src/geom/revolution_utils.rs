//! OCCT Geom_RevolutionUtils (TKG3d/Geom/Geom_RevolutionUtils.pxx) — the
//! surface-of-revolution evaluation helpers (pxx L43-265, the
//! `CalculateD0/D1/D2/DN` bodies), the `D0` / `D1` / `D2` / `DN` template
//! wrappers (pxx L275-442), and the rcad re-hosts of
//! `Geom_SurfaceOfRevolution::EvalD0/EvalD1/EvalD2/EvalDN`
//! (Geom_SurfaceOfRevolution.cxx L236-248 / L252-270 / L274-296 / L331-360)
//! plus the `GeomAdaptor_Surface::EvalDN` revolution arm
//! (GeomAdaptor_Surface.cxx L1768-1780).
//!
//! Revolution surface: `P(U, V) = Rotation(Axis, U) * BasisCurve(V)` where U is
//! the rotation angle and V is the parameter along the basis curve.
//!
//! Architecture differences:
//!   - `gp_Trsf::SetRotation(gp_Ax1, angle)` + `P.Transform(T)` map to the
//!     kernel [`Trsf`] ([`Trsf::set_rotation`] / [`Trsf::apply`] for points and
//!     [`Trsf::transform_vec`] for the `gp_Vec` outputs) — the same gp_Trsf
//!     translation already used by the rest of rcad;
//!   - `occ::handle(Geom_Curve) basisCurve` maps to the rcad [`Curve3`] value
//!     and the `gp_Ax1(loc, direction)` member pair maps to the
//!     [`RevolutionSurface`] `axis_origin` / `axis_dir` fields;
//!   - `theBasis.EvalDN(V, Nv)` (the `Geom_Curve` virtual) is the rcad
//!     [`crate::geom::curve_dn::curve_dn`] union, and `theBasis.EvalD0(V)` is
//!     [`CurveEval::point_at`].

use glam::DVec3;

use super::curve_dn::curve_dn;
use super::offset_surface_utils::{ResD1, ResD2};
use crate::core::precision::SQUARE_CONFUSION;
use crate::geom::{Curve3, CurveEval, RevolutionSurface};
use crate::math::gp::{Ax1, Trsf};

// =========================================================================
// The Geom_RevolutionUtils Calculate* leaves (pxx L43-265)
// =========================================================================

/// OCCT `Geom_RevolutionUtils::CalculateD0` (pxx L43-52).
pub fn calculate_d0(the_curve_pt: DVec3, the_u: f64, the_axis: &Ax1) -> DVec3 {
    let mut a_rotation = Trsf::identity();
    a_rotation.set_rotation(the_axis, the_u);
    a_rotation.apply(the_curve_pt)
}

/// OCCT `Geom_RevolutionUtils::CalculateD1` (pxx L62-88).
pub fn calculate_d1(
    the_curve_pt: DVec3,
    the_curve_d1: DVec3,
    the_u: f64,
    the_axis: &Ax1,
) -> ResD1 {
    // Vector from center of rotation to the point on rotated curve.
    let a_cq = the_curve_pt - the_axis.location;
    let mut the_d1u = the_axis.direction.cross(a_cq);
    // If the point is placed on the axis of revolution then derivatives on U
    // are undefined.  Manually set them to zero.
    if the_d1u.length_squared() < SQUARE_CONFUSION {
        the_d1u = DVec3::ZERO;
    }

    let mut a_rotation = Trsf::identity();
    a_rotation.set_rotation(the_axis, the_u);
    let the_p = a_rotation.apply(the_curve_pt);
    let the_d1u = a_rotation.transform_vec(the_d1u);
    let the_d1v = a_rotation.transform_vec(the_curve_d1);
    ResD1 {
        point: the_p,
        d1u: the_d1u,
        d1v: the_d1v,
    }
}

/// OCCT `Geom_RevolutionUtils::CalculateD2` (pxx L103-140).
pub fn calculate_d2(
    the_curve_pt: DVec3,
    the_curve_d1: DVec3,
    the_curve_d2: DVec3,
    the_u: f64,
    the_axis: &Ax1,
) -> ResD2 {
    // Vector from center of rotation to the point on rotated curve.
    let a_cq = the_curve_pt - the_axis.location;
    let a_dir = the_axis.direction;
    let mut the_d1u = a_dir.cross(a_cq);
    // If the point is placed on the axis of revolution then derivatives on U
    // are undefined.  Manually set them to zero.
    if the_d1u.length_squared() < SQUARE_CONFUSION {
        the_d1u = DVec3::ZERO;
    }
    let the_d2u = a_dir.dot(a_cq) * a_dir - a_cq;
    let the_d2uv = a_dir.cross(the_curve_d1);

    let mut a_rotation = Trsf::identity();
    a_rotation.set_rotation(the_axis, the_u);
    let the_p = a_rotation.apply(the_curve_pt);
    let the_d1u = a_rotation.transform_vec(the_d1u);
    let the_d1v = a_rotation.transform_vec(the_curve_d1);
    let the_d2u = a_rotation.transform_vec(the_d2u);
    let the_d2v = a_rotation.transform_vec(the_curve_d2);
    let the_d2uv = a_rotation.transform_vec(the_d2uv);
    ResD2 {
        point: the_p,
        d1u: the_d1u,
        d1v: the_d1v,
        d2u: the_d2u,
        d2v: the_d2v,
        d2uv: the_d2uv,
    }
}

/// OCCT `Geom_RevolutionUtils::CalculateDN` (pxx L222-265).
///
/// `the_der_v` is part of the interface contract (pxx L228-229) — the caller
/// provides different data based on its value.
pub fn calculate_dn(
    the_curve_pt_or_dn: DVec3,
    the_u: f64,
    the_axis: &Ax1,
    the_der_u: i32,
    the_der_v: i32,
) -> DVec3 {
    let _ = the_der_v;

    let mut a_rotation = Trsf::identity();
    a_rotation.set_rotation(the_axis, the_u);

    let a_result = if the_der_u == 0 {
        // Pure V derivative: just rotate the curve derivative.
        the_curve_pt_or_dn
    } else {
        // For theDerV == 0: theCurvePtOrDN is (P - AxisLocation) as a vector.
        // For theDerV > 0: theCurvePtOrDN is the curve derivative.
        let a_dir = the_axis.direction;
        if the_der_u % 4 == 1 {
            a_dir.cross(the_curve_pt_or_dn)
        } else if the_der_u % 4 == 2 {
            a_dir.dot(the_curve_pt_or_dn) * a_dir - the_curve_pt_or_dn
        } else if the_der_u % 4 == 3 {
            a_dir.cross(the_curve_pt_or_dn) * -1.0
        } else {
            the_curve_pt_or_dn - a_dir.dot(the_curve_pt_or_dn) * a_dir
        }
    };

    a_rotation.transform_vec(a_result)
}

// =========================================================================
// The Geom_RevolutionUtils D0/D1/D2/DN wrappers (pxx L275-442)
// =========================================================================

/// OCCT `Geom_RevolutionUtils::D0` (pxx L275-285) — the basis-curve
/// `EvalD0(theV)` + `CalculateD0`.
pub fn d0(the_u: f64, the_v: f64, the_basis: &Curve3, the_axis: &Ax1) -> DVec3 {
    calculate_d0(the_basis.point_at(the_v), the_u, the_axis)
}

/// OCCT `Geom_RevolutionUtils::D1` (pxx L297-309) — the basis-curve
/// `EvalD1(theV)` + `CalculateD1`.
pub fn d1(the_u: f64, the_v: f64, the_basis: &Curve3, the_axis: &Ax1) -> ResD1 {
    calculate_d1(
        the_basis.point_at(the_v),
        the_basis.derivative_at(the_v),
        the_u,
        the_axis,
    )
}

/// OCCT `Geom_RevolutionUtils::D2` (pxx L325-350) — the basis-curve
/// `EvalD2(theV)` + `CalculateD2`.
pub fn d2(the_u: f64, the_v: f64, the_basis: &Curve3, the_axis: &Ax1) -> ResD2 {
    calculate_d2(
        the_basis.point_at(the_v),
        the_basis.derivative_at(the_v),
        the_basis.derivative2_at(the_v),
        the_u,
        the_axis,
    )
}

/// OCCT `Geom_RevolutionUtils::DN` (pxx L416-442) — the basis-curve
/// `EvalD0`/`EvalDN` selection + `CalculateDN`.
pub fn dn(
    the_u: f64,
    the_v: f64,
    the_basis: &Curve3,
    the_axis: &Ax1,
    the_der_u: i32,
    the_der_v: i32,
) -> DVec3 {
    let a_curve_pt_or_dn = if the_der_u == 0 {
        curve_dn(the_basis, the_v, the_der_v)
    } else if the_der_v == 0 {
        let a_p = the_basis.point_at(the_v);
        a_p - the_axis.location
    } else {
        curve_dn(the_basis, the_v, the_der_v)
    };

    calculate_dn(a_curve_pt_or_dn, the_u, the_axis, the_der_u, the_der_v)
}

// =========================================================================
// Geom_SurfaceOfRevolution::EvalD0/EvalD1/EvalD2/EvalDN re-hosts
// =========================================================================

/// The `gp_Ax1(loc, direction)` of a rcad [`RevolutionSurface`] payload.
fn revolution_axis(rev: &RevolutionSurface) -> Ax1 {
    Ax1::new(rev.axis_origin, rev.axis_dir)
}

/// OCCT `Geom_SurfaceOfRevolution::EvalD0` (cxx L236-248) — the
/// `GeomEval_RepUtils::TryEvalSurfaceD0` short circuit is the rcad payload's
/// own concern, not reproduced here.
pub fn revolution_eval_d0(rev: &RevolutionSurface, u: f64, v: f64) -> DVec3 {
    d0(u, v, &rev.profile, &revolution_axis(rev))
}

/// OCCT `Geom_SurfaceOfRevolution::EvalD1` (cxx L252-270).
pub fn revolution_eval_d1(rev: &RevolutionSurface, u: f64, v: f64) -> ResD1 {
    d1(u, v, &rev.profile, &revolution_axis(rev))
}

/// OCCT `Geom_SurfaceOfRevolution::EvalD2` (cxx L274-296).
pub fn revolution_eval_d2(rev: &RevolutionSurface, u: f64, v: f64) -> ResD2 {
    d2(u, v, &rev.profile, &revolution_axis(rev))
}

/// OCCT `Geom_SurfaceOfRevolution::EvalDN` (cxx L331-360).
pub fn revolution_eval_dn(rev: &RevolutionSurface, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    // OCCT: if (Nu + Nv < 1 || Nu < 0 || Nv < 0) throw Geom_UndefinedDerivative.
    assert!(
        nu + nv >= 1 && nu >= 0 && nv >= 0,
        "Geom_UndefinedDerivative: Geom_SurfaceOfRevolution::EvalDN"
    );
    dn(u, v, &rev.profile, &revolution_axis(rev), nu, nv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle3, Line3, SurfaceEval};

    /// A line `P(v) = (d, 0, v)` parallel to Z at distance `d` from the Z axis,
    /// revolved about Z: the surface is the cylinder of radius `d` with
    /// `S(u, v) = (d cos u, d sin u, v)`.
    ///
    /// Exact derivatives: `Su = d(-sin u, cos u, 0)`, `Sv = (0, 0, 1)`,
    /// `Suu = -d(cos u, sin u, 0)`, `Suv = Suv = 0`, `Svv = 0`,
    /// `Suuu = d(sin u, -cos u, 0)`.
    fn revolved_line(d: f64) -> RevolutionSurface {
        RevolutionSurface {
            profile: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(d, 0.0, 0.0),
                direction: DVec3::Z,
            })),
            axis_origin: DVec3::ZERO,
            axis_dir: DVec3::Z,
        }
    }

    #[test]
    fn revolution_derivatives_of_a_revolved_line() {
        let d = 1.25f64;
        let rev = revolved_line(d);
        let (u, v) = (0.8f64, 0.6f64);
        let (su, cu) = u.sin_cos();
        let want_p = DVec3::new(d * cu, d * su, v);
        let want_d1u = DVec3::new(-d * su, d * cu, 0.0);
        let want_d1v = DVec3::Z;
        let want_d2u = DVec3::new(-d * cu, -d * su, 0.0);
        let want_d3u = DVec3::new(d * su, -d * cu, 0.0);

        assert!(
            (revolution_eval_d0(&rev, u, v) - want_p).length() < 1e-14,
            "d0={:?}",
            revolution_eval_d0(&rev, u, v)
        );
        let d1 = revolution_eval_d1(&rev, u, v);
        assert!((d1.point - want_p).length() < 1e-14);
        assert!((d1.d1u - want_d1u).length() < 1e-14, "d1u={:?}", d1.d1u);
        assert!((d1.d1v - want_d1v).length() < 1e-14, "d1v={:?}", d1.d1v);
        let d2 = revolution_eval_d2(&rev, u, v);
        assert!((d2.point - want_p).length() < 1e-14);
        assert!((d2.d1u - want_d1u).length() < 1e-14);
        assert!((d2.d1v - want_d1v).length() < 1e-14);
        assert!((d2.d2u - want_d2u).length() < 1e-14, "d2u={:?}", d2.d2u);
        assert!(d2.d2uv.length() < 1e-14, "d2uv={:?}", d2.d2uv);
        assert!(d2.d2v.length() < 1e-14, "d2v={:?}", d2.d2v);
        assert!((revolution_eval_dn(&rev, u, v, 1, 0) - want_d1u).length() < 1e-14);
        assert!((revolution_eval_dn(&rev, u, v, 0, 1) - want_d1v).length() < 1e-14);
        assert!((revolution_eval_dn(&rev, u, v, 2, 0) - want_d2u).length() < 1e-14);
        assert!((revolution_eval_dn(&rev, u, v, 3, 0) - want_d3u).length() < 1e-13);
        assert!(revolution_eval_dn(&rev, u, v, 0, 2).length() < 1e-14);
    }

    /// A point on the axis (`d = 0`) is the singular case of
    /// `CalculateD1` (pxx L75-80): the U derivative is forced to zero rather
    /// than left at the (zero-magnitude) cross product.
    #[test]
    fn revolution_singular_point_on_the_axis() {
        let rev = revolved_line(0.0);
        let d1 = revolution_eval_d1(&rev, 0.7, 0.4);
        assert!(d1.d1u.length() < 1e-15, "d1u={:?}", d1.d1u);
        let d2 = revolution_eval_d2(&rev, 0.7, 0.4);
        assert!(d2.d1u.length() < 1e-15);
    }

    /// A circle profile in the XZ plane revolved about Z: the sphere of radius
    /// `R`.  With `C(v) = (R sin v, 0, R cos v)` (a circle of radius R about
    /// the Y axis in the XZ plane) the revolution reduces to
    /// `S(u, v) = R (sin v cos u, sin v sin u, cos v)` — a second, independent
    /// oracle for the D2/DN chain (the mixed term is non-zero here).
    #[test]
    fn revolution_of_a_profile_circle_reduces_to_the_sphere() {
        let r = 2.0f64;
        let rev = RevolutionSurface {
            profile: Box::new(Curve3::Circle(Circle3 {
                center: DVec3::ZERO,
                normal: DVec3::Y,
                radius: r,
                x_dir: DVec3::Z,
                y_dir: DVec3::X,
            })),
            axis_origin: DVec3::ZERO,
            axis_dir: DVec3::Z,
        };
        let (u, v) = (0.35f64, 0.9f64);
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let want_p = DVec3::new(r * sv * cu, r * sv * su, r * cv);
        let want_d1u = DVec3::new(-r * sv * su, r * sv * cu, 0.0);
        let want_d1v = DVec3::new(r * cv * cu, r * cv * su, -r * sv);
        // The U second derivative is the rotation's own second derivative:
        // `-(P - P_axis)`, i.e. only the component perpendicular to the axis
        // flips sign (`Geom_RevolutionUtils::CalculateD2` builds
        // `aDir.Dot(aCQ)*aDir - aCQ`).
        let want_d2u = -DVec3::new(r * sv * cu, r * sv * su, 0.0);
        let want_d2uv = DVec3::new(-r * cv * su, r * cv * cu, 0.0);
        let want_d2v = DVec3::new(-r * sv * cu, -r * sv * su, -r * cv);

        // The payload's own point formula is the sphere.
        assert!(
            (SurfaceEval::point_at(&rev, u, v) - want_p).length() < 1e-14,
            "point_at={:?} want={want_p:?}",
            SurfaceEval::point_at(&rev, u, v)
        );
        let d1 = revolution_eval_d1(&rev, u, v);
        assert!((d1.point - want_p).length() < 1e-13, "p={:?}", d1.point);
        assert!((d1.d1u - want_d1u).length() < 1e-13, "d1u={:?}", d1.d1u);
        assert!((d1.d1v - want_d1v).length() < 1e-13, "d1v={:?}", d1.d1v);
        let d2 = revolution_eval_d2(&rev, u, v);
        assert!((d2.d2u - want_d2u).length() < 1e-13, "d2u={:?}", d2.d2u);
        assert!((d2.d2uv - want_d2uv).length() < 1e-13, "d2uv={:?}", d2.d2uv);
        assert!((d2.d2v - want_d2v).length() < 1e-13, "d2v={:?}", d2.d2v);
        assert!((revolution_eval_dn(&rev, u, v, 1, 1) - want_d2uv).length() < 1e-12);
        assert!((revolution_eval_dn(&rev, u, v, 0, 2) - want_d2v).length() < 1e-12);
    }
}
