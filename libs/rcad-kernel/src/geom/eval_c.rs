//! The evaluation leaves of the OCCT `GeomEval_*` surface family that the
//! rcad surface payloads re-host with a different parameterisation or that
//! have no other rcad home:
//!
//!   - `GeomEval_EllipsoidSurface::EvalD1/EvalD2/EvalDN`
//!     (TKG3d/GeomEval/GeomEval_EllipsoidSurface.cxx L197-221 / L225-259 /
//!     L316-346) for the rcad [`EllipsoidalSurface`] payload.  The V
//!     parameterisation differs between the two (see
//!     [`ellipsoid_eval_d1`]) — the OCCT body is translated literally in the
//!     OCCT parameterisation and evaluated through the documented
//!     `V_occt = Pi/2 - V_rcad` substitution;
//!   - `GeomEval_CircularHelicoidSurface::EvalD1/EvalD2/EvalDN`
//!     (GeomEval_CircularHelicoidSurface.cxx L160-183 / L187-217 / L265-315)
//!     for the rcad [`HelicoidSurface`] payload, whose parameterisation is
//!     identical (`P(U, V) = O + V*cos(U)*XD + V*sin(U)*YD + P*U/(2*Pi)*ZD`).
//!
//! Architecture differences: the OCCT bodies read the `gp_Ax3 pos` frame
//! members (`pos.XDirection()` / `YDirection()` / `Direction()` /
//! `Location()`) while the rcad payloads carry `axis` + `ref_dir` and derive
//! the frame with `geom::orthonormal_frame` (X = `ref_dir` projected, Y = Z ^ X
//! — the gp_Ax3 right-handed convention, gp_Ax3.hxx L57-58); the
//! `GeomEval_RepUtils::TryEvalSurfaceD*` short circuits do not exist for these
//! payloads (they are the surface's own evaluation).
//!
//! `EvalD3` (cxx L263-312 / L221-261) is not translated: rcad has no consumer.

use glam::DVec3;

use super::offset_surface_utils::{ResD1, ResD2};
use super::orthonormal_frame;
use crate::geom::{EllipsoidalSurface, HelicoidSurface};
use std::f64::consts::PI;

// =========================================================================
// GeomEval_EllipsoidSurface (GeomEval_EllipsoidSurface.cxx)
// =========================================================================

/// The [`EllipsoidalSurface`] payload frame as the OCCT `gp_Ax3 pos` members
/// `(XD, YD, ZD)`; the location is `center`.
fn ellipsoid_frame(el: &EllipsoidalSurface) -> (DVec3, DVec3, DVec3) {
    let (axis, x_axis, y_axis) = orthonormal_frame(el.axis, el.ref_dir);
    (x_axis, y_axis, axis)
}

/// OCCT `GeomEval_EllipsoidSurface::EvalD1` (cxx L197-221) — the literal body
/// in the OCCT parameterisation (V = the latitude, `Bounds` = `[-Pi/2, Pi/2]`).
#[allow(clippy::too_many_arguments)]
fn ellipsoid_eval_d1_occt(
    o: DVec3,
    xd: DVec3,
    yd: DVec3,
    zd: DVec3,
    my_a: f64,
    my_b: f64,
    my_c: f64,
    u: f64,
    v: f64,
) -> ResD1 {
    let (a_sin_u, a_cos_u) = u.sin_cos();
    let (a_sin_v, a_cos_v) = v.sin_cos();

    ResD1 {
        point: o + my_a * a_cos_v * a_cos_u * xd
            + my_b * a_cos_v * a_sin_u * yd
            + my_c * a_sin_v * zd,
        // dP/du = -A*cv*su*XD + B*cv*cu*YD
        d1u: my_a * a_cos_v * (-a_sin_u) * xd + my_b * a_cos_v * a_cos_u * yd,
        // dP/dv = -A*sv*cu*XD - B*sv*su*YD + C*cv*ZD
        d1v: my_a * (-a_sin_v) * a_cos_u * xd
            + my_b * (-a_sin_v) * a_sin_u * yd
            + my_c * a_cos_v * zd,
    }
}

/// OCCT `GeomEval_EllipsoidSurface::EvalD2` (cxx L225-259) — the literal body
/// in the OCCT parameterisation.
#[allow(clippy::too_many_arguments)]
fn ellipsoid_eval_d2_occt(
    o: DVec3,
    xd: DVec3,
    yd: DVec3,
    zd: DVec3,
    my_a: f64,
    my_b: f64,
    my_c: f64,
    u: f64,
    v: f64,
) -> ResD2 {
    let (a_sin_u, a_cos_u) = u.sin_cos();
    let (a_sin_v, a_cos_v) = v.sin_cos();

    ResD2 {
        point: o + my_a * a_cos_v * a_cos_u * xd
            + my_b * a_cos_v * a_sin_u * yd
            + my_c * a_sin_v * zd,
        // dP/du
        d1u: my_a * a_cos_v * (-a_sin_u) * xd + my_b * a_cos_v * a_cos_u * yd,
        // dP/dv
        d1v: my_a * (-a_sin_v) * a_cos_u * xd
            + my_b * (-a_sin_v) * a_sin_u * yd
            + my_c * a_cos_v * zd,
        // d2P/du2 = -A*cv*cu*XD - B*cv*su*YD
        d2u: my_a * a_cos_v * (-a_cos_u) * xd + my_b * a_cos_v * (-a_sin_u) * yd,
        // d2P/dv2 = -A*cv*cu*XD - B*cv*su*YD - C*sv*ZD
        d2v: my_a * (-a_cos_v) * a_cos_u * xd
            + my_b * (-a_cos_v) * a_sin_u * yd
            + my_c * (-a_sin_v) * zd,
        // d2P/dudv = A*sv*su*XD - B*sv*cu*YD
        d2uv: my_a * a_sin_v * a_sin_u * xd + my_b * (-a_sin_v) * a_cos_u * yd,
    }
}

/// OCCT `GeomEval_EllipsoidSurface::EvalDN` (cxx L316-346) — the literal body
/// in the OCCT parameterisation.  The `XDir` derivative coefficients use the
/// phase-shift identities `d^n/du^n[cos(u)] = cos(u + n*Pi/2)` and
/// `d^n/du^n[sin(u)] = sin(u + n*Pi/2)` (cxx L330-334).
#[allow(clippy::too_many_arguments)]
fn ellipsoid_eval_dn_occt(
    xd: DVec3,
    yd: DVec3,
    zd: DVec3,
    my_a: f64,
    my_b: f64,
    my_c: f64,
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
) -> DVec3 {
    let a_phase_u = nu as f64 * PI / 2.0;
    let a_phase_v = nv as f64 * PI / 2.0;

    // XD coefficient: A * d^Nv[cos(v)] * d^Nu[cos(u)]
    let a_coeff_x = my_a * (v + a_phase_v).cos() * (u + a_phase_u).cos();
    // YD coefficient: B * d^Nv[cos(v)] * d^Nu[sin(u)]
    let a_coeff_y = my_b * (v + a_phase_v).cos() * (u + a_phase_u).sin();
    // ZD coefficient: C * d^Nv[sin(v)] * (Nu == 0 ? 1 : 0)
    let a_coeff_z = my_c * (v + a_phase_v).sin() * if nu == 0 { 1.0 } else { 0.0 };

    a_coeff_x * xd + a_coeff_y * yd + a_coeff_z * zd
}

/// OCCT `GeomEval_EllipsoidSurface::EvalD1` (cxx L197-221) over the rcad
/// [`EllipsoidalSurface`] payload.
///
/// Parameterisation difference (verified against both sources, not papered
/// over): `GeomEval_EllipsoidSurface` uses `V` as the **latitude** —
/// `P = O + A*cos(V)*cos(U)*XD + B*cos(V)*sin(U)*YD + C*sin(V)*ZD` with
/// `Bounds` `V in [-Pi/2, Pi/2]` (cxx L192 / `.hxx` L52-57), i.e. `V = 0` on
/// the equator and `V = Pi/2` at `+C*ZD`.  The rcad [`EllipsoidalSurface`] uses
/// `V` as the **colatitude** (0 at the `+axis` pole, `geom/mod.rs` L642-646)
/// and its consumers rely on that (`shhealing::shape_analysis` L492-493:
/// "Ellipsoid has two poles at v=0 and v=PI").
///
/// The two point formulas are the same function of `V` under
/// `V_occt = Pi/2 - V_rcad`: substituting `sin(V_occt) = cos(V_rcad)` and
/// `cos(V_occt) = sin(V_rcad)` into the OCCT `EvalD0` reproduces
/// `EllipsoidalSurface::point_at` exactly.  The translated OCCT body is
/// therefore evaluated at `V_occt` and every `V`-partial of order `n` picks up
/// the chain factor `(dV_occt/dV_rcad)^n = (-1)^n`.
pub fn ellipsoid_eval_d1(el: &EllipsoidalSurface, u: f64, v: f64) -> ResD1 {
    let (xd, yd, zd) = ellipsoid_frame(el);
    let v_occt = PI / 2.0 - v;
    let mut a_result = ellipsoid_eval_d1_occt(
        el.center,
        xd,
        yd,
        zd,
        el.radius_x,
        el.radius_y,
        el.radius_z,
        u,
        v_occt,
    );
    a_result.d1v = -a_result.d1v;
    a_result
}

/// OCCT `GeomEval_EllipsoidSurface::EvalD2` (cxx L225-259) over the rcad
/// [`EllipsoidalSurface`] payload — the `V_occt = Pi/2 - V_rcad` substitution
/// of [`ellipsoid_eval_d1`], with the chain factors `(-1)^n` on the
/// `V`-partials.
pub fn ellipsoid_eval_d2(el: &EllipsoidalSurface, u: f64, v: f64) -> ResD2 {
    let (xd, yd, zd) = ellipsoid_frame(el);
    let v_occt = PI / 2.0 - v;
    let mut a_result = ellipsoid_eval_d2_occt(
        el.center,
        xd,
        yd,
        zd,
        el.radius_x,
        el.radius_y,
        el.radius_z,
        u,
        v_occt,
    );
    a_result.d1v = -a_result.d1v;
    a_result.d2uv = -a_result.d2uv;
    // d2v carries (-1)^2 = +1.
    a_result
}

/// OCCT `GeomEval_EllipsoidSurface::EvalDN` (cxx L316-346) over the rcad
/// [`EllipsoidalSurface`] payload — the `V_occt = Pi/2 - V_rcad` substitution
/// of [`ellipsoid_eval_d1`], with the chain factor `(-1)^Nv`.
pub fn ellipsoid_eval_dn(el: &EllipsoidalSurface, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    // OCCT: if (Nu + Nv < 1 || Nu < 0 || Nv < 0) throw Geom_UndefinedDerivative
    // ("GeomEval_EllipsoidSurface::EvalDN: invalid derivative order").
    assert!(
        nu + nv >= 1 && nu >= 0 && nv >= 0,
        "Geom_UndefinedDerivative: GeomEval_EllipsoidSurface::EvalDN: invalid derivative order"
    );
    let (xd, yd, zd) = ellipsoid_frame(el);
    let v_occt = PI / 2.0 - v;
    let a_result = ellipsoid_eval_dn_occt(
        xd,
        yd,
        zd,
        el.radius_x,
        el.radius_y,
        el.radius_z,
        u,
        v_occt,
        nu,
        nv,
    );
    if nv % 2 == 1 {
        -a_result
    } else {
        a_result
    }
}

// =========================================================================
// GeomEval_CircularHelicoidSurface (GeomEval_CircularHelicoidSurface.cxx)
// =========================================================================

/// The [`HelicoidSurface`] payload frame as the OCCT `gp_Ax3 pos` members
/// `(XD, YD, ZD)`; the location is `origin`.
fn helicoid_frame(h: &HelicoidSurface) -> (DVec3, DVec3, DVec3) {
    let (axis, x_axis, y_axis) = orthonormal_frame(h.axis, h.ref_dir);
    (x_axis, y_axis, axis)
}

/// OCCT `GeomEval_CircularHelicoidSurface::EvalD1` (cxx L160-183) over the
/// rcad [`HelicoidSurface`] payload.  The parameterisations are identical:
/// `S(u, v) = O + v*cos(u)*XD + v*sin(u)*YD + (Pitch/(2*Pi))*u*ZD`
/// (cxx L143-156 / `geom/mod.rs` L1147-1157).
pub fn helicoid_eval_d1(h: &HelicoidSurface, u: f64, v: f64) -> ResD1 {
    let (xd, yd, zd) = helicoid_frame(h);
    let (a_sin_u, a_cos_u) = u.sin_cos();

    let a_z_rate = h.pitch / (2.0 * PI);

    ResD1 {
        point: h.origin + v * a_cos_u * xd + v * a_sin_u * yd + (h.pitch * u / (2.0 * PI)) * zd,
        // dS/du = -v*sin(u)*XDir + v*cos(u)*YDir + P/(2*Pi)*ZDir
        d1u: v * (-a_sin_u) * xd + v * a_cos_u * yd + a_z_rate * zd,
        // dS/dv = cos(u)*XDir + sin(u)*YDir
        d1v: a_cos_u * xd + a_sin_u * yd,
    }
}

/// OCCT `GeomEval_CircularHelicoidSurface::EvalD2` (cxx L187-217) over the
/// rcad [`HelicoidSurface`] payload.
pub fn helicoid_eval_d2(h: &HelicoidSurface, u: f64, v: f64) -> ResD2 {
    let (xd, yd, zd) = helicoid_frame(h);
    let (a_sin_u, a_cos_u) = u.sin_cos();

    let a_z_rate = h.pitch / (2.0 * PI);

    ResD2 {
        point: h.origin + v * a_cos_u * xd + v * a_sin_u * yd + (h.pitch * u / (2.0 * PI)) * zd,
        d1u: v * (-a_sin_u) * xd + v * a_cos_u * yd + a_z_rate * zd,
        d1v: a_cos_u * xd + a_sin_u * yd,
        // d2S/dudu = -v*cos(u)*XDir - v*sin(u)*YDir
        d2u: v * (-a_cos_u) * xd + v * (-a_sin_u) * yd,
        // d2S/dvdv = 0 (ruled surface)
        d2v: DVec3::ZERO,
        // d2S/dudv = -sin(u)*XDir + cos(u)*YDir
        d2uv: (-a_sin_u) * xd + a_cos_u * yd,
    }
}

/// OCCT `GeomEval_CircularHelicoidSurface::EvalDN` (cxx L265-315) over the
/// rcad [`HelicoidSurface`] payload.  The v-dependence is linear, so
/// `Nv >= 2` yields zero, and the u-dependence of the Z term is linear, so
/// only `Nu == 1` carries a Z contribution (cxx L276-304).
pub fn helicoid_eval_dn(h: &HelicoidSurface, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    // OCCT: if (Nu + Nv < 1 || Nu < 0 || Nv < 0) throw Geom_UndefinedDerivative.
    assert!(
        nu + nv >= 1 && nu >= 0 && nv >= 0,
        "Geom_UndefinedDerivative: GeomEval_CircularHelicoidSurface::EvalDN"
    );
    // OCCT (cxx L276-282): the v-dependence is linear, so Nv >= 2 vanishes
    // before the phase-shift coefficients are built.
    if nv >= 2 {
        return DVec3::ZERO;
    }

    let (xd, yd, zd) = helicoid_frame(h);

    let a_phase_u = nu as f64 * PI / 2.0;

    if nv == 0 {
        // d^Nu/du^Nu of [v*cos(u)*XD + v*sin(u)*YD + P*u/(2*Pi)*ZD].
        let a_coeff_x = v * (u + a_phase_u).cos();
        let a_coeff_y = v * (u + a_phase_u).sin();
        let mut a_result = a_coeff_x * xd + a_coeff_y * yd;
        if nu == 1 {
            a_result += (h.pitch / (2.0 * PI)) * zd;
        }
        return a_result;
    }

    // Nv == 1: d/dv of d^Nu/du^Nu [v*cos(u)*XD + v*sin(u)*YD]
    // = d^Nu/du^Nu[cos(u)]*XD + d^Nu/du^Nu[sin(u)]*YD.
    let a_coeff_x = (u + a_phase_u).cos();
    let a_coeff_y = (u + a_phase_u).sin();
    a_coeff_x * xd + a_coeff_y * yd
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::SurfaceEval;

    /// The ellipsoid payload with `A = 1.5`, `B = 0.75`, `C = 2.5` on the
    /// identity frame: in the rcad (colatitude) parameterisation
    ///   `P(u, v)  = (A sin v cos u, B sin v sin u, C cos v)`
    ///   `Pu       = (-A sin v sin u, B sin v cos u, 0)`
    ///   `Pv       = (A cos v cos u, B cos v sin u, -C sin v)`
    ///   `Puu      = (-A sin v cos u, -B sin v sin u, 0)`
    ///   `Puv      = (-A cos v sin u, B cos v cos u, 0)`
    ///   `Pvv      = (-A sin v cos u, -B sin v sin u, -C cos v)`
    fn triaxial_ellipsoid() -> EllipsoidalSurface {
        EllipsoidalSurface {
            center: DVec3::new(0.5, -1.0, 2.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 1.5,
            radius_y: 0.75,
            radius_z: 2.5,
        }
    }

    /// `GeomEval_EllipsoidSurface::EvalD1/EvalD2/EvalDN` through the
    /// `V_occt = Pi/2 - V_rcad` conversion: the hand-derived closed form of the
    /// rcad parameterisation must be reproduced exactly, and the point must
    /// agree with `EllipsoidalSurface::point_at` (the parameterisation anchor).
    #[test]
    fn ellipsoid_derivatives_match_the_colatitude_closed_form() {
        let (a, b, c) = (1.5f64, 0.75f64, 2.5f64);
        let el = triaxial_ellipsoid();
        let (u, v) = (0.9f64, 1.1f64);
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let o = el.center;
        let want_p = o + DVec3::new(a * sv * cu, b * sv * su, c * cv);
        let want_d1u = DVec3::new(-a * sv * su, b * sv * cu, 0.0);
        let want_d1v = DVec3::new(a * cv * cu, b * cv * su, -c * sv);
        let want_d2u = DVec3::new(-a * sv * cu, -b * sv * su, 0.0);
        let want_d2uv = DVec3::new(-a * cv * su, b * cv * cu, 0.0);
        let want_d2v = DVec3::new(-a * sv * cu, -b * sv * su, -c * cv);

        // The conversion must land on the payload's own point formula.
        assert!(
            (SurfaceEval::point_at(&el, u, v) - want_p).length() < 1e-14,
            "point_at={:?} want={want_p:?}",
            SurfaceEval::point_at(&el, u, v)
        );

        let d1 = ellipsoid_eval_d1(&el, u, v);
        assert!((d1.point - want_p).length() < 1e-13, "p={:?}", d1.point);
        assert!((d1.d1u - want_d1u).length() < 1e-13, "d1u={:?}", d1.d1u);
        assert!((d1.d1v - want_d1v).length() < 1e-13, "d1v={:?}", d1.d1v);

        let d2 = ellipsoid_eval_d2(&el, u, v);
        assert!((d2.point - want_p).length() < 1e-13);
        assert!((d2.d1u - want_d1u).length() < 1e-13);
        assert!((d2.d1v - want_d1v).length() < 1e-13);
        assert!((d2.d2u - want_d2u).length() < 1e-13, "d2u={:?}", d2.d2u);
        assert!((d2.d2uv - want_d2uv).length() < 1e-13, "d2uv={:?}", d2.d2uv);
        assert!((d2.d2v - want_d2v).length() < 1e-13, "d2v={:?}", d2.d2v);

        assert!((ellipsoid_eval_dn(&el, u, v, 1, 0) - want_d1u).length() < 1e-13);
        assert!((ellipsoid_eval_dn(&el, u, v, 0, 1) - want_d1v).length() < 1e-13);
        assert!((ellipsoid_eval_dn(&el, u, v, 2, 0) - want_d2u).length() < 1e-12);
        assert!((ellipsoid_eval_dn(&el, u, v, 1, 1) - want_d2uv).length() < 1e-12);
        assert!((ellipsoid_eval_dn(&el, u, v, 0, 2) - want_d2v).length() < 1e-12);
    }

    /// The a-sphere limit `A = B = C = R`: the rcad colatitude parameterisation
    /// degenerates to the sphere `P = O + R(sin v cos u, sin v sin u, cos v)`,
    /// whose third-order terms are known in closed form (an oracle independent
    /// of the phase-shift identities the OCCT body uses).
    #[test]
    fn ellipsoid_derivatives_in_the_sphere_limit() {
        let r = 2.0f64;
        let el = EllipsoidalSurface {
            center: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: r,
            radius_y: r,
            radius_z: r,
        };
        let (u, v) = (0.4f64, 0.8f64);
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let want_d3u = DVec3::new(r * sv * su, -r * sv * cu, 0.0);
        let want_d3v = DVec3::new(-r * cv * cu, -r * cv * su, r * sv);
        assert!((ellipsoid_eval_dn(&el, u, v, 3, 0) - want_d3u).length() < 1e-12);
        assert!((ellipsoid_eval_dn(&el, u, v, 0, 3) - want_d3v).length() < 1e-12);
    }

    /// `GeomEval_CircularHelicoidSurface::EvalD1/EvalD2/EvalDN` against the
    /// closed form of `S(u, v) = (v cos u, v sin u, lead*u)` with
    /// `lead = pitch/(2 Pi)`.
    #[test]
    fn helicoid_derivatives_match_the_closed_form() {
        let pitch = 4.0f64;
        let h = HelicoidSurface {
            origin: DVec3::new(0.25, -0.5, 1.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch,
        };
        let lead = pitch / (2.0 * PI);
        let (u, v) = (0.6f64, 1.7f64);
        let (su, cu) = u.sin_cos();
        let o = h.origin;
        let want_p = o + DVec3::new(v * cu, v * su, lead * u);
        let want_d1u = DVec3::new(-v * su, v * cu, lead);
        let want_d1v = DVec3::new(cu, su, 0.0);
        let want_d2u = DVec3::new(-v * cu, -v * su, 0.0);
        let want_d2uv = DVec3::new(-su, cu, 0.0);

        assert!((SurfaceEval::point_at(&h, u, v) - want_p).length() < 1e-14);
        let d1 = helicoid_eval_d1(&h, u, v);
        assert!((d1.point - want_p).length() < 1e-14);
        assert!((d1.d1u - want_d1u).length() < 1e-14, "d1u={:?}", d1.d1u);
        assert!((d1.d1v - want_d1v).length() < 1e-14, "d1v={:?}", d1.d1v);

        let d2 = helicoid_eval_d2(&h, u, v);
        assert!((d2.d1u - want_d1u).length() < 1e-14);
        assert!((d2.d1v - want_d1v).length() < 1e-14);
        assert!((d2.d2u - want_d2u).length() < 1e-14, "d2u={:?}", d2.d2u);
        assert!(d2.d2v.length() < 1e-15);
        assert!((d2.d2uv - want_d2uv).length() < 1e-14, "d2uv={:?}", d2.d2uv);

        assert!((helicoid_eval_dn(&h, u, v, 1, 0) - want_d1u).length() < 1e-14);
        assert!((helicoid_eval_dn(&h, u, v, 0, 1) - want_d1v).length() < 1e-14);
        assert!((helicoid_eval_dn(&h, u, v, 2, 0) - want_d2u).length() < 1e-13);
        assert!((helicoid_eval_dn(&h, u, v, 1, 1) - want_d2uv).length() < 1e-13);
        // The v-dependence is linear: every Nv >= 2 term vanishes.
        assert!(helicoid_eval_dn(&h, u, v, 0, 2).length() < 1e-15);
        assert!(helicoid_eval_dn(&h, u, v, 1, 2).length() < 1e-15);
        // The Z term is linear in u: only Nu == 1 contributes.
        let want_d3u = DVec3::new(v * su, -v * cu, 0.0);
        assert!((helicoid_eval_dn(&h, u, v, 3, 0) - want_d3u).length() < 1e-13);
    }
}
