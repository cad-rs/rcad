//! The analytic Geom2d_OffsetCurve D0/D1/D2/D3/DN evaluation over the rcad
//! [`OffsetCurve2d`] carrier — 1:1 translation:
//!
//! - `Geom2d_OffsetCurve::EvalD0/EvalD1/EvalD2/EvalD3/EvalDN`
//!   (TKG2d/Geom2d/Geom2d_OffsetCurve.cxx L214-357) — the basis-derivative
//!   preamble, the `AdjustDerivative` singular-point guard and the
//!   `Geom2d_OffsetCurveUtils` dispatch;
//! - `Geom2d_OffsetCurveUtils::CalculateD0/CalculateD1/CalculateD2/CalculateD3`
//!   (TKG2d/Geom2d/Geom2d_OffsetCurveUtils.pxx L43-280) — the IICURV
//!   (EUCLID-IS) offset formulas in terms of the rotated basis derivatives
//!   `Ndir = (D1.Y, -D1.X)`, `DNdir = (D2.Y, -D2.X)`, ...;
//! - `Geom2d_OffsetCurveUtils::AdjustDerivative`
//!   (Geom2d_OffsetCurveUtils.pxx L295-364) — the Taylor-series fallback at
//!   singular points where the basis first derivative vanishes.
//!
//! Architecture notes (rcad carrier vs OCCT object):
//! 1. OCCT bundles (Point, D1[, D2[, D3]]) into `Geom2d_Curve::ResD1/2/3`
//!    returns of the basis `EvalD1/EvalD2/EvalD3`; the rcad
//!    [`Curve2dEval`] trait splits them into `point_at` /
//!    `derivative_at` / `derivative2_at` / `derivative3_at`.  The BSpline
//!    basis answers all three exactly (`bspline2d_dn`), the conics
//!    analytically, so the split carries the same exact data.
//! 2. The OCCT exceptions `Geom2d_UndefinedValue` /
//!    `Geom2d_UndefinedDerivative` raised when the offset direction is
//!    undefined (zero tangent) map to panics carrying the exception names
//!    (the same mapping as `curve_dn`).
//! 3. `AdjustDerivative` calls the basis `EvalDN` up to order 6; the rcad
//!    [`Curve2dEval::derivative_n_at`] default answers N > 3 by a finite
//!    difference of the exact D3 (the 2D per-type DN engines beyond D3 are
//!    not re-hosted yet).  The path is reached only at singular points
//!    where the basis first derivative is exactly zero.
//! 4. OCCT quirk preserved: the published higher-order terms of
//!    `CalculateD3` (pxx L252-256 else-branch: the `- 3 * DNdir * Dr^2/R5`
//!    grouping and the unscaled `- D3r` in the Ndir coefficient) deviate
//!    from the exact third derivative of the offset whenever
//!    `Dr = Ndir.DNdir != 0` — OCCT's own EvalD2, differenced centrally,
//!    yields the exact D3, not the OCCT EvalD3 value.  The formulas are
//!    translated verbatim; the deviation vanishes where Dr = D2r = D3r = 0
//!    (e.g. at the quarter-arc symmetry point).

use glam::DVec2;

use crate::core::precision::{REAL_FIRST, REAL_LAST};
use crate::geom::bspline2d_dn::{ResD1, ResD2, ResD3};
use crate::geom::{Curve2d, Curve2dEval, OffsetCurve2d};

/// OCCT `gp::Resolution()` = `RealSmall()` = `DBL_MIN`
/// (gp.hxx L59-60) — the zero-tangent guard of the offset formulas.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT Geom2d_OffsetCurveUtils::CalculateD0 (Geom2d_OffsetCurveUtils.pxx
/// L43-53) — the offset point from the basis D1; `false` is the OCCT
/// `return false` on the zero tangent (raises Geom2d_UndefinedValue at the
/// caller).
pub fn calculate_d0(the_value: &mut DVec2, the_d1: DVec2, the_offset: f64) -> bool {
    if the_d1.length_squared() <= GP_RESOLUTION {
        return false;
    }
    // gp_Dir2d aNormal(theD1.Y(), -theD1.X()) — the Dir ctor normalizes.
    let a_normal = DVec2::new(the_d1.y, -the_d1.x).normalize();
    *the_value = *the_value + a_normal * the_offset;
    true
}

/// OCCT Geom2d_OffsetCurveUtils::CalculateD1 (Geom2d_OffsetCurveUtils.pxx
/// L61-101) — the offset point and first derivative; `false` is the OCCT
/// failure return (R2 degenerate).
pub fn calculate_d1(
    the_value: &mut DVec2,
    the_d1: &mut DVec2,
    the_d2: DVec2,
    the_offset: f64,
) -> bool {
    let mut n_dir = DVec2::new(the_d1.y, -the_d1.x);
    let mut dn_dir = DVec2::new(the_d2.y, -the_d2.x);
    let r2 = n_dir.length_squared();
    let r = r2.sqrt();
    let r3 = r * r2;
    let dr = n_dir.dot(dn_dir);
    if r3 <= GP_RESOLUTION {
        if r2 <= GP_RESOLUTION {
            return false;
        }
        // We try another computation but the stability is not very good.
        dn_dir = dn_dir * r;
        dn_dir = dn_dir - n_dir * (dr / r);
        dn_dir = dn_dir * (the_offset / r2);
    } else {
        // Same computation as IICURV in EUCLID-IS because the stability is
        // better.
        dn_dir = dn_dir * (the_offset / r);
        dn_dir = dn_dir - n_dir * (the_offset * dr / r3);
    }

    n_dir = n_dir * (the_offset / r);
    // P(u)
    *the_value = *the_value + n_dir;
    // P'(u)
    *the_d1 = *the_d1 + dn_dir;
    true
}

/// OCCT Geom2d_OffsetCurveUtils::CalculateD2 (Geom2d_OffsetCurveUtils.pxx
/// L111-178) — the offset point and first/second derivatives.
#[allow(clippy::too_many_arguments)]
pub fn calculate_d2(
    the_value: &mut DVec2,
    the_d1: &mut DVec2,
    the_d2: &mut DVec2,
    the_d3: DVec2,
    the_is_dir_change: bool,
    the_offset: f64,
) -> bool {
    let mut n_dir = DVec2::new(the_d1.y, -the_d1.x);
    let mut dn_dir = DVec2::new(the_d2.y, -the_d2.x);
    let mut d2n_dir = DVec2::new(the_d3.y, -the_d3.x);
    let r2 = n_dir.length_squared();
    let r = r2.sqrt();
    let r3 = r2 * r;
    let r4 = r2 * r2;
    let r5 = r3 * r2;
    let dr = n_dir.dot(dn_dir);
    let d2r = n_dir.dot(d2n_dir) + dn_dir.dot(dn_dir);
    if r5 <= GP_RESOLUTION {
        if r4 <= GP_RESOLUTION {
            return false;
        }
        // We try another computation but the stability is not very good
        // dixit ISG.
        // V2 = P" (U) :
        d2n_dir = d2n_dir - dn_dir * (2.0 * dr / r2);
        d2n_dir = d2n_dir + n_dir * ((3.0 * dr * dr) / r4 - d2r / r2);
        d2n_dir = d2n_dir * (the_offset / r);

        // V1 = P' (U) :
        dn_dir = dn_dir * r;
        dn_dir = dn_dir - n_dir * (dr / r);
        dn_dir = dn_dir * (the_offset / r2);
    } else {
        // Same computation as IICURV in EUCLID-IS because the stability is
        // better.
        // V2 = P" (U) :
        d2n_dir = d2n_dir * (the_offset / r);
        d2n_dir = d2n_dir - dn_dir * (2.0 * the_offset * dr / r3);
        d2n_dir = d2n_dir + n_dir * (the_offset * ((3.0 * dr * dr) / r5 - d2r / r3));

        // V1 = P' (U)
        dn_dir = dn_dir * (the_offset / r);
        dn_dir = dn_dir - n_dir * (the_offset * dr / r3);
    }

    n_dir = n_dir * (the_offset / r);
    // P(u)
    *the_value = *the_value + n_dir;
    // P'(u) :
    *the_d1 = *the_d1 + dn_dir;
    // P"(u) :
    if the_is_dir_change {
        *the_d2 = -*the_d2;
    }
    *the_d2 = *the_d2 + d2n_dir;
    true
}

/// OCCT Geom2d_OffsetCurveUtils::CalculateD3 (Geom2d_OffsetCurveUtils.pxx
/// L189-280) — the offset point and first/second/third derivatives.  The
/// degenerate-branch `R4 = R2 * R2` re-assignment and the sign quirks of
/// the published formulas are kept verbatim.
#[allow(clippy::too_many_arguments)]
pub fn calculate_d3(
    the_value: &mut DVec2,
    the_d1: &mut DVec2,
    the_d2: &mut DVec2,
    the_d3: &mut DVec2,
    the_d4: DVec2,
    the_is_dir_change: bool,
    the_offset: f64,
) -> bool {
    let mut n_dir = DVec2::new(the_d1.y, -the_d1.x);
    let mut dn_dir = DVec2::new(the_d2.y, -the_d2.x);
    let mut d2n_dir = DVec2::new(the_d3.y, -the_d3.x);
    let mut d3n_dir = DVec2::new(the_d4.y, -the_d4.x);
    let r2 = n_dir.length_squared();
    let r = r2.sqrt();
    let r3 = r2 * r;
    let mut r4 = r2 * r2;
    let r5 = r3 * r2;
    let r6 = r3 * r3;
    let r7 = r5 * r2;
    let dr = n_dir.dot(dn_dir);
    let d2r = n_dir.dot(d2n_dir) + dn_dir.dot(dn_dir);
    let d3r = n_dir.dot(d3n_dir) + 3.0 * dn_dir.dot(d2n_dir);

    if r7 <= GP_RESOLUTION {
        if r6 <= GP_RESOLUTION {
            return false;
        }
        // We try another computation but the stability is not very good
        // dixit ISG.
        // V3 = P"' (U) :
        d3n_dir = d3n_dir - d2n_dir * (3.0 * dr / r2);
        d3n_dir = d3n_dir - dn_dir * (3.0 * (d2r / r2 + dr * dr / r4));
        d3n_dir = d3n_dir
            + n_dir * (6.0 * dr * dr / r4 + 6.0 * dr * d2r / r4 - 15.0 * dr * dr * dr / r6 - d3r);
        d3n_dir = d3n_dir * (the_offset / r);
        // V2 = P" (U) :
        r4 = r2 * r2;
        d2n_dir = d2n_dir - dn_dir * (2.0 * dr / r2);
        d2n_dir = d2n_dir - n_dir * (3.0 * dr * dr / r4 - d2r / r2);
        d2n_dir = d2n_dir * (the_offset / r);
        // V1 = P' (U) :
        dn_dir = dn_dir * r;
        dn_dir = dn_dir - n_dir * (dr / r);
        dn_dir = dn_dir * (the_offset / r2);
    } else {
        // Same computation as IICURV in EUCLID-IS because the stability is
        // better.
        // V3 = P"' (U) :
        d3n_dir = d3n_dir * (the_offset / r);
        d3n_dir = d3n_dir - d2n_dir * (3.0 * the_offset * dr / r3);
        d3n_dir = d3n_dir - dn_dir * (3.0 * the_offset * (d2r / r3 + dr * dr / r5));
        d3n_dir = d3n_dir
            + n_dir * (the_offset * (6.0 * dr * dr / r5 + 6.0 * dr * d2r / r5 - 15.0 * dr * dr * dr / r7 - d3r));
        // V2 = P" (U) :
        d2n_dir = d2n_dir * (the_offset / r);
        d2n_dir = d2n_dir - dn_dir * (2.0 * the_offset * dr / r3);
        d2n_dir = d2n_dir - n_dir * (the_offset * ((3.0 * dr * dr) / r5 - d2r / r3));
        // V1 = P' (U) :
        dn_dir = dn_dir * (the_offset / r);
        dn_dir = dn_dir - n_dir * (the_offset * dr / r3);
    }

    n_dir = n_dir * (the_offset / r);
    // P(u)
    *the_value = *the_value + n_dir;
    // P'(u) :
    *the_d1 = *the_d1 + dn_dir;
    // P"(u)
    *the_d2 = *the_d2 + d2n_dir;
    // P"'(u)
    if the_is_dir_change {
        *the_d3 = -*the_d3;
    }
    *the_d3 = *the_d3 + d3n_dir;
    true
}

/// OCCT Geom2d_OffsetCurveUtils::AdjustDerivative
/// (Geom2d_OffsetCurveUtils.pxx L295-364) — the Taylor-series fallback at
/// singular points where the basis first derivative is nearly zero; the
/// basis curve is the rcad [`Curve2d`] value behind the
/// `Geom2d_Curve::FirstParameter/LastParameter/EvalD0/EvalDN` calls.
#[allow(clippy::too_many_arguments)]
fn adjust_derivative(
    the_curve: &Curve2d,
    the_max_derivative: i32,
    the_u: f64,
    the_d1: &mut DVec2,
    the_d2: &mut DVec2,
    the_d3: &mut DVec2,
    the_d4: &mut DVec2,
    the_is_direction_change: &mut bool,
) -> bool {
    // OCCT L305-307: static consts.
    let a_tol = GP_RESOLUTION;
    let a_min_step = 1.0e-7;
    let a_max_deriv_order = 3;

    *the_is_direction_change = false;
    // OCCT L310-311: FirstParameter / LastParameter of the basis curve.
    let domain = the_curve.default_domain();
    let (an_uinfium, an_usupremum) = (domain[0], domain[1]);

    // OCCT L313-324: DivisionFactor and the finite-domain guard.
    const DIVISION_FACTOR: f64 = 1.0e-3;
    let du;
    if (an_usupremum >= REAL_LAST) || (an_uinfium <= REAL_FIRST) {
        du = 0.0;
    } else {
        du = an_usupremum - an_uinfium;
    }
    let a_delta = f64::max(du * DIVISION_FACTOR, a_min_step);

    // OCCT L327-333: the do/while Taylor climb — EvalDN with a growing
    // order until a non-degenerate derivative is found (capped at 3).
    let mut an_index = 1; // Derivative order
    let mut v;
    loop {
        an_index += 1;
        v = the_curve.derivative_n_at(the_u, an_index);
        if !(v.length_squared() <= a_tol && an_index < a_max_deriv_order) {
            break;
        }
    }

    // OCCT L335-344: the evaluation step side.
    let u = if the_u - an_uinfium < a_delta {
        the_u + a_delta
    } else {
        the_u - a_delta
    };

    // OCCT L346-352: EvalD0 on both sides, the direction-change probe.
    let a_p1 = the_curve.point_at(f64::min(the_u, u));
    let a_p2 = the_curve.point_at(f64::max(the_u, u));

    let v1 = a_p2 - a_p1;
    *the_is_direction_change = v.dot(v1) < 0.0;
    let a_sign = if *the_is_direction_change { -1.0 } else { 1.0 };

    // OCCT L355-361: the sign-corrected derivatives; aDeriv[3] maps to
    // (the_d2, the_d3, the_d4) in call order.
    *the_d1 = v * a_sign;
    for i in 1..the_max_derivative {
        let a_dn = the_curve.derivative_n_at(the_u, an_index + i);
        let a_deriv: &mut DVec2 = match i {
            1 => &mut *the_d2,
            2 => &mut *the_d3,
            _ => &mut *the_d4,
        };
        *a_deriv = a_dn * a_sign;
    }

    true
}

/// OCCT Geom2d_OffsetCurve::EvalD0 (Geom2d_OffsetCurve.cxx L214-229) — the
/// basis point plus `CalculateD0` over the basis D1 (the rcad trait splits
/// the OCCT `EvalD1` bundle into `point_at` + `derivative_at`; architecture
/// note 1).
pub fn eval_d0(the_c: &OffsetCurve2d, the_u: f64) -> DVec2 {
    // OCCT L222-223: basisCurve->EvalD1(theU), aValue = aBasisD1.Point.
    let a_basis_d1 = the_c.basis.derivative_at(the_u);
    let mut a_value = the_c.basis.point_at(the_u);
    // OCCT L224-227.
    if !calculate_d0(&mut a_value, a_basis_d1, the_c.offset_distance) {
        panic!("Geom2d_UndefinedValue: Geom2d_OffsetCurve::EvalD0");
    }
    a_value
}

/// OCCT Geom2d_OffsetCurve::EvalD1 (Geom2d_OffsetCurve.cxx L233-249) — the
/// basis D2 bundle through `CalculateD1`.
pub fn eval_d1(the_c: &OffsetCurve2d, the_u: f64) -> ResD1 {
    // OCCT L241-243: basisCurve->EvalD2(theU); aValue = aBasisD2.Point;
    // aD1 = aBasisD2.D1.
    let mut a_value = the_c.basis.point_at(the_u);
    let mut a_d1 = the_c.basis.derivative_at(the_u);
    let a_d2 = the_c.basis.derivative2_at(the_u);
    // OCCT L244-247.
    if !calculate_d1(&mut a_value, &mut a_d1, a_d2, the_c.offset_distance) {
        panic!("Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalD1");
    }
    ResD1 {
        point: a_value,
        d1: a_d1,
    }
}

/// OCCT Geom2d_OffsetCurve::EvalD2 (Geom2d_OffsetCurve.cxx L253-285) — the
/// basis D3 bundle, the `AdjustDerivative(3)` guard on the zero basis D1,
/// then `CalculateD2`.
pub fn eval_d2(the_c: &OffsetCurve2d, the_u: f64) -> ResD2 {
    // OCCT L261-263: basisCurve->EvalD3(theU); the point/derivative split.
    let mut a_value = the_c.basis.point_at(the_u);
    let mut a_d1 = the_c.basis.derivative_at(the_u);
    let mut a_d2 = the_c.basis.derivative2_at(the_u);
    let mut a_d3 = the_c.basis.derivative3_at(the_u);
    // OCCT L264: bool isDirectionChange = false.
    let mut is_direction_change = false;
    // OCCT L265-279: the singular-point guard.
    if a_d1.length_squared() <= GP_RESOLUTION {
        let mut a_dummy_d4 = DVec2::ZERO;
        if !adjust_derivative(
            &the_c.basis,
            3,
            the_u,
            &mut a_d1,
            &mut a_d2,
            &mut a_d3,
            &mut a_dummy_d4,
            &mut is_direction_change,
        ) {
            panic!("Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalD2");
        }
    }
    // OCCT L280-283.
    if !calculate_d2(
        &mut a_value,
        &mut a_d1,
        &mut a_d2,
        a_d3,
        is_direction_change,
        the_c.offset_distance,
    ) {
        panic!("Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalD2");
    }
    ResD2 {
        point: a_value,
        d1: a_d1,
        d2: a_d2,
    }
}

/// OCCT Geom2d_OffsetCurve::EvalD3 (Geom2d_OffsetCurve.cxx L289-328) — the
/// basis D3 bundle plus the basis `EvalDN(U, 4)`, the
/// `AdjustDerivative(4)` guard, then `CalculateD3`.
pub fn eval_d3(the_c: &OffsetCurve2d, the_u: f64) -> ResD3 {
    // OCCT L297-301: basisCurve->EvalD3(theU) and basisCurve->EvalDN(U, 4).
    let mut a_value = the_c.basis.point_at(the_u);
    let mut a_d1 = the_c.basis.derivative_at(the_u);
    let mut a_d2 = the_c.basis.derivative2_at(the_u);
    let mut a_d3 = the_c.basis.derivative3_at(the_u);
    let mut a_d4 = the_c.basis.derivative_n_at(the_u, 4);
    // OCCT L302: bool isDirectionChange = false.
    let mut is_direction_change = false;
    // OCCT L303-316: the singular-point guard.
    if a_d1.length_squared() <= GP_RESOLUTION {
        if !adjust_derivative(
            &the_c.basis,
            4,
            the_u,
            &mut a_d1,
            &mut a_d2,
            &mut a_d3,
            &mut a_d4,
            &mut is_direction_change,
        ) {
            panic!("Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalD3");
        }
    }
    // OCCT L317-326.
    if !calculate_d3(
        &mut a_value,
        &mut a_d1,
        &mut a_d2,
        &mut a_d3,
        a_d4,
        is_direction_change,
        the_c.offset_distance,
    ) {
        panic!("Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalD3");
    }
    ResD3 {
        point: a_value,
        d1: a_d1,
        d2: a_d2,
        d3: a_d3,
    }
}

/// OCCT Geom2d_OffsetCurve::EvalDN (Geom2d_OffsetCurve.cxx L332-357) —
/// N = 1/2/3 delegate to the D1/D2/D3 bodies, higher orders return the
/// basis `EvalDN` (the offset contribution is negligible there).
pub fn eval_dn(the_c: &OffsetCurve2d, the_u: f64, the_n: i32) -> DVec2 {
    if the_n < 1 {
        panic!("Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalDN");
    }
    match the_n {
        1 => eval_d1(the_c, the_u).d1,
        2 => eval_d2(the_c, the_u).d2,
        3 => eval_d3(the_c, the_u).d3,
        _ => the_c.basis.derivative_n_at(the_u, the_n),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// The rational quadratic quarter circle (radius 1, from (1,0) to
    /// (0,1), CCW) on knots [0,0,0,1,1,1] with weights (1, sqrt2/2, 1).
    fn make_quarter_circle() -> crate::geom::BSplineCurve2 {
        let s2 = std::f64::consts::SQRT_2;
        crate::geom::BSplineCurve2 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec2::new(1.0, 0.0),
                DVec2::new(1.0, 1.0),
                DVec2::new(0.0, 1.0),
            ],
            weights: vec![1.0, s2 / 2.0, 1.0],
            is_periodic: false,
        }
    }

    fn offset_of(basis: crate::geom::BSplineCurve2, d: f64) -> OffsetCurve2d {
        OffsetCurve2d {
            basis: Box::new(Curve2d::BSpline(basis)),
            offset_distance: d,
        }
    }

    /// (a) Offset of the straight degree-1 BSpline line P(t) = (4t, 0) by
    /// d — the offset is the parallel line (4t, -d): D1 = (4, 0) is exactly
    /// the direction, D2 and D3 carry no offset term and vanish exactly.
    /// The old finite-difference D2 produced cancellation noise far above
    /// the 1e-9 tolerance asserted here.
    #[test]
    fn straight_bspline_line_offset_d1_exact_d2_zero() {
        let d = 2.5;
        let line = crate::geom::BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(4.0, 0.0)],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        };
        let off = offset_of(line, d);
        for k in 0..=5 {
            let t = k as f64 / 5.0;
            // D0 through the analytic eval_d0 (not the 2-point tangent).
            let p = off.point_at(t);
            let y_want = -d;
            assert!(
                (p - DVec2::new(4.0 * t, y_want)).length() < 1.0e-12,
                "P({t}) = {:?}, want ({}, {y_want})",
                p,
                4.0 * t
            );
            // D1 exactly the direction.
            let d1 = off.derivative_at(t);
            assert!(
                (d1 - DVec2::new(4.0, 0.0)).length() < 1.0e-12,
                "D1({t}) = {:?}, want (4, 0)",
                d1
            );
            // The offset contributes exactly zero to D2/D3 on a line.
            let d2 = off.derivative2_at(t);
            assert!(
                d2.length() < 1.0e-9,
                "D2({t}) = {:?}, want zero",
                d2
            );
            let d3 = off.derivative3_at(t);
            assert!(
                d3.length() < 1.0e-9,
                "D3({t}) = {:?}, want zero",
                d3
            );
        }
    }

    /// Closed-form reference of the unit rational quarter circle and of the
    /// angle derivatives (theta', theta'', theta''') of its parametrization.
    ///
    /// With the homogeneous power forms X = 1 + p t + q t^2,
    /// Y = a t + q t^2, W = 1 + p t + r t^2 (a = sqrt2, p = a - 2,
    /// q = 1 - a, r = 2 - a), the curve lies on the unit circle so
    /// X^2 + Y^2 = W^2, and
    /// theta' = (X Y' - Y X') / W^2 = (a + 2 q (t - t^2)) / W^2;
    /// theta'' = (g' W - 2 g W') / W^3 and
    /// theta''' = (h' W - 3 h W') / W^4 follow by the quotient rule.
    fn circle_theta_derivs(t: f64) -> (DVec2, f64, f64, f64) {
        let a = std::f64::consts::SQRT_2;
        let p = a - 2.0;
        let q = 1.0 - a;
        let r = 2.0 - a;
        let big_x = 1.0 + p * t + q * t * t;
        let big_y = a * t + q * t * t;
        let w = 1.0 + p * t + r * t * t;
        let point = DVec2::new(big_x / w, big_y / w);
        let g = a + 2.0 * q * (t - t * t);
        let g1 = 2.0 * q * (1.0 - 2.0 * t);
        let g2 = -4.0 * q;
        let w1 = p + 2.0 * r * t;
        let w2 = 2.0 * r;
        let th1 = g / (w * w);
        let h = g1 * w - 2.0 * g * w1;
        let h1 = g2 * w - g1 * w1 - 2.0 * g * w2;
        let th2 = h / (w * w * w);
        let th3 = (h1 * w - 3.0 * h * w1) / (w * w * w * w);
        (point, th1, th2, th3)
    }

    /// (b) Offset of the rational quarter-circle BSpline by d — the offset
    /// is the circle of radius 1 + d under the same angle parametrization,
    /// so with A = (-y, x), B = (-x, -y):
    ///   P_off = (1+d)(x, y),
    ///   D1_off = (1+d) th' A,
    ///   D2_off = (1+d)(th'' A + th'^2 B),
    ///   D3_off = (1+d)((th''' - th'^3) A + 3 th' th'' B).
    ///
    /// The D3 closed form is asserted only at t = 0.5, where the symmetry
    /// of the quarter arc gives Dr = D2r = D3r = 0 and the preserved OCCT
    /// quirk of `CalculateD3` (architecture note 4) is inactive; away from
    /// it the OCCT 8.0 formula deviates from its own D2 (verified by a
    /// central difference of the translated EvalD2), so the faithful
    /// engine output there is the published OCCT result, not the calculus.
    #[test]
    fn rational_circle_offset_matches_closed_form() {
        for &d in &[0.7_f64, -0.3] {
            let off = offset_of(make_quarter_circle(), d);
            let radius = 1.0 + d;
            for k in 0..=4 {
                let t = k as f64 / 4.0;
                let (c, th1, th2, th3) = circle_theta_derivs(t);
                let a_vec = DVec2::new(-c.y, c.x);
                let b_vec = DVec2::new(-c.x, -c.y);
                let p_ref = c * radius;
                let d1_ref = a_vec * (th1 * radius);
                let d2_ref = a_vec * (th2 * radius) + b_vec * (th1 * th1 * radius);

                let p = off.point_at(t);
                assert!(
                    (p - p_ref).length() < 1.0e-12,
                    "P({t}, d={d}) = {:?}, want {:?}",
                    p,
                    p_ref
                );
                let d1 = off.derivative_at(t);
                assert!(
                    (d1 - d1_ref).length() < 1.0e-12,
                    "D1({t}, d={d}) = {:?}, want {:?}",
                    d1,
                    d1_ref
                );
                let d2 = off.derivative2_at(t);
                assert!(
                    (d2 - d2_ref).length() < 1.0e-10,
                    "D2({t}, d={d}) = {:?}, want {:?}",
                    d2,
                    d2_ref
                );
                // The offset of the unit circle is the radius-(1+d) circle.
                assert!((p.length() - radius).abs() < 1.0e-12, "radius at {t}");

                if k == 2 {
                    // t = 0.5: Dr = D2r = D3r = 0 — the OCCT formula is the
                    // exact calculus and the closed form applies.  The
                    // residual tolerance covers the basis DN(4) finite
                    // difference (architecture note 3) feeding D3Ndir/D3r.
                    let d3_ref = a_vec * ((th3 - th1 * th1 * th1) * radius)
                        + b_vec * (3.0 * th1 * th2 * radius);
                    let d3 = off.derivative3_at(t);
                    assert!(
                        (d3 - d3_ref).length() < 1.0e-6,
                        "D3({t}, d={d}) = {:?}, want {:?}",
                        d3,
                        d3_ref
                    );
                }
            }
        }
    }

    /// The singular-point guard: the cubic with poles P0 = P1 = P2 = (0,0),
    /// P3 = (1,0) has x(t) = t^3, so at U = 0 the basis D1 = D2 = 0 exactly
    /// and D3 = (6, 0); AdjustDerivative climbs to the third derivative and
    /// the offset evaluation continues through the Taylor form: the offset
    /// of (t^3, 0) by d along the adjusted normal keeps D1 = (6, 0),
    /// D2 = D3 = 0 at the cusp.
    #[test]
    fn singular_cusp_adjust_derivative_path() {
        let basis = crate::geom::BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(0.0, 0.0),
                DVec2::new(0.0, 0.0),
                DVec2::new(1.0, 0.0),
            ],
            weights: vec![1.0; 4],
            is_periodic: false,
        };
        let off = offset_of(basis, 1.0);
        // EvalD2 through the singular point: the adjusted D1 and the
        // point offset along the adjusted normal (0, -1).
        let r2 = eval_d2(&off, 0.0);
        assert!((r2.point - DVec2::new(0.0, -1.0)).length() < 1.0e-12, "P(0)");
        assert!((r2.d1 - DVec2::new(6.0, 0.0)).length() < 1.0e-12, "D1(0)");
        assert!((r2.d2 - DVec2::new(0.0, 0.0)).length() < 1.0e-12, "D2(0)");
        // EvalD3 through the singular point (the order-4 basis DN).
        let r3 = eval_d3(&off, 0.0);
        assert!((r3.d1 - DVec2::new(6.0, 0.0)).length() < 1.0e-12, "D1(0)");
        assert!((r3.d2 - DVec2::new(0.0, 0.0)).length() < 1.0e-12, "D2(0)");
        assert!((r3.d3 - DVec2::new(0.0, 0.0)).length() < 1.0e-12, "D3(0)");
    }

    /// OCCT guard asymmetry (cxx L233-249): EvalD1 has no AdjustDerivative
    /// fallback — at the cusp the zero basis tangent raises
    /// Geom2d_UndefinedDerivative.
    #[test]
    #[should_panic(expected = "Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalD1")]
    fn eval_d1_raises_at_the_cusp() {
        let basis = crate::geom::BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(0.0, 0.0),
                DVec2::new(0.0, 0.0),
                DVec2::new(1.0, 0.0),
            ],
            weights: vec![1.0; 4],
            is_periodic: false,
        };
        let off = offset_of(basis, 1.0);
        off.derivative_at(0.0);
    }

    /// EvalDN (cxx L332-357): N = 1/2/3 return the D1/D2/D3 bodies, N = 4
    /// returns the basis fourth derivative (zero on the straight line),
    /// N < 1 raises Geom2d_UndefinedDerivative.
    #[test]
    fn eval_dn_dispatch() {
        let line = crate::geom::BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(4.0, 0.0)],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        };
        let off = offset_of(line, 2.5);
        let t = 0.3;
        assert!((eval_dn(&off, t, 1) - off.derivative_at(t)).length() < 1.0e-12);
        assert!((eval_dn(&off, t, 2) - off.derivative2_at(t)).length() < 1.0e-12);
        assert!((eval_dn(&off, t, 3) - off.derivative3_at(t)).length() < 1.0e-12);
        assert!(eval_dn(&off, t, 4).length() < 1.0e-12);
    }

    #[test]
    #[should_panic(expected = "Geom2d_UndefinedDerivative: Geom2d_OffsetCurve::EvalDN")]
    fn eval_dn_rejects_n_below_one() {
        let line = crate::geom::BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(4.0, 0.0)],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        };
        let off = offset_of(line, 2.5);
        eval_dn(&off, 0.5, 0);
    }

    /// The zero-tangent guard of EvalD0: the degenerate constant basis
    /// (D1 = 0 everywhere) raises Geom2d_UndefinedValue.
    #[test]
    #[should_panic(expected = "Geom2d_UndefinedValue: Geom2d_OffsetCurve::EvalD0")]
    fn eval_d0_rejects_zero_tangent() {
        let basis = crate::geom::BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![DVec2::new(1.0, 1.0); 4],
            weights: vec![1.0; 4],
            is_periodic: false,
        };
        let off = offset_of(basis, 1.0);
        off.point_at(0.5);
    }
}
