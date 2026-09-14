//! OCCT Convert package conic conversions (TKMath/Convert) — 1:1
//! statement-mapped translations of the ellipse / hyperbola / parabola
//! conversions and the shared `BuildCosAndSin` machinery.
//!
//! Sources (1:1 statement-mapped):
//! - `Convert_ConicToBSplineCurve.cxx` L217-785 (+ .hxx L42-161): the base
//!   class.  The trimmed-parameterisation `BuildCosAndSin` overload
//!   (L359-622) and its evaluators (L235-355) are the existing
//!   [`super::build_cos_and_sin`] translation — reused, not duplicated.
//!   This file adds the periodic overload (L626-785).
//! - `Convert_EllipseToBSplineCurve.cxx` L47-174 (+ .hxx L38-59)
//! - `Convert_HyperbolaToBSplineCurve.cxx` L32-79 (+ .hxx L37-48)
//! - `Convert_ParabolaToBSplineCurve.cxx` L32-69 (+ .hxx L36-47)
//!
//! Architecture notes (same adaptation as `super::mod.rs`):
//! - The base-class data model (`myPoles 2d / myWeights / myKnots / myMults
//!   / myDegree / myIsPeriodic`) is the existing [`ConvertConicToBspline`]
//!   carrier.
//! - Every `GeomConvert` / `Geom2dConvert` caller constructs the conic in
//!   the CANONICAL `gp::OX2d()` frame (`gp_Elips2d(gp::OX2d(), R, r)` etc.).
//!   For that conic the OCCT handedness factor is `+1` (`value = +r`,
//!   `S = +1`) and the internal `gp_Trsf2d` placement is the identity; the
//!   real frame placement happens in the callers' `BSplineCurveBuilder`.
//!   The statements are kept in canonical form and annotated per function.
//! - `Geom_BSplineCurve::SetPeriodic` (needed by the full-ellipse caller
//!   statement, `GeomConvert.cxx` L377) is carried over the legacy
//!   `BSplineCurve3` carrier at the end of this file.

use std::f64::consts::PI;

use glam::{DVec2, DVec3};

use crate::geom::CurveEval;

use super::{build_cos_and_sin, ConvertConicToBspline, ConvertParameterisation};

// ---------------------------------------------------------------------------
// OCCT Convert_ConicToBSplineCurve.cxx L626-785 — the periodic
// BuildCosAndSin(Parameterisation, CosNumerator, SinNumerator,
// Denominator, Degree, Knots, Mults) overload.
// ---------------------------------------------------------------------------

/// OCCT `BSplCLib::D0(Parameter, UDerivativeRequest = 0, Degree, Periodic =
/// false, Poles, Weights, Knots, Mults, Point)` — the flat-array overload
/// used on the RationalC1 temp curves (Convert_ConicToBSplineCurve.cxx
/// L737-764).  The 1-dim rational form (numerator poles + denominator
/// weights) and the non-rational form (denominator alone, `NoWeights`) are
/// carried through the legacy `BSplineCurve3` carrier (the
/// `cos_and_sin_rational_c1` precedent in `super`): the value is lifted to
/// `(v, 0, 0)` and read back from the x coordinate.
fn temp_curve_d0(
    parameter: f64,
    degree: i32,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
) -> f64 {
    let mut flat = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            flat.push(*k);
        }
    }
    let (control_points, carrier_weights) = match weights {
        Some(w) => (
            poles
                .iter()
                .zip(w.iter())
                .map(|(&p, &wi)| DVec3::new(p * wi, 0.0, 0.0))
                .collect::<Vec<DVec3>>(),
            w.to_vec(),
        ),
        None => (
            poles
                .iter()
                .map(|&p| DVec3::new(p, 0.0, 0.0))
                .collect::<Vec<DVec3>>(),
            vec![1.0; poles.len()],
        ),
    };
    let curve = crate::geom::BSplineCurve3 {
        degree: degree as usize,
        knots: flat,
        control_points,
        weights: carrier_weights,
        is_periodic: false,
    };
    curve.point_at(parameter).x
}

/// OCCT Convert_ConicToBSplineCurve::BuildCosAndSin(Parameterisation,
/// CosNumerator, SinNumerator, Denominator, Degree, Knots, Mults)
/// (Convert_ConicToBSplineCurve.cxx L626-785) — the periodic overload.
/// Returns `(cos_numerator, sin_numerator, denominator, degree, knots,
/// mults)`; the denominator array is the base-class `myWeights` payload
/// (the OCCT callers pass `myWeights` as the `Denominator` output).
#[allow(clippy::type_complexity)]
pub fn build_cos_and_sin_periodic(
    parameterisation: ConvertParameterisation,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, i32, Vec<f64>, Vec<i32>) {
    // OCCT L640-643.
    if parameterisation != ConvertParameterisation::TgtThetaOver2
        && parameterisation != ConvertParameterisation::RationalC1
    {
        panic!("Standard_ConstructionError: Convert_ConicToBSplineCurve::BuildCosAndSin");
    }
    if parameterisation == ConvertParameterisation::TgtThetaOver2 {
        // OCCT L646-656: rebuild on [0, 2*PI] with Convert_TgtThetaOver2_3.
        let (temp_cos, temp_sin, temp_denominator, degree, knots, mut mults) =
            build_cos_and_sin(ConvertParameterisation::TgtThetaOver2_3, 0.0, 2.0 * PI);
        // OCCT L657-665: drop the last (closing) pole triple.
        let num_poles = temp_cos.len() - 1;
        let mut cos_numerator = vec![0.0f64; num_poles];
        let mut sin_numerator = vec![0.0f64; num_poles];
        let mut denominator = vec![0.0f64; num_poles];
        for ii in 0..num_poles {
            cos_numerator[ii] = temp_cos[ii];
            sin_numerator[ii] = temp_sin[ii];
            denominator[ii] = temp_denominator[ii];
        }
        // OCCT L666-669.
        for m in mults.iter_mut() {
            *m = degree;
        }
        (cos_numerator, sin_numerator, denominator, degree, knots, mults)
    } else {
        // OCCT L671-683: Convert_RationalC1 — rebuild on [0, PI]; the temp
        // arrays are the poles/knots of the degree-4 rational curve that
        // the interpolation below re-samples.
        let first_param = 0.0;
        let last_param = PI;
        let (temp_cos, temp_sin, temp_denominator, temp_degree, temp_knots, temp_mults) =
            build_cos_and_sin(ConvertParameterisation::RationalC1, first_param, last_param);

        // OCCT L685-693.
        let degree = 4i32;
        let num_knots = 5i32;
        let num_flat_knots = (degree - 1) * num_knots + 2 * 2;
        let num_poles = num_flat_knots - degree - 1;
        let num_periodic_poles = num_poles - 2;
        let mut flat_knots = vec![0.0f64; num_flat_knots as usize];
        let mut cos_numerator = vec![0.0f64; num_periodic_poles as usize];
        let mut sin_numerator = vec![0.0f64; num_periodic_poles as usize];
        let mut denominator = vec![0.0f64; num_periodic_poles as usize];

        // OCCT L695-715.
        let half_pi = PI * 0.5;
        let mut index = 1i32;
        for _jj in 1..=2 {
            flat_knots[(index - 1) as usize] = -half_pi;
            index += 1;
        }
        for ii in 1..=num_knots {
            for _jj in 1..=(degree - 1) {
                flat_knots[(index - 1) as usize] = (ii - 1) as f64 * half_pi;
                index += 1;
            }
        }
        for _jj in 1..=2 {
            flat_knots[(index - 1) as usize] = 2.0 * PI + half_pi;
            index += 1;
        }
        // OCCT L716-722.
        let mut knots = vec![0.0f64; num_knots as usize];
        let mut mults = vec![0i32; num_knots as usize];
        for ii in 1..=num_knots {
            knots[(ii - 1) as usize] = (ii - 1) as f64 * half_pi;
            mults[(ii - 1) as usize] = degree - 1;
        }

        // OCCT L724-727: Schoenberg points of the flat knot sequence.
        let mut parameters = vec![0.0f64; num_poles as usize];
        crate::math::bspl_lib::build_schoenberg_points(
            degree as usize,
            &flat_knots,
            &mut parameters,
        );
        // OCCT L728: inverse = 1.0 (kept across iterations, as in OCCT).
        let mut inverse = 1.0f64;
        let mut poles_array = vec![0.0f64; (num_poles * 3) as usize];
        let mut contact_order_array = vec![0i32; num_poles as usize];
        for (ii, &param0) in parameters.iter().enumerate() {
            let mut param = param0;
            if param > PI {
                inverse = -1.0;
                param -= PI;
            }
            // OCCT L737-764: the three BSplCLib::D0 evaluations of the temp
            // rational cos / sin curves and the non-rational denominator.
            let value1 = temp_curve_d0(
                param,
                temp_degree,
                &temp_cos,
                Some(&temp_denominator),
                &temp_knots,
                &temp_mults,
            );
            let value2 = temp_curve_d0(
                param,
                temp_degree,
                &temp_sin,
                Some(&temp_denominator),
                &temp_knots,
                &temp_mults,
            );
            let value3 = temp_curve_d0(
                param,
                temp_degree,
                &temp_denominator,
                None,
                &temp_knots,
                &temp_mults,
            );
            contact_order_array[ii] = 0;
            poles_array[ii * 3] = value1 * value3 * inverse;
            poles_array[ii * 3 + 1] = value2 * value3 * inverse;
            poles_array[ii * 3 + 2] = value3;
        }
        // OCCT L771-776.
        let pivot_index_problem = crate::math::bspl_lib::interpolate(
            degree as usize,
            &flat_knots,
            &parameters,
            &contact_order_array,
            3,
            &mut poles_array,
        );
        let _ = pivot_index_problem;
        // OCCT L777-783.
        for ii in 0..num_periodic_poles as usize {
            let inverse = 1.0 / poles_array[ii * 3 + 2];
            cos_numerator[ii] = poles_array[ii * 3] * inverse;
            sin_numerator[ii] = poles_array[ii * 3 + 1] * inverse;
            denominator[ii] = poles_array[ii * 3 + 2];
        }
        (cos_numerator, sin_numerator, denominator, degree, knots, mults)
    }
}

// ---------------------------------------------------------------------------
// OCCT Convert_EllipseToBSplineCurve.cxx
// ---------------------------------------------------------------------------

/// OCCT Convert_EllipseToBSplineCurve(E, Parameterisation)
/// (Convert_EllipseToBSplineCurve.cxx L47-112) — the periodic full ellipse
/// on the canonical `gp::OX2d()` frame (R = theMajorRadius,
/// r = theMinorRadius).
pub fn convert_ellipse_to_bspline_periodic(
    major_radius: f64,
    minor_radius: f64,
    parameterisation: ConvertParameterisation,
) -> ConvertConicToBspline {
    let r = minor_radius;
    let is_periodic;
    // OCCT L55-56: NCollection_Array1<double> CosNumerator, SinNumerator.
    let (cos_numerator, sin_numerator, weights, degree, knots, mults);
    if parameterisation != ConvertParameterisation::TgtThetaOver2
        && parameterisation != ConvertParameterisation::RationalC1
    {
        // OCCT L61-75: if BuildCosAndSin cannot manage the periodicity
        // => trim on 0, 2*PI.
        is_periodic = false;
        let built = build_cos_and_sin(parameterisation, 0.0, 2.0 * PI);
        cos_numerator = built.0;
        sin_numerator = built.1;
        weights = built.2;
        degree = built.3;
        knots = built.4;
        mults = built.5;
    } else {
        // OCCT L76-86.
        is_periodic = true;
        let built = build_cos_and_sin_periodic(parameterisation);
        cos_numerator = built.0;
        sin_numerator = built.1;
        weights = built.2;
        degree = built.3;
        knots = built.4;
        mults = built.5;
    }

    // OCCT L90-101: Ox / Oy are the canonical gp::OX2d() directions
    // (1, 0) / (0, 1): the cross product is +1 > 0, so value = +r, and the
    // Trsf2d built from E.XAxis() (L93) is the identity.  The real frame
    // placement happens in the caller's BSplineCurveBuilder.
    let value = r;

    // OCCT L106-111: poles in the canonical frame, weights from the
    // Denominator output.
    let poles_2d = cos_numerator
        .iter()
        .zip(sin_numerator.iter())
        .map(|(&c, &s)| DVec2::new(major_radius * c, value * s))
        .collect::<Vec<DVec2>>();

    ConvertConicToBspline {
        poles_2d,
        weights,
        knots,
        mults,
        degree,
        is_periodic,
    }
}

/// OCCT Convert_EllipseToBSplineCurve(E, UFirst, ULast, Parameterisation)
/// (Convert_EllipseToBSplineCurve.cxx L119-174) — the non-periodic ellipse
/// arc on the canonical `gp::OX2d()` frame.
pub fn convert_ellipse_arc_to_bspline(
    major_radius: f64,
    minor_radius: f64,
    u_first: f64,
    u_last: f64,
    parameterisation: ConvertParameterisation,
) -> ConvertConicToBspline {
    // OCCT L126-131: Tol = Precision::PConfusion(); Standard_DomainError
    // when delta > 2*PI + Tol or delta <= 0.
    let tol = crate::core::precision::PCONFUSION;
    let delta = u_last - u_first;
    if delta > (2.0 * PI + tol) || delta <= 0.0 {
        panic!("Standard_DomainError: Convert_EllipseToBSplineCurve");
    }

    let r = minor_radius;
    // OCCT L139: myIsPeriodic = false.
    let is_periodic = false;
    // OCCT L140-148.
    let (cos_numerator, sin_numerator, weights, degree, knots, mults) =
        build_cos_and_sin(parameterisation, u_first, u_last);

    // OCCT L152-163: canonical frame — value = +r, identity Trsf2d (see the
    // periodic ctor note above).
    let value = r;

    // OCCT L168-173.
    let poles_2d = cos_numerator
        .iter()
        .zip(sin_numerator.iter())
        .map(|(&c, &s)| DVec2::new(major_radius * c, value * s))
        .collect::<Vec<DVec2>>();

    ConvertConicToBspline {
        poles_2d,
        weights,
        knots,
        mults,
        degree,
        is_periodic,
    }
}

// ---------------------------------------------------------------------------
// OCCT Convert_HyperbolaToBSplineCurve.cxx L26-79
// ---------------------------------------------------------------------------

/// OCCT Convert_HyperbolaToBSplineCurve(H, U1, U2)
/// (Convert_HyperbolaToBSplineCurve.cxx L32-79) — the hyperbola arc on the
/// canonical `gp::OX2d()` frame (R = theMajorRadius, r = theMinorRadius).
pub fn convert_hyperbola_to_bspline(
    major_radius: f64,
    minor_radius: f64,
    u1: f64,
    u2: f64,
) -> ConvertConicToBspline {
    // OCCT L38: Standard_DomainError when |U2 - U1| < Epsilon(0.).
    if (u2 - u1).abs() < crate::base::extrema_ext_elc::epsilon_of(0.0) {
        panic!("Standard_DomainError: Convert_HyperbolaToBSplineCurve");
    }

    // OCCT L40-41.
    let uf = u1.min(u2);
    let ul = u2.max(u1);

    // OCCT L26-28: TheDegree = 2, MaxNbKnots = 2, MaxNbPoles = 3; the base
    // ctor sizes myPoles/myWeights/myKnots/myMults (L36).
    let degree = 2i32;
    let mut knots = vec![0.0f64; 2];
    let mut mults = vec![0i32; 2];

    // OCCT L43: myIsPeriodic = false; L44-47.
    let is_periodic = false;
    knots[0] = uf;
    mults[0] = 3;
    knots[1] = ul;
    mults[1] = 3;

    // OCCT L49-55: construction of the hyperbola in the reference xOy.  For
    // the canonical frame Ox = (1, 0), Oy = (0, 1): the cross product is
    // +1 > 0, so S = +1, and the Trsf2d built from H.Axis().XAxis() (L75)
    // is the identity — the real frame placement happens in the caller's
    // BSplineCurveBuilder.
    let s = 1.0f64;

    let big_r = major_radius;
    let r = minor_radius;

    // OCCT L62-64: the 2nd pole is at the intersection of the 2 tangents to
    // the curve at P(UF), P(UL); its weight is cosh((UL - UF) / 2).
    let weights = vec![1.0, ((ul - uf) / 2.0).cosh(), 1.0];

    // OCCT L66-71: poles expressed in the reference mark.
    let delta = (ul - uf).sinh();
    let x = big_r * (ul.sinh() - uf.sinh()) / delta;
    let y = s * r * (ul.cosh() - uf.cosh()) / delta;
    let poles_2d = vec![
        DVec2::new(big_r * uf.cosh(), s * r * uf.sinh()),
        DVec2::new(x, y),
        DVec2::new(big_r * ul.cosh(), s * r * ul.sinh()),
    ];

    ConvertConicToBspline {
        poles_2d,
        weights,
        knots,
        mults,
        degree,
        is_periodic,
    }
}

// ---------------------------------------------------------------------------
// OCCT Convert_ParabolaToBSplineCurve.cxx L26-69
// ---------------------------------------------------------------------------

/// OCCT Convert_ParabolaToBSplineCurve(Prb, U1, U2)
/// (Convert_ParabolaToBSplineCurve.cxx L32-69) — the parabola arc on the
/// canonical `gp::OX2d()` frame.  `parameter_p` is `Prb.Parameter()` (= 2 x
/// the OCCT focal length; the rcad `Parabola2d::focal_param` /
/// `Parabola3::focal_param` carry the same value).
pub fn convert_parabola_to_bspline(parameter_p: f64, u1: f64, u2: f64) -> ConvertConicToBspline {
    // OCCT L37: Standard_DomainError when |U2 - U1| < Epsilon(0.).
    if (u2 - u1).abs() < crate::base::extrema_ext_elc::epsilon_of(0.0) {
        panic!("Standard_DomainError: Convert_ParabolaToBSplineCurve");
    }

    // OCCT L39-40.
    let uf = u1.min(u2);
    let ul = u2.max(u1);

    // OCCT L42: p = Prb.Parameter() — passed in as `parameter_p`.

    // OCCT L26-28 + L35: TheDegree = 2, MaxNbKnots = 2, MaxNbPoles = 3.
    let degree = 2i32;
    let mut knots = vec![0.0f64; 2];
    let mut mults = vec![0i32; 2];

    // OCCT L44-48.
    let is_periodic = false;
    knots[0] = uf;
    mults[0] = 3;
    knots[1] = ul;
    mults[1] = 3;

    // OCCT L50-52: the parabola is non-rational.
    let weights = vec![1.0, 1.0, 1.0];

    // OCCT L54-56: canonical frame — Ox = (1, 0), Oy = (0, 1), so S = +1,
    // and the Trsf2d built from Prb.Axis().XAxis() (L65) is the identity
    // (the real frame placement happens in the caller's BSplineCurveBuilder).
    let s = 1.0f64;

    // OCCT L58-61: poles expressed in the reference mark.
    let poles_2d = vec![
        DVec2::new((uf * uf) / (2.0 * parameter_p), s * uf),
        DVec2::new((uf * ul) / (2.0 * parameter_p), s * (uf + ul) / 2.0),
        DVec2::new((ul * ul) / (2.0 * parameter_p), s * ul),
    ];

    ConvertConicToBspline {
        poles_2d,
        weights,
        knots,
        mults,
        degree,
        is_periodic,
    }
}

// ---------------------------------------------------------------------------
// OCCT Geom_BSplineCurve::SetPeriodic (Geom_BSplineCurve.cxx L777-815) over
// the legacy BSplineCurve3 carrier — the full-ellipse caller statement
// (GeomConvert.cxx L377 / Geom2dConvert.cxx L386).
// ---------------------------------------------------------------------------

/// OCCT Geom_BSplineCurve::SetPeriodic (Geom_BSplineCurve.cxx L777-815).
/// The knot/mult copy spans `FirstUKnotIndex() .. LastUKnotIndex()` — the
/// curve-level accessors return `1` / `myKnots.Length()` for the carrier's
/// non-periodic data — the end multiplicities are clamped to
/// `min(Degree, max(Mults(1), Mults(n)))`, the poles are truncated to
/// `BSplCLib::NbPoles(Degree, true, Mults)` and the periodic flag is set.
pub(crate) fn bspline3_set_periodic(curve: &mut crate::geom::BSplineCurve3) {
    // OCCT L780-788: first = FirstUKnotIndex(); last = LastUKnotIndex();
    // cknots over [first, last].
    let (knots, mults) = curve.knots_mults();
    let first = 1i32;
    let last = knots.len() as i32;
    let cknots: Vec<f64> = ((first - 1)..last).map(|k| knots[k as usize]).collect();
    // OCCT L790-796: cmults over [first, last] with the end clamp.
    let mut cmults: Vec<i32> = ((first - 1)..last).map(|k| mults[k as usize]).collect();
    let nb = cknots.len();
    let end_mult = (curve.degree as i32).min(cmults[0].max(cmults[nb - 1]));
    cmults[0] = end_mult;
    cmults[nb - 1] = end_mult;

    // OCCT L798-799: compute the new number of poles.
    let nbp = crate::math::bspl_lib::nb_poles(curve.degree, true, &cmults);

    // OCCT L801-809: resize the poles (and weights) keeping the leading
    // nbp entries (UnitWeights for the non-rational case — the carrier
    // always materialises the weight array).
    curve.control_points.truncate(nbp);
    curve.weights.truncate(nbp);

    // OCCT L811: myPeriodic = true (L813-814 bookkeeping omitted: no eval
    // representation on the carrier).
    curve.is_periodic = true;
    curve.knots = cknots
        .iter()
        .zip(cmults.iter())
        .flat_map(|(&k, &m)| std::iter::repeat(k).take(m as usize))
        .collect();
}

// ---------------------------------------------------------------------------
// Tests — expectations hand-derived from the OCCT formulas (never from the
// code output).
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Quarter ellipse, TgtThetaOver2 (delta = PI/2 <= 0.9999*PI: one span,
    /// degree 2, 3 poles).  Closed forms of the L450-469 loop with
    /// UFirst = 0, alpha = PI/4:
    ///   pole 1: (cos 0, sin 0) / 1            = (1, 0), w = 1
    ///   pole 2: (cos PI/4, sin PI/4) / cos PI/4 = (1, 1), w = cos PI/4
    ///   pole 3: (cos PI/2, sin PI/2) / 1      = (0, 1), w = 1
    ///   knots [0, PI/2], mults [3, 3].
    /// Placed on the ellipse R = 3, r = 1: poles (3,0), (3,1), (0,1).
    #[test]
    fn quarter_ellipse_tgt_theta_rational_quadratic() {
        let conv = convert_ellipse_arc_to_bspline(
            3.0,
            1.0,
            0.0,
            std::f64::consts::FRAC_PI_2,
            ConvertParameterisation::TgtThetaOver2,
        );
        assert_eq!(conv.degree, 2);
        assert_eq!(conv.poles_2d.len(), 3);
        assert_eq!(conv.knots.len(), 2);
        let c = std::f64::consts::FRAC_PI_4.cos();
        assert!(conv.poles_2d[0].distance(DVec2::new(3.0, 0.0)) < 1e-14);
        assert!(conv.poles_2d[1].distance(DVec2::new(3.0, 1.0)) < 1e-14);
        assert!(conv.poles_2d[2].distance(DVec2::new(0.0, 1.0)) < 1e-14);
        assert!((conv.weights[1] - c).abs() < 1e-14, "w2 = cos(PI/4)");
        assert_eq!(conv.weights[0], 1.0);
        assert_eq!(conv.weights[2], 1.0);
        assert!((conv.knots[0]).abs() < 1e-14);
        assert!((conv.knots[1] - std::f64::consts::FRAC_PI_2).abs() < 1e-14);
        assert_eq!(conv.mults, vec![3, 3]);
        assert!(!conv.is_periodic);

        // The rational quadratic evaluates exactly on the ellipse
        // (x/3)^2 + y^2 = 1.  At u = PI/4 the Bernstein parameter is
        // t = 1/2 and with N = (1/4, 1/2, 1/4), w = (1, c, 1):
        //   x = 3(1/2 + c)/(1 + c) = 3/sqrt(2),  y = (1/2 + c)/(1 + c)
        //     = 1/sqrt(2).
        let s2 = 2.0f64.sqrt();
        let (px, py) = (3.0 / s2, 1.0 / s2);
        let conv2 = crate::geom::BSplineCurve2 {
            degree: conv.degree as usize,
            knots: crate::geom::BSplineCurve3::from_knots_mults(
                conv.degree as usize,
                conv.knots.clone(),
                conv.mults.clone(),
                Vec::new(),
            )
            .knots,
            control_points: conv.poles_2d.clone(),
            weights: conv.weights.clone(),
        };
        let pt = crate::math::bspl::de_boor_2d(
            conv2.degree,
            &conv2.knots,
            &conv2.control_points,
            &conv2.weights,
            std::f64::consts::FRAC_PI_4,
        );
        assert!((pt.x - px).abs() < 1e-12 && (pt.y - py).abs() < 1e-12);
        assert!(((pt.x / 3.0).powi(2) + pt.y * pt.y - 1.0).abs() < 1e-12);
    }

    /// Full ellipse, TgtThetaOver2 (periodic overload L646-669): 6 poles
    /// (the 7th closing pole dropped), 4 knots over [0, 2PI], all mults
    /// forced to Degree = 2, is_periodic = true.  Three 120-degree spans,
    /// alpha = PI/3: weight-1 poles on the arc at 0/120/240 degrees,
    /// interior poles at 60/180/300 degrees scaled by 1/cos(PI/3) = 2 with
    /// weight cos(PI/3) = 0.5.  Ellipse R = 2, r = 5.
    #[test]
    fn full_ellipse_periodic_tgt_theta() {
        let conv =
            convert_ellipse_to_bspline_periodic(2.0, 5.0, ConvertParameterisation::TgtThetaOver2);
        assert!(conv.is_periodic);
        assert_eq!(conv.degree, 2);
        assert_eq!(conv.poles_2d.len(), 6);
        assert_eq!(conv.knots.len(), 4);
        assert_eq!(conv.mults, vec![2, 2, 2, 2]);
        // Knots: 3 equal spans over [0, 2PI].
        for (i, k) in conv.knots.iter().enumerate() {
            assert!((k - i as f64 * 2.0 * PI / 3.0).abs() < 1e-14, "knot {i}");
        }
        let sq3 = 3.0f64.sqrt();
        // Weight-1 poles (CosNumerator odd entries, Denominator = 1):
        //   1: (cos 0, sin 0)                  -> (2, 0)
        //   3: (cos 120, sin 120)              -> (-1, 5*sqrt3/2)
        //   5: (cos 240, sin 240)              -> (-1, -5*sqrt3/2)
        assert!(conv.poles_2d[0].distance(DVec2::new(2.0, 0.0)) < 1e-14);
        assert!(conv.poles_2d[2].distance(DVec2::new(-1.0, 5.0 * sq3 / 2.0)) < 1e-14);
        assert!(conv.poles_2d[4].distance(DVec2::new(-1.0, -5.0 * sq3 / 2.0)) < 1e-14);
        assert_eq!(conv.weights[0], 1.0);
        assert_eq!(conv.weights[2], 1.0);
        assert_eq!(conv.weights[4], 1.0);
        // Interior poles (scaled by inverse = 2, weight = direct = 0.5):
        //   2: 2*(cos 60, sin 60)              -> (2, 5*sqrt3)
        //   4: 2*(cos 180, sin 180)            -> (-4, 0)
        //   6: 2*(cos 300, sin 300)            -> (2, -5*sqrt3)
        assert!(conv.poles_2d[1].distance(DVec2::new(2.0, 5.0 * sq3)) < 1e-14);
        assert!(conv.poles_2d[3].distance(DVec2::new(-4.0, 0.0)) < 1e-14);
        assert!(conv.poles_2d[5].distance(DVec2::new(2.0, -5.0 * sq3)) < 1e-14);
        assert!((conv.weights[1] - 0.5).abs() < 1e-14);
        assert!((conv.weights[3] - 0.5).abs() < 1e-14);
        assert!((conv.weights[5] - 0.5).abs() < 1e-14);
    }

    /// Parabola p = 2, [-1, 3] (L59-61 closed forms; UF = min, UL = max):
    ///   P1 = (1/(2*2), -1) = (0.25, -1), P2 = (-3/(2*2), 1) = (-0.75, 1),
    ///   P3 = (9/(2*2), 3) = (2.25, 3), all weights 1, degree 2,
    ///   knots [-1, 3] mults [3, 3].  The deg-2 polynomial interpolates the
    ///   analytic parabola at the three Bezier ordinates: at t = 1
    ///   (midpoint of [-1, 3] in the Bernstein sense) the quadratic through
    ///   (P1, P2, P3) gives (t^2/(2p), t) with t = 1.
    #[test]
    fn parabola_deg2_poles() {
        let conv = convert_parabola_to_bspline(2.0, 3.0, -1.0);
        assert_eq!(conv.degree, 2);
        assert!(conv.poles_2d[0].distance(DVec2::new(0.25, -1.0)) < 1e-14);
        assert!(conv.poles_2d[1].distance(DVec2::new(-0.75, 1.0)) < 1e-14);
        assert!(conv.poles_2d[2].distance(DVec2::new(2.25, 3.0)) < 1e-14);
        assert_eq!(conv.weights, vec![1.0, 1.0, 1.0]);
        assert!((conv.knots[0] - (-1.0)).abs() < 1e-14);
        assert!((conv.knots[1] - 3.0).abs() < 1e-14);
        assert_eq!(conv.mults, vec![3, 3]);
        assert!(!conv.is_periodic);
    }

    /// Hyperbola R = 3, r = 2, [-1, 1] (L62-71 closed forms):
    ///   w2 = cosh(1), P1 = (3 cosh 1, -2 sinh 1), P3 = (3 cosh 1, 2 sinh 1),
    ///   P2 = (3(sinh 1 - sinh(-1))/sinh 2, 2(cosh 1 - cosh(-1))/sinh 2)
    ///       = (3*2*sinh 1 / (2 sinh 1 cosh 1), 0) = (3 / cosh 1, 0).
    #[test]
    fn hyperbola_rational_quadratic_poles() {
        let conv = convert_hyperbola_to_bspline(3.0, 2.0, -1.0, 1.0);
        assert_eq!(conv.degree, 2);
        let ch = 1.0f64.cosh();
        let sh = 1.0f64.sinh();
        assert!((conv.weights[1] - ch).abs() < 1e-14);
        assert!(conv.poles_2d[0].distance(DVec2::new(3.0 * ch, -2.0 * sh)) < 1e-14);
        assert!(conv.poles_2d[1].distance(DVec2::new(3.0 / ch, 0.0)) < 1e-14);
        assert!(conv.poles_2d[2].distance(DVec2::new(3.0 * ch, 2.0 * sh)) < 1e-14);
        assert_eq!(conv.mults, vec![3, 3]);
    }

    /// SetPeriodic over the carrier (Geom_BSplineCurve.cxx L777-815): a
    /// clamped degree-6 curve with 7 poles and mults [7, 7] becomes a
    /// periodic curve with the end mults clamped to min(6, 7) = 6 and the
    /// poles truncated to NbPoles(6, true, [6, 6]) = 6.
    #[test]
    fn set_periodic_clamps_end_mults_and_truncates_poles() {
        let mut curve = crate::geom::BSplineCurve3 {
            degree: 6,
            knots: vec![0.0; 7].into_iter().chain(vec![1.0; 7]).collect(),
            control_points: (0..7).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect(),
            weights: vec![1.0; 7],
            is_periodic: false,
        };
        bspline3_set_periodic(&mut curve);
        assert!(curve.is_periodic);
        assert_eq!(curve.control_points.len(), 6);
        assert_eq!(curve.knots.len(), 12);
        assert!(curve.knots.iter().take(6).all(|&k| k == 0.0));
        assert!(curve.knots.iter().skip(6).all(|&k| k == 1.0));
    }

    // -- Caller-level checks through the retired Geom2dConvert /
    //    GeomConvert branches (expectations hand-derived from the OCCT
    //    caller + builder statements).

    /// Geom2dConvert.cxx L314-321 + BSplineCurveBuilder L70-96: the trimmed
    /// parabola p = 2 on [-1, 3], apex (10, 20), axis (1, 0).  Canonical
    /// poles (U^2/(2p), U): (0.25, -1), (-0.75, 1), (2.25, 3); the frame
    /// (axis, turn_2d(axis)) is right-handed so no mirror; placed poles:
    /// (10.25, 19), (9.25, 21), (12.25, 23).  At the trim start t = -1 the
    /// curve evaluates to the apex-relative point (t^2/(2p), t) + origin =
    /// (10.25, 19); at t = 0 (Bernstein parameter 1/4, weights all 1) the
    /// quadratic gives exactly the apex (10, 20).
    #[test]
    fn caller_2d_trimmed_parabola_placed_poles() {
        use crate::base::geom2d_convert::curve_to_bspline_curve_2d;
        use crate::geom::{Curve2d, Parabola2d, TrimmedCurve2};
        let prb = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Parabola(Parabola2d {
                origin: DVec2::new(10.0, 20.0),
                axis_dir: DVec2::new(1.0, 0.0),
                focal_param: 2.0,
            })),
            t_min: -1.0,
            t_max: 3.0,
        });
        let bs = curve_to_bspline_curve_2d(&prb, ConvertParameterisation::TgtThetaOver2);
        assert_eq!(bs.degree, 2);
        assert_eq!(bs.control_points.len(), 3);
        assert!(bs.control_points[0].distance(DVec2::new(10.25, 19.0)) < 1e-14);
        assert!(bs.control_points[1].distance(DVec2::new(9.25, 21.0)) < 1e-14);
        assert!(bs.control_points[2].distance(DVec2::new(12.25, 23.0)) < 1e-14);
        assert_eq!(bs.weights, vec![1.0, 1.0, 1.0]);
        // Flat knots: [-1 x3, 3 x3].
        assert_eq!(bs.knots, vec![-1.0, -1.0, -1.0, 3.0, 3.0, 3.0]);
    }

    /// Geom2dConvert.cxx L266-303 + BSplineCurveBuilder L70-96: the
    /// canonical quarter-ellipse poles (R,0), (R,r) w = cos(PI/4), (0,r)
    /// placed in a ROTATED right-handed conic frame X = (3/5, 4/5),
    /// Y = (-4/5, 3/5) on center (1, 2): Loc + x*X + y*Y (the composition
    /// of the handedness mirror and T.SetTransformation(XAxis, OX2d),
    /// gp_Trsf2d.cxx L48-70, collapses to the frame placement).
    /// P1 = (1+9/5, 2+12/5) = (2.8, 4.4); P2 = (1+9/5-4/5, 2+12/5+3/5) =
    /// (2, 5); P3 = (1-4/5, 2+3/5) = (0.2, 2.6).
    #[test]
    fn caller_2d_trimmed_ellipse_rotated_frame() {
        use crate::base::geom2d_convert::curve_to_bspline_curve_2d;
        use crate::geom::{Curve2d, Ellipse2d, TrimmedCurve2};
        let el = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Ellipse(Ellipse2d {
                center: DVec2::new(1.0, 2.0),
                major_dir: DVec2::new(0.6, 0.8),
                minor_dir: DVec2::new(-0.8, 0.6),
                major_radius: 3.0,
                minor_radius: 1.0,
            })),
            t_min: 0.0,
            t_max: std::f64::consts::FRAC_PI_2,
        });
        let bs = curve_to_bspline_curve_2d(&el, ConvertParameterisation::TgtThetaOver2);
        assert_eq!(bs.degree, 2);
        assert!(bs.control_points[0].distance(DVec2::new(2.8, 4.4)) < 1e-14);
        assert!(bs.control_points[1].distance(DVec2::new(2.0, 5.0)) < 1e-14);
        assert!(bs.control_points[2].distance(DVec2::new(0.2, 2.6)) < 1e-14);
        assert!((bs.weights[1] - std::f64::consts::FRAC_PI_4.cos()).abs() < 1e-14);
        // Flat knots: [0 x3, PI/2 x3].
        let hp = std::f64::consts::FRAC_PI_2;
        assert_eq!(bs.knots, vec![0.0, 0.0, 0.0, hp, hp, hp]);
    }

    /// GeomConvert.cxx L284-290 + BSplineCurveBuilder L56-81: the trimmed
    /// hyperbola R = 3, r = 2 on [-1, 1], center (1,2,3), frame (X, Z x X).
    /// Canonical poles (R cosh U, r sinh U) with the middle pole
    /// (R (sinh UL - sinh UF) / sinh(UL-UF), 0) = (3 / cosh 1, 0) and
    /// weight cosh((UL-UF)/2) = cosh 1.
    #[test]
    fn caller_3d_trimmed_hyperbola_placed_poles() {
        use crate::base::convert::geom_convert_curve_to_bspline_curve;
        use crate::geom::{Curve3, Hyperbola3, TrimmedCurve3};
        let hyp = Curve3::Trimmed(TrimmedCurve3::new(
            Curve3::Hyperbola(Hyperbola3 {
                center: DVec3::new(1.0, 2.0, 3.0),
                normal: DVec3::Z,
                major_dir: DVec3::X,
                semi_major: 3.0,
                semi_minor: 2.0,
            }),
            -1.0,
            1.0,
        ));
        let bs = geom_convert_curve_to_bspline_curve(&hyp, ConvertParameterisation::TgtThetaOver2);
        assert_eq!(bs.degree, 2);
        let ch = 1.0f64.cosh();
        let sh = 1.0f64.sinh();
        assert!(bs.control_points[0].distance(DVec3::new(1.0 + 3.0 * ch, 2.0 - 2.0 * sh, 3.0)) < 1e-14);
        assert!(bs.control_points[1].distance(DVec3::new(1.0 + 3.0 / ch, 2.0, 3.0)) < 1e-14);
        assert!(bs.control_points[2].distance(DVec3::new(1.0 + 3.0 * ch, 2.0 + 2.0 * sh, 3.0)) < 1e-14);
        assert!((bs.weights[1] - ch).abs() < 1e-14);
        assert!(!bs.is_periodic);
    }

    /// GeomConvert.cxx L363-385: the FULL ellipse with TgtThetaOver2 —
    /// the periodic 6-pole conversion (see full_ellipse_periodic_tgt_theta)
    /// placed in the conic plane, then SetPeriodic (no-op recomputation on
    /// the already-periodic data: knots [0, 2PI/3, 4PI/3, 2PI] mults 2,
    /// NbPoles(2, true, [2,2,2,2]) = 6).
    #[test]
    fn caller_3d_full_ellipse_periodic() {
        use crate::base::convert::geom_convert_curve_to_bspline_curve;
        use crate::geom::{Curve3, Ellipse3};
        let el = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 2.0,
            minor_radius: 5.0,
        });
        let bs = geom_convert_curve_to_bspline_curve(&el, ConvertParameterisation::TgtThetaOver2);
        assert!(bs.is_periodic);
        assert_eq!(bs.degree, 2);
        assert_eq!(bs.control_points.len(), 6);
        let sq3 = 3.0f64.sqrt();
        let k = 2.0 * PI / 3.0;
        assert_eq!(
            bs.knots,
            vec![0.0, 0.0, k, k, 2.0 * k, 2.0 * k, 2.0 * PI, 2.0 * PI]
        );
        assert!(bs.control_points[0].distance(DVec3::new(2.0, 0.0, 0.0)) < 1e-14);
        assert!(bs.control_points[2].distance(DVec3::new(-1.0, 5.0 * sq3 / 2.0, 0.0)) < 1e-14);
        assert!(bs.control_points[3].distance(DVec3::new(-4.0, 0.0, 0.0)) < 1e-14);
        assert!((bs.weights[1] - 0.5).abs() < 1e-14);
    }
}
