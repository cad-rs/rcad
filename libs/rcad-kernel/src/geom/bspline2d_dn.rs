//! The exact Geom2d_BSplineCurve D0/D1/D2/D3 evaluation over the rcad
//! [`BSplineCurve2`] carrier — the BSplCLib statement chain, translated
//! 1:1:
//!
//! - `Geom2d_BSplineCurve::EvalD0/EvalD1/EvalD2/EvalD3`
//!   (TKG2d/Geom2d/Geom2d_BSplineCurve_1.cxx L175-304) — the span location
//!   preamble plus the `BSplCLib::D0/D1/D2/D3` (knots, mults) dispatch;
//! - `Geom2d_BSplineCurve::PeriodicNormalization`
//!   (Geom2d_BSplineCurve.cxx L1328-1346);
//! - `BSplCLib_D0/D1/D2/D3<gp_Pnt2d, gp_Vec2d, ..., 2>`
//!   (TKMath/BSplCLib/BSplCLib_CurveComputation.pxx L843-1084) with the
//!   inlined `PrepareEval_T<..., 2>` (L777-852) — LocateParameter,
//!   BuildKnots, PoleIndex, the span-local IsRational probe, BuildEval,
//!   then BSplCLib::Bohm and the rational PLib::RationalDerivative pass;
//! - `BSplCLib::LocateParameter(Degree, Knots, Mults, U, Periodic, ...)`
//!   (BSplCLib.cxx L321-356);
//! - `BSplCLib::IsRational` (BSplCLib.cxx L842-859).
//!
//! Architecture notes (rcad carrier vs OCCT object):
//! 1. `BSplineCurve2` stores the FLAT (multiplicity-expanded) knot vector
//!    while the OCCT D-kernels read the (Knots, Mults) pair; the pair is
//!    reconstructed by run-length compression (the same bijection as
//!    `geom::bspline_ops` / `geom2d_convert::comp_curve_to_bspline_2d`).
//! 2. `BSplineCurve2` carries no periodic flag — every OCCT `myPeriodic`
//!    site below passes `false` (the same architecture annotation as
//!    `geom2d_convert::comp_curve_to_bspline_2d`).
//! 3. OCCT passes `Weights()` (a null handle for a non-rational curve);
//!    the rcad carrier always stores a weight vector (all 1.0 when
//!    non-rational).  The span-local `IsRational` probe inside
//!    `PrepareEval_T` makes the two forms take identical branches.

use glam::DVec2;

use crate::geom::BSplineCurve2;
use crate::math::bspl_lib::{
    at, ati, bohm, bspl_clib_eval_inplace, build_knots_local, first_uknot_index_mults,
    last_uknot_index_mults, locate_parameter_main, pole_index,
};
use crate::math::plib::rational_derivative;

/// OCCT Geom2d_Curve::ResD1 (Geom2d_Curve.hxx L63-67) — point and first
/// derivative.
pub struct ResD1 {
    pub point: DVec2,
    pub d1: DVec2,
}

/// OCCT Geom2d_Curve::ResD2 (Geom2d_Curve.hxx L70-75).
pub struct ResD2 {
    pub point: DVec2,
    pub d1: DVec2,
    pub d2: DVec2,
}

/// OCCT Geom2d_Curve::ResD3 (Geom2d_Curve.hxx L78-84).
pub struct ResD3 {
    pub point: DVec2,
    pub d1: DVec2,
    pub d2: DVec2,
    pub d3: DVec2,
}

/// The (Knots, Multiplicities) pair of an rcad flat knot vector
/// (run-length compression; OCCT stores the pair directly).
fn knots_mults_of(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots = Vec::new();
    let mut mults = Vec::new();
    for (i, k) in flat.iter().enumerate() {
        if i > 0 && *k == knots.last().copied().unwrap_or(f64::NAN) {
            *mults.last_mut().expect("non-empty") += 1;
        } else {
            knots.push(*k);
            mults.push(1);
        }
    }
    (knots, mults)
}

/// OCCT BSplCLib::IsRational(Weights, I1, I2) (BSplCLib.cxx L842-859) —
/// the consecutive-weight probe over [I1, I2]; a single differing
/// neighbouring pair marks the span rational.  `weights` is 1-based
/// addressed through `at`-style indexing (Lower == 1).
fn bsplclib_is_rational(weights: &[f64], i1: i32, i2: i32) -> bool {
    // OCCT L844-848: f = Weights.Lower() (1), l = Weights.Length(),
    // I3 = I2 - f.
    let f = 1i32;
    let l = weights.len() as i32;
    let i3 = i2 - f;
    let mut i = i1 - f;
    while i < i3 {
        // OCCT L851: WG[f + (i % l)] != WG[f + ((i + 1) % l)] — with the
        // WG pointer folding the 1-based Lower offset, the element read is
        // Weights(f + (i % l)) == 0-based weights[i % l].
        if at(weights, f + (i % l)) != at(weights, f + ((i + 1) % l)) {
            return true;
        }
        i += 1;
    }
    false
}

/// OCCT BSplCLib::LocateParameter(Degree, Knots, Mults, U, Periodic,
/// KnotIndex, NewU) (BSplCLib.cxx L321-356) — the (knots, mults) form:
/// first/last from FirstUKnotIndex/LastUKnotIndex, the inner locate runs
/// only when the incoming index is outside [first, last], otherwise
/// NewU = U.
fn locate_parameter_bspline(
    degree: i32,
    knots: &[f64],
    mults: &[i32],
    u: f64,
    is_periodic: bool,
    knot_index: &mut i32,
    new_u: &mut f64,
) {
    let first = first_uknot_index_mults(degree as usize, mults);
    let last = last_uknot_index_mults(degree as usize, mults);
    if *knot_index < first || *knot_index > last {
        // OCCT L351-355: LocateParameter(Knots, U, Periodic, first, last,
        // KnotIndex, NewU, Knots(first), Knots(last)).
        locate_parameter_main(
            knots,
            u,
            is_periodic,
            first,
            last,
            knot_index,
            new_u,
            at(knots, first),
            at(knots, last),
        );
    } else {
        *new_u = u;
    }
}

/// OCCT PrepareEval_T<gp_Pnt2d, gp_Vec2d, NCollection_Array1<gp_Pnt2d>, 2>
/// (BSplCLib_CurveComputation.pxx L777-852) — normalizes/locates the span,
/// builds the local knot window, converts the span index to the pole
/// cursor, probes the span-local rationality and copies the (Degree+1)
/// span poles (homogenized when rational) into the local array.
/// Returns `(dim, rational, dc_knots, dc_poles)`.
#[allow(clippy::too_many_arguments)]
fn prepare_eval_2d(
    u: &mut f64,
    index: &mut i32,
    degree: i32,
    is_periodic: bool,
    poles: &[DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
) -> (usize, bool, Vec<f64>, Vec<f64>) {
    // OCCT L798: LocateParameter(Degree, Knots, Mults, u, Periodic, index, u).
    locate_parameter_bspline(degree, knots, mults, *u, is_periodic, index, u);

    // OCCT L801: BuildKnots(Degree, index, Periodic, Knots, Mults, *dc.knots).
    let deg = degree as usize;
    let mut dc_knots = vec![0.0f64; 2 * deg];
    build_knots_local(deg, *index, is_periodic, knots, Some(mults), &mut dc_knots);

    // OCCT L806: Mults != nullptr — index = PoleIndex(Degree, index,
    // Periodic, *Mults).
    *index = pole_index(deg, *index, is_periodic, mults);

    // OCCT L809-816: rational = (Weights != nullptr) and the span-local
    // IsRational probe with WLower = Weights.Lower() + index.
    let mut rational = weights.is_some();
    if rational {
        let w_lower = 1 + *index; // Weights.Lower() == 1.
        rational = bsplclib_is_rational(weights.expect("weights"), w_lower, w_lower + degree);
    }

    // OCCT L818-838: BSplCLib_BuildEval<..., 2> (L719-754) — the
    // (Degree+1) local poles, homogeneous (stride 3) when rational.  The
    // cursor `ip = Poles.Lower() + Index - 1` (Lower == 1) with the
    // wrap `if (ip > PUpper) ip = PLower` is kept.
    let dim = if rational { 3 } else { 2 };
    let mut dc_poles = vec![0.0f64; (deg + 1) * dim];
    let poles_upper = poles.len() as i32;
    let mut ip = 1 + *index - 1;
    let mut ptr = 0usize;
    for _i in 0..=deg {
        ip += 1;
        if ip > poles_upper {
            ip = 1;
        }
        let p = poles[(ip - 1) as usize];
        if rational {
            let w = weights.expect("weights")[(ip - 1) as usize];
            dc_poles[ptr] = p.x * w;
            dc_poles[ptr + 1] = p.y * w;
            dc_poles[ptr + 2] = w;
            ptr += 3;
        } else {
            dc_poles[ptr] = p.x;
            dc_poles[ptr + 1] = p.y;
            ptr += 2;
        }
    }
    (dim, rational, dc_knots, dc_poles)
}

/// OCCT validateBSplineDegree (BSplCLib_CurveComputation.pxx L254-258) —
/// THE_MAX_DEGREE == BSplCLib::MaxDegree() == 25.
fn validate_bspline_degree(degree: i32) {
    if degree > 25 {
        panic!("Standard_OutOfRange: BSplCLib: bspline degree is greater than maximum supported");
    }
}

/// OCCT BSplCLib::D0<gp_Pnt2d, ..., 2> (BSplCLib_CurveComputation.pxx
/// L843-892) — Eval on the local window, then the homogeneous division for
/// the rational case.
#[allow(clippy::too_many_arguments)]
fn bsplclib_d0(
    u: f64,
    index: i32,
    degree: i32,
    is_periodic: bool,
    poles: &[DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
) -> DVec2 {
    let mut index = index;
    let mut u = u;
    validate_bspline_degree(degree);
    let (dim, rational, dc_knots, mut dc_poles) =
        prepare_eval_2d(&mut u, &mut index, degree, is_periodic, poles, weights, knots, mults);
    // OCCT: BSplCLib::Eval(u, Degree, *dc.knots, dim, *dc.poles).
    bspl_clib_eval_inplace(u, degree, &dc_knots, dim, &mut dc_poles);

    if rational {
        // OCCT: double w = dc.poles[Dimension]; CoordsToPointScaled — the
        // homogeneous division.
        let w = dc_poles[2];
        DVec2::new(dc_poles[0] / w, dc_poles[1] / w)
    } else {
        DVec2::new(dc_poles[0], dc_poles[1])
    }
}

/// OCCT BSplCLib::D1<gp_Pnt2d, ..., 2> (BSplCLib_CurveComputation.pxx
/// L896-948) — Bohm with N = 1, the rational PLib::RationalDerivative pass,
/// then P and V from the (DerivativeRequest+1) blocks of stride Dimension.
#[allow(clippy::too_many_arguments)]
fn bsplclib_d1(
    u: f64,
    index: i32,
    degree: i32,
    is_periodic: bool,
    poles: &[DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
) -> ResD1 {
    let mut index = index;
    let mut u = u;
    validate_bspline_degree(degree);
    let (dim, rational, dc_knots, mut dc_poles) =
        prepare_eval_2d(&mut u, &mut index, degree, is_periodic, poles, weights, knots, mults);
    // OCCT: BSplCLib::Bohm(u, Degree, 1, *dc.knots, dim, *dc.poles).
    bohm(u, degree, 1, &dc_knots, dim, &mut dc_poles);
    let mut dc_ders;
    // OCCT: double* result = dc.poles.
    let mut result: &[f64] = &dc_poles;
    if rational {
        // OCCT: PLib::RationalDerivative(Degree, 1, Dimension, *dc.poles,
        // *dc.ders) — All defaulted true.
        dc_ders = vec![0.0f64; 2 * 2];
        rational_derivative(degree, 1, 2, &dc_poles, &mut dc_ders, true);
        result = &dc_ders;
    }

    ResD1 {
        point: DVec2::new(result[0], result[1]),
        d1: DVec2::new(result[2], result[3]),
    }
}

/// OCCT BSplCLib::D2<gp_Pnt2d, ..., 2> (BSplCLib_CurveComputation.pxx
/// L950-1012) — Bohm with N = 2; the non-rational Degree < 2 nullification
/// of V2 is kept verbatim.
#[allow(clippy::too_many_arguments)]
fn bsplclib_d2(
    u: f64,
    index: i32,
    degree: i32,
    is_periodic: bool,
    poles: &[DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
) -> ResD2 {
    let mut index = index;
    let mut u = u;
    validate_bspline_degree(degree);
    let (dim, rational, dc_knots, mut dc_poles) =
        prepare_eval_2d(&mut u, &mut index, degree, is_periodic, poles, weights, knots, mults);
    // OCCT: BSplCLib::Bohm(u, Degree, 2, *dc.knots, dim, *dc.poles).
    bohm(u, degree, 2, &dc_knots, dim, &mut dc_poles);
    let mut dc_ders;
    let mut result: &[f64] = &dc_poles;
    if rational {
        // OCCT: PLib::RationalDerivative(Degree, 2, Dimension, ...).
        dc_ders = vec![0.0f64; 3 * 2];
        rational_derivative(degree, 2, 2, &dc_poles, &mut dc_ders, true);
        result = &dc_ders;
    }

    let point = DVec2::new(result[0], result[1]);
    let v1 = DVec2::new(result[2], result[3]);
    let v2 = if !rational && degree < 2 {
        DVec2::ZERO
    } else {
        DVec2::new(result[4], result[5])
    };
    ResD2 { point, d1: v1, d2: v2 }
}

/// OCCT BSplCLib::D3<gp_Pnt2d, ..., 2> (BSplCLib_CurveComputation.pxx
/// L1014-1083) — Bohm with N = 3; the non-rational Degree < 2 / Degree < 3
/// nullifications of V2 / V3 are kept verbatim.
#[allow(clippy::too_many_arguments)]
fn bsplclib_d3(
    u: f64,
    index: i32,
    degree: i32,
    is_periodic: bool,
    poles: &[DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
) -> ResD3 {
    let mut index = index;
    let mut u = u;
    validate_bspline_degree(degree);
    let (dim, rational, dc_knots, mut dc_poles) =
        prepare_eval_2d(&mut u, &mut index, degree, is_periodic, poles, weights, knots, mults);
    // OCCT: BSplCLib::Bohm(u, Degree, 3, *dc.knots, dim, *dc.poles).
    bohm(u, degree, 3, &dc_knots, dim, &mut dc_poles);
    let mut dc_ders;
    let mut result: &[f64] = &dc_poles;
    if rational {
        // OCCT: PLib::RationalDerivative(Degree, 3, Dimension, ...).
        dc_ders = vec![0.0f64; 4 * 2];
        rational_derivative(degree, 3, 2, &dc_poles, &mut dc_ders, true);
        result = &dc_ders;
    }

    let point = DVec2::new(result[0], result[1]);
    let v1 = DVec2::new(result[2], result[3]);
    let v2 = if !rational && degree < 2 {
        DVec2::ZERO
    } else {
        DVec2::new(result[4], result[5])
    };
    let v3 = if !rational && degree < 3 {
        DVec2::ZERO
    } else {
        DVec2::new(result[6], result[7])
    };
    ResD3 {
        point,
        d1: v1,
        d2: v2,
        d3: v3,
    }
}

/// OCCT Geom2d_BSplineCurve::PeriodicNormalization(Parameter)
/// (Geom2d_BSplineCurve.cxx L1328-1346) — wraps the parameter into the
/// flat-knot period.  The carrier carries no periodic flag (architecture
/// note): `is_periodic` is false from every caller, matching the OCCT
/// `if (myPeriodic)` guard on the non-periodic curves the carrier holds.
fn periodic_normalization(bs: &BSplineCurve2, is_periodic: bool, parameter: &mut f64) {
    if is_periodic {
        let flat = &bs.knots;
        let deg = bs.degree as i32;
        let k_upper = flat.len() as i32;
        // OCCT L1332: Period = myFlatKnots.Value(Upper - myDeg)
        // - myFlatKnots.Value(myDeg + 1).
        let period = at(flat, k_upper - deg) - at(flat, deg + 1);
        // OCCT L1333-1339.
        while *parameter > at(flat, k_upper - deg) {
            *parameter -= period;
        }
        // OCCT L1340-1345.
        while *parameter < at(flat, deg + 1) {
            *parameter += period;
        }
    }
}

/// The carrier weights in the OCCT `Weights()` position — always present
/// on the rcad carrier (unit weights when non-rational; see the module
/// architecture note 3).
fn curve_weights(bs: &BSplineCurve2) -> Option<&[f64]> {
    Some(&bs.weights[..])
}

/// OCCT Geom2d_BSplineCurve::EvalD0(U) (Geom2d_BSplineCurve_1.cxx
/// L175-201) — PeriodicNormalization, LocateParameter, the span adjust and
/// the BSplCLib::D0 dispatch.
pub fn eval_d0(bs: &BSplineCurve2, u: f64) -> DVec2 {
    let mut a_span_index = 0i32;
    let mut a_new_u = u;
    // OCCT L192: PeriodicNormalization(aNewU).
    periodic_normalization(bs, false, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L193: LocateParameter(myDeg, myKnots, &myMults, U, myPeriodic,
    // aSpanIndex, aNewU).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, false, &mut a_span_index, &mut a_new_u);
    // OCCT L195-198.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L200: BSplCLib::D0(...).
    bsplclib_d0(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        false,
        &bs.control_points,
        curve_weights(bs),
        &knots,
        &mults,
    )
}

/// OCCT Geom2d_BSplineCurve::EvalD1(U) (Geom2d_BSplineCurve_1.cxx
/// L199-226).
pub fn eval_d1(bs: &BSplineCurve2, u: f64) -> ResD1 {
    let mut a_span_index = 0i32;
    let mut a_new_u = u;
    // OCCT L207: PeriodicNormalization(aNewU).
    periodic_normalization(bs, false, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L208: LocateParameter(...).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, false, &mut a_span_index, &mut a_new_u);
    // OCCT L210-213.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L215-224: BSplCLib::D1(...).
    bsplclib_d1(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        false,
        &bs.control_points,
        curve_weights(bs),
        &knots,
        &mults,
    )
}

/// OCCT Geom2d_BSplineCurve::EvalD2(U) (Geom2d_BSplineCurve_1.cxx
/// L232-264).
pub fn eval_d2(bs: &BSplineCurve2, u: f64) -> ResD2 {
    let mut a_span_index = 0i32;
    let mut a_new_u = u;
    // OCCT L240: PeriodicNormalization(aNewU).
    periodic_normalization(bs, false, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L241: LocateParameter(...).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, false, &mut a_span_index, &mut a_new_u);
    // OCCT L243-246.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L248-257: BSplCLib::D2(...).
    bsplclib_d2(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        false,
        &bs.control_points,
        curve_weights(bs),
        &knots,
        &mults,
    )
}

/// OCCT Geom2d_BSplineCurve::EvalD3(U) (Geom2d_BSplineCurve_1.cxx
/// L266-299).
pub fn eval_d3(bs: &BSplineCurve2, u: f64) -> ResD3 {
    let mut a_span_index = 0i32;
    let mut a_new_u = u;
    // OCCT L274: PeriodicNormalization(aNewU).
    periodic_normalization(bs, false, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L275: LocateParameter(...).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, false, &mut a_span_index, &mut a_new_u);
    // OCCT L277-280.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L282-291: BSplCLib::D3(...).
    bsplclib_d3(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        false,
        &bs.control_points,
        curve_weights(bs),
        &knots,
        &mults,
    )
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// A degree-2 non-rational BSpline of two parabola segments joined C1 at
    /// u = 1:
    ///   segment [0, 1]: poles (0,0), (1,2), (2,0) — P(u) = (2u, 4u(1-u)),
    ///     D1 = (2, 4-8u), D2 = (0, -8) constant;
    ///   segment [1, 2]: poles (2,0), (3,-2), (4,0) — with v = u-1,
    ///     P = (2+2v, -4v(1-v)), D1 = (2, -4+8v), D2 = (0, +8) constant.
    /// The D2 jump at the interior knot makes the finite-difference
    /// default non-discriminative there and at the clamped ends.
    fn two_segment_parabola() -> BSplineCurve2 {
        BSplineCurve2 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0],
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(1.0, 2.0),
                DVec2::new(2.0, 0.0),
                DVec2::new(3.0, -2.0),
                DVec2::new(4.0, 0.0),
            ],
            weights: vec![1.0; 5],
        }
    }

    #[test]
    fn non_rational_d2_constant_per_segment() {
        let bs = two_segment_parabola();
        // Interior of segment 1.
        for &t in &[0.0_f64, 0.1, 0.5, 0.9] {
            let r2 = eval_d2(&bs, t);
            assert!(
                (r2.d2 - DVec2::new(0.0, -8.0)).length() < 1.0e-12,
                "D2({t}) = {:?}, want (0, -8)",
                r2.d2
            );
            let r1 = eval_d1(&bs, t);
            assert!((r1.d1 - DVec2::new(2.0, 4.0 - 8.0 * t)).length() < 1.0e-12);
            let r3 = eval_d3(&bs, t);
            assert!(r3.d3.length() < 1.0e-12, "degree-2 D3 must be zero");
        }
        // Interior of segment 2.
        for &t in &[1.1, 1.5, 1.9, 2.0] {
            let r2 = eval_d2(&bs, t);
            assert!(
                (r2.d2 - DVec2::new(0.0, 8.0)).length() < 1.0e-12,
                "D2({t}) = {:?}, want (0, +8)",
                r2.d2
            );
        }
    }

    #[test]
    fn c1_continuity_and_right_span_at_the_knot() {
        let bs = two_segment_parabola();
        // C1: D1 matches from both sides.
        let r_left = eval_d1(&bs, 1.0 - 1.0e-9);
        let r_at = eval_d1(&bs, 1.0);
        let r_right = eval_d1(&bs, 1.0 + 1.0e-9);
        assert!((r_at.d1 - r_left.d1).length() < 1.0e-8);
        assert!((r_at.d1 - r_right.d1).length() < 1.0e-8);
        // OCCT semantics: U = knot locates the RIGHT span — D2(1) is the
        // second segment's constant (the finite-difference default would
        // return the jump average).
        let r2 = eval_d2(&bs, 1.0);
        assert!(
            (r2.d2 - DVec2::new(0.0, 8.0)).length() < 1.0e-12,
            "D2 at the knot takes the right span: {:?}",
            r2.d2
        );
    }

    /// The rational quarter circle (degree 2, poles (1,0), (1,1), (0,1),
    /// weights (1, sqrt(2)/2, 1)) and the exact closed form of the rational
    /// quadratic parameterization, evaluated by the quotient rule over the
    /// polynomial (X, Y, W) with
    ///   W = 1 - b*t + b*t^2 (b = 2 - sqrt(2)),
    ///   X = W - b/2 * 2*t*(1-t) ... computed directly from the homogeneous
    /// de Boor coefficients.
    fn quarter_circle() -> BSplineCurve2 {
        let s2 = std::f64::consts::SQRT_2;
        BSplineCurve2 {
            degree: 2,
            knots: vec![
                0.0,
                0.0,
                0.0,
                std::f64::consts::FRAC_PI_2,
                std::f64::consts::FRAC_PI_2,
                std::f64::consts::FRAC_PI_2,
            ],
            control_points: vec![
                DVec2::new(1.0, 0.0),
                DVec2::new(1.0, 1.0),
                DVec2::new(0.0, 1.0),
            ],
            weights: vec![1.0, s2 / 2.0, 1.0],
        }
    }

    /// The closed-form (point, D1, D2) of the rational quadratic
    /// parameterization at parameter t, from the quotient rule over the
    /// Bernstein numerator (X, Y) and denominator W (polynomial evaluation
    /// in the power basis — machine-precision accurate).
    fn circle_reference(t: f64) -> (DVec2, DVec2, DVec2) {
        let s2 = std::f64::consts::SQRT_2;
        let b = 2.0 - s2; // W(t) = 1 - b t + b t^2.
        let c = s2 / 2.0; // the middle homogeneous coordinate pair.
        // Power-basis coefficients:
        // W  = 1 - b t + b t^2
        // X  = 1 + (c*2 - 2) t + (2 - 2c) t^2 ... expand:
        //   X = 1*(1-t)^2 + 2c*t(1-t)*1 + 0*t^2 (numerator x)
        //     = 1 + (2c - 2) t + (1 - 2c) t^2
        //   Y = 2c*t(1-t)*1 + t^2
        //     = 2c t + (1 - 2c) t^2
        let w0 = 1.0;
        let w1 = -b;
        let w2 = b;
        let x0 = 1.0;
        let x1 = 2.0 * c - 2.0;
        let x2 = 1.0 - 2.0 * c;
        let y0 = 0.0;
        let y1 = 2.0 * c;
        let y2 = 1.0 - 2.0 * c;
        // Evaluate the polynomials and their derivatives (Horner).
        let (x, xd) = (x0 + t * (x1 + t * x2), x1 + 2.0 * t * x2);
        let (y, yd) = (y0 + t * (y1 + t * y2), y1 + 2.0 * t * y2);
        let (w, wd) = (w0 + t * (w1 + t * w2), w1 + 2.0 * t * w2);
        // f = N/W: f' = (N'W - N W') / W^2; f'' = (N''W - N W'')/W^2 - 2 W' f' / W.
        let wdd = 2.0 * w2;
        let xdd = 2.0 * x2;
        let ydd = 2.0 * y2;
        let p = DVec2::new(x / w, y / w);
        let d1 = DVec2::new((xd * w - x * wd) / (w * w), (yd * w - y * wd) / (w * w));
        let d2 = DVec2::new(
            (xdd * w - x * wdd) / (w * w) - 2.0 * wd * d1.x / w,
            (ydd * w - y * wdd) / (w * w) - 2.0 * wd * d1.y / w,
        );
        (p, d1, d2)
    }

    #[test]
    fn rational_circle_d1_d2_matches_closed_form() {
        let bs = quarter_circle();
        for k in 0..=8 {
            let t = std::f64::consts::FRAC_PI_2 * (k as f64) / 8.0;
            let (p_ref, d1_ref, d2_ref) = circle_reference(t);
            let r1 = eval_d1(&bs, t);
            assert!((r1.point - p_ref).length() < 1.0e-12, "P({t})");
            assert!(
                (r1.d1 - d1_ref).length() < 1.0e-10,
                "D1({t}): got {:?}, want {:?}",
                r1.d1,
                d1_ref
            );
            let r2 = eval_d2(&bs, t);
            assert!(
                (r2.d2 - d2_ref).length() < 1.0e-9,
                "D2({t}): got {:?}, want {:?}",
                r2.d2,
                d2_ref
            );
        }
    }

    /// Endpoint semantics — OCCT LocateParameter clamps U into the first /
    /// last span, so D1/D2/D3 at U = First/Last evaluate the clamped span
    /// boundary exactly.  On the circle the endpoint derivatives are the
    /// exact axis-aligned vectors.
    #[test]
    fn endpoint_semantics() {
        let bs = quarter_circle();
        let half_pi = std::f64::consts::FRAC_PI_2;
        // At t = 0: P = (1, 0); D1 = (0, 2c)/W(0) = (0, sqrt(2)); D2 from
        // the reference.
        let (_, _, d2_ref0) = circle_reference(0.0);
        let r1 = eval_d1(&bs, 0.0);
        assert!((r1.point - DVec2::new(1.0, 0.0)).length() < 1.0e-12);
        assert!(r1.d1.x.abs() < 1.0e-12, "D1(0).x = {}", r1.d1.x);
        assert!(r1.d1.y > 0.0);
        let r2 = eval_d2(&bs, 0.0);
        assert!(
            (r2.d2 - d2_ref0).length() < 1.0e-9,
            "D2(0): got {:?}, want {:?} (the finite-difference default would \
             blow up on the clamped end)",
            r2.d2,
            d2_ref0
        );
        let r3 = eval_d3(&bs, 0.0);
        assert!(r3.d3.length() < 1.0e-7, "D3(0) = {:?}", r3.d3);
        // At t = half_pi: P = (0, 1), D1 = (-sqrt(2), 0).
        let r1e = eval_d1(&bs, half_pi);
        assert!((r1e.point - DVec2::new(0.0, 1.0)).length() < 1.0e-12);
        assert!(r1e.d1.y.abs() < 1.0e-12 && r1e.d1.x < 0.0);
        // Slightly outside the domain clamps into the end spans (OCCT
        // LocateParameter clamps NewU into [Knots(first), Knots(last)]).
        let r_out = eval_d1(&bs, half_pi + 0.01);
        assert!((r_out.d1 - r1e.d1).length() < 1.0e-9, "clamped end");
        // The non-rational curve D2 at its ends is the segment constant.
        let bs2 = two_segment_parabola();
        let r2a = eval_d2(&bs2, 0.0);
        assert!((r2a.d2 - DVec2::new(0.0, -8.0)).length() < 1.0e-12);
        let r2b = eval_d2(&bs2, 2.0);
        assert!((r2b.d2 - DVec2::new(0.0, 8.0)).length() < 1.0e-12);
    }
}
