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
//! 2. `BSplineCurve2` carries the `is_periodic` flag (the carrier mirror of
//!    OCCT `myPeriodic`); every OCCT `myPeriodic` site below reads it.  The
//!    carrier stores the PLAIN (non-wrapped) flat knot expansion, while the
//!    OCCT object keeps the wrapped `myFlatKnots` sequence — where that
//!    difference changes an index (`PeriodicNormalization` endpoints) the
//!    carrier equivalent is annotated at the statement.
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
/// periodic curves locate over the full knot range (first = Knots.Lower(),
/// last = Knots.Upper(), L332-337), non-periodic curves over the
/// FirstUKnotIndex/LastUKnotIndex pair (L338-341); the inner locate runs
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
    // OCCT L331-345.
    let first;
    let last;
    if is_periodic {
        first = 1i32;
        last = knots.len() as i32;
    } else {
        first = first_uknot_index_mults(degree as usize, mults);
        last = last_uknot_index_mults(degree as usize, mults);
    }
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
    // The rcad BSplineCurve2 carries an EMPTY weights vec for the
    // non-rational encoding — that is the OCCT null Weights handle
    // (`Weights != nullptr`, pxx L805): an empty slice must read as None,
    // or the span-local IsRational probe divides by a zero length.
    let weights = weights.filter(|the_w| !the_w.is_empty());

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
/// flat-knot period.  `is_periodic` is the carrier mirror of OCCT
/// `myPeriodic`; the carrier stores the PLAIN (non-wrapped) flat knot
/// expansion (module architecture note 2), whose first/last entries are
/// exactly the wrapped-sequence endpoints the OCCT statement reads:
/// myFlatKnots.Value(myDeg + 1) == Knots(1) == flat[0] and
/// myFlatKnots.Value(Upper - myDeg) == Knots(NbKnots) == flat[Upper]
/// (the end multiplicities of a periodic layout equal the degree).
fn periodic_normalization(bs: &BSplineCurve2, is_periodic: bool, parameter: &mut f64) {
    if is_periodic {
        let flat = &bs.knots;
        let k_upper = flat.len() as i32;
        // OCCT L1332: Period = myFlatKnots.Value(Upper - myDeg)
        // - myFlatKnots.Value(myDeg + 1).
        let period = at(flat, k_upper) - at(flat, 1);
        // OCCT L1333-1339.
        while *parameter > at(flat, k_upper) {
            *parameter -= period;
        }
        // OCCT L1340-1345.
        while *parameter < at(flat, 1) {
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
    periodic_normalization(bs, bs.is_periodic, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L193: LocateParameter(myDeg, myKnots, &myMults, U, myPeriodic,
    // aSpanIndex, aNewU).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, bs.is_periodic, &mut a_span_index, &mut a_new_u);
    // OCCT L195-198.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L200: BSplCLib::D0(...).
    bsplclib_d0(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        bs.is_periodic,
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
    periodic_normalization(bs, bs.is_periodic, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L208: LocateParameter(...).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, bs.is_periodic, &mut a_span_index, &mut a_new_u);
    // OCCT L210-213.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L215-224: BSplCLib::D1(...).
    bsplclib_d1(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        bs.is_periodic,
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
    periodic_normalization(bs, bs.is_periodic, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L241: LocateParameter(...).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, bs.is_periodic, &mut a_span_index, &mut a_new_u);
    // OCCT L243-246.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L248-257: BSplCLib::D2(...).
    bsplclib_d2(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        bs.is_periodic,
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
    periodic_normalization(bs, bs.is_periodic, &mut a_new_u);
    let (knots, mults) = knots_mults_of(&bs.knots);
    // OCCT L275: LocateParameter(...).
    locate_parameter_bspline(bs.degree as i32, &knots, &mults, u, bs.is_periodic, &mut a_span_index, &mut a_new_u);
    // OCCT L277-280.
    if a_new_u < at(&knots, a_span_index) {
        a_span_index -= 1;
    }
    // OCCT L282-291: BSplCLib::D3(...).
    bsplclib_d3(
        a_new_u,
        a_span_index,
        bs.degree as i32,
        bs.is_periodic,
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
    // The tests drive the Curve2dEval trait methods so they exercise the
    // override dispatch (not the module functions directly).
    use crate::geom::Curve2dEval;

/// Trait-dispatched D1 (the eval.rs override of BSplineCurve2).
fn d1_of(bs: &BSplineCurve2, t: f64) -> (DVec2, DVec2) {
    (bs.point_at(t), bs.derivative_at(t))
}

/// Trait-dispatched D2.
fn d2_of(bs: &BSplineCurve2, t: f64) -> DVec2 {
    bs.derivative2_at(t)
}

/// Trait-dispatched D3.
fn d3_of(bs: &BSplineCurve2, t: f64) -> DVec2 {
    bs.derivative3_at(t)
}

fn make_two_segment() -> BSplineCurve2 {
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
        is_periodic: false,
    }
}

fn make_quarter_circle() -> BSplineCurve2 {
    let s2 = std::f64::consts::SQRT_2;
    BSplineCurve2 {
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

/// Closed-form (P, D1, D2, D3) of the rational quarter circle over the
/// clamped knots [0,0,0,1,1,1] (the rational quadratic Bezier in the
/// plain [0,1] parameter), by the quotient rule over the degree-2 power
/// forms  X = 1 + (s2-2)t + (1-s2)t^2,  Y = s2 t + (1-s2)t^2,
/// W = 1 - (2-s2) t + (2-s2) t^2.  The third numerator derivative is
/// zero, so the third quotient derivative is (-3 f2 W1 - 3 f1 W2) / W.
fn circle_ref(t: f64) -> (DVec2, DVec2, DVec2, DVec2) {
    let s2 = std::f64::consts::SQRT_2;
    let b = 2.0 - s2;
    let c = s2 / 2.0;
    let (x0, x1, x2) = (1.0, 2.0 * c - 2.0, 1.0 - 2.0 * c);
    let (y0, y1, y2) = (0.0, 2.0 * c, 1.0 - 2.0 * c);
    let (w0, w1, w2) = (1.0, -b, b);
    let (x, xd) = (x0 + t * (x1 + t * x2), x1 + 2.0 * t * x2);
    let (y, yd) = (y0 + t * (y1 + t * y2), y1 + 2.0 * t * y2);
    let (w, wd) = (w0 + t * (w1 + t * w2), w1 + 2.0 * t * w2);
    let wdd = 2.0 * w2;
    let xdd = 2.0 * x2;
    let ydd = 2.0 * y2;
    let p = DVec2::new(x / w, y / w);
    let d1 = DVec2::new((xd * w - x * wd) / (w * w), (yd * w - y * wd) / (w * w));
    let d2 = DVec2::new(
        (xdd * w - x * wdd) / (w * w) - 2.0 * wd * d1.x / w,
        (ydd * w - y * wdd) / (w * w) - 2.0 * wd * d1.y / w,
    );
    let d3 = DVec2::new(
        (-3.0 * d2.x * wd - 3.0 * d1.x * wdd) / w,
        (-3.0 * d2.y * wd - 3.0 * d1.y * wdd) / w,
    );
    (p, d1, d2, d3)
}

/// (a) degree-2 non-rational: D2 is constant per segment; the constants
/// and the linear D1 come from the derivative-pole recursion
///   D_i = 2 (P_{i+1} - P_i) / (T_{i+3} - T_{i+1}):
/// D_1 = (2,4), D_2 = (1,-2), D_3 = (2,-4)  =>
/// on (0,1): D1(u) = (2-u, 4-6u), D2 = (-1,-6);
/// on (1,2): D1(u) = (1+v, -2-2v), v = u-1, D2 = (1,-2).
#[test]
fn non_rational_d2_constant_per_segment() {
    let bs = make_two_segment();
    for &t in &[0.0_f64, 0.1, 0.3, 0.5, 0.9] {
        let (p, d1) = d1_of(&bs, t);
        let p_ref = bs.point_at(t);
        assert!((p - p_ref).length() < 1.0e-12);
        assert!(
            (d1 - DVec2::new(2.0 - t, 4.0 - 6.0 * t)).length() < 1.0e-12,
            "D1({t}) = {:?}",
            d1
        );
        let d2 = d2_of(&bs, t);
        assert!(
            (d2 - DVec2::new(-1.0, -6.0)).length() < 1.0e-12,
            "D2({t}) = {:?}, want (-1, -6)",
            d2
        );
        let d3 = d3_of(&bs, t);
        assert!(d3.length() < 1.0e-12, "degree-2 D3 must be zero");
    }
    for &v in &[0.1_f64, 0.5, 0.9, 1.0] {
        let t = 1.0 + v;
        let (_, d1) = d1_of(&bs, t);
        assert!(
            (d1 - DVec2::new(1.0 + v, -2.0 - 2.0 * v)).length() < 1.0e-12,
            "D1({t}) = {:?}",
            d1
        );
        let d2 = d2_of(&bs, t);
        assert!(
            (d2 - DVec2::new(1.0, -2.0)).length() < 1.0e-12,
            "D2({t}) = {:?}, want (1, -2)",
            d2
        );
    }
}

/// C1 continuity across the knot and the OCCT right-span semantics for
/// D2 at the knot (the stencil default would average the jump).
#[test]
fn c1_continuity_and_right_span_at_the_knot() {
    let bs = make_two_segment();
    let d_left = bs.derivative_at(1.0 - 1.0e-9);
    let d_at = bs.derivative_at(1.0);
    let d_right = bs.derivative_at(1.0 + 1.0e-9);
    assert!((d_at - d_left).length() < 1.0e-8);
    assert!((d_at - d_right).length() < 1.0e-8);
    assert!((d_at - DVec2::new(1.0, -2.0)).length() < 1.0e-12);
    let d2 = d2_of(&bs, 1.0);
    assert!(
        (d2 - DVec2::new(1.0, -2.0)).length() < 1.0e-12,
        "D2 at the knot takes the right span: {:?}",
        d2
    );
}

/// (b) rational quarter circle: D1/D2 match the closed form at sampled
/// parameters to 1e-10 / 1e-9.
#[test]
fn rational_circle_d1_d2_matches_closed_form() {
    let bs = make_quarter_circle();
    for k in 0..=8 {
        let t = k as f64 / 8.0;
        let (p_ref, d1_ref, d2_ref, _) = circle_ref(t);
        let (p, d1) = d1_of(&bs, t);
        assert!((p - p_ref).length() < 1.0e-12, "P({t})");
        assert!(
            (d1 - d1_ref).length() < 1.0e-10,
            "D1({t}): got {:?}, want {:?}",
            d1,
            d1_ref
        );
        let d2 = d2_of(&bs, t);
        assert!(
            (d2 - d2_ref).length() < 1.0e-9,
            "D2({t}): got {:?}, want {:?}",
            d2,
            d2_ref
        );
    }
}

/// (c) endpoint behavior: exact one-sided boundary derivatives at
/// U = First / Last (the point-clamped stencil default returns a
/// half-length D1 and O(1/h) noise on D2 there); beyond the domain the
/// D-kernels clamp the located SPAN INDEX, not the parameter
/// (LocateParameter assigns NewU = U), so the evaluation follows the
/// polynomial continuation.
#[test]
fn endpoint_semantics() {
    let bs = make_quarter_circle();
    let (_, d1_ref0, d2_ref0, d3_ref0) = circle_ref(0.0);
    let (p0, d1_0) = d1_of(&bs, 0.0);
    assert!((p0 - DVec2::new(1.0, 0.0)).length() < 1.0e-12);
    assert!((d1_0 - d1_ref0).length() < 1.0e-10, "D1(0) = {:?}", d1_0);
    let d2_0 = d2_of(&bs, 0.0);
    assert!((d2_0 - d2_ref0).length() < 1.0e-9, "D2(0) = {:?}", d2_0);
    let d3_0 = d3_of(&bs, 0.0);
    assert!((d3_0 - d3_ref0).length() < 1.0e-8, "D3(0) = {:?}", d3_0);

    let (_, d1_ref1, _, _) = circle_ref(1.0);
    let (p1, d1_1) = d1_of(&bs, 1.0);
    assert!((p1 - DVec2::new(0.0, 1.0)).length() < 1.0e-12);
    assert!((d1_1 - d1_ref1).length() < 1.0e-10, "D1(1) = {:?}", d1_1);

    let d1_out = bs.derivative_at(1.0 + 0.01);
    let (_, d1_cont, _, _) = circle_ref(1.01);
    assert!((d1_out - d1_cont).length() < 1.0e-9, "continuation");

    // The non-rational curve D2 at its ends are the segment constants.
    let bs2 = make_two_segment();
    assert!((d2_of(&bs2, 0.0) - DVec2::new(-1.0, -6.0)).length() < 1.0e-12);
    assert!((d2_of(&bs2, 2.0) - DVec2::new(1.0, -2.0)).length() < 1.0e-12);
}

/// Periodic DN evaluation (the `myPeriodic` arms now reachable through
/// `BSplineCurve2::is_periodic`): a degree-1 CLOSED bspline over knots
/// [0, 1, 2], mults [1, 1, 1] (end multiplicities = Degree), poles
/// (0,0), (1,1) — the closed polyline P(u) = (u, u) on [0, 1] and
/// P(u) = (2-u, 2-u) on [1, 2], period 2.  The wrapped-span samples:
///   P(2.5) = P(0.5) = (0.5, 0.5)   (PeriodicNormalization -1 x period)
///   P(-0.5) = P(1.5) = (0.5, 0.5)  (PeriodicNormalization +1 x period)
///   P(0) = P(2) = (0, 0)           (the seam), D1 wraps the same way.
/// Without the flag the evaluation follows the clamped continuation
/// (span index clamped to the last span) and none of these hold.
#[test]
fn periodic_wrapped_span_evaluation() {
    let bs = BSplineCurve2 {
        degree: 1,
        knots: vec![0.0, 1.0, 2.0],
        control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(1.0, 1.0)],
        weights: vec![1.0, 1.0],
        is_periodic: true,
    };
    // Wrapped samples (hand-computed closed form above).
    let (p, d1) = d1_of(&bs, 2.5);
    assert!(p.distance(DVec2::new(0.5, 0.5)) < 1.0e-12, "P(2.5) = {:?}", p);
    assert!(d1.distance(DVec2::new(1.0, 1.0)) < 1.0e-12, "D1(2.5) = {:?}", d1);
    let (p, d1) = d1_of(&bs, -0.5);
    assert!(p.distance(DVec2::new(0.5, 0.5)) < 1.0e-12, "P(-0.5) = {:?}", p);
    assert!(d1.distance(DVec2::new(-1.0, -1.0)) < 1.0e-12, "D1(-0.5) = {:?}", d1);
    // The seam: first and last knot evaluate identically.
    let p_first = bs.point_at(0.0);
    let p_last = bs.point_at(2.0);
    assert!(p_first.distance(DVec2::new(0.0, 0.0)) < 1.0e-12);
    assert!(p_last.distance(p_first) < 1.0e-12, "seam wrap");
    // D2/D3 stay zero (degree 1) across the wrap.
    assert!(d2_of(&bs, 2.5).length() < 1.0e-12);
    assert!(d3_of(&bs, 2.5).length() < 1.0e-12);
}
}
