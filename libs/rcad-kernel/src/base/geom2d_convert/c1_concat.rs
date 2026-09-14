//! OCCT Geom2dConvert C1 concatenation family (TKGeomBase/Geom2dConvert/
//! Geom2dConvert.cxx L430-1542): the `ConcatC1` statics, the
//! `MultNumandDenom`/`Continuity`/... helpers and `C0BSplineToC1BSplineCurve`
//! consumers, plus the `Geom2d_BSplineCurve` member functions this batch
//! needs that are not hosted in `bspline_curve.rs` (whose edit surface is
//! fixed to D0/D1/SetPeriodic by this batch).
//!
//! Architecture note (value semantics): OCCT mutates the curve handle in
//! place (`SetPole`, `SetWeight`, `RemoveKnot`); rcad's `Geom2dBSplineCurve`
//! keeps its pole/knot storage private, so the mutation methods below are
//! expressed as functions taking the curve by reference and returning the
//! updated curve (read-out through the public accessors, mutate, rebuild
//! through the constructors — the constructors run the same
//! `updateKnots()`/rationality fixups the OCCT methods end with).

use glam::DVec2;

use crate::base::geom2d_convert::bspline_curve::Geom2dBSplineCurve;
use crate::math::gp::GP_RESOLUTION;
use crate::math::bspl_lib::{at, knot_sequence, knot_sequence_length, remove_knot, BSPLIB_MAX_DEGREE};

// ===========================================================================
// Geom2d_BSplineCurve member functions consumed by the concatenation family
// ===========================================================================

/// OCCT Geom2d_BSplineCurve::Knots() accessor
/// (Geom2d_BSplineCurve_1.cxx L381-386) — the distinct-knot array.
pub fn curve_knots(bs: &Geom2dBSplineCurve) -> Vec<f64> {
    (1..=bs.nb_knots()).map(|i| bs.knot(i)).collect()
}

/// OCCT Geom2d_BSplineCurve::Multiplicities() accessor
/// (Geom2d_BSplineCurve_1.cxx L591-596) — the distinct multiplicities.
pub fn curve_mults(bs: &Geom2dBSplineCurve) -> Vec<i32> {
    (1..=bs.nb_knots()).map(|i| bs.multiplicity(i)).collect()
}

/// OCCT Geom2d_BSplineCurve::WeightsArray() — the per-pole weights (unit
/// weights for a non-rational curve, matching the OCCT storage kept by
/// updateKnots / UnitWeights).
pub fn curve_weights(bs: &Geom2dBSplineCurve) -> Vec<f64> {
    (1..=bs.nb_poles_curve()).map(|i| bs.weight(i)).collect()
}

/// OCCT Geom2d_BSplineCurve::KnotSequence()
/// (Geom2d_BSplineCurve_1.cxx L398-403) — the flat knot sequence myFlatKnots,
/// rebuilt through BSplCLib::KnotSequence (the identical kernel updateKnots
/// runs, deterministic in knots/mults/degree/periodic).
pub fn curve_flat_knots(bs: &Geom2dBSplineCurve) -> Vec<f64> {
    let knots = curve_knots(bs);
    let mults = curve_mults(bs);
    let len = knot_sequence_length(&mults, bs.degree(), bs.is_periodic());
    let mut flat = vec![0.0f64; len];
    knot_sequence(&knots, &mults, bs.degree(), bs.is_periodic(), &mut flat);
    flat
}

/// OCCT Geom2d_BSplineCurve::SetPole(Index, P) (Geom2d_BSplineCurve.cxx
/// L1134-1142) in value semantics: returns the curve with pole `index`
/// replaced by `p`.
pub fn set_pole(bs: &Geom2dBSplineCurve, index: i32, p: DVec2) -> Geom2dBSplineCurve {
    if index < 1 || index > bs.nb_poles_curve() {
        panic!("Standard_OutOfRange: BSpline curve: SetPole: index and #pole mismatch");
    }
    let mut poles: Vec<DVec2> = (1..=bs.nb_poles_curve()).map(|i| bs.pole(i)).collect();
    poles[(index - 1) as usize] = p;
    set_pole_value(bs, poles)
}

/// OCCT Geom2d_BSplineCurve::SetWeight(Index, W) (Geom2d_BSplineCurve.cxx
/// L1154-1191) in value semantics: returns the curve with weight `index` set
/// to `w` (with the OCCT rationality fixup: a uniform weight set collapses
/// back to non-rational unit weights).
pub fn set_weight(bs: &Geom2dBSplineCurve, index: i32, w: f64) -> Geom2dBSplineCurve {
    if index < 1 || index > bs.nb_poles_curve() {
        panic!("Standard_OutOfRange: BSpline curve: SetWeight: Index and #pole mismatch");
    }
    if w <= GP_RESOLUTION {
        panic!("Standard_ConstructionError: BSpline curve: SetWeight: Weight too small");
    }
    // rat = IsRational() || (std::abs(W - 1.) > gp::Resolution()); if (rat)
    // { myWeights.SetValue(Index, W); ... }
    let rat = bs.is_rational() || ((w - 1.0).abs() > GP_RESOLUTION);
    if !rat {
        return bs.clone();
    }
    let mut weights = curve_weights(bs);
    weights[(index - 1) as usize] = w;
    // Geom2d_BSplineCurve(Poles, Weights, Knots, Mults, Degree, Periodic)
    // applies Rational(myWeights) and resets to unit weights when uniform —
    // the OCCT SetWeight tail (L1178-1188).
    let poles: Vec<DVec2> = (1..=bs.nb_poles_curve()).map(|i| bs.pole(i)).collect();
    Geom2dBSplineCurve::new_rational(
        poles,
        weights,
        curve_knots(bs),
        curve_mults(bs),
        bs.degree(),
        bs.is_periodic(),
    )
}

/// OCCT Geom2d_BSplineCurve::RemoveKnot(Index, M, Tolerance)
/// (Geom2d_BSplineCurve.cxx L410-476) in value semantics: `None` is the
/// OCCT `false` return (knot not removed).
pub fn remove_knot_curve(
    bs: &Geom2dBSplineCurve,
    index: i32,
    m: i32,
    tolerance: f64,
) -> Option<Geom2dBSplineCurve> {
    if m < 0 {
        return Some(bs.clone());
    }
    let i1 = bs.first_uknot_index();
    let i2 = bs.last_uknot_index();
    if index < i1 || index > i2 {
        panic!("Standard_OutOfRange: BSpline curve: RemoveKnot: index out of range");
    }
    let step = bs.multiplicity(index) - m;
    if step <= 0 {
        return Some(bs.clone());
    }

    let nbpoles = (bs.nb_poles_curve() - step) as usize;
    let nbknots = (bs.nb_knots() - (if m == 0 { 1 } else { 0 })) as usize;

    let rational = bs.is_rational();
    let weights = curve_weights(bs);
    // PLib::SetPoles: homogeneous dim 3 when rational, dim 2 otherwise.
    let dim = if rational { 3usize } else { 2usize };
    let poles: Vec<DVec2> = (1..=bs.nb_poles_curve()).map(|i| bs.pole(i)).collect();
    let mut flat = Vec::with_capacity(poles.len() * dim);
    for (k, p) in poles.iter().enumerate() {
        match rational {
            true => {
                let w = weights[k];
                flat.push(p.x * w);
                flat.push(p.y * w);
                flat.push(w);
            }
            false => {
                flat.push(p.x);
                flat.push(p.y);
            }
        }
    }
    let mut npoles_flat = vec![0.0f64; nbpoles * dim];
    let mut nknots = vec![0.0f64; nbknots];
    let mut nmults = vec![0i32; nbknots];

    // BSplCLib::RemoveKnot(Index, M, myDeg, myPeriodic, myPoles, Weights(),
    //                      myKnots, myMults, npoles, nweights, nknots, nmults,
    //                      Tolerance)
    let ok = remove_knot(
        index as usize,
        m,
        bs.degree(),
        bs.is_periodic(),
        dim,
        &flat,
        &curve_knots(bs),
        &curve_mults(bs),
        &mut npoles_flat,
        &mut nknots,
        &mut nmults,
        tolerance,
    );
    if !ok {
        return None;
    }

    // PLib::GetPoles / dehomogenize (OCCT L461-471).
    let mut npoles: Vec<DVec2> = Vec::with_capacity(nbpoles);
    let mut nweights: Vec<f64> = Vec::new();
    if rational {
        nweights.resize(nbpoles, 0.0);
        for (k, c) in npoles_flat.chunks_exact(dim).enumerate() {
            npoles.push(DVec2::new(c[0] / c[2], c[1] / c[2]));
            nweights[k] = c[2];
        }
        Some(Geom2dBSplineCurve::new_rational(
            npoles, nweights, nknots, nmults, bs.degree(), bs.is_periodic(),
        ))
    } else {
        for c in npoles_flat.chunks_exact(2) {
            npoles.push(DVec2::new(c[0], c[1]));
        }
        Some(Geom2dBSplineCurve::new(
            npoles, nknots, nmults, bs.degree(), bs.is_periodic(),
        ))
    }
}

/// OCCT BSplCLib::Resolution — the `case 2:` arm (BSplCLib.cxx L4316-4445)
/// used by Geom2d_BSplineCurve::Resolution (Geom2d_BSplineCurve_1.cxx L764-806
/// with the myMaxDerivInv cache elided — the kernel struct carries no cache).
/// Returns `UTolerance = ToleranceUV * myMaxDerivInv`.
pub fn resolution_curve(bs: &Geom2dBSplineCurve, tolerance_uv: f64) -> f64 {
    // The periodic arm (L768-792) unperiodizes the poles first.
    let (flat, weights, num_poles, flat_knots) = if bs.is_periodic() {
        // BSplCLib::PrepareUnperiodize sizes; the poles repeat modulo the
        // non-periodic count (L774-784).
        let mults = curve_mults(bs);
        let degree = bs.degree();
        let knots = curve_knots(bs);
        let mut nb_knots = 0i32;
        let mut nb_poles_count = 0i32;
        crate::math::bspl_lib::prepare_unperiodize(degree, &mults, &mut nb_knots, &mut nb_poles_count);
        let nb = nb_poles_count as usize;
        let np = bs.nb_poles_curve();
        let raw_weights = curve_weights(bs);
        let rational = bs.is_rational();
        let dim = if rational { 3usize } else { 2usize };
        let mut new_flat = Vec::with_capacity(nb * dim);
        let mut new_weights = vec![0.0f64; nb];
        for ii in 1..=nb as i32 {
            let p = bs.pole(((ii - 1) % np) + 1);
            if rational {
                let w = raw_weights[(((ii - 1) % np) + 1 - 1) as usize];
                new_weights[(ii - 1) as usize] = w;
                new_flat.push(p.x * w);
                new_flat.push(p.y * w);
                new_flat.push(w);
            } else {
                new_flat.push(p.x);
                new_flat.push(p.y);
            }
        }
        let len = knot_sequence_length(&mults, degree, true);
        let mut fk = vec![0.0f64; len];
        knot_sequence(&knots, &mults, degree, true, &mut fk);
        (new_flat, if rational { Some(new_weights) } else { None }, nb, fk)
    } else {
        let rational = bs.is_rational();
        let dim = if rational { 3usize } else { 2usize };
        let weights = curve_weights(bs);
        let poles: Vec<DVec2> = (1..=bs.nb_poles_curve()).map(|i| bs.pole(i)).collect();
        let mut flat = Vec::with_capacity(poles.len() * dim);
        for (k, p) in poles.iter().enumerate() {
            match rational {
                true => {
                    let w = weights[k];
                    flat.push(p.x * w);
                    flat.push(p.y * w);
                    flat.push(w);
                }
                false => {
                    flat.push(p.x);
                    flat.push(p.y);
                }
            }
        }
        (
            flat,
            if rational { Some(weights) } else { None },
            poles.len(),
            curve_flat_knots(bs),
        )
    };

    // UTolerance = ToleranceUV * myMaxDerivInv with Tolerance3D = 1.
    let max_deriv_inv = resolution_2d(&flat, weights.as_deref(), num_poles, &flat_knots, bs.degree());
    tolerance_uv * max_deriv_inv
}

/// OCCT BSplCLib::Resolution `case 2:` body (BSplCLib.cxx L4340-4445) —
/// returns myMaxDerivInv for a 2-component pole array (`weights == None` is
/// BSplCLib::NoWeights()).
fn resolution_2d(
    pa: &[f64],
    weights: Option<&[f64]>,
    num_poles_in: usize,
    flat_knots: &[f64],
    degree: usize,
) -> f64 {
    let deg1 = degree as i32 + 1;
    let deg2 = ((degree as i32) << 1) + 1;
    let mut max_derivative = 0.0f64;
    let num_poles = (flat_knots.len() as i32 - deg1) as usize;
    let _ = num_poles_in;
    let fk = |i: i32| at(flat_knots, i);
    match weights {
        Some(wg) => {
            // min_weights over WG[0..NumPoles-1] (L4343-4353).
            let mut min_weights = wg[0];
            for ii in 1..num_poles {
                let w = wg[ii];
                if w < min_weights {
                    min_weights = w;
                }
            }
            for ii in 1..num_poles as i32 {
                let ii_index = ii % num_poles as i32;
                let mut ii_in_dim = ii_index << 1;
                let ii_minus = (ii - 1) % num_poles as i32;
                let mut ii_mi_dim = ii_minus << 1;
                let pa_ii_in_dim_0 = pa[ii_in_dim as usize];
                ii_in_dim += 1;
                let pa_ii_in_dim_1 = pa[ii_in_dim as usize];
                let pa_ii_mi_dim_0 = pa[ii_mi_dim as usize];
                ii_mi_dim += 1;
                let pa_ii_mi_dim_1 = pa[ii_mi_dim as usize];
                let wg_ii_index = wg[ii_index as usize];
                let wg_ii_minus = wg[ii_minus as usize];
                let mut inverse = fk(ii + degree as i32) - fk(ii);
                inverse = 1.0 / inverse;
                let mut lower = ii - deg1;
                if lower < 0 {
                    lower = 0;
                }
                let mut upper = deg2 + ii;
                if upper > num_poles as i32 {
                    upper = num_poles as i32;
                }

                let mut jj = lower;
                while jj < upper {
                    let mut jj_index = jj % num_poles as i32;
                    jj_index <<= 1;
                    let mut value = 0.0f64;
                    let mut factor = ((pa[jj_index as usize] - pa_ii_in_dim_0) * wg_ii_index)
                        - ((pa[jj_index as usize] - pa_ii_mi_dim_0) * wg_ii_minus);
                    if factor < 0.0 {
                        factor = -factor;
                    }
                    value += factor;
                    jj_index += 1;
                    factor = ((pa[jj_index as usize] - pa_ii_in_dim_1) * wg_ii_index)
                        - ((pa[jj_index as usize] - pa_ii_mi_dim_1) * wg_ii_minus);
                    if factor < 0.0 {
                        factor = -factor;
                    }
                    value += factor;
                    value *= inverse;
                    if max_derivative < value {
                        max_derivative = value;
                    }
                    jj += 1;
                }
            }
            max_derivative /= min_weights;
        }
        None => {
            for ii in 1..num_poles as i32 {
                let mut ii_index = ii % num_poles as i32;
                ii_index <<= 1;
                let ii_minus = ((ii - 1) % num_poles as i32) << 1;
                let mut inverse = fk(ii + degree as i32) - fk(ii);
                inverse = 1.0 / inverse;
                let mut value = 0.0f64;
                let mut factor = pa[ii_index as usize] - pa[ii_minus as usize];
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                ii_index += 1;
                let ii_minus1 = ii_minus + 1;
                factor = pa[ii_index as usize] - pa[ii_minus1 as usize];
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                value *= inverse;
                if max_derivative < value {
                    max_derivative = value;
                }
            }
        }
    }
    // OCCT tail: max_derivative *= Degree; UTol = Tol / max (or / RealSmall).
    max_derivative *= degree as f64;
    if max_derivative > f64::MIN_POSITIVE {
        1.0 / max_derivative
    } else {
        1.0 / f64::MIN_POSITIVE
    }
}

/// OCCT Geom2d_BSplineCurve::MaxDegree() (Geom2d_BSplineCurve.cxx L229-231).
pub fn max_degree() -> usize {
    BSPLIB_MAX_DEGREE
}

/// Rebuild a curve from new poles keeping knots/mults/weights/degree (the
/// OCCT methods mutate one member then run updateKnots(); the constructors
/// run the identical fixups).
pub(crate) fn set_pole_value(bs: &Geom2dBSplineCurve, poles: Vec<DVec2>) -> Geom2dBSplineCurve {
    if bs.is_rational() {
        let weights: Vec<f64> = (1..=bs.nb_poles_curve()).map(|i| bs.weight(i)).collect();
        Geom2dBSplineCurve::new_rational(
            poles,
            weights,
            curve_knots(bs),
            curve_mults(bs),
            bs.degree(),
            bs.is_periodic(),
        )
    } else {
        Geom2dBSplineCurve::new(
            poles,
            curve_knots(bs),
            curve_mults(bs),
            bs.degree(),
            bs.is_periodic(),
        )
    }
}

// ===========================================================================
// Geom2dConvert statics (Geom2dConvert.cxx L430-1542)
// ===========================================================================

use crate::base::convert::ConvertParameterisation;
use crate::base::extrema_ext_elc::epsilon_of;
use crate::base::geom2d_convert::comp_curve_to_bspline_2d::Geom2dConvertCompCurveToBSplineCurve;
use crate::core::precision::{ANGULAR, CONFUSION};
use crate::geom::{BSplineCurve2, Curve2d};
use crate::math::bspl_lib::{build_schoenberg_points, eval_flat, interpolate, BSplCLibEvaluatorFunction};
use crate::math::bspl_lib::{nb_poles, reparametrize};

/// OCCT GeomAbs_Shape — the enumerators in declaration order (C0, G1, C1,
/// G2, C2, C3, CN); the derived Ord carries the OCCT comparisons
/// (`< GeomAbs_C0`, `>= GeomAbs_G1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum GeomAbsConcatShape {
    C0,
    G1,
    C1,
    G2,
    C2,
    C3,
    CN,
}

impl From<crate::base::geom2d_convert::bspline_curve::SmoothShape> for GeomAbsConcatShape {
    fn from(s: crate::base::geom2d_convert::bspline_curve::SmoothShape) -> Self {
        use crate::base::geom2d_convert::bspline_curve::SmoothShape as S;
        match s {
            S::CN => GeomAbsConcatShape::CN,
            S::C0 => GeomAbsConcatShape::C0,
            S::C1 => GeomAbsConcatShape::C1,
            S::C2 => GeomAbsConcatShape::C2,
            S::C3 => GeomAbsConcatShape::C3,
        }
    }
}

/// OCCT static GeomAbsToInteger (Geom2dConvert.cxx L716-744).
pub fn geom_abs_to_integer(gcont: GeomAbsConcatShape) -> i32 {
    match gcont {
        GeomAbsConcatShape::C0 => 0,
        GeomAbsConcatShape::G1 => 1,
        GeomAbsConcatShape::C1 => 2,
        GeomAbsConcatShape::G2 => 3,
        GeomAbsConcatShape::C2 => 4,
        GeomAbsConcatShape::C3 => 5,
        GeomAbsConcatShape::CN => 6,
    }
}

/// OCCT static Continuity (Geom2dConvert.cxx L748-852) — continuity at the
/// junction of two curves.  Architecture note: OCCT receives Geom2d_Curve
/// handles, extracts Geom2d_TrimmedCurve basis curves (L771-780) and
/// downcasts to Geom2d_BSplineCurve (L781-810); ConcatC1 only ever passes
/// Geom2d_BSplineCurve handles, so the arguments are the BSpline curves
/// directly and the downcast branches are total.
#[allow(clippy::too_many_arguments, clippy::many_single_char_names)]
pub fn continuity(
    c1: &Geom2dBSplineCurve,
    c2: &Geom2dBSplineCurve,
    u1: f64,
    u2: f64,
    r1: bool,
    r2: bool,
    tl: f64,
    ta: f64,
) -> GeomAbsConcatShape {
    let mut cont = GeomAbsConcatShape::C0;
    let mut cont1 = geom_abs_to_integer(c1.continuity().into());
    let mut cont2 = geom_abs_to_integer(c2.continuity().into());

    // BSplineCurve branch (L781-795).
    let tolerance = resolution_curve(c1, tl);
    let mut index1 = 0i32;
    let mut index2 = 0i32;
    c1.locate_u(u1, tolerance, &mut index1, &mut index2, false);
    if index1 > 1 && index2 < c1.nb_knots() && index1 == index2 {
        cont1 = c1.degree() as i32 - c1.multiplicity(index1);
    } else {
        cont1 = 5;
    }

    // BSplineCurve branch (L796-810).
    let tolerance = resolution_curve(c2, tl);
    let mut index1 = 0i32;
    let mut index2 = 0i32;
    c2.locate_u(u2, tolerance, &mut index1, &mut index2, false);
    if index1 > 1 && index2 < c2.nb_knots() && index1 == index2 {
        cont2 = c2.degree() as i32 - c2.multiplicity(index1);
    } else {
        cont2 = 5;
    }

    let (point1, mut d1) = c1.eval_d1(u1);
    let (point2, mut d2) = c2.eval_d1(u2);
    if point1.distance_squared(point2) <= tl * tl {
        if cont1 != 0 && cont2 != 0 {
            if d1.length_squared() >= tl * tl && d2.length_squared() >= tl * tl {
                if r1 {
                    d1 = DVec2::new(-d1.x, -d1.y);
                }
                if r2 {
                    d2 = DVec2::new(-d2.x, -d2.y);
                }
                let value = d1.dot(d2);
                if (d1.length() <= (d2.length() + tl))
                    && (d1.length() >= (d2.length() - tl))
                    && (value / (d1.length() * d2.length()) >= 1.0 - ta * ta)
                {
                    cont = GeomAbsConcatShape::C1;
                } else {
                    let n1 = d1.normalize();
                    let n2 = d2.normalize();
                    let value = n1.dot(n2).abs();
                    if value >= 1.0 - ta * ta {
                        cont = GeomAbsConcatShape::G1;
                    }
                }
            }
        }
    } else {
        panic!("Standard_Failure: Courbes non jointives");
    }
    let _ = &mut cont1;
    let _ = &mut cont2;
    cont
}

/// OCCT static Continuity — the tl/ta default overload (Geom2dConvert.cxx
/// L856-864) with Precision::Confusion() and Precision::Angular().
pub fn continuity_default(
    c1: &Geom2dBSplineCurve,
    c2: &Geom2dBSplineCurve,
    u1: f64,
    u2: f64,
    r1: bool,
    r2: bool,
) -> GeomAbsConcatShape {
    continuity(c1, c2, u1, u2, r1, r2, CONFUSION, ANGULAR)
}

/// OCCT class Geom2dConvert_reparameterise_evaluator
/// (Geom2dConvert.cxx L868-895): evaluates the degree-2 reparameterizing
/// polynomial through PLib::EvalPolynomial(Parameter, DerivativeRequest, 2,
/// 1, myPolynomialCoefficient, Result).
pub struct Geom2dConvertReparameteriseEvaluator {
    my_polynomial_coefficient: [f64; 3],
}

impl Geom2dConvertReparameteriseEvaluator {
    /// OCCT ctor (L872-875).
    pub fn new(the_polynomial_coefficient: &[f64; 3]) -> Self {
        Geom2dConvertReparameteriseEvaluator {
            my_polynomial_coefficient: *the_polynomial_coefficient,
        }
    }
}

impl BSplCLibEvaluatorFunction for Geom2dConvertReparameteriseEvaluator {
    fn evaluate(
        &self,
        derivative_request: i32,
        _start_end: &[f64],
        parameter: f64,
        result: &mut f64,
        error_code: &mut i32,
    ) {
        // OCCT L883-891: PLib::EvalPolynomial(Parameter,
        // DerivativeRequest, 2, 1, coeffs, Result) — the degree-2,
        // dimension-1 polynomial and its derivative chain: c0 + c1*p +
        // c2*p^2, c1 + 2*c2*p, 2*c2.
        *error_code = 0;
        let c = &self.my_polynomial_coefficient;
        match derivative_request {
            0 => *result = c[0] + c[1] * parameter + c[2] * parameter * parameter,
            1 => *result = c[1] + 2.0 * c[2] * parameter,
            2 => *result = 2.0 * c[2],
            _ => *result = 0.0,
        }
    }
}

/// OCCT class Geom2dConvert_law_evaluator (Geom2dConvert.cxx L453-484): the
/// y-coordinate law of the anchor curve.
struct Geom2dConvertLawEvaluator<'a> {
    my_ancore: &'a Geom2dBSplineCurve,
}

impl<'a> BSplCLibEvaluatorFunction for Geom2dConvertLawEvaluator<'a> {
    fn evaluate(
        &self,
        derivative_request: i32,
        start_end: &[f64],
        parameter: f64,
        result: &mut f64,
        error_code: &mut i32,
    ) {
        *error_code = 0;
        // !myAncore.IsNull() is total on the reference.
        if parameter >= start_end[0] && parameter <= start_end[1] && derivative_request == 0 {
            let a_point = self.my_ancore.eval_d0(parameter);
            *result = a_point.y; // Coord(2)
        } else {
            *error_code = 1;
        }
    }
}

/// OCCT BSplCLib::MergeBSplineKnots (BSplCLib_2.cxx L1063-1173), hosted
/// locally (the math/bspl_lib.rs edit surface is fixed for this batch).
/// Returns (NumPoles, NewKnots, NewMults).
#[allow(clippy::too_many_arguments)]
fn merge_bspline_knots(
    tolerance: f64,
    start_value: f64,
    end_value: f64,
    degree1: usize,
    knots1: &[f64],
    mults1: &[i32],
    degree2: usize,
    knots2: &[f64],
    mults2: &[i32],
) -> (usize, Vec<f64>, Vec<i32>) {
    if start_value < end_value - tolerance {
        let degree = degree1 + degree2;
        let mut knots1 = knots1.to_vec();
        let mut knots2 = knots2.to_vec();
        // BSplCLib::Reparametrize(StartValue, EndValue, knots1/2).
        reparametrize(start_value, end_value, &mut knots1);
        reparametrize(start_value, end_value, &mut knots2);
        // Count pass (L1099-1116).
        let mut num_knots = 0i32;
        let mut jj = 1i32;
        for ii in 1..=knots1.len() as i32 {
            while jj <= knots2.len() as i32 && at(&knots2, jj) <= at(&knots1, ii) - tolerance {
                jj += 1;
                num_knots += 1;
            }
            while jj <= knots2.len() as i32 && at(&knots2, jj) <= at(&knots1, ii) + tolerance {
                jj += 1;
            }
            num_knots += 1;
        }
        let mut new_knots = vec![0.0f64; num_knots as usize];
        let mut new_mults = vec![0i32; num_knots as usize];
        // Fill pass (L1119-1148).
        num_knots = 1;
        jj = 1;
        for ii in 1..=knots1.len() as i32 {
            while jj <= knots2.len() as i32 && at(&knots2, jj) <= at(&knots1, ii) - tolerance {
                new_knots[(num_knots - 1) as usize] = at(&knots2, jj);
                new_mults[(num_knots - 1) as usize] = mults2[(jj - 1) as usize] + degree1 as i32;
                jj += 1;
                num_knots += 1;
            }
            let mut set_mults_flag = 0;
            while jj <= knots2.len() as i32 && at(&knots2, jj) <= at(&knots1, ii) + tolerance {
                let cont = (degree1 as i32 - mults1[(ii - 1) as usize])
                    .min(degree2 as i32 - mults2[(jj - 1) as usize]);
                set_mults_flag = 1;
                new_mults[(num_knots - 1) as usize] = degree as i32 - cont;
                jj += 1;
            }
            new_knots[(num_knots - 1) as usize] = at(&knots1, ii);
            if set_mults_flag == 0 {
                new_mults[(num_knots - 1) as usize] = mults1[(ii - 1) as usize] + degree2 as i32;
            }
            num_knots += 1;
        }
        num_knots -= 1;
        new_mults[0] = degree as i32 + 1;
        new_mults[(num_knots - 1) as usize] = degree as i32 + 1;
        let mut index = 0i32;
        for ii in 1..=num_knots {
            index += new_mults[(ii - 1) as usize];
        }
        let num_poles = (index - degree as i32 - 1) as usize;
        (num_poles, new_knots, new_mults)
    } else {
        // OCCT L1160-1172.
        let degree = degree1 + degree2;
        let new_knots = vec![start_value, end_value];
        let new_mults = vec![degree as i32 + 1, degree as i32 + 1];
        let num_poles = nb_poles(degree, false, &new_mults);
        (num_poles, new_knots, new_mults)
    }
}

/// OCCT BSplCLib::FunctionMultiply (BSplCLib_2.cxx L861-934), the flat
/// overload, hosted locally for the MultNumandDenom consumers.
#[allow(clippy::too_many_arguments)]
fn function_multiply(
    function: &dyn BSplCLibEvaluatorFunction,
    bspline_degree: usize,
    bspline_flat_knots: &[f64],
    poles_dimension: usize,
    poles: &[f64],
    flat_knots: &[f64],
    new_degree: usize,
    new_poles: &mut [f64],
    the_status: &mut i32,
) {
    let a_num_new_poles = flat_knots.len() as i32 - new_degree as i32 - 1;
    let a_start_end = [
        at(flat_knots, new_degree as i32 + 1),
        at(flat_knots, a_num_new_poles + 1),
    ];

    let mut a_parameters = vec![0.0f64; a_num_new_poles as usize];
    let mut a_contact_order_array = vec![0i32; a_num_new_poles as usize];
    let mut a_new_poles_array = vec![0.0f64; (a_num_new_poles * poles_dimension as i32) as usize];

    build_schoenberg_points(new_degree, flat_knots, &mut a_parameters);

    if a_parameters[0] < a_start_end[0] {
        a_parameters[0] = a_start_end[0];
    }
    if a_parameters[(a_num_new_poles - 1) as usize] > a_start_end[1] {
        a_parameters[(a_num_new_poles - 1) as usize] = a_start_end[1];
    }

    // int anExtrapMode = BSplineDegree (same ExtrapMode model as
    // function_reparameterise: Eval reads entries [0] and [1]).
    let mut an_extrap_mode = [bspline_degree as i32, bspline_degree as i32];
    let mut an_index = 0usize;
    for i in 1..=a_num_new_poles {
        a_contact_order_array[(i - 1) as usize] = 0;
        let mut a_result = 0.0f64;
        let mut an_error_code = 0i32;
        function.evaluate(
            a_contact_order_array[(i - 1) as usize],
            &a_start_end,
            a_parameters[(i - 1) as usize],
            &mut a_result,
            &mut an_error_code,
        );
        if an_error_code != 0 {
            *the_status = 1;
            return;
        }

        let dim = poles_dimension;
        eval_flat(
            a_parameters[(i - 1) as usize],
            false,
            0,
            &mut an_extrap_mode,
            bspline_degree,
            bspline_flat_knots,
            dim,
            poles,
            &mut a_new_poles_array[an_index..an_index + dim],
        );

        for _j in 0..poles_dimension {
            a_new_poles_array[an_index] *= a_result;
            an_index += 1;
        }
    }

    *the_status = interpolate(
        new_degree,
        flat_knots,
        &a_parameters,
        &a_contact_order_array,
        poles_dimension,
        &mut a_new_poles_array,
    );

    let total = (a_num_new_poles * poles_dimension as i32) as usize;
    new_poles[..total].copy_from_slice(&a_new_poles_array[..total]);
}

/// OCCT static MultNumandDenom (Geom2dConvert.cxx L488-567): multiplies the
/// denominator `bs` by the law `a` and returns the rational product curve.
pub fn mult_numand_denom(a: &Geom2dBSplineCurve, bs: &Geom2dBSplineCurve) -> Geom2dBSplineCurve {
    let bs_knots = curve_knots(bs);
    let bs_mults = curve_mults(bs);
    let mut bs_poles: Vec<DVec2> = (1..=bs.nb_poles_curve()).map(|i| bs.pole(i)).collect();
    let bs_weights = curve_weights(bs);
    let bs_flat_knots = curve_flat_knots(bs);
    let start_value = bs_knots[0];
    let end_value = bs_knots[(bs.nb_knots() - 1) as usize];
    let tolerance = 10.0 * epsilon_of(end_value.abs());

    let mut a_knots = curve_knots(a);
    let a_poles: Vec<DVec2> = (1..=a.nb_poles_curve()).map(|i| a.pole(i)).collect();
    let a_mults = curve_mults(a);
    // BSplCLib::Reparametrize(BS->FirstParameter(), BS->LastParameter(), aKnots).
    reparametrize(bs.first_parameter(), bs.last_parameter(), &mut a_knots);
    let an_ancore =
        Geom2dBSplineCurve::new(a_poles, a_knots.clone(), a_mults.clone(), a.degree(), false);

    let (res_nb_poles, res_knots, res_mults) = merge_bspline_knots(
        tolerance,
        start_value,
        end_value,
        a.degree(),
        &a_knots,
        &a_mults,
        bs.degree(),
        &bs_knots,
        &bs_mults,
    );
    let degree = bs.degree() + a.degree();
    let mut res_flat_knots = vec![0.0f64; res_nb_poles + degree + 1];
    knot_sequence(&res_knots, &res_mults, degree, false, &mut res_flat_knots);
    // Homogenize the BS poles by the BS weights (L532-538).
    for ii in 0..bs_poles.len() {
        bs_poles[ii].x *= bs_weights[ii];
        bs_poles[ii].y *= bs_weights[ii];
    }
    // POP for NT
    let ev = Geom2dConvertLawEvaluator { my_ancore: &an_ancore };
    let mut a_status = 0i32;
    let mut res_num_poles = vec![0.0f64; res_nb_poles * 2];
    let mut bs_poles_flat = Vec::with_capacity(bs_poles.len() * 2);
    for p in &bs_poles {
        bs_poles_flat.push(p.x);
        bs_poles_flat.push(p.y);
    }
    function_multiply(
        &ev,
        bs.degree(),
        &bs_flat_knots,
        2,
        &bs_poles_flat,
        &res_flat_knots,
        degree,
        &mut res_num_poles,
        &mut a_status,
    );
    let mut res_den_poles = vec![0.0f64; res_nb_poles];
    function_multiply(
        &ev,
        bs.degree(),
        &bs_flat_knots,
        1,
        &bs_weights,
        &res_flat_knots,
        degree,
        &mut res_den_poles,
        &mut a_status,
    );
    let mut res_poles = Vec::with_capacity(res_nb_poles);
    for ii in 0..res_nb_poles {
        res_poles.push(DVec2::new(
            res_num_poles[2 * ii] / res_den_poles[ii],
            res_num_poles[2 * ii + 1] / res_den_poles[ii],
        ));
    }
    Geom2dBSplineCurve::new_rational(res_poles, res_den_poles, res_knots, res_mults, degree, false)
}

/// OCCT static Pretreatment (Geom2dConvert.cxx L571-593): normalizes away
/// uniform non-unit end-anchored weights (value semantics: each SetWeight
/// rebuilds the curve).
pub fn pretreatment(tab: &mut [Geom2dBSplineCurve]) {
    for curve in tab.iter_mut() {
        if curve.is_rational() {
            let a = curve.weight(1);
            if (curve.weight(2) == a)
                && (curve.weight(curve.nb_poles_curve() - 1) == a)
                && (curve.weight(curve.nb_poles_curve()) == a)
            {
                for j in 1..=curve.nb_poles_curve() {
                    let w = curve.weight(j) / a;
                    *curve = set_weight(curve, j, w);
                }
            }
        }
    }
}

/// OCCT static NeedToBeTreated (Geom2dConvert.cxx L597-617).
pub fn need_to_be_treated(bs: &Geom2dBSplineCurve) -> bool {
    if bs.is_rational() {
        let tab_weights = curve_weights(bs);
        crate::math::hermit::bspl_eval_kernels::is_rational_window(
            &tab_weights,
            1,
            bs.nb_poles_curve(),
        ) && ((bs.weight(1) < (1.0 - CONFUSION))
            || (bs.weight(1) > (1.0 + CONFUSION))
            || (bs.weight(2) < (1.0 - CONFUSION))
            || (bs.weight(2) > (1.0 + CONFUSION))
            || (bs.weight(bs.nb_poles_curve() - 1) < (1.0 - CONFUSION))
            || (bs.weight(bs.nb_poles_curve() - 1) > (1.0 + CONFUSION))
            || (bs.weight(bs.nb_poles_curve()) < (1.0 - CONFUSION))
            || (bs.weight(bs.nb_poles_curve()) > (1.0 + CONFUSION)))
    } else {
        false
    }
}

/// OCCT static Need2DegRepara (Geom2dConvert.cxx L621-637).
pub fn need_2deg_repara(tab: &[Geom2dBSplineCurve]) -> bool {
    let mut rapport = 1.0f64;
    for i in 0..tab.len().saturating_sub(1) {
        let (_, vec1) = tab[i + 1].eval_d1(tab[i + 1].first_parameter());
        let (_, vec2) = tab[i].eval_d1(tab[i].last_parameter());
        rapport *= vec2.length() / vec1.length();
    }
    (rapport > (1.0 + CONFUSION)) || (rapport < (1.0 - CONFUSION))
}

/// OCCT static Indexmin (Geom2dConvert.cxx L641-655) — the LAST index
/// carrying the minimal degree (`<=` keeps it).
pub fn indexmin(tab: &[Geom2dBSplineCurve]) -> usize {
    let mut index = 0usize;
    let mut degree = tab[0].degree();
    for (i, curve) in tab.iter().enumerate() {
        if curve.degree() <= degree {
            degree = curve.degree();
            index = i;
        }
    }
    index
}

/// OCCT static ReorderArrayOfG1 (Geom2dConvert.cxx L659-712) — rotates the
/// closed junction to the front.  All arrays are 0-based as in OCCT; the
/// last bis entries of toler/g1 are left unset (OCCT default init) and
/// modeled with 0/false.
pub fn reorder_array_of_g1(
    array_of_curves: &mut [Geom2dBSplineCurve],
    array_of_toler: &mut [f64],
    tab_g1: &mut [bool],
    start_index: i32,
    closed_tolerance: f64,
) {
    let len = array_of_curves.len();
    let arraybis_of_curves: Vec<Geom2dBSplineCurve> = array_of_curves.to_vec();
    // The bis copies carry every entry but the last (OCCT leaves the last
    // toler/g1 slots at their default-init values).
    let mut arraybis_of_toler = vec![0.0f64; len];
    let mut tabbis_g1 = vec![false; len];
    for i in 0..len {
        if i != len - 1 {
            arraybis_of_toler[i] = array_of_toler[i];
            tabbis_g1[i] = tab_g1[i];
        }
    }
    let arraybis_of_toler = &arraybis_of_toler;
    let tabbis_g1 = &tabbis_g1;

    for i in 0..=(len as i32 - (start_index + 2)) {
        let i = i as usize;
        array_of_curves[i] = arraybis_of_curves[i + start_index as usize + 1].clone();
        if i as i32 != (len as i32 - (start_index + 2)) {
            array_of_toler[i] = arraybis_of_toler[i + start_index as usize + 1];
            tab_g1[i] = tabbis_g1[i + start_index as usize + 1];
        }
    }

    array_of_toler[(len as i32 - (start_index + 2)) as usize] = closed_tolerance;
    tab_g1[(len as i32 - (start_index + 2)) as usize] = true;

    for i in (len as i32 - (start_index + 1))..=(len as i32 - 1) {
        let i = i as usize;
        if i != len - 1 {
            array_of_curves[i] = arraybis_of_curves[i - (len - (start_index as usize + 1))].clone();
            array_of_toler[i] = arraybis_of_toler[i - (len - (start_index as usize + 1))];
            tab_g1[i] = tabbis_g1[i - (len - (start_index as usize + 1))];
        } else {
            array_of_curves[i] = arraybis_of_curves[i - (len - (start_index as usize + 1))].clone();
        }
    }
}

/// gp_Vec2d::Angle (gp_Vec2d.cxx L47-85) — the acos/asin branch pair.
fn gp_vec2d_angle(coord: DVec2, other: DVec2) -> f64 {
    let a_norm = coord.length();
    let an_other_norm = other.length();
    if a_norm <= GP_RESOLUTION || an_other_norm <= GP_RESOLUTION {
        panic!("gp_VectorWithNullMagnitude");
    }
    let a_d = a_norm * an_other_norm;
    let a_cosinus = coord.dot(other) / a_d;
    let a_sinus = (coord.x * other.y - coord.y * other.x) / a_d;
    let a_cos_45_deg = std::f64::consts::FRAC_1_SQRT_2;
    if a_cosinus > -a_cos_45_deg && a_cosinus < a_cos_45_deg {
        if a_sinus > 0.0 {
            a_cosinus.acos()
        } else {
            -a_cosinus.acos()
        }
    } else if a_cosinus > 0.0 {
        a_sinus.asin()
    } else if a_sinus > 0.0 {
        std::f64::consts::PI - a_sinus.asin()
    } else {
        -std::f64::consts::PI - a_sinus.asin()
    }
}

/// gp_Vec2d::IsParallel (gp_Vec2d.hxx L372-376).
pub fn gp_vec2d_is_parallel(v: DVec2, other: DVec2, angular_tolerance: f64) -> bool {
    let an_ang = gp_vec2d_angle(v, other).abs();
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// Architecture bridge (value semantics): pack the kernel curve into the
/// legacy flat representation consumed by the existing
/// Geom2dConvert_CompCurveToBSplineCurve port (knots expanded per
/// multiplicity; weights carried).
fn to_bspline2(bs: &Geom2dBSplineCurve) -> BSplineCurve2 {
    let mut knots = Vec::new();
    for i in 1..=bs.nb_knots() {
        let k = bs.knot(i);
        let m = bs.multiplicity(i);
        for _ in 0..m {
            knots.push(k);
        }
    }
    BSplineCurve2 {
        degree: bs.degree(),
        knots,
        control_points: (1..=bs.nb_poles_curve()).map(|i| bs.pole(i)).collect(),
        weights: (1..=bs.nb_poles_curve()).map(|i| bs.weight(i)).collect(),
    }
}

/// OCCT Geom2dConvert_CompCurveToBSplineCurve C(Curve2);
/// C.Add(Curve1, Tolerance[, After]); Curve2 = C.BSplineCurve()
/// (CompCurveToBSplineCurve.cxx L57-244, via the existing port) with the
/// OCCT failure throw.  The legacy BSplineCurve2 result is unpacked with the
/// faithful from_bspline2 (non-periodic: the OCCT accumulator is only made
/// periodic by the explicit SetPeriodic call on the ConcatC1 result).
pub fn comp_curve_add(
    first: &Geom2dBSplineCurve,
    new_curve: &Geom2dBSplineCurve,
    tolerance: f64,
    after: bool,
) -> Geom2dBSplineCurve {
    let mut comp = Geom2dConvertCompCurveToBSplineCurve::with_basis_curve(
        &Curve2d::BSpline(to_bspline2(first)),
        ConvertParameterisation::TgtThetaOver2,
    );
    let fusion = comp.add(&Curve2d::BSpline(to_bspline2(new_curve)), tolerance);
    if !fusion {
        panic!("Standard_ConstructionError: Geom2dConvert Concatenation Error");
    }
    let merged = comp.bspline_curve().expect("Geom2dConvert myCurve");
    Geom2dBSplineCurve::from_bspline2(&merged, false)
}

/// OCCT Geom2dConvert::ConcatC1 — the default AngularTolerance overload
/// (Geom2dConvert.cxx L1149-1164, Precision::Angular()).
pub fn concat_c1_default(
    array_of_curves: &mut [Geom2dBSplineCurve],
    array_of_toler: &[f64],
    closed_flag: &mut bool,
    closed_tolerance: f64,
) -> (Vec<i32>, Vec<Geom2dBSplineCurve>) {
    concat_c1(
        array_of_curves,
        array_of_toler,
        closed_flag,
        closed_tolerance,
        ANGULAR,
    )
}

/// OCCT Geom2dConvert::ConcatC1 (Geom2dConvert.cxx L1168-1455) —
/// concatenates the G1 array of curves into the C1 array of concatenated
/// curves; returns (ArrayOfIndices, ArrayOfConcatenated) (the OCCT out
/// handles).
#[allow(clippy::too_many_arguments)]
pub fn concat_c1(
    array_of_curves: &mut [Geom2dBSplineCurve],
    array_of_toler: &[f64],
    closed_flag: &mut bool,
    closed_tolerance: f64,
    angular_tolerance: f64,
) -> (Vec<i32>, Vec<Geom2dBSplineCurve>) {
    let nb_curve = array_of_curves.len() as i32;
    let mut nb_vertex_g1;
    let mut nb_group = 0i32;
    let mut index = 0i32;
    let mut nb_vertex_group0 = 0i32;
    let mut pre_last = 0.0f64;
    let mut tab_g1 = vec![false; (nb_curve - 1).max(0) as usize]; // G1 continuity table at junctions
    let mut local_tolerance = array_of_toler.to_vec();

    for i in 0..array_of_toler.len() {
        local_tolerance[i] = array_of_toler[i];
    }
    for i in 0..nb_curve {
        if i >= 1 {
            let first = array_of_curves[i as usize].first_parameter();
            if continuity(
                &array_of_curves[(i - 1) as usize],
                &array_of_curves[i as usize],
                pre_last,
                first,
                true,
                true,
                array_of_toler[(i - 1) as usize],
                angular_tolerance,
            ) < GeomAbsConcatShape::C0
            {
                panic!("Standard_ConstructionError: Geom2dConvert curves not C0"); // renvoi d'une erreur
            } else {
                if continuity(
                    &array_of_curves[(i - 1) as usize],
                    &array_of_curves[i as usize],
                    pre_last,
                    first,
                    true,
                    true,
                    array_of_toler[(i - 1) as usize],
                    angular_tolerance,
                ) >= GeomAbsConcatShape::G1
                {
                    tab_g1[(i - 1) as usize] = true; // True=Continuite G1
                } else {
                    tab_g1[(i - 1) as usize] = false;
                }
            }
        }
        pre_last = array_of_curves[i as usize].last_parameter();
    }

    // Determination of Wire characteristics.
    while index <= nb_curve - 1 {
        nb_vertex_g1 = 0;
        while ((index + nb_vertex_g1) <= nb_curve - 2)
            && tab_g1[(index + nb_vertex_g1) as usize]
        {
            nb_vertex_g1 += 1;
        }
        nb_group += 1;
        if index == 0 {
            nb_vertex_group0 = nb_vertex_g1;
        }
        index = index + 1 + nb_vertex_g1;
    }

    if *closed_flag && nb_group != 1 {
        // rearrangement du tableau
        nb_group -= 1;
        reorder_array_of_g1(
            array_of_curves,
            &mut local_tolerance,
            &mut tab_g1,
            nb_vertex_group0,
            closed_tolerance,
        );
    }

    let mut array_of_indices = vec![0i32; (nb_group + 1) as usize];
    let mut array_of_concatenated: Vec<Option<Geom2dBSplineCurve>> =
        (0..nb_group).map(|_| None).collect();

    let mut k = 0i32;
    index = 0;
    pretreatment(array_of_curves);
    let mut a_polynomial_coefficient = [0.0f64; 3];

    let need_double_deg_repara = need_2deg_repara(array_of_curves);
    if nb_group == 1 && *closed_flag && need_double_deg_repara {
        let curve1 = &array_of_curves[(nb_curve - 1) as usize];
        if curve1.degree() > max_degree() / 2 {
            *closed_flag = false;
        }
    }

    if nb_group == 1 && *closed_flag {
        // treatment of a particular case
        array_of_indices[0] = 0;
        array_of_indices[1] = 0;
        let indexmin = indexmin(array_of_curves);
        if indexmin != (array_of_curves.len() - 1) {
            reorder_array_of_g1(
                array_of_curves,
                &mut local_tolerance,
                &mut tab_g1,
                indexmin as i32,
                closed_tolerance,
            );
        }
        let mut curve2_opt: Option<Geom2dBSplineCurve> = None;
        for j in 0..nb_curve {
            // secondary loop inside each group
            let j = j as usize;
            let mut curve1 = if need_to_be_treated(&array_of_curves[j]) {
                let solution = crate::math::hermit::hermit_solution_2d(
                    &array_of_curves[j],
                    1.0e-6,
                    1.0e-6,
                );
                mult_numand_denom(&solution, &array_of_curves[j])
            } else {
                array_of_curves[j].clone()
            };

            let a_new_curve_degree = 2 * curve1.degree();

            if j == 0 {
                // initialisation en debut de groupe
                curve2_opt = Some(curve1);
            } else {
                if (j as i32 == (nb_curve - 1)) && need_double_deg_repara {
                    let curve2_ref = curve2_opt.as_ref().expect("Curve2");
                    let (_, vec1) = curve2_ref.eval_d1(curve2_ref.last_parameter());
                    let (_, vec2) = curve1.eval_d1(curve1.first_parameter());
                    let lambda = vec2.length() / vec1.length();
                    let mut knot_c1 = curve_knots(&curve1);
                    let (_, vec2) = curve1.eval_d1(curve1.last_parameter());
                    let (_, vec1) =
                        array_of_curves[0].eval_d1(array_of_curves[0].first_parameter());
                    let lambda2 = vec1.length() / vec2.length();
                    let umin = curve1.first_parameter();
                    let umax = curve1.last_parameter();
                    let tmax = 2.0 * lambda * (umax - umin) / (1.0 + lambda * lambda2);
                    let a_coef = (lambda * lambda2 - 1.0) / (2.0 * lambda * tmax);
                    a_polynomial_coefficient[2] = a_coef;
                    let b = 1.0 / lambda;
                    a_polynomial_coefficient[1] = b;
                    let c = umin;
                    a_polynomial_coefficient[0] = c;
                    let mut curve1_flat_knots =
                        vec![0.0f64; curve1.nb_poles_curve() as usize + curve1.degree() + 1];
                    let mut knot_c1_mults = curve_mults(&curve1);
                    knot_sequence(
                        &knot_c1,
                        &knot_c1_mults,
                        curve1.degree(),
                        false,
                        &mut curve1_flat_knots,
                    );
                    knot_c1[0] = 0.0;
                    for ii in 1..knot_c1.len() {
                        knot_c1[ii] =
                            (-b + (b * b - 4.0 * a_coef * (c - knot_c1[ii])).sqrt())
                                / (2.0 * a_coef); // ifv 17.05.00 buc60667
                    }
                    let mut curve1_poles: Vec<DVec2> =
                        (1..=curve1.nb_poles_curve()).map(|i| curve1.pole(i)).collect();
                    let curve1_weights = curve_weights(&curve1);
                    for ii in 0..curve1_poles.len() {
                        curve1_poles[ii].x *= curve1_weights[ii];
                        curve1_poles[ii].y *= curve1_weights[ii];
                    }
                    for ii in 0..knot_c1_mults.len() {
                        knot_c1_mults[ii] = curve1.degree() as i32 + knot_c1_mults[ii];
                    }
                    let mut flat_knots = vec![
                        0.0f64;
                        curve1_flat_knots.len()
                            + (curve1.degree() * curve1.nb_knots() as usize)
                    ];
                    knot_sequence(
                        &knot_c1,
                        &knot_c1_mults,
                        a_new_curve_degree,
                        false,
                        &mut flat_knots,
                    );
                    let new_poles_count = flat_knots.len() - (a_new_curve_degree + 1);
                    let mut new_poles = vec![0.0f64; 2 * new_poles_count];
                    let mut a_status = 0i32;
                    // Size checks from BSplCLib_FunctionReparameterise
                    // (CurveComputation.pxx L2028-2031).
                    assert_eq!(
                        curve1_poles.len(),
                        curve1_flat_knots.len() - curve1.degree() - 1,
                        "FunctionReparameterise: poles count mismatch"
                    );
                    let mut curve1_poles_hom = Vec::with_capacity(curve1_poles.len() * 2);
                    for p in &curve1_poles {
                        curve1_poles_hom.push(p.x);
                        curve1_poles_hom.push(p.y);
                    }
                    let ev = Geom2dConvertReparameteriseEvaluator::new(&a_polynomial_coefficient);
                    crate::math::bspl_lib::function_reparameterise(
                        &ev,
                        curve1.degree() as i32,
                        &curve1_flat_knots,
                        2,
                        &curve1_poles_hom,
                        &flat_knots,
                        a_new_curve_degree as i32,
                        &mut new_poles,
                        &mut a_status,
                    );
                    let mut new_weights = vec![0.0f64; new_poles_count];
                    crate::math::bspl_lib::function_reparameterise(
                        &ev,
                        curve1.degree() as i32,
                        &curve1_flat_knots,
                        1,
                        &curve1_weights,
                        &flat_knots,
                        a_new_curve_degree as i32,
                        &mut new_weights,
                        &mut a_status,
                    );
                    let mut new_poles_vec: Vec<DVec2> = Vec::with_capacity(new_poles_count);
                    for ii in 0..new_poles_count {
                        new_poles_vec.push(DVec2::new(
                            new_poles[2 * ii] / new_weights[ii],
                            new_poles[2 * ii + 1] / new_weights[ii],
                        ));
                    }
                    curve1 = Geom2dBSplineCurve::new_rational(
                        new_poles_vec,
                        new_weights,
                        knot_c1,
                        knot_c1_mults,
                        a_new_curve_degree,
                        false,
                    );
                }
                let curve2_ref = curve2_opt.as_ref().expect("Curve2");
                let merged = comp_curve_add(curve2_ref, &curve1, local_tolerance[j - 1], false);
                curve2_opt = Some(merged);
            }
        }
        let mut curve2 = curve2_opt.expect("Curve2");
        curve2.set_periodic(); // single C1 curve
        let last = curve2.last_uknot_index();
        let m = curve2.multiplicity(last) - 1;
        if let Some(reduced) = remove_knot_curve(&curve2, last, m, CONFUSION) {
            curve2 = reduced;
        }
        array_of_concatenated[0] = Some(curve2);
    } else {
        // boucle principale sur chaque groupe de continuite interne G1
        for i in 0..nb_group {
            nb_vertex_g1 = 0;
            while ((index + nb_vertex_g1) <= nb_curve - 2)
                && tab_g1[(index + nb_vertex_g1) as usize]
            {
                nb_vertex_g1 += 1;
            }

            if (!*closed_flag) || (nb_group == 1) {
                // Filling the array of preserved indices
                k += 1;
                array_of_indices[(k - 1) as usize] = index;
                if k == nb_group {
                    array_of_indices[k as usize] = 0;
                }
            } else {
                k += 1;
                array_of_indices[(k - 1) as usize] = index + nb_vertex_group0 + 1;
                if k == nb_group {
                    array_of_indices[k as usize] = nb_vertex_group0 + 1;
                }
            }

            for j in index..=(index + nb_vertex_g1) {
                let j = j as usize;
                let curve1 = if need_to_be_treated(&array_of_curves[j]) {
                    let solution = crate::math::hermit::hermit_solution_2d(
                        &array_of_curves[j],
                        1.0e-6,
                        1.0e-6,
                    );
                    mult_numand_denom(&solution, &array_of_curves[j])
                } else {
                    array_of_curves[j].clone()
                };

                if index == j as i32 {
                    // initialisation en debut de groupe
                    array_of_concatenated[i as usize] = Some(curve1);
                } else {
                    let current =
                        array_of_concatenated[i as usize].as_ref().expect("concatenated");
                    let merged = comp_curve_add(current, &curve1, array_of_toler[j - 1], false);
                    array_of_concatenated[i as usize] = Some(merged);
                }
            }
            index = index + 1 + nb_vertex_g1;
        }
    }
    (
        array_of_indices,
        array_of_concatenated
            .into_iter()
            .map(|c| c.expect("ArrayOfConcatenated entry"))
            .collect(),
    )
}

// ===========================================================================
// Tests
// ===========================================================================
#[cfg(test)]
mod c1_concat_tests {
    use super::*;

    fn line(p0: DVec2, p1: DVec2) -> Geom2dBSplineCurve {
        Geom2dBSplineCurve::new(vec![p0, p1], vec![0.0, 1.0], vec![2, 2], 1, false)
    }

    /// OCCT GeomAbsToInteger table (L721-742).
    #[test]
    fn geom_abs_to_integer_table() {
        assert_eq!(geom_abs_to_integer(GeomAbsConcatShape::C0), 0);
        assert_eq!(geom_abs_to_integer(GeomAbsConcatShape::G1), 1);
        assert_eq!(geom_abs_to_integer(GeomAbsConcatShape::C1), 2);
        assert_eq!(geom_abs_to_integer(GeomAbsConcatShape::G2), 3);
        assert_eq!(geom_abs_to_integer(GeomAbsConcatShape::C2), 4);
        assert_eq!(geom_abs_to_integer(GeomAbsConcatShape::C3), 5);
        assert_eq!(geom_abs_to_integer(GeomAbsConcatShape::CN), 6);
    }

    /// Indexmin keeps the LAST minimal degree (`<=` in OCCT L648): degrees
    /// [3, 2, 2] -> 2; [2, 3, 4] -> 0; [4, 2, 3] -> 1.
    #[test]
    fn indexmin_last_minimum() {
        let mk = |d: usize| {
            let n = d + 1;
            Geom2dBSplineCurve::new(
                (0..n).map(|i| DVec2::new(i as f64, 0.0)).collect(),
                vec![0.0, 1.0],
                vec![n as i32, n as i32],
                d,
                false,
            )
        };
        let tab = vec![mk(3), mk(2), mk(2)];
        assert_eq!(indexmin(&tab), 2);
        let tab = vec![mk(2), mk(3), mk(4)];
        assert_eq!(indexmin(&tab), 0);
        let tab = vec![mk(4), mk(2), mk(3)];
        assert_eq!(indexmin(&tab), 1);
    }

    /// ReorderArrayOfG1 with 4 curves and StartIndex=1 (hand-derived from
    /// OCCT L659-712): curves rotate to [c2, c3, c0, c1]; the junction slot
    /// len-(Start+2) = 1 receives ClosedTolerance with tabG1[1] = true; the
    /// remaining tolerances/g1 permute as [t2, T, t0, t1] / [g2, true, g0,
    /// g1] (t3/g3 unset in the bis copy).
    #[test]
    fn reorder_array_of_g1_permutation() {
        let mk = |x: f64| line(DVec2::new(x, 0.0), DVec2::new(x + 1.0, 0.0));
        let mut curves = vec![mk(0.0), mk(10.0), mk(20.0), mk(30.0)];
        let mut toler = vec![1.0, 2.0, 3.0, 0.0];
        let mut g1 = vec![false, false, true, false];
        reorder_array_of_g1(&mut curves, &mut toler, &mut g1, 1, 7.5);
        let xs: Vec<f64> = curves.iter().map(|c| c.pole(1).x).collect();
        assert_eq!(xs, vec![20.0, 30.0, 0.0, 10.0]);
        assert_eq!(toler, vec![3.0, 7.5, 1.0, 0.0]);
        assert_eq!(g1, vec![true, true, false, false]);
    }

    /// Need2DegRepara on two collinear unit lines is false (Rapport = 1);
    /// with the second line twice as fast (tangent magnitude 2) Rapport =
    /// 1/2 and the answer is true.
    #[test]
    fn need_2deg_repara_ratio() {
        let a = line(DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0));
        let b = line(DVec2::new(1.0, 0.0), DVec2::new(2.0, 0.0));
        assert!(!need_2deg_repara(&[a.clone(), b.clone()]));
        let fast_b = line(DVec2::new(1.0, 0.0), DVec2::new(3.0, 0.0));
        assert!(need_2deg_repara(&[a, fast_b]));
    }

    /// Continuity at the junction: two collinear unit lines join C1
    /// (parallel equal-length tangents); a right-angle corner stays C0.
    #[test]
    fn continuity_junction_classification() {
        let a = line(DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0));
        let b = line(DVec2::new(1.0, 0.0), DVec2::new(2.0, 0.0));
        let tol = 1.0e-7;
        let cont = continuity(&a, &b, 1.0, 0.0, true, true, tol, ANGULAR);
        assert_eq!(cont, GeomAbsConcatShape::C1);
        let corner = line(DVec2::new(1.0, 0.0), DVec2::new(1.0, 1.0));
        let cont = continuity(&a, &corner, 1.0, 0.0, true, true, tol, ANGULAR);
        assert_eq!(cont, GeomAbsConcatShape::C0);
    }

    /// MultNumandDenom of the constant law a = 1/2 (line y = 1/2) with the
    /// unit line bs = (0,0)-(1,0), hand-derived: the merged knot table is
    /// [0,1] x [3,3] (a degree-2 Bezier), the law evaluates to 1/2 at the
    /// Schoenberg points {1/3, 2/3, 1}, the numerator interpolates to x
    /// poles [0, 1/4, 1/2] and the denominator to [1/2, 1/2, 1/2], so the
    /// quotient recovers bs raised to degree 2: poles (0,0), (1/2,0), (1,0)
    /// with uniform weights (non-rational after the ctor rationality fixup).
    #[test]
    fn mult_numand_denom_constant_law() {
        let a = line(DVec2::new(0.0, 0.5), DVec2::new(1.0, 0.5));
        let bs = line(DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0));
        let res = mult_numand_denom(&a, &bs);
        assert_eq!(res.degree(), 2);
        assert_eq!(res.nb_knots(), 2);
        assert_eq!((res.knot(1), res.knot(2)), (0.0, 1.0));
        assert_eq!((res.multiplicity(1), res.multiplicity(2)), (3, 3));
        assert_eq!(res.nb_poles_curve(), 3);
        assert!(!res.is_rational());
        let expected = [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0)];
        for (i, &(x, y)) in expected.iter().enumerate() {
            let p = res.pole(i as i32 + 1);
            assert!((p.x - x).abs() < 1e-13 && (p.y - y).abs() < 1e-13);
        }
    }

    /// ConcatC1 of two C1-connected unit lines [0,1] and [1,2]:
    /// tabG1 = [true] (parallel equal tangents), one group, non-closed, no
    /// rational treatment, and the accumulator merges into the degree-1
    /// curve on [0,2] with poles (0,0), (2,0) (junction knot removed down to
    /// multiplicity 0: the collinear removal deviation is 0).
    /// ArrayOfIndices = [0, 0] (single preserved group).
    #[test]
    fn concat_c1_two_lines() {
        let a = line(DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0));
        let b = line(DVec2::new(1.0, 0.0), DVec2::new(2.0, 0.0));
        let mut curves = vec![a, b];
        let toler = vec![1.0e-7];
        let mut closed = false;
        let (indices, concatenated) = concat_c1_default(&mut curves, &toler, &mut closed, 0.0);
        assert!(!closed);
        assert_eq!(indices, vec![0, 0]);
        assert_eq!(concatenated.len(), 1);
        let c = &concatenated[0];
        assert_eq!(c.degree(), 1);
        assert_eq!(c.nb_knots(), 2);
        assert_eq!((c.knot(1), c.knot(2)), (0.0, 2.0));
        assert_eq!(c.nb_poles_curve(), 2);
        assert!((c.pole(1).x - 0.0).abs() < 1e-13 && c.pole(1).y.abs() < 1e-13);
        assert!((c.pole(2).x - 2.0).abs() < 1e-13 && c.pole(2).y.abs() < 1e-13);
    }
}
