//! OCCT ChFi3d_PerformElSpine (TKFillet/ChFi3d, Builder_0.cxx L4797-5327)
//! and the ChFi3d static helpers it consumes (CurveCleaner L3659-3681,
//! GoodExt L4773-4793, ChFi3d_ApproxByC2 L5835-5851, ChFi3d_IsSmooth
//! L5855-5944).
//!
//! The Geom/GeomConvert/GeomLib primitives the function needs live here as
//! local translations because the owning kernel batches have not landed:
//!   - Geom_BSplineCurve methods (TKG3d/Geom/Geom_BSplineCurve.cxx):
//!     RemoveKnot L420-492, Reverse L496-516, Segment L527-715,
//!     SetPeriodic L777-815, SetOrigin(Index) L819-909,
//!     SetOrigin(U, Tol) L913-970, SetNotPeriodic L974-1016,
//!     InsertKnots L351-416.
//!   - GeomConvert_CompCurveToBSplineCurve (TKGeomBase/GeomConvert/
//!     GeomConvert_CompCurveToBSplineCurve.cxx L32-274).
//!   - GeomLib::ExtendCurveToPoint (TKGeomBase/GeomLib/GeomLib.cxx
//!     L1269-1413) with its dependencies PLib::HermiteCoefficients
//!     (TKMath/PLib/PLib.cxx L1404-1478), PLib::CoefficientsPoles
//!     (L1482-1605) and the static GeomLib::ComputeLambda (GeomLib.cxx
//!     L126-253).
//!   - GeomLib::AdjustExtremity / GeomConvert::C0BSplineToC1BSplineCurve /
//!     GeomAPI_PointsToBSpline — GAP carriers (annotated at the call
//!     sites, OCCT failure paths preserved).
//!
//! Architecture mappings: `occ::handle<Geom_Curve>` -> `Curve3`;
//! `Geom_BSplineCurve` -> `BSplineCurve3` (rcad stores the flat knot
//! vector; the (knots, mults) form is derived and rebuilt by the helpers);
//! `occ::handle<ChFiDS_ElSpine>&` -> `&mut ChFiDSElSpine`;
//! `occ::handle<ChFiDS_Spine>&` -> `&mut ChFiDSSpineHandle` (the spine
//! Parameter/D1 queries take &mut in the rcad ChFiDS translation);
//! `GeomAbs_Shape` continuity -> the rcad GeomAbsShape enum.

use glam::DVec3;
use rcad_kernel::base::convert::{geom_convert_curve_to_bspline_curve, ConvertParameterisation};
use rcad_kernel::base::gcpnts::abscissa_point::abscissa_point_parameter;
use rcad_kernel::geom::{
    BSplineCurve3, Circle3, Curve3, CurveEval as _, Ellipse3, Line3, TrimmedCurve3,
};
use rcad_kernel::math::bspl_lib;
use rcad_kernel::math::bspl_lib::{ati, at};
use rcad_kernel::math::bspl::de_boor;
use rcad_kernel::math::el as elclib;
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::topo::topods::{Orientation, Shape};

use super::chfi3d_builder_0::topexp_common_vertex;
use super::chfi3d_geom_lib::{geom_lib_adjust_extremity, geom_lib_extend_curve_to_point};
use super::chfi_ds::{ChFiDSElSpine, ChFiDSSpineHandle};
use crate::geomalgo::gtests_stubs::GeomAbsShape as ChFiDSGeomAbsShape;

// =========================================================================
// OCCT Standard_Real.hxx L242-248 — Epsilon(theValue): one ULP toward the
// same-sign infinity.
// =========================================================================
fn epsilon(the_value: f64) -> f64 {
    if the_value >= 0.0 {
        next_after_up(the_value) - the_value
    } else {
        the_value - next_after_down(the_value)
    }
}

fn next_after_up(x: f64) -> f64 {
    if x == 0.0 || x.is_nan() {
        return f64::MIN_POSITIVE;
    }
    let bits = x.to_bits();
    if x > 0.0 {
        f64::from_bits(bits + 1)
    } else {
        f64::from_bits(bits - 1)
    }
}

fn next_after_down(x: f64) -> f64 {
    if x == 0.0 {
        return -f64::MIN_POSITIVE;
    }
    let bits = x.to_bits();
    if x > 0.0 {
        f64::from_bits(bits - 1)
    } else {
        f64::from_bits(bits + 1)
    }
}

// =========================================================================
// Geom_BSplineCurve accessors over the rcad flat-knot BSplineCurve3.
// =========================================================================

/// OCCT Geom_BSplineCurve::NbKnots().
fn bspline_nb_knots(bs: &BSplineCurve3) -> usize {
    bs.knots_mults().0.len()
}

/// OCCT Geom_BSplineCurve::Knot(Index) (1-based).
fn bspline_knot(bs: &BSplineCurve3, index: i32) -> f64 {
    bs.knots_mults().0[(index - 1) as usize]
}

/// OCCT Geom_BSplineCurve::Multiplicity(Index) (1-based).
fn bspline_multiplicity(bs: &BSplineCurve3, index: i32) -> i32 {
    bs.knots_mults().1[(index - 1) as usize]
}

/// OCCT Geom_BSplineCurve::NbPoles().
fn bspline_nb_poles(bs: &BSplineCurve3) -> usize {
    bs.control_points.len()
}

/// OCCT Geom_BSplineCurve::Pole(Index) (1-based).
fn bspline_pole(bs: &BSplineCurve3, index: i32) -> DVec3 {
    bs.control_points[(index - 1) as usize]
}

// OCCT Geom_BSplineCurve::FirstUKnotIndex()/LastUKnotIndex() route through
// bspl_lib::first_uknot_index_mults / last_uknot_index_mults at their call
// sites (RemoveKnot L427-428, SetPeriodic L780-781, SetOrigin L829-830).

/// OCCT Geom_BSplineCurve::Value(U) — the pole evaluation.
pub(crate) fn bspline_value(bs: &BSplineCurve3, u: f64) -> DVec3 {
    de_boor(bs.degree, &bs.knots, &bs.control_points, &bs.weights, u)
}

/// OCCT Geom_BSplineCurve::D1 magnitude query via DN(U, 1).
fn bspline_d1(bs: &BSplineCurve3, u: f64) -> DVec3 {
    bs.dn(u, 1)
}

/// Rebuild the flat knot vector from OCCT (knots, mults).
pub(crate) fn flat_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut flat = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            flat.push(*k);
        }
    }
    flat
}

fn poles_flat(control_points: &[DVec3]) -> Vec<f64> {
    control_points.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

fn unflatten_poles(flat: &[f64]) -> Vec<DVec3> {
    (0..flat.len() / 3)
        .map(|i| DVec3::new(flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]))
        .collect()
}

// =========================================================================
// OCCT Geom_BSplineCurve.cxx L420-492 — RemoveKnot(Index, M, Tolerance).
// GAP: the kernel BSplCLib::RemoveKnot carries the non-rational 3D path
// (weights arrays pending the TKMath batch); rational inputs fall through
// with unit weights (the OCCT non-rational arithmetic is exact).
// =========================================================================
pub(crate) fn bspline_remove_knot(bs: &mut BSplineCurve3, index: i32, m: i32, tolerance: f64) -> bool {
    if m < 0 {
        return true;
    }
    let (knots, mults) = bs.knots_mults();
    let i1 = bspl_lib::first_uknot_index_mults(bs.degree, &mults);
    let i2 = bspl_lib::last_uknot_index_mults(bs.degree, &mults);
    if !bs.is_periodic && (index <= i1 || index >= i2) {
        panic!("Standard_OutOfRange: BSpline curve: RemoveKnot: index out of range");
    } else if bs.is_periodic && (index < i1 || index > i2) {
        panic!("Standard_OutOfRange: BSpline curve: RemoveKnot: index out of range");
    }

    // OCCT L441: int step = myMults.Value(Index) - M;
    let step = bspl_lib::ati(&mults, index) - m;
    if step <= 0 {
        return true;
    }

    let old_poles_len = bspline_nb_poles(bs);
    let new_poles_len = old_poles_len - step as usize;
    let new_knots_len = knots.len() - if m == 0 { 1 } else { 0 };
    let poles = poles_flat(&bs.control_points);
    // OCCT passes Weights() — the kernel helper models the unit-weight path.
    let mut new_poles = vec![0.0f64; 3 * new_poles_len];
    let mut new_knots = vec![0.0f64; new_knots_len];
    let mut new_mults = vec![0i32; new_knots_len];
    let ok = bspl_lib::remove_knot(
        index as usize,
        m,
        bs.degree,
        bs.is_periodic,
        3,
        &poles,
        &knots,
        &mults,
        &mut new_poles,
        &mut new_knots,
        &mut new_mults,
        tolerance,
    );
    if !ok {
        return false;
    }
    bs.control_points = unflatten_poles(&new_poles);
    bs.weights = vec![1.0; new_poles_len];
    bs.knots = flat_knots(&new_knots, &new_mults);
    true
}

// =========================================================================
// OCCT Geom_BSplineCurve.cxx L351-416 — InsertKnots(Knots, Mults, Epsilon,
// Add).  GAP: the kernel BSplCLib::InsertKnots carries the non-rational 3D
// path (weights arrays pending the TKMath batch).
// =========================================================================
fn bspline_insert_knots(bs: &mut BSplineCurve3, knots_add: &[f64], mults_add: &[i32], eps: f64) {
    let (knots, mults) = bs.knots_mults();
    let mut nb_poles = 0i32;
    let mut nb_knots = 0i32;
    if !bspl_lib::prepare_insert_knots(
        bs.degree,
        bs.is_periodic,
        &knots,
        &mults,
        knots_add,
        Some(mults_add),
        &mut nb_poles,
        &mut nb_knots,
        eps,
        false,
    ) {
        panic!("Standard_ConstructionError: Geom_BSplineCurve::InsertKnots");
    }
    if nb_poles as usize == bspline_nb_poles(bs) {
        return;
    }
    let poles = poles_flat(&bs.control_points);
    let mut new_poles = vec![0.0f64; 3 * nb_poles as usize];
    let mut new_knots = vec![0.0f64; nb_knots as usize];
    let mut new_mults = vec![0i32; nb_knots as usize];
    bspl_lib::insert_knots(
        bs.degree,
        bs.is_periodic,
        3,
        &poles,
        &knots,
        &mults,
        knots_add,
        Some(mults_add),
        &mut new_poles,
        &mut new_knots,
        &mut new_mults,
        eps,
        false,
    );
    bs.control_points = unflatten_poles(&new_poles);
    bs.weights = vec![1.0; nb_poles as usize];
    bs.knots = flat_knots(&new_knots, &new_mults);
}

// =========================================================================
// OCCT Geom_BSplineCurve.cxx L777-815 — SetPeriodic().
// =========================================================================
pub(crate) fn bspline_set_periodic(bs: &mut BSplineCurve3) {
    let (knots, mults) = bs.knots_mults();
    let first = bspl_lib::first_uknot_index_mults(bs.degree, &mults);
    let last = bspl_lib::last_uknot_index_mults(bs.degree, &mults);

    let cknots: Vec<f64> = (first..=last).map(|k| at(&knots, k)).collect();
    let mut cmults: Vec<i32> = (first..=last).map(|k| ati(&mults, k)).collect();
    // OCCT L795: cmults(1) = cmults(Upper) = min(Deg, max(cmults(1), cmults(Upper))).
    let end_mult = (cmults[0].max(cmults[cmults.len() - 1])).min(bs.degree as i32);
    let cmults_upper = cmults.len() - 1;
    cmults[cmults_upper] = end_mult;
    cmults[0] = end_mult;

    // OCCT L799: compute new number of poles.
    let nbp = bspl_lib::nb_poles(bs.degree, true, &cmults);
    bs.control_points.truncate(nbp);
    if !bs.is_rational() {
        bs.weights = vec![1.0; nbp];
    } else {
        bs.weights.truncate(nbp);
    }
    bs.knots = flat_knots(&cknots, &cmults);
    bs.is_periodic = true;
}

// =========================================================================
// OCCT Geom_BSplineCurve.cxx L819-909 — SetOrigin(Index) (periodic only).
// =========================================================================
pub(crate) fn bspline_set_origin_index(bs: &mut BSplineCurve3, index: i32) {
    if !bs.is_periodic {
        panic!("Standard_NoSuchObject: Geom_BSplineCurve::SetOrigin");
    }
    let (knots, mults) = bs.knots_mults();
    let first = bspl_lib::first_uknot_index_mults(bs.degree, &mults);
    let last = bspl_lib::last_uknot_index_mults(bs.degree, &mults);
    if index < first || index > last {
        panic!("Standard_DomainError: Geom_BSplineCurve::SetOrigin");
    }
    let nbknots = knots.len() as i32;

    // OCCT L840-858: rotate the knots and mults.
    let mut newknots = Vec::with_capacity(nbknots as usize);
    let mut newmults = Vec::with_capacity(nbknots as usize);
    let period = at(&knots, last) - at(&knots, first);
    let mut i = index;
    while i <= last {
        newknots.push(at(&knots, i));
        newmults.push(ati(&mults, i));
        i += 1;
    }
    let mut i = first + 1;
    while i <= index {
        newknots.push(at(&knots, i) + period);
        newmults.push(ati(&mults, i));
        i += 1;
    }

    // OCCT L860-864: the pole index of the new origin.
    let mut pindex = 1i32;
    let mut i = first + 1;
    while i <= index {
        pindex += ati(&mults, i);
        i += 1;
    }

    // OCCT L867-884: rotate the poles (and weights).
    let nbpoles = bspline_nb_poles(bs) as i32;
    let mut newpoles = Vec::with_capacity(nbpoles as usize);
    let mut newweights = Vec::with_capacity(nbpoles as usize);
    let mut i = pindex;
    while i <= nbpoles {
        newpoles.push(bs.control_points[(i - 1) as usize]);
        newweights.push(bs.weights[(i - 1) as usize]);
        i += 1;
    }
    let mut i = 1i32;
    while i < pindex {
        newpoles.push(bs.control_points[(i - 1) as usize]);
        newweights.push(bs.weights[(i - 1) as usize]);
        i += 1;
    }
    bs.control_points = newpoles;
    bs.weights = newweights;
    bs.knots = flat_knots(&newknots, &newmults);
}

// =========================================================================
// OCCT Geom_BSplineCurve.cxx L913-970 — SetOrigin(U, Tol) (periodic only;
// called by ChFiDS_ElSpine::SetOrigin with Tol = Precision::PConfusion()).
// =========================================================================
pub(crate) fn bspline_set_origin_u_tol(bs: &mut BSplineCurve3, u: f64, tol: f64) {
    if !bs.is_periodic {
        panic!("Standard_NoSuchObject: Geom_BSplineCurve::SetOrigin");
    }
    // Is U within the period? (OCCT L919-929)
    let uf0 = bs.first_parameter();
    let ul0 = bs.last_parameter();
    let period = ul0 - uf0;
    let mut uu = u;
    while tol < uf0 - uu {
        uu += period;
    }
    while tol > ul0 - uu {
        uu -= period;
    }

    let mut uf = uf0;
    let mut ul = ul0;
    if (u - uu).abs() > tol {
        // Reparametrize the curve (OCCT L931-943).
        let delta = u - uu;
        uf += delta;
        ul += delta;
        for k in bs.knots.iter_mut() {
            *k += delta;
        }
    }
    // For periodic curve, uf and ul represent the same point (OCCT L944-948).
    if (u - uf).abs() < tol || (u - ul).abs() < tol {
        return;
    }

    // OCCT L950-960: locate the knot nearest to U.
    let nbknots = bspline_nb_knots(bs) as i32;
    let mut ik = 0i32;
    let mut delta = f64::MAX; // OCCT RealLast()
    for i in 1..=nbknots {
        let dki = bspline_knot(bs, i) - u;
        if dki.abs() < delta.abs() {
            ik = i;
            delta = dki;
        }
    }
    if delta.abs() > tol {
        // OCCT L963: InsertKnot(U) — the (M = 1, ParametricTolerance = 0,
        // Add = false) defaults.
        bspline_insert_knots(bs, &[u], &[1], 0.0);
        if delta < 0.0 {
            ik += 1;
        }
    }
    bspline_set_origin_index(bs, ik);
}

// =========================================================================
// OCCT Geom_BSplineCurve.cxx L974-1016 — SetNotPeriodic().
// GAP: the kernel BSplCLib::Unperiodize carries the non-rational 3D path
// (weights arrays pending the TKMath batch).
// =========================================================================
fn bspline_set_not_periodic(bs: &mut BSplineCurve3) {
    if bs.is_periodic {
        let (_, mults) = bs.knots_mults();
        let mut nb_knots = 0i32;
        let mut nb_poles = 0i32;
        bspl_lib::prepare_unperiodize(bs.degree, &mults, &mut nb_knots, &mut nb_poles);
        let poles = poles_flat(&bs.control_points);
        let (knots, mults) = bs.knots_mults();
        let mut new_mults = vec![0i32; nb_knots as usize];
        let mut new_knots = vec![0.0f64; nb_knots as usize];
        let mut new_poles = vec![0.0f64; 3 * nb_poles as usize];
        bspl_lib::unperiodize(
            bs.degree,
            &mults,
            &knots,
            &poles,
            &mut new_mults,
            &mut new_knots,
            &mut new_poles,
        );
        bs.control_points = unflatten_poles(&new_poles);
        bs.weights = vec![1.0; nb_poles as usize];
        bs.knots = flat_knots(&new_knots, &new_mults);
        bs.is_periodic = false;
    }
}

// =========================================================================
// OCCT Geom_BSplineCurve.cxx L527-715 — Segment(U1, U2, theTolerance).
// =========================================================================
pub(crate) fn bspline_segment(bs: &mut BSplineCurve3, u1: f64, u2: f64, the_tolerance: f64) {
    if u2 < u1 {
        panic!("Standard_DomainError: Geom_BSplineCurve::Segment");
    }

    let was_periodic = bs.is_periodic;
    let mut du = 0.0;
    let mut a_ddu = 0.0;
    if bs.is_periodic {
        // OCCT L545-558: define param distance to keep.
        let period = bs.last_parameter() - bs.first_parameter();
        du = u2 - u1;
        if du - period > rcad_kernel::core::precision::p_confusion() {
            panic!("Standard_DomainError: Geom_BSplineCurve::Segment");
        }
        if du > period {
            du = period;
        }
        a_ddu = du;
    }

    let mut new_u1 = 0.0;
    let mut new_u2 = 0.0;
    {
        let (knots, mults) = bs.knots_mults();
        let mut index = 0i32;
        let mut _u = 0.0f64;
        bspl_lib::locate_parameter_knots_mults(
            bs.degree, &knots, &mults, u1, bs.is_periodic, 1, knots.len() as i32,
            &mut index, &mut new_u1,
        );
        let mut index_b = 0i32;
        bspl_lib::locate_parameter_knots_mults(
            bs.degree, &knots, &mults, u2, bs.is_periodic, 1, knots.len() as i32,
            &mut index_b, &mut new_u2,
        );
    }

    let a_nu2 = new_u2;

    let seg_knots = [new_u1.min(new_u2), new_u1.max(new_u2)];
    let seg_mults = [bs.degree as i32, bs.degree as i32];
    let mut abs_umax = new_u1.abs().max(new_u2.abs());
    abs_umax = abs_umax.max(bs.first_parameter().abs().max(bs.last_parameter().abs()));
    let eps = epsilon(abs_umax).max(the_tolerance);

    bspline_insert_knots(bs, &seg_knots, &seg_mults, eps);

    if bs.is_periodic {
        // set the origine at NewU1 (OCCT L599-619).
        let (knots, mults) = bs.knots_mults();
        let mut index = 0i32;
        let mut u = 0.0f64;
        bspl_lib::locate_parameter_knots_mults(
            bs.degree, &knots, &mults, u1, bs.is_periodic, 1, knots.len() as i32,
            &mut index, &mut u,
        );
        if (at(&knots, index + 1) - u).abs() <= eps {
            index += 1;
        }
        bspline_set_origin_index(bs, index);
        bspline_set_not_periodic(bs);
        new_u2 = new_u1 + du;
    }

    // compute index1 and index2 to set the new knots and mults (OCCT L621-635).
    let (knots, mults) = bs.knots_mults();
    let mut index1 = 0i32;
    let mut index2 = 0i32;
    let mut u = 0.0f64;
    bspl_lib::locate_parameter_knots_mults(
        bs.degree, &knots, &mults, new_u1, bs.is_periodic, 1, knots.len() as i32,
        &mut index1, &mut u,
    );
    if (at(&knots, index1 + 1) - u).abs() <= eps {
        index1 += 1;
    }
    bspl_lib::locate_parameter_knots_mults(
        bs.degree, &knots, &mults, new_u2, bs.is_periodic, 1, knots.len() as i32,
        &mut index2, &mut u,
    );
    if (at(&knots, index2 + 1) - u).abs() <= eps || index2 == index1 {
        index2 += 1;
    }
    let nbknots = index2 - index1 + 1;

    // to restore changed U1 (OCCT L642-646).
    if du > 0.0 {
        du = new_u1 - u1;
    }

    let mut nknots = Vec::with_capacity(nbknots as usize);
    let mut nmults = Vec::with_capacity(nbknots as usize);
    let mut i = index1;
    while i <= index2 {
        nknots.push(at(&knots, i) - du);
        nmults.push(ati(&mults, i));
        i += 1;
    }
    nmults[0] = bs.degree as i32 + 1;
    nmults[(nbknots - 1) as usize] = bs.degree as i32 + 1;

    // compute index1 and index2 to set the new poles and weights (OCCT L658-665).
    let mut pindex1 = bspl_lib::pole_index(bs.degree, index1, bs.is_periodic, &mults);
    let mut pindex2 = bspl_lib::pole_index(bs.degree, index2, bs.is_periodic, &mults);
    pindex1 += 1;
    pindex2 = (pindex2 + 1).min(bspline_nb_poles(bs) as i32);
    let nbpoles = pindex2 - pindex1 + 1;

    let mut npoles = Vec::with_capacity(nbpoles as usize);
    let mut i = pindex1;
    while i <= pindex2 {
        npoles.push(bs.control_points[(i - 1) as usize]);
        i += 1;
    }

    // OCCT L690-698 (DBB).
    if was_periodic {
        nknots[0] = u1;
        if a_nu2 < u2 {
            nknots[(nbknots - 1) as usize] = u1 + a_ddu;
        }
    }

    bs.knots = flat_knots(&nknots, &nmults);
    bs.control_points = npoles;
    bs.weights = vec![1.0; nbpoles as usize];
    bs.is_periodic = false;
}

// =========================================================================
// OCCT TKGeomBase/GeomConvert/GeomConvert_CompCurveToBSplineCurve.cxx
// L32-274 — the curve concatenator (1:1 translation).
// =========================================================================
pub(crate) struct GeomConvertCompCurveToBSplineCurve {
    /// OCCT: occ::handle<Geom_BSplineCurve> myCurve.
    my_curve: Option<BSplineCurve3>,
    /// OCCT: double myTol.
    my_tol: f64,
    /// OCCT: Convert_ParameterisationType myType.
    my_type: ConvertParameterisation,
}

impl GeomConvertCompCurveToBSplineCurve {
    /// OCCT L32-37 — ctor(Parameterisation); myTol = Precision::Confusion().
    pub fn new(the_parameterisation: ConvertParameterisation) -> Self {
        GeomConvertCompCurveToBSplineCurve {
            my_curve: None,
            my_tol: rcad_kernel::core::precision::CONFUSION,
            my_type: the_parameterisation,
        }
    }

    /// OCCT L41-56 — ctor(BasisCurve, Parameterisation).
    pub fn with_basis_curve(basis_curve: &Curve3, parameterisation: ConvertParameterisation) -> Self {
        let mut out = Self::new(parameterisation);
        // OCCT: down_cast<Geom_BSplineCurve>(BasisCurve); non-null -> Copy(),
        // else GeomConvert::CurveToBSplineCurve(BasisCurve, myType).
        match basis_curve {
            Curve3::BSpline(bs) => {
                out.my_curve = Some(bs.clone());
            }
            _ => {
                out.my_curve = Some(geom_convert_curve_to_bspline_curve(basis_curve, out.my_type));
            }
        }
        out
    }

    /// OCCT L60-131 — Add(NewCurve, Tolerance, After = true, WithRatio =
    /// true, MinM = 2).
    pub fn add(&mut self, new_curve: &Curve3, tolerance: f64, after: bool) -> bool {
        self.add_full(new_curve, tolerance, after, true, 2)
    }

    /// OCCT L60-131 — Add with the explicit defaults.
    fn add_full(
        &mut self,
        new_curve: &Curve3,
        tolerance: f64,
        after: bool,
        with_ratio: bool,
        min_m: i32,
    ) -> bool {
        // conversion (OCCT L66-75).
        let mut bs: BSplineCurve3 = match new_curve {
            Curve3::BSpline(b) => b.clone(),
            _ => geom_convert_curve_to_bspline_curve(new_curve, self.my_type),
        };
        if self.my_curve.is_none() {
            self.my_curve = Some(bs);
            return true;
        }
        let my_curve = self.my_curve.as_mut().expect("myCurve");

        self.my_tol = tolerance;

        // Use actual curve endpoints instead of poles for proper G0
        // continuity check (OCCT L85-94).
        let a_curve_start = bspline_value(my_curve, my_curve.first_parameter());
        let a_curve_end = bspline_value(my_curve, my_curve.last_parameter());
        let a_bs_start = bspline_value(&bs, bs.first_parameter());
        let a_bs_end = bspline_value(&bs, bs.last_parameter());

        let mut avant = a_curve_start.distance(a_bs_start) < self.my_tol
            || a_curve_start.distance(a_bs_end) < self.my_tol;
        let mut apres = a_curve_end.distance(a_bs_start) < self.my_tol
            || a_curve_end.distance(a_bs_end) < self.my_tol;

        // Will myCurve be (or become) closed? (OCCT L96-107)
        if avant && apres {
            if after {
                avant = false;
            } else {
                apres = false;
            }
        }

        // Append after? (OCCT L109-118)
        if apres {
            if a_curve_end.distance(a_bs_end) < self.my_tol {
                bs = bs.reversed();
            }
            let mut first_curve = self.my_curve.take().expect("myCurve");
            self.my_curve = Some(add_sweep(
                &mut first_curve,
                bs,
                true,
                with_ratio,
                min_m,
                self.my_tol,
            ));
            return true;
        }
        // Prepend before? (OCCT L119-128)
        else if avant {
            if a_curve_start.distance(a_bs_start) < self.my_tol {
                bs = bs.reversed();
            }
            let mut first_curve = bs;
            let second_curve = self.my_curve.take().expect("myCurve");
            self.my_curve = Some(add_sweep(
                &mut first_curve,
                second_curve,
                false,
                with_ratio,
                min_m,
                self.my_tol,
            ));
            return true;
        }

        false
    }

    /// OCCT L264-267 — BSplineCurve().
    pub fn bspline_curve(&self) -> Option<BSplineCurve3> {
        self.my_curve.clone()
    }

    /// OCCT L271-274 — Clear().  Pending consumer: no translated ChFi3d
    /// code calls it yet (kept for the 1:1 class surface).
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.my_curve = None;
    }
}

/// OCCT GeomConvert_CompCurveToBSplineCurve.cxx L135-260 — the private
/// Add(FirstCurve, SecondCurve, After, WithRatio, MinM).  Returns the new
/// merged curve (the OCCT callee stores it into myCurve); FirstCurve is
/// increased in degree in place.
fn add_sweep(
    first_curve: &mut BSplineCurve3,
    mut second_curve: BSplineCurve3,
    after: bool,
    with_ratio: bool,
    min_m: i32,
    my_tol: f64,
) -> BSplineCurve3 {
    // Harmonize the degrees (OCCT L141-150).
    let deg = first_curve.degree.max(second_curve.degree);
    if first_curve.degree < deg {
        // GAP: the kernel IncreaseDegree carries the non-rational path only
        // (OCCT Geom_BSplineCurve.cxx L243-287; weights pending TKMath).
        first_curve.increase_degree(deg);
    }
    if second_curve.degree < deg {
        // GAP: rational weights — see above.
        second_curve.increase_degree(deg);
    }

    // OCCT L152-161: declarations (1-based OCCT arrays -> Vec + index - 1).
    let mut ratio = 1.0f64;
    let nb_p1 = bspline_nb_poles(first_curve) as i32;
    let nb_p2 = bspline_nb_poles(&second_curve) as i32;
    let nb_k1 = bspline_nb_knots(first_curve) as i32;
    let nb_k2 = bspline_nb_knots(&second_curve) as i32;
    let mut noeuds = vec![0.0f64; (nb_k1 + nb_k2 - 1) as usize];
    let mut poles = vec![DVec3::ZERO; (nb_p1 + nb_p2 - 1) as usize];
    let mut poids = vec![0.0f64; (nb_p1 + nb_p2 - 1) as usize];
    let mut mults = vec![0i32; (nb_k1 + nb_k2 - 1) as usize];

    // Reparameterization ratio (C1 if possible) (OCCT L163-177).
    if with_ratio {
        let l1 = bspline_d1(first_curve, first_curve.last_parameter()).length();
        let l2 = bspline_d1(&second_curve, second_curve.first_parameter()).length();
        if l1 > rcad_kernel::core::precision::CONFUSION
            && l2 > rcad_kernel::core::precision::CONFUSION
        {
            ratio = l1 / l2;
        }
        if ratio < rcad_kernel::core::precision::CONFUSION
            || ratio > 1.0 / rcad_kernel::core::precision::CONFUSION
        {
            ratio = 1.0;
        }
    }

    let (ratio1, delta1, ratio2, delta2);
    if after {
        // Do not move the first curve (OCCT L179-186).
        ratio1 = 1.0;
        delta1 = 0.0;
        ratio2 = 1.0 / ratio;
        delta2 = ratio2 * bspline_knot(&second_curve, 1) - bspline_knot(first_curve, nb_k1);
    } else {
        // Do not move the second curve (OCCT L187-194).
        ratio1 = ratio;
        delta1 = ratio1 * bspline_knot(first_curve, nb_k1) - bspline_knot(&second_curve, 1);
        ratio2 = 1.0;
        delta2 = 0.0;
    }

    // The Knots (OCCT L196-229).
    for ii in 1..=nb_k1 {
        let v = ratio1 * bspline_knot(first_curve, ii) - delta1;
        noeuds[(ii - 1) as usize] = v;
        if ii > 1 {
            let mut eps = epsilon(noeuds[(ii - 2) as usize].abs());
            if eps < 5.0e-10 {
                eps = 5.0e-10;
            }
            if noeuds[(ii - 1) as usize] - noeuds[(ii - 2) as usize] <= eps {
                noeuds[(ii - 1) as usize] += eps;
            }
        }
        mults[(ii - 1) as usize] = bspline_multiplicity(first_curve, ii);
    }
    mults[(nb_k1 - 1) as usize] = first_curve.degree as i32;
    let mut jj = nb_k1 + 1;
    for ii in 2..=nb_k2 {
        let v = ratio2 * bspline_knot(&second_curve, ii) - delta2;
        noeuds[(jj - 1) as usize] = v;
        let mut eps = epsilon(noeuds[(jj - 2) as usize].abs());
        if eps < 5.0e-10 {
            eps = 5.0e-10;
        }
        if noeuds[(jj - 1) as usize] - noeuds[(jj - 2) as usize] <= eps {
            noeuds[(jj - 1) as usize] += eps;
        }
        mults[(jj - 1) as usize] = bspline_multiplicity(&second_curve, ii);
        jj += 1;
    }

    // OCCT L231-232: Ratio = Weight(NbP1) / Weight(1) of the second curve.
    ratio = first_curve.weights[(nb_p1 - 1) as usize] / second_curve.weights[0];
    // The Poles and Weights (OCCT L233-247).
    for ii in 1..nb_p1 {
        poles[(ii - 1) as usize] = bspline_pole(first_curve, ii);
        poids[(ii - 1) as usize] = first_curve.weights[(ii - 1) as usize];
    }
    let mut jj = nb_p1;
    for ii in 1..=nb_p2 {
        poles[(jj - 1) as usize] = bspline_pole(&second_curve, ii);
        poids[(jj - 1) as usize] = ratio * second_curve.weights[(ii - 1) as usize];
        jj += 1;
    }

    // Create the BSpline (OCCT L249-250).
    let mut my_curve = BSplineCurve3 {
        degree: deg,
        knots: flat_knots(&noeuds, &mults),
        control_points: poles,
        weights: poids,
        is_periodic: false,
    };

    // Optionally reduce multiplicity down to MinM (OCCT L252-259).
    let mut ok = true;
    let mut m = mults[(nb_k1 - 1) as usize];
    while m > min_m && ok {
        m -= 1;
        ok = bspline_remove_knot(&mut my_curve, nb_k1, m, my_tol);
    }
    my_curve
}

// =========================================================================
// OCCT TKGeomBase/GeomLib/GeomLib.cxx L126-253 — the static ComputeLambda
// (Constraint, Hermit, Length, Lambda): the Gauss-quadrature minimization of
// the integrated square deviation |C'(t)|^2 - |C'(0)|^2 over [0, 1],
// followed by the extremum refinement.
// GAP: the extremum refinement runs OCCT math_FunctionAllRoots over a
// GeomLib_LogSample (log-spaced parameters); the kernel math_FunctionAllRoots
// binds the concrete linear math_FunctionSample (no virtual GetParameter),
// so the refinement samples the same interval linearly — pending the TKMath
// batch.  The OCCT fallback (keep Lambda when no better extremum is found)
// is preserved.
// =========================================================================
// =========================================================================
// OCCT TKGeomBase/GeomConvert/GeomConvert.cxx L1313-1339 —
// C0BSplineToC1BSplineCurve(BS, tolerance, AngularTol).
// GAP: pending the TKGeomBase batch (C0BSplineToArrayOfC1BSplineCurve +
// ConcatC1).  The carrier keeps the input curve; the caller's
// FirstParameter identity check (Builder_0.cxx L5205) then keeps the
// input, preserving the OCCT failure fallback (L5209-5212).
// =========================================================================
pub(crate) fn geom_convert_c0_bspline_to_c1_bspline_curve(
    _bs: &mut BSplineCurve3,
    _tolerance: f64,
    _angular_tol: f64,
) {
}

// =========================================================================
// OCCT TKTopAlgo/BRepLib/BRepLib.cxx L301-456 — BuildCurve3d(AnEdge,
// Tolerance, Continuity, MaxDegree, MaxSegment).  "if the edge has a 3d
// curve returns true" (L319-325).  The pcurve reconstruction branches
// (L330-454) route through GeomLib::BuildCurve3d — the GeomLib body has
// landed (`geomalgo::geom_lib::build_curve3d`); this BRepLib-level carrier is
// still a stand-in and returns false when the 3d curve cannot be produced
// (L363-366 / L436-439 / L452).
// =========================================================================
pub(crate) fn brep_lib_build_curve3d(an_edge: &mut Shape) -> bool {
    an_edge
        .as_edge()
        .map_or(false, |ed| ed.curve.is_some())
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L3659-3681 — CurveCleaner(BS, Tol, MultMin):
// makes a BSpline as much continued as possible at a given tolerance.
// =========================================================================
pub(crate) fn curve_cleaner(bs: &mut BSplineCurve3, tol: f64, mult_min: i32) {
    let mut tol = tol;
    let nb_k = bspline_nb_knots(bs) as i32;
    let mut mult = bs.degree as i32;
    while mult > mult_min {
        tol *= 0.5; // Progressive reduction
        let mut ii = nb_k;
        while ii > 1 {
            if bspline_multiplicity(bs, ii) == mult {
                bspline_remove_knot(bs, ii, mult - 1, tol);
            }
            ii -= 1;
        }
        mult -= 1;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L4773-4793 — GoodExt(C, V, f, l, a).
// =========================================================================
fn good_ext(c: &Curve3, v: DVec3, f: f64, l: f64, a: f64) -> bool {
    for i in 0..6i32 {
        let t = i as f64 * 0.2;
        let par = (1.0 - t) * f + t * l;
        // OCCT: C->D1(par, d0, d1); d1.Angle(V).
        let d1 = c.derivative_at(par);
        let ang = gp_vec_angle(d1, v);
        let angref = a * t + 0.002;
        if ang > angref {
            return false;
        }
    }
    true
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5835-5851 — ChFi3d_ApproxByC2(C).
// GAP: GeomAPI_PointsToBSpline(Points, Approx_ChordLength, 3, 8, GeomAbs_C2,
// 1.000001e-3) (TKGeomBase/GeomAPI) is pending; the kernel
// GeomAPI_Interpolate carrier (exact pass through the sample points,
// chord-length parameterization, degree 3) stands in.
// =========================================================================
pub(crate) fn chfi3d_approx_by_c2(c: &Curve3) -> BSplineCurve3 {
    let first = c.default_domain()[0];
    let last = c.default_domain()[1];
    let nb_points = 101usize;

    let mut points: Vec<DVec3> = Vec::with_capacity(nb_points);
    let delta = (last - first) / (nb_points as f64 - 1.0);
    for i in 1..=nb_points - 1 {
        points.push(c.point_at(first + (i as f64 - 1.0) * delta));
    }
    points.push(c.point_at(last));

    match rcad_kernel::base::geom_api::interpolate::interpolate_points(&points) {
        Ok(bs) => bs,
        Err(_) => panic!("GAP carrier: GeomAPI_PointsToBSpline (ChFi3d_Builder_0.cxx L5848)"),
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5855-5944 — ChFi3d_IsSmooth(C).
// GAP: GeomLProp_CLProps(C, 2, Resolution) — the rcad CLProps is the 2D
// variant (base/geom_lprop/cl_props_base.rs); the 3D tangent/curvature/
// centre queries below use the OCCT GeomLProp_CLProps.cxx formulas on the
// curve D1/D2 (IsTangentDefined <=> |D1| nonzero; Curvature = |D1^D2|/|D1|^3;
// CentreOfCurvature = P + normal/curvature).
// =========================================================================
pub(crate) fn chfi3d_is_smooth(c: &Curve3) -> bool {
    let nbintv = c.nb_intervals(GeomAbsShape::CN);
    let mut ti: Vec<f64> = Vec::new();
    c.intervals(&mut ti, GeomAbsShape::CN);
    let resolution = bspl_lib::GP_RESOLUTION; // OCCT gp::Resolution()

    let is_tangent_defined = |t: f64| -> bool { c.derivative_at(t).length() > resolution };
    let curvature = |t: f64| -> f64 {
        let d1 = c.derivative_at(t);
        let d2 = c.derivative2_at(t);
        let n = d1.cross(d2).length();
        let m3 = d1.length().powi(3);
        if m3 <= 0.0 {
            0.0
        } else {
            n / m3
        }
    };
    let centre_of_curvature = |t: f64| -> DVec3 {
        let p = c.point_at(t);
        let d1 = c.derivative_at(t);
        let d2 = c.derivative2_at(t);
        let n = d1.cross(d2);
        let k = curvature(t);
        if k > resolution && n.length() > 0.0 {
            p + n.normalize() / k
        } else {
            p
        }
    };

    let mut prev_vec = DVec3::ZERO;
    let mut prev_vec_found = false;
    let mut intrv_found = 0usize;
    for intrv in 1..=nbintv {
        let mut t = ti[intrv - 1];
        let step = (ti[intrv] - t) / 30.0; // Discretisation = 30
        for _ii in 1..=30i32 {
            if !is_tangent_defined(t) {
                return false;
            }
            let curv = curvature(t).abs();
            if curv > resolution {
                let p1 = c.point_at(t);
                let p2 = centre_of_curvature(t);
                prev_vec = p2 - p1; // OCCT: gp_Vec(P1, P2)
                prev_vec_found = true;
                break;
            }
            t += step;
        }
        if prev_vec_found {
            intrv_found = intrv;
            break;
        }
    }

    if !prev_vec_found {
        return true;
    }

    for intrv in intrv_found..=nbintv {
        let mut t = ti[intrv - 1];
        let step = (ti[intrv] - t) / 30.0;
        for ii in 1..=30i32 {
            if !is_tangent_defined(t) {
                return false;
            }
            let curv = curvature(t).abs();
            if curv > resolution {
                let p1 = c.point_at(t);
                let p2 = centre_of_curvature(t);
                let vec = p2 - p1;
                let angle = gp_vec_angle(prev_vec, vec);
                if angle > std::f64::consts::PI / 3.0 {
                    return false;
                }
                let mut ratio = vec.length() / prev_vec.length();
                if ratio < 1.0 {
                    ratio = 1.0 / ratio;
                }
                if ratio > 2.0 && (intrv != nbintv || ii != 30) {
                    return false;
                }
                prev_vec = vec;
            }
            t += step;
        }
    }

    true
}

// =========================================================================
// OCCT Geom_Curve::Reversed — R(t) = C(ReversedParameter(t)).
// Geom_Line::Reverse negates the direction; Geom_Conic::Reverse mirrors the
// conic frame (the rcad eval derives the minor axis from normal x
// major_dir, so flipping `normal` mirrors the parameterization over the
// same point set); Geom_BSplineCurve::Reverse L496-516 (kernel
// BSplineCurve3::reversed); Geom_TrimmedCurve::Reverse swaps first/last
// through the basis ReversedParameter.
// =========================================================================
fn geom_curve_reversed(cv: &Curve3) -> Curve3 {
    match cv {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin,
            direction: -l.direction,
        }),
        Curve3::Circle(ci) => Curve3::Circle(Circle3 {
            y_dir: -ci.y_dir,
            ..ci.clone()
        }),
        Curve3::Ellipse(el) => Curve3::Ellipse(Ellipse3 {
            normal: -el.normal,
            ..el.clone()
        }),
        Curve3::BSpline(bs) => Curve3::BSpline(bs.reversed()),
        Curve3::Trimmed(ct) => {
            let basis = geom_curve_reversed(ct.basis_curve());
            let first = basis.reversed_parameter(ct.last);
            let last = basis.reversed_parameter(ct.first);
            Curve3::Trimmed(TrimmedCurve3::new(basis, first, last))
        }
        _ => panic!(
            "Geom_Curve::Reversed: curve variant not carried by the rcad edge geometry"
        ),
    }
}

// =========================================================================
// OCCT BRep_Tool::Curve(E, First, Last) — the geometry curve and range
// (orientation-independent).
// =========================================================================
fn brep_tool_curve(e: &Shape) -> (Curve3, f64, f64) {
    let ed = e.as_edge().expect("not an edge");
    (
        ed.curve.clone().expect("edge curve"),
        ed.range[0],
        ed.range[1],
    )
}

// OCCT BRep_Tool::Tolerance(V).
fn brep_tool_tolerance(v: &Shape) -> f64 {
    v.as_vertex().map(|vd| vd.tolerance).unwrap_or(0.0)
}

// OCCT BRepAdaptor_Curve::Resolution(R3d) -> GeomAdaptor_Curve::Resolution.
fn brep_adaptor_curve_resolution(e: &Shape, r3d: f64) -> f64 {
    let ed = e.as_edge().expect("not an edge");
    ed.curve.as_ref().map(|c| c.resolution(r3d)).unwrap_or(r3d)
}

// OCCT gp_Pnt::IsEqual(theOther, theLinTol) — Distance <= theLinTol.
fn gp_pnt_is_equal(p: DVec3, other: DVec3, tol: f64) -> bool {
    p.distance(other) <= tol
}

// OCCT gp_Vec::Angle(theOther) — the angle in [0, PI] between the vectors
// (acos of the normalized dot product; OCCT raises on null vectors, the
// carrier reports PI through the clamped zero dot).
fn gp_vec_angle(a: DVec3, b: DVec3) -> f64 {
    let dot = a.normalize_or_zero().dot(b.normalize_or_zero());
    dot.clamp(-1.0, 1.0).acos()
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L4797-5327 — ChFi3d_PerformElSpine(HES, Spine,
// continuity, tol, IsOffset).
// =========================================================================
pub fn chfi3d_perform_elspine(
    hes: &mut ChFiDSElSpine,
    spine: &mut ChFiDSSpineHandle,
    continuity: ChFiDSGeomAbsShape,
    tol: f64,
    is_offset: bool,
) {
    let periodic;
    let mut bof;
    let checkdeb;
    let mut cepadur;
    let mut i_edge;
    let if_;
    let mut il;
    let nbed;
    let mut i_to_approx_by_c2;
    let mut wf;
    let mut wl;
    let mut wrefdeb;
    let mut wreffin;
    let mut nwf;
    let mut nwl;
    let mut pared = 0.0;
    let mut first;
    let mut last;
    let mut eps_v;
    let mut a_continuity;
    let p_deb;
    let p_fin;
    let mut vref_deb;
    let mut vref_fin;
    let mut cv;
    let mut tc: Curve3;
    let mut bs;
    let mut bspline;
    let mut e;
    let mut eold;
    let mut v = Shape::null();

    // OCCT L4817-4825: setup.
    wf = hes.first_parameter();
    wl = hes.last_parameter();
    wrefdeb = wf;
    wreffin = wl;
    nwf = wf;
    nwl = wl;
    nbed = spine.base().nb_edges();
    periodic = spine.base().is_periodic();
    if periodic {
        let period = spine.base().period();
        nwf = elclib::in_period(nwf, -tol, period - tol);
        if_ = spine.base().index_of_param(nwf, true);
        nwl = elclib::in_period(nwl, tol, period + tol);
        il = spine.base().index_of_param(nwl, false);
        if nwl < nwf + tol {
            il += nbed;
        }
    } else {
        if_ = spine.base().index_of_param(wf, true);
        il = spine.base().index_of_param(wl, false);
        wrefdeb = spine.base().first_parameter_of(if_).max(wf);
        wreffin = spine.base().last_parameter_of(il).min(wl);
    }
    // OCCT L4846-4849.
    let d1_deb = spine.base_mut().d1(wf);
    p_deb = d1_deb.0;
    vref_deb = d1_deb.1;
    let d1_fin = spine.base_mut().d1(wl);
    p_fin = d1_fin.0;
    vref_fin = d1_fin.1;
    vref_deb = vref_deb.normalize();
    vref_fin = vref_fin.normalize();

    // OCCT L4851-4853: ExtrapPole/ExtraCoeffs/Cont (1,5) arrays — declared
    // by OCCT for the extremity treatment; the GeomLib carriers carry their
    // own scratch.

    // Attention on segmente eventuellement la premiere et la derniere arete.
    // Traitment de la premiere arete (OCCT L4857-4862).
    cepadur = false;
    e = if is_offset {
        spine.base().offset_edges(if_).clone()
    } else {
        spine.base().edges(if_).clone()
    };
    bof = brep_lib_build_curve3d(&mut e);
    let _ = bof;
    let tolpared = brep_adaptor_curve_resolution(spine.base().edges(if_), tol);
    let cv0 = brep_tool_curve(&e);
    cv = cv0.0;
    first = cv0.1;
    last = cv0.2;
    // Add vertex with tangent (OCCT L4864-4871).
    if hes.is_periodic() {
        let par_for_elspine = if e.orientation == Orientation::Forward {
            first
        } else {
            last
        };
        let pnt_for_elspine = cv.point_at(par_for_elspine);
        let dir_for_elspine = cv.derivative_at(par_for_elspine);
        hes.add_vertex_with_tangent(Ax1::new(pnt_for_elspine, dir_for_elspine));
    }

    // OCCT L4873-4878.
    let urefdeb = spine.base().first_parameter_of(if_);
    checkdeb = nwf > urefdeb;
    if checkdeb {
        pared = spine.base_mut().parameter_on(nwf, if_, false);
    }

    // OCCT L4880-4944: the first-edge trimming (REVERSED edge vs #1).
    if e.orientation == Orientation::Reversed {
        let sov = first;
        first = cv.reversed_parameter(last);
        last = cv.reversed_parameter(sov);
        if checkdeb {
            pared = cv.reversed_parameter(pared);
        } else {
            pared = first;
        }
        if first < pared {
            first = pared;
        }
        if il == if_ {
            let ureffin = spine.base().last_parameter_of(il);
            let checkfin = nwl < ureffin;
            if checkfin {
                pared = spine.base_mut().parameter_on(nwl, il, false);
                pared = cv.reversed_parameter(pared);
            } else {
                pared = last;
            }
            if pared < last {
                last = pared;
            }
        }
        cv = geom_curve_reversed(&cv);
    } else {
        if !checkdeb {
            pared = first;
        }
        if first < pared {
            first = pared;
        }
        if il == if_ {
            let ureffin = spine.base().last_parameter_of(il);
            let checkfin = nwl < ureffin;
            if checkfin {
                pared = spine.base_mut().parameter_on(nwl, il, false);
            } else {
                pared = last;
            }
            if pared < last {
                last = pared;
            }
        }
    }

    // OCCT L4946-4949.
    if (last - first).abs() < tolpared {
        cepadur = true;
    }

    // Petite veru pour les cas ou un KPart a bouffe l arete sans parvenir a
    // terminer. On tire une droite. (OCCT L4951-4977)
    if cepadur {
        if wl < spine.base().first_parameter_of(1) + tol {
            let (_ptemp, vtemp) = hes.last_point_and_tgt();
            let d = vtemp.normalize();
            // OCCT: olin.ChangeCoord().SetLinearForm(-WL, d.XYZ(), PFin.XYZ()).
            let olin = p_fin - d * wl;
            let l = Curve3::Line(Line3::new(olin, d));
            hes.set_curve(l);
        } else if wf > spine.base().last_parameter_of(nbed) - tol {
            let (_ptemp, vtemp) = hes.first_point_and_tgt();
            let d = vtemp.normalize();
            let olin = p_deb - d * wf;
            let l = Curve3::Line(Line3::new(olin, d));
            hes.set_curve(l);
        }
        return; // =>
    }

    // OCCT L4979-4992.
    tc = Curve3::Trimmed(TrimmedCurve3::new(cv, first, last));
    bs = geom_convert_curve_to_bspline_curve(&tc, ConvertParameterisation::TgtThetaOver2);
    curve_cleaner(&mut bs, (wl - wf).abs() * 1.0e-4, 0);

    // Smoothing of the curve (OCCT L4984-4992).
    i_to_approx_by_c2 = 0;
    a_continuity = tc.continuity();
    let b_is_smooth = chfi3d_is_smooth(&tc);
    if a_continuity < GeomAbsShape::C2 && !b_is_smooth {
        i_to_approx_by_c2 += 1;
        bs = chfi3d_approx_by_c2(&tc);
        tc = Curve3::BSpline(bs.clone());
    }

    // Concatenation des aretes suivantes (OCCT L4995).
    let mut concat =
        GeomConvertCompCurveToBSplineCurve::with_basis_curve(&tc, ConvertParameterisation::QuasiAngular);

    eold = e.clone();
    i_edge = if_ + 1;
    while i_edge <= il {
        // OCCT L5000-5004: the periodic spine indexing.
        let mut iloc = i_edge;
        if periodic {
            iloc = (i_edge - 1) % nbed + 1;
        }
        // OCCT L5006-5010.
        e = if is_offset {
            spine.base().offset_edges(iloc).clone()
        } else {
            spine.base().edges(iloc).clone()
        };
        if e.as_edge().map(|ed| ed.degenerated).unwrap_or(false) {
            i_edge += 1;
            continue;
        }

        // OCCT L5012-5017.
        eps_v = tol;
        let v_opt = topexp_common_vertex(&eold, &e);
        bof = v_opt.is_some();
        if let Some(vv) = v_opt {
            v = vv;
            eps_v = brep_tool_tolerance(&v);
        }

        // OCCT L5019-5023.
        bof = brep_lib_build_curve3d(&mut e);
        if !bof {
            panic!("Standard_ConstructionError: PerformElSpine : BuildCurve3d error");
        }

        // OCCT L5025-5031.
        let cv1 = brep_tool_curve(&e);
        cv = cv1.0;
        first = cv1.1;
        last = cv1.2;
        // Add vertex with tangent.
        let par_for_elspine = if e.orientation == Orientation::Forward {
            first
        } else {
            last
        };
        let pnt_for_elspine = cv.point_at(par_for_elspine);
        let dir_for_elspine = cv.derivative_at(par_for_elspine);
        hes.add_vertex_with_tangent(Ax1::new(pnt_for_elspine, dir_for_elspine));

        // OCCT L5033-5064: the last-edge trimming.
        if i_edge == il {
            let ureffin = spine.base().last_parameter_of(iloc);
            let checkfin = nwl < ureffin;
            if checkfin {
                pared = spine.base_mut().parameter_on(nwl, iloc, false);
            } else {
                pared = last;
            }
            if e.orientation == Orientation::Reversed {
                let sov = first;
                first = cv.reversed_parameter(last);
                last = cv.reversed_parameter(sov);
                if checkfin {
                    pared = cv.reversed_parameter(pared);
                } else {
                    pared = last;
                }
                cv = geom_curve_reversed(&cv);
            }
            if pared < last {
                last = pared;
            }
        }

        // OCCT L5066-5078.
        tc = Curve3::Trimmed(TrimmedCurve3::new(cv, first, last));
        bs = geom_convert_curve_to_bspline_curve(&tc, ConvertParameterisation::TgtThetaOver2);
        curve_cleaner(&mut bs, (wl - wf).abs() * 1.0e-4, 0);

        // Smoothing of the curve.
        a_continuity = tc.continuity();
        let b_is_smooth = chfi3d_is_smooth(&tc);
        if a_continuity < GeomAbsShape::C2 && !b_is_smooth {
            i_to_approx_by_c2 += 1;
            bs = chfi3d_approx_by_c2(&tc);
            tc = Curve3::BSpline(bs.clone());
        }

        // OCCT L5080-5094.
        let tolrac = tol.min(eps_v);
        bof = concat.add(&tc, 2.0 * tolrac, true);
        // si l'ajout ne s'est pas bien passe on essai d'augmenter la tolerance
        if !bof {
            bof = concat.add(&tc, 2.0 * eps_v, true);
        }
        if !bof {
            bof = concat.add(&tc, 200.0 * eps_v, true);
            if !bof {
                panic!("Standard_ConstructionError: PerformElSpine: spine merged error");
            }
        }
        eold = e.clone();
        i_edge += 1;
    } // for (IEdge=IF+1; IEdge<=IL; ++IEdge)

    // On a la portion d elspine calculee sans prolongements sur la partie
    // valide des aretes du chemin. (OCCT L5098-5104)
    let mut bspline_ = concat
        .bspline_curve()
        .expect("PerformElSpine: null concatenated curve");
    // There is a reparametrisation to maximally connect the abscissas of edges.
    let (mut bs_noeuds, bs_mults) = bspline_.knots_mults();
    bspl_lib::reparametrize(wrefdeb, wreffin, &mut bs_noeuds);
    bspline_.set_knots(&bs_noeuds, &bs_mults);
    bspline = bspline_;

    // Traitement des Extremites (OCCT L5106-5180).
    let mut caredeb;
    let mut carefin;
    let angle = std::f64::consts::PI * 0.75;
    let mut local_wl = wl;
    let mut local_wf = wf;
    caredeb = 0;
    carefin = 0;
    if !hes.is_periodic() && !gp_pnt_is_equal(p_deb, bspline_pole(&bspline, 1), tol) {
        // Prolongement C3 au debut : afin d'eviter des pts d'inflexions dans
        // la partie utile de la spine le prolongement se fait jusqu'a un
        // point eloigne. (OCCT L5119-5149)
        if bspline.is_rational() {
            caredeb = 1;
        }
        let rabdist = wrefdeb - wf;
        let bout = p_deb - vref_deb * (20.0 * rabdist);
        let mut goodext = false;
        let mut newc = bspline.clone();
        let mut icont = 3i32;
        while icont >= 1 && !goodext {
            let mut an_ext_curve = bspline.clone();
            geom_lib_extend_curve_to_point(&mut an_ext_curve, bout, icont, false);
            newc = an_ext_curve;
            // OCCT: gacurve.Load(newc); GCPnts_AbscissaPoint GCP(gacurve,
            // -rabdist, Wrefdeb, WF); — the kernel helper merges IsDone
            // (the seed is returned when the iteration does not converge).
            wf = abscissa_point_parameter(&Curve3::BSpline(newc.clone()), wrefdeb, wf, -rabdist, wrefdeb);
            goodext = good_ext(&Curve3::BSpline(newc.clone()), vref_deb, wrefdeb, wf, angle);
            icont -= 1;
        }
        if caredeb != 0 {
            caredeb = bspline_nb_knots(&newc) as i32 - bspline_nb_knots(&bspline) as i32;
        }
        bspline = newc;
        local_wf = bspline.first_parameter();
    }

    if !hes.is_periodic() && !gp_pnt_is_equal(p_fin, bspline_pole(&bspline, bspline_nb_poles(&bspline) as i32), tol) {
        // Prolongement C3 en fin (OCCT L5153-5180).
        if bspline.is_rational() {
            carefin = 1;
        }
        let rabdist = wl - wreffin;
        let bout = p_fin + vref_fin * (20.0 * rabdist);
        let mut goodext = false;
        let mut newc = bspline.clone();
        let mut icont = 3i32;
        while icont >= 1 && !goodext {
            let mut an_ext_curve = bspline.clone();
            geom_lib_extend_curve_to_point(&mut an_ext_curve, bout, icont, true);
            newc = an_ext_curve;
            let wl_gcp =
                abscissa_point_parameter(&Curve3::BSpline(newc.clone()), wreffin, wl, rabdist, wreffin);
            wl = wl_gcp;
            goodext = good_ext(&Curve3::BSpline(newc.clone()), vref_fin, wreffin, wl, angle);
            icont -= 1;
        }
        if carefin != 0 {
            carefin = bspline_nb_knots(&newc) as i32 - bspline_nb_knots(&bspline) as i32;
        }
        bspline = newc;
        local_wl = bspline.last_parameter();
    }

    // Reparametrisation et segmentation sur le domaine de la Spine
    // (OCCT L5182-5197).
    if (bspline.first_parameter() - wf).abs() < tol {
        wf = bspline.first_parameter();
    }
    if (bspline.last_parameter() - wl).abs() < tol {
        wl = bspline.last_parameter();
    }

    if local_wf < wf || local_wl > wl {
        // pour eviter des pb avec segment!
        bspline_segment(&mut bspline, wf, wl, 0.0);
        hes.set_first_parameter(wf);
        hes.set_last_parameter(wl);
    }

    // OCCT L5199-5213: the C0 -> C1 split on rational curves.
    if bspline.is_rational() {
        let mut c1 = bspline.clone();
        geom_convert_c0_bspline_to_c1_bspline_curve(&mut c1, tol, 0.1);
        // Il faut s'assurer que l'origine n'a pas bouge (cts21158).
        if c1.first_parameter() == bspline.first_parameter() {
            bspline = c1;
        }
    }

    // OCCT L5215-5221: deformation eventuelle pour rendre la spine C2
    // (ou C3 pour des approx C2).
    if (caredeb != 0 || carefin != 0) && bspline.degree < 8 {
        // GAP: the kernel IncreaseDegree carries the non-rational path only
        // (OCCT Geom_BSplineCurve.cxx L243-287; weights pending TKMath).
        bspline.increase_degree(8);
    }

    let mut fk = 2i32;
    let mut lk = bspline_nb_knots(&bspline) as i32 - 1;
    if bspline.is_periodic {
        fk = 1;
    }
    if caredeb != 0 {
        fk += caredeb;
    }
    if carefin != 0 {
        lk -= carefin;
    }

    let mult_max;
    if continuity == ChFiDSGeomAbsShape::C3 {
        if bspline.degree < 7 {
            bspline.increase_degree(7);
        }
        mult_max = bspline.degree as i32 - 3;
    } else {
        if bspline.degree < 5 {
            bspline.increase_degree(5);
        }
        mult_max = bspline.degree as i32 - 2;
    }
    // correction C2 or C3 (if possible) (OCCT L5254-5256).
    curve_cleaner(&mut bspline, (wl - wf).abs() * 1.0e-4, 1);
    curve_cleaner(&mut bspline, (wl - wf).abs() * 1.0e-2, mult_max);
    let mult_min = (bspline.degree as i32 - 4).max(1);
    let mut ii = fk;
    while ii <= lk {
        if bspline_multiplicity(&bspline, ii) > mult_max {
            bof = bspline_remove_knot(&mut bspline, ii, mult_max, (wl - wf).abs() / 10.0);
        }
        // See C4 (OCCT L5265-5268).
        if bspline_multiplicity(&bspline, ii) > mult_min {
            bof = bspline_remove_knot(&mut bspline, ii, mult_min, (wl - wf).abs() * 1.0e-4);
        }
        ii += 1;
    }

    // elspine periodic => BSpline Periodic (OCCT L5270-5284).
    if hes.is_periodic() {
        if !bspline.is_periodic {
            bspline_set_periodic(&mut bspline);
            // modified by NIZNHY-PKV Fri Dec 10 12:20:22 2010ft
            if i_to_approx_by_c2 != 0 {
                bof = bspline_remove_knot(&mut bspline, 1, mult_max, (wl - wf).abs() / 10.0);
            }
        }
    } else {
        // Otherwise is it necessary to move the poles to adapt them to new
        // tangents ? (OCCT L5285-5318)
        let mut adjust = false;
        let p1 = bspline_value(&bspline, wf);
        let mut v1 = bspline_d1(&bspline, wf);
        v1 = v1.normalize();
        let (pdeb_t, vrefdeb_t) = hes.first_point_and_tgt();
        let scaldeb = vrefdeb_t.dot(v1);
        let disdeb = pdeb_t.distance(p1);
        if (wf - local_wf).abs() < 1.0e-12 && (scaldeb <= 0.9999999 || disdeb >= tol) {
            // Yes if there was no extension and the tangent is not the good one.
            adjust = true;
        }
        let p2 = bspline_value(&bspline, wl);
        let mut v2 = bspline_d1(&bspline, wl);
        v2 = v2.normalize();
        let (pfin_t, vreffin_t) = hes.last_point_and_tgt();
        let scalfin = vreffin_t.dot(v2);
        let disfin = pfin_t.distance(p2);
        if (wl - local_wl).abs() < 1.0e-12 && (scalfin <= 0.9999999 || disfin >= tol) {
            // the same at the end
            adjust = true;
        }
        if adjust {
            geom_lib_adjust_extremity(&mut bspline, pdeb_t, pfin_t, vrefdeb_t, vreffin_t);
        }
    }

    // Le Resultat (OCCT L5320-5321).
    hes.set_curve(Curve3::BSpline(bspline));
}
