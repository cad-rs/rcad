//! OCCT GeomConvert_CompCurveToBSplineCurve (TKGeomBase/GeomConvert) — the
//! incremental pole-by-pole curve concatenator (3D twin of
//! `base::geom2d_convert::comp_curve_to_bspline_2d`).
//!
//! Sources (1:1 statement-mapped):
//! - `GeomConvert_CompCurveToBSplineCurve.cxx` L32-274 (+ .hxx L30-78)
//! - `GeomConvert.cxx` L163-380 (CurveToBSplineCurve — the trimmed-BSpline
//!   clamp/segment branch lives here; the trimmed line/circle branches are
//!   re-used from the existing `geom_convert_curve_to_bspline_curve`)
//! - `Geom_BSplineCurve.cxx` L243-287 (IncreaseDegree), L420-492
//!   (RemoveKnot), L527-715 (Segment) — the curve-level operations the
//!   concatenator mutates.
//!
//! Architecture note (same adaptation as `geom::bspline_ops`): the rcad
//! legacy `BSplineCurve3` carries the flat (expanded) knot vector instead of
//! OCCT's (knots, mults) pairs; the two are bijective and the conversions
//! here preserve OCCT semantics exactly.  Rational curves are homogenized to
//! the flat dim-4 kernel form (PLib::SetPoles) exactly as the OCCT
//! gp_Pnt overload does.  The legacy struct keeps a periodic flag but the
//! concatenation result is always built non-periodic in OCCT
//! (`GeomConvert_CompCurveToBSplineCurve.cxx` L250 uses the non-periodic
//! constructor); the periodic-only Segment branch carries the architecture
//! annotation.

use glam::DVec3;

use crate::base::convert::ConvertParameterisation;
use crate::base::extrema_ext_elc::epsilon_of;
use crate::geom::{BezierCurve3, BSplineCurve3, Curve3, CurveEval};
use crate::math::bspl_lib::{
    at, ati, build_cache_3d, first_uknot_index_mults, flat_bezier_knots, insert_knots,
    increase_degree, increase_degree_count_knots, last_uknot_index_mults,
    locate_parameter_knots_mults, pole_index, prepare_insert_knots, remove_knot,
};
use crate::math::plib::{coefficients_poles_3d, trimming_3d};

/// OCCT BSplCLib::MaxDegree() == 25.
const BSPLIB_MAX_DEGREE: usize = 25;

// ---------------------------------------------------------------------------
// Legacy-carrier helpers.
// ---------------------------------------------------------------------------

/// Split the flat knot vector into OCCT (knots, mults).
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

/// Build the flat knot vector from OCCT (knots, mults).
fn flat_knots_of(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut flat = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            flat.push(*k);
        }
    }
    flat
}

/// OCCT Geom_BSplineCurve::IsRational — the legacy carrier observes
/// rationality as "any weight != 1" (`geom::bspline_ops` precedent).
fn bspline3_is_rational(c: &BSplineCurve3) -> bool {
    c.weights.iter().any(|&w| w != 1.0)
}

/// OCCT PLib::SetPoles — flatten 3D poles for the BSplCLib kernels.
/// Non-rational: (x, y, z) per pole, dim 3.  Rational
/// (PLib::SetPoles(Poles, Weights, FPoles)): homogeneous (x*w, y*w, z*w, w)
/// per pole, dim 4.
fn poles_flat(c: &BSplineCurve3) -> (Vec<f64>, usize) {
    if bspline3_is_rational(c) {
        let mut flat = Vec::with_capacity(c.control_points.len() * 4);
        for (p, &w) in c.control_points.iter().zip(c.weights.iter()) {
            flat.push(p.x * w);
            flat.push(p.y * w);
            flat.push(p.z * w);
            flat.push(w);
        }
        (flat, 4)
    } else {
        let mut flat = Vec::with_capacity(c.control_points.len() * 3);
        for p in &c.control_points {
            flat.push(p.x);
            flat.push(p.y);
            flat.push(p.z);
        }
        (flat, 3)
    }
}

/// OCCT PLib::GetPoles — restore poles (and weights when homogeneous).
fn poles_unflatten(flat: &[f64], dim: usize, rational: bool) -> (Vec<DVec3>, Vec<f64>) {
    if dim == 4 {
        let n = flat.len() / 4;
        let mut poles = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);
        for chunk in flat.chunks_exact(4) {
            weights.push(chunk[3]);
            poles.push(DVec3::new(
                chunk[0] / chunk[3],
                chunk[1] / chunk[3],
                chunk[2] / chunk[3],
            ));
        }
        (poles, weights)
    } else {
        let poles: Vec<DVec3> = flat
            .chunks_exact(3)
            .map(|ch| DVec3::new(ch[0], ch[1], ch[2]))
            .collect();
        let weights = if rational {
            Vec::new()
        } else {
            vec![1.0; poles.len()]
        };
        (poles, weights)
    }
}

// ---------------------------------------------------------------------------
// Geom_BSplineCurve operations over the legacy carrier.
// ---------------------------------------------------------------------------

/// OCCT Geom_BSplineCurve::IncreaseDegree(Degree)
/// (Geom_BSplineCurve.cxx L243-287) — rational curves homogenized to dim 4.
fn bspline3_increase_degree(c: &mut BSplineCurve3, degree: usize) {
    if degree == c.degree {
        return;
    }
    if degree < c.degree || degree > BSPLIB_MAX_DEGREE {
        panic!("Standard_ConstructionError: BSpline curve: IncreaseDegree: bad degree value");
    }

    let (knots, mults) = knots_mults_of(&c.knots);
    let from_k1 = first_uknot_index_mults(c.degree, &mults);
    let to_k2 = last_uknot_index_mults(c.degree, &mults);
    let step = (degree - c.degree) as i32;

    // OCCT: npoles(1, myPoles.Length() + Step * (ToK2 - FromK1)).
    let npoles_len = c.control_points.len() + (step * (to_k2 - from_k1)) as usize;
    let nbknots = increase_degree_count_knots(c.degree, degree, false, &mults);

    let (poles_flat_v, dim) = poles_flat(c);
    let mut npoles_flat = vec![0.0f64; npoles_len * dim];
    let mut nknots = vec![0.0f64; nbknots];
    let mut nmults = vec![0i32; nbknots];

    increase_degree(
        c.degree,
        degree,
        false,
        dim,
        &poles_flat_v,
        &knots,
        &mults,
        &mut npoles_flat,
        &mut nknots,
        &mut nmults,
    );

    let rational = bspline3_is_rational(c);
    c.degree = degree;
    let (poles, weights) = poles_unflatten(&npoles_flat, dim, rational);
    c.control_points = poles;
    c.weights = weights;
    c.knots = flat_knots_of(&nknots, &nmults);
}

/// OCCT Geom_BSplineCurve::RemoveKnot(Index, M, Tolerance)
/// (Geom_BSplineCurve.cxx L420-492) — rational curves homogenized to dim 4.
fn bspline3_remove_knot(c: &mut BSplineCurve3, index: i32, m: i32, tolerance: f64) -> bool {
    if m < 0 {
        return true;
    }

    let (knots, mults) = knots_mults_of(&c.knots);
    let i1 = first_uknot_index_mults(c.degree, &mults);
    let i2 = last_uknot_index_mults(c.degree, &mults);
    if index < i1 || index > i2 {
        panic!("Standard_OutOfRange: BSpline curve: RemoveKnot: index out of range");
    }

    // OCCT L441: int step = myMults.Value(Index) - M;
    let step = ati(&mults, index) - m;
    if step <= 0 {
        return true;
    }

    let old_poles_len = c.control_points.len();
    let new_poles_len = old_poles_len - step as usize;
    let new_knots_len = knots.len() - if m == 0 { 1 } else { 0 };
    let rational = bspline3_is_rational(c);
    let (poles_flat_v, dim) = poles_flat(c);
    let mut npoles_flat = vec![0.0f64; new_poles_len * dim];
    let mut nknots = vec![0.0f64; new_knots_len];
    let mut nmults = vec![0i32; new_knots_len];

    if !remove_knot(
        index as usize,
        m,
        c.degree,
        false,
        dim,
        &poles_flat_v,
        &knots,
        &mults,
        &mut npoles_flat,
        &mut nknots,
        &mut nmults,
        tolerance,
    ) {
        return false;
    }

    let (poles, weights) = poles_unflatten(&npoles_flat, dim, rational);
    c.control_points = poles;
    c.weights = weights;
    c.knots = flat_knots_of(&nknots, &nmults);
    true
}

/// OCCT Geom_BSplineCurve::InsertKnots(Knots, Mults, Epsilon, Add)
/// (Geom_BSplineCurve.cxx L351-416) — the Segment support (Add = false).
fn bspline3_insert_knots(c: &mut BSplineCurve3, knots_add: &[f64], mults_add: &[i32], eps: f64) {
    let (knots, mults) = knots_mults_of(&c.knots);
    let mut nbpoles = 0i32;
    let mut nbknots = 0i32;
    if !prepare_insert_knots(
        c.degree,
        false,
        &knots,
        &mults,
        knots_add,
        Some(mults_add),
        &mut nbpoles,
        &mut nbknots,
        eps,
        false,
    ) {
        panic!("Standard_ConstructionError: Geom_BSplineCurve::InsertKnots");
    }
    if nbpoles as usize == c.control_points.len() {
        return;
    }

    let rational = bspline3_is_rational(c);
    let (poles_flat_v, dim) = poles_flat(c);
    let mut npoles_flat = vec![0.0f64; nbpoles as usize * dim];
    let mut nknots = vec![0.0f64; nbknots as usize];
    let mut nmults = vec![0i32; nbknots as usize];

    insert_knots(
        c.degree,
        false,
        dim,
        &poles_flat_v,
        &knots,
        &mults,
        knots_add,
        Some(mults_add),
        &mut npoles_flat,
        &mut nknots,
        &mut nmults,
        eps,
        false,
    );

    let (poles, weights) = poles_unflatten(&npoles_flat, dim, rational);
    c.control_points = poles;
    c.weights = weights;
    c.knots = flat_knots_of(&nknots, &nmults);
}

/// OCCT Geom_BSplineCurve::Segment(aU1, aU2, theTolerance)
/// (Geom_BSplineCurve.cxx L527-715) — the non-periodic mainline.  The
/// periodic-only branch (SetOrigin/SetNotPeriodic, L599-628) carries the
/// architecture annotation (see module docs).
fn bspline3_segment(c: &mut BSplineCurve3, au1: f64, au2: f64, the_tolerance: f64) {
    if au2 < au1 {
        panic!("Standard_DomainError: Geom_BSplineCurve::Segment");
    }

    let u1 = au1;
    let u2 = au2;

    let mut new_u1 = 0.0f64;
    let mut new_u2 = 0.0f64;
    {
        let (knots, mults) = knots_mults_of(&c.knots);
        let knots_upper = knots.len() as i32;
        let mut index = 0i32;
        locate_parameter_knots_mults(
            c.degree,
            &knots,
            &mults,
            u1,
            false,
            1,
            knots_upper,
            &mut index,
            &mut new_u1,
        );
        let mut index2 = 0i32;
        locate_parameter_knots_mults(
            c.degree,
            &knots,
            &mults,
            u2,
            false,
            1,
            knots_upper,
            &mut index2,
            &mut new_u2,
        );
    }

    let mut abs_umax = new_u1.abs().max(new_u2.abs());
    abs_umax = abs_umax
        .max(c.first_parameter().abs())
        .max(c.last_parameter().abs());
    let eps = epsilon_of(abs_umax).max(the_tolerance);

    let knots2 = [new_u1.min(new_u2), new_u1.max(new_u2)];
    let mults2 = [c.degree as i32, c.degree as i32];
    bspline3_insert_knots(c, &knots2, &mults2, eps);

    // compute index1 and index2 to set the new knots and mults (L621-635).
    let (knots, mults) = knots_mults_of(&c.knots);
    let from_u1 = 1i32;
    let to_u2 = knots.len() as i32;
    let mut index1 = 0i32;
    let mut index2 = 0i32;
    let mut u = 0.0f64;
    locate_parameter_knots_mults(
        c.degree,
        &knots,
        &mults,
        new_u1,
        false,
        from_u1,
        to_u2,
        &mut index1,
        &mut u,
    );
    if (at(&knots, index1 + 1) - u).abs() <= eps {
        index1 += 1;
    }
    locate_parameter_knots_mults(
        c.degree,
        &knots,
        &mults,
        new_u2,
        false,
        from_u1,
        to_u2,
        &mut index2,
        &mut u,
    );
    if (at(&knots, index2 + 1) - u).abs() <= eps || index2 == index1 {
        index2 += 1;
    }

    let nbknots = (index2 - index1 + 1) as usize;
    let mut nknots = vec![0.0f64; nbknots];
    let mut nmults = vec![0i32; nbknots];
    let mut k = 1i32;
    let mut i = index1;
    while i <= index2 {
        nknots[(k - 1) as usize] = at(&knots, i);
        nmults[(k - 1) as usize] = ati(&mults, i);
        k += 1;
        i += 1;
    }
    nmults[0] = c.degree as i32 + 1;
    nmults[nbknots - 1] = c.degree as i32 + 1;
    // compute index1 and index2 to set the new poles and weights (L658-665).
    let mut pindex1 = pole_index(c.degree, index1, false, &mults);
    let mut pindex2 = pole_index(c.degree, index2, false, &mults);
    pindex1 += 1;
    pindex2 = (pindex2 + 1).min(c.control_points.len() as i32);
    let nbpoles = (pindex2 - pindex1 + 1) as usize;

    let mut npoles = Vec::with_capacity(nbpoles);
    let mut nweights = Vec::with_capacity(nbpoles);
    let mut i = pindex1;
    while i <= pindex2 {
        npoles.push(c.control_points[(i - 1) as usize]);
        nweights.push(c.weights[(i - 1) as usize]);
        i += 1;
    }

    c.knots = flat_knots_of(&nknots, &nmults);
    c.control_points = npoles;
    c.weights = nweights;
}

// ---------------------------------------------------------------------------
// OCCT GeomConvert.cxx L163-380 — CurveToBSplineCurve (the branch the
// concatenator needs beyond the existing `geom_convert_curve_to_bspline_curve`).
// ---------------------------------------------------------------------------

/// OCCT GeomConvert::CurveToBSplineCurve(C, Parameterisation) — the
/// trimmed-BSpline branch (GeomConvert.cxx L336-340: Copy + Segment over the
/// clamped trim parameters L177-188); every other non-BSpline input goes
/// through the existing OCCT-anchored kernel conversion (trimmed line /
/// trimmed circle exact, further branches staged there).
pub fn curve_to_bspline_curve_3d(
    c: &Curve3,
    parameterisation: ConvertParameterisation,
) -> BSplineCurve3 {
    if let Curve3::Trimmed(ctrim) = c {
        let curv: &Curve3 = &ctrim.curve;
        if let Curve3::BSpline(bs) = curv {
            // OCCT L177-188: if the curve is not truly restricted, there is
            // no risk of a Raise in BS->Segment (clamp U1/U2).
            let mut u1 = ctrim.first;
            let mut u2 = ctrim.last;
            if !curv.is_periodic() {
                if u1 < bs.first_parameter() {
                    u1 = bs.first_parameter();
                }
                if u2 > bs.last_parameter() {
                    u2 = bs.last_parameter();
                }
            }
            // OCCT L336-340.
            let mut the_curve = bs.clone();
            bspline3_segment(&mut the_curve, u1, u2, 0.0);
            return the_curve;
        }
        if let Curve3::Bezier(cbez) = curv {
            // OCCT L177-188: the clamp (a Geom_BezierCurve is never
            // periodic; FirstParameter() == 0, LastParameter() == 1).
            let mut u1 = ctrim.first;
            let mut u2 = ctrim.last;
            if !curv.is_periodic() {
                if u1 < 0.0 {
                    u1 = 0.0;
                }
                if u2 > 1.0 {
                    u2 = 1.0;
                }
            }
            // OCCT L300-321: Copy + Geom_BezierCurve::Segment, then the
            // clamped-knot BSpline build (knots 0/1, mults Degree+1).
            let mut the_bez = cbez.clone();
            bezier3_segment(&mut the_bez, u1, u2);
            let degree = the_bez.control_points.len() - 1;
            let kts = [0.0f64, 1.0];
            let mults = [degree as i32 + 1, degree as i32 + 1];
            return BSplineCurve3 {
                degree,
                knots: flat_knots_of(&kts, &mults),
                control_points: the_bez.control_points,
                weights: the_bez.weights,
                is_periodic: false,
            };
        }
    }
    crate::base::convert::geom_convert_curve_to_bspline_curve(c, parameterisation)
}

/// OCCT Geom_BezierCurve::Segment(U1, U2) (Geom_BezierCurve.cxx L388-425):
/// reparameterizes the Bezier onto the sub-interval [U1, U2] of [0, 1] —
/// BSplCLib::BuildCache(0, 1, false, aDeg, KnotSequence(), ...) produces
/// the Taylor (power) coefficients at 0, PLib::Trimming performs the
/// substitution u = U1 + v*(U2-U1), PLib::CoefficientsPoles converts the
/// power coefficients back to Bernstein poles.
fn bezier3_segment(c: &mut BezierCurve3, u1: f64, u2: f64) {
    // OCCT L390: myClosed = (|Value(U1).Distance(Value(U2))| <=
    // Precision::Confusion()) — the legacy BezierCurve3 carrier carries no
    // closed flag (architecture note).
    // OCCT L392: const int aDeg = myPoles.Length() - 1.
    let degree = c.control_points.len() as i32 - 1;
    // KnotSequence() — BSplCLib::FlatBezierKnots(Degree) (BSplCLib.cxx
    // L4971-4977).
    let knot_sequence = flat_bezier_knots(degree);
    // OCCT L394: NCollection_Array1<gp_Pnt> coeffs(1, myPoles.Length()).
    let mut coeffs = c.control_points.clone();
    if c.weights.iter().any(|&w| w != 1.0) {
        // OCCT L395-409: the rational arm.
        let mut wcoeffs = c.weights.clone();
        build_cache_3d(
            0.0,
            1.0,
            false,
            degree,
            &knot_sequence,
            &c.control_points,
            Some(&c.weights),
            &mut coeffs,
            Some(&mut wcoeffs),
        );
        trimming_3d(u1, u2, &mut coeffs, Some(&mut wcoeffs));
        coefficients_poles_3d(
            &coeffs,
            Some(&wcoeffs),
            &mut c.control_points,
            Some(&mut c.weights),
        );
    } else {
        // OCCT L410-423: the non-rational arm — NoWeights throughout.
        build_cache_3d(
            0.0,
            1.0,
            false,
            degree,
            &knot_sequence,
            &c.control_points,
            None,
            &mut coeffs,
            None,
        );
        trimming_3d(u1, u2, &mut coeffs, None);
        coefficients_poles_3d(&coeffs, None, &mut c.control_points, None);
    }
    // OCCT L424: myMaxDerivInvOk = false — no carrier field
    // (architecture note).
}

// ---------------------------------------------------------------------------
// OCCT GeomConvert_CompCurveToBSplineCurve.cxx L32-274 — the concatenator.
// ---------------------------------------------------------------------------

/// OCCT GeomConvert_CompCurveToBSplineCurve — converts and concatenates
/// several curves into one BSplineCurve (CompCurveToBSplineCurve.hxx L30-78).
pub struct GeomConvertCompCurveToBSplineCurve {
    /// OCCT: occ::handle<Geom_BSplineCurve> myCurve.
    my_curve: Option<BSplineCurve3>,
    /// OCCT: double myTol.
    my_tol: f64,
    /// OCCT: Convert_ParameterisationType myType.
    my_type: ConvertParameterisation,
}

impl GeomConvertCompCurveToBSplineCurve {
    /// OCCT L32-37 — ctor(Parameterisation = Convert_TgtThetaOver2);
    /// myTol = Precision::Confusion().
    pub fn new(the_parameterisation: ConvertParameterisation) -> Self {
        GeomConvertCompCurveToBSplineCurve {
            my_curve: None,
            my_tol: crate::core::precision::CONFUSION,
            my_type: the_parameterisation,
        }
    }

    /// OCCT L41-56 — ctor(BasisCurve, Parameterisation): the BSpline basis
    /// is copied, every other type goes through CurveToBSplineCurve.
    pub fn with_basis_curve(basis_curve: &Curve3, parameterisation: ConvertParameterisation) -> Self {
        let mut out = Self::new(parameterisation);
        match basis_curve {
            // OCCT L47-51: down_cast<Geom_BSplineCurve> non-null -> Copy().
            Curve3::BSpline(bs) => {
                out.my_curve = Some(bs.clone());
            }
            // OCCT L52-55: GeomConvert::CurveToBSplineCurve(BasisCurve, myType).
            _ => {
                out.my_curve = Some(curve_to_bspline_curve_3d(basis_curve, out.my_type));
            }
        }
        out
    }

    /// OCCT L60-131 — Add(NewCurve, Tolerance, After = false, WithRatio =
    /// true, MinM = 0) with the OCCT defaults.
    pub fn add(&mut self, new_curve: &Curve3, tolerance: f64) -> bool {
        self.add_full(new_curve, tolerance, false, true, 0)
    }

    /// OCCT L60-131 — Add with the explicit After argument (the
    /// BRepOffset_Tool::Glue call form, BRepOffset_Tool.cxx L1217).
    pub fn add_after(&mut self, new_curve: &Curve3, tolerance: f64, after: bool) -> bool {
        self.add_full(new_curve, tolerance, after, true, 0)
    }

    /// OCCT L60-131 — Add(NewCurve, Tolerance, After, WithRatio, MinM).
    pub fn add_full(
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
            _ => curve_to_bspline_curve_3d(new_curve, self.my_type),
        };
        if self.my_curve.is_none() {
            // OCCT L76-80.
            self.my_curve = Some(bs);
            return true;
        }

        self.my_tol = tolerance;

        // Use actual curve endpoints instead of poles for proper G0
        // continuity check (OCCT L85-91).
        let my = self.my_curve.as_ref().expect("myCurve");
        let a_curve_start = my.point_at(my.first_parameter());
        let a_curve_end = my.point_at(my.last_parameter());
        let a_bs_start = bs.point_at(bs.first_parameter());
        let a_bs_end = bs.point_at(bs.last_parameter());

        // OCCT L93-94.
        let mut avant = a_curve_start.distance(a_bs_start) < self.my_tol
            || a_curve_start.distance(a_bs_end) < self.my_tol;
        let mut apres = a_curve_end.distance(a_bs_start) < self.my_tol
            || a_curve_end.distance(a_bs_end) < self.my_tol;

        // Will myCurve be (or become) closed? (OCCT L96-107).
        if avant && apres {
            if after {
                avant = false;
            } else {
                apres = false;
            }
        }

        // Append after? (OCCT L109-118).
        if apres {
            if a_curve_end.distance(a_bs_end) < self.my_tol {
                bs = bs.reversed();
            }
            let mut first_curve = self.my_curve.take().expect("myCurve");
            let merged =
                add_concat(&mut first_curve, &mut bs, true, with_ratio, min_m, self.my_tol);
            self.my_curve = Some(merged);
            return true;
        }
        // Prepend before? (OCCT L119-128).
        else if avant {
            if a_curve_start.distance(a_bs_start) < self.my_tol {
                bs = bs.reversed();
            }
            let mut first_curve = bs;
            let mut second_curve = self.my_curve.take().expect("myCurve");
            let merged = add_concat(
                &mut first_curve,
                &mut second_curve,
                false,
                with_ratio,
                min_m,
                self.my_tol,
            );
            self.my_curve = Some(merged);
            return true;
        }

        // OCCT L130.
        false
    }

    /// OCCT L264-267 — BSplineCurve().
    pub fn bspline_curve(&self) -> Option<BSplineCurve3> {
        self.my_curve.clone()
    }

    /// OCCT L271-274 — Clear().
    pub fn clear(&mut self) {
        self.my_curve = None;
    }
}

/// OCCT GeomConvert_CompCurveToBSplineCurve.cxx L135-260 — the private
/// Add(FirstCurve, SecondCurve, After, WithRatio, MinM): pole-by-pole
/// concatenation of two BSplines at a degree-multiplicity junction with the
/// knot epsilon clamp, followed by the RemoveKnot reduction down to MinM.
/// FirstCurve is raised in degree in place; the merged curve is returned
/// (the OCCT callee stores it into myCurve).
fn add_concat(
    first_curve: &mut BSplineCurve3,
    second_curve: &mut BSplineCurve3,
    after: bool,
    with_ratio: bool,
    min_m: i32,
    my_tol: f64,
) -> BSplineCurve3 {
    // Harmonize the degrees (OCCT L141-150).
    let deg = first_curve.degree.max(second_curve.degree);
    if first_curve.degree < deg {
        bspline3_increase_degree(first_curve, deg);
    }
    if second_curve.degree < deg {
        bspline3_increase_degree(second_curve, deg);
    }

    // Declarations (OCCT L152-161): 1-based OCCT arrays -> Vec + index - 1.
    let mut ratio = 1.0f64;
    let nb_p1 = first_curve.control_points.len() as i32;
    let nb_p2 = second_curve.control_points.len() as i32;
    let (knots1, mults1) = knots_mults_of(&first_curve.knots);
    let (knots2, mults2) = knots_mults_of(&second_curve.knots);
    let nb_k1 = knots1.len() as i32;
    let nb_k2 = knots2.len() as i32;
    let mut noeuds = vec![0.0f64; (nb_k1 + nb_k2 - 1) as usize];
    let mut poles = vec![DVec3::ZERO; (nb_p1 + nb_p2 - 1) as usize];
    let mut poids = vec![0.0f64; (nb_p1 + nb_p2 - 1) as usize];
    let mut mults = vec![0i32; (nb_k1 + nb_k2 - 1) as usize];

    // Reparameterization ratio (C1 if possible) (OCCT L163-177).
    if with_ratio {
        let l1 = first_curve.dn(first_curve.last_parameter(), 1).length();
        let l2 = second_curve.dn(second_curve.first_parameter(), 1).length();
        if l1 > crate::core::precision::CONFUSION && l2 > crate::core::precision::CONFUSION {
            ratio = l1 / l2;
        }
        if ratio < crate::core::precision::CONFUSION
            || ratio > 1.0 / crate::core::precision::CONFUSION
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
        delta2 = ratio2 * at(&knots2, 1) - at(&knots1, nb_k1);
    } else {
        // Do not move the second curve (OCCT L187-194).
        ratio1 = ratio;
        delta1 = ratio1 * at(&knots1, nb_k1) - at(&knots2, 1);
        ratio2 = 1.0;
        delta2 = 0.0;
    }

    // The Knots (OCCT L196-229).
    let mut ii = 1i32;
    while ii <= nb_k1 {
        noeuds[(ii - 1) as usize] = ratio1 * at(&knots1, ii) - delta1;
        if ii > 1 {
            let mut eps = epsilon_of(noeuds[(ii - 2) as usize].abs());
            if eps < 5.0e-10 {
                eps = 5.0e-10;
            }
            if noeuds[(ii - 1) as usize] - noeuds[(ii - 2) as usize] <= eps {
                noeuds[(ii - 1) as usize] += eps;
            }
        }
        mults[(ii - 1) as usize] = ati(&mults1, ii);
        ii += 1;
    }
    mults[(nb_k1 - 1) as usize] = first_curve.degree as i32;
    let mut jj = nb_k1 + 1;
    ii = 2;
    while ii <= nb_k2 {
        noeuds[(jj - 1) as usize] = ratio2 * at(&knots2, ii) - delta2;
        let mut eps = epsilon_of(noeuds[(jj - 2) as usize].abs());
        if eps < 5.0e-10 {
            eps = 5.0e-10;
        }
        if noeuds[(jj - 1) as usize] - noeuds[(jj - 2) as usize] <= eps {
            noeuds[(jj - 1) as usize] += eps;
        }
        mults[(jj - 1) as usize] = ati(&mults2, ii);
        ii += 1;
        jj += 1;
    }

    // OCCT L231-232: Ratio = Weight(NbP1) / Weight(1).
    ratio = first_curve.weights[(nb_p1 - 1) as usize] / second_curve.weights[0];
    // The Poles and Weights (OCCT L233-247).
    let mut ii = 1i32;
    while ii < nb_p1 {
        poles[(ii - 1) as usize] = first_curve.control_points[(ii - 1) as usize];
        poids[(ii - 1) as usize] = first_curve.weights[(ii - 1) as usize];
        ii += 1;
    }
    let mut jj = nb_p1;
    let mut ii = 1i32;
    while ii <= nb_p2 {
        poles[(jj - 1) as usize] = second_curve.control_points[(ii - 1) as usize];
        // Note: the weights may not necessarily connect with C0 continuity,
        // hence the use of Ratio (OCCT L242-246).
        poids[(jj - 1) as usize] = ratio * second_curve.weights[(ii - 1) as usize];
        ii += 1;
        jj += 1;
    }

    // Create the BSpline (OCCT L249-250).
    let mut my_curve = BSplineCurve3 {
        degree: deg,
        knots: flat_knots_of(&noeuds, &mults),
        control_points: poles,
        weights: poids,
        is_periodic: false,
    };

    // Optionally reduce multiplicity down to MinM (OCCT L252-259).
    let mut ok = true;
    let mut m = mults[(nb_k1 - 1) as usize];
    while m > min_m && ok {
        m -= 1;
        ok = bspline3_remove_knot(&mut my_curve, nb_k1, m, my_tol);
    }
    my_curve
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Two collinear trimmed lines concatenate into one degree-1 BSpline
    /// (the BRepFill_SectionPlacement call form: Add(TC, tol) with the OCCT
    /// defaults).
    #[test]
    fn concat_two_lines() {
        use crate::geom::{Line3, TrimmedCurve3};
        let l1 = Curve3::Trimmed(TrimmedCurve3::new(
            Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X)),
            0.0,
            1.0,
        ));
        let l2 = Curve3::Trimmed(TrimmedCurve3::new(
            Curve3::Line(Line3::new(DVec3::new(1.0, 0.0, 0.0), DVec3::Y)),
            0.0,
            1.0,
        ));
        let mut comp =
            GeomConvertCompCurveToBSplineCurve::new(ConvertParameterisation::TgtThetaOver2);
        assert!(comp.add(&l1, 1.0e-7));
        assert!(comp.add(&l2, 1.0e-7));
        let bs = comp.bspline_curve().expect("merged");
        assert_eq!(bs.control_points.len(), 3);
        assert!(
            (bs.control_points[2] - DVec3::new(1.0, 1.0, 0.0)).length() < 1.0e-9,
            "end pole at the second line end"
        );
    }

    /// The trimmed-BSpline branch segments the basis (OCCT L336-340).
    #[test]
    fn trimmed_bspline_branch_segments() {
        let bs = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        };
        let out = curve_to_bspline_curve_3d(
            &Curve3::Trimmed(crate::geom::TrimmedCurve3::new(
                Curve3::BSpline(bs),
                0.25,
                0.75,
            )),
            ConvertParameterisation::TgtThetaOver2,
        );
        assert!((out.control_points[0] - DVec3::new(0.25, 0.0, 0.0)).length() < 1.0e-9);
        assert!((out.control_points[1] - DVec3::new(0.75, 0.0, 0.0)).length() < 1.0e-9);
    }

    /// The trimmed-Bezier branch (GeomConvert.cxx L300-321): Copy +
    /// Geom_BezierCurve::Segment — the quadratic (0,0,0),(2,4,0),(4,0,0)
    /// (x(u) = 4u, y(u) = 8u(1-u)) trimmed to [0.25, 0.75]
    /// (u = (1 + 2v)/4):
    ///   x(v) = 1 + 2v                    -> power coeffs (1, 2, 0)
    ///   y(v) = 8(1+2v)(3-2v)/16 = 1.5 + 2v - 2v^2
    ///                                    -> power coeffs (1.5, 2, -2)
    /// and the PLib::CoefficientsPoles Pascal recombination over both passes
    /// gives the poles (1, 1.5, 0), (2, 2.5, 0), (3, 1.5, 0) over the
    /// clamped knots [0, 1], mults 3.
    #[test]
    fn trimmed_bezier_branch_segments() {
        let cbez = BezierCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(2.0, 4.0, 0.0),
                DVec3::new(4.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let out = curve_to_bspline_curve_3d(
            &Curve3::Trimmed(crate::geom::TrimmedCurve3::new(
                Curve3::Bezier(cbez),
                0.25,
                0.75,
            )),
            ConvertParameterisation::TgtThetaOver2,
        );
        assert_eq!(out.degree, 2);
        assert_eq!(out.knots, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let want_poles = [
            DVec3::new(1.0, 1.5, 0.0),
            DVec3::new(2.0, 2.5, 0.0),
            DVec3::new(3.0, 1.5, 0.0),
        ];
        for (p, w) in out.control_points.iter().zip(want_poles.iter()) {
            assert!((p - w).length() < 1.0e-12, "pole {p:?} vs {w:?}");
        }
        // On-curve check at v = 1/2: u = 1/2, x = 2, y = 8*(1/2)*(1/2) = 2.
        let mid = CurveEval::point_at(&out, 0.5);
        assert!((mid - DVec3::new(2.0, 2.0, 0.0)).length() < 1.0e-12, "mid {mid:?}");
    }
}
