//! OCCT Geom2dConvert_CompCurveToBSplineCurve (TKGeomBase/Geom2dConvert)
//! plus the `Geom2dConvert::CurveToBSplineCurve` / `SplitBSplineCurve`
//! static entries consumed by the fillet/offset pipelines.
//!
//! Sources (1:1 statement-mapped):
//! - `Geom2dConvert_CompCurveToBSplineCurve.cxx` L29-244 (+ .hxx)
//! - `Geom2dConvert.cxx` L70-98 (BSplineCurveBuilder), L102-177
//!   (SplitBSplineCurve overloads), L181-449 (CurveToBSplineCurve)
//! - `Geom2d_BSplineCurve.cxx` L236-293 (IncreaseDegree), L410-460
//!   (RemoveKnot), L677-694 (Reverse), L707-888 (Segment) — the curve-level
//!   operations the concatenator mutates.
//!
//! Architecture note (same adaptation as `geom::bspline_ops`): the rcad
//! legacy `BSplineCurve2` carries the flat (expanded) knot vector instead of
//! OCCT's (knots, mults) pairs; the two are bijective and the conversions
//! here preserve OCCT semantics exactly.  The legacy struct carries no
//! periodic flag — the concatenation result is always built non-periodic in
//! OCCT (`Geom2dConvert_CompCurveToBSplineCurve.cxx` L220 uses the
//! non-periodic constructor) and the pipeline feeds clamped/trimmed inputs,
//! so the periodic-only branches (SetOrigin/SetNotPeriodic inside Segment)
//! carry the architecture annotation where OCCT would take them.

use glam::{DVec2, DVec3};

use crate::base::convert::{
    convert_circle_arc_to_bspline, convert_ellipse_arc_to_bspline, convert_hyperbola_to_bspline,
    convert_parabola_to_bspline, ConvertConicToBspline, ConvertParameterisation,
};
use crate::base::extrema_ext_elc::epsilon_of;
use crate::geom::{
    turn_2d, BezierCurve2, BSplineCurve2, Circle2d, Curve2d, Curve2dEval, Ellipse2d,
};
use crate::math::bspl_lib::{
    at, ati, build_cache_2d, first_uknot_index_mults, flat_bezier_knots, insert_knots,
    increase_degree, increase_degree_count_knots, last_uknot_index_mults,
    locate_parameter_knots_mults, pole_index, prepare_insert_knots, remove_knot,
};
use crate::math::plib::{coefficients_poles_2d, trimming_2d};

/// OCCT BSplCLib::MaxDegree() == 25.
const BSPLIB_MAX_DEGREE: usize = 25;

// ---------------------------------------------------------------------------
// Legacy-carrier helpers (flat knots <-> (knots, mults), poles <-> flat
// homogeneous kernels).  Same bijection as `geom::bspline_ops`.
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

/// OCCT Geom2d_BSplineCurve::IsRational — the legacy carrier observes
/// rationality as "any weight != 1" (`geom::bspline_ops` precedent).
fn bspline2_is_rational(c: &BSplineCurve2) -> bool {
    c.weights.iter().any(|&w| w != 1.0)
}

/// OCCT Geom2d_BSplineCurve::FirstParameter (_1.cxx L349-352).
fn bspline2_first_parameter(c: &BSplineCurve2) -> f64 {
    c.knots[c.degree]
}

/// OCCT Geom2d_BSplineCurve::LastParameter (_1.cxx L419-422).
fn bspline2_last_parameter(c: &BSplineCurve2) -> f64 {
    c.knots[c.knots.len() - c.degree - 1]
}

/// OCCT PLib::SetPoles — flatten 2D poles for the BSplCLib kernels.
/// Non-rational: (x, y) per pole, dim 2.  Rational (PLib::SetPoles(Poles,
/// Weights, FPoles)): homogeneous (x*w, y*w, w) per pole, dim 3.
fn poles_flat(c: &BSplineCurve2) -> (Vec<f64>, usize) {
    if bspline2_is_rational(c) {
        let mut flat = Vec::with_capacity(c.control_points.len() * 3);
        for (p, &w) in c.control_points.iter().zip(c.weights.iter()) {
            flat.push(p.x * w);
            flat.push(p.y * w);
            flat.push(w);
        }
        (flat, 3)
    } else {
        let mut flat = Vec::with_capacity(c.control_points.len() * 2);
        for p in &c.control_points {
            flat.push(p.x);
            flat.push(p.y);
        }
        (flat, 2)
    }
}

/// OCCT PLib::GetPoles — restore poles (and weights when homogeneous).
fn poles_unflatten(flat: &[f64], dim: usize, rational: bool) -> (Vec<DVec2>, Vec<f64>) {
    if dim == 3 {
        let n = flat.len() / 3;
        let mut poles = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);
        for chunk in flat.chunks_exact(3) {
            weights.push(chunk[2]);
            poles.push(DVec2::new(chunk[0] / chunk[2], chunk[1] / chunk[2]));
        }
        (poles, weights)
    } else {
        let poles: Vec<DVec2> = flat
            .chunks_exact(2)
            .map(|ch| DVec2::new(ch[0], ch[1]))
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
// Geom2d_BSplineCurve operations over the legacy carrier.
// ---------------------------------------------------------------------------

/// OCCT Geom2d_BSplineCurve::IncreaseDegree(Degree)
/// (Geom2d_BSplineCurve.cxx L236-293).  The BSplCLib gp_Pnt2d overload
/// dispatches through `BSplCLib_IncreaseDegree`: rational curves are
/// homogenized to dim 3, the flat Dimension kernel runs, the result is
/// dehomogenized.
fn bspline2_increase_degree(c: &mut BSplineCurve2, degree: usize) {
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

    let rational = bspline2_is_rational(c);
    c.degree = degree;
    let (poles, weights) = poles_unflatten(&npoles_flat, dim, rational);
    c.control_points = poles;
    c.weights = weights;
    c.knots = flat_knots_of(&nknots, &nmults);
}

/// OCCT Geom2d_BSplineCurve::RemoveKnot(Index, M, Tolerance)
/// (Geom2d_BSplineCurve.cxx L410-460).
fn bspline2_remove_knot(c: &mut BSplineCurve2, index: i32, m: i32, tolerance: f64) -> bool {
    if m < 0 {
        return true;
    }

    let (knots, mults) = knots_mults_of(&c.knots);
    let i1 = first_uknot_index_mults(c.degree, &mults);
    let i2 = last_uknot_index_mults(c.degree, &mults);
    if index < i1 || index > i2 {
        panic!("Standard_OutOfRange: BSpline curve: RemoveKnot: index out of range");
    }

    // OCCT L434: int step = myMults.Value(Index) - M;
    let step = ati(&mults, index) - m;
    if step <= 0 {
        return true;
    }

    let old_poles_len = c.control_points.len();
    let new_poles_len = old_poles_len - step as usize;
    let new_knots_len = knots.len() - if m == 0 { 1 } else { 0 };
    let rational = bspline2_is_rational(c);
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

/// OCCT Geom2d_BSplineCurve::Reverse (Geom2d_BSplineCurve.cxx L677-694):
/// BSplCLib::Reverse(myKnots) + Reverse(myMults) + Reverse(myPoles, last) +
/// Reverse(myWeights, last).  On the flat knot vector, reversing the
/// (knots, mults) arrays and re-expanding equals order-reversing the flat
/// array while reflecting each value about (kfirst + klast).
fn bspline2_reverse(c: &mut BSplineCurve2) {
    let kfirst = c.knots[0];
    let klast = c.knots[c.knots.len() - 1];
    c.knots = c.knots.iter().rev().map(|k| kfirst + klast - k).collect();
    c.control_points.reverse();
    if bspline2_is_rational(c) {
        c.weights.reverse();
    }
}

/// OCCT Geom2d_BSplineCurve::InsertKnots(Knots, Mults, Epsilon, Add)
/// (Geom2d_BSplineCurve.cxx L343-406) — the Segment support (Add = false).
fn bspline2_insert_knots(c: &mut BSplineCurve2, knots_add: &[f64], mults_add: &[i32], eps: f64) {
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
        panic!("Standard_ConstructionError: Geom2d_BSplineCurve::InsertKnots");
    }
    if nbpoles as usize == c.control_points.len() {
        return;
    }

    let rational = bspline2_is_rational(c);
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

/// OCCT Geom2d_BSplineCurve::Segment(aU1, aU2, theTolerance)
/// (Geom2d_BSplineCurve.cxx L707-888) — the non-periodic mainline.  The
/// periodic-only branch (SetOrigin/SetNotPeriodic, L795-828) cannot be
/// reached on the legacy carrier (no periodic flag; architecture note in the
/// module docs).
fn bspline2_segment(c: &mut BSplineCurve2, au1: f64, au2: f64, tolerance: f64) {
    if au2 < au1 {
        panic!("Standard_DomainError: Geom2d_BSplineCurve::Segment");
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
        .max(bspline2_first_parameter(c).abs())
        .max(bspline2_last_parameter(c).abs());
    let eps = epsilon_of(abs_umax).max(tolerance);

    let knots2 = [new_u1.min(new_u2), new_u1.max(new_u2)];
    let mults2 = [c.degree as i32, c.degree as i32];
    bspline2_insert_knots(c, &knots2, &mults2, eps);

    // compute index1 and index2 to set the new knots and mults (L830-856).
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

    // compute index1 and index2 to set the new poles and weights (L861-869).
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
// OCCT Geom2dConvert.cxx L70-98 — BSplineCurveBuilder (file static): builds
// the BSpline from the canonical conic conversion and places it on the
// conic position (handedness mirror + frame transformation).
// ---------------------------------------------------------------------------

/// OCCT BSplineCurveBuilder(TheConic, Convert) reduced to its 2D data flow
/// (Geom2dConvert.cxx L70-98): the canonical poles are mirrored about OX2d
/// when the conic frame is left-handed (L84-91), then mapped by
/// `T.SetTransformation(TheConic->XAxis(), gp::OX2d())` (L92-96).
///
/// The two gp_Trsf2d statements collapse (gp_Trsf2d.cxx L48-84): the OX
/// mirror is (x, y) -> (x, -y) and the XAxis transform with ToA2 = OX2d is
/// P = Loc + x*XDir + y*(-XDir.y, XDir.x) (matrix columns (XDir, perp(XDir)),
/// translation Loc; gp_XY::Multiply is the column-vector convention).  Over
/// the conic's own frame the net effect is P = Loc + x*XDir + y*YDir — the
/// left-handed YDir already carries the mirror, so the y-offset is NOT
/// flipped a second time.
fn bspline_curve_builder_2d(
    loc: DVec2,
    x_dir: DVec2,
    y_dir: DVec2,
    conv: &ConvertConicToBspline,
) -> BSplineCurve2 {
    // OCCT L84-91: (Axis.XDirection() ^ Axis.YDirection()) < 0 -> mirror.
    let cross = x_dir.x * y_dir.y - x_dir.y * y_dir.x;
    let mirror = cross < 0.0;
    // OCCT L93-96: the SetTransformation(XAxis, OX2d) collapse — the
    // transform maps through perp(XDir), not YDir; YDir enters only the
    // handedness test above.
    let perp_x = DVec2::new(-x_dir.y, x_dir.x);
    let control_points = conv
        .poles_2d
        .iter()
        .map(|p| {
            // OCCT L87-89: Sym.SetMirror(gp::OX2d()) + Transform —
            // gp_Trsf2d.cxx L48-62: (x, y) -> (x, -y).
            let y = if mirror { -p.y } else { p.y };
            // gp_Trsf2d.cxx L64-84: P = Loc + x*XDir + y*perp(XDir).
            loc + p.x * x_dir + y * perp_x
        })
        .collect();
    BSplineCurve2 {
        degree: conv.degree as usize,
        knots: flat_knots_of(&conv.knots, &conv.mults),
        control_points,
        weights: conv.weights.clone(),
    }
}

/// The canonical gp::OX2d() circle of radius R used by
/// `Geom2dConvert::CurveToBSplineCurve` for every conic (Geom2dConvert.cxx
/// L231/L270/L309/L318/L383/L393): `Circ2d C2d(gp::OX2d(), R)`.
fn canonical_circle(radius: f64) -> crate::geom::Circle3 {
    crate::geom::Circle3 {
        center: DVec3::ZERO,
        normal: DVec3::Z,
        x_dir: DVec3::X,
        y_dir: DVec3::Y,
        radius,
    }
}

// ---------------------------------------------------------------------------
// OCCT Geom2dConvert.cxx L181-449 — CurveToBSplineCurve(C, Parameterisation).
// ---------------------------------------------------------------------------

/// OCCT Geom2dConvert::CurveToBSplineCurve (Geom2dConvert.cxx L181-449) —
/// per-type exact conversion.  Branches whose underlying `Convert_*` /
/// `ApproxCurve` engines are not yet translated keep the OCCT failure path
/// (staged panics with the source anchor); no sampling approximation is
/// substituted where OCCT is exact.
pub fn curve_to_bspline_curve_2d(
    c: &Curve2d,
    parameterisation: ConvertParameterisation,
) -> BSplineCurve2 {
    if let Curve2d::Trimmed(ctrim) = c {
        // OCCT L189-209: unwrap the trim, clamp U1/U2 to the (non-periodic)
        // basis parameters.
        let curv: &Curve2d = &ctrim.curve;
        let mut u1 = ctrim.t_min;
        let mut u2 = ctrim.t_max;
        if !curv.is_periodic() {
            let dom = curv.default_domain();
            if u1 < dom[0] {
                u1 = dom[0];
            }
            if u2 > dom[1] {
                u2 = dom[1];
            }
        }

        match curv {
            // OCCT L211-226: the trimmed line — two poles at the trim ends
            // (Ctrim->StartPoint()/EndPoint(), the unclamped trim bounds)
            // over the trim-parameter knots.
            Curve2d::Line(_) => {
                let pdeb = curv.point_at(ctrim.t_min);
                let pfin = curv.point_at(ctrim.t_max);
                BSplineCurve2 {
                    degree: 1,
                    knots: vec![ctrim.t_min, ctrim.t_min, ctrim.t_max, ctrim.t_max],
                    control_points: vec![pdeb, pfin],
                    weights: vec![1.0, 1.0],
                }
            }
            // OCCT L228-264: the trimmed circle.
            Curve2d::Circle(the_conic) => {
                circle_arc_to_bspline_2d(the_conic, u1, u2, parameterisation)
            }
            // OCCT L266-303: the trimmed ellipse — Convert_EllipseToBSplineCurve
            // (kernel `convert_ellipse_arc_to_bspline`) with the RationalC1
            // over-length split through the concatenator.
            Curve2d::Ellipse(the_conic) => {
                ellipse_arc_to_bspline_2d(the_conic, u1, u2, parameterisation)
            }
            // OCCT L305-312: Convert_HyperbolaToBSplineCurve
            // (kernel `convert_hyperbola_to_bspline`), placed by the
            // BSplineCurveBuilder on the conic frame (y = turn_2d(major_dir),
            // the Conic2dEval YDirection convention).
            Curve2d::Hyperbola(the_conic) => {
                let conv =
                    convert_hyperbola_to_bspline(the_conic.semi_major, the_conic.semi_minor, u1, u2);
                bspline_curve_builder_2d(
                    the_conic.center,
                    the_conic.major_dir,
                    turn_2d(the_conic.major_dir),
                    &conv,
                )
            }
            // OCCT L314-321: Convert_ParabolaToBSplineCurve
            // (kernel `convert_parabola_to_bspline`); the rcad focal_param
            // is Prb.Parameter().
            Curve2d::Parabola(the_conic) => {
                let conv = convert_parabola_to_bspline(the_conic.focal_param, u1, u2);
                bspline_curve_builder_2d(
                    the_conic.origin,
                    the_conic.axis_dir,
                    turn_2d(the_conic.axis_dir),
                    &conv,
                )
            }
            // OCCT L323-345: the trimmed Bezier — Copy +
            // Geom2d_BezierCurve::Segment (BuildCache + PLib::Trimming +
            // PLib::CoefficientsPoles), then the clamped-knot BSpline build
            // (knots 0/1, mults Degree+1, over the segmented poles).
            Curve2d::Bezier(cbez) => {
                let mut the_bez = cbez.clone(); // OCCT L325: Curv->Copy().
                bezier2_segment(&mut the_bez, u1, u2); // OCCT L326.
                bezier_to_bspline_2d(&the_bez) // OCCT L327-345.
            }
            // OCCT L347-351: the trimmed BSpline — Copy + Segment.
            Curve2d::BSpline(bs) => {
                let mut the_curve = bs.clone();
                bspline2_segment(&mut the_curve, u1, u2, 0.0);
                the_curve
            }
            // OCCT L353-368: the offset branch — Geom2dConvert_ApproxCurve is
            // pending (staged).
            Curve2d::Offset(_) => {
                panic!(
                    "Staged: Geom2dConvert::CurveToBSplineCurve trimmed offset branch \
                     (Geom2dConvert.cxx L353-368, Geom2dConvert_ApproxCurve pending)"
                );
            }
            // OCCT L370-373: throw Standard_DomainError("No such curve").
            _ => panic!("Standard_DomainError: Geom2dConvert::CurveToBSplineCurve: No such curve"),
        }
    } else {
        match c {
            // OCCT L379-387: the full ellipse — the Convert engine
            // (`convert_ellipse_to_bspline_periodic`) is translated, but
            // the caller statement `TheCurve->SetPeriodic()` (L386) has no
            // equivalent on the legacy BSplineCurve2 carrier: the struct
            // carries no periodic flag, so the periodic knot layout the
            // engine produces cannot be represented (architecture
            // difference; BSplineCurve3 carries is_periodic and the 3D
            // twin branch is wired).
            Curve2d::Ellipse(_) => {
                panic!(
                    "Staged: Geom2dConvert::CurveToBSplineCurve periodic ellipse branch \
                     (Geom2dConvert.cxx L379-387, SetPeriodic unrepresentable on \
                     BSplineCurve2 — engine available in convert_conic_to_bspline)"
                );
            }
            // OCCT L389-397: the full circle -> SetPeriodic — pending (staged).
            Curve2d::Circle(_) => {
                panic!(
                    "Staged: Geom2dConvert::CurveToBSplineCurve periodic circle branch \
                     (Geom2dConvert.cxx L389-397, full-circle Convert pending)"
                );
            }
            // OCCT L399-419: the full Bezier — clamped knots of the poles.
            Curve2d::Bezier(cbez) => bezier_to_bspline_2d(cbez),
            // OCCT L420-423: the full BSpline — Copy.
            Curve2d::BSpline(bs) => bs.clone(),
            // OCCT L425-440: the offset branch — ApproxCurve pending (staged).
            Curve2d::Offset(_) => {
                panic!(
                    "Staged: Geom2dConvert::CurveToBSplineCurve offset branch \
                     (Geom2dConvert.cxx L425-440, Geom2dConvert_ApproxCurve pending)"
                );
            }
            // OCCT L442-445: throw Standard_DomainError.
            _ => panic!("Standard_DomainError: Geom2dConvert::CurveToBSplineCurve"),
        }
    }
}

/// OCCT L228-264 — the trimmed circle branch.  The canonical circle arc is
/// converted by `Convert_CircleToBSplineCurve` (kernel
/// `convert_circle_arc_to_bspline`) and placed by BSplineCurveBuilder; the
/// RationalC1 over-length split (L239-263) joins the two halves through the
/// concatenator itself.
fn circle_arc_to_bspline_2d(
    the_conic: &Circle2d,
    u1: f64,
    u2: f64,
    parameterisation: ConvertParameterisation,
) -> BSplineCurve2 {
    let canon = canonical_circle(the_conic.radius);
    if parameterisation != ConvertParameterisation::RationalC1 {
        // OCCT L232-236.
        let conv = convert_circle_arc_to_bspline(&canon, u1, u2, parameterisation);
        bspline_curve_builder_2d(the_conic.center, the_conic.x_dir, the_conic.y_dir, &conv)
    } else if u2 - u1 < 6.0 {
        // OCCT L239-243.
        let conv = convert_circle_arc_to_bspline(&canon, u1, u2, parameterisation);
        bspline_curve_builder_2d(the_conic.center, the_conic.x_dir, the_conic.y_dir, &conv)
    } else {
        // OCCT L244-263: split the circle to avoid numerical overflow when
        // U2 - U1 =~ 2*PI.
        let umed = (u1 + u2) * 0.5;
        let conv1 = convert_circle_arc_to_bspline(&canon, u1, umed, parameterisation);
        let the_curve1 =
            bspline_curve_builder_2d(the_conic.center, the_conic.x_dir, the_conic.y_dir, &conv1);
        let conv2 = convert_circle_arc_to_bspline(&canon, umed, u2, parameterisation);
        let the_curve2 =
            bspline_curve_builder_2d(the_conic.center, the_conic.x_dir, the_conic.y_dir, &conv2);

        let mut cctbspl = Geom2dConvertCompCurveToBSplineCurve::with_basis_curve(
            &Curve2d::BSpline(the_curve1),
            parameterisation,
        );
        cctbspl.add_after(&Curve2d::BSpline(the_curve2), crate::core::precision::PCONFUSION, true);
        cctbspl
            .bspline_curve()
            .expect("Geom2dConvert_CompCurveToBSplineCurve result after Add")
    }
}

/// OCCT L266-303 — the trimmed ellipse branch.  The canonical ellipse arc
/// is converted by `Convert_EllipseToBSplineCurve` (kernel
/// `convert_ellipse_arc_to_bspline`) and placed by BSplineCurveBuilder; the
/// RationalC1 over-length split (L283-302) joins the two halves through the
/// concatenator itself.
fn ellipse_arc_to_bspline_2d(
    the_conic: &Ellipse2d,
    u1: f64,
    u2: f64,
    parameterisation: ConvertParameterisation,
) -> BSplineCurve2 {
    if parameterisation != ConvertParameterisation::RationalC1 {
        // OCCT L271-275.
        let conv = convert_ellipse_arc_to_bspline(
            the_conic.major_radius,
            the_conic.minor_radius,
            u1,
            u2,
            parameterisation,
        );
        bspline_curve_builder_2d(
            the_conic.center,
            the_conic.major_dir,
            the_conic.minor_dir,
            &conv,
        )
    } else if u2 - u1 < 6.0 {
        // OCCT L278-282.
        let conv = convert_ellipse_arc_to_bspline(
            the_conic.major_radius,
            the_conic.minor_radius,
            u1,
            u2,
            parameterisation,
        );
        bspline_curve_builder_2d(
            the_conic.center,
            the_conic.major_dir,
            the_conic.minor_dir,
            &conv,
        )
    } else {
        // OCCT L283-302: split the ellipse to avoid numerical overflow when
        // U2 - U1 =~ 2*PI.
        let umed = (u1 + u2) * 0.5;
        let conv1 = convert_ellipse_arc_to_bspline(
            the_conic.major_radius,
            the_conic.minor_radius,
            u1,
            umed,
            parameterisation,
        );
        let the_curve1 = bspline_curve_builder_2d(
            the_conic.center,
            the_conic.major_dir,
            the_conic.minor_dir,
            &conv1,
        );
        let conv2 = convert_ellipse_arc_to_bspline(
            the_conic.major_radius,
            the_conic.minor_radius,
            umed,
            u2,
            parameterisation,
        );
        let the_curve2 = bspline_curve_builder_2d(
            the_conic.center,
            the_conic.major_dir,
            the_conic.minor_dir,
            &conv2,
        );

        let mut cctbspl = Geom2dConvertCompCurveToBSplineCurve::with_basis_curve(
            &Curve2d::BSpline(the_curve1),
            parameterisation,
        );
        cctbspl.add_after(&Curve2d::BSpline(the_curve2), crate::core::precision::PCONFUSION, true);
        cctbspl
            .bspline_curve()
            .expect("Geom2dConvert_CompCurveToBSplineCurve result after Add")
    }
}

/// OCCT L399-419 (the full Bezier branch): clamped knots 0/1 with
/// multiplicity Degree+1 over the Bezier poles.
fn bezier_to_bspline_2d(cbez: &BezierCurve2) -> BSplineCurve2 {
    let degree = cbez.control_points.len() - 1;
    let kts = [0.0f64, 1.0];
    let mults = [degree as i32 + 1, degree as i32 + 1];
    BSplineCurve2 {
        degree,
        knots: flat_knots_of(&kts, &mults),
        control_points: cbez.control_points.clone(),
        weights: cbez.weights.clone(),
    }
}

/// OCCT Geom2d_BezierCurve::Segment(U1, U2) (Geom2d_BezierCurve.cxx
/// L356-387): reparameterizes the Bezier onto the sub-interval [U1, U2] of
/// [0, 1] — BSplCLib::BuildCache(0, 1, false, Degree, KnotSequence(), ...)
/// produces the Taylor (power) coefficients at 0, PLib::Trimming performs
/// the substitution u = U1 + v*(U2-U1), PLib::CoefficientsPoles converts
/// the power coefficients back to Bernstein poles.
fn bezier2_segment(c: &mut BezierCurve2, u1: f64, u2: f64) {
    // OCCT L357: myClosed = (|Value(U1).Distance(Value(U2))| <=
    // gp::Resolution()) — the legacy BezierCurve2 carrier carries no closed
    // flag (architecture note); the statement has no other observable
    // effect here.
    // OCCT L359: NCollection_Array1<gp_Pnt2d> coeffs(1, myPoles.Length()).
    let degree = c.control_points.len() as i32 - 1;
    // KnotSequence() — BSplCLib::FlatBezierKnots(Degree) (BSplCLib.cxx
    // L4971-4977): Degree+1 zeros followed by Degree+1 ones.
    let knot_sequence = flat_bezier_knots(degree);
    let mut coeffs = c.control_points.clone();
    if c.weights.iter().any(|&w| w != 1.0) {
        // OCCT L360-375: the rational arm — BuildCache with &myWeights,
        // PLib::Trimming with &wcoeffs, PLib::CoefficientsPoles with
        // &myWeights.
        let mut wcoeffs = c.weights.clone();
        build_cache_2d(
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
        trimming_2d(u1, u2, &mut coeffs, Some(&mut wcoeffs));
        coefficients_poles_2d(
            &coeffs,
            Some(&wcoeffs),
            &mut c.control_points,
            Some(&mut c.weights),
        );
    } else {
        // OCCT L376-386: the non-rational arm — NoWeights throughout.
        build_cache_2d(
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
        trimming_2d(u1, u2, &mut coeffs, None);
        coefficients_poles_2d(&coeffs, None, &mut c.control_points, None);
    }
    // OCCT L387: myMaxDerivInvOk = false — no carrier field
    // (architecture note).
}

// ---------------------------------------------------------------------------
// OCCT Geom2dConvert.cxx L102-177 — SplitBSplineCurve overloads.
// ---------------------------------------------------------------------------

/// OCCT Geom2dConvert::SplitBSplineCurve(C, FromK1, ToK2, SameOrientation)
/// (Geom2dConvert.cxx L102-142) — the knot-index form.
pub fn split_bspline_curve_2d_knots(
    c: &BSplineCurve2,
    from_k1: i32,
    to_k2: i32,
    same_orientation: bool,
) -> BSplineCurve2 {
    let (knots, mults) = knots_mults_of(&c.knots);
    let the_first = first_uknot_index_mults(c.degree, &mults);
    let the_last = last_uknot_index_mults(c.degree, &mults);
    if from_k1 == to_k2 {
        panic!("Standard_DomainError: Geom2dConvert::SplitBSplineCurve");
    }
    let first_k = from_k1.min(to_k2);
    let last_k = from_k1.max(to_k2);
    if first_k < the_first || last_k > the_last {
        panic!("Standard_OutOfRange: Geom2dConvert::SplitBSplineCurve");
    }

    let mut new_curve = c.clone();
    bspline2_segment(&mut new_curve, at(&knots, first_k), at(&knots, last_k), 0.0);

    // The legacy carrier is never periodic (architecture note); the
    // non-periodic branch (L135-140) applies.
    if from_k1 > to_k2 {
        bspline2_reverse(&mut new_curve);
    }
    let _ = same_orientation;
    new_curve
}

/// OCCT Geom2dConvert::SplitBSplineCurve(C, FromU1, ToU2,
/// ParametricTolerance, SameOrientation) (Geom2dConvert.cxx L146-177) — the
/// parameter form (the ParametricTolerance argument is unused by the OCCT
/// body, L151).
pub fn split_bspline_curve_2d(
    c: &BSplineCurve2,
    from_u1: f64,
    to_u2: f64,
    _parametric_tolerance: f64,
    same_orientation: bool,
) -> BSplineCurve2 {
    let first_u = from_u1.min(to_u2);
    let last_u = from_u1.max(to_u2);

    let mut c1 = c.clone();
    bspline2_segment(&mut c1, first_u, last_u, 0.0);

    // The legacy carrier is never periodic (architecture note); the
    // non-periodic branch (L169-174) applies.
    if from_u1 > to_u2 {
        bspline2_reverse(&mut c1);
    }
    let _ = same_orientation;
    c1
}

// ---------------------------------------------------------------------------
// OCCT Geom2dConvert_CompCurveToBSplineCurve.cxx L29-244 — the incremental
// concatenator.
// ---------------------------------------------------------------------------

/// OCCT Geom2dConvert_CompCurveToBSplineCurve — converts and concatenates
/// several curves into one BSplineCurve (CompCurveToBSplineCurve.hxx L29-68).
pub struct Geom2dConvertCompCurveToBSplineCurve {
    /// OCCT: occ::handle<Geom2d_BSplineCurve> myCurve.
    my_curve: Option<BSplineCurve2>,
    /// OCCT: double myTol.
    my_tol: f64,
    /// OCCT: Convert_ParameterisationType myType.
    my_type: ConvertParameterisation,
}

impl Geom2dConvertCompCurveToBSplineCurve {
    /// OCCT L29-34 — ctor(Parameterisation = Convert_TgtThetaOver2);
    /// myTol = Precision::Confusion().
    pub fn new(the_parameterisation: ConvertParameterisation) -> Self {
        Geom2dConvertCompCurveToBSplineCurve {
            my_curve: None,
            my_tol: crate::core::precision::CONFUSION,
            my_type: the_parameterisation,
        }
    }

    /// OCCT L38-53 — ctor(BasisCurve, Parameterisation): the BSpline basis
    /// is copied, every other type goes through CurveToBSplineCurve.
    pub fn with_basis_curve(basis_curve: &Curve2d, parameterisation: ConvertParameterisation) -> Self {
        let mut out = Self::new(parameterisation);
        match basis_curve {
            // OCCT L44-48: down_cast<Geom2d_BSplineCurve> non-null -> Copy().
            Curve2d::BSpline(bs) => {
                out.my_curve = Some(bs.clone());
            }
            // OCCT L49-52: Geom2dConvert::CurveToBSplineCurve(BasisCurve, myType).
            _ => {
                out.my_curve = Some(curve_to_bspline_curve_2d(basis_curve, out.my_type));
            }
        }
        out
    }

    /// OCCT L57-125 — Add(NewCurve, Tolerance, After = false).
    pub fn add(&mut self, new_curve: &Curve2d, tolerance: f64) -> bool {
        self.add_after(new_curve, tolerance, false)
    }

    /// OCCT L57-125 — Add(NewCurve, Tolerance, After) with the explicit
    /// After argument.
    pub fn add_after(&mut self, new_curve: &Curve2d, tolerance: f64, after: bool) -> bool {
        // conversion (OCCT L61-70).
        let mut bs: BSplineCurve2 = match new_curve {
            Curve2d::BSpline(b) => b.clone(),
            _ => curve_to_bspline_curve_2d(new_curve, self.my_type),
        };
        if self.my_curve.is_none() {
            // OCCT L71-75.
            self.my_curve = Some(bs);
            return true;
        }

        self.my_tol = tolerance;
        let a_sq_tol = tolerance * tolerance;

        let my = self.my_curve.as_ref().expect("myCurve");
        let l_bs = bs.control_points.len() as i32;
        let l_cb = my.control_points.len() as i32;

        // OCCT L80-85.
        let d1 = my.control_points[0].distance_squared(bs.control_points[0]);
        let d2 = my.control_points[0].distance_squared(bs.control_points[(l_bs - 1) as usize]);
        let is_before_reversed =
            (my.control_points[0].distance_squared(bs.control_points[0]) < a_sq_tol) && (d1 < d2);
        let mut is_before = (my.control_points[0]
            .distance_squared(bs.control_points[(l_bs - 1) as usize])
            < a_sq_tol)
            || is_before_reversed;

        // OCCT L87-91.
        let d1 = my.control_points[(l_cb - 1) as usize].distance_squared(bs.control_points[0]);
        let d2 = my.control_points[(l_cb - 1) as usize]
            .distance_squared(bs.control_points[(l_bs - 1) as usize]);
        let is_after_reversed = (my.control_points[(l_cb - 1) as usize]
            .distance_squared(bs.control_points[(l_bs - 1) as usize])
            < a_sq_tol)
            && (d2 < d1);
        let mut is_after = (my.control_points[(l_cb - 1) as usize]
            .distance_squared(bs.control_points[0])
            < a_sq_tol)
            || is_after_reversed;

        // myCurve and NewCurve together form a closed curve (OCCT L93-104).
        if is_before && is_after {
            if after {
                is_before = false;
            } else {
                is_after = false;
            }
        }
        if is_after {
            // OCCT L105-113.
            if is_after_reversed {
                bspline2_reverse(&mut bs);
            }
            let mut first_curve = self.my_curve.take().expect("myCurve");
            let merged = add_concat(&mut first_curve, &mut bs, true, self.my_tol);
            self.my_curve = Some(merged);
            return true;
        } else if is_before {
            // OCCT L114-122.
            if is_before_reversed {
                bspline2_reverse(&mut bs);
            }
            let mut first_curve = bs;
            let mut second_curve = self.my_curve.take().expect("myCurve");
            let merged = add_concat(&mut first_curve, &mut second_curve, false, self.my_tol);
            self.my_curve = Some(merged);
            return true;
        }

        // OCCT L124.
        false
    }

    /// OCCT L234-237 — BSplineCurve().
    pub fn bspline_curve(&self) -> Option<BSplineCurve2> {
        self.my_curve.clone()
    }

    /// OCCT L241-244 — Clear().
    pub fn clear(&mut self) {
        self.my_curve = None;
    }
}

/// OCCT Geom2dConvert_CompCurveToBSplineCurve.cxx L129-230 — the private
/// Add(FirstCurve, SecondCurve, After): pole-by-pole concatenation of two
/// BSplines at a degree-multiplicity junction, followed by the
/// RemoveKnot reduction loop.  FirstCurve is raised in degree in place; the
/// merged curve is returned (the OCCT callee stores it into myCurve).
fn add_concat(
    first_curve: &mut BSplineCurve2,
    second_curve: &mut BSplineCurve2,
    after: bool,
    my_tol: f64,
) -> BSplineCurve2 {
    // Harmonize the degrees (OCCT L133-142).
    let deg = first_curve.degree.max(second_curve.degree);
    if first_curve.degree < deg {
        bspline2_increase_degree(first_curve, deg);
    }
    if second_curve.degree < deg {
        bspline2_increase_degree(second_curve, deg);
    }

    // Declarations (OCCT L144-153): 1-based OCCT arrays -> Vec + index - 1.
    let mut ratio = 1.0f64;
    let nb_p1 = first_curve.control_points.len() as i32;
    let nb_p2 = second_curve.control_points.len() as i32;
    let (knots1, mults1) = knots_mults_of(&first_curve.knots);
    let (knots2, mults2) = knots_mults_of(&second_curve.knots);
    let nb_k1 = knots1.len() as i32;
    let nb_k2 = knots2.len() as i32;
    let mut noeuds = vec![0.0f64; (nb_k1 + nb_k2 - 1) as usize];
    let mut poles = vec![DVec2::ZERO; (nb_p1 + nb_p2 - 1) as usize];
    let mut poids = vec![0.0f64; (nb_p1 + nb_p2 - 1) as usize];
    let mut mults = vec![0i32; (nb_k1 + nb_k2 - 1) as usize];

    // Reparameterization ratio (C1 if possible) (OCCT L155-166).
    let l1 = first_curve
        .derivative_at(bspline2_last_parameter(first_curve))
        .length();
    let l2 = second_curve
        .derivative_at(bspline2_first_parameter(second_curve))
        .length();
    if l1 > crate::core::precision::CONFUSION && l2 > crate::core::precision::CONFUSION {
        ratio = l1 / l2;
    }
    if ratio < crate::core::precision::CONFUSION || ratio > 1.0 / crate::core::precision::CONFUSION
    {
        ratio = 1.0;
    }

    let (ratio1, delta1, ratio2, delta2, u_de_raccord);
    if after {
        // Do not move the first curve (OCCT L168-176).
        ratio1 = 1.0;
        delta1 = 0.0;
        ratio2 = 1.0 / ratio;
        delta2 = ratio2 * at(&knots2, 1) - at(&knots1, nb_k1);
        u_de_raccord = bspline2_last_parameter(first_curve);
    } else {
        // Do not move the second curve (OCCT L177-185).
        ratio1 = ratio;
        delta1 = ratio1 * at(&knots1, nb_k1) - at(&knots2, 1);
        ratio2 = 1.0;
        delta2 = 0.0;
        u_de_raccord = bspline2_first_parameter(second_curve);
    }

    // The Knots (OCCT L187-200).
    let mut ii = 1i32;
    while ii < nb_k1 {
        noeuds[(ii - 1) as usize] = ratio1 * at(&knots1, ii) - delta1;
        mults[(ii - 1) as usize] = ati(&mults1, ii);
        ii += 1;
    }
    noeuds[(nb_k1 - 1) as usize] = u_de_raccord;
    mults[(nb_k1 - 1) as usize] = first_curve.degree as i32;
    let mut jj = nb_k1 + 1;
    ii = 2;
    while ii <= nb_k2 {
        noeuds[(jj - 1) as usize] = ratio2 * at(&knots2, ii) - delta2;
        mults[(jj - 1) as usize] = ati(&mults2, ii);
        ii += 1;
        jj += 1;
    }

    // OCCT L201-202: Ratio = Weight(NbP1) / Weight(1).
    ratio = first_curve.weights[(nb_p1 - 1) as usize] / second_curve.weights[0];
    // The Poles and Weights (OCCT L203-217).
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
        // hence the use of Ratio (OCCT L212-216).
        poids[(jj - 1) as usize] = ratio * second_curve.weights[(ii - 1) as usize];
        ii += 1;
        jj += 1;
    }

    // Create the BSpline (OCCT L219-220).
    let mut my_curve = BSplineCurve2 {
        degree: deg,
        knots: flat_knots_of(&noeuds, &mults),
        control_points: poles,
        weights: poids,
    };

    // Optionally reduce multiplicity (OCCT L222-229).
    let mut ok = true;
    let mut m = mults[(nb_k1 - 1) as usize];
    while m > 0 && ok {
        m -= 1;
        ok = bspline2_remove_knot(&mut my_curve, nb_k1, m, my_tol);
    }
    my_curve
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Line2d, TrimmedCurve2};

    /// Two trimmed lines meeting at a corner concatenate into one degree-1
    /// BSpline whose junction knot keeps multiplicity = Degree (the
    /// RemoveKnot reduction rejects the corner: deviation |(1,0) from the
    /// chord| = 1/sqrt(2) >> tol per BSplCLib::RemoveKnot anti-Boor check).
    ///
    /// Expected poles per CompCurveToBSplineCurve.cxx L203-217: the input
    /// poles concatenated (unit weights); expected knots per L187-200 with
    /// Ratio = |D1(end of L1)| / |D1(start of L2)| = 1/1 = 1 (unit line
    /// directions): Noeuds = [0, 1(U_de_raccord), 2] (Ratio2*K2(2) - Delta2
    /// = 1*1 + 1), Mults = [2, Degree, 2].
    ///
    /// Note: Line2d::new takes a gp_Dir2d-style direction (normalized), so
    /// the second line from (1,0) reaching (1,1) at t=1 carries the unit
    /// direction (0,1).
    #[test]
    fn concat_two_lines() {
        let l1 = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Line(Line2d::new(
                DVec2::new(0.0, 0.0),
                DVec2::new(1.0, 0.0),
            ))),
            t_min: 0.0,
            t_max: 1.0,
        });
        let l2 = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Line(Line2d::new(
                DVec2::new(1.0, 0.0),
                DVec2::new(0.0, 1.0),
            ))),
            t_min: 0.0,
            t_max: 1.0,
        });
        let mut comp =
            Geom2dConvertCompCurveToBSplineCurve::new(ConvertParameterisation::TgtThetaOver2);
        assert!(comp.add(&l1, 1.0e-7));
        assert!(comp.add(&l2, 1.0e-7));
        let bs = comp.bspline_curve().expect("merged");
        let merged: &BSplineCurve2 = &bs;
        assert_eq!(merged.degree, 1);
        assert_eq!(merged.control_points.len(), 3);
        assert!((merged.control_points[0] - DVec2::new(0.0, 0.0)).length() < 1.0e-12);
        assert!((merged.control_points[1] - DVec2::new(1.0, 0.0)).length() < 1.0e-12);
        assert!(
            (merged.control_points[2] - DVec2::new(1.0, 1.0)).length() < 1.0e-12,
            "end pole at the second line end"
        );
        // Knots: [0, 0, 1, 2, 2] (junction knot at multiplicity Degree — the
        // RemoveKnot loop failed on the corner and kept it).
        let want_knots = [0.0, 0.0, 1.0, 2.0, 2.0];
        assert_eq!(merged.knots.len(), want_knots.len());
        for (k, w) in merged.knots.iter().zip(want_knots.iter()) {
            assert!((k - w).abs() < 1.0e-12, "knot {k} vs {w}");
        }
        for w in merged.weights.iter() {
            assert!((w - 1.0).abs() < 1.0e-15, "unit weights");
        }
    }

    /// The trimmed-circle conversion reproduces the exact OCCT TgtThetaOver2
    /// construction (Geom2dConvert.cxx L228-264 + BSplineCurveBuilder L70-98):
    /// canonical numerator poles (R,0), (R,R), (0,R) — BuildCosAndSin with
    /// alpha = (U2-U1)/2 = PI/4 — placed on the circle frame; weights
    /// (1, cos(alpha), 1); clamped knots [0, PI/2].
    ///
    /// The TgtThetaOver2 parameterization is NOT the angle parameterization:
    /// at u = PI/8 (t = 1/4 in the rational Bezier) the exact evaluation is
    ///   x = 2 + (65 + 21*sqrt(2)) / (10 + 3*sqrt(2))
    ///   y = 3 + (35 + 24*sqrt(2)) / (10 + 3*sqrt(2))
    /// (W = (10+3*sqrt2)/16, xW = (65+21*sqrt2)/16, yW = (35+24*sqrt2)/16),
    /// and the point lies exactly on the circle (exact rational form).
    #[test]
    fn trimmed_circle_to_bspline_evaluates() {
        let sq2 = std::f64::consts::SQRT_2;
        let circ = Circle2d {
            center: DVec2::new(2.0, 3.0),
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: 5.0,
        };
        let trc = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Circle(circ)),
            t_min: 0.0,
            t_max: std::f64::consts::FRAC_PI_2,
        });
        let bs = curve_to_bspline_curve_2d(&trc, ConvertParameterisation::TgtThetaOver2);

        // The exact OCCT pole data.
        assert_eq!(bs.degree, 2);
        let want_poles = [
            DVec2::new(7.0, 3.0),
            DVec2::new(7.0, 8.0),
            DVec2::new(2.0, 8.0),
        ];
        assert_eq!(bs.control_points.len(), 3);
        for (p, w) in bs.control_points.iter().zip(want_poles.iter()) {
            assert!((p - w).length() < 1.0e-12, "pole {p:?} vs {w:?}");
        }
        let want_weights = [1.0, sq2 / 2.0, 1.0];
        for (w, want) in bs.weights.iter().zip(want_weights.iter()) {
            assert!((w - want).abs() < 1.0e-15);
        }
        let half_pi = std::f64::consts::FRAC_PI_2;
        let want_knots = [0.0, 0.0, 0.0, half_pi, half_pi, half_pi];
        assert_eq!(bs.knots.len(), want_knots.len());
        for (k, w) in bs.knots.iter().zip(want_knots.iter()) {
            assert!((k - w).abs() < 1.0e-15);
        }

        // Endpoints exact.
        assert!(Curve2dEval::point_at(&bs, 0.0).distance(DVec2::new(7.0, 3.0)) < 1.0e-12);
        assert!(Curve2dEval::point_at(&bs, half_pi).distance(DVec2::new(2.0, 8.0)) < 1.0e-12);

        // Exact rational evaluation at u = PI/8.  With t = 1/4 in the Bezier:
        // W = (10+3*sqrt2)/16, xW = (65+21*sqrt2)/16, yW = (35+24*sqrt2)/16
        // (the center is baked into the poles), plus on-circle exactness.
        let want_x = (65.0 + 21.0 * sq2) / (10.0 + 3.0 * sq2);
        let want_y = (35.0 + 24.0 * sq2) / (10.0 + 3.0 * sq2);
        let got = Curve2dEval::point_at(&bs, std::f64::consts::FRAC_PI_8);
        assert!(
            (got.x - want_x).abs() < 1.0e-12 && (got.y - want_y).abs() < 1.0e-12,
            "got {got:?}, want ({want_x}, {want_y})"
        );
        for k in 0..=8 {
            let t = half_pi * (k as f64) / 8.0;
            let p = Curve2dEval::point_at(&bs, t);
            let r = (p - DVec2::new(2.0, 3.0)).length();
            assert!((r - 5.0).abs() < 1.0e-12, "u={t}: radius deviation {}", r - 5.0);
        }
    }

    /// SplitBSplineCurve (parameter form): the segment is extracted on
    /// [min, max] (Geom2dConvert.cxx L154-159) and, for a non-periodic
    /// basis, reversed when FromU1 > ToU2 (L169-174; Geom2dConvert.hxx
    /// L66-68: "C is oriented fromU1 toU2").
    #[test]
    fn split_bspline_u_form() {
        let bs = BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::ZERO, DVec2::new(1.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        // FromU1 = 0.75 > ToU2 = 0.25: the result is oriented from 0.75 to
        // 0.25 (reversed).
        let seg = split_bspline_curve_2d(&bs, 0.75, 0.25, 1.0e-9, true);
        assert!((seg.control_points[0] - DVec2::new(0.75, 0.0)).length() < 1.0e-9);
        assert!((seg.control_points[1] - DVec2::new(0.25, 0.0)).length() < 1.0e-9);
        assert!((seg.knots[0] - 0.25).abs() < 1.0e-12 && (seg.knots[3] - 0.75).abs() < 1.0e-12);

        // Ascending form: FromU1 < ToU2 keeps the orientation.
        let seg_asc = split_bspline_curve_2d(&bs, 0.25, 0.75, 1.0e-9, true);
        assert!((seg_asc.control_points[0] - DVec2::new(0.25, 0.0)).length() < 1.0e-9);
        assert!((seg_asc.control_points[1] - DVec2::new(0.75, 0.0)).length() < 1.0e-9);
    }

    /// LEFT-HANDED frame (mirror placement) — the discriminating test for
    /// the BSplineCurveBuilder transform collapse.
    ///
    /// OCCT statement chain (Geom2dConvert.cxx L84-96): the left-handed
    /// conic first takes Sym.SetMirror(gp::OX2d()) — gp_Trsf2d.cxx L48-62:
    /// (x, y) -> (x, -y) — then T.SetTransformation(XAxis, OX2d) —
    /// gp_Trsf2d.cxx L64-84 with ToA2 = OX2d: matrix columns (XDir,
    /// perp(XDir)), translation Loc, column-vector application
    /// (gp_XY::Multiply = M*x).  The collapse is therefore
    ///   P = Loc + x*XDir + y*perp(XDir)
    /// and with the LEFT-HANDED frame's own YDir = -perp(XDir):
    ///   P = Loc + x*XDir + y*YDir  (no second y flip).
    ///
    /// Hand-derived expectation: canonical TgtThetaOver2 quarter arc of the
    /// radius-2 circle has numerator poles (2,0), (2,2), (0,2) with weights
    /// (1, sqrt(2)/2, 1).  With Loc = (5,0), X = (1,0), Y = (0,-1):
    ///   (2,0) -> (5,0) + (2,0)            = (7, 0)
    ///   (2,2) -> (5,0) + (2,0) + (0,-2)   = (7,-2)
    ///   (0,2) -> (5,0) + 0     + (0,-2)   = (5,-2)
    /// The pre-fix body (flipping the y-offset through YDir) produced
    /// (7,0), (7,2), (5,2) — the mirrored arc.
    #[test]
    fn left_handed_circle_placement_matches_occt_collapse() {
        let sq2 = std::f64::consts::SQRT_2;
        let circ = Circle2d {
            center: DVec2::new(5.0, 0.0),
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, -1.0), // left-handed frame.
            radius: 2.0,
        };
        let trc = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Circle(circ)),
            t_min: 0.0,
            t_max: std::f64::consts::FRAC_PI_2,
        });
        let bs = curve_to_bspline_curve_2d(&trc, ConvertParameterisation::TgtThetaOver2);

        assert_eq!(bs.degree, 2);
        let want_poles = [
            DVec2::new(7.0, 0.0),
            DVec2::new(7.0, -2.0),
            DVec2::new(5.0, -2.0),
        ];
        assert_eq!(bs.control_points.len(), 3);
        for (p, w) in bs.control_points.iter().zip(want_poles.iter()) {
            assert!((p - w).length() < 1.0e-12, "pole {p:?} vs {w:?}");
        }
        let want_weights = [1.0, sq2 / 2.0, 1.0];
        for (w, want) in bs.weights.iter().zip(want_weights.iter()) {
            assert!((w - want).abs() < 1.0e-15);
        }
        // Endpoints: Value(0) = Loc + R*XDir = (7,0); the trim end is the
        // collapse-mapped canonical (0, 2) = (5,-2).
        assert!(Curve2dEval::point_at(&bs, 0.0).distance(DVec2::new(7.0, 0.0)) < 1.0e-12);
        assert!(
            Curve2dEval::point_at(&bs, std::f64::consts::FRAC_PI_2).distance(DVec2::new(5.0, -2.0))
                < 1.0e-12
        );
        // The arc lies on the placed circle (exact rational form).
        for k in 0..=8 {
            let t = std::f64::consts::FRAC_PI_2 * (k as f64) / 8.0;
            let p = Curve2dEval::point_at(&bs, t);
            let r = (p - DVec2::new(5.0, 0.0)).length();
            assert!((r - 2.0).abs() < 1.0e-12, "u={t}: radius deviation {}", r - 2.0);
        }
    }

    /// The trimmed-Bezier branch (Geom2dConvert.cxx L323-345): Copy +
    /// Geom2d_BezierCurve::Segment.  For the quadratic (0,0),(2,4),(4,0)
    /// (x(u) = 4u, y(u) = 8u(1-u)) trimmed to [0.5, 1] (u = 1/2 + v/2):
    ///   x(v) = 2 + 2v      -> power coeffs (2, 2)
    ///   y(v) = 2 - 2v^2    -> power coeffs (2, 0, -2)
    /// and the PLib::CoefficientsPoles Pascal recombination over both
    /// passes turns the power form into the Bezier poles
    /// (2,2), (3,2), (4,0) over the clamped knots [0, 1], mults 3.
    #[test]
    fn trimmed_bezier_segment_poles() {
        use crate::geom::BezierCurve2 as Bez2;
        let cbez = Bez2 {
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(2.0, 4.0), DVec2::new(4.0, 0.0)],
            weights: vec![1.0, 1.0, 1.0],
        };
        let trc = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Bezier(cbez)),
            t_min: 0.5,
            t_max: 1.0,
        });
        let bs = curve_to_bspline_curve_2d(&trc, ConvertParameterisation::TgtThetaOver2);
        assert_eq!(bs.degree, 2);
        // OCCT L327-344: knots 0/1, mults Degree+1 = 3.
        assert_eq!(bs.knots, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let want_poles = [
            DVec2::new(2.0, 2.0),
            DVec2::new(3.0, 2.0),
            DVec2::new(4.0, 0.0),
        ];
        assert_eq!(bs.control_points.len(), 3);
        for (p, w) in bs.control_points.iter().zip(want_poles.iter()) {
            assert!((p - w).length() < 1.0e-12, "pole {p:?} vs {w:?}");
        }
        for w in bs.weights.iter() {
            assert!((w - 1.0).abs() < 1.0e-15);
        }
        // Endpoints exact; on-curve check at v = 1/2: u = 3/4,
        // x(3/4) = 3, y(3/4) = 8*(3/4)*(1/4) = 1.5.
        let mid = Curve2dEval::point_at(&bs, 0.5);
        assert!((mid - DVec2::new(3.0, 1.5)).length() < 1.0e-12, "mid {mid:?}");
    }

    /// The RATIONAL Segment arm (Geom2d_BezierCurve.cxx L360-375):
    /// homogeneous BuildCache -> PLib::Trimming with wcoeffs ->
    /// PLib::CoefficientsPoles with weights.  Quadratic (0,0),(2,4),(4,0)
    /// with weights (1,2,1) — homogeneous poles (0,0,1),(4,8,2),(4,0,1),
    /// Taylor coefficients at 0: (0,0,1),(8,16,2),(-4,-16,-2) — trimmed to
    /// [1/2, 1] (u = 1/2 + v/2) gives the power form
    ///   xW = 3 + 2v - v^2,  yW = 4 - 4v^2,  W = 1.5 - 0.5v^2
    /// and the Pascal recombination + weight normalization yields the poles
    /// (2, 8/3), (8/3, 8/3), (4, 0) with weights (1.5, 1.5, 1).
    #[test]
    fn trimmed_rational_bezier_segment_poles() {
        use crate::geom::BezierCurve2 as Bez2;
        let cbez = Bez2 {
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(2.0, 4.0), DVec2::new(4.0, 0.0)],
            weights: vec![1.0, 2.0, 1.0],
        };
        let trc = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Bezier(cbez)),
            t_min: 0.5,
            t_max: 1.0,
        });
        let bs = curve_to_bspline_curve_2d(&trc, ConvertParameterisation::TgtThetaOver2);
        assert_eq!(bs.degree, 2);
        assert_eq!(bs.knots, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let want_poles = [
            DVec2::new(2.0, 8.0 / 3.0),
            DVec2::new(8.0 / 3.0, 8.0 / 3.0),
            DVec2::new(4.0, 0.0),
        ];
        for (p, w) in bs.control_points.iter().zip(want_poles.iter()) {
            assert!((p - w).length() < 1.0e-12, "pole {p:?} vs {w:?}");
        }
        let want_weights = [1.5, 1.5, 1.0];
        for (w, want) in bs.weights.iter().zip(want_weights.iter()) {
            assert!((w - want).abs() < 1.0e-12, "weight {w} vs {want}");
        }
        // Endpoint values: y(u) = 16u(1-u) / (1 + 2u - 2u^2); y(1/2) = 8/3,
        // y(1) = 0.
        let p0 = Curve2dEval::point_at(&bs, 0.0);
        assert!((p0 - DVec2::new(2.0, 8.0 / 3.0)).length() < 1.0e-12);
        let p1 = Curve2dEval::point_at(&bs, 1.0);
        assert!((p1 - DVec2::new(4.0, 0.0)).length() < 1.0e-12);
    }
}
