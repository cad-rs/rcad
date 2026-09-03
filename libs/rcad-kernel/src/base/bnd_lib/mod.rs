//! Bounding-box computation for geometric entities (BndLib).
//!
//! OCCT TKGeomBase BndLib package: BndLib_Add3dCurve, BndLib_AddSurface.
//!
//! In OCCT, `Bnd_Box` (the data structure) lives in TKMath while `BndLib`
//! (algorithms that add `Geom_*` curves/surfaces to a `Bnd_Box`) lives in
//! TKGeomBase.  The functions below correspond to `BndLib_Add3dCurve` and
//! `BndLib_AddSurface`.
//!
//! OCCT references:
//!   BndLib → BndLib_Add3dCurve::Add(Adaptor3d_Curve, U1, U2, Tol, B)
//!         → GeomBndLib_Curve(C).Add(U1, U2, Tol, B)
//!         → per-type dispatch:
//!             Circle  → GeomBndLib_Circle::Box(gp_Circ, U1, U2, Tol)
//!             Ellipse → GeomBndLib_Ellipse::Box(gp_Elips, U1, U2, Tol)
//!             Line    → GeomBndLib_Line::Box(gp_Lin, U1, U2, Tol)
//!             BSpline → GeomBndLib_BSplineCurve::Box(U1, U2, Tol)
//!             Bezier  → GeomBndLib_BezierCurve::Box(U1, U2, Tol)
//!             Other   → GeomBndLib_OtherCurve::Box(U1, U2, Tol) — 33pt sampling
//!
//! Each per-type Box applies `Enlarge(Tol)` internally, matching OCCT.

use glam::{DVec2, DVec3};
use crate::core::precision::{PCONFUSION, parametric_default};
use crate::geom::{
    BSplineCurve3, BezierCurve3, Circle2d, Curve2d, Curve2dEval, Curve3, CurveEval, Ellipse2d,
    Line2d, Surface3,
};

// ---------------------------------------------------------------------------
// OCCT-aligned helpers
// ---------------------------------------------------------------------------

/// OCCT ElCLib::InPeriod — adjusts a parameter to the interval [theFirst, theLast).
fn in_period(par: f64, first: f64, last: f64) -> f64 {
    let period = last - first;
    if period <= 0.0 {
        return par;
    }
    (par - first).rem_euclid(period) + first
}

/// 3D point on a circle at parameter t.
/// OCCT: ElCLib::CircleValue(t, gp_Ax2(pos), R) = pos.Location + R*(cos(t)*XDir + sin(t)*YDir)
fn circle_point(center: DVec3, x_dir: DVec3, y_dir: DVec3, radius: f64, t: f64) -> DVec3 {
    center + radius * (t.cos() * x_dir + t.sin() * y_dir)
}

/// 3D point on an ellipse at parameter t.
/// OCCT: ElCLib::EllipseValue(t, gp_Ax2(pos), aMajR, aMinR)
fn ellipse_point(
    center: DVec3,
    major_dir: DVec3,
    minor_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
    t: f64,
) -> DVec3 {
    center + major_radius * t.cos() * major_dir + minor_radius * t.sin() * minor_dir
}

/// OCCT GeomBndLib_SplineHelpers::FillBox (also GeomBndLib_OtherCurve FillBox).
///
/// Samples [theFirst, theLast] with 2*theN+1 points, adds them to the running
/// bounding box (mn, mx), and returns the max deflection: the maximum distance
/// from the midpoint between two consecutive samples to the curve at the
/// mid-parameter.
fn fill_box_into<E: CurveEval + ?Sized>(
    curve: &E,
    first: f64,
    last: f64,
    n: usize,
    mn: &mut DVec3,
    mx: &mut DVec3,
) -> f64 {
    let p1 = curve.point_at(first);
    *mn = mn.min(p1);
    *mx = mx.max(p1);

    let mut max_tol: f64 = 0.0;
    let diff = last - first;
    if diff.abs() > PCONFUSION {
        let dp = diff / (2.0 * n as f64);
        let mut p = first;
        let mut a_p1 = p1;
        for _ in 1..=n {
            p += dp;
            let a_p2 = curve.point_at(p);
            *mn = mn.min(a_p2);
            *mx = mx.max(a_p2);
            p += dp;
            let a_p3 = curve.point_at(p);
            *mn = mn.min(a_p3);
            *mx = mx.max(a_p3);
            let a_pc = (a_p1 + a_p3) * 0.5;
            max_tol = max_tol.max(a_pc.distance(a_p2));
            a_p1 = a_p3;
        }
    } else {
        // OCCT degenerate branch: add the last point (== first point for a
        // zero-width span; OtherCurve re-adds first — functionally identical).
        let p3 = curve.point_at(last);
        *mn = mn.min(p3);
        *mx = mx.max(p3);
    }
    max_tol
}

/// OCCT GeomBndLib_SplineHelpers::ReduceSplineBox (non-indexed overload).
///
/// Intersects the sampled box with the convex hull of the poles:
///   reduced.min = max(sampled.min, poles.min) per coordinate
///   reduced.max = min(sampled.max, poles.max) per coordinate
/// If the poles box is void, the sampled box is copied unchanged.
fn reduce_spline_box(poles: &[DVec3], sampled: [DVec3; 2]) -> [DVec3; 2] {
    let mut p_mn = DVec3::splat(f64::INFINITY);
    let mut p_mx = DVec3::splat(f64::NEG_INFINITY);
    for &p in poles {
        p_mn = p_mn.min(p);
        p_mx = p_mx.max(p);
    }
    if !p_mn.is_finite() {
        return sampled;
    }
    [
        DVec3::new(
            sampled[0].x.max(p_mn.x),
            sampled[0].y.max(p_mn.y),
            sampled[0].z.max(p_mn.z),
        ),
        DVec3::new(
            sampled[1].x.min(p_mx.x),
            sampled[1].y.min(p_mx.y),
            sampled[1].z.min(p_mx.z),
        ),
    ]
}

/// Apply OCCT `Bnd_Box::Enlarge(Tol)` to a [min, max] pair (adds Tol to all faces).
fn enlarge_box(mn: DVec3, mx: DVec3, tol: f64) -> [DVec3; 2] {
    let en = DVec3::splat(tol);
    [mn - en, mx + en]
}

// ---------------------------------------------------------------------------
// Full curve bounding box (conservative, no range awareness)
// ---------------------------------------------------------------------------

/// Conservative bounding box for an analytic curve.
/// OCCT: BndLib_Add3dCurve::Add (full curve, no sub-range).
pub fn curve_bounding_box(curve: &Curve3) -> Option<[DVec3; 2]> {
    match curve {
        Curve3::Circle(c) => {
            // OCCT GeomBndLib_Circle::Box(gp_Circ, Tol) — full circle
            // Per-coordinate: amp_k = R * sqrt(x_dir_k^2 + y_dir_k^2)
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for k in 0..3 {
                let a_xk = c.x_dir[k];
                let a_yk = c.y_dir[k];
                let a_amp =
                    (c.radius * c.radius * a_xk * a_xk + c.radius * c.radius * a_yk * a_yk)
                        .sqrt();
                mn[k] = c.center[k] - a_amp;
                mx[k] = c.center[k] + a_amp;
            }
            Some([mn, mx])
        }
        Curve3::Ellipse(e) => {
            // OCCT GeomBndLib_Ellipse::Box(gp_Elips, Tol)
            let minor_dir = e.normal.cross(e.major_dir);
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for k in 0..3 {
                let a_xk = e.major_dir[k];
                let a_yk = minor_dir[k];
                let a_amp =
                    (e.major_radius * e.major_radius * a_xk * a_xk
                        + e.minor_radius * e.minor_radius * a_yk * a_yk)
                    .sqrt();
                mn[k] = e.center[k] - a_amp;
                mx[k] = e.center[k] + a_amp;
            }
            Some([mn, mx])
        }
        Curve3::Line(_) => None,
        Curve3::BSpline(b) => {
            let (mut mn, mut mx) =
                (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for &p in &b.control_points {
                mn = mn.min(p);
                mx = mx.max(p);
            }
            if mn.is_finite() {
                Some([mn, mx])
            } else {
                None
            }
        }
        Curve3::Bezier(b) => {
            let (mut mn, mut mx) =
                (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for &p in &b.control_points {
                mn = mn.min(p);
                mx = mx.max(p);
            }
            if mn.is_finite() {
                Some([mn, mx])
            } else {
                None
            }
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Range-aware bounding box (OCCT GeomBndLib per-type dispatch)
// ---------------------------------------------------------------------------

/// OCCT-aligned 1:1: Compute bounding box for a 3D curve over parameter range [u1, u2]
/// with tolerance `tol`.
///
/// Maps directly to GeomBndLib_Curve(C).Box(U1, U2, Tol). Every per-type Box
/// applies `Enlarge(Tol)` internally, exactly like OCCT.
pub fn curve_bounding_box_range(curve: &Curve3, u1: f64, u2: f64, tol: f64) -> Option<[DVec3; 2]> {
    match curve {
        Curve3::Circle(c) => circle_arc_box(c, u1, u2, tol),
        Curve3::Ellipse(e) => ellipse_arc_box(e, u1, u2, tol),
        Curve3::Line(l) => line_box(l, u1, u2, tol),
        Curve3::BSpline(b) => bspline_curve_box(b, u1, u2, tol),
        Curve3::Bezier(b) => bezier_curve_box(b, u1, u2, tol),
        // OCCT GeomBndLib_Curve unwraps Geom_TrimmedCurve and dispatches on the
        // basis curve type with the trimmed parameter range.
        Curve3::Trimmed(t) => curve_bounding_box_range(&t.curve, t.first, t.last, tol),
        // Remainder: OCCT GeomBndLib_OtherCurve::Box — 33-point sampling
        _ => other_curve_box(curve, u1, u2, tol),
    }
}

// ---------------------------------------------------------------------------
// Line box (OCCT GeomBndLib_Line::Box(gp_Lin, U1, U2, Tol))
// ---------------------------------------------------------------------------

/// OCCT GeomBndLib_Line::Box — endpoints + Enlarge(Tol).
///
/// Infinite parameters are handled by the edge's infinite-vertex creation in
/// prepare_edges; for an infinite line edge the box is subsumed by those
/// vertices (rcad BndBox cannot represent an open/infinite box — architecture
/// difference vs OCCT GeomBndLib_InfiniteHelpers::OpenMin/OpenMax).
fn line_box(l: &crate::geom::Line3, u1: f64, u2: f64, tol: f64) -> Option<[DVec3; 2]> {
    use crate::core::precision::{is_negative_infinite_value, is_positive_infinite_value};

    // Finite finite: add both endpoints.
    if !is_negative_infinite_value(u1) && !is_positive_infinite_value(u1)
        && !is_negative_infinite_value(u2) && !is_positive_infinite_value(u2)
    {
        let p1 = l.point_at(u1);
        let p2 = l.point_at(u2);
        let mn = p1.min(p2);
        let mx = p1.max(p2);
        return Some(enlarge_box(mn, mx, tol));
    }
    // Infinite range: return None — the infinite vertices created by
    // prepare_edges (point_at(±inf)) carry the box. OCCT opens the box in the
    // line direction (OpenMin/OpenMax); rcad relies on the vertex boxes.
    None
}

// ---------------------------------------------------------------------------
// Circle arc box (OCCT GeomBndLib_Circle::Box(gp_Circ, U1, U2, Tol))
// ---------------------------------------------------------------------------

/// Full circle bounding box (no range clipping).
fn circle_full_box(c: &crate::geom::Circle3) -> [DVec3; 2] {
    let mut mn = DVec3::splat(f64::INFINITY);
    let mut mx = DVec3::splat(f64::NEG_INFINITY);
    for k in 0..3 {
        let a_xk = c.x_dir[k];
        let a_yk = c.y_dir[k];
        let a_amp = (c.radius * c.radius * a_xk * a_xk + c.radius * c.radius * a_yk * a_yk).sqrt();
        mn[k] = c.center[k] - a_amp;
        mx[k] = c.center[k] + a_amp;
    }
    [mn, mx]
}

/// OCCT GeomBndLib_Circle::Box(gp_Circ, U1, U2, Tol) — arc-aware analytic bounding box.
fn circle_arc_box(c: &crate::geom::Circle3, u1: f64, u2: f64, tol: f64) -> Option<[DVec3; 2]> {
    use std::f64::consts::{PI, TAU};

    let period = TAU - PCONFUSION;

    // OCCT L58-61: if arc spans full period, return full circle box
    if u2 - u1 >= period {
        let [mn, mx] = circle_full_box(c);
        return Some(enlarge_box(mn, mx, tol));
    }

    // Normalize parameter range so u2' > u1' (handle wrapping through 0/2*PI)
    let a_u1 = u1.rem_euclid(TAU); // [0, 2*PI)
    let a_u2 = u2.rem_euclid(TAU); // [0, 2*PI)
    let a_u2 = if a_u2 > a_u1 {
        a_u2
    } else {
        a_u2 + TAU
    }; // u2' ∈ (u1', u1' + 2*PI)

    let mut mn = DVec3::splat(f64::INFINITY);
    let mut mx = DVec3::splat(f64::NEG_INFINITY);

    // OCCT L69-70: Add arc endpoints.
    let p1 = circle_point(c.center, c.x_dir, c.y_dir, c.radius, a_u1);
    let p2 = circle_point(c.center, c.x_dir, c.y_dir, c.radius, a_u2);
    mn = mn.min(p1).min(p2);
    mx = mx.max(p1).max(p2);

    // OCCT L72-109: For each coordinate, check if the extremal parameter lies within the arc.
    for k in 0..3 {
        let a_xk = c.x_dir[k];
        let a_yk = c.y_dir[k];

        // OCCT L79-84: extremal parameter t where d(P_k)/dt = 0
        // P_k(t) = center_k + R*cos(t)*x_dir_k + R*sin(t)*y_dir_k
        // dP_k/dt = -R*sin(t)*x_dir_k + R*cos(t)*y_dir_k = 0 => tan(t) = y_dir_k / x_dir_k
        // (OCCT uses gp::Resolution()=1e-15 here; PCONFUSION=1e-9 is a form
        //  difference — numerically equivalent for unit-length x_dir/y_dir.)
        let a_t_extr_min = if a_xk.abs() > PCONFUSION {
            let t = (a_yk / a_xk).atan();
            in_period(t, 0.0, TAU)
        } else {
            PI / 2.0
        };

        // OCCT L89: opposite extremal parameter (PI apart)
        let a_t_extr_max = if a_t_extr_min <= PI {
            a_t_extr_min + PI
        } else {
            a_t_extr_min - PI
        };

        // OCCT L91-97: evaluate and swap so a_t_extr_min gives the smaller value
        let a_val_min = c.radius * a_t_extr_min.cos() * a_xk
            + c.radius * a_t_extr_min.sin() * a_yk
            + c.center[k];
        let a_val_max = c.radius * a_t_extr_max.cos() * a_xk
            + c.radius * a_t_extr_max.sin() * a_yk
            + c.center[k];

        let (a_t_min, a_t_max) = if a_val_min > a_val_max {
            (a_t_extr_max, a_t_extr_min)
        } else {
            (a_t_extr_min, a_t_extr_max)
        };

        // OCCT L99-108: check min-extremal parameter in range
        let a_tk = in_period(a_t_min, a_u1, a_u1 + TAU);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            let p = circle_point(c.center, c.x_dir, c.y_dir, c.radius, a_t_min);
            mn = mn.min(p);
            mx = mx.max(p);
        }

        // OCCT L104-108: check max-extremal parameter in range
        let a_tk = in_period(a_t_max, a_u1, a_u1 + TAU);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            let p = circle_point(c.center, c.x_dir, c.y_dir, c.radius, a_t_max);
            mn = mn.min(p);
            mx = mx.max(p);
        }
    }

    if mn.is_finite() {
        Some(enlarge_box(mn, mx, tol))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Ellipse arc box (OCCT GeomBndLib_Ellipse::Box(gp_Elips, U1, U2, Tol))
// ---------------------------------------------------------------------------

/// OCCT GeomBndLib_Ellipse::Box(gp_Elips, U1, U2, Tol) — arc-aware analytic bounding box.
fn ellipse_arc_box(e: &crate::geom::Ellipse3, u1: f64, u2: f64, tol: f64) -> Option<[DVec3; 2]> {
    use std::f64::consts::{PI, TAU};

    let period = TAU - PCONFUSION;

    let minor_dir = e.normal.cross(e.major_dir);

    // Full ellipse check
    if u2 - u1 >= period {
        let mut mn = DVec3::splat(f64::INFINITY);
        let mut mx = DVec3::splat(f64::NEG_INFINITY);
        for k in 0..3 {
            let a_xk = e.major_dir[k];
            let a_yk = minor_dir[k];
            let a_amp = (e.major_radius * e.major_radius * a_xk * a_xk
                + e.minor_radius * e.minor_radius * a_yk * a_yk)
                .sqrt();
            mn[k] = e.center[k] - a_amp;
            mx[k] = e.center[k] + a_amp;
        }
        return Some(enlarge_box(mn, mx, tol));
    }

    // Normalize parameter range
    let a_u1 = u1.rem_euclid(TAU);
    let a_u2 = u2.rem_euclid(TAU);
    let a_u2 = if a_u2 > a_u1 {
        a_u2
    } else {
        a_u2 + TAU
    };

    let mut mn = DVec3::splat(f64::INFINITY);
    let mut mx = DVec3::splat(f64::NEG_INFINITY);

    // OCCT L70-72: Add arc endpoints
    let p1 = ellipse_point(e.center, e.major_dir, minor_dir, e.major_radius, e.minor_radius, a_u1);
    let p2 = ellipse_point(e.center, e.major_dir, minor_dir, e.major_radius, e.minor_radius, a_u2);
    mn = mn.min(p1).min(p2);
    mx = mx.max(p1).max(p2);

    // OCCT L74-111: per-coordinate extremal parameter check
    for k in 0..3 {
        let a_xk = e.major_dir[k];
        let a_yk = minor_dir[k];

        // OCCT L80-84: Ellipse extremal: tan(t) = (aMinR * Yk) / (aMajR * Xk)
        let a_t_extr_min = if a_xk.abs() > PCONFUSION {
            let t = ((e.minor_radius * a_yk) / (e.major_radius * a_xk)).atan();
            in_period(t, 0.0, TAU)
        } else {
            PI / 2.0
        };

        // OCCT L90: opposite extremal parameter
        let a_t_extr_max = if a_t_extr_min <= PI {
            a_t_extr_min + PI
        } else {
            a_t_extr_min - PI
        };

        // OCCT L92-98: evaluate and swap
        let a_val_min = e.major_radius * a_t_extr_min.cos() * a_xk
            + e.minor_radius * a_t_extr_min.sin() * a_yk
            + e.center[k];
        let a_val_max = e.major_radius * a_t_extr_max.cos() * a_xk
            + e.minor_radius * a_t_extr_max.sin() * a_yk
            + e.center[k];

        let (a_t_min, a_t_max) = if a_val_min > a_val_max {
            (a_t_extr_max, a_t_extr_min)
        } else {
            (a_t_extr_min, a_t_extr_max)
        };

        // OCCT L101-103: check min extremal in range
        let a_tk = in_period(a_t_min, a_u1, a_u1 + TAU);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            let p = ellipse_point(
                e.center, e.major_dir, minor_dir, e.major_radius, e.minor_radius, a_t_min,
            );
            mn = mn.min(p);
            mx = mx.max(p);
        }

        // OCCT L106-108: check max extremal in range
        let a_tk = in_period(a_t_max, a_u1, a_u1 + TAU);
        if a_tk >= a_u1 && a_tk <= a_u2 {
            let p = ellipse_point(
                e.center, e.major_dir, minor_dir, e.major_radius, e.minor_radius, a_t_max,
            );
            mn = mn.min(p);
            mx = mx.max(p);
        }
    }

    if mn.is_finite() {
        Some(enlarge_box(mn, mx, tol))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Bezier box (OCCT GeomBndLib_BezierCurve::Box(U1, U2, Tol))
// ---------------------------------------------------------------------------

/// OCCT GeomBndLib_BezierCurve::Box(U1, U2, Tol).
///
///   aDeflection = FillBox(aSampledBox, aGACurve, U1, U2, Degree);
///   aSampledBox.Enlarge(aWeakness * aDeflection);
///   ReduceSplineBox(Poles, aSampledBox, aBox);
///   aBox.Enlarge(theTol);
fn bezier_curve_box(b: &BezierCurve3, u1: f64, u2: f64, tol: f64) -> Option<[DVec3; 2]> {
    const WEAKNESS: f64 = 1.5; // OCCT GeomBndLib_BezierCurve.cxx L31

    // OCCT: myGeom->Degree() = NbPoles - 1
    let degree = b.control_points.len().saturating_sub(1);

    let mut smn = DVec3::splat(f64::INFINITY);
    let mut smx = DVec3::splat(f64::NEG_INFINITY);
    let a_deflection = fill_box_into(b, u1, u2, degree, &mut smn, &mut smx);

    if !smn.is_finite() {
        // OCCT: aSampledBox.IsVoid() → empty box
        return None;
    }

    // OCCT L40: aSampledBox.Enlarge(aWeakness * aDeflection);
    let en = DVec3::splat(WEAKNESS * a_deflection);
    smn -= en;
    smx += en;

    // OCCT L43: ReduceSplineBox(myGeom->Poles(), aSampledBox, aBox);
    let reduced = reduce_spline_box(&b.control_points, [smn, smx]);

    // OCCT L45: aBox.Enlarge(theTol);
    Some(enlarge_box(reduced[0], reduced[1], tol))
}

// ---------------------------------------------------------------------------
// BSpline box (OCCT GeomBndLib_BSplineCurve::Box(U1, U2, Tol))
// ---------------------------------------------------------------------------

/// Distinct (compressed) knot values of an expanded knot vector.
fn unique_knots(knots: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        match out.last() {
            Some(&last) if (k - last).abs() < 1e-15 => {}
            _ => out.push(k),
        }
    }
    out
}

/// OCCT GeomBndLib_BSplineCurve::Box(U1, U2, Tol) — knot-span FillBox + ReduceSplineBox.
///
///   1. If [U1,U2] is a strict sub-range, the curve is segmented via real knot
///      insertion (Boehm, BSplCLib::Segment) with aSegmentTol.
///   2. Walk the (possibly segmented) curve's knot spans, FillBox-sampling each.
///   3. Enlarge by weakness*maxDeflection, then intersect with the poles hull.
///   4. Enlarge by theTol.
///
/// Periodic curves (is_periodic == true) use the OCCT periodic branch:
/// AdjustPeriodic normalizes [U1,U2] into one period, and the effective domain
/// is [Knots(FirstUKnotIndex), Knots(LastUKnotIndex)].
fn bspline_curve_box(b: &BSplineCurve3, u1: f64, u2: f64, tol: f64) -> Option<[DVec3; 2]> {
    const WEAKNESS: f64 = 1.5; // OCCT GeomBndLib_BSplineCurve.cxx L33
    use crate::math::bspl::{first_uknot_index, last_uknot_index, segment_bspline_curve};

    let degree = b.degree;
    let n = b.knots.len();
    if n < 2 * degree + 2 || degree == 0 {
        return None;
    }

    // OCCT: aCurve->FirstParameter()/LastParameter() = Knots(FirstUKnotIndex) /
    // Knots(LastUKnotIndex). Derived from the knot multiplicities — correct for
    // both clamped and periodic knot vectors (OCCT BSplCLib).
    let uk0 = unique_knots(&b.knots);
    let fk = first_uknot_index(&b.knots, degree);
    let lk = last_uknot_index(&b.knots, degree);
    let first_param = uk0[fk];
    let last_param = uk0[lk];

    // OCCT L38-39: sub-range check. Precision::Parametric(theTol) = theTol * 0.01.
    let parametric_tol = parametric_default(tol);

    // Curve used for knot-span walking — the original, or the segmented sub-curve.
    let walk_curve: BSplineCurve3;

    if (first_param - u1).abs() > parametric_tol || (last_param - u2).abs() > parametric_tol {
        // OCCT L41-82: Copy + Segment.
        let mut au1 = u1;
        let mut au2 = u2;
        let period = last_param - first_param;

        if b.is_periodic {
            // OCCT L43-50: ElCLib::AdjustPeriodic(First, Last, PConfusion, U1, U2).
            au1 = in_period(au1, first_param, first_param + period);
            if (au2 - au1).abs() > period * 0.5 {
                au2 = in_period(au2, first_param, first_param + period);
            } else {
                au2 = in_period(au2, au1, au1 + period);
            }
        } else {
            // OCCT L53-56: clamp to curve bounds.
            if first_param > au1 {
                au1 = first_param;
            }
            if last_param < au2 {
                au2 = last_param;
            }
        }

        // OCCT L63-79: aSegmentTol = 2.0 * PConfusion();
        let mut a_segment_tol = 2.0 * PCONFUSION;
        if b.is_periodic {
            let a_direct_diff = (au2 - au1).abs();
            let a_cross1 = (au2 - period - au1).abs();
            let a_cross2 = (au1 - period - au2).abs();
            let a_min_diff = a_direct_diff.min(a_cross1.min(a_cross2));
            if a_min_diff < a_segment_tol {
                a_segment_tol = a_min_diff * 0.01;
            }
        } else if (au2 - au1).abs() < a_segment_tol {
            a_segment_tol = (au2 - au1).abs() * 0.01;
        }

        // OCCT L81: aCurve = aSegmentedCurve->Segment(au1, au2, aSegmentTol).
        let seg = segment_bspline_curve(b, au1, au2, a_segment_tol)?;
        walk_curve = seg;
    } else {
        // Full range — walk the original curve.
        walk_curve = b.clone();
    }

    // OCCT L85-104: walk the (possibly segmented) curve's knot spans.
    // A segmented curve is clamped (walk all spans); a periodic full-range
    // curve is still periodic (walk only the interior spans between
    // FirstUKnotIndex and LastUKnotIndex, excluding the seam spans).
    let uk = unique_knots(&walk_curve.knots);
    let (span_lo, span_hi) = if walk_curve.is_periodic {
        (fk, lk)
    } else {
        (0, uk.len() - 1)
    };

    let mut smn = DVec3::splat(f64::INFINITY);
    let mut smx = DVec3::splat(f64::NEG_INFINITY);
    let mut a_max_deflection: f64 = 0.0;

    // OCCT L96-102: for (aKnot = aKnotFirst+1; aKnot <= aKnotLast; ++aKnot)
    let mut a_first = uk[span_lo];
    for i in (span_lo + 1)..=span_hi {
        let a_last = uk[i];
        let defl = fill_box_into(&walk_curve, a_first, a_last, degree, &mut smn, &mut smx);
        a_max_deflection = a_max_deflection.max(defl);
        a_first = a_last;
    }

    if !smn.is_finite() {
        // OCCT: aSampledBox.IsVoid() → empty box
        return None;
    }

    // OCCT L109: aSampledBox.Enlarge(aWeakness * aMaxDeflection);
    let en = DVec3::splat(WEAKNESS * a_max_deflection);
    smn -= en;
    smx += en;

    // OCCT L110: ReduceSplineBox(myGeom->Poles(), aSampledBox, aBox);
    // NOTE: uses the ORIGINAL curve's poles, not the segmented curve's.
    let reduced = reduce_spline_box(&b.control_points, [smn, smx]);

    // OCCT L111: aBox.Enlarge(theTol);
    Some(enlarge_box(reduced[0], reduced[1], tol))
}

// ---------------------------------------------------------------------------
// OtherCurve box (OCCT GeomBndLib_OtherCurve::Box — 33pt sampling)
// ---------------------------------------------------------------------------

/// OCCT GeomBndLib_OtherCurve::Box(U1, U2, Tol) — 33-point sampling with deflection.
///
/// OCCT FillBox + Enlarge(weakness * max_deflection) + Enlarge(Tol).
fn other_curve_box(curve: &Curve3, u1: f64, u2: f64, tol: f64) -> Option<[DVec3; 2]> {
    const WEAKNESS: f64 = 1.5; // OCCT GeomBndLib_OtherCurve.cxx L75
    const N: usize = 33; // OCCT GeomBndLib_OtherCurve.cxx L76

    let mut mn = DVec3::splat(f64::INFINITY);
    let mut mx = DVec3::splat(f64::NEG_INFINITY);
    let max_tol = fill_box_into(curve, u1, u2, N, &mut mn, &mut mx);

    if !mn.is_finite() {
        return None;
    }

    // OCCT L79: aB1.Enlarge(weakness * tol);
    let en = DVec3::splat(WEAKNESS * max_tol);
    mn -= en;
    mx += en;

    // OCCT L83: aBox.Enlarge(theTol);
    Some(enlarge_box(mn, mx, tol))
}

// ---------------------------------------------------------------------------
// Surface bounding box (unchanged)
// ---------------------------------------------------------------------------

/// Conservative bounding box for an analytic surface.
/// OCCT: BndLib_AddSurface::Add.
pub fn surface_bounding_box(
    surface: &Surface3,
    vertices: &[crate::topo::topology::Vertex],
) -> Option<[DVec3; 2]> {
    match surface {
        Surface3::Cylinder(cyl) => {
            let r = cyl.radius;
            let axis = cyl.axis.normalize_or_zero();
            if axis.length_squared() < 0.5 {
                return None;
            }
            let (mut min_axial, mut max_axial) = (f64::INFINITY, f64::NEG_INFINITY);
            for v in vertices {
                let proj = (v.point - cyl.origin).dot(axis);
                min_axial = min_axial.min(proj);
                max_axial = max_axial.max(proj);
            }
            if !min_axial.is_finite() {
                return None;
            }
            let p_lo = cyl.origin + axis * min_axial;
            let p_hi = cyl.origin + axis * max_axial;
            // OCCT BndLib_AddSurface::Add(CylindricalSurface) samples the
            // surface grid: the axial extent is bounded by the v-isoline
            // circles (the face vertices), only the radial (perpendicular to
            // the axis) extent grows by the radius.
            let radial = DVec3::new(
                r * (1.0 - axis.x * axis.x).max(0.0).sqrt(),
                r * (1.0 - axis.y * axis.y).max(0.0).sqrt(),
                r * (1.0 - axis.z * axis.z).max(0.0).sqrt(),
            );
            Some([p_lo.min(p_hi) - radial, p_lo.max(p_hi) + radial])
        }
        Surface3::Sphere(sph) => Some(
            [sph.center - DVec3::splat(sph.radius), sph.center + DVec3::splat(sph.radius)],
        ),
        Surface3::Torus(tor) => {
            let r = DVec3::splat(tor.major_radius + tor.minor_radius);
            Some([tor.center - r, tor.center + r])
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis_dir();
            let apex = cone.apex_point();
            let (mut min_axial, mut max_axial) = (f64::INFINITY, f64::NEG_INFINITY);
            for v in vertices {
                let proj = (v.point - apex).dot(axis);
                min_axial = min_axial.min(proj);
                max_axial = max_axial.max(proj);
            }
            if !min_axial.is_finite() {
                return None;
            }
            let max_r = cone.radius_at_axial(min_axial).max(cone.radius_at_axial(max_axial));
            let r_eff = max_r.max(cone.radius);
            let p_lo = apex + axis * min_axial;
            let p_hi = apex + axis * max_axial;
            // OCCT BndLib_AddSurface::Add(Cone): the axial extent is bounded
            // by the v-isoline circles at the face vertices, only the radial
            // (perpendicular to the axis) extent grows by the local radius.
            let radial = DVec3::new(
                r_eff * (1.0 - axis.x * axis.x).max(0.0).sqrt(),
                r_eff * (1.0 - axis.y * axis.y).max(0.0).sqrt(),
                r_eff * (1.0 - axis.z * axis.z).max(0.0).sqrt(),
            );
            Some([p_lo.min(p_hi) - radial, p_lo.max(p_hi) + radial])
        }
        Surface3::Ellipsoid(e) => {
            let max_r = e.radius_x.max(e.radius_y).max(e.radius_z);
            Some([e.center - DVec3::splat(max_r), e.center + DVec3::splat(max_r)])
        }
        Surface3::BSpline(b) => {
            let (mut mn, mut mx) =
                (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for row in &b.control_points {
                for p in row {
                    mn = mn.min(*p);
                    mx = mx.max(*p);
                }
            }
            if mn.is_finite() {
                Some([mn, mx])
            } else {
                None
            }
        }
        Surface3::Bezier(b) => {
            let (mut mn, mut mx) =
                (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
            for row in &b.control_points {
                for &p in row {
                    mn = mn.min(p);
                    mx = mx.max(p);
                }
            }
            if mn.is_finite() {
                Some([mn, mx])
            } else {
                None
            }
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// 2D curve bounding box — OCCT GeomBndLib_Curve2d::Box (used by
// BRepTools::AddUVBounds through BndLib_Add2dCurve::Add)
// ---------------------------------------------------------------------------

/// Add a 2D point into a [u_min, u_max, v_min, v_max] box.
fn box2d_add_point(b: &mut [f64; 4], p: DVec2) {
    b[0] = b[0].min(p.x);
    b[1] = b[1].max(p.x);
    b[2] = b[2].min(p.y);
    b[3] = b[3].max(p.y);
}

/// OCCT GeomBndLib_Line2d::Box (GeomBndLib_Line2d.hxx L70-127) — segment
/// endpoints. Infinite parameters open the corresponding side of the box
/// (GeomBndLib_InfiniteHelpers OpenMin/OpenMax); rcad uses the fully open box
/// as a conservative superset.
fn line2d_box_uv(l: &Line2d, u1: f64, u2: f64, tol: f64) -> [f64; 4] {
    let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    if !u1.is_infinite() {
        box2d_add_point(&mut b, l.origin + l.direction * u1);
    }
    if !u2.is_infinite() {
        box2d_add_point(&mut b, l.origin + l.direction * u2);
    }
    // OpenMin/OpenMax: the missing endpoint leaves the box open on that side.
    b[0] -= tol;
    b[1] += tol;
    b[2] -= tol;
    b[3] += tol;
    b
}

/// OCCT GeomBndLib_Circle2d::Box (GeomBndLib_Circle2d.cxx L23-105) — exact
/// bounds of a circular arc; a full circle uses the analytical extrema.
fn circle2d_box_uv(c: &Circle2d, u1: f64, u2: f64, tol: f64) -> [f64; 4] {
    let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    let r = c.radius;
    let o = c.center;
    let xd = c.x_dir;
    let yd = c.y_dir;
    let eval = |t: f64| o + r * t.cos() * xd + r * t.sin() * yd;
    let two_pi = std::f64::consts::TAU;
    if u2 - u1 >= two_pi - PCONFUSION {
        // OCCT L23-44: full circle — analytical extrema per coordinate.
        for k in 0..2 {
            let a_xk = if k == 0 { xd.x } else { xd.y };
            let a_yk = if k == 0 { yd.x } else { yd.y };
            let a_amp = (r * r * a_xk * a_xk + r * r * a_yk * a_yk).sqrt();
            if k == 0 {
                b[0] = o.x - a_amp;
                b[1] = o.x + a_amp;
            } else {
                b[2] = o.y - a_amp;
                b[3] = o.y + a_amp;
            }
        }
    } else {
        // OCCT L48-100: arc endpoints + in-range extrema.
        box2d_add_point(&mut b, eval(u1));
        box2d_add_point(&mut b, eval(u2));
        for k in 0..2 {
            let a_xk = if k == 0 { xd.x } else { xd.y };
            let a_yk = if k == 0 { yd.x } else { yd.y };
            // OCCT: gp::Resolution() == 1e-15.
            let a_t_extr_min = if a_xk.abs() > 1e-15 {
                in_period((a_yk / a_xk).atan(), 0.0, two_pi)
            } else {
                std::f64::consts::FRAC_PI_2
            };
            let a_t_extr_max = if a_t_extr_min <= std::f64::consts::PI {
                a_t_extr_min + std::f64::consts::PI
            } else {
                a_t_extr_min - std::f64::consts::PI
            };
            let a_tk = in_period(a_t_extr_min, u1, u1 + two_pi);
            if a_tk >= u1 && a_tk <= u2 {
                box2d_add_point(&mut b, eval(a_t_extr_min));
            }
            let a_tk = in_period(a_t_extr_max, u1, u1 + two_pi);
            if a_tk >= u1 && a_tk <= u2 {
                box2d_add_point(&mut b, eval(a_t_extr_max));
            }
        }
    }
    // OCCT: aBox.Enlarge(theTol).
    b[0] -= tol;
    b[1] += tol;
    b[2] -= tol;
    b[3] += tol;
    b
}

/// OCCT GeomBndLib_Ellipse2d::Box (GeomBndLib_Ellipse2d.cxx L23-105) — exact
/// bounds of an elliptical arc; a full ellipse uses the analytical extrema.
fn ellipse2d_box_uv(c: &Ellipse2d, u1: f64, u2: f64, tol: f64) -> [f64; 4] {
    let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    let a_maj_r = c.major_radius;
    let a_min_r = c.minor_radius;
    let o = c.center;
    let xd = c.major_dir;
    // OCCT gp_Elips2d: YAxis = Rot90(XAxis).
    let yd = DVec2::new(-xd.y, xd.x);
    let eval = |t: f64| o + a_maj_r * t.cos() * xd + a_min_r * t.sin() * yd;
    let two_pi = std::f64::consts::TAU;
    if u2 - u1 >= two_pi - PCONFUSION {
        // OCCT L23-44: full ellipse — analytical extrema per coordinate.
        for k in 0..2 {
            let a_xk = if k == 0 { xd.x } else { xd.y };
            let a_yk = if k == 0 { yd.x } else { yd.y };
            let a_amp = (a_maj_r * a_maj_r * a_xk * a_xk + a_min_r * a_min_r * a_yk * a_yk).sqrt();
            if k == 0 {
                b[0] = o.x - a_amp;
                b[1] = o.x + a_amp;
            } else {
                b[2] = o.y - a_amp;
                b[3] = o.y + a_amp;
            }
        }
    } else {
        // OCCT L48-100: arc endpoints + in-range extrema.
        box2d_add_point(&mut b, eval(u1));
        box2d_add_point(&mut b, eval(u2));
        for k in 0..2 {
            let a_xk = if k == 0 { xd.x } else { xd.y };
            let a_yk = if k == 0 { yd.x } else { yd.y };
            let a_t_extr_min = if a_xk.abs() > 1e-15 {
                in_period(((a_min_r * a_yk) / (a_maj_r * a_xk)).atan(), 0.0, two_pi)
            } else {
                std::f64::consts::FRAC_PI_2
            };
            let a_t_extr_max = if a_t_extr_min <= std::f64::consts::PI {
                a_t_extr_min + std::f64::consts::PI
            } else {
                a_t_extr_min - std::f64::consts::PI
            };
            let a_tk = in_period(a_t_extr_min, u1, u1 + two_pi);
            if a_tk >= u1 && a_tk <= u2 {
                box2d_add_point(&mut b, eval(a_t_extr_min));
            }
            let a_tk = in_period(a_t_extr_max, u1, u1 + two_pi);
            if a_tk >= u1 && a_tk <= u2 {
                box2d_add_point(&mut b, eval(a_t_extr_max));
            }
        }
    }
    // OCCT: aBox.Enlarge(theTol).
    b[0] -= tol;
    b[1] += tol;
    b[2] -= tol;
    b[3] += tol;
    b
}

/// OCCT GeomBndLib_OtherCurve2d::Box (GeomBndLib_OtherCurve2d.cxx L213-228) —
/// 33-point sampling with deflection-based enlargement
/// (Enlarge(1.5 * deflection), then Enlarge(theTol)).
fn other_curve2d_box_uv(c: &Curve2d, u1: f64, u2: f64, tol: f64) -> [f64; 4] {
    const N: usize = 33;
    const WEAKNESS: f64 = 1.5;
    let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    let p1 = c.point_at(u1);
    box2d_add_point(&mut b, p1);
    let mut max_tol: f64 = 0.0;
    let diff = u2 - u1;
    if diff.abs() > PCONFUSION {
        let dp = diff / (2.0 * N as f64);
        let mut p = u1;
        let mut a_p1 = p1;
        for _ in 1..=N {
            p += dp;
            let a_p2 = c.point_at(p);
            box2d_add_point(&mut b, a_p2);
            p += dp;
            let a_p3 = c.point_at(p);
            box2d_add_point(&mut b, a_p3);
            let a_pc = (a_p1 + a_p3) * 0.5;
            max_tol = max_tol.max(a_pc.distance(a_p2));
            a_p1 = a_p3;
        }
    } else {
        // OCCT degenerate branch: add the last point.
        box2d_add_point(&mut b, c.point_at(u2));
    }
    let en = WEAKNESS * max_tol + tol;
    b[0] -= en;
    b[1] += en;
    b[2] -= en;
    b[3] += en;
    b
}

/// OCCT GeomBndLib_Curve2d::Box (GeomBndLib_Curve2d.cxx L223-238) — exact 2D
/// bounding box of a curve sub-range [u1, u2] (the box already enlarged by
/// `tol`). Returns [u_min, u_max, v_min, v_max]. Used by
/// BRepTools::AddUVBounds via BndLib_Add2dCurve::Add (BRepTools.cxx L185).
pub fn curve2d_bounding_box(c: &Curve2d, u1: f64, u2: f64, tol: f64) -> [f64; 4] {
    match c {
        Curve2d::Line(l) => line2d_box_uv(l, u1, u2, tol),
        Curve2d::Circle(cir) => circle2d_box_uv(cir, u1, u2, tol),
        Curve2d::Ellipse(el) => ellipse2d_box_uv(el, u1, u2, tol),
        _ => other_curve2d_box_uv(c, u1, u2, tol),
    }
}

/// Function-version of [`curve_bounding_box_range`] (OCCT
/// `BndLib_Add3dCurve::Add(Adaptor3d_Curve, U1, U2, Tol, Box)` with an
/// arbitrary point-evaluation adaptor, e.g. GeomFill_SnglrFunc used as a
/// curve): sample `eval` over [u1, u2], pad by `tol`, and return the box
/// together with the gap (= `tol`, as `Bnd_Box::GetGap()`).
pub fn curve_box_range_fn(
    eval: &dyn Fn(f64) -> DVec3,
    u1: f64,
    u2: f64,
    tol: f64,
) -> ([DVec3; 2], f64) {
    let mut mn = DVec3::splat(f64::INFINITY);
    let mut mx = DVec3::splat(f64::NEG_INFINITY);
    const N_GRID: usize = 64;
    for i in 0..=N_GRID {
        let u = u1 + (u2 - u1) * (i as f64) / (N_GRID as f64);
        let p = eval(u);
        mn = mn.min(p);
        mx = mx.max(p);
    }
    mn -= DVec3::splat(tol);
    mx += DVec3::splat(tol);
    ([mn, mx], tol)
}
