//! The rcad re-host of `Geom_Curve::EvalDN(U, N)` over the [`Curve3`] value —
//! the union of the per-type bodies.  It is the `basisCurve->EvalDN(...)`
//! leaf of the swept-surface derivative translations
//! (`Geom_SurfaceOfLinearExtrusion::EvalDN` L245-271,
//! `Geom_SurfaceOfRevolution::EvalDN` L331-360) and of the
//! `Geom_ExtrusionUtils::DN` / `Geom_RevolutionUtils::DN` helpers
//! (`GeomAdaptor_Surface::EvalDN` L1754-1780).
//!
//! Per-type anchors:
//!   - `Geom_Line::EvalDN` (Geom_Line.cxx L207-221) / `ElCLib::LineDN`
//!     (ElCLib.cxx L911-918);
//!   - `Geom_Circle::EvalDN` (Geom_Circle.cxx L184-191) / `ElCLib::CircleDN`
//!     (ElCLib.cxx L922-955);
//!   - `Geom_Ellipse::EvalDN` (Geom_Ellipse.cxx L230-237) / `ElCLib::EllipseDN`
//!     (ElCLib.cxx L957-994);
//!   - `Geom_Hyperbola::EvalDN` (Geom_Hyperbola.cxx L269-276) /
//!     `ElCLib::HyperbolaDN` (ElCLib.cxx L996-1018);
//!   - `Geom_Parabola::EvalDN` (Geom_Parabola.cxx L205-212) /
//!     `ElCLib::ParabolaDN` (ElCLib.cxx L1020-1047);
//!   - `Geom_BSplineCurve::EvalDN` (Geom_BSplineCurve_1.cxx L300-316) —
//!     `BSplCLib::DN`, the rcad engine already hosted by
//!     [`BSplineCurve3::dn`];
//!   - `Geom_BezierCurve::EvalDN` (Geom_BezierCurve.cxx L601-620) — the same
//!     `BSplCLib::DN` over the Bezier's own knots/multiplicities;
//!   - `Geom_TrimmedCurve::EvalDN` (Geom_TrimmedCurve.cxx L240-243) — the
//!     basis curve.
//!
//! Architecture differences: `ElCLib::DN(U, gp_Circ/gp_Elips/gp_Hypr/gp_Parab)`
//! reads the `gp_Ax2` frame of the conic while the rcad payloads carry the
//! frame axes directly, and the rcad `Parabola3::focal_param` is the OCCT
//! `p = 2 * Focal()` (see `geom/mod.rs` L264/L1147), so the `U / (2*Focal)` and
//! `1 / (2*Focal)` coefficients of `ElCLib::ParabolaDN` are `U / focal_param`
//! and `1 / focal_param` here.
//!
//! The rcad-only curve kinds (`OffsetCurve` / `CircularHelix` / `SineWave`)
//! have no translated `Geom_Curve::EvalDN` counterpart yet: they raise an
//! explicit gap rather than a silent zero.

use glam::DVec3;

use crate::geom::Curve3;

/// OCCT `gp::Resolution()` = `RealSmall()` = `DBL_MIN`
/// (gp.hxx L59-60, Standard_Real.hxx L132-135) — the `Focal` degenerate guard
/// of `ElCLib::ParabolaDN` (ElCLib.cxx L1032/L1043).
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT `Geom_Line::EvalDN(U, N)` (Geom_Line.cxx L207-221) — the line's
/// direction for `N == 1`, the null vector above.
///
/// `gp_Lin::Direction()` is a unit `gp_Dir` by construction; the rcad payload
/// carries the direction vector as stored (`Line3::point_at` uses it raw too),
/// so the two agree on the unit-direction invariant every rcad line payload
/// holds.
fn line_dn(l: &crate::geom::Line3, n: i32) -> DVec3 {
    if n == 1 {
        l.direction
    } else {
        DVec3::ZERO
    }
}

/// OCCT `ElCLib::CircleDN(U, Pos, Radius, N)` (ElCLib.cxx L922-955).
fn circle_dn(c: &crate::geom::Circle3, u: f64, n: i32) -> DVec3 {
    let mut xc = 0.0;
    let mut yc = 0.0;
    if n == 1 {
        xc = c.radius * -u.sin();
        yc = c.radius * u.cos();
    } else if (n + 2) % 4 == 0 {
        xc = c.radius * -u.cos();
        yc = c.radius * -u.sin();
    } else if (n + 1) % 4 == 0 {
        xc = c.radius * u.sin();
        yc = c.radius * -u.cos();
    } else if n % 4 == 0 {
        xc = c.radius * u.cos();
        yc = c.radius * u.sin();
    } else if (n - 1) % 4 == 0 {
        xc = c.radius * -u.sin();
        yc = c.radius * u.cos();
    }
    // Coord1 = Pos.XDirection(); Coord1.SetLinearForm(Xc, Coord1, Yc, YDir).
    xc * c.x_dir + yc * c.y_dir
}

/// OCCT `ElCLib::EllipseDN(U, Pos, MajorRadius, MinorRadius, N)`
/// (ElCLib.cxx L957-994).
fn ellipse_dn(e: &crate::geom::Ellipse3, u: f64, n: i32) -> DVec3 {
    let x_dir = e.major_dir;
    let y_dir = e.normal.cross(x_dir).normalize_or_zero();
    let mut xc = 0.0;
    let mut yc = 0.0;
    if n == 1 {
        xc = e.major_radius * -u.sin();
        yc = e.minor_radius * u.cos();
    } else if (n + 2) % 4 == 0 {
        xc = e.major_radius * -u.cos();
        yc = e.minor_radius * -u.sin();
    } else if (n + 1) % 4 == 0 {
        xc = e.major_radius * u.sin();
        yc = e.minor_radius * -u.cos();
    } else if n % 4 == 0 {
        xc = e.major_radius * u.cos();
        yc = e.minor_radius * u.sin();
    } else if (n - 1) % 4 == 0 {
        xc = e.major_radius * -u.sin();
        yc = e.minor_radius * u.cos();
    }
    xc * x_dir + yc * y_dir
}

/// OCCT `ElCLib::HyperbolaDN(U, Pos, MajorRadius, MinorRadius, N)`
/// (ElCLib.cxx L996-1018) — the `IsOdd(N)` / `IsEven(N)` branches.
fn hyperbola_dn(h: &crate::geom::Hyperbola3, u: f64, n: i32) -> DVec3 {
    let x_dir = h.major_dir;
    let y_dir = h.normal.cross(x_dir).normalize_or_zero();
    let mut xc = 0.0;
    let mut yc = 0.0;
    if n % 2 == 1 {
        // IsOdd(N)
        xc = h.semi_major * u.sinh();
        yc = h.semi_minor * u.cosh();
    } else if n % 2 == 0 {
        // IsEven(N)
        xc = h.semi_major * u.cosh();
        yc = h.semi_minor * u.sinh();
    }
    xc * x_dir + yc * y_dir
}

/// OCCT `ElCLib::ParabolaDN(U, Pos, Focal, N)` (ElCLib.cxx L1020-1047).
fn parabola_dn(p: &crate::geom::Parabola3, u: f64, n: i32) -> DVec3 {
    if n > 2 || n <= 0 {
        return DVec3::ZERO;
    }

    // Coord1 = Pos.XDirection() (the focal axis).
    let coord1 = p.axis_dir;
    if n == 1 {
        if p.focal_param.abs() <= GP_RESOLUTION {
            return coord1;
        }
        // SetLinearForm(U / (2 * Focal), Coord1, YDir) with
        // focal_param = 2 * Focal.
        return (u / p.focal_param) * coord1 + p.normal.cross(p.axis_dir).normalize_or_zero();
    }

    if p.focal_param.abs() <= GP_RESOLUTION {
        return DVec3::ZERO;
    }

    coord1 / p.focal_param
}

/// OCCT `Geom_BSplineCurve::EvalDN(U, N)` (Geom_BSplineCurve_1.cxx L300-316) —
/// `BSplCLib::DN` through the rcad [`BSplineCurve3::dn`] engine.
fn bspline_dn(bs: &crate::geom::BSplineCurve3, u: f64, n: i32) -> DVec3 {
    if n < 1 {
        panic!("Geom_UndefinedDerivative: Geom_BSplineCurve::EvalDN");
    }
    bs.dn(u, n as usize)
}

/// OCCT `Geom_BezierCurve::EvalDN(U, N)` (Geom_BezierCurve.cxx L601-620) —
/// `BSplCLib::DN(U, N, 0, aDeg, false, myPoles, Weights(), Knots(),
/// &Multiplicities(), V)`.
///
/// The Bezier's `Knots()` / `Multiplicities()` are `{0, 1}` / `{Deg+1, Deg+1}`
/// (Geom_BezierCurve.cxx L382-392), so the flat knot vector of the equivalent
/// BSpline payload is `[0] * (Deg+1) ++ [1] * (Deg+1)`.
fn bezier_dn(bz: &crate::geom::BezierCurve3, u: f64, n: i32) -> DVec3 {
    if n < 1 {
        panic!("Geom_UndefinedDerivative: Geom_BezierCurve::EvalDN");
    }
    let deg = bz.control_points.len().saturating_sub(1);
    let mut knots = vec![0.0f64; deg + 1];
    knots.extend(std::iter::repeat(1.0f64).take(deg + 1));
    let as_bspline = crate::geom::BSplineCurve3 {
        degree: deg,
        knots,
        control_points: bz.control_points.clone(),
        weights: bz.weights.clone(),
        is_periodic: false,
    };
    as_bspline.dn(u, n as usize)
}

/// OCCT `Geom_Curve::EvalDN(U, N)` over the rcad [`Curve3`] value — the
/// union of the per-type `EvalDN` bodies listed in the module header.
pub fn curve_dn(the_curve: &Curve3, u: f64, n: i32) -> DVec3 {
    match the_curve {
        Curve3::Line(l) => line_dn(l, n),
        Curve3::Circle(c) => circle_dn(c, u, n),
        Curve3::Ellipse(e) => ellipse_dn(e, u, n),
        Curve3::Hyperbola(h) => hyperbola_dn(h, u, n),
        Curve3::Parabola(p) => parabola_dn(p, u, n),
        Curve3::BSpline(b) => bspline_dn(b, u, n),
        Curve3::Bezier(b) => bezier_dn(b, u, n),
        // OCCT Geom_TrimmedCurve::EvalDN (cxx L240-243).
        Curve3::Trimmed(tc) => curve_dn(&tc.curve, u, n),
        _ => panic!(
            "GAP: Geom_Curve::EvalDN (TKG3d/Geom) is not translated for this curve kind \
             (the rcad union covers Geom_Line, Geom_Circle, Geom_Ellipse, Geom_Hyperbola, \
             Geom_Parabola, Geom_BSplineCurve, Geom_BezierCurve and Geom_TrimmedCurve) — \
             Geom_SurfaceOfLinearExtrusion::EvalDN / Geom_SurfaceOfRevolution::EvalDN"
        ),
    }
}
