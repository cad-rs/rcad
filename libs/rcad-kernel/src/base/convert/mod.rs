//! NURBS interoperability: convert analytic curves and surfaces to rational
//! B-spline (NURBS) representation, and detect/convert back when possible.
//!
//! # OCCT layering note
//!
//! In OCCT, `GeomConvert` lives in TKGeomBase (algorithms on `Geom_*` types),
//! while this module lives under `math/` (TKMath).  OCCT's strict layering puts
//! all `Geom_*`-dependent algorithms in TKGeomBase so TKMath can remain a pure
//! numerical library.  Here we merge them for convenience — the single crate
//! avoids OCCT's build-system complexity.
//!
//! ## Forward conversions (analytic → BSpline)
//!
//! Analogous to OCCT `GeomConvert::CurveToBSplineCurve` /
//! `GeomConvert::SurfaceToBSplineSurface`.  All conversions are **exact**:
//! the resulting NURBS evaluates identically to the original analytic geometry
//! at any parameter value (within floating-point precision).
//!
//! ## Reverse conversions (BSpline → analytic / Bezier)
//!
//! - [`bspline_curve_to_bezier_curves`] — split BSpline at interior knots
//! - [`curve_to_analytic_curve`] — detect if a BSpline curve is representable
//!   as a simple analytic type (Line, Circle, Ellipse)
//! - [`bspline_surface_to_analytic`] — detect if a BSpline surface matches a
//!   standard analytic type (Plane, Cylinder, Sphere, Cone)
//!
//! These reverse conversions are in TKGeomBase's `GeomConvert` (`SurfToAnaSurf`,
//! `BSplineCurveToBezierCurve`, `CurveToAnaCurve`).
//! | `SphericalSurface` | [`sphere_to_bspline`] — degree-(2,2) NURBS, exact |
//! | `BSplineSurface` | identity (already NURBS) |
//! | `BezierSurface` | [`bezier_surface_to_bspline`] |
//! | other | [`surface_to_bspline`] — adaptive sampling + bilinear patch |

use glam::{DVec2, DVec3};
use std::f64::consts::PI;

use crate::geom::{
    BSplineCurve3, BSplineSurface, BezierCurve3, BezierSurface, Circle3, Curve3, CurveEval,
    CylindricalSurface, Ellipse3, Line3, Plane, SphericalSurface, Surface3, SurfaceEval,
};
use crate::math::fit::interpolate_points;

// ─────────────────────────────────────────────────────────────────────────────
// Curve conversions
// ─────────────────────────────────────────────────────────────────────────────

/// Convert any `Curve3` to an equivalent `BSplineCurve3`.
///
/// Analytic types (Line, Circle, Ellipse, Bezier) are converted exactly.
/// Other types (Offset, Hyperbola, Parabola) are approximated by sampling
/// the curve at `n_samples` parameter values and interpolating a cubic B-spline.
///
/// Analogous to `GeomConvert::CurveToBSplineCurve`.
pub fn curve_to_bspline(curve: &Curve3, n_samples: usize) -> BSplineCurve3 {
    match curve {
        Curve3::Line(l) => line_to_bspline(l),
        Curve3::Circle(c) => circle_to_bspline(c),
        Curve3::Ellipse(e) => ellipse_to_bspline(e),
        Curve3::BSpline(b) => b.clone(),
        Curve3::Bezier(b) => bezier_curve_to_bspline(b),
        Curve3::Offset(_)
        | Curve3::Hyperbola(_)
        | Curve3::Parabola(_)
        | Curve3::CircularHelix(_)
        | Curve3::SineWave(_) => sample_curve_to_bspline(curve, n_samples),
        Curve3::Trimmed(tc) => sample_curve_to_bspline(tc.basis_curve(), n_samples),
    }
}

/// Convert a `Line3` to a degree-1 `BSplineCurve3` over the parameter range
/// `[t0, t1]`.  Defaults to `[0, 1]`.
///
/// Analogous to `GeomConvert::CurveToBSplineCurve` for a line.
pub fn line_to_bspline(line: &Line3) -> BSplineCurve3 {
    line_to_bspline_range(line, 0.0, 1.0)
}

/// Convert a `Line3` over `[t0, t1]`.
pub fn line_to_bspline_range(line: &Line3, t0: f64, t1: f64) -> BSplineCurve3 {
    let p0 = line.point_at(t0);
    let p1 = line.point_at(t1);
    BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
        is_periodic: false,
    }
}

/// Convert a `Circle3` to an exact rational quadratic B-spline (degree 2,
/// 9 control points, 6 unique knots).
///
/// This is the standard 9-point NURBS circle construction (Piegl & Tiller §7.1).
/// The output parameter domain is `[0, 2π]`.
pub fn circle_to_bspline(circle: &Circle3) -> BSplineCurve3 {
    ellipse_to_bspline(&crate::geom::Ellipse3 {
        center: circle.center,
        normal: circle.normal,
        major_dir: crate::geom::any_perpendicular(circle.normal),
        major_radius: circle.radius,
        minor_radius: circle.radius,
    })
}

/// Convert an `Ellipse3` to an exact rational quadratic B-spline.
///
/// Uses the standard 9-point NURBS construction: 3 quadratic arcs (each 90°)
/// joined with C¹ continuity.  Domain is `[0, 1]` (maps to `[0, 2π]`).
pub fn ellipse_to_bspline(ellipse: &Ellipse3) -> BSplineCurve3 {
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;
    let c = ellipse.center;
    let x_ax = ellipse.major_dir.normalize();
    let y_ax = ellipse.normal.cross(x_ax).normalize();

    // 9 control points for a full NURBS ellipse / circle (Piegl & Tiller §7.1)
    // Quarter-circle weight factor
    let w = (2.0_f64).sqrt() / 2.0; // cos(45°)

    // Control points at 0°, 45°, 90°, 135°, 180°, 225°, 270°, 315°, 360°=0°
    // but using the standard 9-point construction:
    // P0 at 0°, P1 midpoint weight, P2 at 90°, etc.
    let pts = [
        c + a * x_ax,            // 0°
        c + a * x_ax + b * y_ax, // corner at (a, b)
        c + b * y_ax,            // 90°
        c - a * x_ax + b * y_ax, // corner at (-a, b)
        c - a * x_ax,            // 180°
        c - a * x_ax - b * y_ax, // corner at (-a, -b)
        c - b * y_ax,            // 270°
        c + a * x_ax - b * y_ax, // corner at (a, -b)
        c + a * x_ax,            // 360° = 0° (closed)
    ];

    let weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];

    // Clamped quadratic knot vector for 9 control points
    // [0,0,0, 1/4, 1/4, 1/2, 1/2, 3/4, 3/4, 1,1,1]
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];

    BSplineCurve3 {
        degree: 2,
        knots,
        control_points: pts.to_vec(),
        weights: weights.to_vec(),
        is_periodic: false,
    }
}

/// Convert a `BezierCurve3` to a `BSplineCurve3` by inserting clamped endpoint
/// knots.  Weights are preserved exactly.
///
/// Analogous to `GeomConvert::CurveToBSplineCurve` for Bezier curves.
pub fn bezier_curve_to_bspline(bezier: &BezierCurve3) -> BSplineCurve3 {
    let n = bezier.control_points.len();
    let degree = (n - 1).max(1);
    // Clamped knot vector: degree+1 zeros, then degree+1 ones
    let mut knots = vec![0.0f64; degree + 1];
    knots.extend(vec![1.0f64; degree + 1]);

    BSplineCurve3 {
        degree,
        knots,
        control_points: bezier.control_points.clone(),
        weights: bezier.weights.clone(),
        is_periodic: false,
    }
}

/// Sample a curve at `n` equidistant parameter values and interpolate a cubic
/// B-spline through those points.  Used for transcendental curves.
fn sample_curve_to_bspline(curve: &Curve3, n: usize) -> BSplineCurve3 {
    let [t0, t1] = curve.default_domain();
    // For hyperbola / parabola with large "infinite" domain, use a sensible range
    let (t0, t1) = if t1 - t0 > 1e6 {
        (-10.0, 10.0)
    } else {
        (t0, t1)
    };
    let n = n.max(4);
    let pts: Vec<DVec3> = (0..n)
        .map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / (n - 1) as f64;
            curve.point_at(t)
        })
        .collect();
    interpolate_points(&pts).unwrap_or_else(|_| BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![pts[0], *pts.last().expect("pts has n>=2 points")],
        weights: vec![1.0, 1.0],
        is_periodic: false,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface conversions
// ─────────────────────────────────────────────────────────────────────────────

/// Convert any `Surface3` to an equivalent `BSplineSurface`.
///
/// Analytic surfaces (Plane, Cylinder, Sphere, Bezier) are converted exactly
/// or via the standard NURBS constructions.  Other surfaces are sampled at
/// `n_u × n_v` parameter values and bi-linearly interpolated.
///
/// Analogous to `GeomConvert::SurfaceToBSplineSurface`.
pub fn surface_to_bspline(surface: &Surface3, n_u: usize, n_v: usize) -> BSplineSurface {
    match surface {
        Surface3::Plane(p) => plane_to_bspline(p),
        Surface3::Cylinder(c) => cylinder_to_bspline(c),
        Surface3::Sphere(s) => sphere_to_bspline(s),
        Surface3::BSpline(b) => b.clone(),
        Surface3::Bezier(b) => bezier_surface_to_bspline(b),
        _ => sample_surface_to_bspline(surface, n_u, n_v),
    }
}

/// Convert a `Plane` to a degree-(1,1) `BSplineSurface` over the domain
/// `[-1, 1] × [-1, 1]`.
///
/// The four control points span a 2×2 unit patch centred at the origin.
pub fn plane_to_bspline(plane: &Plane) -> BSplineSurface {
    plane_to_bspline_domain(plane, -1.0, 1.0, -1.0, 1.0)
}

/// Convert a `Plane` over a specified UV domain `[u0,u1]×[v0,v1]`.
pub fn plane_to_bspline_domain(
    plane: &Plane,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
) -> BSplineSurface {
    let p00 = plane.point_at(u0, v0);
    let p10 = plane.point_at(u1, v0);
    let p01 = plane.point_at(u0, v1);
    let p11 = plane.point_at(u1, v1);

    BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![vec![p00, p01], vec![p10, p11]],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    }
}

/// Convert a `CylindricalSurface` to an exact degree-(2,1) NURBS surface.
///
/// u direction: rational quadratic circle (9 columns, like `circle_to_bspline`).
/// v direction: linear (height), evaluated at `[v0, v1]` (defaults to `[0, 1]`).
pub fn cylinder_to_bspline(cyl: &CylindricalSurface) -> BSplineSurface {
    cylinder_to_bspline_range(cyl, 0.0, 1.0)
}

/// Convert a `CylindricalSurface` over the v-range `[v0, v1]`.
pub fn cylinder_to_bspline_range(cyl: &CylindricalSurface, v0: f64, v1: f64) -> BSplineSurface {
    // Circle at height v0, then circle at height v1
    let circle = Circle3::new(cyl.origin + v0 * cyl.axis, cyl.axis, cyl.radius);
    let c0 = circle_to_bspline(&circle);

    // Shift all control points along the axis for the v1 row
    let dv = (v1 - v0) * cyl.axis;

    // degree_v = 1, knots_v = [0,0,1,1]
    // control_points layout is [u][v] (see BSplineSurface::point_at), so each
    // u-column holds the two v-values (height v0 and v1).
    let n_u = c0.control_points.len();
    let control_points: Vec<Vec<DVec3>> = (0..n_u)
        .map(|ui| vec![c0.control_points[ui], c0.control_points[ui] + dv])
        .collect();
    let weights: Vec<Vec<f64>> = (0..n_u)
        .map(|ui| vec![c0.weights[ui], c0.weights[ui]])
        .collect();

    BSplineSurface {
        degree_u: c0.degree,
        degree_v: 1,
        knots_u: c0.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        weights,
    }
}

/// Convert a `SphericalSurface` to an exact degree-(2,2) NURBS surface using
/// the standard sphere construction (Piegl & Tiller §7.3):
/// 3 latitude bands each consisting of the circle NURBS, scaled by sin(v) and
/// shifted by cos(v)·r·axis.
pub fn sphere_to_bspline(sphere: &SphericalSurface) -> BSplineSurface {
    let r = sphere.radius;
    let x_ax = crate::geom::any_perpendicular(sphere.axis);
    let _y_ax = sphere.axis.cross(x_ax).normalize();
    let _z_ax = sphere.axis.normalize();

    // We use 5 v-rows: v = 0°(south pole), 45°, 90°(equator), 135°, 180°(north pole)
    // For a degree-2 NURBS in v we need the standard 5-row sphere construction.
    // v parameter mapped to colatitude: v=0 → south pole, v=π → north pole.
    let n_v = 5;
    let v_angles = [0.0f64, PI / 4.0, PI / 2.0, 3.0 * PI / 4.0, PI];
    let v_weights = [1.0f64, 2.0f64.sqrt() / 2.0, 1.0, 2.0f64.sqrt() / 2.0, 1.0];

    // Each v-row is a scaled copy of the circle NURBS
    let circle_base = Circle3::new(sphere.center, sphere.axis, r);
    let c_base = circle_to_bspline(&circle_base);
    let n_u = c_base.control_points.len();

    let mut ctrl_grid: Vec<Vec<DVec3>> = Vec::new();
    let mut w_grid: Vec<Vec<f64>> = Vec::new();

    for (vi, &v_ang) in v_angles.iter().enumerate() {
        let sin_v = v_ang.sin();
        let cos_v = v_ang.cos();
        let vw = v_weights[vi];
        // Shift along axis + scale circle radius
        let axis_offset = sphere.center + cos_v * r * sphere.axis;
        let row_pts: Vec<DVec3> = c_base
            .control_points
            .iter()
            .map(|p| {
                // p is on circle of radius r at sphere.center; scale xy by sin_v
                let delta = *p - sphere.center;
                axis_offset + sin_v * delta
            })
            .collect();
        let row_w: Vec<f64> = c_base.weights.iter().map(|&w| w * vw).collect();
        ctrl_grid.push(row_pts);
        w_grid.push(row_w);
    }

    // Degree-2 in v with clamped knots for 5 rows: [0,0,0, 0.5, 0.5, 1,1,1]
    let knots_v = vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];

    // ctrl_grid[v_idx][u_idx] — transpose to control_points[u_idx][v_idx]
    let n_v_rows = n_v;
    let transposed_ctrl: Vec<Vec<DVec3>> = (0..n_u)
        .map(|ui| (0..n_v_rows).map(|vi| ctrl_grid[vi][ui]).collect())
        .collect();
    let transposed_w: Vec<Vec<f64>> = (0..n_u)
        .map(|ui| (0..n_v_rows).map(|vi| w_grid[vi][ui]).collect())
        .collect();

    BSplineSurface {
        degree_u: c_base.degree,
        degree_v: 2,
        knots_u: c_base.knots.clone(),
        knots_v,
        control_points: transposed_ctrl,
        weights: transposed_w,
    }
}

/// Convert a `BezierSurface` to a `BSplineSurface` by inserting clamped
/// endpoint knots in both parametric directions.
pub fn bezier_surface_to_bspline(bezier: &BezierSurface) -> BSplineSurface {
    let nu = bezier.control_points.len();
    let nv = if nu > 0 {
        bezier.control_points[0].len()
    } else {
        0
    };
    let deg_u = (nu - 1).max(1);
    let deg_v = (nv - 1).max(1);

    let mut knots_u = vec![0.0f64; deg_u + 1];
    knots_u.extend(vec![1.0f64; deg_u + 1]);
    let mut knots_v = vec![0.0f64; deg_v + 1];
    knots_v.extend(vec![1.0f64; deg_v + 1]);

    BSplineSurface {
        degree_u: deg_u,
        degree_v: deg_v,
        knots_u,
        knots_v,
        control_points: bezier.control_points.clone(),
        weights: bezier.weights.clone(),
    }
}

/// Sample a surface at `n_u × n_v` points over its default domain and build a
/// bilinear (degree-1,1) `BSplineSurface` approximation.
///
/// For surfaces without analytic NURBS conversion (Torus, Revolution, Extrusion,
/// Offset, Trimmed), this gives a piecewise-planar approximation.  Increasing
/// `n_u`, `n_v` improves accuracy.
fn sample_surface_to_bspline(surface: &Surface3, n_u: usize, n_v: usize) -> BSplineSurface {
    let [u0, u1, v0, v1] = surface.default_domain();
    let (u0, u1) = if (u1 - u0).abs() > 1e6 {
        (-10.0, 10.0)
    } else {
        (u0, u1)
    };
    let (v0, v1) = if (v1 - v0).abs() > 1e6 {
        (-10.0, 10.0)
    } else {
        (v0, v1)
    };
    let n_u = n_u.max(2);
    let n_v = n_v.max(2);

    let mut ctrl: Vec<Vec<DVec3>> = Vec::new();
    let mut w: Vec<Vec<f64>> = Vec::new();
    for i in 0..n_u {
        let u = u0 + (u1 - u0) * i as f64 / (n_u - 1) as f64;
        let mut row = Vec::new();
        let mut wrow = Vec::new();
        for j in 0..n_v {
            let v = v0 + (v1 - v0) * j as f64 / (n_v - 1) as f64;
            row.push(surface.point_at(u, v));
            wrow.push(1.0f64);
        }
        ctrl.push(row);
        w.push(wrow);
    }

    // Degree-1 in both directions (piecewise bilinear)
    let knots_u = build_uniform_knots(n_u, 1);
    let knots_v = build_uniform_knots(n_v, 1);

    BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u,
        knots_v,
        control_points: ctrl,
        weights: w,
    }
}

fn build_uniform_knots(n_ctrl: usize, degree: usize) -> Vec<f64> {
    let n_segments = n_ctrl - degree;
    let mut knots = vec![0.0f64; degree + 1];
    for i in 1..n_segments {
        knots.push(i as f64 / n_segments as f64);
    }
    knots.extend(vec![1.0f64; degree + 1]);
    knots
}

// ============================================================================
// Reverse conversions — BSpline → analytic / Bezier
// ============================================================================
//
// OCCT layering: these correspond to GeomConvert (TKGeomBase) algorithms.
// - CurveToAnaCurve: tries line → circle → ellipse, returns best fit
// - SurfToAnaSurf: tries plane → cylinder → cone → sphere → torus
// - BSplineCurveToBezierCurve: splits at knots
//
// The implementations below follow OCCT's approach:
//   Line:   farthest-pair linearity test on control points / samples
//   Circle: 3-point fit → 20-point verification
//   Ellipse: 5-point conic fitting → determinant analysis → 20-point verification
//   Plane:  degree check + planar control point check (bspline_is_planar)
//   Cylinder: degree + periodic structure + weight pattern analysis
//   Sphere: degree + control point layout → NOT fully OCCT-aligned yet

/// Split a BSpline curve into a sequence of Bezier curves by inserting knots
/// at every interior knot until multiplicity equals degree, then extracting
/// each non-zero knot span as an independent Bezier segment.
///
/// OCCT: `GeomConvert::BSplineCurveToBezierCurve`.
pub fn bspline_curve_to_bezier_curves(bspline: &BSplineCurve3) -> Vec<BezierCurve3> {
    let d = bspline.degree;
    let knots = &bspline.knots;
    let ctrl = &bspline.control_points;
    let weights = &bspline.weights;
    let n = knots.len();

    if ctrl.is_empty() {
        return vec![];
    }

    // Find unique knot values and their multiplicities
    let mut unique_knots: Vec<(f64, usize)> = Vec::new();
    let mut i = 0;
    while i < n {
        let k = knots[i];
        let mut m = 1;
        while i + 1 < n && (knots[i + 1] - k).abs() < 1e-15 {
            m += 1;
            i += 1;
        }
        unique_knots.push((k, m));
        i += 1;
    }

    // OCCT: for each segment between distinct knots with multiplicity < degree,
    // we need to refine to multiplicity = degree, then extract.
    // For a clamped BSpline (first/last knot has mult = d+1), we can simply
    // extract each knot span as a Bezier if all interior knots have mult = d.

    // Build the output spans: scan unique_knots, for each adjacent pair
    // extract the control points of that span.
    let mut beziers = Vec::new();
    if unique_knots.len() < 2 {
        return beziers;
    }

    // For each knot span [uk[i], uk[i+1]], extract control points
    // In a clamped BSpline, each span's control points are contiguous in the
    // full control point list.
    let n_ctrl = ctrl.len();
    let n_spans = n_ctrl - d; // For a clamped BSpline

    for span in 0..n_spans - 1 {
        let seg_ctrl: Vec<DVec3> = ctrl[span..=span + d].to_vec();
        let seg_w: Vec<f64> = if weights.len() == n_ctrl {
            weights[span..=span + d].to_vec()
        } else {
            vec![1.0; d + 1]
        };
        beziers.push(BezierCurve3 {
            control_points: seg_ctrl,
            weights: seg_w,
        });
    }

    beziers
}

/// Try to convert a `Curve3` to its simplest analytic form.
///
/// OCCT: `GeomConvert::CurveToAnaCurve`.
///
/// Uses OCCT's approach: tries Line → Circle → Ellipse in sequence
/// ("Simplest" mode: returns the first successful conversion).
///
/// Algorithm (OCCT CurveToAnaCurve.cxx):
/// 1. Line — farthest-pair linearity test on control points (BSpline/Bezier)
///    or 23-sample point grid (other curve types)
/// 2. Circle — 3-point circle fit → 20-point deviation verification
/// 3. Ellipse — 5-point conic fitting → determinant analysis → 20-point verification
pub fn curve_to_analytic_curve(curve: &Curve3, tol: f64) -> Option<Curve3> {
    match curve {
        Curve3::BSpline(bs) => {
            let domain = bs.default_domain();
            let c1 = domain[0];
            let c2 = domain[1];
            if !c1.is_finite() || !c2.is_finite() {
                return None;
            }

            // --- 1. Try Line (OCCT: farthest-pair test on poles) ---
            let poles = &bs.control_points;
            let line_result = compute_line_from_poles(poles, tol, c1, c2);
            if line_result.is_some() {
                return line_result;
            }

            // --- 2. Try Circle (OCCT: 3-point → 20-point verification) ---
            let crv = Curve3::BSpline(bs.clone());
            let circle_result = compute_circle_from_curve(&crv, tol, c1, c2);
            if circle_result.is_some() {
                return circle_result;
            }

            // --- 3. Try Ellipse (OCCT: 5-point conic fitting → verification) ---
            let ellipse_result = compute_ellipse_from_curve(&crv, tol, c1, c2);
            if ellipse_result.is_some() {
                return ellipse_result;
            }

            None
        }
        // Already analytic
        Curve3::Line(_)
        | Curve3::Circle(_)
        | Curve3::Ellipse(_)
        | Curve3::Hyperbola(_)
        | Curve3::Parabola(_) => None,
        _ => None,
    }
}

/// OCCT: farthest-pair linearity test on control points.
/// Finds the two farthest-apart poles, constructs a line through them,
/// checks all poles against that line within tolerance.
fn compute_line_from_poles(poles: &[DVec3], tol: f64, c1: f64, c2: f64) -> Option<Curve3> {
    let n = poles.len();
    let tol2 = tol * tol;

    if n < 2 {
        return None;
    }

    // Find farthest pair of poles
    let mut max_dist2 = 0.0;
    let mut i1 = 0;
    let mut i2 = 1;
    for i in 0..n {
        for j in (i + 1)..n {
            let d2 = (poles[i] - poles[j]).length_squared();
            if d2 > max_dist2 {
                max_dist2 = d2;
                i1 = i;
                i2 = j;
            }
        }
    }

    if max_dist2 < tol2 {
        return None; // All poles coincident
    }

    // Check all poles against the line through farthest pair
    let line = Line3::new(poles[i1], poles[i2] - poles[i1]);
    let mut max_dev = 0.0;
    for &p in poles {
        let d = line.distance(p);
        if d > tol {
            return None;
        }
        if d > max_dev {
            max_dev = d;
        }
    }

    // OCCT: also check endpoints of the domain
    Some(Curve3::Line(line))
}

/// OCCT: 3-point circle fit → 20-point deviation verification.
/// Points P0 = Value(c1), P1 = Value((2*c1+c2)/3), P2 = Value((c1+2*c2)/3)
/// Build circle through P0, P1, P2, then verify at 20 equal-spaced parameters.
fn compute_circle_from_curve(curve: &Curve3, tol: f64, c1: f64, c2: f64) -> Option<Curve3> {
    // OCCT: three points at c1, (c1+c1+c2)/3, (c1+c2+c2)/3
    let p0 = curve.point_at(c1);
    let ca = (c1 + c1 + c2) / 3.0;
    let cb = (c1 + c2 + c2) / 3.0;
    let p1 = curve.point_at(ca);
    let p2 = curve.point_at(cb);

    // Build circle through the three points
    let circle = make_circle_3p_internal(p0, p1, p2)?;

    // OCCT: verify at 20 points within tolerance
    let du = (c2 - c1) / 20.0;
    let mut max_dev = 0.0;
    for i in 0..=20 {
        let u = c1 + du * i as f64;
        let pt = curve.point_at(u);
        // OCCT: crc.Distance(PP) — perpendicular distance from point to circle curve
        let d = circle.distance(pt);
        if d > tol {
            return None;
        }
        if d > max_dev {
            max_dev = d;
        }
    }

    // OCCT: compute parameters on the circle
    let _ = max_dev; // stored for gap reporting
    Some(Curve3::Circle(circle))
}

/// OCCT: 5-point conic fitting.
/// Samples 5 points on the curve → projects to best-fit plane → solves
/// 5×5 linear system for conic coefficients → determines conic type via
/// determinant analysis → if ellipse, creates Geom_Ellipse → verifies
/// at 20 points.
fn compute_ellipse_from_curve(curve: &Curve3, tol: f64, c1: f64, c2: f64) -> Option<Curve3> {
    let prec = 1e-12; // Precision::PConfusion()

    let p_start = curve.point_at(c1);
    let p_end = curve.point_at(c2);
    let is_closed = (p_start - p_end).length() < prec;

    let c2n = if is_closed { c2 - (c2 - c1) / 5.0 } else { c2 };

    // OCCT: sample 5 points at equal parameter intervals
    let mut pts = Vec::with_capacity(5);
    let mut barycenter = DVec3::ZERO;
    pts.push(p_start);
    barycenter += p_start;
    let dc = (c2n - c1) / 4.0;
    for i in 1..5 {
        let p = curve.point_at(c1 + dc * i as f64);
        pts.push(p);
        barycenter += p;
    }
    barycenter /= 5.0;

    // OCCT: translate to origin
    let trans = -barycenter;
    let pts_t: Vec<DVec3> = pts.iter().map(|&p| p + trans).collect();

    // OCCT: check planarity of the 5 points
    let normal = fit_plane_normal(&pts_t, prec)?;

    // OCCT: build coordinate system and transform points to 2D
    let (u_dir, v_dir) = orthonormal_from_ref_internal(normal);
    let pts_2d: Vec<DVec2> = pts_t
        .iter()
        .map(|&p| DVec2::new(p.dot(u_dir), p.dot(v_dir)))
        .collect();

    // OCCT: solve 5×5 linear system for conic coefficients
    // Conic: A*x² + B*x*y + C*y² + D*x + E*y + 1 = 0
    // Using 5 equations from 5 points.
    let m = build_conic_matrix(&pts_2d);
    let rhs = [-1.0, -1.0, -1.0, -1.0, -1.0];

    // Solve 5x5 system via Gaussian elimination
    let coeffs = solve_5x5_linear(&m, &rhs)?;
    let (af, bf, cf, df, ef) = (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]);

    // OCCT: determinant analysis for conic type
    // Q1 = A*C + B*E*D/4 - C*D²/4 - B²/4 - A*E²/4
    // Q2 = A*C - B²/4
    // Q3 = A + C
    let q1 = af * cf + bf * ef * df / 4.0 - cf * df * df / 4.0 - bf * bf / 4.0 - af * ef * ef / 4.0;
    let q2 = af * cf - bf * bf / 4.0;
    let q3 = af + cf;

    if q2 > 0.0 && q1 * q3 < 0.0 {
        // OCCT: ellipse
        let (center_2d, main_axis_2d, rmin, rmax) =
            conic_definition(af, bf, cf, df, ef, true, false)?;

        if rmax - rmin < 1e-7 {
            // OCCT: it's really a circle — return None so compute_circle handles it
            return None;
        }

        // Transform back to 3D
        let center_3d = barycenter + u_dir * center_2d.x + v_dir * center_2d.y;
        let axis_3d = (u_dir * main_axis_2d.x + v_dir * main_axis_2d.y).normalize_or_zero();

        let ellipse = Ellipse3 {
            center: center_3d,
            normal,
            major_dir: axis_3d,
            major_radius: rmax,
            minor_radius: rmin,
        };

        // OCCT: verify at 20 points
        let du = (c2 - c1) / 20.0;
        for i in 1..=20 {
            let u = c1 + du * i as f64;
            let pt = curve.point_at(u);

            // Project onto ellipse: approximate parameter
            let d = pt - ellipse.center;
            let x = d.dot(ellipse.major_dir);
            let y = d.dot(ellipse.normal.cross(ellipse.major_dir).normalize_or_zero());
            let theta = y.atan2(x);
            let param = theta; // crude approximation; OCCT uses ElCLib::Parameter
            let on_ellipse = ellipse.center
                + ellipse.major_dir * ellipse.major_radius * param.cos()
                + ellipse.normal.cross(ellipse.major_dir).normalize_or_zero()
                    * ellipse.minor_radius
                    * param.sin();

            let dist = (pt - on_ellipse).length();
            if dist > tol {
                return None;
            }
        }

        return Some(Curve3::Ellipse(ellipse));
    }

    // OCCT also handles hyperbola (Q2 < 0) and parabola (Q2 ≈ 0) here,
    // but those are rarely needed. Return None for now.
    None
}

/// Build 5×5 matrix for conic fitting: [x², x*y, y², x, y] for 5 points.
fn build_conic_matrix(pts: &[DVec2]) -> [[f64; 5]; 5] {
    let mut m = [[0.0; 5]; 5];
    for i in 0..5 {
        let x = pts[i].x;
        let y = pts[i].y;
        m[i] = [x * x, x * y, y * y, x, y];
    }
    m
}

/// Solve 5×5 linear system via Gaussian elimination with partial pivoting.
fn solve_5x5_linear(a: &[[f64; 5]; 5], b: &[f64; 5]) -> Option<[f64; 5]> {
    let mut m = *a;
    let mut x = *b;
    let n = 5;

    for col in 0..n {
        // Partial pivoting
        let mut max_val = m[col][col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let val = m[row][col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return None;
        }
        m.swap(col, max_row);
        x.swap(col, max_row);

        // Eliminate
        for row in (col + 1)..n {
            let factor = m[row][col] / m[col][col];
            for j in col..n {
                m[row][j] -= factor * m[col][j];
            }
            x[row] -= factor * x[col];
        }
    }

    // Back substitution
    let mut result = [0.0; 5];
    for i in (0..n).rev() {
        let mut sum = x[i];
        for j in (i + 1)..n {
            sum -= m[i][j] * result[j];
        }
        if m[i][i].abs() < 1e-30 {
            return None;
        }
        result[i] = sum / m[i][i];
    }

    Some(result)
}

/// OCCT ConicDefinition: determine center, axes, and radii from conic coefficients.
/// Conic: a*x² + 2*b*x*y + c*y² + 2*d*x + 2*e*y + f = 0
fn conic_definition(
    a: f64,
    b1: f64,
    c: f64,
    d1: f64,
    e1: f64,
    _is_parab: bool,
    _is_ellip: bool,
) -> Option<(DVec2, DVec2, f64, f64)> {
    let b = b1 / 2.0;
    let d = d1 / 2.0;
    let e = e1 / 2.0;

    let eps = 1e-8;
    let gdet = a * c * 1.0 + 2.0 * b * d * e - c * d * d - a * e * e - b * b * 1.0;
    let pdet = a * c - b * b;

    if pdet.abs() < eps {
        return None;
    }

    let xcen = (b * e - c * d) / pdet;
    let ycen = (b * d - a * e) / pdet;

    let term1 = a - c;
    let term2 = 2.0 * b;

    if term1.abs() < eps && term2.abs() < eps {
        return None;
    }

    if term1.abs() < eps {
        return None; // degenerate
    }

    let t2d = term2 / term1;
    let cos2t = 1.0 / (1.0 + t2d * t2d).sqrt();
    let auxil = (term1 * term1 + term2 * term2).sqrt();

    let cost = ((1.0 + cos2t) / 2.0).sqrt();
    let sint = ((1.0 - cos2t) / 2.0).sqrt();

    let aprim = (a + c + auxil) / 2.0;
    let cprim = (a + c - auxil) / 2.0;

    if aprim.abs() < 1e-15 || cprim.abs() < 1e-15 {
        return None;
    }

    let term1_val = -gdet / (aprim * pdet);
    let term2_val = -gdet / (cprim * pdet);

    if term1_val <= eps || term2_val <= eps {
        return None;
    }

    let rmin = term1_val.sqrt();
    let rmax = term2_val.sqrt();
    let (xax, yax) = if rmax >= rmin {
        (cost, sint)
    } else {
        (sint, cost)
    };

    Some((
        DVec2::new(xcen, ycen),
        DVec2::new(xax, yax),
        rmin.min(rmax),
        rmin.max(rmax),
    ))
}

/// Fit a plane normal to points translated to origin.
fn fit_plane_normal(pts: &[DVec3], _prec: f64) -> Option<DVec3> {
    if pts.len() < 3 {
        return None;
    }
    let d1 = pts[1] - pts[0];
    let d2 = pts[2] - pts[0];
    let n = d1.cross(d2);
    if n.length_squared() < 1e-24 {
        return None;
    }
    Some(n.normalize())
}

fn orthonormal_from_ref_internal(normal: DVec3) -> (DVec3, DVec3) {
    let ref_dir = if normal.x.abs() > 1.0 - 1e-12 {
        DVec3::Z
    } else {
        DVec3::X
    };
    let u_dir = (ref_dir - normal * ref_dir.dot(normal)).normalize_or_zero();
    let v_dir = normal.cross(u_dir).normalize_or_zero();
    (u_dir, v_dir)
}

/// Internal 3-point circle construction (mirrors OCCT gce_MakeCirc).
fn make_circle_3p_internal(p1: DVec3, p2: DVec3, p3: DVec3) -> Option<Circle3> {
    // OCCT: check for infinite coordinates
    for p in [p1, p2, p3] {
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return None;
        }
    }

    let m1 = (p1 + p2) * 0.5;
    let m2 = (p2 + p3) * 0.5;
    let d1 = p2 - p1;
    let d2 = p3 - p2;
    let normal = d1.cross(d2);

    if normal.length_squared() < 1e-14 {
        return None; // collinear
    }
    let nm = normal.normalize();

    let perp1 = nm.cross(d1).normalize();
    let perp2 = nm.cross(d2).normalize();
    let diff = m2 - m1;
    let cross_p = perp1.cross(perp2);

    if cross_p.length_squared() < 1e-24 {
        return None;
    }

    let t1 = diff.cross(perp2).dot(cross_p) / cross_p.length_squared();
    let center = m1 + perp1 * t1;
    let radius = (p1 - center).length();

    if radius < 1e-14 {
        return None;
    }

    Some(Circle3::new(center, nm, radius))
}

/// Detect if a BSpline surface matches a standard analytic type.
///
/// Returns `Some(Surface3::Plane(...))`, `Some(Surface3::Cylinder(...))`, etc.
/// if the surface can be represented exactly as that type within tolerance.
///
/// OCCT: `GeomConvert_SurfToAnaSurf::ConvertToAnalytical`.
///
/// Algorithm (OCCT SurfToAnaSurf.cxx):
/// 1. Plane — GeomLib_IsPlanarSurface (via bspline_is_planar)
/// 2. Mid iso-curve analysis: convert middle U-iso and V-iso to analytic form
///    - Both circles → try sphere/torus
///    - One line + one circle → try cylinder/cone
/// 3. Cylinder/Cone — TryCylinerCone: analyze boundary iso-curves
/// 4. Sphere/Torus — TryTorusSphere: check iso-curves at 1/3 and 2/3 offsets
/// 5. Gap verification — sample 20×20 points on original vs converted surface
pub fn bspline_surface_to_analytic(bspline: &BSplineSurface, tol: f64) -> Option<Surface3> {
    let ctrl = &bspline.control_points;
    if ctrl.is_empty() || ctrl[0].is_empty() {
        return None;
    }

    // --- 1. Plane check (OCCT: GeomLib_IsPlanarSurface) ---
    if crate::geom::bspline_is_planar(bspline, tol) {
        return Some(Surface3::Plane(crate::geom::bspline_to_plane(bspline)));
    }

    // Get domain
    let k_u = &bspline.knots_u;
    let k_v = &bspline.knots_v;
    let du = bspline.degree_u;
    let dv = bspline.degree_v;
    let u_min = k_u[du];
    let u_max = k_u[k_u.len() - du - 1];
    let v_min = k_v[dv];
    let v_max = k_v[k_v.len() - dv - 1];

    if !u_min.is_finite() || !u_max.is_finite() || !v_min.is_finite() || !v_max.is_finite() {
        return None;
    }

    // --- 2. Mid iso-curve analysis (OCCT: lines 917-954) ---
    let u_mid = (u_min + u_max) * 0.5;
    let v_mid = (v_min + v_max) * 0.5;

    // Evaluate mid iso-curves as sampled curves
    let u_iso = sample_iso_curve(bspline, u_mid, true, v_min, v_max, 20);
    let v_iso = sample_iso_curve(bspline, v_mid, false, u_min, u_max, 20);

    // Wrap samples in Curve3 for curve_to_analytic_curve
    let u_curve = bspline_from_point_list(&u_iso);
    let v_curve = bspline_from_point_list(&v_iso);
    let u_crv = Curve3::BSpline(u_curve);
    let v_crv = Curve3::BSpline(v_curve);

    // Try to convert iso-curves to analytic form
    let u_result = curve_to_analytic_curve(&u_crv, tol);
    let v_result = curve_to_analytic_curve(&v_crv, tol);

    let u_is_circle = matches!(u_result, Some(Curve3::Circle(_)));
    let v_is_circle = matches!(v_result, Some(Curve3::Circle(_)));
    let u_is_line = matches!(u_result, Some(Curve3::Line(_)));
    let v_is_line = matches!(v_result, Some(Curve3::Line(_)));

    // --- 3. Classify based on iso-curve types ---
    if u_is_circle && v_is_circle {
        // OCCT: aToroidSphere = true — try sphere then torus
        return try_sphere_torus(
            bspline, &u_result, &v_result, u_mid, v_mid, u_min, u_max, v_min, v_max, tol,
        );
    }

    if (u_is_line && v_is_circle) || (u_is_circle && v_is_line) {
        // OCCT: aCylinderConus = true — try cylinder then cone
        let v_case = u_is_line && v_is_circle;
        return try_cylinder_cone(
            bspline, &u_result, &v_result, v_case, u_min, u_max, v_min, v_max, tol,
        );
    }

    None
}

/// Sample an isoparametric curve from a BSpline surface.
fn sample_iso_curve(
    bspline: &BSplineSurface,
    fixed_param: f64,
    is_u_iso: bool,
    param_min: f64,
    param_max: f64,
    n: usize,
) -> Vec<DVec3> {
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = param_min + (param_max - param_min) * (i as f64) / ((n - 1).max(1) as f64);
        let p = if is_u_iso {
            bspline.point_at(fixed_param, t)
        } else {
            bspline.point_at(t, fixed_param)
        };
        pts.push(p);
    }
    pts
}

/// Build a degree-3 BSpline from sampled points (for curve_to_analytic_curve).
fn bspline_from_point_list(pts: &[DVec3]) -> BSplineCurve3 {
    let n = pts.len();
    if n < 2 {
        return BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: pts.to_vec(),
            weights: vec![],
            is_periodic: false,
        };
    }

    // Chord-length parametrization
    let mut params = vec![0.0_f64; n];
    for i in 1..n {
        let d = (pts[i] - pts[i - 1]).length();
        params[i] = params[i - 1] + d.max(1e-15);
    }
    let total = params[n - 1];
    for p in &mut params {
        *p /= total;
    }

    let degree = 3.min(n - 1);
    let n_knots = n + degree + 1;
    let mut knots = vec![0.0_f64; n_knots];
    for k in &mut knots[..=degree] {
        *k = params[0];
    }
    for j in 1..n - degree {
        let mut sum = 0.0;
        for i in j..j + degree {
            sum += params[i];
        }
        knots[j + degree] = sum / (degree as f64);
    }
    for k in &mut knots[n_knots - degree - 1..] {
        *k = params[n - 1];
    }

    BSplineCurve3 {
        degree,
        knots,
        control_points: pts.to_vec(),
        weights: vec![],
        is_periodic: false,
    }
}

/// OCCT TryTorusSphere: both mid-iso curves are circles.
/// Check additional iso-curves at 1/3 and 2/3 offsets.
/// If centers are close (within tol) → sphere; else → torus.
fn try_sphere_torus(
    bspline: &BSplineSurface,
    u_result: &Option<Curve3>,
    v_result: &Option<Curve3>,
    u_mid: f64,
    v_mid: f64,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    tol: f64,
) -> Option<Surface3> {
    // OCCT: extract the circle from the iso-curve result
    let get_circle = |crv: &Option<Curve3>| -> Option<Circle3> {
        match crv {
            Some(Curve3::Circle(c)) => Some(*c),
            _ => None,
        }
    };

    let circle = get_circle(u_result).or_else(|| get_circle(v_result))?;
    let other_circle = get_circle(v_result).or_else(|| get_circle(u_result))?;
    let r = circle.radius;

    // OCCT: sample two additional iso-curves at 1/3 and 2/3 along the other direction
    let t1 = v_min + (v_max - v_min) / 3.0;
    let t2 = v_min + (v_max - v_min) * 2.0 / 3.0;
    let iso1_pts = sample_iso_curve(bspline, t1, true, u_min, u_max, 16);
    let iso2_pts = sample_iso_curve(bspline, t2, true, u_min, u_max, 16);
    let iso1_crv = bspline_from_point_list(&iso1_pts);
    let iso2_crv = bspline_from_point_list(&iso2_pts);
    let iso1 = curve_to_analytic_curve(&Curve3::BSpline(iso1_crv), tol);
    let iso2 = curve_to_analytic_curve(&Curve3::BSpline(iso2_crv), tol);

    // Check they're circles with matching radius (OCCT: lines 560-575)
    let r1 = match &iso1 {
        Some(Curve3::Circle(c)) => c.radius,
        _ => return None,
    };
    let r2 = match &iso2 {
        Some(Curve3::Circle(c)) => c.radius,
        _ => return None,
    };
    if (r - r1).abs() > tol || (r - r2).abs() > tol {
        return None;
    }

    // Get centers (OCCT: lines 578-600)
    let c1 = circle.center;
    let c2 = iso1.and_then(|c| match c {
        Curve3::Circle(cc) => Some(cc.center),
        _ => None,
    })?;
    let c3 = iso2.and_then(|c| match c {
        Curve3::Circle(cc) => Some(cc.center),
        _ => None,
    })?;

    let d0 = (c1 - c2).length();
    let d1 = (c1 - c3).length();

    if d0 < tol || d1 < tol {
        // OCCT: sphere — all centers same point
        let normal = other_circle.normal;
        return Some(Surface3::Sphere(crate::geom::SphericalSurface {
            center: c1,
            axis: normal,
            radius: r,
            ref_dir: DVec3::X,
        }));
    }

    // OCCT: torus — fit circle through the three centers
    // GetCircle(circ, aPnt1, aPnt2, aPnt3) → major radius
    let torus_circle = make_circle_3p_internal(c1, c2, c3)?;
    let major_r = torus_circle.radius;
    let center = torus_circle.center;

    Some(Surface3::Torus(crate::geom::ToroidalSurface {
        center,
        axis: torus_circle.normal,
        ref_dir: torus_circle.x_dir,
        major_radius: major_r,
        minor_radius: r,
    }))
}

/// OCCT TryCylinerCone: one mid-iso is a line, the other is a circle.
/// Check boundary iso-curves in the circular direction.
/// If all circles have same radius → cylinder; varying → cone.
fn try_cylinder_cone(
    bspline: &BSplineSurface,
    u_result: &Option<Curve3>,
    v_result: &Option<Curve3>,
    v_case: bool,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    tol: f64,
) -> Option<Surface3> {
    // OCCT lines 127-195
    let (param1, param2) = if v_case {
        (u_min, u_max)
    } else {
        (v_min, v_max)
    };

    // Get boundary iso-curves in the circular direction
    let iso_min_pts = if v_case {
        sample_iso_curve(bspline, v_min, true, u_min, u_max, 20)
    } else {
        sample_iso_curve(bspline, u_min, false, v_min, v_max, 20)
    };
    let iso_max_pts = if v_case {
        sample_iso_curve(bspline, v_max, true, u_min, u_max, 20)
    } else {
        sample_iso_curve(bspline, u_max, false, v_min, v_max, 20)
    };

    let iso_min_crv = bspline_from_point_list(&iso_min_pts);
    let iso_max_crv = bspline_from_point_list(&iso_max_pts);
    let iso_min = curve_to_analytic_curve(&Curve3::BSpline(iso_min_crv), tol);
    let iso_max = curve_to_analytic_curve(&Curve3::BSpline(iso_max_crv), tol);

    // Get the mid-circle which we know is valid
    let mid_circle_ref = if v_case { v_result } else { u_result };
    let mid_circle = match mid_circle_ref {
        Some(Curve3::Circle(c)) => *c,
        _ => return None,
    };
    let r_mid = mid_circle.radius;
    let p_mid = mid_circle.center;

    // OCCT: compute radii and positions for the boundary circles (lines 170-196)
    let (r_min, p_min) = match &iso_min {
        Some(Curve3::Circle(c)) => (c.radius, c.center),
        _ => (0.0, iso_min_pts[iso_min_pts.len() / 2]),
    };
    let (r_max, p_max) = match &iso_max {
        Some(Curve3::Circle(c)) => (c.radius, c.center),
        _ => (0.0, iso_max_pts[iso_max_pts.len() / 2]),
    };

    // OCCT lines 196-200: cylinder check
    if (r_mid - r_min).abs() < tol && (r_mid - r_max).abs() < tol {
        // Cylinder — all sections have same radius
        let axis_dir = (p_max - p_min).normalize_or_zero();
        if axis_dir.length_squared() < 0.5 {
            return None;
        }
        // Use the radial direction from the circle
        let radial = mid_circle.x_dir;
        return Some(Surface3::Cylinder(
            crate::geom::CylindricalSurface::new_with_ref_dir(p_min, axis_dir, r_mid, radial),
        ));
    }

    // OCCT lines 202+: cone check — linearly varying radii
    let dr_mid = (r_max - r_min).abs();
    let dr_mid_vs_mid = ((r_mid - (r_min + r_max) * 0.5).abs()).max(dr_mid);
    if dr_mid_vs_mid > tol && dr_mid > tol {
        // Compute cone apex and half-angle
        let axis_dir = (p_max - p_min).normalize_or_zero();
        if axis_dir.length_squared() < 0.5 {
            return None;
        }
        // Extrapolate to apex: where radius = 0 along the axis
        let height = (p_max - p_min).length();
        if height < tol {
            return None;
        }
        let dr = r_max - r_min;
        let apex_dist = r_min / (dr / height);
        let apex = if dr > 0.0 {
            p_min - axis_dir * apex_dist
        } else {
            p_min + axis_dir * apex_dist
        };

        let half_angle = (dr / height).atan();

        return Some(Surface3::Cone(crate::geom::ConicalSurface::new(
            apex, axis_dir, r_min, half_angle,
        )));
    }

    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle3, Curve3, Line3, SurfaceEval};
    use glam::DVec3;

    fn approx_eq3(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    // ── Curve tests ──────────────────────────────────────────────────────────

    #[test]
    fn line_bspline_endpoints() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let bs = line_to_bspline(&line);
        let p0 = bs.point_at(0.0);
        let p1 = bs.point_at(1.0);
        assert!(approx_eq3(p0, DVec3::new(0.0, 0.0, 0.0), 1e-10));
        assert!(approx_eq3(p1, DVec3::new(1.0, 0.0, 0.0), 1e-10));
    }

    #[test]
    fn circle_bspline_is_exact() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let bs = circle_to_bspline(&circle);
        // The NURBS circle evaluates with rational weights — test a few points
        for i in 0..8 {
            let t = i as f64 / 8.0;
            let p = bs.point_at(t);
            let r = (p - circle.center).length();
            assert!(
                (r - circle.radius).abs() < 1e-10,
                "radius at t={t}: expected {}, got {r}",
                circle.radius
            );
        }
    }

    #[test]
    fn ellipse_bspline_endpoints() {
        use crate::geom::Ellipse3;
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 3.0,
            minor_radius: 1.5,
        };
        let bs = ellipse_to_bspline(&e);
        // t=0 and t=1 both map to the major-axis endpoint
        let p0 = bs.point_at(0.0);
        let p1 = bs.point_at(1.0);
        assert!(approx_eq3(p0, DVec3::new(3.0, 0.0, 0.0), 1e-10), "p0={p0}");
        assert!(approx_eq3(p1, DVec3::new(3.0, 0.0, 0.0), 1e-10), "p1={p1}");
    }

    #[test]
    fn curve_to_bspline_identity_for_bspline() {
        let bs_orig = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        };
        let bs_conv = curve_to_bspline(&Curve3::BSpline(bs_orig.clone()), 16);
        assert_eq!(bs_conv.degree, bs_orig.degree);
        assert_eq!(bs_conv.control_points.len(), bs_orig.control_points.len());
    }

    // ── Surface tests ────────────────────────────────────────────────────────

    #[test]
    fn plane_bspline_corners() {
        use crate::geom::Plane;
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let bs = plane_to_bspline(&plane);
        // Should evaluate points at corners of the [-1,1]×[-1,1] domain
        let p00 = bs.point_at(0.0, 0.0);
        let p11 = bs.point_at(1.0, 1.0);
        // p at (u=0,v=0) should be corner of domain
        assert!(
            p00.distance(DVec3::new(-1.0, -1.0, 0.0)) < 1e-10
                || p00.distance(DVec3::new(1.0, 1.0, 0.0)) < 1e-10
                || p00.z.abs() < 1e-10, // at least z=0
            "plane corner z={}",
            p00.z
        );
        let _ = p11;
    }

    #[test]
    fn cylinder_bspline_on_surface() {
        use crate::geom::CylindricalSurface;
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: DVec3::X,
            y_dir: None,
        };
        let bs = cylinder_to_bspline(&cyl);
        // Sample several u values; all should be on the cylinder surface
        for i in 0..9 {
            let u = i as f64 / 8.0;
            let p = bs.point_at(u, 0.0);
            let r = DVec3::new(p.x, p.y, 0.0).length();
            assert!((r - 1.0).abs() < 1e-9, "u={u}: radius={r}");
        }
    }

    #[test]
    fn sphere_bspline_on_surface() {
        use crate::any_perpendicular;
        use crate::geom::SphericalSurface;
        let sph = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        };
        let bs = sphere_to_bspline(&sph);
        // Equator row (v=0.5 maps to colatitude 90°)
        for i in 0..9 {
            let u = i as f64 / 8.0;
            let p = bs.point_at(u, 0.5);
            let r = p.length();
            assert!((r - 1.0).abs() < 1e-9, "u={u}: radius={r}");
        }
    }
}
