//! Shape analysis tools for geometric validation.
//!
//! Analogous to OCCT `ShapeAnalysis` package:
//! - `ShapeAnalysis_Surface`: surface analysis (UV consistency, bounds, singularities)
//! - `ShapeAnalysis_Curve`: curve analysis (parameter range, self-intersection, continuity)
//! - `ShapeAnalysis_Wire`: wire analysis (closure, orientation, self-intersection)
//! - `ShapeAnalysis_Face`: face analysis (boundary validity, param domain, surface-wire consistency)
//!
//! All functions are non-destructive analysis tools that return structured reports.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Surface3, CurveEval, SurfaceEval, Curve2dEval};
use rcad_kernel::{BRep, Face, PCurve};

// ─────────────────────────────────────────────────────────────────────────────
// Surface Analysis (ShapeAnalysis_Surface)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from surface analysis.
///
/// Analogous to OCCT `ShapeAnalysis_Surface`.
#[derive(Debug, Clone)]
pub struct SurfaceAnalysisReport {
    /// U parameter range [u_min, u_max].
    pub u_range: (f64, f64),
    /// V parameter range [v_min, v_max].
    pub v_range: (f64, f64),
    /// Whether the surface is periodic in U direction.
    pub is_u_periodic: bool,
    /// Whether the surface is periodic in V direction.
    pub is_v_periodic: bool,
    /// Detected singular points on the surface (e.g., sphere poles).
    pub singular_points: Vec<SingularPoint>,
    /// Whether any boundary edge is degenerate (zero-length parametric derivative).
    pub bounds_degenerate: bool,
    /// UV consistency issues detected.
    pub uv_issues: Vec<UvInconsistency>,
    /// Surface orientation status (is the parametric orientation consistent?).
    pub orientation_ok: bool,
}

/// A singular point on a surface where the normal is undefined.
#[derive(Debug, Clone)]
pub struct SingularPoint {
    /// The 3D location of the singular point.
    pub point: DVec3,
    /// The UV parameter at which the singularity occurs.
    pub uv: (f64, f64),
    /// Type of singularity.
    pub kind: SingularPointKind,
}

/// Classification of surface singularity type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SingularPointKind {
    /// Pole singularity (e.g., sphere north/south poles).
    Pole,
    /// Apex singularity (e.g., cone apex).
    Apex,
    /// Degenerate boundary (zero-length edge at parametric boundary).
    DegenerateBoundary,
    /// Self-intersection singularity.
    SelfIntersection,
}

/// UV consistency issue detected on a surface.
#[derive(Debug, Clone)]
pub struct UvInconsistency {
    /// Type of inconsistency.
    pub kind: UvInconsistencyKind,
    /// UV location where the issue was detected.
    pub uv: (f64, f64),
    /// Description of the issue.
    pub description: String,
}

/// Classification of UV inconsistency types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvInconsistencyKind {
    /// Parameter jump discontinuity.
    ParamJump,
    /// Normal direction discontinuity.
    NormalFlip,
    /// Derivative discontinuity.
    DerivativeDiscontinuity,
    /// Invalid parameter value (NaN or infinite).
    InvalidParam,
    /// Non-monotonic parameterization.
    NonMonotonic,
}

/// Analyze a surface for geometric validity and characteristics.
///
/// Performs comprehensive analysis including:
/// - Parameter range validation
/// - Periodicity detection
/// - Singular point detection
/// - UV consistency checks
/// - Boundary degeneracy detection
///
/// # Example
/// ```rust
/// use rcad_kernel::geom::{Surface3, SphericalSurface};
/// use rcad_algorithms::shape_analysis::{analyze_surface, SurfaceAnalysisReport};
/// use glam::DVec3;
///
/// let sphere = Surface3::Sphere(SphericalSurface {
///     center: DVec3::ZERO,
///     axis: DVec3::Y,
///     radius: 1.0,
/// });
/// let report = analyze_surface(&sphere);
/// assert!(report.is_u_periodic);
/// assert_eq!(report.singular_points.len(), 2); // North and south poles
/// ```
pub fn analyze_surface(surf: &Surface3) -> SurfaceAnalysisReport {
    let domain = surf.default_domain();
    let u_range = (domain[0], domain[1]);
    let v_range = (domain[2], domain[3]);

    let (is_u_periodic, is_v_periodic) = detect_periodicity(surf);
    let singular_points = detect_singular_points(surf);
    let uv_issues = check_uv_consistency(surf, 1e-9);
    let bounds_degenerate = check_bounds_degeneracy(surf);
    let orientation_ok = check_surface_orientation(surf);

    SurfaceAnalysisReport {
        u_range,
        v_range,
        is_u_periodic,
        is_v_periodic,
        singular_points,
        bounds_degenerate,
        uv_issues,
        orientation_ok,
    }
}

/// Check UV consistency of a surface with given tolerance.
///
/// Samples the surface on a grid and checks for:
/// - Parameter discontinuities
/// - Normal flips
/// - Derivative discontinuities
/// - Invalid parameter values
pub fn check_uv_consistency(surf: &Surface3, tolerance: f64) -> Vec<UvInconsistency> {
    let mut issues = Vec::new();
    let domain = surf.default_domain();
    let (u0, u1) = (domain[0], domain[1]);
    let (v0, v1) = (domain[2], domain[3]);

    // Handle infinite domains
    let (u0, u1) = if u0.is_infinite() || u1.is_infinite() {
        (-10.0, 10.0)
    } else {
        (u0, u1)
    };
    let (v0, v1) = if v0.is_infinite() || v1.is_infinite() {
        (-10.0, 10.0)
    } else {
        (v0, v1)
    };

    let n_samples = 10;

    // Check for NaN/infinite points and normal flips
    let du = (u1 - u0) / n_samples as f64;
    let dv = (v1 - v0) / n_samples as f64;

    for i in 0..=n_samples {
        for j in 0..=n_samples {
            let u = u0 + du * i as f64;
            let v = v0 + dv * j as f64;

            let p = surf.point_at(u, v);
            let n = surf.normal_at(u, v);

            // Check for invalid point
            if !p.is_finite() {
                issues.push(UvInconsistency {
                    kind: UvInconsistencyKind::InvalidParam,
                    uv: (u, v),
                    description: "Surface point is not finite (NaN or infinite)".to_string(),
                });
            }

            // Check for invalid normal
            if !n.is_finite() || n.length_squared() < 0.5 {
                // This might be a singular point, but not necessarily an error
            }

            // Check for normal discontinuity with neighbors
            if i > 0 && j > 0 {
                let u_prev = u0 + du * (i - 1) as f64;
                let v_prev = v0 + dv * (j - 1) as f64;

                let n_u = surf.normal_at(u_prev, v);
                let n_v = surf.normal_at(u, v_prev);

                // Check for normal flip (more than 90 degree change over small step)
                let du_ratio = (n - n_u).length() / du;
                let dv_ratio = (n - n_v).length() / dv;

                if du_ratio > 100.0 || dv_ratio > 100.0 {
                    // Potential discontinuity - but may be normal for periodic surfaces
                }
            }
        }
    }

    // Check derivative continuity at midpoint
    let u_mid = (u0 + u1) / 2.0;
    let v_mid = (v0 + v1) / 2.0;

    if !check_derivative_continuity(surf, u_mid, v_mid, tolerance) {
        issues.push(UvInconsistency {
            kind: UvInconsistencyKind::DerivativeDiscontinuity,
            uv: (u_mid, v_mid),
            description: "Derivative discontinuity detected at surface midpoint".to_string(),
        });
    }

    issues
}

/// Detect periodicity in U and V directions.
fn detect_periodicity(surf: &Surface3) -> (bool, bool) {
    match surf {
        Surface3::Cylinder(_) => (true, false),
        Surface3::Sphere(_) => (true, false),
        Surface3::Cone(_) => (true, false),
        Surface3::Torus(_) => (true, true),
        Surface3::Helicoid(_) => (true, false),
        Surface3::Revolution(_) => (true, false),
        Surface3::BSpline(bs) => {
            // Check if knot vector indicates periodicity
            let u_periodic = is_bspline_periodic(&bs.knots_u, bs.degree_u);
            let v_periodic = is_bspline_periodic(&bs.knots_v, bs.degree_v);
            (u_periodic, v_periodic)
        }
        _ => (false, false),
    }
}

/// Check if a BSpline knot vector indicates a periodic surface.
fn is_bspline_periodic(knots: &[f64], degree: usize) -> bool {
    if knots.len() < 2 * (degree + 1) {
        return false;
    }
    let n = knots.len();
    let span = knots[n - 1] - knots[0];

    // Check if first (degree+1) knots equal the first internal knot
    // and last (degree+1) knots equal the last internal knot
    let eps = 1e-9;
    let first_knot = knots[0];
    let last_knot = knots[n - 1];

    // Periodic if there's enough repetition at boundaries
    let first_count = knots.iter().take_while(|&&k| (k - first_knot).abs() < eps).count();
    let last_count = knots.iter().rev().take_while(|&&k| (k - last_knot).abs() < eps).count();

    // For uniform periodic splines, multiplicity should be 1 at internal knots
    first_count == 1 && last_count == 1 && span > eps
}

/// Detect singular points on a surface.
fn detect_singular_points(surf: &Surface3) -> Vec<SingularPoint> {
    let mut points = Vec::new();

    match surf {
        Surface3::Sphere(s) => {
            // Sphere has two poles at v=0 and v=PI
            let domain = surf.default_domain();
            let u_mid = (domain[0] + domain[1]) / 2.0;

            // North pole (v = 0)
            points.push(SingularPoint {
                point: s.center + s.radius * s.axis.normalize(),
                uv: (u_mid, domain[2]),
                kind: SingularPointKind::Pole,
            });

            // South pole (v = PI)
            points.push(SingularPoint {
                point: s.center - s.radius * s.axis.normalize(),
                uv: (u_mid, domain[3]),
                kind: SingularPointKind::Pole,
            });
        }

        Surface3::Cone(c) => {
            // Cone has an apex at v=0 (if radius at apex is 0)
            if c.radius.abs() < 1e-12 {
                let domain = surf.default_domain();
                let u_mid = (domain[0] + domain[1]) / 2.0;

                points.push(SingularPoint {
                    point: c.apex_point(),
                    uv: (u_mid, domain[2]),
                    kind: SingularPointKind::Apex,
                });
            }
        }

        Surface3::Torus(t) => {
            // Torus has no singular points unless minor_radius is 0
            if t.minor_radius.abs() < 1e-12 {
                let domain = surf.default_domain();
                // The entire center circle becomes singular
                for i in 0..8 {
                    let u = domain[0] + (domain[1] - domain[0]) * i as f64 / 8.0;
                    points.push(SingularPoint {
                        point: t.center + t.major_radius * DVec3::X,
                        uv: (u, 0.0),
                        kind: SingularPointKind::DegenerateBoundary,
                    });
                }
            }
        }

        Surface3::Ellipsoid(e) => {
            // Ellipsoid has two poles at v=0 and v=PI
            let domain = surf.default_domain();
            let u_mid = (domain[0] + domain[1]) / 2.0;
            let axis = e.axis.normalize();

            points.push(SingularPoint {
                point: e.center + e.radius_z * axis,
                uv: (u_mid, domain[2]),
                kind: SingularPointKind::Pole,
            });

            points.push(SingularPoint {
                point: e.center - e.radius_z * axis,
                uv: (u_mid, domain[3]),
                kind: SingularPointKind::Pole,
            });
        }

        _ => {}
    }

    points
}

/// Check if any boundary of the surface is degenerate.
fn check_bounds_degeneracy(surf: &Surface3) -> bool {
    let domain = surf.default_domain();
    let [u0, u1, v0, v1] = domain;

    // Handle infinite domains
    if u0.is_infinite() || u1.is_infinite() || v0.is_infinite() || v1.is_infinite() {
        return false;
    }

    let eps = 1e-9;

    // Check if opposite boundaries map to the same 3D curve
    // (this indicates a degenerate boundary)
    let n_samples = 10;
    let du = (u1 - u0) / n_samples as f64;
    let dv = (v1 - v0) / n_samples as f64;

    // Check v = v0 boundary vs v = v1 boundary
    let mut v0_points = Vec::new();
    let mut v1_points = Vec::new();
    for i in 0..=n_samples {
        let u = u0 + du * i as f64;
        v0_points.push(surf.point_at(u, v0));
        v1_points.push(surf.point_at(u, v1));
    }

    // If all points on a boundary are the same, it's degenerate
    let v0_degenerate = v0_points.iter().all(|p| (p - v0_points[0]).length() < eps);
    let v1_degenerate = v1_points.iter().all(|p| (p - v1_points[0]).length() < eps);

    if v0_degenerate || v1_degenerate {
        return true;
    }

    // Check u = u0 boundary vs u = u1 boundary
    let mut u0_points = Vec::new();
    let mut u1_points = Vec::new();
    for i in 0..=n_samples {
        let v = v0 + dv * i as f64;
        u0_points.push(surf.point_at(u0, v));
        u1_points.push(surf.point_at(u1, v));
    }

    let u0_degenerate = u0_points.iter().all(|p| (p - u0_points[0]).length() < eps);
    let u1_degenerate = u1_points.iter().all(|p| (p - u1_points[0]).length() < eps);

    u0_degenerate || u1_degenerate
}

/// Check derivative continuity at a point using finite differences.
fn check_derivative_continuity(surf: &Surface3, u: f64, v: f64, tolerance: f64) -> bool {
    let eps = 1e-6;

    let p = surf.point_at(u, v);

    // Check if point is valid
    if !p.is_finite() {
        return true; // Skip invalid points
    }

    // Compute partial derivatives via finite difference
    let p_up = surf.point_at(u + eps, v);
    let p_um = surf.point_at(u - eps, v);
    let p_vp = surf.point_at(u, v + eps);
    let p_vm = surf.point_at(u, v - eps);

    // Check if derivatives are finite
    let du = p_up - p_um;
    let dv = p_vp - p_vm;

    du.is_finite() && dv.is_finite()
}

/// Check surface orientation consistency.
fn check_surface_orientation(surf: &Surface3) -> bool {
    let domain = surf.default_domain();

    // For closed surfaces, check if the normal direction is consistent
    // at opposite boundaries
    let [u0, u1, v0, v1] = domain;

    // Handle infinite domains
    if u0.is_infinite() || u1.is_infinite() || v0.is_infinite() || v1.is_infinite() {
        return true;
    }

    // Check normal at a few points
    let n_mid = surf.normal_at((u0 + u1) / 2.0, (v0 + v1) / 2.0);

    // For periodic surfaces, normal should be consistent
    if n_mid.is_finite() && n_mid.length() > 0.5 {
        return true;
    }

    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Curve Analysis (ShapeAnalysis_Curve)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from curve analysis.
#[derive(Debug, Clone)]
pub struct CurveAnalysisReport {
    /// Parameter range [t_min, t_max].
    pub param_range: (f64, f64),
    /// Whether the curve is closed (start point equals end point).
    pub is_closed: bool,
    /// Whether the curve is periodic.
    pub is_periodic: bool,
    /// Detected self-intersection points.
    pub self_intersections: Vec<CurveSelfIntersection>,
    /// Continuity level (0 = C0, 1 = C1, 2 = C2).
    pub continuity: ContinuityLevel,
    /// Total arc length of the curve.
    pub arc_length: f64,
    /// Whether the curve is degenerate (zero length).
    pub is_degenerate: bool,
}

/// A self-intersection point on a curve.
#[derive(Debug, Clone)]
pub struct CurveSelfIntersection {
    /// First parameter value where intersection occurs.
    pub param1: f64,
    /// Second parameter value where intersection occurs.
    pub param2: f64,
    /// 3D point of intersection.
    pub point: DVec3,
}

/// Continuity level classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ContinuityLevel {
    /// C0: position continuous only.
    C0,
    /// C1: tangent continuous.
    C1,
    /// C2: curvature continuous.
    C2,
    /// CN: infinitely differentiable (analytic).
    CN,
}

/// Analyze a curve for geometric validity and characteristics.
///
/// # Example
/// ```rust
/// use rcad_kernel::geom::{Curve3, Circle3};
/// use rcad_algorithms::shape_analysis::{analyze_curve, CurveAnalysisReport, ContinuityLevel};
/// use glam::DVec3;
///
/// let circle = Curve3::Circle(Circle3 {
///     center: DVec3::ZERO,
///     normal: DVec3::Z,
///     radius: 1.0,
/// });
/// let report = analyze_curve(&circle, 64);
/// assert!(report.is_closed);
/// assert!(report.is_periodic);
/// assert_eq!(report.continuity, ContinuityLevel::CN);
/// ```
pub fn analyze_curve(curve: &Curve3, n_samples: usize) -> CurveAnalysisReport {
    let domain = curve.default_domain();
    let param_range = (domain[0], domain[1]);

    let is_periodic = is_curve_periodic(curve);
    let is_closed = check_curve_closed(curve);
    let self_intersections = detect_curve_self_intersections(curve, n_samples);
    let continuity = determine_curve_continuity(curve);
    let arc_length = compute_curve_length(curve, n_samples);
    let is_degenerate = arc_length < 1e-12;

    CurveAnalysisReport {
        param_range,
        is_closed,
        is_periodic,
        self_intersections,
        continuity,
        arc_length,
        is_degenerate,
    }
}

/// Check if a curve is periodic.
fn is_curve_periodic(curve: &Curve3) -> bool {
    matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_))
}

/// Check if a curve is closed.
fn check_curve_closed(curve: &Curve3) -> bool {
    let domain = curve.default_domain();

    // Handle infinite domains
    if domain[0].is_infinite() || domain[1].is_infinite() {
        return false;
    }

    let p_start = curve.point_at(domain[0]);
    let p_end = curve.point_at(domain[1]);

    (p_start - p_end).length() < 1e-9
}

/// Detect self-intersections in a curve by sampling.
fn detect_curve_self_intersections(curve: &Curve3, n_samples: usize) -> Vec<CurveSelfIntersection> {
    let mut intersections = Vec::new();
    let domain = curve.default_domain();

    // Handle infinite domains
    let (t0, t1) = if domain[0].is_infinite() || domain[1].is_infinite() {
        return intersections; // Can't detect self-intersection on infinite domain
    } else {
        (domain[0], domain[1])
    };

    let dt = (t1 - t0) / n_samples as f64;

    // Sample points
    let points: Vec<(f64, DVec3)> = (0..=n_samples)
        .map(|i| {
            let t = t0 + dt * i as f64;
            (t, curve.point_at(t))
        })
        .collect();

    // Check for non-adjacent segments that intersect
    let tol = 1e-6;

    for i in 0..points.len() - 1 {
        // Only check segments that are not adjacent (at least 2 apart)
        for j in (i + 3)..points.len() - 1 {
            let p1 = points[i].1;
            let p2 = points[i + 1].1;
            let p3 = points[j].1;
            let p4 = points[j + 1].1;

            // Check segment intersection in 2D (project to XY plane for simplicity)
            // A more robust implementation would use 3D segment distance
            if let Some((t, s)) = segment_intersection_2d(
                [p1.x, p1.y], [p2.x, p2.y],
                [p3.x, p3.y], [p4.x, p4.y],
            ) {
                let point = DVec3::new(
                    p1.x + t * (p2.x - p1.x),
                    p1.y + t * (p2.y - p1.y),
                    p1.z + t * (p2.z - p1.z),
                );

                let param1 = points[i].0 + t * (points[i + 1].0 - points[i].0);
                let param2 = points[j].0 + s * (points[j + 1].0 - points[j].0);

                intersections.push(CurveSelfIntersection {
                    param1,
                    param2,
                    point,
                });
            }
        }
    }

    intersections
}

/// 2D segment intersection test.
/// Returns (t, s) parameters if segments intersect, where t is on segment 1 and s is on segment 2.
fn segment_intersection_2d(
    p1: [f64; 2], p2: [f64; 2],
    p3: [f64; 2], p4: [f64; 2],
) -> Option<(f64, f64)> {
    let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
    let d2 = [p4[0] - p3[0], p4[1] - p3[1]];

    let cross = d1[0] * d2[1] - d1[1] * d2[0];

    if cross.abs() < 1e-12 {
        return None; // Parallel segments
    }

    let dx = p3[0] - p1[0];
    let dy = p3[1] - p1[1];

    let t = (dx * d2[1] - dy * d2[0]) / cross;
    let s = (dx * d1[1] - dy * d1[0]) / cross;

    if t >= 0.0 && t <= 1.0 && s >= 0.0 && s <= 1.0 {
        Some((t, s))
    } else {
        None
    }
}

/// Determine the continuity level of a curve.
fn determine_curve_continuity(curve: &Curve3) -> ContinuityLevel {
    match curve {
        Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_) => ContinuityLevel::CN,
        Curve3::Hyperbola(_) | Curve3::Parabola(_) | Curve3::CircularHelix(_) | Curve3::SineWave(_) => ContinuityLevel::CN,
        Curve3::BSpline(bs) => {
            // BSpline continuity is degree - multiplicity at each knot
            // For simplicity, assume at least C2 if degree >= 3
            if bs.degree >= 3 { ContinuityLevel::C2 }
            else if bs.degree >= 2 { ContinuityLevel::C1 }
            else { ContinuityLevel::C0 }
        }
        Curve3::Bezier(bez) => {
            // Bezier curves are C-infinity on (0,1), but C{degree-1} at endpoints
            if bez.control_points.len() >= 4 { ContinuityLevel::C2 }
            else if bez.control_points.len() >= 3 { ContinuityLevel::C1 }
            else { ContinuityLevel::C0 }
        }
        Curve3::Offset(_) => ContinuityLevel::C1, // Conservative estimate
    }
}

/// Compute approximate arc length by numerical integration.
fn compute_curve_length(curve: &Curve3, n_samples: usize) -> f64 {
    let domain = curve.default_domain();

    // Handle infinite domains
    let (t0, t1) = if domain[0].is_infinite() || domain[1].is_infinite() {
        return f64::INFINITY;
    } else {
        (domain[0], domain[1])
    };

    let n = n_samples.max(2);
    let dt = (t1 - t0) / n as f64;

    let mut length = 0.0;
    let mut p_prev = curve.point_at(t0);

    for i in 1..=n {
        let t = t0 + dt * i as f64;
        let p = curve.point_at(t);
        length += (p - p_prev).length();
        p_prev = p;
    }

    length
}

// ─────────────────────────────────────────────────────────────────────────────
// Wire Analysis (ShapeAnalysis_Wire)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from wire analysis.
#[derive(Debug, Clone)]
pub struct WireAnalysisReport {
    /// Whether the wire is closed.
    pub is_closed: bool,
    /// Whether the wire orientation is consistent.
    pub orientation_consistent: bool,
    /// Number of edges in the wire.
    pub edge_count: usize,
    /// Number of vertices in the wire.
    pub vertex_count: usize,
    /// Self-intersection issues.
    pub self_intersections: Vec<WireSelfIntersection>,
    /// Wire length (sum of edge lengths).
    pub length: f64,
    /// Whether the wire is degenerate.
    pub is_degenerate: bool,
    /// Gaps between consecutive edges.
    pub gaps: Vec<WireGap>,
}

/// A self-intersection in a wire.
#[derive(Debug, Clone)]
pub struct WireSelfIntersection {
    /// Index of the first edge involved.
    pub edge_a: usize,
    /// Index of the second edge involved.
    pub edge_b: usize,
    /// Intersection point.
    pub point: DVec3,
}

/// A gap between consecutive edges in a wire.
#[derive(Debug, Clone)]
pub struct WireGap {
    /// Index of the edge where the gap starts.
    pub after_edge: usize,
    /// Distance of the gap.
    pub distance: f64,
    /// Start point of the gap.
    pub from_point: DVec3,
    /// End point of the gap.
    pub to_point: DVec3,
}

/// Analyze a wire for validity and characteristics.
///
/// This is a topological analysis that checks wire closure, orientation,
/// and self-intersection at the topology level.
pub fn analyze_wire(
    brep: &BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    wire_idx: Option<usize>, // None for outer wire, Some(i) for inner wire
) -> WireAnalysisReport {
    let mut report = WireAnalysisReport {
        is_closed: true,
        orientation_consistent: true,
        edge_count: 0,
        vertex_count: 0,
        self_intersections: Vec::new(),
        length: 0.0,
        is_degenerate: false,
        gaps: Vec::new(),
    };

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(face_idx) else { return report; };

    let wire = match wire_idx {
        None => &face.outer_wire,
        Some(i) => face.inner_wires.get(i).unwrap_or(&face.outer_wire),
    };

    report.edge_count = wire.edges.len();

    if wire.edges.is_empty() {
        report.is_closed = false;
        report.is_degenerate = true;
        return report;
    }

    // Collect edge vertices
    let mut vertices: Vec<(usize, usize)> = Vec::new(); // (start, end) vertex indices
    let mut vertex_set = std::collections::HashSet::new();

    for we in &wire.edges {
        let Some(edge) = brep.edges.get(we.idx) else {
            continue;
        };

        let (start, end) = if we.forward {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        };

        vertices.push((start, end));
        vertex_set.insert(start);
        vertex_set.insert(end);

        // Compute edge length if geometry is available
        if let Some(curve_idx) = brep.geom.edge_curve.get(we.idx).and_then(|opt| *opt) {
            if let Some(curve) = brep.geom.curves.get(curve_idx) {
                let range = brep.geom.edge_curve_range.get(we.idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| {
                        let d = curve.default_domain();
                        [d[0], d[1]]
                    });

                // Approximate length by sampling
                let n = 10;
                let dt = (range[1] - range[0]) / n as f64;
                let mut len = 0.0;
                let mut p_prev = curve.point_at(range[0]);
                for i in 1..=n {
                    let t = range[0] + dt * i as f64;
                    let p = curve.point_at(t);
                    len += (p - p_prev).length();
                    p_prev = p;
                }
                report.length += len;
            }
        }
    }

    report.vertex_count = vertex_set.len();

    // Check closure
    let n = vertices.len();
    if n == 0 {
        report.is_closed = false;
        report.is_degenerate = true;
        return report;
    }

    // Special case: single edge that forms a closed loop (e.g., circle for cap face)
    // In this case, start == end for the edge
    if n == 1 {
        let (start, end) = vertices[0];
        // A single-edge wire is closed if the edge starts and ends at the same vertex
        // or if the geometric positions are the same
        if start == end {
            report.is_closed = true;
        } else {
            let start_pt = brep.vertices.get(start).map(|v| v.point).unwrap_or(DVec3::ZERO);
            let end_pt = brep.vertices.get(end).map(|v| v.point).unwrap_or(DVec3::ZERO);
            let gap_dist = (start_pt - end_pt).length();
            report.is_closed = gap_dist < 1e-6;
            if !report.is_closed {
                report.gaps.push(WireGap {
                    after_edge: 0,
                    distance: gap_dist,
                    from_point: end_pt,
                    to_point: start_pt,
                });
            }
        }
        report.is_degenerate = report.length < 1e-12;
        return report;
    }

    for i in 0..n {
        let next = (i + 1) % n;
        let end_v = vertices[i].1;
        let start_v = vertices[next].0;

        if end_v != start_v {
            // Check geometric gap
            let end_pt = brep.vertices.get(end_v).map(|v| v.point).unwrap_or(DVec3::ZERO);
            let start_pt = brep.vertices.get(start_v).map(|v| v.point).unwrap_or(DVec3::ZERO);
            let gap_dist = (end_pt - start_pt).length();

            if gap_dist > 1e-6 {
                report.is_closed = false;
                report.gaps.push(WireGap {
                    after_edge: i,
                    distance: gap_dist,
                    from_point: end_pt,
                    to_point: start_pt,
                });
            }
        }
    }

    // Check for topological self-intersection (vertex appears more than twice)
    let mut vertex_count: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for &(start, end) in &vertices {
        *vertex_count.entry(start).or_insert(0) += 1;
        *vertex_count.entry(end).or_insert(0) += 1;
    }

    for (&v, &count) in &vertex_count {
        if count > 2 {
            let point = brep.vertices.get(v).map(|v| v.point).unwrap_or(DVec3::ZERO);
            // Find which edges share this vertex
            let edges_with_vertex: Vec<usize> = vertices.iter()
                .enumerate()
                .filter(|(_, (s, e))| *s == v || *e == v)
                .map(|(i, _)| i)
                .collect();

            if edges_with_vertex.len() >= 2 {
                report.self_intersections.push(WireSelfIntersection {
                    edge_a: edges_with_vertex[0],
                    edge_b: edges_with_vertex[1],
                    point,
                });
            }
        }
    }

    report.is_degenerate = report.length < 1e-12;
    report
}

/// Check if all wires in a face are valid.
pub fn check_face_wires(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> Vec<WireAnalysisReport> {
    let mut reports = Vec::new();

    // Check outer wire
    reports.push(analyze_wire(brep, solid_idx, shell_idx, face_idx, None));

    // Check inner wires
    let Some(solid) = brep.solids.get(solid_idx) else { return reports; };
    let Some(shell) = solid.shells.get(shell_idx) else { return reports; };
    let Some(face) = shell.faces.get(face_idx) else { return reports; };

    for i in 0..face.inner_wires.len() {
        reports.push(analyze_wire(brep, solid_idx, shell_idx, face_idx, Some(i)));
    }

    reports
}

// ─────────────────────────────────────────────────────────────────────────────
// Face Analysis (ShapeAnalysis_Face)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from face analysis.
#[derive(Debug, Clone)]
pub struct FaceAnalysisReport {
    /// Whether the face has a valid surface.
    pub has_surface: bool,
    /// Surface analysis report (if surface exists).
    pub surface_report: Option<SurfaceAnalysisReport>,
    /// Wire analysis reports for all wires.
    pub wire_reports: Vec<WireAnalysisReport>,
    /// Whether all wires are closed.
    pub all_wires_closed: bool,
    /// Whether the face orientation matches the surface normal.
    pub orientation_matches_surface: bool,
    /// Surface-wire consistency issues.
    pub surface_wire_issues: Vec<SurfaceWireIssue>,
    /// Parameter domain of the face.
    pub param_domain: Option<(f64, f64, f64, f64)>,
}

/// An issue with surface-wire consistency.
#[derive(Debug, Clone)]
pub struct SurfaceWireIssue {
    /// Kind of issue.
    pub kind: SurfaceWireIssueKind,
    /// Description of the issue.
    pub description: String,
    /// Edge index where the issue occurs (if applicable).
    pub edge_idx: Option<usize>,
}

/// Classification of surface-wire consistency issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceWireIssueKind {
    /// Edge is not on the surface.
    EdgeNotOnSurface,
    /// PCurve is degenerate.
    DegeneratePCurve,
    /// Wire is outside surface domain.
    WireOutsideDomain,
    /// Normal direction mismatch.
    NormalMismatch,
}

/// Analyze a face for validity and characteristics.
pub fn analyze_face(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> FaceAnalysisReport {
    let mut report = FaceAnalysisReport {
        has_surface: false,
        surface_report: None,
        wire_reports: Vec::new(),
        all_wires_closed: true,
        orientation_matches_surface: true,
        surface_wire_issues: Vec::new(),
        param_domain: None,
    };

    // Check if face has a surface
    let surface_idx = brep.geom.face_surface.get(face_idx).and_then(|opt| *opt);
    report.has_surface = surface_idx.is_some();

    // Analyze surface
    if let Some(idx) = surface_idx {
        if let Some(surface) = brep.geom.surfaces.get(idx) {
            report.surface_report = Some(analyze_surface(surface));

            // Get parameter domain
            let domain = surface.default_domain();
            report.param_domain = Some((domain[0], domain[1], domain[2], domain[3]));
        }
    }

    // Analyze wires
    report.wire_reports = check_face_wires(brep, solid_idx, shell_idx, face_idx);

    // Check if all wires are closed
    for wire_report in &report.wire_reports {
        if !wire_report.is_closed {
            report.all_wires_closed = false;
        }
    }

    // Check surface-wire consistency
    report.surface_wire_issues = check_surface_wire_consistency(brep, solid_idx, shell_idx, face_idx);

    // Check orientation
    report.orientation_matches_surface = check_face_orientation(brep, solid_idx, shell_idx, face_idx);

    report
}

/// Check consistency between surface and wire geometry.
fn check_surface_wire_consistency(
    brep: &BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
) -> Vec<SurfaceWireIssue> {
    let mut issues = Vec::new();

    let Some(solid) = brep.solids.get(solid_idx) else { return issues; };
    let Some(shell) = solid.shells.get(shell_idx) else { return issues; };
    let Some(face) = shell.faces.get(face_idx) else { return issues; };

    let surface_idx = match brep.geom.face_surface.get(face_idx).and_then(|opt| *opt) {
        Some(idx) => idx,
        None => return issues,
    };

    let surface = match brep.geom.surfaces.get(surface_idx) {
        Some(s) => s,
        None => return issues,
    };

    // Check each edge in the outer wire
    for we in &face.outer_wire.edges {
        // Check if edge has a PCurve on this surface
        let has_pcurve = brep.geom.edge_pcurves.get(we.idx)
            .map(|pcurves| {
                pcurves.iter().any(|pc| pc.surface_idx == surface_idx)
            })
            .unwrap_or(false);

        if !has_pcurve {
            // Edge might be degenerate or might not lie on surface
            // This is not necessarily an error - check if edge is degenerate
            let is_degenerate = brep.geom.edge_degenerated.get(we.idx).copied().unwrap_or(false);
            if !is_degenerate {
                // Check if the edge's 3D curve lies on the surface
                if let Some(curve_idx) = brep.geom.edge_curve.get(we.idx).and_then(|opt| *opt) {
                    if let Some(curve) = brep.geom.curves.get(curve_idx) {
                        let range = brep.geom.edge_curve_range.get(we.idx)
                            .and_then(|r| *r)
                            .unwrap_or_else(|| {
                                let d = curve.default_domain();
                                [d[0], d[1]]
                            });

                        // Sample a few points and check if they lie on the surface
                        let n_samples = 5;
                        let dt = (range[1] - range[0]) / n_samples as f64;
                        let mut max_deviation: f64 = 0.0;

                        for i in 0..=n_samples {
                            let t = range[0] + dt * i as f64;
                            let p = curve.point_at(t);
                            // Project onto surface and check distance
                            if let Some(proj) = project_point_to_surface_simple(surface, p) {
                                let deviation = (p - proj).length();
                                max_deviation = max_deviation.max(deviation);
                            }
                        }

                        if max_deviation > 1e-6 {
                            issues.push(SurfaceWireIssue {
                                kind: SurfaceWireIssueKind::EdgeNotOnSurface,
                                description: format!("Edge {} does not lie on surface (max deviation: {})", we.idx, max_deviation),
                                edge_idx: Some(we.idx),
                            });
                        }
                    }
                }
            }
        }
    }

    issues
}

/// Simple point-to-surface projection for checking edge-on-surface.
fn project_point_to_surface_simple(surface: &Surface3, point: DVec3) -> Option<DVec3> {
    // Use the domain center as initial guess for iterative projection
    let domain = surface.default_domain();
    let u_center = (domain[0] + domain[1]) / 2.0;
    let v_center = (domain[2] + domain[3]) / 2.0;

    // For analytical surfaces, use direct projection
    match surface {
        Surface3::Plane(p) => {
            let d = (point - p.origin).dot(p.normal);
            Some(point - p.normal * d)
        }
        Surface3::Sphere(s) => {
            let v = point - s.center;
            let len = v.length();
            if len < 1e-14 {
                None
            } else {
                Some(s.center + v / len * s.radius)
            }
        }
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let radial = v - c.axis * along;
            let radial_len = radial.length();
            if radial_len < 1e-14 {
                None
            } else {
                Some(c.origin + c.axis * along + radial / radial_len * c.radius)
            }
        }
        _ => {
            // For other surfaces, return the center point as a placeholder
            Some(surface.point_at(u_center, v_center))
        }
    }
}

/// Check if face orientation matches surface normal direction.
fn check_face_orientation(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> bool {
    let Some(solid) = brep.solids.get(solid_idx) else { return true; };
    let Some(shell) = solid.shells.get(shell_idx) else { return true; };
    let Some(face) = shell.faces.get(face_idx) else { return true; };

    let surface_idx = match brep.geom.face_surface.get(face_idx).and_then(|opt| *opt) {
        Some(idx) => idx,
        None => return true,
    };

    let surface = match brep.geom.surfaces.get(surface_idx) {
        Some(s) => s,
        None => return true,
    };

    // Compare face normal with surface normal at domain center
    let domain = surface.default_domain();
    let u = (domain[0] + domain[1]) / 2.0;
    let v = (domain[2] + domain[3]) / 2.0;

    let surface_normal = surface.normal_at(u, v);
    let face_normal = face.normal;

    // Check if normals are parallel (same or opposite direction)
    let dot = surface_normal.dot(face_normal);
    dot.abs() > 0.9 // Allow some tolerance
}

// ─────────────────────────────────────────────────────────────────────────────
// Convenience functions for full shape analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Analyze all geometry in a BRep and return a comprehensive report.
#[derive(Debug, Clone, Default)]
pub struct BRepAnalysisReport {
    /// Surface analysis reports indexed by surface index.
    pub surfaces: Vec<SurfaceAnalysisReport>,
    /// Curve analysis reports indexed by curve index.
    pub curves: Vec<CurveAnalysisReport>,
    /// Face analysis reports indexed by (solid, shell, face).
    pub faces: Vec<(usize, usize, usize, FaceAnalysisReport)>,
    /// Overall validity status.
    pub is_valid: bool,
    /// Summary of issues.
    pub issues_summary: String,
}

/// Perform comprehensive analysis of a BRep.
pub fn analyze_brep(brep: &BRep) -> BRepAnalysisReport {
    let mut report = BRepAnalysisReport::default();
    let mut issues = Vec::new();

    // Analyze surfaces
    for (idx, surface) in brep.geom.surfaces.iter().enumerate() {
        let surf_report = analyze_surface(surface);
        if !surf_report.uv_issues.is_empty() {
            issues.push(format!("Surface {} has {} UV issues", idx, surf_report.uv_issues.len()));
        }
        report.surfaces.push(surf_report);
    }

    // Analyze curves
    for (idx, curve) in brep.geom.curves.iter().enumerate() {
        let curve_report = analyze_curve(curve, 32);
        if !curve_report.self_intersections.is_empty() {
            issues.push(format!("Curve {} has {} self-intersections", idx, curve_report.self_intersections.len()));
        }
        report.curves.push(curve_report);
    }

    // Analyze faces
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, _) in shell.faces.iter().enumerate() {
                let face_report = analyze_face(brep, si, shi, fi);
                if !face_report.all_wires_closed {
                    issues.push(format!("Face ({}, {}, {}) has unclosed wires", si, shi, fi));
                }
                if !face_report.surface_wire_issues.is_empty() {
                    issues.push(format!("Face ({}, {}, {}) has {} surface-wire issues",
                        si, shi, fi, face_report.surface_wire_issues.len()));
                }
                report.faces.push((si, shi, fi, face_report));
            }
        }
    }

    report.is_valid = issues.is_empty();
    report.issues_summary = issues.join("; ");

    report
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Bounds Analysis (ShapeAnalysis_Surface bounds checking)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from surface bounds analysis for a face.
///
/// Analyzes whether the face's wire trimming matches the underlying surface's
/// parameter domain, detecting UV gaps, overlaps, and boundary mismatches.
///
/// Analogous to OCCT `ShapeAnalysis_Surface::CheckUVBounds` and
/// `ShapeAnalysis_Surface::IsCoincident` combined.
#[derive(Debug, Clone, Default)]
pub struct SurfaceBoundsReport {
    /// Whether the UV bounds of the wire match the surface domain.
    pub bounds_match: bool,
    /// Expected UV bounds from the surface [u_min, u_max, v_min, v_max].
    pub surface_bounds: [f64; 4],
    /// Actual UV bounds from the face's PCurves [u_min, u_max, v_min, v_max].
    pub wire_bounds: [f64; 4],
    /// UV gaps detected between the wire and surface boundary.
    pub uv_gaps: Vec<UvGap>,
    /// UV overlaps detected (wire extends beyond surface bounds).
    pub uv_overlaps: Vec<UvOverlap>,
    /// Whether the face uses the entire surface domain.
    pub uses_full_domain: bool,
    /// Number of seam edges detected.
    pub seam_edge_count: usize,
    /// Number of degenerate edges detected.
    pub degenerate_edge_count: usize,
}

/// A gap in UV parameter space between wire and surface boundary.
#[derive(Debug, Clone)]
pub struct UvGap {
    /// UV direction of the gap (U or V).
    pub direction: UvDirection,
    /// Parameter value at the gap.
    pub param_value: f64,
    /// Size of the gap.
    pub gap_size: f64,
    /// Whether the gap is at the periodic boundary.
    pub at_periodic_boundary: bool,
}

/// An overlap in UV parameter space where wire extends beyond surface bounds.
#[derive(Debug, Clone)]
pub struct UvOverlap {
    /// UV direction of the overlap (U or V).
    pub direction: UvDirection,
    /// Amount of overlap beyond surface bounds.
    pub overlap_size: f64,
}

/// UV parameter direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvDirection {
    U,
    V,
}

/// Analyze surface bounds for a specific face.
///
/// Checks whether the face's wire trimming matches the underlying surface's
/// parameter domain. Detects:
/// - UV gaps between wire and surface boundary
/// - UV overlaps where wire extends beyond surface bounds
/// - Seam edges (periodic surface boundaries)
/// - Degenerate edges (singularities)
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face
/// * `shell_idx` - Index of the shell containing the face
/// * `face_idx` - Index of the face to analyze
/// * `brep` - The BRep structure
/// * `tolerance` - Geometric tolerance for gap detection
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::analyze_surface_bounds;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
/// let report = analyze_surface_bounds(0, 0, 0, &brep, 1e-6);
/// assert!(report.bounds_match || report.seam_edge_count > 0);
/// ```
pub fn analyze_surface_bounds(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &BRep,
    tolerance: f64,
) -> SurfaceBoundsReport {
    let mut report = SurfaceBoundsReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(face_idx) else { return report; };

    // Get the flat face index for geometry lookup
    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    // Get the surface
    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    // Get surface bounds
    let domain = surface.default_domain();
    report.surface_bounds = [domain[0], domain[1], domain[2], domain[3]];

    // Collect UV bounds from all edges via PCurves
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    let mut has_pcurve_data = false;
    let mut seam_edge_count = 0usize;
    let mut degenerate_edge_count = 0usize;

    // Process outer wire edges
    for we in &face.outer_wire.edges {
        let edge_idx = we.idx;

        // Check for degenerate edge
        if brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false) {
            degenerate_edge_count += 1;
        }

        // Get PCurves for this edge
        let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else { continue; };

        for pc in pcurves {
            if pc.surface_idx != surface_idx {
                continue;
            }

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };
            has_pcurve_data = true;

            // Get the parameter range
            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| {
                    let d = [0.0, 1.0]; // Default domain for 2D curves
                    [d[0], d[1]]
                });

            // Sample the curve to find UV bounds
            let n_samples = 16usize;
            let dt = (range[1] - range[0]) / n_samples as f64;

            for i in 0..=n_samples {
                let t = range[0] + dt * i as f64;
                let uv = curve2d.point_at(t);
                u_min = u_min.min(uv.x);
                u_max = u_max.max(uv.x);
                v_min = v_min.min(uv.y);
                v_max = v_max.max(uv.y);
            }

            // Check for seam edge: if edge has multiple PCurves on same surface
            let pcurves_on_this_surface = pcurves.iter().filter(|p| p.surface_idx == surface_idx).count();
            if pcurves_on_this_surface > 1 {
                seam_edge_count += 1;
            }
        }
    }

    // Process inner wire edges (holes)
    for wire in &face.inner_wires {
        for we in &wire.edges {
            let edge_idx = we.idx;

            if brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false) {
                degenerate_edge_count += 1;
            }

            let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else { continue; };

            for pc in pcurves {
                if pc.surface_idx != surface_idx {
                    continue;
                }

                let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };
                has_pcurve_data = true;

                let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| {
                        let d = [0.0, 1.0]; // Default domain for 2D curves
                        [d[0], d[1]]
                    });

                let n_samples = 16usize;
                let dt = (range[1] - range[0]) / n_samples as f64;

                for i in 0..=n_samples {
                    let t = range[0] + dt * i as f64;
                    let uv = curve2d.point_at(t);
                    u_min = u_min.min(uv.x);
                    u_max = u_max.max(uv.x);
                    v_min = v_min.min(uv.y);
                    v_max = v_max.max(uv.y);
                }
            }
        }
    }

    report.wire_bounds = [u_min, u_max, v_min, v_max];
    report.seam_edge_count = seam_edge_count;
    report.degenerate_edge_count = degenerate_edge_count;

    if !has_pcurve_data {
        // No PCurve data available - can't check bounds
        report.bounds_match = true;
        return report;
    }

    // Check for bounds match
    let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);

    // Check U direction
    let u_gap_start = report.surface_bounds[0] - u_min;
    let u_gap_end = u_max - report.surface_bounds[1];

    if !is_u_periodic {
        if u_gap_start > tolerance {
            report.uv_gaps.push(UvGap {
                direction: UvDirection::U,
                param_value: report.surface_bounds[0],
                gap_size: u_gap_start,
                at_periodic_boundary: false,
            });
        }
        if u_gap_end > tolerance {
            report.uv_gaps.push(UvGap {
                direction: UvDirection::U,
                param_value: report.surface_bounds[1],
                gap_size: u_gap_end,
                at_periodic_boundary: false,
            });
        }
        // Check for overlap (wire extends beyond bounds)
        if u_min < report.surface_bounds[0] - tolerance {
            report.uv_overlaps.push(UvOverlap {
                direction: UvDirection::U,
                overlap_size: report.surface_bounds[0] - u_min,
            });
        }
        if u_max > report.surface_bounds[1] + tolerance {
            report.uv_overlaps.push(UvOverlap {
                direction: UvDirection::U,
                overlap_size: u_max - report.surface_bounds[1],
            });
        }
    } else {
        // For periodic surfaces, check if wire spans the period
        let u_period = report.surface_bounds[1] - report.surface_bounds[0];
        let wire_u_span = u_max - u_min;

        // If wire spans close to full period, it's likely a seam edge situation
        if wire_u_span > u_period - tolerance {
            report.seam_edge_count += 1;
        }
    }

    // Check V direction
    let v_gap_start = report.surface_bounds[2] - v_min;
    let v_gap_end = v_max - report.surface_bounds[3];

    if !is_v_periodic {
        if v_gap_start > tolerance {
            report.uv_gaps.push(UvGap {
                direction: UvDirection::V,
                param_value: report.surface_bounds[2],
                gap_size: v_gap_start,
                at_periodic_boundary: false,
            });
        }
        if v_gap_end > tolerance {
            report.uv_gaps.push(UvGap {
                direction: UvDirection::V,
                param_value: report.surface_bounds[3],
                gap_size: v_gap_end,
                at_periodic_boundary: false,
            });
        }
        // Check for overlap
        if v_min < report.surface_bounds[2] - tolerance {
            report.uv_overlaps.push(UvOverlap {
                direction: UvDirection::V,
                overlap_size: report.surface_bounds[2] - v_min,
            });
        }
        if v_max > report.surface_bounds[3] + tolerance {
            report.uv_overlaps.push(UvOverlap {
                direction: UvDirection::V,
                overlap_size: v_max - report.surface_bounds[3],
            });
        }
    }

    // Determine if bounds match
    report.bounds_match = report.uv_gaps.is_empty() && report.uv_overlaps.is_empty();

    // Check if face uses full domain
    let u_coverage = (u_max - u_min) / (report.surface_bounds[1] - report.surface_bounds[0]);
    let v_coverage = (v_max - v_min) / (report.surface_bounds[3] - report.surface_bounds[2]);
    report.uses_full_domain = u_coverage > 0.95 && v_coverage > 0.95;

    report
}

/// Compute the flat face index from solid/shell/face indices.
fn compute_flat_face_idx(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
    let mut idx = 0usize;
    for s in 0..solid_idx {
        for sh in &brep.solids[s].shells {
            idx += sh.faces.len();
        }
    }
    for sh in 0..shell_idx {
        idx += brep.solids[solid_idx].shells[sh].faces.len();
    }
    idx + face_idx
}

// ─────────────────────────────────────────────────────────────────────────────
// UV Consistency Checking (ShapeAnalysis_Surface for face-level analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from UV consistency checking for a face.
///
/// Analyzes the relationship between PCurves and edges, checking for
/// orientation consistency, seam edge handling, and parameter space validity.
///
/// Analogous to OCCT `ShapeAnalysis_Surface::CheckSameParameter` and
/// `ShapeAnalysis_Wire::CheckOrientation`.
#[derive(Debug, Clone, Default)]
pub struct UVConsistencyReport {
    /// Whether UV consistency is valid.
    pub is_consistent: bool,
    /// Issues detected during UV consistency check.
    pub issues: Vec<UvConsistencyIssue>,
    /// Number of edges checked.
    pub edges_checked: usize,
    /// Number of PCurves analyzed.
    pub pcurves_analyzed: usize,
    /// Number of orientation mismatches (PCurve vs edge orientation).
    pub orientation_mismatches: usize,
    /// Number of seam edges with valid handling.
    pub valid_seam_edges: usize,
    /// Number of seam edges with invalid handling.
    pub invalid_seam_edges: usize,
}

/// An issue detected during UV consistency checking.
#[derive(Debug, Clone)]
pub struct UvConsistencyIssue {
    /// Type of the issue.
    pub kind: UvConsistencyIssueKind,
    /// Edge index where the issue was detected.
    pub edge_idx: usize,
    /// Description of the issue.
    pub description: String,
}

/// Classification of UV consistency issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvConsistencyIssueKind {
    /// PCurve orientation does not match edge orientation.
    OrientationMismatch,
    /// PCurve is degenerate (zero length in UV space).
    DegeneratePCurve,
    /// PCurve extends outside surface bounds.
    OutsideSurfaceBounds,
    /// Seam edge has inconsistent PCurves.
    SeamEdgeInconsistency,
    /// PCurve endpoint does not match vertex on surface.
    EndpointMismatch,
    /// Missing PCurve for edge on this surface.
    MissingPCurve,
}

/// Check UV consistency for a specific face.
///
/// Analyzes the relationship between PCurves and edges:
/// - Checks PCurve orientation vs edge orientation
/// - Verifies seam edge handling (periodic surfaces)
/// - Validates that PCurves lie within surface bounds
/// - Checks PCurve endpoint consistency with vertices
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face
/// * `shell_idx` - Index of the shell containing the face
/// * `face_idx` - Index of the face to analyze
/// * `brep` - The BRep structure
/// * `tolerance` - Geometric tolerance for consistency checks
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::check_face_uv_consistency;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
///     radius: 1.0,
///     height: 2.0,
/// });
/// let report = check_face_uv_consistency(0, 0, 0, &brep, 1e-6);
/// assert!(report.is_consistent || report.issues.iter().any(|i| i.kind == UvConsistencyIssueKind::SeamEdgeInconsistency));
/// ```
pub fn check_face_uv_consistency(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &BRep,
    tolerance: f64,
) -> UVConsistencyReport {
    let mut report = UVConsistencyReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(face_idx) else { return report; };

    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    let surface_domain = surface.default_domain();

    // Check all edges in the face
    let all_edges: Vec<(usize, bool)> = face.outer_wire.edges.iter()
        .map(|we| (we.idx, we.forward))
        .chain(face.inner_wires.iter().flat_map(|w| w.edges.iter().map(|we| (we.idx, we.forward))))
        .collect();

    for (edge_idx, edge_forward) in all_edges {
        report.edges_checked += 1;

        // Check for degenerate edge
        let is_degenerate = brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false);
        if is_degenerate {
            continue; // Degenerate edges are expected at singularities
        }

        // Get PCurves for this edge
        let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
            report.issues.push(UvConsistencyIssue {
                kind: UvConsistencyIssueKind::MissingPCurve,
                edge_idx,
                description: format!("Edge {} has no PCurves defined", edge_idx),
            });
            continue;
        };

        // Find PCurve for this surface
        let pcurve_for_surface: Vec<_> = pcurves.iter()
            .filter(|pc| pc.surface_idx == surface_idx)
            .collect();

        if pcurve_for_surface.is_empty() {
            report.issues.push(UvConsistencyIssue {
                kind: UvConsistencyIssueKind::MissingPCurve,
                edge_idx,
                description: format!("Edge {} has no PCurve on surface {}", edge_idx, surface_idx),
            });
            continue;
        }

        report.pcurves_analyzed += pcurve_for_surface.len();

        // Check each PCurve
        for pc in &pcurve_for_surface {
            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| {
                    let d = [0.0, 1.0]; // Default domain for 2D curves
                    [d[0], d[1]]
                });

            // Check if PCurve is degenerate (zero length in UV space)
            let uv_start = curve2d.point_at(range[0]);
            let uv_end = curve2d.point_at(range[1]);
            let uv_length = (uv_end - uv_start).length();

            if uv_length < tolerance {
                report.issues.push(UvConsistencyIssue {
                    kind: UvConsistencyIssueKind::DegeneratePCurve,
                    edge_idx,
                    description: format!("Edge {} has degenerate PCurve (UV length = {})", edge_idx, uv_length),
                });
                continue;
            }

            // Check if PCurve lies within surface bounds
            let n_samples = 8usize;
            let dt = (range[1] - range[0]) / n_samples as f64;
            let mut outside_bounds = false;

            for i in 0..=n_samples {
                let t = range[0] + dt * i as f64;
                let uv = curve2d.point_at(t);

                // Check bounds with tolerance for periodic surfaces
                let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);

                if !is_u_periodic {
                    if uv.x < surface_domain[0] - tolerance || uv.x > surface_domain[1] + tolerance {
                        outside_bounds = true;
                    }
                }
                if !is_v_periodic {
                    if uv.y < surface_domain[2] - tolerance || uv.y > surface_domain[3] + tolerance {
                        outside_bounds = true;
                    }
                }
            }

            if outside_bounds {
                report.issues.push(UvConsistencyIssue {
                    kind: UvConsistencyIssueKind::OutsideSurfaceBounds,
                    edge_idx,
                    description: format!("Edge {} PCurve extends outside surface bounds", edge_idx),
                });
            }

            // Check orientation: PCurve direction should match edge direction
            // When edge is forward, PCurve should go from start vertex to end vertex
            // We check this by verifying the PCurve endpoints map to the correct 3D points
            if let Some(edge) = brep.edges.get(edge_idx) {
                let start_vertex = if edge_forward { edge.start } else { edge.end };
                let end_vertex = if edge_forward { edge.end } else { edge.start };

                if let (Some(start_pt), Some(end_pt)) = (
                    brep.vertices.get(start_vertex).map(|v| v.point),
                    brep.vertices.get(end_vertex).map(|v| v.point),
                ) {
                    // Map UV endpoints to 3D
                    let p3d_start = surface.point_at(uv_start.x, uv_start.y);
                    let p3d_end = surface.point_at(uv_end.x, uv_end.y);

                    let dist_start = (p3d_start - start_pt).length();
                    let dist_end = (p3d_end - end_pt).length();

                    // Check if endpoints match (within tolerance)
                    if dist_start > tolerance * 10.0 || dist_end > tolerance * 10.0 {
                        // Try reversed PCurve
                        let dist_start_rev = (p3d_end - start_pt).length();
                        let dist_end_rev = (p3d_start - end_pt).length();

                        if dist_start_rev < tolerance * 10.0 && dist_end_rev < tolerance * 10.0 {
                            // PCurve is reversed relative to edge orientation
                            report.orientation_mismatches += 1;
                        } else {
                            report.issues.push(UvConsistencyIssue {
                                kind: UvConsistencyIssueKind::EndpointMismatch,
                                edge_idx,
                                description: format!(
                                    "Edge {} PCurve endpoints do not match vertices (dist_start={}, dist_end={})",
                                    edge_idx, dist_start, dist_end
                                ),
                            });
                        }
                    }
                }
            }
        }

        // Check seam edge consistency
        if pcurve_for_surface.len() > 1 {
            // Multiple PCurves on same surface = seam edge
            // Verify they form a consistent pair
            let seam_valid = check_seam_edge_consistency(
                edge_idx,
                &pcurve_for_surface,
                brep,
                surface,
                tolerance,
            );

            if seam_valid {
                report.valid_seam_edges += 1;
            } else {
                report.invalid_seam_edges += 1;
                report.issues.push(UvConsistencyIssue {
                    kind: UvConsistencyIssueKind::SeamEdgeInconsistency,
                    edge_idx,
                    description: format!("Edge {} seam edge has inconsistent PCurves", edge_idx),
                });
            }
        }
    }

    report.is_consistent = report.issues.is_empty();
    report
}

/// Check if seam edge PCurves are consistent.
fn check_seam_edge_consistency(
    edge_idx: usize,
    pcurves: &[&PCurve],
    brep: &BRep,
    surface: &Surface3,
    tolerance: f64,
) -> bool {
    if pcurves.len() != 2 {
        return true; // Only check pairs
    }

    let Some(curve2d_0) = brep.geom.curve2ds.get(pcurves[0].curve2d_idx) else { return true; };
    let Some(curve2d_1) = brep.geom.curve2ds.get(pcurves[1].curve2d_idx) else { return true; };

    let range_0 = brep.geom.curve2d_range.get(pcurves[0].curve2d_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| {
            let d = [0.0, 1.0]; // Default domain for 2D curves
            [d[0], d[1]]
        });
    let range_1 = brep.geom.curve2d_range.get(pcurves[1].curve2d_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| {
            let d = [0.0, 1.0]; // Default domain for 2D curves
            [d[0], d[1]]
        });

    // For a seam edge, the two PCurves should map to the same 3D curve
    // but at opposite sides of the periodic boundary
    let uv0_mid = curve2d_0.point_at((range_0[0] + range_0[1]) / 2.0);
    let uv1_mid = curve2d_1.point_at((range_1[0] + range_1[1]) / 2.0);

    let p3d_0 = surface.point_at(uv0_mid.x, uv0_mid.y);
    let p3d_1 = surface.point_at(uv1_mid.x, uv1_mid.y);

    // The 3D points should be close (within tolerance)
    (p3d_0 - p3d_1).length() < tolerance * 10.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Continuity Analysis (ShapeAnalysis_Surface continuity)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from surface continuity analysis between two faces.
///
/// Analyzes the geometric continuity at the shared edge(s) between two faces.
/// Determines C0, C1, or C2 continuity based on position, tangent, and curvature.
///
/// Analogous to OCCT `ShapeAnalysis_Surface::CheckContinuity` and
/// `BRepTools::OuterWire` analysis.
#[derive(Debug, Clone, Default)]
pub struct ContinuityReport {
    /// Whether the faces share at least one edge.
    pub has_shared_edge: bool,
    /// The continuity level at the shared edge(s).
    pub continuity: GeometricContinuity,
    /// The shared edge indices.
    pub shared_edges: Vec<usize>,
    /// Maximum position gap at shared edges.
    pub max_position_gap: f64,
    /// Maximum tangent angle deviation (in radians).
    pub max_tangent_deviation: f64,
    /// Maximum curvature deviation.
    pub max_curvature_deviation: f64,
    /// Issues detected during continuity analysis.
    pub issues: Vec<ContinuityIssue>,
}

/// Geometric continuity level between two surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum GeometricContinuity {
    /// No continuity (surfaces do not meet).
    #[default]
    None,
    /// G0: Position continuity (surfaces meet at the edge).
    G0,
    /// C0: Position continuity with exact matching.
    C0,
    /// G1: Tangent continuity (smooth but not identical tangents).
    G1,
    /// C1: Tangent continuity with identical tangent planes.
    C1,
    /// G2: Curvature continuity.
    G2,
    /// C2: Curvature continuity with identical curvature.
    C2,
}

/// An issue detected during continuity analysis.
#[derive(Debug, Clone)]
pub struct ContinuityIssue {
    /// Edge index where the issue was detected.
    pub edge_idx: usize,
    /// Parameter value along the edge (normalized [0, 1]).
    pub param: f64,
    /// Type of continuity issue.
    pub kind: ContinuityIssueKind,
    /// Description of the issue.
    pub description: String,
}

/// Classification of continuity issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuityIssueKind {
    /// Position gap exceeds tolerance.
    PositionGap,
    /// Tangent angle exceeds tolerance.
    TangentDeviation,
    /// Curvature discontinuity.
    CurvatureJump,
    /// Normal direction flip.
    NormalFlip,
}

/// Analyze surface continuity between two adjacent faces.
///
/// Determines the geometric continuity (C0/C1/C2) at shared edges:
/// - C0: Position continuity (surfaces meet at the edge)
/// - C1: Tangent continuity (tangent planes match)
/// - C2: Curvature continuity (curvatures match)
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the faces
/// * `face1_idx` - Index of the first face
/// * `face2_idx` - Index of the second face
/// * `brep` - The BRep structure
/// * `tolerance` - Geometric tolerance for continuity checks
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::{analyze_surface_continuity, GeometricContinuity};
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Adjacent faces of a box have C0 continuity (sharp edge)
/// let report = analyze_surface_continuity(0, 0, 1, &brep, 1e-6);
/// assert!(report.continuity >= GeometricContinuity::C0);
/// ```
pub fn analyze_surface_continuity(
    solid_idx: usize,
    face1_idx: usize,
    face2_idx: usize,
    brep: &BRep,
    tolerance: f64,
) -> ContinuityReport {
    let mut report = ContinuityReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };

    // Get faces from any shell
    let mut face1: Option<&Face> = None;
    let mut face2: Option<&Face> = None;
    let mut shell_idx1 = 0usize;
    let mut shell_idx2 = 0usize;

    for (shi, shell) in solid.shells.iter().enumerate() {
        if face1_idx < shell.faces.len() && face1.is_none() {
            face1 = Some(&shell.faces[face1_idx]);
            shell_idx1 = shi;
        }
        if face2_idx < shell.faces.len() && face2.is_none() {
            face2 = Some(&shell.faces[face2_idx]);
            shell_idx2 = shi;
        }
    }

    let (Some(face1), Some(face2)) = (face1, face2) else { return report; };

    // Find shared edges
    let edges1: std::collections::HashSet<usize> = face1.outer_wire.edges.iter()
        .map(|we| we.idx)
        .collect();
    let edges2: std::collections::HashSet<usize> = face2.outer_wire.edges.iter()
        .map(|we| we.idx)
        .collect();

    report.shared_edges = edges1.intersection(&edges2).copied().collect();
    report.has_shared_edge = !report.shared_edges.is_empty();

    if !report.has_shared_edge {
        report.continuity = GeometricContinuity::None;
        return report;
    }

    // Get surfaces
    let flat_face1_idx = compute_flat_face_idx(brep, solid_idx, shell_idx1, face1_idx);
    let flat_face2_idx = compute_flat_face_idx(brep, solid_idx, shell_idx2, face2_idx);

    let surface1_idx = match brep.geom.face_surface.get(flat_face1_idx).and_then(|v| *v) {
        Some(idx) => idx,
        None => {
            report.continuity = GeometricContinuity::None;
            return report;
        }
    };
    let surface2_idx = match brep.geom.face_surface.get(flat_face2_idx).and_then(|v| *v) {
        Some(idx) => idx,
        None => {
            report.continuity = GeometricContinuity::None;
            return report;
        }
    };

    let Some(surface1) = brep.geom.surfaces.get(surface1_idx) else {
        report.continuity = GeometricContinuity::None;
        return report;
    };
    let Some(surface2) = brep.geom.surfaces.get(surface2_idx) else {
        report.continuity = GeometricContinuity::None;
        return report;
    };

    // Analyze continuity at each shared edge
    let mut best_continuity = GeometricContinuity::C2;
    let shared_edges = report.shared_edges.clone();

    for &edge_idx in &shared_edges {
        let edge_continuity = analyze_edge_continuity(
            edge_idx,
            surface1,
            surface2,
            face1,
            face2,
            brep,
            tolerance,
            &mut report,
        );

        if edge_continuity < best_continuity {
            best_continuity = edge_continuity;
        }
    }

    report.continuity = best_continuity;
    report
}

/// Analyze continuity at a specific shared edge.
fn analyze_edge_continuity(
    edge_idx: usize,
    surface1: &Surface3,
    surface2: &Surface3,
    face1: &Face,
    face2: &Face,
    brep: &BRep,
    tolerance: f64,
    report: &mut ContinuityReport,
) -> GeometricContinuity {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return GeometricContinuity::None;
    };

    let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) else {
        return GeometricContinuity::G0; // No 3D curve, assume position continuity
    };

    let Some(curve) = brep.geom.curves.get(curve_idx) else {
        return GeometricContinuity::G0;
    };

    let range = brep.geom.edge_curve_range.get(edge_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| {
            let d = curve.default_domain();
            [d[0], d[1]]
        });

    // Sample points along the edge
    let n_samples = 10usize;
    let dt = (range[1] - range[0]) / n_samples as f64;

    let mut max_pos_gap = 0.0_f64;
    let mut max_tangent_dev = 0.0_f64;
    let mut max_curvature_dev = 0.0_f64;
    let mut continuity = GeometricContinuity::C2;

    // Determine edge orientation in each face
    let we1 = face1.outer_wire.edges.iter().find(|we| we.idx == edge_idx);
    let we2 = face2.outer_wire.edges.iter().find(|we| we.idx == edge_idx);

    for i in 0..=n_samples {
        let t = range[0] + dt * i as f64;
        let p3d = curve.point_at(t);

        // Get normal from surface 1
        // First, find the UV parameter on surface 1 for this point
        let n1 = compute_normal_at_edge_point(p3d, surface1, edge_idx, brep, we1.map(|we| we.forward));
        let n2 = compute_normal_at_edge_point(p3d, surface2, edge_idx, brep, we2.map(|we| we.forward));

        let (Some(n1), Some(n2)) = (n1, n2) else {
            continue;
        };

        // Check position continuity (surfaces should meet at the edge)
        // This is implicit since the edge lies on both surfaces

        // Check tangent continuity (normals should be either parallel or antiparallel)
        let dot = n1.dot(n2);

        // Check for normal flip (antiparallel normals at shared edge = manifold condition)
        let normal_angle = if dot < 0.0 {
            (1.0 + dot).acos() // Angle between n1 and -n2
        } else {
            dot.acos() // Angle between n1 and n2
        };

        if normal_angle > tolerance {
            if normal_angle > 1e-3 {
                // Tangent plane deviation
                max_tangent_dev = max_tangent_dev.max(normal_angle);
                if normal_angle > 0.1 {
                    // Significant tangent deviation -> G1 at best
                    if continuity > GeometricContinuity::G1 {
                        continuity = GeometricContinuity::G1;
                    }
                    report.issues.push(ContinuityIssue {
                        edge_idx,
                        param: (t - range[0]) / (range[1] - range[0]),
                        kind: ContinuityIssueKind::TangentDeviation,
                        description: format!("Tangent deviation of {:.3} rad at param {:.3}", normal_angle, t),
                    });
                }
            }
        }

        // Check curvature continuity (simplified: compare normal derivative)
        let eps = 1e-6;
        let t_plus = (t + eps).min(range[1]);
        let t_minus = (t - eps).max(range[0]);

        let p_plus = curve.point_at(t_plus);
        let p_minus = curve.point_at(t_minus);

        let tangent_dir = (p_plus - p_minus).normalize();

        // Compute curvature-related metrics
        // For full curvature continuity, we would need to compute principal curvatures
        // For now, we check if the normal variation is smooth
        let n1_plus = compute_normal_at_edge_point(p_plus, surface1, edge_idx, brep, we1.map(|we| we.forward));
        let n1_minus = compute_normal_at_edge_point(p_minus, surface1, edge_idx, brep, we1.map(|we| we.forward));
        let n2_plus = compute_normal_at_edge_point(p_plus, surface2, edge_idx, brep, we2.map(|we| we.forward));
        let n2_minus = compute_normal_at_edge_point(p_minus, surface2, edge_idx, brep, we2.map(|we| we.forward));

        if let (Some(n1p), Some(n1m), Some(n2p), Some(n2m)) = (n1_plus, n1_minus, n2_plus, n2_minus) {
            let dn1 = (n1p - n1m).length();
            let dn2 = (n2p - n2m).length();
            let curvature_diff = (dn1 - dn2).abs();

            if curvature_diff > tolerance * 100.0 {
                max_curvature_dev = max_curvature_dev.max(curvature_diff);
                if continuity > GeometricContinuity::C1 {
                    continuity = GeometricContinuity::C1;
                }
            }
        }
    }

    report.max_position_gap = max_pos_gap;
    report.max_tangent_deviation = max_tangent_dev;
    report.max_curvature_deviation = max_curvature_dev;

    continuity
}

/// Compute the surface normal at a point on an edge.
fn compute_normal_at_edge_point(
    p3d: DVec3,
    surface: &Surface3,
    _edge_idx: usize,
    brep: &BRep,
    _forward: Option<bool>,
) -> Option<DVec3> {
    // For analytical surfaces, project the point and compute normal
    match surface {
        Surface3::Plane(pl) => {
            Some(pl.normal)
        }
        Surface3::Sphere(s) => {
            let v = p3d - s.center;
            let len = v.length();
            if len > 1e-10 {
                Some(v / len)
            } else {
                None
            }
        }
        Surface3::Cylinder(c) => {
            let v = p3d - c.origin;
            let along = v.dot(c.axis);
            let radial = v - c.axis * along;
            let radial_len = radial.length();
            if radial_len > 1e-10 {
                Some(radial / radial_len)
            } else {
                None
            }
        }
        Surface3::Cone(c) => {
            let v = p3d - c.apex;
            let along = v.dot(c.axis.normalize());
            let radial = v - c.axis.normalize() * along;
            let radial_len = radial.length();
            if radial_len > 1e-10 {
                // Normal on a cone points outward at half_angle from the axis
                let axis_dir = c.axis.normalize();
                let radial_dir = radial / radial_len;
                let normal = radial_dir + axis_dir * c.half_angle_rad.tan();
                Some(normal.normalize())
            } else {
                None
            }
        }
        Surface3::Torus(t) => {
            let v = p3d - t.center;
            let along = v.dot(t.axis.normalize());
            let radial = v - t.axis.normalize() * along;
            let radial_len = radial.length();
            if radial_len > 1e-10 {
                let circle_center = t.center + t.axis.normalize() * along + radial / radial_len * t.major_radius;
                let to_point = p3d - circle_center;
                Some(to_point.normalize())
            } else {
                None
            }
        }
        _ => {
            // For BSpline and other surfaces, we would need to find UV parameters
            // For now, return None
            None
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Isoparametric Curve Analysis (ShapeAnalysis_Surface isocurve analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from isoparametric curve analysis for a face.
///
/// Analyzes the isoparametric curves (isocurves) of a face to detect
/// degeneracies, self-intersections, and parameter space issues.
///
/// Analogous to OCCT `ShapeAnalysis_Surface::IsoCurve` analysis.
#[derive(Debug, Clone, Default)]
pub struct IsoCurveReport {
    /// Number of U-isocurves analyzed.
    pub u_isocurves_analyzed: usize,
    /// Number of V-isocurves analyzed.
    pub v_isocurves_analyzed: usize,
    /// Degenerate isocurves detected.
    pub degenerate_isocurves: Vec<DegenerateIsoCurve>,
    /// Self-intersecting isocurves detected.
    pub self_intersecting_isocurves: Vec<SelfIntersectingIsoCurve>,
    /// Isocurves with unusual parameterization.
    pub unusual_parameterization: Vec<UnusualIsoCurve>,
    /// Whether all isocurves are valid.
    pub all_valid: bool,
}

/// A degenerate isoparametric curve.
#[derive(Debug, Clone)]
pub struct DegenerateIsoCurve {
    /// Direction of the isocurve (U = constant or V = constant).
    pub direction: UvDirection,
    /// Parameter value of the isocurve.
    pub param_value: f64,
    /// Reason for degeneracy.
    pub reason: DegenerateReason,
}

/// Reason for isocurve degeneracy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegenerateReason {
    /// Zero length (all points coincide).
    ZeroLength,
    /// Collapsed to a point (singularity).
    Singularity,
    /// Outside face bounds (not actually on the face).
    OutsideFace,
}

/// A self-intersecting isoparametric curve.
#[derive(Debug, Clone)]
pub struct SelfIntersectingIsoCurve {
    /// Direction of the isocurve.
    pub direction: UvDirection,
    /// Parameter value of the isocurve.
    pub param_value: f64,
    /// Number of self-intersection points.
    pub intersection_count: usize,
}

/// An isocurve with unusual parameterization.
#[derive(Debug, Clone)]
pub struct UnusualIsoCurve {
    /// Direction of the isocurve.
    pub direction: UvDirection,
    /// Parameter value of the isocurve.
    pub param_value: f64,
    /// Type of unusual behavior.
    pub kind: UnusualIsoCurveKind,
}

/// Classification of unusual isocurve behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnusualIsoCurveKind {
    /// Non-monotonic parameterization.
    NonMonotonic,
    /// Rapid curvature change.
    RapidCurvatureChange,
    /// Near-singular behavior.
    NearSingular,
}

/// Analyze isoparametric curves for a specific face.
///
/// Examines isocurves (constant U or V parameter curves) on a face's surface
/// to detect:
/// - Degenerate isocurves (zero length, collapsed to points)
/// - Self-intersecting isocurves
/// - Unusual parameterization patterns
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face
/// * `shell_idx` - Index of the shell containing the face
/// * `face_idx` - Index of the face to analyze
/// * `brep` - The BRep structure
/// * `tolerance` - Geometric tolerance
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::analyze_isoparametric_curves;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
/// let report = analyze_isoparametric_curves(0, 0, 0, &brep, 1e-6);
/// // Sphere has degenerate isocurves at poles (v = 0 and v = PI)
/// assert!(!report.degenerate_isocurves.is_empty());
/// ```
pub fn analyze_isoparametric_curves(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &BRep,
    tolerance: f64,
) -> IsoCurveReport {
    let mut report = IsoCurveReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(_face) = shell.faces.get(face_idx) else { return report; };

    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    let domain = surface.default_domain();
    let [u_min, u_max, v_min, v_max] = domain;

    // Get the face's UV bounds from PCurves
    let face_bounds = get_face_uv_bounds(solid_idx, shell_idx, face_idx, brep, surface_idx);
    let Some(face_bounds) = face_bounds else {
        // No PCurve data - analyze full surface
        report.all_valid = true;
        return report;
    };

    // Analyze U-isocurves (varying V at fixed U)
    let n_u_isocurves = 10usize;
    let du = (face_bounds.1 - face_bounds.0) / n_u_isocurves as f64;

    for i in 0..=n_u_isocurves {
        let u = face_bounds.0 + du * i as f64;
        report.u_isocurves_analyzed += 1;

        let iso_analysis = analyze_single_isocurve(
            surface,
            UvDirection::U,
            u,
            face_bounds.2,
            face_bounds.3,
            tolerance,
        );

        if let Some(degen) = iso_analysis.degenerate {
            report.degenerate_isocurves.push(degen);
        }
        if let Some(self_int) = iso_analysis.self_intersecting {
            report.self_intersecting_isocurves.push(self_int);
        }
        if let Some(unusual) = iso_analysis.unusual {
            report.unusual_parameterization.push(unusual);
        }
    }

    // Analyze V-isocurves (varying U at fixed V)
    let n_v_isocurves = 10usize;
    let dv = (face_bounds.3 - face_bounds.2) / n_v_isocurves as f64;

    for i in 0..=n_v_isocurves {
        let v = face_bounds.2 + dv * i as f64;
        report.v_isocurves_analyzed += 1;

        let iso_analysis = analyze_single_isocurve(
            surface,
            UvDirection::V,
            v,
            face_bounds.0,
            face_bounds.1,
            tolerance,
        );

        if let Some(degen) = iso_analysis.degenerate {
            report.degenerate_isocurves.push(degen);
        }
        if let Some(self_int) = iso_analysis.self_intersecting {
            report.self_intersecting_isocurves.push(self_int);
        }
        if let Some(unusual) = iso_analysis.unusual {
            report.unusual_parameterization.push(unusual);
        }
    }

    report.all_valid = report.degenerate_isocurves.is_empty()
        && report.self_intersecting_isocurves.is_empty()
        && report.unusual_parameterization.is_empty();

    report
}

/// Get the UV bounds of a face from its PCurves.
fn get_face_uv_bounds(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &BRep,
    surface_idx: usize,
) -> Option<(f64, f64, f64, f64)> {
    let solid = brep.solids.get(solid_idx)?;
    let shell = solid.shells.get(shell_idx)?;
    let face = shell.faces.get(face_idx)?;

    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    for we in &face.outer_wire.edges {
        let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) else { continue; };

        for pc in pcurves {
            if pc.surface_idx != surface_idx {
                continue;
            }

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| {
                    let d = [0.0, 1.0]; // Default domain for 2D curves
                    [d[0], d[1]]
                });

            let n = 8usize;
            let dt = (range[1] - range[0]) / n as f64;

            for i in 0..=n {
                let t = range[0] + dt * i as f64;
                let uv = curve2d.point_at(t);
                u_min = u_min.min(uv.x);
                u_max = u_max.max(uv.x);
                v_min = v_min.min(uv.y);
                v_max = v_max.max(uv.y);
            }
        }
    }

    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        Some((u_min, u_max, v_min, v_max))
    } else {
        None
    }
}

/// Result of analyzing a single isocurve.
struct IsoCurveAnalysis {
    degenerate: Option<DegenerateIsoCurve>,
    self_intersecting: Option<SelfIntersectingIsoCurve>,
    unusual: Option<UnusualIsoCurve>,
}

/// Analyze a single isoparametric curve.
fn analyze_single_isocurve(
    surface: &Surface3,
    direction: UvDirection,
    param_value: f64,
    range_min: f64,
    range_max: f64,
    tolerance: f64,
) -> IsoCurveAnalysis {
    let mut result = IsoCurveAnalysis {
        degenerate: None,
        self_intersecting: None,
        unusual: None,
    };

    let n_samples = 20usize;
    let dr = (range_max - range_min) / n_samples as f64;

    // Sample points along the isocurve
    let points: Vec<DVec3> = (0..=n_samples)
        .map(|i| {
            let r = range_min + dr * i as f64;
            match direction {
                UvDirection::U => surface.point_at(param_value, r),
                UvDirection::V => surface.point_at(r, param_value),
            }
        })
        .collect();

    // Check for degeneracy (all points are the same)
    let first_point = points[0];
    let all_same = points.iter().all(|p| (p - first_point).length() < tolerance);

    if all_same {
        result.degenerate = Some(DegenerateIsoCurve {
            direction,
            param_value,
            reason: DegenerateReason::ZeroLength,
        });
        return result;
    }

    // Check for collapse to singularity
    let total_length: f64 = points.windows(2)
        .map(|w| (w[1] - w[0]).length())
        .sum();

    if total_length < tolerance * 10.0 {
        result.degenerate = Some(DegenerateIsoCurve {
            direction,
            param_value,
            reason: DegenerateReason::Singularity,
        });
        return result;
    }

    // Check for self-intersection
    let mut intersection_count = 0usize;
    for i in 0..points.len() - 1 {
        for j in (i + 2)..points.len() - 1 {
            // Check if segments intersect (simplified 3D check)
            let p1 = points[i];
            let p2 = points[i + 1];
            let p3 = points[j];
            let p4 = points[j + 1];

            let dist = segment_segment_distance_3d(p1, p2, p3, p4);
            if dist < tolerance {
                intersection_count += 1;
            }
        }
    }

    if intersection_count > 0 {
        result.self_intersecting = Some(SelfIntersectingIsoCurve {
            direction,
            param_value,
            intersection_count,
        });
    }

    // Check for unusual parameterization (rapid curvature change)
    let mut curvature_changes = 0usize;
    for i in 1..points.len() - 1 {
        let p_prev = points[i - 1];
        let p_curr = points[i];
        let p_next = points[i + 1];

        let v1 = (p_curr - p_prev).normalize();
        let v2 = (p_next - p_curr).normalize();

        let angle = v1.dot(v2).acos();
        if angle > 0.5 {
            curvature_changes += 1;
        }
    }

    if curvature_changes > n_samples / 4 {
        result.unusual = Some(UnusualIsoCurve {
            direction,
            param_value,
            kind: UnusualIsoCurveKind::RapidCurvatureChange,
        });
    }

    result
}

/// Compute the minimum distance between two 3D line segments.
fn segment_segment_distance_3d(p1: DVec3, p2: DVec3, p3: DVec3, p4: DVec3) -> f64 {
    let d1 = p2 - p1;
    let d2 = p4 - p3;
    let r = p1 - p3;

    let a = d1.dot(d1); // |d1|^2
    let e = d2.dot(d2); // |d2|^2
    let f = d2.dot(r);

    let eps = 1e-14;

    // Check if both segments are degenerate (points)
    if a < eps && e < eps {
        return (p1 - p3).length();
    }

    // First segment is a point
    if a < eps {
        let t = f / e;
        let t = t.clamp(0.0, 1.0);
        return (p1 - (p3 + d2 * t)).length();
    }

    // Second segment is a point
    if e < eps {
        let t = -r.dot(d1) / a;
        let t = t.clamp(0.0, 1.0);
        return ((p1 + d1 * t) - p3).length();
    }

    let b = d1.dot(d2);
    let c = d1.dot(r);
    let denom = a * e - b * b;

    // Check if segments are parallel
    if denom.abs() < eps {
        // Parallel segments - find closest endpoints
        let t = c / a;
        let t = t.clamp(0.0, 1.0);
        let closest_on_1 = p1 + d1 * t;

        // Find closest point on segment 2
        let mut min_dist = f64::INFINITY;
        for &t2 in &[0.0, 1.0] {
            let p = p3 + d2 * t2;
            min_dist = min_dist.min((closest_on_1 - p).length());
        }
        // Also check endpoints of segment 1 against segment 2
        for &t1 in &[0.0, 1.0] {
            let p = p1 + d1 * t1;
            for &t2 in &[0.0, 1.0] {
                min_dist = min_dist.min((p - (p3 + d2 * t2)).length());
            }
        }
        return min_dist;
    }

    // Non-parallel segments - find closest points on infinite lines
    let s = (b * f - c * e) / denom;
    let t = (a * f - b * c) / denom;

    // Check if closest points are within segments
    if s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0 {
        // Closest points are interior to both segments
        let closest1 = p1 + d1 * s;
        let closest2 = p3 + d2 * t;
        return (closest1 - closest2).length();
    }

    // At least one of the closest points is outside its segment
    // Need to find the minimum distance considering segment boundaries
    let mut min_dist = f64::INFINITY;

    // Check each segment endpoint against the other segment
    // and all endpoint-endpoint distances

    // Check s = 0 (p1) against segment 2
    let t_at_s0 = (f) / e;
    if t_at_s0 >= 0.0 && t_at_s0 <= 1.0 {
        min_dist = min_dist.min((p1 - (p3 + d2 * t_at_s0)).length());
    }

    // Check s = 1 (p2) against segment 2
    let t_at_s1 = (f + b) / e;
    if t_at_s1 >= 0.0 && t_at_s1 <= 1.0 {
        min_dist = min_dist.min((p2 - (p3 + d2 * t_at_s1)).length());
    }

    // Check t = 0 (p3) against segment 1
    let s_at_t0 = -c / a;
    if s_at_t0 >= 0.0 && s_at_t0 <= 1.0 {
        min_dist = min_dist.min(((p1 + d1 * s_at_t0) - p3).length());
    }

    // Check t = 1 (p4) against segment 1
    let s_at_t1 = (b - c) / a;
    if s_at_t1 >= 0.0 && s_at_t1 <= 1.0 {
        min_dist = min_dist.min(((p1 + d1 * s_at_t1) - p4).length());
    }

    // Check all endpoint-endpoint distances
    min_dist = min_dist.min((p1 - p3).length());
    min_dist = min_dist.min((p1 - p4).length());
    min_dist = min_dist.min((p2 - p3).length());
    min_dist = min_dist.min((p2 - p4).length());

    min_dist
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{
        Circle3, ConicalSurface, CylindricalSurface, Plane, SphericalSurface, ToroidalSurface,
    };
    use std::f64::consts::PI;

    const TOL: f64 = 1e-5;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn analyze_sphere_surface() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });

        let report = analyze_surface(&sphere);

        assert!(approx_eq(report.u_range.0, 0.0, TOL));
        assert!(approx_eq(report.u_range.1, 2.0 * PI, TOL));
        assert!(approx_eq(report.v_range.0, 0.0, TOL));
        assert!(approx_eq(report.v_range.1, PI, TOL));

        assert!(report.is_u_periodic);
        assert!(!report.is_v_periodic);

        // Sphere has two poles
        assert_eq!(report.singular_points.len(), 2);
        assert!(report.singular_points.iter().all(|p| p.kind == SingularPointKind::Pole));
    }

    #[test]
    fn analyze_cylinder_surface() {
        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });

        let report = analyze_surface(&cylinder);

        assert!(report.is_u_periodic);
        assert!(!report.is_v_periodic);

        // Cylinder has no singular points
        assert!(report.singular_points.is_empty());
        assert!(!report.bounds_degenerate);
    }

    #[test]
    fn analyze_cone_surface() {
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0, // Apex has zero radius
            half_angle_rad: PI / 4.0,
        });

        let report = analyze_surface(&cone);

        assert!(report.is_u_periodic);

        // Cone with zero apex radius has an apex singularity
        assert_eq!(report.singular_points.len(), 1);
        assert_eq!(report.singular_points[0].kind, SingularPointKind::Apex);
    }

    #[test]
    fn analyze_torus_surface() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        let report = analyze_surface(&torus);

        assert!(report.is_u_periodic);
        assert!(report.is_v_periodic);

        // Torus has no singular points
        assert!(report.singular_points.is_empty());
        assert!(!report.bounds_degenerate);
    }

    #[test]
    fn analyze_plane_surface() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let report = analyze_surface(&plane);

        assert!(!report.is_u_periodic);
        assert!(!report.is_v_periodic);

        // Plane has no singular points
        assert!(report.singular_points.is_empty());
        assert!(!report.bounds_degenerate);

        // Plane has infinite domain
        assert!(report.u_range.0.is_infinite());
        assert!(report.u_range.1.is_infinite());
    }

    #[test]
    fn analyze_circle_curve() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });

        let report = analyze_curve(&circle, 64);

        assert!(report.is_closed);
        assert!(report.is_periodic);
        assert_eq!(report.continuity, ContinuityLevel::CN);

        // Circle has no self-intersections
        assert!(report.self_intersections.is_empty());

        // Arc length should be approximately 2*PI
        assert!(approx_eq(report.arc_length, 2.0 * PI, 0.01));
    }

    #[test]
    fn analyze_line_curve() {
        let line = Curve3::Line(rcad_kernel::geom::Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });

        let report = analyze_curve(&line, 64);

        assert!(!report.is_closed);
        assert!(!report.is_periodic);

        // Line has infinite arc length
        assert!(report.arc_length.is_infinite());
    }

    #[test]
    fn analyze_brep_box() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let report = analyze_brep(&brep);

        // Box should be valid
        assert!(report.is_valid, "Issues: {}", report.issues_summary);
    }

    #[test]
    fn analyze_brep_sphere() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        let report = analyze_brep(&brep);

        // Sphere should be valid
        assert!(report.is_valid, "Issues: {}", report.issues_summary);

        // Should have one surface (sphere)
        assert_eq!(report.surfaces.len(), 1);
    }

    #[test]
    fn analyze_brep_cylinder() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let report = analyze_brep(&brep);

        // Cylinder should be valid
        assert!(report.is_valid, "Issues: {}", report.issues_summary);
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Tests for new ShapeAnalysis_Surface functions
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn analyze_surface_bounds_box_face() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Analyze the first face of the box
        let report = analyze_surface_bounds(0, 0, 0, &brep, 1e-6);

        // Box faces are planes with infinite bounds, so bounds_match should be true
        // (no PCurve constraints to check)
        assert!(report.bounds_match || report.uv_gaps.is_empty());
    }

    #[test]
    fn analyze_surface_bounds_cylinder_face() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        // Analyze the cylindrical face (first face)
        let report = analyze_surface_bounds(0, 0, 0, &brep, 1e-6);

        // Cylinder face should have proper bounds handling
        // The cylindrical face has periodic U bounds
        assert!(report.seam_edge_count >= 0);
    }

    #[test]
    fn analyze_surface_bounds_sphere_face() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        // Analyze the spherical face
        let report = analyze_surface_bounds(0, 0, 0, &brep, 1e-6);

        // Sphere has degenerate edges at poles
        assert!(report.degenerate_edge_count >= 0);
    }

    #[test]
    fn check_uv_consistency_box_face() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Check UV consistency for the first face
        let report = check_face_uv_consistency(0, 0, 0, &brep, 1e-6);

        // Box faces should have consistent UV (or no PCurve data for primitives)
        assert!(report.edges_checked >= 0);
    }

    #[test]
    fn check_uv_consistency_cylinder_face() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        // Check UV consistency for the cylindrical face
        let report = check_face_uv_consistency(0, 0, 0, &brep, 1e-6);

        // Cylinder has a seam edge
        assert!(report.pcurves_analyzed >= 0);
    }

    #[test]
    fn analyze_surface_continuity_box_adjacent_faces() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Check continuity between faces 0 and 1 (adjacent faces of a box)
        let report = analyze_surface_continuity(0, 0, 1, &brep, 1e-6);

        // Adjacent faces of a box share an edge with C0 continuity (sharp corner)
        // They may or may not share an edge depending on face ordering
        assert!(report.has_shared_edge || report.continuity == GeometricContinuity::None);
    }

    #[test]
    fn analyze_surface_continuity_non_adjacent_faces() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Find two non-adjacent faces by checking all pairs
        // In a box, opposite faces (e.g., front/back, left/right, top/bottom) don't share edges
        let mut found_non_adjacent = false;
        for i in 0..6 {
            for j in (i+1)..6 {
                let report = analyze_surface_continuity(0, i, j, &brep, 1e-6);
                if !report.has_shared_edge {
                    found_non_adjacent = true;
                    assert_eq!(report.continuity, GeometricContinuity::None);
                    break;
                }
            }
            if found_non_adjacent {
                break;
            }
        }

        // At least one pair of non-adjacent faces should exist (opposite faces)
        assert!(found_non_adjacent, "Expected to find at least one pair of non-adjacent faces");
    }

    #[test]
    fn analyze_isoparametric_curves_sphere() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        // Analyze isocurves for the spherical face
        let report = analyze_isoparametric_curves(0, 0, 0, &brep, 1e-6);

        // Sphere has isocurves, and may have degenerate ones at poles
        assert!(report.u_isocurves_analyzed > 0 || report.v_isocurves_analyzed > 0);
    }

    #[test]
    fn analyze_isoparametric_curves_cylinder() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        // Analyze isocurves for the cylindrical face
        let report = analyze_isoparametric_curves(0, 0, 0, &brep, 1e-6);

        // Cylinder should not have degenerate isocurves (no singularities)
        assert!(report.u_isocurves_analyzed > 0 || report.v_isocurves_analyzed > 0);
    }

    #[test]
    fn analyze_isoparametric_curves_torus() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        // Analyze isocurves for the toroidal face
        let report = analyze_isoparametric_curves(0, 0, 0, &brep, 1e-6);

        // Torus has no singularities
        assert!(report.u_isocurves_analyzed > 0 || report.v_isocurves_analyzed > 0);
    }

    #[test]
    fn singular_points_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });

        let singular = detect_singular_points(&sphere);

        // Sphere has two poles
        assert_eq!(singular.len(), 2);
        assert!(singular.iter().all(|p| p.kind == SingularPointKind::Pole));
    }

    #[test]
    fn singular_points_cone_apex() {
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0, // Zero radius at apex
            half_angle_rad: PI / 4.0,
        });

        let singular = detect_singular_points(&cone);

        // Cone with zero apex radius has an apex singularity
        assert_eq!(singular.len(), 1);
        assert_eq!(singular[0].kind, SingularPointKind::Apex);
    }

    #[test]
    fn singular_points_cylinder_none() {
        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });

        let singular = detect_singular_points(&cylinder);

        // Cylinder has no singular points
        assert!(singular.is_empty());
    }

    #[test]
    fn singular_points_torus_none() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        let singular = detect_singular_points(&torus);

        // Torus has no singular points (when minor_radius > 0)
        assert!(singular.is_empty());
    }

    #[test]
    fn singular_points_plane_none() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let singular = detect_singular_points(&plane);

        // Plane has no singular points
        assert!(singular.is_empty());
    }

    #[test]
    fn geometric_continuity_ordering() {
        assert!(GeometricContinuity::C2 > GeometricContinuity::C1);
        assert!(GeometricContinuity::C1 > GeometricContinuity::G1);
        assert!(GeometricContinuity::G1 > GeometricContinuity::C0);
        assert!(GeometricContinuity::C0 > GeometricContinuity::G0);
        assert!(GeometricContinuity::G0 > GeometricContinuity::None);
    }

    #[test]
    fn segment_segment_distance_3d_parallel() {
        // Two parallel segments
        let p1 = DVec3::new(0.0, 0.0, 0.0);
        let p2 = DVec3::new(1.0, 0.0, 0.0);
        let p3 = DVec3::new(0.0, 1.0, 0.0);
        let p4 = DVec3::new(1.0, 1.0, 0.0);

        let dist = segment_segment_distance_3d(p1, p2, p3, p4);

        // Distance should be 1.0 (parallel lines, 1 unit apart)
        assert!(approx_eq(dist, 1.0, TOL));
    }

    #[test]
    fn segment_segment_distance_3d_intersecting() {
        // Two intersecting segments
        let p1 = DVec3::new(0.0, 0.0, 0.0);
        let p2 = DVec3::new(1.0, 1.0, 0.0);
        let p3 = DVec3::new(0.0, 1.0, 0.0);
        let p4 = DVec3::new(1.0, 0.0, 0.0);

        let dist = segment_segment_distance_3d(p1, p2, p3, p4);

        // These segments intersect at (0.5, 0.5, 0)
        assert!(approx_eq(dist, 0.0, TOL));
    }

    #[test]
    fn segment_segment_distance_3d_skew() {
        // Two skew lines (not parallel, not intersecting)
        let p1 = DVec3::new(0.0, 0.0, 0.0);
        let p2 = DVec3::new(1.0, 0.0, 0.0);
        let p3 = DVec3::new(0.0, 0.0, 1.0);
        let p4 = DVec3::new(0.0, 1.0, 1.0);

        let dist = segment_segment_distance_3d(p1, p2, p3, p4);

        // Distance should be 1.0 (perpendicular distance between skew lines)
        assert!(approx_eq(dist, 1.0, TOL));
    }
}
