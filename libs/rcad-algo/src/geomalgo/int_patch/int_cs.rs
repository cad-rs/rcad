// OCCT IntCurveSurface_HInter (canonic-curve vs quadric exact path) — the
// exact curve-surface intersection used by BoundedArc inside
// IntStart_SearchOnBoundaries (IntStart_SearchOnBoundaries.gxx L339-472).
//
// 1:1 Rust translation of the canonic conic-quadric exact path
// (PerformBounds -> PerformConicSurf* -> IntAna_IntConicQuad -> AppendIntAna).
// rcad data-model notes:
//   - Adaptor3d_CurveOnSurface (a 2D arc lifted to 3D on a quadric surface)
//     maps to a Curve3 whose parameterization equals the 2D arc parameter.
//   - The generic polygon/polyhedron path (for non-canonic curves and Torus
//     quadrics) is not reachable from ImpImp FF: quadric boundary arcs are
//     always lines/circles, and IntCS is only invoked when TypeQuad != Torus.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};

/// OCCT IntCurveSurface_TransitionOnCurve.hxx.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionOnCurve {
    In,
    Out,
    Tangent,
}

/// OCCT IntCurveSurface_IntersectionPoint.hxx.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionPoint {
    pub p: DVec3,
    pub u_surf: f64,
    pub v_surf: f64,
    pub w_curve: f64,
    pub tr_curv: TransitionOnCurve,
}

impl IntersectionPoint {
    pub fn pnt(&self) -> DVec3 {
        self.p
    }
    pub fn u(&self) -> f64 {
        self.u_surf
    }
    pub fn v(&self) -> f64 {
        self.v_surf
    }
    pub fn w(&self) -> f64 {
        self.w_curve
    }
    pub fn transition(&self) -> TransitionOnCurve {
        self.tr_curv
    }
}

/// OCCT IntCurveSurface_IntersectionSegment.hxx.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionSegment {
    pub p1: IntersectionPoint,
    pub p2: IntersectionPoint,
}

impl IntersectionSegment {
    pub fn first_point(&self) -> IntersectionPoint {
        self.p1
    }
    pub fn second_point(&self) -> IntersectionPoint {
        self.p2
    }
}

/// OCCT IntCurveSurface_HInter — the exact curve-surface intersection for a
/// canonic (Line/Circle/Ellipse/Parabola/Hyperbola) curve against a real
/// quadric surface.  Equivalent to the IntCS.Perform(HConS, GAHsurf) call in
/// IntStart_SearchOnBoundaries BoundedArc.
pub struct IntCurveSurface {
    done: bool,
    points: Vec<IntersectionPoint>,
    segments: Vec<IntersectionSegment>,
    is_parallel: bool,
}

/// OCCT IntCurveSurface_Inter.pxx constants.
const THE_TOLTANGENCY: f64 = 0.00000001;
const THE_TOLERANCE_ANGULAIRE: f64 = 1.0e-12;
const THE_TOLERANCE: f64 = 0.00000001;

impl IntCurveSurface {
    pub fn new() -> Self {
        IntCurveSurface {
            done: false,
            points: Vec::new(),
            segments: Vec::new(),
            is_parallel: false,
        }
    }

    /// OCCT IntCurveSurface_HInter::Perform(curve, surface).
    ///
    /// curve: a canonic 3D curve whose parameter equals the 2D arc parameter
    /// (the CurveOnSurface parameterization).  surface: the quadric surface.
    pub fn perform(&mut self, curve: &Curve3, surface: &Surface3) {
        self.points.clear();
        self.segments.clear();
        self.is_parallel = false;
        self.done = false;

        let quad = crate::geomalgo::int_surf::quadric::Quadric::from_surface3(surface);
        let Some(quad) = quad else { return };
        let stype = quad.type_quadric();

        // OCCT: if the conic lies IN the quadric (coincident) or is parallel,
        // IntAna_IntConicQuad reports IsInQuadric/IsParallel; AppendIntAna sets
        // myIsParallel and appends nothing.  No points -> Nbp=0 -> the caller
        // falls back to math_FunctionAllRoots (which detects the whole-arc
        // null interval as a restriction segment).
        let mut ana_points: Vec<(DVec3, f64)> = Vec::new(); // (3D point, W)
        match curve {
            Curve3::Line(l) => {
                let pts = intersect_line_quadric(l, &quad);
                if pts.is_none() {
                    self.done = false;
                    return;
                }
                ana_points = pts.unwrap();
            }
            Curve3::Circle(c) => {
                let pts = intersect_circle_quadric(c, &quad);
                if pts.is_none() {
                    self.done = false;
                    return;
                }
                ana_points = pts.unwrap();
            }
            _ => {
                // Non-canonic curve: the generic polygon path in OCCT would be
                // used; not reachable for quadric-quadric FF (boundary arcs are
                // lines/circles).  Mark not-done so BoundedArc falls back to
                // math_FunctionAllRoots (OCCT does the same when IntCS fails).
                self.done = false;
                return;
            }
        }

        self.done = true;

        // OCCT AppendIntAna -> ProcessIntAna: for each IntAna point, compute
        // surface params (ComputeParamsOnQuadric) and validate+transition via
        // ComputeAppendPoint.
        for (p, w) in ana_points {
            let (u, v) = compute_params_on_quadric(&quad, p);
            if let Some(pt) = compute_append_point(curve, w, surface, u, v) {
                self.points.push(pt);
            }
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }
    /// OCCT Point(Index) — 1-based.
    pub fn point(&self, index: usize) -> IntersectionPoint {
        self.points[index - 1]
    }
    /// OCCT NbSegments().
    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }
    /// OCCT Segment(Index) — 1-based.
    pub fn segment(&self, index: usize) -> IntersectionSegment {
        self.segments[index - 1]
    }
    /// OCCT myIsParallel.
    pub fn is_parallel(&self) -> bool {
        self.is_parallel
    }
}

/// OCCT IntCurveSurface_InterUtils::ComputeParamsOnQuadric
/// (IntCurveSurface_InterUtils.pxx L901-925).
fn compute_params_on_quadric(
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
    p: DVec3,
) -> (f64, f64) {
    quad.parameters(p)
}

/// OCCT IntCurveSurface_InterUtils::ComputeAppendPoint
/// (IntCurveSurface_InterUtils.pxx L1117-1176).
fn compute_append_point(
    curve: &Curve3,
    lw: f64,
    surface: &Surface3,
    su: f64,
    sv: f64,
) -> Option<IntersectionPoint> {
    let dom = curve.default_domain();
    let (w0, w1) = (dom[0], dom[1]);
    let surf_dom = surface.default_domain();
    let (u0, u1) = (surf_dom[0], surf_dom[1]);
    let (v0, v1) = (surf_dom[2], surf_dom[3]);

    let mut w = lw;
    let mut u = su;
    let mut v = sv;

    let is_circle = matches!(curve, Curve3::Circle(_));
    if curve.is_periodic() || is_circle {
        w = elclib_in_period(w, w0, w0 + curve_period(curve));
    }

    if (w0 - w) >= THE_TOLTANGENCY || (w - w1) >= THE_TOLTANGENCY {
        return None;
    }

    let s_u_periodic = surface.is_u_periodic()
        || matches!(surface, Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_));
    if s_u_periodic {
        u = elclib_in_period(u, u0, u0 + surface_period(surface, true));
    }
    if surface.is_v_periodic() {
        v = elclib_in_period(v, v0, v0 + surface_period(surface, false));
    }

    if (u0 - u) >= THE_TOLTANGENCY || (u - u1) >= THE_TOLTANGENCY {
        return None;
    }
    if (v0 - v) >= THE_TOLTANGENCY || (v - v1) >= THE_TOLTANGENCY {
        return None;
    }

    let tr = compute_transitions(curve, w, surface, u, v);
    let p = curve.point_at(w);
    Some(IntersectionPoint {
        p,
        u_surf: u,
        v_surf: v,
        w_curve: w,
        tr_curv: tr,
    })
}

/// OCCT IntCurveSurface_InterUtils::ComputeTransitions
/// (IntCurveSurface_InterUtils.pxx L856-895).
fn compute_transitions(
    curve: &Curve3,
    w: f64,
    surface: &Surface3,
    u: f64,
    v: f64,
) -> TransitionOnCurve {
    let (psurf, d1u, d1v) = surface.derivatives(u, v);
    let n_surf = d1u.cross(d1v);
    let d1u_curve = curve.derivative_at(w);
    let norm = n_surf.length();
    if norm > THE_TOLERANCE_ANGULAIRE
        && d1u_curve.length_squared() > THE_TOLERANCE_ANGULAIRE
    {
        let d1u_n = d1u_curve.normalize_or_zero();
        let cos_dir = n_surf.dot(d1u_n);
        let cos_dir = cos_dir / norm;
        if -cos_dir > THE_TOLERANCE_ANGULAIRE {
            // --Curve--->    <----Surface----
            TransitionOnCurve::In
        } else if cos_dir > THE_TOLERANCE_ANGULAIRE {
            // --Curve--->  ----Surface-->
            TransitionOnCurve::Out
        } else {
            TransitionOnCurve::Tangent
        }
    } else {
        TransitionOnCurve::Tangent
    }
}

/// OCCT ElCLib::InPeriod(X, A, B).
fn elclib_in_period(x: f64, a: f64, b: f64) -> f64 {
    let p = b - a;
    let mut x = x;
    while x < a {
        x += p;
    }
    while x > b {
        x -= p;
    }
    x
}

fn curve_period(curve: &Curve3) -> f64 {
    match curve {
        Curve3::Circle(_) => std::f64::consts::TAU,
        _ => curve.default_domain()[1] - curve.default_domain()[0],
    }
}

fn surface_period(surface: &Surface3, is_u: bool) -> f64 {
    if is_u {
        if surface.is_u_periodic() {
            std::f64::consts::TAU
        } else {
            0.0
        }
    } else if surface.is_v_periodic() {
        std::f64::consts::TAU
    } else {
        0.0
    }
}

/// OCCT IntAna_IntConicQuad (line vs quadric) — the exact analytic
/// intersection.  Returns None when the conic is parallel/coincident (the
/// isAnaProcessed=false path in OCCT also falls back); otherwise the list of
/// (point, W) where W is the parameter on the curve.
fn intersect_line_quadric(
    line: &rcad_kernel::geom::Line3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<Vec<(DVec3, f64)>> {
    use super::curve_surface as cs;
    let tr = [f64::NEG_INFINITY, f64::INFINITY];
    let hits: Vec<(DVec3, f64)> = match quad.type_quadric() {
        crate::geomalgo::int_surf::quadric::QuadricType::Plane => {
            // IntAna_IntConicQuad(Line, Plane) — direct analytic intersection.
            let pl = quad.plane();
            let d = pl.normal.dot(line.origin) - pl.normal.dot(pl.origin);
            let dn = pl.normal.dot(line.direction);
            if dn.abs() < 1e-12 {
                Vec::new() // parallel
            } else {
                let t = -d / dn;
                vec![(line.origin + line.direction * t, t)]
            }
        }
        crate::geomalgo::int_surf::quadric::QuadricType::Cylinder => {
            let cyl = quad.cylinder();
            cs::intersect_line_cylinder_with_tol(line, tr, &cyl, 1e-7)
                .iter()
                .map(|h| (h.point, h.curve_param))
                .collect()
        }
        crate::geomalgo::int_surf::quadric::QuadricType::Sphere => {
            let sph = quad.sphere();
            cs::intersect_line_sphere_with_tol(line, tr, &sph, 1e-7)
                .iter()
                .map(|h| (h.point, h.curve_param))
                .collect()
        }
        crate::geomalgo::int_surf::quadric::QuadricType::Cone => {
            let con = quad.cone();
            cs::intersect_line_cone_with_tol(line, tr, &con, 1e-7)
                .iter()
                .map(|h| (h.point, h.curve_param))
                .collect()
        }
        _ => return None,
    };
    Some(hits)
}

/// OCCT IntAna_IntConicQuad (circle vs quadric).
fn intersect_circle_quadric(
    circle: &rcad_kernel::geom::Circle3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<Vec<(DVec3, f64)>> {
    // OCCT IntCurveSurface_HInter -> IntAna_IntConicQuad(Circle, Quadric):
    // analytic conic-quadric intersection (no sampling/Newton).
    // The Plane surface uses the (Circ, Pln, Tolang, Tol) overload with
    // THE_TOLERANCE_ANGULAIRE / THE_TOLERANCE (IntCurveSurface_Inter.pxx
    // L731-735); the Cylinder/Cone/Sphere quadrics use the (Circ, Quadric)
    // overload which takes no tolerance (Eps is internal, 1.5e-12).
    let quad_type = quad.type_quadric();
    let res = if quad_type == crate::geomalgo::int_surf::quadric::QuadricType::Plane {
        let pl = quad.plane();
        let (_, in_quadric, pts) = super::int_conic_quad::intersect_circle_plane(
            circle,
            &pl,
            THE_TOLERANCE_ANGULAIRE,
            THE_TOLERANCE,
        );
        if in_quadric {
            // The circle lies entirely on the plane: no isolated points
            // (BoundedArc treats it as an all-arc solution segment).
            Some((true, Vec::new()))
        } else {
            Some((false, pts))
        }
    } else {
        super::int_conic_quad::intersect_circle_quadric(circle, quad)
    }?;
    let (in_quadric, pts) = res;
    let _ = in_quadric;
    Some(pts)
}
