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

impl Default for IntersectionPoint {
    /// OCCT IntCurveSurface_IntersectionPoint() (IntCurveSurface_IntersectionPoint.cxx
    /// L21-28): zero parameters and the Transition initialized to Tangent.
    fn default() -> Self {
        IntersectionPoint {
            p: DVec3::ZERO,
            u_surf: 0.0,
            v_surf: 0.0,
            w_curve: 0.0,
            tr_curv: TransitionOnCurve::Tangent,
        }
    }
}

impl IntersectionPoint {
    /// OCCT IntCurveSurface_IntersectionPoint(P, USurf, VSurf, UCurv,
    /// TrOnCurv) (L30-42).
    pub fn new(
        p: DVec3,
        u_surf: f64,
        v_surf: f64,
        w_curve: f64,
        tr_curv: TransitionOnCurve,
    ) -> Self {
        IntersectionPoint {
            p,
            u_surf,
            v_surf,
            w_curve,
            tr_curv,
        }
    }

    /// OCCT IntCurveSurface_IntersectionPoint::SetValues (L44-55).
    pub fn set_values(
        &mut self,
        p: DVec3,
        u_surf: f64,
        v_surf: f64,
        w_curve: f64,
        tr_curv: TransitionOnCurve,
    ) {
        self.p = p;
        self.u_surf = u_surf;
        self.v_surf = v_surf;
        self.w_curve = w_curve;
        self.tr_curv = tr_curv;
    }

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
    /// OCCT Adaptor3d_CurveOnSurface::FirstParameter/LastParameter — the 2D
    /// boundary arc parameter domain.  The IntCS W() must be validated and
    /// period-wrapped against this domain (IntCurveSurface_InterUtils
    /// ComputeAppendPoint uses CurveTool::FirstParameter/LastParameter),
    /// NOT against the 3D curve's geometric domain.
    curve_domain: Option<[f64; 2]>,
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
            curve_domain: None,
        }
    }

    /// OCCT IntCurveSurface_HInter::Perform(curve, surface).
    ///
    /// curve: a canonic 3D curve whose parameter equals the 2D arc parameter
    /// (the CurveOnSurface parameterization).  surface: the quadric surface.
    pub fn perform(&mut self, curve: &Curve3, surface: &Surface3, curve_domain: [f64; 2]) {
        // OCCT IntCurveSurface_HInter::Perform(Adaptor3d_CurveOnSurface, ...):
        // W0/W1 come from CurveTool::FirstParameter/LastParameter on the
        // curve-on-surface, i.e. the 2D boundary arc parameter domain.
        self.curve_domain = Some(curve_domain);
        self.points.clear();
        self.segments.clear();
        self.is_parallel = false;
        self.done = false;

        let quad = crate::geomalgo::int_surf::quadric::Quadric::from_surface3(surface);
        let Some(quad) = quad else { return };
        let stype = quad.type_quadric();

        // OCCT ProcessIntAna (IntCurveSurface_InterUtils.pxx L1188-1228): if
        // the conic lies IN the quadric (coincident) or is parallel,
        // IntAna_IntConicQuad reports IsInQuadric/IsParallel; the caller sets
        // myIsParallel and appends nothing.  No points -> Nbp=0 -> the caller
        // falls back to math_FunctionAllRoots (which detects the whole-arc
        // null interval as a restriction segment).
        let mut ana_points: Vec<(DVec3, f64)> = Vec::new(); // (3D point, W)
        self.is_parallel = false;
        match curve {
            Curve3::Line(l) => {
                let (in_quadric, pts) = match intersect_line_quadric(l, &quad) {
                    None => {
                        self.done = false;
                        return;
                    }
                    Some(r) => r,
                };
                if in_quadric {
                    self.is_parallel = true;
                }
                ana_points = pts;
            }
            Curve3::Circle(c) => {
                let (in_quadric, pts) = match intersect_circle_quadric(c, &quad) {
                    None => {
                        self.done = false;
                        return;
                    }
                    Some(r) => r,
                };
                if in_quadric {
                    self.is_parallel = true;
                }
                ana_points = pts;
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
        let w0_w1 = self.curve_domain.unwrap_or(curve.default_domain());
        for (p, w) in ana_points {
            let (u, v) = compute_params_on_quadric(&quad, p);
            if let Some(pt) = compute_append_point(curve, w, surface, u, v, w0_w1) {
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
    curve_domain: [f64; 2],
) -> Option<IntersectionPoint> {
    // OCCT IntCurveSurface_InterUtils::ComputeAppendPoint: W0/W1 =
    // CurveTool::FirstParameter/LastParameter — on the curve-on-surface these
    // are the 2D boundary arc parameter bounds (NOT the 3D curve geometric
    // domain).  The W range check and period wrap use this domain.
    let (w0, w1) = (curve_domain[0], curve_domain[1]);
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

/// OCCT IntAna_IntConicQuad::Perform(const gp_Lin&, const IntAna_Quadric&)
/// (IntAna_IntConicQuad.cxx L67-127) — the exact analytic line-quadric
/// intersection.  The line direction/location are substituted into the
/// absolute-frame quadric coefficients, giving a quadratic in t solved by
/// math_DirectPolynomialRoots.  Returns None when not done, otherwise
/// (in_quadric, (point, W) pairs) where W is the line parameter.
pub(crate) fn intersect_line_quadric(
    line: &rcad_kernel::geom::Line3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    use rcad_kernel::math::direct_polynomial_roots::DirectPolynomialRoots;

    let co = super::int_conic_quad::quadric_frame_coefs(quad)?;
    let (qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, qcte) = (
        co.xx, co.yy, co.zz, co.xy, co.xz, co.yz, co.x, co.y, co.z, co.cte,
    );
    let (lx0, ly0, lz0) = (line.origin.x, line.origin.y, line.origin.z);
    let (lx, ly, lz) = (line.direction.x, line.direction.y, line.direction.z);

    // OCCT L95-96: A0.
    let a0 = qcte + qxx * lx0 * lx0 + qyy * ly0 * ly0 + qzz * lz0 * lz0
        + 2.0 * (lx0 * (qx + qxy * ly0 + qxz * lz0) + ly0 * (qy + qyz * lz0) + qz * lz0);

    // OCCT L98-101: A1.
    let a1 = 2.0
        * (lx * (qx + qxx * lx0 + qxy * ly0 + qxz * lz0)
            + ly * (qy + qxy * lx0 + qyy * ly0 + qyz * lz0)
            + lz * (qz + qxz * lx0 + qyz * ly0 + qzz * lz0));

    // OCCT L103-104: A2.
    let a2 = qxx * lx * lx + qyy * ly * ly + qzz * lz * lz
        + 2.0 * (lx * (qxy * ly + qxz * lz) + qyz * ly * lz);

    let lin_quad_pol = DirectPolynomialRoots::new_quadratic(a2, a1, a0);
    if !lin_quad_pol.is_done() {
        return None;
    }
    if lin_quad_pol.infinite_roots() {
        // OCCT L111-115: the line lies in the quadric (inquadric = true) — no
        // isolated points; the caller falls back to the function-root path.
        return Some((true, Vec::new()));
    }

    let mut pts = Vec::new();
    for i in 1..=lin_quad_pol.nb_solutions() {
        let t = lin_quad_pol.value(i);
        pts.push((line.origin + line.direction * t, t));
    }
    Some((false, pts))
}

/// OCCT IntAna_IntConicQuad (circle vs quadric).  Returns
/// (in_quadric, points): the Plane surface uses the (Circ, Pln, Tolang, Tol)
/// overload; the Cylinder/Cone/Sphere quadrics use the (Circ, Quadric)
/// overload.
fn intersect_circle_quadric(
    circle: &rcad_kernel::geom::Circle3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    // OCCT IntCurveSurface_HInter -> IntAna_IntConicQuad(Circle, Quadric):
    // analytic conic-quadric intersection (no sampling/Newton).
    // The Plane surface uses the (Circ, Pln, Tolang, Tol) overload with
    // THE_TOLERANCE_ANGULAIRE / THE_TOLERANCE (IntCurveSurface_Inter.pxx
    // L731-735); the Cylinder/Cone/Sphere quadrics use the (Circ, Quadric)
    // overload which takes no tolerance (Eps is internal, 1.5e-12).
    let quad_type = quad.type_quadric();
    if quad_type == crate::geomalgo::int_surf::quadric::QuadricType::Plane {
        let pl = quad.plane();
        let (parallel, in_quadric, pts) = super::int_conic_quad::intersect_circle_plane(
            circle,
            &pl,
            THE_TOLERANCE_ANGULAIRE,
            THE_TOLERANCE,
        );
        if parallel || in_quadric {
            // The circle lies entirely on the plane or is parallel to it: no
            // isolated points (BoundedArc treats it as an all-arc solution
            // segment).
            Some((true, Vec::new()))
        } else {
            Some((false, pts))
        }
    } else {
        super::int_conic_quad::intersect_circle_quadric(circle, quad)
    }
}
