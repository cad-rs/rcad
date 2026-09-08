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
//   - The ChFi3d corner-tail flow (ChFi3d_Builder_C1.cxx Update, L191-280)
//     performs the same IntCurveSurface_HInter::Perform(ct, fb) with
//     Ellipse/Parabola/Hyperbola interference curves against Plane/Cylinder/
//     Cone/Sphere prolongation surfaces, so the full PerformBounds conic
//     switch (Inter.pxx L121-139) is dispatched here.

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
            Curve3::Ellipse(e) => {
                // OCCT PerformBounds case GeomAbs_Ellipse (Inter.pxx L129-131)
                // -> PerformConicSurfEllipse.
                let (in_quadric, pts) = match perform_conic_surf_ellipse(e, &quad) {
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
            Curve3::Parabola(p) => {
                // OCCT PerformBounds case GeomAbs_Parabola (Inter.pxx L132-134)
                // -> PerformConicSurfParabola.
                let (in_quadric, pts) = match perform_conic_surf_parabola(p, &quad) {
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
            Curve3::Hyperbola(h) => {
                // OCCT PerformBounds case GeomAbs_Hyperbola (Inter.pxx L135-137)
                // -> PerformConicSurfHyperbola.
                let (in_quadric, pts) = match perform_conic_surf_hyperbola(h, &quad) {
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
                // Non-canonic curve (Bezier/BSpline/...): the OCCT PerformBounds
                // default branch (IntCurveSurface_Inter.pxx L139-180) runs the
                // polygonal/polyhedron interference, or the exact quadric solver
                // (TheQuadCurvExactHInter) for Plane/Cylinder/Cone/Sphere — both
                // carried by geomalgo::int_curve_surface over adaptors.  This
                // concrete IntCurveSurface carrier reports not-done there, so the
                // caller takes its documented fallback.
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
    /// OCCT IntCurveSurface_Intersection::NbPoints
    /// (IntCurveSurface_Intersection.cxx L45-52): StdFail_NotDone when not
    /// done.
    pub fn nb_points(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: IntCurveSurface_Intersection::NbPoints");
        }
        self.points.len()
    }
    /// OCCT IntCurveSurface_Intersection::Point(N)
    /// (IntCurveSurface_Intersection.cxx L65-72): StdFail_NotDone when not
    /// done; lpnt.Value(N) is the 1-based sequence access and raises
    /// Standard_OutOfRange outside [1, NbPoints()].
    pub fn point(&self, index: usize) -> IntersectionPoint {
        if !self.done {
            panic!("StdFail_NotDone: IntCurveSurface_Intersection::Point");
        }
        if index < 1 || index > self.points.len() {
            panic!(
                "Standard_OutOfRange: IntCurveSurface_Intersection::Point({}) with NbPoints()={}",
                index,
                self.points.len()
            );
        }
        self.points[index - 1]
    }
    /// OCCT IntCurveSurface_Intersection::NbSegments
    /// (IntCurveSurface_Intersection.cxx L55-62): StdFail_NotDone when not
    /// done.
    pub fn nb_segments(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: IntCurveSurface_Intersection::NbSegments");
        }
        self.segments.len()
    }
    /// OCCT IntCurveSurface_Intersection::Segment(N)
    /// (IntCurveSurface_Intersection.cxx L75-82): StdFail_NotDone when not
    /// done; lseg.Value(N) is the 1-based sequence access and raises
    /// Standard_OutOfRange outside [1, NbSegments()].
    pub fn segment(&self, index: usize) -> IntersectionSegment {
        if !self.done {
            panic!("StdFail_NotDone: IntCurveSurface_Intersection::Segment");
        }
        if index < 1 || index > self.segments.len() {
            panic!(
                "Standard_OutOfRange: IntCurveSurface_Intersection::Segment({}) with NbSegments()={}",
                index,
                self.segments.len()
            );
        }
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

    // OCCT: CurveTool::IsPeriodic || aCType == GeomAbs_Circle ||
    // aCType == GeomAbs_Ellipse (InterUtils.pxx L1142-1145).
    let is_circle_or_ellipse = matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_));
    if curve.is_periodic() || is_circle_or_ellipse {
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
        // OCCT Adaptor3d_Curve::Period() = LastParameter - FirstParameter;
        // Geom_Circle/Geom_Ellipse domains are always [0, 2*PI].
        Curve3::Circle(_) | Curve3::Ellipse(_) => std::f64::consts::TAU,
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

/// OCCT IntCurveSurface_Inter.pxx ProcessIntAna tail (InterUtils.pxx
/// L1188-1231): not-done -> None; IsInQuadric()/IsParallel() -> (true, no
/// points); otherwise the (Point(i), ParamOnConic(i)) pairs (1-based).
fn finish_int_ana(
    ana: &super::int_conic_quad::IntConicQuad,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    if !ana.is_done() {
        return None;
    }
    if ana.is_in_quadric() || ana.is_parallel() {
        return Some((true, Vec::new()));
    }
    let pts = (1..=ana.nb_points())
        .map(|i| (ana.point(i), ana.param_on_conic(i)))
        .collect();
    Some((false, pts))
}

/// OCCT IntCurveSurface_Inter.pxx PerformConicSurfEllipse (L769-820): the
/// Plane surface uses the (Elips, Pln, Tolang, Tol) overload (IntAna
/// IntConicQuad.cxx L562-565 delegates it to the quadric path), the
/// Cylinder/Cone/Sphere quadrics use the (Elips, Quadric) overload.  The
/// default branch (Torus and non-quadrics) falls back to the polygonal
/// interference path (IntCurveSurface_Inter.pxx L816-819, carried by
/// geomalgo::int_curve_surface::inter_impl); this concrete IntCurveSurface
/// carrier reports not-done there.  Returns None when not done, otherwise
/// (in_quadric, (point, W) pairs).
fn perform_conic_surf_ellipse(
    ellipse: &rcad_kernel::geom::Ellipse3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    use super::int_conic_quad::IntConicQuad;
    let quad_type = quad.type_quadric();
    if quad_type == crate::geomalgo::int_surf::quadric::QuadricType::Plane {
        let ana = IntConicQuad::new_ellipse_plane(ellipse, &quad.plane());
        finish_int_ana(&ana)
    } else if matches!(
        quad_type,
        crate::geomalgo::int_surf::quadric::QuadricType::Cylinder
            | crate::geomalgo::int_surf::quadric::QuadricType::Cone
            | crate::geomalgo::int_surf::quadric::QuadricType::Sphere
    ) {
        let ana = IntConicQuad::new_ellipse_quadric(ellipse, quad);
        finish_int_ana(&ana)
    } else {
        None
    }
}

/// OCCT IntCurveSurface_Inter.pxx PerformConicSurfParabola (L822-880): the
/// Plane surface uses the (Parab, Pln, Tolang) overload (IntAna
/// IntConicQuad.cxx L567-570 delegates it to the quadric path), the
/// Cylinder/Cone/Sphere quadrics use the (Parab, Quadric) overload.  The
/// default branch falls back to the polyhedron path (IntCurveSurface_Inter.pxx
/// L852-878); this carrier reports not-done there.
fn perform_conic_surf_parabola(
    parabola: &rcad_kernel::geom::Parabola3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    use super::int_conic_quad::IntConicQuad;
    let quad_type = quad.type_quadric();
    if quad_type == crate::geomalgo::int_surf::quadric::QuadricType::Plane {
        let ana = IntConicQuad::new_parabola_plane(parabola, &quad.plane());
        finish_int_ana(&ana)
    } else if matches!(
        quad_type,
        crate::geomalgo::int_surf::quadric::QuadricType::Cylinder
            | crate::geomalgo::int_surf::quadric::QuadricType::Cone
            | crate::geomalgo::int_surf::quadric::QuadricType::Sphere
    ) {
        let ana = IntConicQuad::new_parabola_quadric(parabola, quad);
        finish_int_ana(&ana)
    } else {
        None
    }
}

/// OCCT IntCurveSurface_Inter.pxx PerformConicSurfHyperbola (L882-940): the
/// Plane surface uses the (Hypr, Pln, Tolang) overload (IntAna
/// IntConicQuad.cxx L572-575 delegates it to the quadric path), the
/// Cylinder/Cone/Sphere quadrics use the (Hypr, Quadric) overload.  The
/// default branch falls back to the polyhedron path (IntCurveSurface_Inter.pxx
/// L912-938); this carrier reports not-done there.
fn perform_conic_surf_hyperbola(
    hyperbola: &rcad_kernel::geom::Hyperbola3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    use super::int_conic_quad::IntConicQuad;
    let quad_type = quad.type_quadric();
    if quad_type == crate::geomalgo::int_surf::quadric::QuadricType::Plane {
        let ana = IntConicQuad::new_hyperbola_plane(hyperbola, &quad.plane());
        finish_int_ana(&ana)
    } else if matches!(
        quad_type,
        crate::geomalgo::int_surf::quadric::QuadricType::Cylinder
            | crate::geomalgo::int_surf::quadric::QuadricType::Cone
            | crate::geomalgo::int_surf::quadric::QuadricType::Sphere
    ) {
        let ana = IntConicQuad::new_hyperbola_quadric(hyperbola, quad);
        finish_int_ana(&ana)
    } else {
        None
    }
}
