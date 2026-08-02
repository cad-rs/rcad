// OCCT IntPatch_TheSOnBounds = IntStart_SearchOnBoundaries + IntPatch_ArcFunction
// + IntPatch_ThePathPointOfTheSOnBounds + IntPatch_TheSegmentOfTheSOnBounds.
//
// IntStart_SearchOnBoundaries.gxx (1232 lines) instantiated as
// IntPatch_TheSOnBounds.  1:1 Rust translation.
//
// rcad data-model notes:
//   - Adaptor3d_TopolTool (the face domain) maps to the corrected FF UV
//     rectangle [u_min, u_max, v_min, v_max]; its boundary arcs are the four
//     edges of the rectangle (2D lines in UV space) and its vertices are the
//     four corners.
//   - Adaptor2d_Curve2d (boundary arc) -> Curve2d.
//   - Adaptor3d_HVertex (domain corner) -> DomainVertex { u, v }.
//   - Adaptor3d_Surface -> Surface3; IntSurf_Quadric -> Quadric.
//   - The generic TreatLC/IsRegularity BRep-edge branches are not reachable:
//     rcad's domain has no BRep edge (aDomain->Edge() == nullptr), so both
//     return their early value, exactly as OCCT does for a null edge.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, SurfaceEval};

use super::math_roots::{BrentMinimum, FunctionAllRoots, FunctionSample, FunctionValue, FunctionWithDerivative};
use super::int_cs::{IntCurveSurface, TransitionOnCurve};
use crate::topalgo::int_surf::quadric::Quadric;

// =====================================================================
// IntPatch_ThePathPointOfTheSOnBounds
// =====================================================================

/// OCCT IntPatch_ThePathPointOfTheSOnBounds.hxx — a solution point on a
/// boundary arc: 3D point, tolerance, the arc, its parameter, and (when not
/// new) the domain vertex it coincides with.
#[derive(Debug, Clone)]
pub struct PathPoint {
    point: DVec3,
    tol: f64,
    isnew: bool,
    vtx: Option<DomainVertex>,
    arc: Curve2d,
    param: f64,
}

impl PathPoint {
    pub fn new() -> Self {
        PathPoint {
            point: DVec3::ZERO,
            tol: 0.0,
            isnew: true,
            vtx: None,
            arc: Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::X,
            }),
            param: 0.0,
        }
    }

    /// OCCT SetValue(P, Tol, V, A, Parameter) — with vertex (isnew = false).
    pub fn set_value_with_vtx(
        &mut self,
        p: DVec3,
        tol: f64,
        v: DomainVertex,
        a: Curve2d,
        param: f64,
    ) {
        self.isnew = false;
        self.point = p;
        self.tol = tol;
        self.vtx = Some(v);
        self.arc = a;
        self.param = param;
    }

    /// OCCT SetValue(P, Tol, A, Parameter) — no vertex (isnew = true).
    pub fn set_value(&mut self, p: DVec3, tol: f64, a: Curve2d, param: f64) {
        self.isnew = true;
        self.point = p;
        self.tol = tol;
        self.arc = a;
        self.param = param;
    }

    /// OCCT Value().
    pub fn value(&self) -> DVec3 {
        self.point
    }
    /// OCCT Tolerance().
    pub fn tolerance(&self) -> f64 {
        self.tol
    }
    /// OCCT IsNew().
    pub fn is_new(&self) -> bool {
        self.isnew
    }
    /// OCCT Vertex().
    pub fn vertex(&self) -> DomainVertex {
        self.vtx.unwrap_or(DomainVertex { u: 0.0, v: 0.0 })
    }
    /// OCCT Arc().
    pub fn arc(&self) -> &Curve2d {
        &self.arc
    }
    /// OCCT Parameter().
    pub fn parameter(&self) -> f64 {
        self.param
    }
}

// =====================================================================
// IntPatch_TheSegmentOfTheSOnBounds
// =====================================================================

/// OCCT IntPatch_TheSegmentOfTheSOnBounds.hxx — a solution segment on a
/// boundary arc (an interval where F is null on the whole arc).
#[derive(Debug, Clone)]
pub struct Segment {
    arc: Curve2d,
    has_fp: bool,
    the_fp: PathPoint,
    has_lp: bool,
    the_lp: PathPoint,
}

impl Segment {
    /// OCCT SetValue(A).
    pub fn set_value(&mut self, a: Curve2d) {
        self.has_fp = false;
        self.has_lp = false;
        self.arc = a;
    }

    /// OCCT SetLimitPoint(V, First).
    pub fn set_limit_point(&mut self, v: PathPoint, first: bool) {
        if first {
            self.has_fp = true;
            self.the_fp = v;
        } else {
            self.has_lp = true;
            self.the_lp = v;
        }
    }

    /// OCCT Curve().
    pub fn curve(&self) -> &Curve2d {
        &self.arc
    }
    /// OCCT HasFirstPoint().
    pub fn has_first_point(&self) -> bool {
        self.has_fp
    }
    /// OCCT FirstPoint().
    pub fn first_point(&self) -> &PathPoint {
        &self.the_fp
    }
    /// OCCT HasLastPoint().
    pub fn has_last_point(&self) -> bool {
        self.has_lp
    }
    /// OCCT LastPoint().
    pub fn last_point(&self) -> &PathPoint {
        &self.the_lp
    }
}

// =====================================================================
// IntPatch_ArcFunction
// =====================================================================

/// OCCT IntPatch_ArcFunction.hxx / .cxx — F(t) = Q(P(C(t))), the algebraic
/// distance from the point C(t) on the parametric surface to the quadric Q.
/// Implements math_FunctionWithDerivative.
pub struct ArcFunction {
    my_arc: Curve2d,
    my_surf: Surface3,
    my_quad: Quadric,
    ptsol: DVec3,
    seqpt: Vec<DVec3>,
}

impl ArcFunction {
    /// OCCT default constructor.
    pub fn new() -> Self {
        ArcFunction {
            my_arc: Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::X,
            }),
            my_surf: Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
            my_quad: Quadric::new(),
            ptsol: DVec3::ZERO,
            seqpt: Vec::new(),
        }
    }

    /// OCCT SetQuadric(Q).
    pub fn set_quadric(&mut self, q: Quadric) {
        self.my_quad = q;
    }
    /// OCCT Set(S) — the parametric surface.
    pub fn set_surface(&mut self, s: Surface3) {
        self.my_surf = s;
    }
    /// OCCT Set(A) — the boundary arc.
    pub fn set_arc(&mut self, a: Curve2d) {
        self.my_arc = a;
    }
    /// OCCT Value(X, F).
    pub fn value(&mut self, x: f64) -> Option<f64> {
        let p2d = self.my_arc.point_at(x);
        let ptsol = self.my_surf.point_at(p2d.x, p2d.y);
        if !ptsol.is_finite() {
            return None;
        }
        self.ptsol = ptsol;
        Some(self.my_quad.distance(ptsol))
    }
    /// OCCT Derivative(X, D).
    pub fn derivative(&mut self, x: f64) -> Option<f64> {
        let p2d = self.my_arc.point_at(x);
        let d2d = self.my_arc.derivative_at(x);
        let (ptsol, d1u, d1v) = self.my_surf.derivatives(p2d.x, p2d.y);
        if !ptsol.is_finite() {
            return None;
        }
        self.ptsol = ptsol;
        let v = d2d.x * d1u + d2d.y * d1v;
        Some(v.dot(self.my_quad.gradient(ptsol)))
    }
    /// OCCT Values(X, F, D).
    pub fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let p2d = self.my_arc.point_at(x);
        let d2d = self.my_arc.derivative_at(x);
        let (ptsol, d1u, d1v) = self.my_surf.derivatives(p2d.x, p2d.y);
        if !ptsol.is_finite() {
            return None;
        }
        self.ptsol = ptsol;
        let v1 = d2d.x * d1u + d2d.y * d1v;
        let (f, v2) = self.my_quad.val_and_grad(ptsol);
        let d = v1.dot(v2);
        Some((f, d))
    }
    /// OCCT GetStateNumber().
    pub fn get_state_number(&mut self) -> i32 {
        self.seqpt.push(self.ptsol);
        self.seqpt.len() as i32
    }
    /// OCCT Valpoint(Index) — 1-based.
    pub fn valpoint(&self, index: i32) -> DVec3 {
        self.seqpt[(index - 1) as usize]
    }
    /// OCCT Quadric().
    pub fn quadric(&self) -> &Quadric {
        &self.my_quad
    }
    /// OCCT Surface().
    pub fn surface(&self) -> &Surface3 {
        &self.my_surf
    }
    /// OCCT Arc().
    pub fn arc(&self) -> &Curve2d {
        &self.my_arc
    }
    /// OCCT LastComputedPoint().
    pub fn last_computed_point(&self) -> DVec3 {
        self.ptsol
    }
    /// OCCT NbSamples().
    pub fn nb_samples(&self) -> i32 {
        let u = nb_samples_u(&self.my_surf);
        let v = nb_samples_v(&self.my_surf);
        let a = nb_samples_on_arc(&self.my_arc);
        u.max(v).max(a)
    }
}

impl FunctionValue for ArcFunction {
    fn value(&mut self, x: f64) -> Option<f64> {
        ArcFunction::value(self, x)
    }
}

impl FunctionWithDerivative for ArcFunction {
    fn derivative(&mut self, x: f64) -> Option<f64> {
        ArcFunction::derivative(self, x)
    }
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        ArcFunction::values(self, x)
    }
    fn get_state_number(&mut self) -> i32 {
        ArcFunction::get_state_number(self)
    }
}

/// OCCT IntPatch_HInterTool::NbSamplesU (IntPatch_HInterTool.cxx L76-112).
fn nb_samples_u(s: &Surface3) -> i32 {
    match super::classify_surface_type(s) {
        super::GeomAbsSurfaceType::Plane => 2,
        super::GeomAbsSurfaceType::Torus => 20,
        super::GeomAbsSurfaceType::Cylinder
        | super::GeomAbsSurfaceType::Cone
        | super::GeomAbsSurfaceType::Sphere
        | super::GeomAbsSurfaceType::BezierSurface
        | super::GeomAbsSurfaceType::BSplineSurface
        | super::GeomAbsSurfaceType::SurfaceOfRevolution
        | super::GeomAbsSurfaceType::SurfaceOfExtrusion
        | super::GeomAbsSurfaceType::OffsetSurface
        | super::GeomAbsSurfaceType::OtherSurface => 10,
    }
}

/// OCCT IntPatch_HInterTool::NbSamplesV (IntPatch_HInterTool.cxx L37-74).
fn nb_samples_v(s: &Surface3) -> i32 {
    match super::classify_surface_type(s) {
        super::GeomAbsSurfaceType::Plane => 2,
        super::GeomAbsSurfaceType::Cylinder
        | super::GeomAbsSurfaceType::Cone
        | super::GeomAbsSurfaceType::Sphere
        | super::GeomAbsSurfaceType::Torus
        | super::GeomAbsSurfaceType::SurfaceOfRevolution
        | super::GeomAbsSurfaceType::SurfaceOfExtrusion => 15,
        super::GeomAbsSurfaceType::BezierSurface
        | super::GeomAbsSurfaceType::BSplineSurface
        | super::GeomAbsSurfaceType::OffsetSurface
        | super::GeomAbsSurfaceType::OtherSurface => 10,
    }
}

/// OCCT IntPatch_HInterTool::NbSamplesOnArc (IntPatch_HInterTool.cxx L232-261).
fn nb_samples_on_arc(a: &Curve2d) -> i32 {
    match a {
        Curve2d::Line(_) => 2,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) | Curve2d::Hyperbola(_) | Curve2d::Parabola(_) => 10,
        Curve2d::Bezier(_) => 10,
        Curve2d::BSpline(_) => 10,
        _ => 10,
    }
}

// =====================================================================
// Domain (Adaptor3d_TopolTool equivalent for the FF UV rectangle)
// =====================================================================

/// A corner of the UV domain rectangle (Adaptor3d_HVertex equivalent).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DomainVertex {
    pub u: f64,
    pub v: f64,
}

/// The domain of restriction of a surface: the corrected FF UV rectangle.
/// Its boundary arcs are the four edges of the rectangle.
pub struct Domain {
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    // Current arc index for Init/More/Value/Next.
    arc_idx: usize,
    // Current vertex iterator state for the current arc.
    vtx_idx: usize,
    // The arc the vertex iterator is attached to (index).
    vtx_arc: usize,
}

impl Domain {
    pub fn new(u_min: f64, u_max: f64, v_min: f64, v_max: f64) -> Self {
        Domain {
            u_min,
            u_max,
            v_min,
            v_max,
            arc_idx: 0,
            vtx_idx: 0,
            vtx_arc: 0,
        }
    }

    /// OCCT TopolTool::Init().
    pub fn init(&mut self) {
        self.arc_idx = 0;
    }
    /// OCCT TopolTool::More().
    pub fn more(&self) -> bool {
        self.arc_idx < 4
    }
    /// OCCT TopolTool::Value() — the current boundary arc (2D curve in UV).
    pub fn value(&self) -> Curve2d {
        let (o, d) = match self.arc_idx {
            0 => (DVec2::new(0.0, self.v_min), DVec2::new(1.0, 0.0)),
            1 => (DVec2::new(0.0, self.v_max), DVec2::new(1.0, 0.0)),
            2 => (DVec2::new(self.u_min, 0.0), DVec2::new(0.0, 1.0)),
            _ => (DVec2::new(self.u_max, 0.0), DVec2::new(0.0, 1.0)),
        };
        Curve2d::Line(rcad_kernel::geom::Line2d { origin: o, direction: d })
    }
    /// OCCT TopolTool::Next().
    pub fn next(&mut self) {
        self.arc_idx += 1;
    }

    /// OCCT IntPatch_HInterTool::Bounds(A, Ufirst, Ulast).
    pub fn bounds(&self, a: &Curve2d) -> (f64, f64) {
        // Arc 0/1: V=const, U in [u_min,u_max]; Arc 2/3: U=const, V in [v_min,v_max].
        if matches!(a, Curve2d::Line(l) if l.direction.y.abs() > 0.5) {
            (self.v_min, self.v_max)
        } else {
            (self.u_min, self.u_max)
        }
    }

    /// OCCT TopolTool::Initialize(A) — attach the vertex iterator to an arc.
    pub fn initialize(&mut self, a: &Curve2d) {
        self.vtx_arc = self.arc_of(a);
        self.vtx_idx = 0;
    }
    /// OCCT TopolTool::InitVertexIterator().
    pub fn init_vertex_iterator(&mut self) {
        self.vtx_idx = 0;
    }
    /// OCCT TopolTool::MoreVertex() — each arc has two endpoint corners.
    pub fn more_vertex(&self) -> bool {
        self.vtx_idx < 2
    }
    /// OCCT TopolTool::Vertex().
    pub fn vertex(&self) -> DomainVertex {
        match self.vtx_arc {
            0 => {
                if self.vtx_idx == 0 {
                    DomainVertex { u: self.u_min, v: self.v_min }
                } else {
                    DomainVertex { u: self.u_max, v: self.v_min }
                }
            }
            1 => {
                if self.vtx_idx == 0 {
                    DomainVertex { u: self.u_min, v: self.v_max }
                } else {
                    DomainVertex { u: self.u_max, v: self.v_max }
                }
            }
            2 => {
                if self.vtx_idx == 0 {
                    DomainVertex { u: self.u_min, v: self.v_min }
                } else {
                    DomainVertex { u: self.u_min, v: self.v_max }
                }
            }
            _ => {
                if self.vtx_idx == 0 {
                    DomainVertex { u: self.u_max, v: self.v_min }
                } else {
                    DomainVertex { u: self.u_max, v: self.v_max }
                }
            }
        }
    }
    /// OCCT TopolTool::NextVertex().
    pub fn next_vertex(&mut self) {
        self.vtx_idx += 1;
    }
    /// OCCT TopolTool::Identical(V1, V2).
    pub fn identical(&self, v1: DomainVertex, v2: DomainVertex) -> bool {
        (v1.u - v2.u).abs() <= rcad_kernel::precision::CONFUSION
            && (v1.v - v2.v).abs() <= rcad_kernel::precision::CONFUSION
    }
    /// OCCT IntPatch_HInterTool::Parameter(V, A).
    pub fn parameter(&self, v: DomainVertex, a: &Curve2d) -> f64 {
        if matches!(a, Curve2d::Line(l) if l.direction.y.abs() > 0.5) {
            v.v
        } else {
            v.u
        }
    }
    /// OCCT IntPatch_HInterTool::Tolerance(V, A) — vertex resolution on the arc.
    pub fn vertex_tolerance(&self, _v: DomainVertex, _a: &Curve2d) -> f64 {
        // rcad: no BRep vertex resolution; use the point-confusion tolerance.
        rcad_kernel::precision::CONFUSION
    }
    /// The index of the arc matching `a`.
    fn arc_of(&self, a: &Curve2d) -> usize {
        if let Curve2d::Line(l) = a {
            if l.direction.y.abs() > 0.5 {
                if (l.origin.y - self.v_min).abs() < (l.origin.y - self.v_max).abs() {
                    2
                } else {
                    3
                }
            } else if (l.origin.y - self.v_min).abs() < (l.origin.y - self.v_max).abs() {
                0
            } else {
                1
            }
        } else {
            0
        }
    }

    /// OCCT TopolTool::Edge() — the BRep edge; rcad's domain has no BRep edge.
    pub fn edge(&self) -> Option<()> {
        None
    }
}

// =====================================================================
// MinFunction — wraps ArcFunction as F -> F^2 (IntStart_SearchOnBoundaries.gxx L105-127)
// =====================================================================

struct MinFunction<'a> {
    func: &'a mut ArcFunction,
}

impl FunctionValue for MinFunction<'_> {
    fn value(&mut self, x: f64) -> Option<f64> {
        let v = self.func.value(x)?;
        Some(v * v)
    }
}

// =====================================================================
// SolInfo — sortable (parameter, index) pair (IntStart_SearchOnBoundaries.gxx L191-227)
// =====================================================================

#[derive(Clone, Copy)]
struct SolInfo {
    math_index: i32,
    value: f64,
}

impl SolInfo {
    fn new() -> Self {
        SolInfo {
            math_index: -1,
            value: f64::MAX,
        }
    }
    fn init_from_all_roots(&mut self, sol: &FunctionAllRoots, index: i32) {
        self.math_index = index;
        self.value = sol.get_point(index as usize);
    }
    fn init_from_int_cs(&mut self, sol: &IntCurveSurface, index: i32) {
        self.math_index = index;
        self.value = sol.point(index as usize).w();
    }
}

// =====================================================================
// IntPatch_TheSOnBounds
// =====================================================================

/// OCCT IntPatch_TheSOnBounds (= IntStart_SearchOnBoundaries).
pub struct SOnBounds {
    done: bool,
    all: bool,
    sseg: Vec<Segment>,
    spnt: Vec<PathPoint>,
}

impl SOnBounds {
    pub fn new() -> Self {
        SOnBounds {
            done: false,
            all: false,
            sseg: Vec::new(),
            spnt: Vec::new(),
        }
    }

    /// OCCT IntStart_SearchOnBoundaries::Perform (L1135-1232).
    pub fn perform(
        &mut self,
        func: &mut ArcFunction,
        domain: &mut Domain,
        tol_boundary: f64,
        tol_tangency: f64,
        recheck_on_regularity: bool,
    ) {
        self.done = false;
        self.spnt.clear();
        self.sseg.clear();

        domain.init();
        self.all = domain.more();

        while domain.more() {
            let a = domain.value();
            // IntPatch_HInterTool::HasBeenSeen always returns FALSE.
            if !has_been_seen(&a) {
                func.set_arc(a.clone());
                find_vertex(a.clone(), domain, func, &mut self.spnt, tol_boundary);
                let (mut p_deb, mut p_fin) = domain.bounds(&a);
                if rcad_kernel::precision::is_negative_infinite_value(p_deb)
                    || rcad_kernel::precision::is_positive_infinite_value(p_fin)
                {
                    let mut nb_echant = 0i32;
                    compute_bounds_from_infinite(func, &mut p_deb, &mut p_fin, &mut nb_echant);
                }
                let mut arcsol = false;
                bounded_arc(
                    a.clone(),
                    domain,
                    p_deb,
                    p_fin,
                    func,
                    &mut self.spnt,
                    &mut self.sseg,
                    tol_boundary,
                    tol_tangency,
                    &mut arcsol,
                    recheck_on_regularity,
                );
                self.all = self.all && arcsol;
            }
            domain.next();
        }
        self.done = true;
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT AllArcSolution().
    pub fn all_arc_solution(&self) -> bool {
        self.all
    }
    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.spnt.len()
    }
    /// OCCT Point(Index) — 1-based.
    pub fn point(&self, index: usize) -> &PathPoint {
        &self.spnt[index - 1]
    }
    /// OCCT NbSegments().
    pub fn nb_segments(&self) -> usize {
        self.sseg.len()
    }
    /// OCCT Segment(Index) — 1-based.
    pub fn segment(&self, index: usize) -> &Segment {
        &self.sseg[index - 1]
    }
}

/// OCCT IntPatch_HInterTool::HasBeenSeen (IntPatch_HInterTool.cxx L318-321).
fn has_been_seen(_a: &Curve2d) -> bool {
    false
}

// =====================================================================
// FindVertex (IntStart_SearchOnBoundaries.gxx L131-165)
// =====================================================================

fn find_vertex(
    a: Curve2d,
    domain: &mut Domain,
    func: &mut ArcFunction,
    pnt: &mut Vec<PathPoint>,
    toler: f64,
) {
    domain.initialize(&a);
    domain.init_vertex_iterator();
    while domain.more_vertex() {
        let vtx = domain.vertex();
        let param = domain.parameter(vtx, &a);
        if let Some(valf) = func.value(param) {
            if valf.abs() <= toler {
                let itemp = func.get_state_number();
                let mut pp = PathPoint::new();
                pp.set_value_with_vtx(func.valpoint(itemp), toler, vtx, a.clone(), param);
                pnt.push(pp);
            }
        }
        domain.next_vertex();
    }
}

// =====================================================================
// ComputeBoundsfromInfinite (IntStart_SearchOnBoundaries.gxx L804-876)
// =====================================================================

fn compute_bounds_from_infinite(
    func: &mut ArcFunction,
    p_deb: &mut f64,
    p_fin: &mut f64,
    _nb_echant: &mut i32,
) {
    let u0 = 0.0;
    let du = 0.001;
    let mut dist0;
    let mut dist1;
    let mut ok = func.value(u0);
    if ok.is_none() {
        return;
    }
    dist0 = ok.unwrap();
    ok = func.value(u0 + du);
    if ok.is_none() {
        return;
    }
    dist1 = ok.unwrap();
    let mut d_dist = dist1 - dist0;
    if d_dist != 0.0 {
        *p_deb = u0 - du * dist0 / d_dist;
        *p_fin = *p_deb;
        let mut u_min = *p_deb - 1e5;
        ok = func.value(u_min);
        if ok.is_none() {
            return;
        }
        dist0 = ok.unwrap();
        ok = func.value(u_min + du);
        if ok.is_none() {
            return;
        }
        dist1 = ok.unwrap();
        d_dist = dist1 - dist0;
        if d_dist != 0.0 {
            u_min -= du * dist0 / d_dist;
        } else {
            u_min -= 10.0;
        }
        let mut u_max = *p_deb + 1e8;
        ok = func.value(u_max);
        if ok.is_none() {
            return;
        }
        dist0 = ok.unwrap();
        ok = func.value(u_max + du);
        if ok.is_none() {
            return;
        }
        dist1 = ok.unwrap();
        d_dist = dist1 - dist0;
        if d_dist != 0.0 {
            u_max -= du * dist0 / d_dist;
        } else {
            u_max += 10.0;
        }
        if u_min > *p_deb {
            u_min = *p_deb - 10.0;
        }
        if u_max < *p_deb {
            u_max = *p_deb + 10.0;
        }
        *p_fin = u_max + 10.0 * (u_max - u_min);
        *p_deb = u_min - 10.0 * (u_max - u_min);
    } else {
        // Possibility of an arc entirely in the quadric.
        *p_deb = 1e10;
        *p_fin = -1e10;
    }
}

// =====================================================================
// BoundedArc (IntStart_SearchOnBoundaries.gxx L229-798)
// =====================================================================

#[allow(clippy::too_many_arguments)]
fn bounded_arc(
    a: Curve2d,
    domain: &mut Domain,
    p_deb: f64,
    p_fin: f64,
    func: &mut ArcFunction,
    pnt: &mut Vec<PathPoint>,
    seg: &mut Vec<Segment>,
    tol_boundary: f64,
    tol_tangency: f64,
    arcsol: &mut bool,
    recheck_on_regularity: bool,
) {
    let mut nbi = 0;
    let mut nbp = 0;

    // EpsX ~ 1e-5 and ResolutionU/V ~ 1e-9.
    let mut eps_x = 1.0e-10;

    let mut nb_echant = func.nb_samples();
    if nb_echant < 100 {
        nb_echant = 100;
    }

    // Adjust tolerances for short arcs.
    let mut n_tol_tangency = tol_tangency;
    if (p_fin - p_deb) < (tol_tangency * 10.0) {
        n_tol_tangency = (p_fin - p_deb) * 0.1;
    }
    if eps_x > (n_tol_tangency + n_tol_tangency) {
        eps_x = n_tol_tangency * 0.1;
    }

    let mut para = 0.0;
    let mut dist: f64;
    let mut maxdist: f64;

    // ---- Rejection test (L292-333).
    let mut rejection = true;
    let mut maxdr;
    let mut maxr;
    let mut minr;
    let mut ur;
    let mut dur;
    minr = f64::MAX;
    maxr = -minr;
    maxdr = -minr;
    dur = (p_fin - p_deb) * 0.2;
    let mut i = 1;
    ur = p_deb;
    while i <= 6 {
        if let Some((f, d)) = func.values(ur) {
            let mut lminr;
            let mut lmaxr;
            let mut d = d;
            if d < 0.0 {
                d = -d;
            }
            d *= dur + dur;
            if d > maxdr {
                maxdr = d;
            }
            lminr = f - d;
            lmaxr = f + d;
            if lminr < minr {
                minr = lminr;
            }
            if lmaxr > maxr {
                maxr = lmaxr;
            }
            if minr < 0.0 && maxr > 0.0 {
                rejection = false;
                break;
            }
        }
        ur += dur;
        i += 1;
    }
    if rejection {
        dur = 0.001 + maxdr + (maxr - minr) * 0.1;
        minr -= dur;
        maxr += dur;
        if minr < 0.0 && maxr > 0.0 {
            rejection = false;
        }
    }

    *arcsol = false;

    if !rejection {
        let quadric = func.quadric().clone();
        let type_quad = quadric.type_quadric();
        let mut type_con_s = CurveType3d::Other;

        // ---- Exact intersection (L339-472).
        let mut int_cs = IntCurveSurface::new();
        let mut is_int_cs_done = false;
        let mut params: Vec<f64> = Vec::new();
        let mut p_sol: Option<FunctionAllRoots> = None;

        let echant = FunctionSample::new(p_deb, p_fin, nb_echant);

        // maxdist: L354-372.
        let mut aelargir = true;
        maxdist = tol_boundary + tol_tangency;
        i = 1;
        while i <= nb_echant && aelargir {
            let u = echant.get_parameter(i);
            if let Some(d) = func.value(u) {
                if d > maxdist || -d > maxdist {
                    aelargir = false;
                }
            }
            i += 1;
        }
        if !(aelargir && maxdist < 0.01) {
            maxdist = tol_boundary;
        }

        if type_quad != crate::topalgo::int_surf::quadric::QuadricType::Other {
            // Build the 3D curve of the boundary arc on the surface.
            let surf = func.surface().clone();
            if let Some((curve3, ctype)) = curve_on_surface(&a, &surf) {
                type_con_s = ctype;
                // Exact solution only for canonic curves, non-Torus quadrics,
                // and non-degenerated entities.
                if is_canonic(type_con_s)
                    && type_quad != crate::topalgo::int_surf::quadric::QuadricType::Torus
                    && !is_degenerated_curve(&curve3)
                    && !is_degenerated_quadric(&quadric)
                {
                    let quad_surf = quadric_to_surface3(&quadric);
                    int_cs.perform(&curve3, &quad_surf);
                }
            }
            is_int_cs_done = int_cs.is_done();
            if is_int_cs_done {
                nbp = int_cs.nb_points();
                nbi = int_cs.nb_segments();
            }
            if nbp == 0 && nbi == 0 {
                is_int_cs_done = false;
            }
        }

        if !is_int_cs_done {
            let p = FunctionAllRoots::new(func, &echant, eps_x, maxdist, maxdist);
            assert!(p.is_done(), "FunctionAllRoots not done");
            p_sol = Some(p);
            nbp = p_sol.as_ref().unwrap().nb_points();
        }

        // ---- RecheckOnRegularity (L491-560) — rcad: the domain has no BRep
        // edge so IsRegularity returns false; branch is not entered.
        if recheck_on_regularity && nbp > 0 && is_regularity(domain) {
            // Ported structure; not reachable without a BRep edge.
        }

        // ---- Solution points (L568-738).
        if nbp != 0 {
            let mut a_si: Vec<SolInfo> = (1..=nbp).map(|_| SolInfo::new()).collect();
            for i in 1..=nbp {
                if is_int_cs_done {
                    a_si[i - 1].init_from_int_cs(&int_cs, i as i32);
                } else {
                    a_si[i - 1].init_from_all_roots(p_sol.as_ref().unwrap(), i as i32);
                }
            }
            a_si.sort_by(|x, y| x.value.partial_cmp(&y.value).unwrap_or(std::cmp::Ordering::Equal));

            // TreatLC: BRep-edge special handling.  rcad has no BRep edge, so
            // TreatLC returns 1 (use the old way), exactly like OCCT for a
            // null edge.
            let ip = treat_lc(&a, domain, &quadric, tol_boundary, pnt);
            if ip != 0 {
                // ---- Tangent/quasi-tangent treatment (L599-681).
                for i in 1..nbp {
                    let mut parap1 = a_si[i].value;
                    para = a_si[i - 1].value;

                    let mut param = (para + parap1) * 0.5;
                    let mut yf = 0.0;
                    let mut ym = 0.0;
                    let mut yl = 0.0;
                    if let Some(vm) = func.value(param) {
                        ym = vm;
                        if ym.abs() < maxdist {
                            let sm = ym.signum();
                            let mut a_tang = match func.value(para) {
                                Some(vf) => {
                                    yf = vf;
                                    match func.value(parap1) {
                                        Some(vl) => {
                                            yl = vl;
                                            true
                                        }
                                        None => false,
                                    }
                                }
                                None => false,
                            };
                            if a_tang {
                                a_tang = a_tang && yf.abs() < maxdist && yl.abs() < maxdist;
                            }
                            if a_tang && is_int_cs_done && type_con_s == CurveType3d::Line {
                                // Interval from exact intersection: tangent only if
                                // all points are on the same side.
                                let sf = yf.signum();
                                let sl = yl.signum();
                                a_tang = a_tang && (sm == sf) && (sm == sl);
                            }
                            if a_tang {
                                // Consider this interval as tangent.  Find the
                                // parameter with the lowest function value.
                                let mut a_tol = tol_boundary * 1000.0;
                                if a_tol > 0.001 {
                                    a_tol = 0.001;
                                }
                                parap1 = if parap1.abs() < 1.0e9 {
                                    parap1
                                } else if parap1 >= 0.0 {
                                    1.0e9
                                } else {
                                    -1.0e9
                                };
                                para = if para.abs() < 1.0e9 {
                                    para
                                } else if para >= 0.0 {
                                    1.0e9
                                } else {
                                    -1.0e9
                                };

                                let a_nb_nodes =
                                    ((parap1 - para) / a_tol).ceil() as i64;
                                let mut a_val = f64::MAX;
                                let mut a_val_max = 0.0;
                                let a_delta = (parap1 - para) / (a_nb_nodes + 1) as f64;
                                let mut ii = 0i64;
                                while ii <= a_nb_nodes + 1 {
                                    let a_cur_par = if ii < a_nb_nodes + 1 {
                                        para + ii as f64 * a_delta
                                    } else {
                                        parap1
                                    };
                                    if let Some(a_cur_val) = func.value(a_cur_par) {
                                        let an_abs_val = a_cur_val.abs();
                                        if an_abs_val < a_val {
                                            a_val = an_abs_val;
                                            param = a_cur_par;
                                        }
                                        if an_abs_val > a_val_max {
                                            a_val_max = an_abs_val;
                                        }
                                    }
                                    ii += 1;
                                }
                                if is_int_cs_done && a_nb_nodes > 1 {
                                    a_tang = (param - para).abs() > eps_x
                                        && (parap1 - param).abs() > eps_x
                                        && 0.01 * a_val_max <= a_val;
                                }
                                if a_tang {
                                    a_si[i - 1].value = p_deb - 1.0;
                                    a_si[i].value = param;
                                }
                            }
                        }
                    }
                }

                // ---- Emit points (L684-736).
                for i in 1..=nbp {
                    para = a_si[i - 1].value;
                    if (para - p_deb) < eps_x || (p_fin - para) < eps_x {
                        continue;
                    }
                    let d = match func.value(para) {
                        Some(v) => v,
                        None => continue,
                    };
                    let dist = d.abs();

                    let mut an_indx = -1;
                    let a_param = a_si[i - 1].value;
                    if dist < maxdist {
                        if !is_int_cs_done
                            && ((a_param - p_deb).abs() <= rcad_kernel::precision::PCONFUSION
                                || (a_param - p_fin).abs() <= rcad_kernel::precision::PCONFUSION)
                        {
                            an_indx = p_sol.as_ref().unwrap().get_point_state(a_si[i - 1].math_index as usize);
                        }
                    }

                    let mut a_pnt = if an_indx < 0 {
                        func.last_computed_point()
                    } else {
                        func.valpoint(an_indx)
                    };

                    if dist > 0.1 * rcad_kernel::precision::CONFUSION {
                        // Precise found points: make the vertex nearer to the
                        // intersection line and merge near vertices.
                        let a_f_par = if i == 1 {
                            p_deb
                        } else {
                            (para + a_si[i - 2].value) / 2.0
                        };
                        let a_l_par = if i == nbp {
                            p_fin
                        } else {
                            (para + a_si[i].value) / 2.0
                        };

                        let mut a_new_func = MinFunction { func };
                        let mut a_min = BrentMinimum::new(
                            rcad_kernel::precision::CONFUSION,
                            100,
                            1.0e-12,
                        );
                        a_min.perform(&mut a_new_func, a_f_par, para, a_l_par);
                        if a_min.is_done() {
                            para = a_min.location();
                            let a_p2d = a.point_at(para);
                            let surf = func.surface().clone();
                            a_pnt = surf.point_at(a_p2d.x, a_p2d.y);
                        }
                    }

                    let mut range = 0usize;
                    point_process(a_pnt, para, &a, domain, pnt, tol_boundary, &mut range);
                }
            } // end if(ip)
        } // end if(nbp)

        // ---- Restriction segments (L744-787).
        if !is_int_cs_done {
            nbi = p_sol.as_ref().unwrap().nb_intervals();
        }

        if !recheck_on_regularity && nbp != 0 {
            nbi = 0;
        }

        let mut pardeb = 0.0;
        let mut parfin = 0.0;
        let mut ptdeb = DVec3::ZERO;
        let mut ptfin = DVec3::ZERO;
        for i in 1..=nbi {
            let mut newseg = Segment {
                arc: a.clone(),
                has_fp: false,
                the_fp: PathPoint::new(),
                has_lp: false,
                the_lp: PathPoint::new(),
            };
            newseg.set_value(a.clone());
            if is_int_cs_done {
                let int_seg = int_cs.segment(i);
                let end1 = int_seg.first_point();
                let end2 = int_seg.second_point();
                pardeb = end1.w();
                parfin = end2.w();
                ptdeb = end1.pnt();
                ptfin = end2.pnt();
            } else {
                let (pd, pf) = p_sol.as_ref().unwrap().get_interval(i);
                pardeb = pd;
                parfin = pf;
                let (ideb, ifin) = p_sol.as_ref().unwrap().get_interval_state(i);
                ptdeb = func.valpoint(ideb);
                ptfin = func.valpoint(ifin);
            }

            let mut ranged = 0usize;
            point_process(ptdeb, pardeb, &a, domain, pnt, tol_boundary, &mut ranged);
            newseg.set_limit_point(pnt[ranged - 1].clone(), true);
            let mut rangef = 0usize;
            point_process(ptfin, parfin, &a, domain, pnt, tol_boundary, &mut rangef);
            newseg.set_limit_point(pnt[rangef - 1].clone(), false);
            seg.push(newseg);
        }

        if nbi == 1 {
            if (pardeb - p_deb).abs() < rcad_kernel::precision::PCONFUSION
                && (parfin - p_fin).abs() < rcad_kernel::precision::PCONFUSION
            {
                *arcsol = true;
            }
        }
    }
}

// =====================================================================
// PointProcess (IntStart_SearchOnBoundaries.gxx L880-1001)
// =====================================================================

fn point_process(
    pt: DVec3,
    para: f64,
    a: &Curve2d,
    domain: &mut Domain,
    pnt: &mut Vec<PathPoint>,
    tol: f64,
    range: &mut usize,
) {
    let nbsol = pnt.len();
    let mut found = false;
    let mut goon;

    domain.initialize(a);
    domain.init_vertex_iterator();
    found = false;
    goon = domain.more_vertex();
    while goon {
        let vtx = domain.vertex();
        let dist = (para - domain.parameter(vtx, a)).abs();
        let toler = domain.vertex_tolerance(vtx, a);
        if dist <= toler {
            // Locate the vertex in the list of solutions.
            let mut k = 1usize;
            found = k > nbsol;
            while !found {
                let ptsol = &pnt[k - 1];
                if !ptsol.is_new() {
                    if domain.identical(ptsol.vertex(), vtx)
                        && curves_same(ptsol.arc(), a)
                        && (ptsol.parameter() - para).abs() <= toler
                    {
                        found = true;
                    } else {
                        k += 1;
                        found = k > nbsol;
                    }
                } else {
                    k += 1;
                    found = k > nbsol;
                }
            }
            if k <= nbsol {
                *range = k;
            } else {
                let mut ptsol = PathPoint::new();
                ptsol.set_value_with_vtx(pt, tol, vtx, a.clone(), para);
                pnt.push(ptsol);
                *range = pnt.len();
            }
            found = true;
            goon = false;
        } else {
            domain.next_vertex();
            goon = domain.more_vertex();
        }
    }

    if !found {
        // No vertex matched: do not add segment's extremities if they already exist.
        let mut found_internal = false;
        for k in 1..=pnt.len() {
            let ptsol = &pnt[k - 1];
            if !curves_same(ptsol.arc(), a) || !ptsol.is_new() {
                continue;
            }
            if (ptsol.parameter() - para).abs() <= rcad_kernel::precision::PCONFUSION {
                found_internal = true;
                *range = k;
            }
        }
        if !found_internal {
            let mut t = tol;
            t *= 1000.0;
            if t > 0.005 {
                t = 0.005;
            }
            let mut ptsol = PathPoint::new();
            ptsol.set_value(pt, t, a.clone(), para);
            pnt.push(ptsol);
            *range = pnt.len();
        }
    }
}

// =====================================================================
// IsRegularity (IntStart_SearchOnBoundaries.gxx L1005-1016)
// =====================================================================

fn is_regularity(domain: &Domain) -> bool {
    // OCCT: BRep_Tool::HasContinuity(*anE) where anE = aDomain->Edge().
    // rcad's domain has no BRep edge -> false.
    domain.edge().is_some()
}

// =====================================================================
// TreatLC (IntStart_SearchOnBoundaries.gxx L1020-1123)
// =====================================================================

fn treat_lc(
    a: &Curve2d,
    domain: &Domain,
    _quadric: &Quadric,
    _tol_boundary: f64,
    _pnt: &mut Vec<PathPoint>,
) -> i32 {
    // OCCT: if (aDomain->Edge() == nullptr) return 1;  — rcad's domain has no
    // BRep edge, so the special line-cylinder tangent treatment never applies.
    let _ = a;
    if domain.edge().is_none() {
        return 1;
    }
    1
}

// =====================================================================
// 3D curve of a boundary arc on the surface (Adaptor3d_CurveOnSurface)
// =====================================================================

/// OCCT GeomAbs_CurveType for the 3D curve of the boundary arc.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveType3d {
    Line,
    Circle,
    Other,
}

fn is_canonic(t: CurveType3d) -> bool {
    matches!(t, CurveType3d::Line | CurveType3d::Circle)
}

/// OCCT IsDegenerated(Handle(Adaptor3d_CurveOnSurface)) (L167-176):
/// a Circle with radius <= Confusion is degenerated.
fn is_degenerated_curve(curve: &rcad_kernel::geom::Curve3) -> bool {
    if let rcad_kernel::geom::Curve3::Circle(c) = curve {
        if c.radius <= rcad_kernel::precision::CONFUSION {
            return true;
        }
    }
    false
}

/// OCCT IsDegenerated(const IntSurf_Quadric&) (L178-189): a Cone with
/// |SemiAngle| < 0.02 or > 1.55 is degenerated.
fn is_degenerated_quadric(quadric: &Quadric) -> bool {
    if quadric.type_quadric() == crate::topalgo::int_surf::quadric::QuadricType::Cone {
        let a = quadric.semi_angle().abs();
        if a < 0.02 || a > 1.55 {
            return true;
        }
    }
    false
}

/// Build the 3D curve (parameterized by the 2D arc parameter) of a boundary
/// arc on a surface.  The boundary arcs of the FF domain are lines in UV
/// space; on an analytic quadric their 3D image is a Line or a Circle.
fn curve_on_surface(
    a: &Curve2d,
    surf: &Surface3,
) -> Option<(rcad_kernel::geom::Curve3, CurveType3d)> {
    let arc_along_u = matches!(a, Curve2d::Line(l) if l.direction.y.abs() <= 0.5);
    match surf {
        Surface3::Plane(p) => {
            // P(t) = p.origin + u(t)*p.u_dir + v(t)*p.v_dir.
            let (u0, v0) = a.point_at(0.0).into();
            let o = p.origin + u0 * p.u_dir + v0 * p.v_dir;
            let dir = if arc_along_u { p.u_dir } else { p.v_dir };
            Some((
                rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 { origin: o, direction: dir }),
                CurveType3d::Line,
            ))
        }
        Surface3::Cylinder(c) => {
            let z = c.axis.normalize_or_zero();
            let x = c.ref_dir.normalize_or_zero();
            let y = z.cross(x).normalize_or_zero();
            if arc_along_u {
                // V = const -> circle at height V.
                let (u0, v0) = a.point_at(0.0).into();
                let center = c.origin + v0 * z;
                let radius = c.radius;
                Some((
                    rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3 {
                        center,
                        normal: z,
                        x_dir: x,
                        y_dir: y,
                        radius,
                    }),
                    CurveType3d::Circle,
                ))
            } else {
                // U = const -> generatrix line.
                let (u0, v0) = a.point_at(0.0).into();
                let o = c.origin + c.radius * (u0.cos() * x + u0.sin() * y) + v0 * z;
                Some((
                    rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 { origin: o, direction: z }),
                    CurveType3d::Line,
                ))
            }
        }
        Surface3::Cone(c) => {
            let z = c.axis.normalize_or_zero();
            let x = rcad_kernel::geom::any_perpendicular(z).normalize_or_zero();
            let y = z.cross(x).normalize_or_zero();
            let apex = c.apex_point();
            let semi = c.half_angle_rad;
            let _ = y;
            if arc_along_u {
                // V = const -> circle at height V.
                let (u0, v0) = a.point_at(0.0).into();
                let r = c.radius + v0 * semi.sin();
                let center = apex + v0 * semi.cos() * z;
                Some((
                    rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3 {
                        center,
                        normal: z,
                        x_dir: x,
                        y_dir: y,
                        radius: r,
                    }),
                    CurveType3d::Circle,
                ))
            } else {
                // U = const -> generatrix line.
                let (u0, v0) = a.point_at(0.0).into();
                let o = apex + (c.radius + v0 * semi.sin()) * (u0.cos() * x + u0.sin() * y)
                    + v0 * semi.cos() * z;
                let dir = semi.sin() * (u0.cos() * x + u0.sin() * y) + semi.cos() * z;
                Some((
                    rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 { origin: o, direction: dir }),
                    CurveType3d::Line,
                ))
            }
        }
        Surface3::Sphere(s) => {
            let z = s.axis.normalize_or_zero();
            let x = s.ref_dir.normalize_or_zero();
            let y = z.cross(x).normalize_or_zero();
            let r = s.radius;
            let (u0, v0) = a.point_at(0.0).into();
            if arc_along_u {
                // V = const -> parallel circle.
                let center = s.center + r * v0.sin() * z;
                let radius = (r * v0.cos()).abs();
                Some((
                    rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3 {
                        center,
                        normal: z,
                        x_dir: x,
                        y_dir: y,
                        radius,
                    }),
                    CurveType3d::Circle,
                ))
            } else {
                // U = const -> meridian circle in the plane spanned by (x,y) at U.
                let ux = u0.cos() * x + u0.sin() * y;
                Some((
                    rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3 {
                        center: s.center,
                        normal: ux,
                        x_dir: z,
                        y_dir: ux.cross(z).normalize_or_zero(),
                        radius: r,
                    }),
                    CurveType3d::Circle,
                ))
            }
        }
        Surface3::Torus(t) => {
            let z = t.axis.normalize_or_zero();
            let x = rcad_kernel::geom::any_perpendicular(z).normalize_or_zero();
            let y = z.cross(x).normalize_or_zero();
            let (u0, v0) = a.point_at(0.0).into();
            if arc_along_u {
                // V = const -> circle of radius R + r*cos(v).
                let center = t.center + (t.minor_radius * v0.sin()) * z;
                let radius = t.major_radius + t.minor_radius * v0.cos();
                Some((
                    rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3 {
                        center,
                        normal: z,
                        x_dir: x,
                        y_dir: y,
                        radius,
                    }),
                    CurveType3d::Circle,
                ))
            } else {
                // U = const -> circle of radius r in the plane at U.
                let ux = u0.cos() * x + u0.sin() * y;
                let center = t.center + t.major_radius * ux;
                Some((
                    rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3 {
                        center,
                        normal: ux,
                        x_dir: z,
                        y_dir: ux.cross(z).normalize_or_zero(),
                        radius: t.minor_radius,
                    }),
                    CurveType3d::Circle,
                ))
            }
        }
        _ => None,
    }
}

/// Reconstruct a Surface3 from a Quadric.
fn quadric_to_surface3(quad: &Quadric) -> Surface3 {
    use crate::topalgo::int_surf::quadric::QuadricType;
    match quad.type_quadric() {
        QuadricType::Plane => Surface3::Plane(quad.plane()),
        QuadricType::Cylinder => Surface3::Cylinder(quad.cylinder()),
        QuadricType::Sphere => Surface3::Sphere(quad.sphere()),
        QuadricType::Cone => Surface3::Cone(quad.cone()),
        QuadricType::Torus => Surface3::Torus(quad.torus()),
        QuadricType::Other => Surface3::Plane(quad.plane()),
    }
}

/// rcad: OCCT compares handles (Arc() == A) by pointer identity.  The rcad
/// boundary arcs are value types; equality is checked structurally.
pub fn curves_same(a: &Curve2d, b: &Curve2d) -> bool {
    if let (Curve2d::Line(l1), Curve2d::Line(l2)) = (a, b) {
        (l1.origin.x - l2.origin.x).abs() < rcad_kernel::precision::CONFUSION
            && (l1.origin.y - l2.origin.y).abs() < rcad_kernel::precision::CONFUSION
            && (l1.direction.x - l2.direction.x).abs() < rcad_kernel::precision::CONFUSION
            && (l1.direction.y - l2.direction.y).abs() < rcad_kernel::precision::CONFUSION
    } else {
        false
    }
}

