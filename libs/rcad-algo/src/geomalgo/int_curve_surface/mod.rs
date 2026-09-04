//! IntCurveSurface — the HInter assembly: intersection of a 3D curve with a
//! parametrised surface through adaptors (OCCT TKGeomAlgo IntCurveSurface
//! package).
//!
//! 1:1 translations:
//! - [`TransitionOnCurve`] — IntCurveSurface_TransitionOnCurve.hxx (L37-42).
//! - [`IntersectionPoint`] — IntCurveSurface_IntersectionPoint.hxx/.cxx/.lxx.
//! - [`IntersectionSegment`] — IntCurveSurface_IntersectionSegment.hxx/.cxx.
//! - [`Intersection`] — IntCurveSurface_Intersection.hxx/.cxx (the base of
//!   [`HInter`], with the PARAMEQUAL(1e-8) append de-duplication).
//! - [`HCurveTool`] — IntCurveSurface_TheHCurveTool.hxx (L41-191) + .cxx
//!   (NbSamples L34-70, SamplePars L72-299) as a static tool trait over the
//!   curve adaptor type `C` (the OCCT template parameter
//!   `TheCurveTool = IntCurveSurface_TheHCurveTool`).
//! - [`HSurfaceTool`] — Adaptor3d_HSurfaceTool.hxx (L40-293) + .cxx
//!   (NbSamplesU/V L27-115) as a static tool trait over the surface adaptor
//!   type `S` (`TheSurfaceTool = Adaptor3d_HSurfaceTool`).  UTrim/VTrim need
//!   Adaptor3d_TopolTool (TKHLR runway 2a-4) and stay documented
//!   unimplemented.
//! - [`inter_utils`] — IntCurveSurface_InterUtils.pxx (L1-1637).
//! - [`inter_impl`] — IntCurveSurface_Inter.pxx (L1-1006).
//! - [`quad_curv_exact`] — IntCurveSurface_TheQuadCurvExactHInter +
//!   TheQuadCurvFuncOfTheQuadCurvExactHInter + QuadricCurveExactInterUtils.pxx.
//! - [`hinter`] — IntCurveSurface_HInter.hxx/.cxx (L1-581): the thin shell
//!   delegating to inter_impl / inter_utils.
//!
//! The OCCT callback lambdas passed from HInter.cxx into the Inter.pxx
//! templates (ResetFunc / PerformBoundsFunc / AppendFunc / ...) are encoded
//! as the [`HInterHost`] trait — one method per OCCT member function the
//! lambdas invoke — implemented by [`HInter`].

pub mod hinter;
pub mod inter_impl;
pub mod inter_utils;
pub mod quad_curv_exact;

pub use hinter::HInter;
pub use inter_utils::{SortedStartPoints, UVBounds};

use rcad_kernel::geom::{
    Circle3, ConicalSurface, CylindricalSurface, Ellipse3, Hyperbola3, Line3, Parabola3, Plane,
    SphericalSurface, ToroidalSurface,
};
use rcad_kernel::geom::CurveEval;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::bnd::BoundSortBox;

use crate::geomalgo::int_patch::int_conic_quad::IntConicQuad;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};

/// OCCT IntCurveSurface_TransitionOnCurve (TransitionOnCurve.hxx L37-42).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionOnCurve {
    Tangent,
    In,
    Out,
}

/// OCCT GeomAbs_SurfaceType — the values distinguished by the assembly; the
/// enum itself lives with the IntPatch translation (same OCCT enumeration).
pub type SurfaceType = crate::geomalgo::int_patch::GeomAbsSurfaceType;

/// OCCT IntCurveSurface_IntersectionPoint (IntersectionPoint.cxx L20-41 +
/// .lxx L19-42) — one curve/surface intersection point.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionPoint {
    my_p: rcad_kernel::geom::Point3,
    my_u_surf: f64,
    my_v_surf: f64,
    my_u_curv: f64,
    my_tr_on_curv: TransitionOnCurve,
}

impl IntersectionPoint {
    /// OCCT IntersectionPoint() (cxx L20-26) — empty.
    pub fn new() -> Self {
        IntersectionPoint {
            my_p: rcad_kernel::geom::Point3::ZERO,
            my_u_surf: 0.0,
            my_v_surf: 0.0,
            my_u_curv: 0.0,
            my_tr_on_curv: TransitionOnCurve::Tangent,
        }
    }

    /// OCCT IntersectionPoint(P, USurf, VSurf, UCurv, TrOnCurv) (cxx L29-41).
    pub fn with_values(
        p: rcad_kernel::geom::Point3,
        u_surf: f64,
        v_surf: f64,
        u_curv: f64,
        tr_on_curv: TransitionOnCurve,
    ) -> Self {
        IntersectionPoint {
            my_p: p,
            my_u_surf: u_surf,
            my_v_surf: v_surf,
            my_u_curv: u_curv,
            my_tr_on_curv: tr_on_curv,
        }
    }

    /// OCCT SetValues(P, USurf, VSurf, UCurv, TrOnCurv) (cxx L44-55).
    pub fn set_values(
        &mut self,
        p: rcad_kernel::geom::Point3,
        u_surf: f64,
        v_surf: f64,
        u_curv: f64,
        tr_on_curv: TransitionOnCurve,
    ) {
        self.my_p = p;
        self.my_u_surf = u_surf;
        self.my_v_surf = v_surf;
        self.my_u_curv = u_curv;
        self.my_tr_on_curv = tr_on_curv;
    }

    /// OCCT Values(P, USurf, VSurf, UCurv, TrOnCurv) (cxx L58-69).
    pub fn values(
        &self,
    ) -> (
        rcad_kernel::geom::Point3,
        f64,
        f64,
        f64,
        TransitionOnCurve,
    ) {
        (
            self.my_p,
            self.my_u_surf,
            self.my_v_surf,
            self.my_u_curv,
            self.my_tr_on_curv,
        )
    }

    /// OCCT Pnt() (.lxx L19-22) — the geometric point.
    pub fn pnt(&self) -> rcad_kernel::geom::Point3 {
        self.my_p
    }

    /// OCCT U() (.lxx L24-27) — the U parameter on the surface.
    pub fn u(&self) -> f64 {
        self.my_u_surf
    }

    /// OCCT V() (.lxx L29-32) — the V parameter on the surface.
    pub fn v(&self) -> f64 {
        self.my_v_surf
    }

    /// OCCT W() (.lxx L34-37) — the parameter on the curve.
    pub fn w(&self) -> f64 {
        self.my_u_curv
    }

    /// OCCT Transition() (.lxx L39-42).
    pub fn transition(&self) -> TransitionOnCurve {
        self.my_tr_on_curv
    }
}

impl Default for IntersectionPoint {
    fn default() -> Self {
        IntersectionPoint::new()
    }
}

/// OCCT IntCurveSurface_IntersectionSegment (IntersectionSegment.cxx L19-77).
#[derive(Debug, Clone, Copy)]
pub struct IntersectionSegment {
    my_p1: IntersectionPoint,
    my_p2: IntersectionPoint,
}

impl IntersectionSegment {
    /// OCCT IntersectionSegment() (cxx L19).
    pub fn new() -> Self {
        IntersectionSegment {
            my_p1: IntersectionPoint::new(),
            my_p2: IntersectionPoint::new(),
        }
    }

    /// OCCT IntersectionSegment(P1, P2) (cxx L22-28).
    pub fn with_values(p1: IntersectionPoint, p2: IntersectionPoint) -> Self {
        IntersectionSegment { my_p1: p1, my_p2: p2 }
    }

    /// OCCT SetValues(P1, P2) (cxx L31-36).
    pub fn set_values(&mut self, p1: IntersectionPoint, p2: IntersectionPoint) {
        self.my_p1 = p1;
        self.my_p2 = p2;
    }

    /// OCCT Values(P1, P2) (cxx L39-44).
    pub fn values(&self) -> (IntersectionPoint, IntersectionPoint) {
        (self.my_p1, self.my_p2)
    }

    /// OCCT FirstPoint(P1) (cxx L47-50).
    pub fn first_point_out(&self) -> IntersectionPoint {
        self.my_p1
    }

    /// OCCT SecondPoint(P2) (cxx L53-56).
    pub fn second_point_out(&self) -> IntersectionPoint {
        self.my_p2
    }

    /// OCCT FirstPoint() (cxx L59-62).
    pub fn first_point(&self) -> &IntersectionPoint {
        &self.my_p1
    }

    /// OCCT SecondPoint() (cxx L65-68).
    pub fn second_point(&self) -> &IntersectionPoint {
        &self.my_p2
    }
}

impl Default for IntersectionSegment {
    fn default() -> Self {
        IntersectionSegment::new()
    }
}

/// OCCT Intersection.cxx L23: `#define PARAMEQUAL(a, b) (std::abs((a) - (b)) < (1e-8))`.
fn paramequal(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-8
}

/// OCCT IntCurveSurface_Intersection (Intersection.hxx L31-109 + .cxx L26-203)
/// — the base of HInter: the done/parallel flags and the point/segment
/// sequences.
#[derive(Debug, Clone, Default)]
pub struct Intersection {
    pub(crate) done: bool,
    /// Curve is "parallel" surface (hxx L102-104).
    pub(crate) my_is_parallel: bool,
    pub(crate) lpnt: Vec<IntersectionPoint>,
    pub(crate) lseg: Vec<IntersectionSegment>,
}

impl Intersection {
    /// OCCT IntCurveSurface_Intersection() (cxx L26-30).
    pub fn new() -> Self {
        Intersection {
            done: false,
            my_is_parallel: false,
            lpnt: Vec::new(),
            lseg: Vec::new(),
        }
    }

    /// OCCT IsDone() (cxx L33-36).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsParallel() (cxx L39-42).
    pub fn is_parallel(&self) -> bool {
        self.my_is_parallel
    }

    /// OCCT NbPoints() (cxx L45-52).
    pub fn nb_points(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        self.lpnt.len()
    }

    /// OCCT NbSegments() (cxx L55-62).
    pub fn nb_segments(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        self.lseg.len()
    }

    /// OCCT Point(N) (cxx L65-72) — 1-based.
    pub fn point(&self, n: usize) -> &IntersectionPoint {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        &self.lpnt[n - 1]
    }

    /// OCCT Segment(N) (cxx L75-82) — 1-based.
    pub fn segment(&self, n: usize) -> &IntersectionSegment {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        &self.lseg[n - 1]
    }

    /// OCCT SetValues(Inter) (cxx L85-110).
    pub fn set_values(&mut self, other: &Intersection) {
        if other.done {
            self.lseg.clear();
            self.lpnt.clear();
            for p in &other.lpnt {
                self.lpnt.push(*p);
            }
            for s in &other.lseg {
                self.lseg.push(*s);
            }
            self.done = true;
            self.my_is_parallel = other.my_is_parallel;
        } else {
            self.done = false;
            self.my_is_parallel = false;
        }
    }

    /// OCCT Append(Inter, FirstParamOnCurve, LastParamOnCurve) (cxx L113-133)
    /// — the parameter arguments are unused in OCCT (commented out there).
    pub fn append_inter(&mut self, other: &Intersection, _first: f64, _last: f64) {
        if other.done {
            for i in 1..=other.lpnt.len() {
                self.append(other.point(i));
            }
            for i in 1..=other.lseg.len() {
                self.append_seg(other.segment(i));
            }
        }
    }

    /// OCCT Append(IntersectionPoint) (cxx L136-162) — de-duplicates against
    /// the already-stored points with the PARAMEQUAL(1e-8) test on (u, v, w,
    /// transition).
    pub fn append(&mut self, other_point: &IntersectionPoint) {
        let (p, u, v, w, tr_on_curve) = other_point.values();
        for existing in &self.lpnt {
            let (an_p, an_u, an_v, an_w, an_tr_on_curve) = existing.values();
            let _ = an_p;
            let _ = p;
            if paramequal(u, an_u)
                && paramequal(v, an_v)
                && paramequal(w, an_w)
                && an_tr_on_curve == tr_on_curve
            {
                return;
            }
        }
        self.lpnt.push(*other_point);
    }

    /// OCCT Append(IntersectionSegment) (cxx L165-168).
    pub fn append_seg(&mut self, other_segment: &IntersectionSegment) {
        self.lseg.push(*other_segment);
    }

    /// OCCT ResetFields() (cxx L171-180).
    pub fn reset_fields(&mut self) {
        if self.done {
            self.lseg.clear();
            self.lpnt.clear();
            self.done = false;
            self.my_is_parallel = false;
        }
    }
}

/// OCCT IntCurveSurface_TheHCurveTool (TheHCurveTool.hxx L41-191 + .cxx
/// L34-299) — the static tool over the curve adaptor `C`, i.e. the OCCT
/// template parameter `TheCurveTool`.  Every accessor forwards to the
/// adaptor, exactly as the OCCT static class forwards to
/// `Adaptor3d_Curve`.  NbSamples / SamplePars carry their OCCT default
/// implementations (TheHCurveTool.cxx).
pub trait HCurveTool {
    type Curve: ?Sized;

    /// OCCT FirstParameter(C) (hxx L46-49).
    fn first_parameter(c: &Self::Curve) -> f64;
    /// OCCT LastParameter(C) (hxx L51).
    fn last_parameter(c: &Self::Curve) -> f64;
    /// OCCT Continuity(C) (hxx L53).
    fn continuity(c: &Self::Curve) -> GeomAbsShape;
    /// OCCT NbIntervals(C, S) (hxx L57-60).
    fn nb_intervals(c: &Self::Curve, s: GeomAbsShape) -> usize;
    /// OCCT Intervals(C, T, S) (hxx L67-72) — fills T with the parameters
    /// bounding the continuity intervals (T is 0-based here; length is
    /// NbIntervals + 1).
    fn intervals(c: &Self::Curve, t: &mut [f64], s: GeomAbsShape);
    /// OCCT IsClosed(C) (hxx L74).
    fn is_closed(c: &Self::Curve) -> bool;
    /// OCCT IsPeriodic(C) (hxx L76).
    fn is_periodic(c: &Self::Curve) -> bool;
    /// OCCT Period(C) (hxx L78).
    fn period(c: &Self::Curve) -> f64;
    /// OCCT Value(C, U) (hxx L81).
    fn value(c: &Self::Curve, u: f64) -> rcad_kernel::geom::Point3;
    /// OCCT D0(C, U, P) (hxx L84).
    fn d0(c: &Self::Curve, u: f64) -> rcad_kernel::geom::Point3;
    /// OCCT D1(C, U, P, V) (hxx L90-93) — P is the curve point itself (the
    /// OCCT contract).
    fn d1(c: &Self::Curve, u: f64) -> (rcad_kernel::geom::Point3, rcad_kernel::geom::Vec3);
    /// OCCT D2(C, U, P, V1, V2) (hxx L99-106).
    fn d2(
        c: &Self::Curve,
        u: f64,
    ) -> (
        rcad_kernel::geom::Point3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
    );
    /// OCCT D3(C, U, P, V1, V2, V3) (hxx L112-120).
    fn d3(
        c: &Self::Curve,
        u: f64,
    ) -> (
        rcad_kernel::geom::Point3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
    );
    /// OCCT DN(C, U, N) (hxx L127-130).
    fn dn(c: &Self::Curve, u: f64, n: usize) -> rcad_kernel::geom::Vec3;
    /// OCCT Resolution(C, R3d) (hxx L134-137).
    fn resolution(c: &Self::Curve, r3d: f64) -> f64;
    /// OCCT GetType(C) (hxx L142).
    fn get_type(c: &Self::Curve) -> rcad_kernel::base::proj_lib::CurveType;
    /// OCCT Line(C) (hxx L144).
    fn line(c: &Self::Curve) -> Line3;
    /// OCCT Circle(C) (hxx L146).
    fn circle(c: &Self::Curve) -> Circle3;
    /// OCCT Ellipse(C) (hxx L148).
    fn ellipse(c: &Self::Curve) -> Ellipse3;
    /// OCCT Hyperbola(C) (hxx L150).
    fn hyperbola(c: &Self::Curve) -> Hyperbola3;
    /// OCCT Parabola(C) (hxx L152).
    fn parabola(c: &Self::Curve) -> Parabola3;
    /// OCCT Bezier(C) (hxx L154-157).
    fn bezier(c: &Self::Curve) -> &rcad_kernel::geom::BezierCurve3;
    /// OCCT BSpline(C) (hxx L159-162).
    fn bspline(c: &Self::Curve) -> &rcad_kernel::geom::BSplineCurve3;

    /// OCCT NbSamples(C, U0, U1) (TheHCurveTool.cxx L34-70).
    fn nb_samples(c: &Self::Curve, u0: f64, u1: f64) -> usize {
        let typ_c = Self::get_type(c);
        const NBS_OTHER: f64 = 10.0;
        let mut nbs = NBS_OTHER;

        if typ_c == rcad_kernel::base::proj_lib::CurveType::Line {
            nbs = 2.0;
        } else if typ_c == rcad_kernel::base::proj_lib::CurveType::Bezier {
            nbs = 3.0 + Self::bezier(c).control_points.len() as f64;
        } else if typ_c == rcad_kernel::base::proj_lib::CurveType::BSpline {
            nbs = distinct_knots(&Self::bspline(c).knots).len() as f64;
            nbs *= Self::bspline(c).degree as f64;
            nbs *= Self::last_parameter(c) - Self::first_parameter(c);
            let a_range = u1 - u0;
            if a_range.abs() > rcad_kernel::precision::PCONFUSION {
                nbs /= a_range;
            }
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        if nbs > 50.0 {
            nbs = 50.0;
        }
        nbs as usize
    }

    /// OCCT SamplePars(C, U0, U1, Defl, NbMin) (TheHCurveTool.cxx L72-299) —
    /// sample parameters within [U0, U1]; the BSpline branch keeps the knot
    /// interval subdivision and the deflection analysis.
    fn sample_pars(c: &Self::Curve, u0: f64, u1: f64, defl: f64, nb_min: usize) -> Vec<f64> {
        use rcad_kernel::base::proj_lib::CurveType;
        let typ_c = Self::get_type(c);
        const NBS_OTHER: f64 = 10.0;
        let mut nbs = NBS_OTHER;

        if typ_c == CurveType::Line {
            nbs = 2.0;
        } else if typ_c == CurveType::Bezier {
            nbs = 3.0 + Self::bezier(c).control_points.len() as f64;
        }

        if typ_c != CurveType::BSpline {
            if nbs > 50.0 {
                nbs = 50.0;
            }
            let nnbs = nbs as usize;

            let mut pars = vec![0.0f64; nnbs];
            let du = (u1 - u0) / (nnbs - 1) as f64;

            pars[0] = u0;
            pars[nnbs - 1] = u1;
            let mut u = u0 + du;
            for p in pars.iter_mut().skip(1).take(nnbs - 2) {
                *p = u;
                u += du;
            }
            return pars;
        }

        let a_bc = Self::bspline(c);
        let knots = distinct_knots(&a_bc.knots);
        // OCCT ui1 = FirstUKnotIndex(), ui2 = LastUKnotIndex() (1-based over
        // the distinct knots).
        let mut ui1 = 1usize;
        let mut ui2 = knots.len();
        for i in ui1..ui2 {
            if u0 >= knots[i - 1] && u0 < knots[i] {
                ui1 = i;
                break;
            }
        }
        for i in (ui1 + 1..=ui2).rev() {
            if u1 <= knots[i - 1] && u1 > knots[i - 2] {
                ui2 = i;
                break;
            }
        }

        let mut nbsu = ui2 - ui1 + 1;
        nbsu += (nbsu - 1) * (a_bc.degree.saturating_sub(1));
        let mut b_uniform = false;
        if nbsu < nb_min {
            nbsu = nb_min;
            b_uniform = true;
        }

        let mut a_pars = vec![0.0f64; nbsu];
        let mut a_flg = vec![false; nbsu];
        // Filling of sample parameters.
        if b_uniform {
            let mut t1 = u0;
            let t2 = u1;
            let dt = (t2 - t1) / (nbsu - 1) as f64;
            a_pars[0] = t1;
            a_flg[0] = false;
            a_pars[nbsu - 1] = t2;
            a_flg[nbsu - 1] = false;
            let mut tt = t1 + dt;
            for p in a_pars.iter_mut().skip(1).take(nbsu - 2) {
                *p = tt;
                tt += dt;
            }
        } else {
            let nbi = a_bc.degree;
            let mut k = 0usize;
            let mut t1 = u0;
            for i in ui1 + 1..=ui2 {
                let t2 = if i == ui2 { u1 } else { knots[i - 1] };
                let dt = (t2 - t1) / nbi as f64;
                let mut j = 1usize;
                loop {
                    k += 1;
                    a_pars[k - 1] = t1;
                    a_flg[k - 1] = false;
                    t1 += dt;
                    j += 1;
                    if j > nbi {
                        break;
                    }
                }
                t1 = t2;
            }
            k += 1;
            a_pars[k - 1] = t1;
        }
        // Analysis of deflection.

        let a_defl2 = (defl * defl).max(1.0e-9);
        let tol = (0.01 * a_defl2).max(1.0e-9);
        let mut l;

        let mut nb_samples = 2usize;
        a_flg[0] = true;
        a_flg[nbsu - 1] = true;
        let mut j = 1usize;
        let mut b_cont = true;
        while j < nbsu - 1 && b_cont {
            if a_flg[j] {
                j += 1;
                continue;
            }

            let mut t2 = a_pars[j - 1];
            let p1 = a_bc.point_at(t2);
            let mut broke = false;
            for k in j + 2..=nbsu {
                t2 = a_pars[k - 1];
                let p2 = a_bc.point_at(t2);

                if p1.distance_squared(p2) <= tol {
                    continue;
                }

                // gce_MakeLin(p1, p2) — the line through p1 and p2; its
                // gp_Lin::SquareDistance(pp) = |(pp-p1) x d|^2 with |d| = 1.
                let d = (p2 - p1).normalize();
                let mut ok = true;
                for l_i in j + 1..k {
                    if a_flg[l_i - 1] {
                        ok = false;
                        break;
                    }
                    let pp = a_bc.point_at(a_pars[l_i - 1]);
                    let cross = d.cross(pp - p1);
                    let dist = cross.length_squared();
                    if dist <= a_defl2 {
                        continue;
                    }
                    ok = false;
                    break;
                }
                l = k;

                if !ok {
                    j = l - 1;
                    a_flg[j - 1] = true;
                    nb_samples += 1;
                    broke = true;
                    break;
                }

                if a_flg[k - 1] {
                    j = k;
                    broke = true;
                    break;
                }
            }

            if !broke {
                b_cont = false;
            }
        }

        const MY_MIN_PNTS: usize = 5;
        if nb_samples < MY_MIN_PNTS {
            // Uniform distribution.
            let nb = MY_MIN_PNTS;
            let mut pars = vec![0.0f64; nb];
            let mut t1 = u0;
            let t2 = u1;
            let dt = (t2 - t1) / (nb - 1) as f64;
            pars[0] = t1;
            pars[nb - 1] = t2;
            for p in pars.iter_mut().skip(1).take(nb - 2) {
                t1 += dt;
                *p = t1;
            }
            return pars;
        }

        let mut pars = Vec::with_capacity(nb_samples);
        for i in 1..=nbsu {
            if a_flg[i - 1] {
                pars.push(a_pars[i - 1]);
            }
        }
        pars
    }
}

/// The distinct knot values of an expanded knot vector (OCCT
/// Geom_BSplineCurve::Knot(i) counts each distinct knot once).
pub(crate) fn distinct_knots(knots: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if out.is_empty() || (k - out[out.len() - 1]).abs() > 1e-15 {
            out.push(k);
        }
    }
    out
}

/// OCCT Adaptor3d_Curve — the direct adaptor virtuals used where the OCCT
/// code calls the basis-curve adaptor without the tool (InterUtils EstLim*).
pub trait Adaptor3dCurveBasis {
    /// OCCT Adaptor3d_Curve::GetType().
    fn get_type(&self) -> rcad_kernel::base::proj_lib::CurveType;
    /// OCCT Adaptor3d_Curve::Value(U).
    fn value(&self, u: f64) -> rcad_kernel::geom::Point3;
    /// OCCT Adaptor3d_Curve::Line().
    fn line(&self) -> Line3;
    /// OCCT Adaptor3d_Curve::Parabola().
    fn parabola(&self) -> Parabola3;
    /// OCCT Adaptor3d_Curve::Hyperbola().
    fn hyperbola(&self) -> Hyperbola3;
}

/// OCCT Adaptor3d_Surface — the direct adaptor virtuals used on the basis
/// surface of an offset surface (InterUtils EstLimForInfOffs).  The curve
/// adaptor type is the trait parameter BCurve (the OCCT
/// handle<Adaptor3d_Curve>).
pub trait Adaptor3dSurfaceBasis<BCurve: Adaptor3dCurveBasis> {
    /// OCCT Adaptor3d_Surface::GetType().
    fn get_type(&self) -> SurfaceType;
    /// OCCT Adaptor3d_Surface::Plane().
    fn plane(&self) -> Plane;
    /// OCCT Adaptor3d_Surface::Cylinder().
    fn cylinder(&self) -> CylindricalSurface;
    /// OCCT Adaptor3d_Surface::Cone().
    fn cone(&self) -> ConicalSurface;
    /// OCCT Adaptor3d_Surface::Direction().
    fn direction(&self) -> rcad_kernel::geom::Vec3;
    /// OCCT Adaptor3d_Surface::BasisCurve().
    fn basis_curve(&self) -> BCurve;
}

/// OCCT Adaptor3d_HSurfaceTool (Adaptor3d_HSurfaceTool.hxx L40-293 + .cxx
/// L27-115) — the static tool over the surface adaptor `S`, i.e. the OCCT
/// template parameter `TheSurfaceTool`.  UTrim/VTrim (hxx L90-105) return
/// trimmed adaptors, which depend on Adaptor3d_TopolTool (runway 2a-4) and
/// stay documented unimplemented.
pub trait HSurfaceTool {
    type Surface: ?Sized;
    /// The adaptor type returned by BasisCurve (OCCT
    /// handle<Adaptor3d_Curve>).
    type BasisCurve: Adaptor3dCurveBasis;
    /// The adaptor type returned by BasisSurface (OCCT
    /// handle<Adaptor3d_Surface>).
    type BasisSurface: Adaptor3dSurfaceBasis<Self::BasisCurve>;

    /// OCCT FirstUParameter(S) (hxx L45-48).
    fn first_u_parameter(s: &Self::Surface) -> f64;
    /// OCCT FirstVParameter(S) (hxx L50-53).
    fn first_v_parameter(s: &Self::Surface) -> f64;
    /// OCCT LastUParameter(S) (hxx L55-58).
    fn last_u_parameter(s: &Self::Surface) -> f64;
    /// OCCT LastVParameter(S) (hxx L60-63).
    fn last_v_parameter(s: &Self::Surface) -> f64;
    /// OCCT NbUIntervals(S, Sh) (hxx L65-68).
    fn nb_u_intervals(s: &Self::Surface, sh: GeomAbsShape) -> usize;
    /// OCCT NbVIntervals(S, Sh) (hxx L70-73).
    fn nb_v_intervals(s: &Self::Surface, sh: GeomAbsShape) -> usize;
    /// OCCT UIntervals(S, Tab, Sh) (hxx L75-80).
    fn u_intervals(s: &Self::Surface, tab: &mut [f64], sh: GeomAbsShape);
    /// OCCT VIntervals(S, Tab, Sh) (hxx L82-87).
    fn v_intervals(s: &Self::Surface, tab: &mut [f64], sh: GeomAbsShape);
    /// OCCT UTrim(S, First, Last, Tol) (hxx L90-96) — needs
    /// Adaptor3d_TopolTool (runway 2a-4).
    fn u_trim(_s: &Self::Surface, _first: f64, _last: f64, _tol: f64) -> ! {
        unimplemented!("Adaptor3d_HSurfaceTool::UTrim — Adaptor3d_TopolTool dependency (runway 2a-4)");
    }
    /// OCCT VTrim(S, First, Last, Tol) (hxx L99-105) — same dependency.
    fn v_trim(_s: &Self::Surface, _first: f64, _last: f64, _tol: f64) -> ! {
        unimplemented!("Adaptor3d_HSurfaceTool::VTrim — Adaptor3d_TopolTool dependency (runway 2a-4)");
    }
    /// OCCT IsUClosed(S) (hxx L107-110).
    fn is_u_closed(s: &Self::Surface) -> bool;
    /// OCCT IsVClosed(S) (hxx L112-115).
    fn is_v_closed(s: &Self::Surface) -> bool;
    /// OCCT IsUPeriodic(S) (hxx L117-120).
    fn is_u_periodic(s: &Self::Surface) -> bool;
    /// OCCT UPeriod(S) (hxx L122-125).
    fn u_period(s: &Self::Surface) -> f64;
    /// OCCT IsVPeriodic(S) (hxx L127-130).
    fn is_v_periodic(s: &Self::Surface) -> bool;
    /// OCCT VPeriod(S) (hxx L132-135).
    fn v_period(s: &Self::Surface) -> f64;
    /// OCCT Value(S, U, V) (hxx L137-142).
    fn value(s: &Self::Surface, u: f64, v: f64) -> rcad_kernel::geom::Point3;
    /// OCCT D0(S, U, V, P) (hxx L144-150).
    fn d0(s: &Self::Surface, u: f64, v: f64) -> rcad_kernel::geom::Point3;
    /// OCCT D1(S, U, V, P, D1U, D1V) (hxx L152-160) — P is the surface point
    /// itself (the OCCT contract).
    fn d1(s: &Self::Surface, u: f64, v: f64) -> (rcad_kernel::geom::Point3, rcad_kernel::geom::Vec3, rcad_kernel::geom::Vec3);
    /// OCCT D2(S, U, V, P, D1U, D1V, D2U, D2V, D2UV) (hxx L162-173).
    #[allow(clippy::type_complexity)]
    fn d2(
        s: &Self::Surface,
        u: f64,
        v: f64,
    ) -> (
        rcad_kernel::geom::Point3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
    );
    /// OCCT DN(S, U, V, NU, NV) (hxx L203-210).
    fn dn(s: &Self::Surface, u: f64, v: f64, nu: usize, nv: usize) -> rcad_kernel::geom::Vec3;
    /// OCCT UResolution(S, R3d) (hxx L212-215).
    fn u_resolution(s: &Self::Surface, r3d: f64) -> f64;
    /// OCCT VResolution(S, R3d) (hxx L217-220).
    fn v_resolution(s: &Self::Surface, r3d: f64) -> f64;
    /// OCCT GetType(S) (hxx L222-225).
    fn get_type(s: &Self::Surface) -> SurfaceType;
    /// OCCT Plane(S) (hxx L227).
    fn plane(s: &Self::Surface) -> Plane;
    /// OCCT Cylinder(S) (hxx L229-232).
    fn cylinder(s: &Self::Surface) -> CylindricalSurface;
    /// OCCT Cone(S) (hxx L234).
    fn cone(s: &Self::Surface) -> ConicalSurface;
    /// OCCT Torus(S) (hxx L236).
    fn torus(s: &Self::Surface) -> ToroidalSurface;
    /// OCCT Sphere(S) (hxx L238-241).
    fn sphere(s: &Self::Surface) -> SphericalSurface;
    /// OCCT AxeOfRevolution(S) (hxx L253-256) — the gp_Ax1 as (Location,
    /// Direction).
    fn axe_of_revolution(s: &Self::Surface) -> (rcad_kernel::geom::Point3, rcad_kernel::geom::Vec3);
    /// OCCT Direction(S) (hxx L258-261).
    fn direction(s: &Self::Surface) -> rcad_kernel::geom::Vec3;
    /// OCCT BasisCurve(S) (hxx L263-266).
    fn basis_curve(s: &Self::Surface) -> Self::BasisCurve;
    /// OCCT BasisSurface(S) (hxx L268-271).
    fn basis_surface(s: &Self::Surface) -> Self::BasisSurface;
    /// OCCT OffsetValue(S) (hxx L273-276).
    fn offset_value(s: &Self::Surface) -> f64;

    /// OCCT Adaptor3d_Surface::NbUPoles() — used by NbSamplesU(S).
    fn nb_u_poles(s: &Self::Surface) -> usize;
    /// OCCT Adaptor3d_Surface::NbVPoles() — used by NbSamplesV(S).
    fn nb_v_poles(s: &Self::Surface) -> usize;
    /// OCCT Adaptor3d_Surface::NbUKnots() (distinct knots).
    fn nb_u_knots(s: &Self::Surface) -> usize;
    /// OCCT Adaptor3d_Surface::NbVKnots() (distinct knots).
    fn nb_v_knots(s: &Self::Surface) -> usize;
    /// OCCT Adaptor3d_Surface::UDegree().
    fn u_degree(s: &Self::Surface) -> usize;
    /// OCCT Adaptor3d_Surface::VDegree().
    fn v_degree(s: &Self::Surface) -> usize;

    /// OCCT NbSamplesU(S, u1, u2) (Adaptor3d_HSurfaceTool.cxx L72-93).
    fn nb_samples_u(s: &Self::Surface, u1: f64, u2: f64) -> usize {
        let nbs = Self::nb_samples_u_total(s);
        let mut n = nbs;
        if nbs > 10 {
            let uf = Self::first_u_parameter(s);
            let ul = Self::last_u_parameter(s);
            n = (nbs as f64 * ((u2 - u1) / (ul - uf))) as usize;
            if n > nbs || n > 50 {
                n = nbs;
            }
            if n < 5 {
                n = 5;
            }
        }
        n
    }

    /// OCCT NbSamplesV(S, v1, v2) (Adaptor3d_HSurfaceTool.cxx L95-115).
    fn nb_samples_v(s: &Self::Surface, v1: f64, v2: f64) -> usize {
        let nbs = Self::nb_samples_v_total(s);
        let mut n = nbs;
        if nbs > 10 {
            let vf = Self::first_v_parameter(s);
            let vl = Self::last_v_parameter(s);
            n = (nbs as f64 * ((v2 - v1) / (vl - vf))) as usize;
            if n > nbs || n > 50 {
                n = nbs;
            }
            if n < 5 {
                n = 5;
            }
        }
        n
    }

    /// OCCT NbSamplesU(S) (Adaptor3d_HSurfaceTool.cxx L27-45).
    fn nb_samples_u_total(s: &Self::Surface) -> usize {
        match Self::get_type(s) {
            SurfaceType::Plane => 2,
            SurfaceType::BezierSurface => 3 + Self::nb_u_poles(s),
            SurfaceType::BSplineSurface => {
                let nbs = Self::nb_u_knots(s) * Self::u_degree(s);
                if nbs < 2 {
                    2
                } else {
                    nbs
                }
            }
            SurfaceType::Torus => 20,
            _ => 10,
        }
    }

    /// OCCT NbSamplesV(S) (Adaptor3d_HSurfaceTool.cxx L47-70).
    fn nb_samples_v_total(s: &Self::Surface) -> usize {
        match Self::get_type(s) {
            SurfaceType::Plane => 2,
            SurfaceType::BezierSurface => 3 + Self::nb_v_poles(s),
            SurfaceType::BSplineSurface => {
                let nbs = Self::nb_v_knots(s) * Self::v_degree(s);
                if nbs < 2 {
                    2
                } else {
                    nbs
                }
            }
            SurfaceType::Cylinder
            | SurfaceType::Cone
            | SurfaceType::Sphere
            | SurfaceType::Torus
            | SurfaceType::SurfaceOfRevolution
            | SurfaceType::SurfaceOfExtrusion => 15,
            _ => 10,
        }
    }
}

/// The OCCT callback bundle of the HInter assembly: one method per member
/// function IntCurveSurface_HInter.cxx passes into the Inter.pxx templates
/// as a lambda ([this](...) { ... }).  Encoding the lambdas as a trait keeps
/// the engine functions 1:1 with the OCCT template bodies while remaining
/// borrow-safe.  The `done` / `myIsParallel` field references the OCCT
/// lambdas capture become [`HInterHost::done_flag`] /
/// [`HInterHost::is_parallel_flag`].
pub trait HInterHost<
    C: ?Sized,
    CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
    S: ?Sized,
    ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
> {
    /// The `done` field (HInter.cxx passes `done` as `bool& theDone`).
    fn done_flag(&mut self) -> &mut bool;
    /// The `myIsParallel` field (passed as `bool&` into AppendIntAna).
    fn is_parallel_flag(&mut self) -> &mut bool;
    /// `[this]() { this->ResetFields(); }`.
    fn reset_fields(&mut self);
    /// `[this](const IntCurveSurface_IntersectionPoint& pt) { this->Append(pt); }`.
    fn append(&mut self, pt: &IntersectionPoint);
    /// `this->Perform(c, s, u0, v0, u1, v1)` — the bounds-protected Perform
    /// (HInter.cxx L120-150).
    fn perform_bounds(&mut self, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64);
    /// `this->PerformConicSurf(conic, c, s, u1, v1, u2, v2)` — the gp_Lin
    /// overload.
    fn perform_conic_line(&mut self, line: &Line3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64);
    /// The gp_Circ overload.
    fn perform_conic_circle(&mut self, circle: &Circle3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64);
    /// The gp_Elips overload.
    fn perform_conic_ellipse(&mut self, ellipse: &Ellipse3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64);
    /// The gp_Parab overload.
    fn perform_conic_parabola(&mut self, parab: &Parabola3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64);
    /// The gp_Hypr overload.
    fn perform_conic_hyperbola(&mut self, hyper: &Hyperbola3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64);
    /// `this->Perform(c, p, s, ph)` — the 4-arg public Perform (HInter.cxx
    /// L194-219).
    fn perform_polygon_polyhedron(
        &mut self,
        c: &C,
        p: &ThePolygonOfHInter,
        s: &S,
        ph: &ThePolyhedronOfHInter,
    );
    /// `this->InternalPerform(c, p, s, ph, u1, v1, u2, v2)` — the 8-arg
    /// version (HInter.cxx L288-315).  The CSFunction/ExactInter template
    /// arguments are fixed by the concrete host (the OCCT _0.cxx
    /// instantiation).
    fn internal_perform(
        &mut self,
        c: &C,
        p: &ThePolygonOfHInter,
        s: &S,
        ph: &ThePolyhedronOfHInter,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    );
    /// `this->InternalPerform(c, p, s, ph, u1, v1, u2, v2, bsb)` — the 9-arg
    /// version with the bound sort box (HInter.cxx L255-284).
    fn internal_perform_bsb(
        &mut self,
        c: &C,
        p: &ThePolygonOfHInter,
        s: &S,
        ph: &ThePolyhedronOfHInter,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        bsb: &mut BoundSortBox,
    );
    /// `this->InternalPerform(c, p, s, U1, V1, U2, V2)` — the 7-arg version
    /// without polyhedron (HInter.cxx L334-365), instantiating
    /// InternalPerformPolygonBounds with TheQuadCurvExactInter.
    fn internal_perform_bounds(&mut self, c: &C, p: &ThePolygonOfHInter, s: &S, u1: f64, v1: f64, u2: f64, v2: f64);
    /// `this->InternalPerformCurveQuadric(c, s)` (HInter.cxx L319-330).
    fn internal_perform_curve_quadric(&mut self, c: &C, s: &S);
    /// `this->AppendIntAna(c, s, ana)` (HInter.cxx L540-550).
    fn append_int_ana(&mut self, c: &C, s: &S, ana: &IntConicQuad);
}
