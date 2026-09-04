//! IntCurveSurface_TheQuadCurvExactHInter + TheQuadCurvFuncOfTheQuadCurvExactHInter
//! + QuadricCurveExactInterUtils.pxx — the exact intersection of a curve
//! with a quadric surface via math_FunctionAllRoots.
//!
//! 1:1 translations:
//! - [`perform_intersection`] — IntCurveSurface_QuadricCurveExactInterUtils.pxx
//!   (L25-134).
//! - [`TheQuadCurvFuncOfTheQuadCurvExactHInter`] — the hxx/.cxx (L24-58):
//!   the signed distance Q(w) and its derivative.
//! - [`TheQuadCurvExactHInter`] — the hxx/.cxx (L28-79).

use rcad_kernel::math::root::{FunctionAllRoots, FunctionSample, FunctionValue, FunctionWithDerivative};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::int_surf::quadric::Quadric;

use super::inter_utils::QuadCurvExactLike;
use super::{HCurveTool, HSurfaceTool, SurfaceType};

/// OCCT IntCurveSurface_QuadricCurveExactInterUtils EPSX (pxx L29).
const EPSX: f64 = 0.00000000000001;
/// OCCT IntCurveSurface_QuadricCurveExactInterUtils EPSDIST (pxx L30).
const EPSDIST: f64 = 0.00000001;
/// OCCT IntCurveSurface_QuadricCurveExactInterUtils EPSNUL (pxx L31).
const EPSNUL: f64 = 0.00000001;

/// OCCT IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter
/// (TheQuadCurvFuncOfTheQuadCurvExactHInter.hxx L27-55 + .cxx L24-58) —
/// the signed distance between the implicit surface and the point at
/// parameter Param on the curve, with its first derivative.
pub struct TheQuadCurvFuncOfTheQuadCurvExactHInter<'a, C: ?Sized, CT: HCurveTool<Curve = C>> {
    my_quadric: Quadric,
    my_curve: &'a C,
    _tool: std::marker::PhantomData<fn(&CT)>,
}

impl<'a, C: ?Sized, CT: HCurveTool<Curve = C>> TheQuadCurvFuncOfTheQuadCurvExactHInter<'a, C, CT> {
    /// OCCT TheQuadCurvFuncOfTheQuadCurvExactHInter(Q, C) (cxx L24-30).
    pub fn new(q: Quadric, c: &'a C) -> Self {
        TheQuadCurvFuncOfTheQuadCurvExactHInter {
            my_quadric: q,
            my_curve: c,
            _tool: std::marker::PhantomData,
        }
    }

    /// OCCT Value(Param, F) (cxx L32-36) — F = the quadric distance at
    /// C(Param); always returns true in OCCT.
    pub fn value(&mut self, param: f64) -> f64 {
        self.my_quadric.distance(<CT as HCurveTool>::value(self.my_curve, param))
    }

    /// OCCT Derivative(Param, D) (cxx L38-46) — D = C'(Param)·grad
    /// Q(C(Param)); always returns true in OCCT.
    pub fn derivative(&mut self, param: f64) -> f64 {
        let (p, t) = <CT as HCurveTool>::d1(self.my_curve, param);
        t.dot(self.my_quadric.gradient(p))
    }

    /// OCCT Values(Param, F, D) (cxx L48-58) — value and derivative from one
    /// quadric evaluation (ValAndGrad).
    pub fn values(&mut self, param: f64) -> (f64, f64) {
        let (p, t) = <CT as HCurveTool>::d1(self.my_curve, param);
        let (f, grad) = self.my_quadric.val_and_grad(p);
        let d = t.dot(grad);
        (f, d)
    }
}

impl<'a, C: ?Sized, CT: HCurveTool<Curve = C>> FunctionValue
    for TheQuadCurvFuncOfTheQuadCurvExactHInter<'a, C, CT>
{
    fn value(&mut self, x: f64) -> Option<f64> {
        Some(TheQuadCurvFuncOfTheQuadCurvExactHInter::value(self, x))
    }
}

impl<'a, C: ?Sized, CT: HCurveTool<Curve = C>> FunctionWithDerivative
    for TheQuadCurvFuncOfTheQuadCurvExactHInter<'a, C, CT>
{
    fn derivative(&mut self, x: f64) -> Option<f64> {
        Some(TheQuadCurvFuncOfTheQuadCurvExactHInter::derivative(self, x))
    }

    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        Some(TheQuadCurvFuncOfTheQuadCurvExactHInter::values(self, x))
    }
}

/// The construction contract of the QuadCurvFuncType template parameter
/// (OCCT: `QuadCurvFuncType aFunction(aQuadric, theCurve)`).
pub trait TheQuadCurvFuncNew<'a, C: ?Sized, CT: HCurveTool<Curve = C>> {
    fn new(q: Quadric, c: &'a C) -> Self;
}

impl<'a, C: ?Sized, CT: HCurveTool<Curve = C>> TheQuadCurvFuncNew<'a, C, CT>
    for TheQuadCurvFuncOfTheQuadCurvExactHInter<'a, C, CT>
{
    fn new(q: Quadric, c: &'a C) -> Self {
        TheQuadCurvFuncOfTheQuadCurvExactHInter::new(q, c)
    }
}

/// OCCT IntCurveSurface_TheQuadCurvExactHInter (TheQuadCurvExactHInter.hxx
/// L28-55 + .cxx L29-78) — the roots of the signed distance function.  The
/// struct stores no adaptor references (like the OCCT value class); the
/// `'a`/tool parameters only carry the instantiation identity.
pub struct TheQuadCurvExactHInter<'a, S: ?Sized, ST: HSurfaceTool<Surface = S>, C: ?Sized, CT: HCurveTool<Curve = C>> {
    nbpnts: i32,
    pnts: Vec<f64>,
    nbintv: i32,
    intv: Vec<f64>,
    _types: std::marker::PhantomData<fn(&'a (), &ST, &CT)>,
}

impl<'a, S: ?Sized, ST: HSurfaceTool<Surface = S>, C: ?Sized, CT: HCurveTool<Curve = C>>
    TheQuadCurvExactHInter<'a, S, ST, C, CT>
{
    /// OCCT TheQuadCurvExactHInter(S, C) (cxx L29-41) — runs
    /// QuadricCurveExactInterUtils::PerformIntersection.
    pub fn new(s: &S, c: &'a C) -> Self {
        let mut nbpnts = -1i32;
        let mut pnts: Vec<f64> = Vec::new();
        let mut nbintv = -1i32;
        let mut intv: Vec<f64> = Vec::new();
        perform_intersection::<S, ST, C, CT, TheQuadCurvFuncOfTheQuadCurvExactHInter<C, CT>>(
            s, c, &mut pnts, &mut intv, &mut nbpnts, &mut nbintv,
        );
        TheQuadCurvExactHInter {
            nbpnts,
            pnts,
            nbintv,
            intv,
            _types: std::marker::PhantomData,
        }
    }

    /// OCCT IsDone() (cxx L45-48).
    pub fn is_done(&self) -> bool {
        self.nbpnts != -1
    }

    /// OCCT NbRoots() (cxx L52-55).
    pub fn nb_roots(&self) -> usize {
        self.nbpnts as usize
    }

    /// OCCT Root(Index) (cxx L66-69) — 1-based.
    pub fn root(&self, index: usize) -> f64 {
        self.pnts[index - 1]
    }

    /// OCCT NbIntervals() (cxx L59-62).
    pub fn nb_intervals(&self) -> i32 {
        self.nbintv
    }

    /// OCCT Intervals(Index, U1, U2) (cxx L73-78) — 1-based; returns (U1,
    /// U2) of the segment on the curve.
    pub fn intervals(&self, index: usize) -> (f64, f64) {
        let index2 = index + index - 1;
        (self.intv[index2 - 1], self.intv[index2])
    }
}

impl<'a, S: ?Sized, ST: HSurfaceTool<Surface = S>, C: ?Sized, CT: HCurveTool<Curve = C>>
    QuadCurvExactLike<C, CT, S, ST> for TheQuadCurvExactHInter<'a, S, ST, C, CT>
{
    fn construct(s: &S, c: &C) -> Self {
        TheQuadCurvExactHInter::new(s, c)
    }
    fn is_done(&self) -> bool {
        TheQuadCurvExactHInter::is_done(self)
    }
    fn nb_roots(&self) -> usize {
        TheQuadCurvExactHInter::nb_roots(self)
    }
    fn root(&self, index: usize) -> f64 {
        TheQuadCurvExactHInter::root(self, index)
    }
}

/// OCCT IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection
/// (QuadricCurveExactInterUtils.pxx L45-132).
pub fn perform_intersection<'a, S: ?Sized, ST, C: ?Sized, CT, QuadCurvFuncType>(
    the_surface: &S,
    the_curve: &'a C,
    the_pnts: &mut Vec<f64>,
    the_intv: &mut Vec<f64>,
    the_nb_pnts: &mut i32,
    the_nb_intv: &mut i32,
) where
    ST: HSurfaceTool<Surface = S>,
    CT: HCurveTool<Curve = C>,
    QuadCurvFuncType: TheQuadCurvFuncNew<'a, C, CT> + FunctionValue + FunctionWithDerivative,
{
    *the_nb_pnts = -1;
    *the_nb_intv = -1;
    the_pnts.clear();
    the_intv.clear();

    let a_quadric_type = <ST as HSurfaceTool>::get_type(the_surface);
    let a_quadric: Quadric = match a_quadric_type {
        SurfaceType::Plane => Quadric::from_plane(&<ST as HSurfaceTool>::plane(the_surface)),
        SurfaceType::Cylinder => Quadric::from_cylinder(&<ST as HSurfaceTool>::cylinder(the_surface)),
        SurfaceType::Cone => Quadric::from_cone(&<ST as HSurfaceTool>::cone(the_surface)),
        SurfaceType::Sphere => Quadric::from_sphere(&<ST as HSurfaceTool>::sphere(the_surface)),
        _ => return,
    };

    let a_nb_intervals = <CT as HCurveTool>::nb_intervals(the_curve, GeomAbsShape::C1);
    let mut an_intervals = vec![0.0f64; a_nb_intervals + 1];

    <CT as HCurveTool>::intervals(the_curve, &mut an_intervals, GeomAbsShape::C1);

    let mut ii = 1usize;
    while ii <= a_nb_intervals {
        let u1 = an_intervals[ii - 1];
        let u2 = an_intervals[ii];

        let a_sample = FunctionSample::new(u1, u2, <CT as HCurveTool>::nb_samples(the_curve, u1, u2) as i32);
        // OCCT passes the quadric by const ref; the function keeps its own
        // copy.
        let mut a_function = QuadCurvFuncType::new(a_quadric.clone(), the_curve);
        let a_roots = FunctionAllRoots::new(&mut a_function, &a_sample, EPSX, EPSDIST, EPSNUL);

        if a_roots.is_done() {
            let a_nb_points = a_roots.nb_points();
            let a_nb_root_intv = a_roots.nb_intervals();

            for i in 1..=a_nb_points {
                the_pnts.push(a_roots.get_point(i));
            }

            for i in 1..=a_nb_root_intv {
                let (a, b) = a_roots.get_interval(i);
                the_intv.push(a);
                the_intv.push(b);
            }
        } else {
            break;
        }
        ii += 1;
    }

    if ii > a_nb_intervals {
        *the_nb_pnts = the_pnts.len() as i32;
        *the_nb_intv = (the_intv.len() / 2) as i32;
    }
}
