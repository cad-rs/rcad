//! OCCT IntImp_ZerCSParFunc (TKGeomAlgo IntImp package).
//!
//! 1:1 translation of `IntImp_ZerCSParFunc.gxx` (L23-116) as an
//! implementation of `FunctionSetWithDerivatives` over the
//! [`PSurfaceTool`]/[`CurveTool3d`] template parameters.

use glam::DVec3;
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::int_cs::ZerCSAccessors;
use super::{CurveTool3d, PSurfaceTool};

/// OCCT IntImp_ZerCSParFunc — the function set F(u, v, w) = S(u,v) - C(w)
/// whose zero is a curve/surface intersection point.  It also caches the
/// squared distance `f` and the mid-point `p` after each evaluation, as
/// OCCT does.
pub struct ZerCSParFunc<'a, S, C, PS: PSurfaceTool<Surface = S>, CT: CurveTool3d<Curve = C>> {
    surface: &'a S,
    curve: &'a C,
    p: DVec3,
    f: f64,
    _tools: std::marker::PhantomData<fn(&S, &C) -> (PS, CT)>,
}

impl<'a, S, C, PS: PSurfaceTool<Surface = S>, CT: CurveTool3d<Curve = C>>
    ZerCSParFunc<'a, S, C, PS, CT>
{
    /// OCCT IntImp_ZerCSParFunc(S, C) — gxx L23-29.
    pub fn new(surface: &'a S, curve: &'a C) -> Self {
        ZerCSParFunc {
            surface,
            curve,
            p: DVec3::ZERO,
            f: 0.0,
            _tools: std::marker::PhantomData,
        }
    }

    /// OCCT Point() — gxx L98-101.
    pub fn point(&self) -> DVec3 {
        self.p
    }

    /// OCCT Root() — gxx L103-106.
    pub fn root(&self) -> f64 {
        self.f
    }

    /// OCCT AuxillarSurface() — gxx L108-111.
    pub fn auxillar_surface(&self) -> &'a S {
        self.surface
    }

    /// OCCT AuxillarCurve() — gxx L113-116.
    pub fn auxillar_curve(&self) -> &'a C {
        self.curve
    }
}

impl<'a, S, C, PS: PSurfaceTool<Surface = S>, CT: CurveTool3d<Curve = C>> ZerCSAccessors<S, C>
    for ZerCSParFunc<'a, S, C, PS, CT>
{
    fn root(&self) -> f64 {
        ZerCSParFunc::root(self)
    }
    fn point(&self) -> DVec3 {
        ZerCSParFunc::point(self)
    }
    fn auxillar_surface(&self) -> &S {
        ZerCSParFunc::auxillar_surface(self)
    }
    fn auxillar_curve(&self) -> &C {
        ZerCSParFunc::auxillar_curve(self)
    }
}

impl<'a, S, C, PS: PSurfaceTool<Surface = S>, CT: CurveTool3d<Curve = C>>
    FunctionSetWithDerivatives for ZerCSParFunc<'a, S, C, PS, CT>
{
    /// OCCT NbVariables() — gxx L31-34.
    fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations() — gxx L36-39.
    fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(X, F) — gxx L41-53.
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let psurf = PS::value(self.surface, x[0], x[1]);
        let pcurv = CT::value(self.curve, x[2]);
        f[0] = psurf.x - pcurv.x;
        f[1] = psurf.y - pcurv.y;
        let f3 = psurf.z - pcurv.z;
        f[2] = f3;
        self.f = f[0] * f[0] + f[1] * f[1] + f3 * f3;
        self.p = (psurf + pcurv) * 0.5;
        true
    }

    /// OCCT Derivatives(X, D) — gxx L55-71.
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        let (_psurf, d1u, d1v) = PS::d1(self.surface, x[0], x[1]);
        let (_pcurv, d1w) = CT::d1(self.curve, x[2]);
        // D(1,1..3) = D1u.X, D1v.X, -D1w.X, etc. — rows are the equations,
        // columns the variables (math_Matrix 1-based layout).
        df[0][0] = d1u.x;
        df[0][1] = d1v.x;
        df[0][2] = -d1w.x;
        df[1][0] = d1u.y;
        df[1][1] = d1v.y;
        df[1][2] = -d1w.y;
        df[2][0] = d1u.z;
        df[2][1] = d1v.z;
        df[2][2] = -d1w.z;
        true
    }

    /// OCCT Values(X, F, D) — gxx L73-96.
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        let (psurf, d1u, d1v) = PS::d1(self.surface, x[0], x[1]);
        let (pcurv, d1w) = CT::d1(self.curve, x[2]);
        df[0][0] = d1u.x;
        df[0][1] = d1v.x;
        df[0][2] = -d1w.x;
        df[1][0] = d1u.y;
        df[1][1] = d1v.y;
        df[1][2] = -d1w.y;
        df[2][0] = d1u.z;
        df[2][1] = d1v.z;
        df[2][2] = -d1w.z;

        f[0] = psurf.x - pcurv.x;
        f[1] = psurf.y - pcurv.y;
        let f3 = psurf.z - pcurv.z;
        f[2] = f3;
        self.f = f[0] * f[0] + f[1] * f[1] + f3 * f3;
        self.p = (psurf + pcurv) * 0.5;
        true
    }
}
