//! OCCT Extrema_GlobOptFuncCC (TKGeomBase/Extrema/Extrema_GlobOptFuncCC.hxx
//! L16-101 and Extrema_GlobOptFuncCC.cxx L24-388) — the objective functions of
//! the Evtushenko curve-curve search:
//! - `Extrema_GlobOptFuncCCC0` (C0 objective, value only),
//! - `Extrema_GlobOptFuncCCC1` (C1 objective, value + gradient),
//! - `Extrema_GlobOptFuncCCC2` (C2 objective, value + gradient + Hessian).
//!
//! Also carries the file-local statics `_NbVariables` (cxx L24-27), the two
//! `_Value` overloads (cxx L30-65), the two `_Gradient` overloads
//! (cxx L70-119) and the two `_Hessian` overloads (cxx L122-185).
//!
//! rcad encoding: the OCCT `myType` discriminator (1 = 3D adaptor, 2 = 2D
//! adaptor) selects which pair of the file-local statics runs; the rcad enum
//! [`GlobOptCurves`] carries the same two cases.

use glam::DVec3;

use crate::base::extrema_math_opt::MultipleVarFunctionWithHessian;
use crate::base::extrema_curve_tool::ExtremaCurveTool;
use crate::base::proj_lib::adaptor::Adaptor2dCurve2d;
use crate::math::math_bfgs::{MultipleVarFunction, MultipleVarFunctionWithGradient};
use crate::math::math_matrix::{Matrix, Vector};

/// OCCT static `_NbVariables()` (cxx L24-27).
fn nb_variables_impl() -> i32 {
    2
}

/// The OCCT `myType` discriminator (cxx L198/L225) — 1 selects the 3D pair of
/// statics, 2 the 2D pair.
pub enum GlobOptCurves<'a> {
    /// OCCT myType == 1 — `const Adaptor3d_Curve* myC1_3d, *myC2_3d`.
    Curves3d(&'a dyn ExtremaCurveTool, &'a dyn ExtremaCurveTool),
    /// OCCT myType == 2 — `const Adaptor2d_Curve2d* myC1_2d, *myC2_2d`.
    Curves2d(&'a dyn Adaptor2dCurve2d, &'a dyn Adaptor2dCurve2d),
}

/// OCCT static `_Value(C1, C2, X, F)` (cxx L30-46) for the 3D pair.
fn value_3d(c1: &dyn ExtremaCurveTool, c2: &dyn ExtremaCurveTool, x: &Vector, f: &mut f64) -> bool {
    let u = x.get(1);
    let v = x.get(2);
    if u < c1.first_parameter() || u > c1.last_parameter() || v < c2.first_parameter()
        || v > c2.last_parameter()
    {
        return false;
    }
    // OCCT: F = C2.Value(v).Distance(C1.Value(u)).
    *f = (c2.value(v) - c1.value(u)).length();
    true
}

/// OCCT static `_Value(C1, C2, X, F)` (cxx L49-65) for the 2D pair.
fn value_2d(c1: &dyn Adaptor2dCurve2d, c2: &dyn Adaptor2dCurve2d, x: &Vector, f: &mut f64) -> bool {
    let u = x.get(1);
    let v = x.get(2);
    if u < c1.first_parameter() || u > c1.last_parameter() || v < c2.first_parameter()
        || v > c2.last_parameter()
    {
        return false;
    }
    *f = (c2.value(v) - c1.value(u)).length();
    true
}

/// OCCT static `_Gradient(C1, C2, X, G)` (cxx L70-93) for the 3D pair.
fn gradient_3d(c1: &dyn ExtremaCurveTool, c2: &dyn ExtremaCurveTool, x: &Vector, g: &mut Vector) -> bool {
    if x.get(1) < c1.first_parameter()
        || x.get(1) > c1.last_parameter()
        || x.get(2) < c2.first_parameter()
        || x.get(2) > c2.last_parameter()
    {
        return false;
    }

    let (c1d0, c1d1) = c1.d1(x.get(1));
    let (c2d0, c2d1) = c2.d1(x.get(2));

    let mut g1 = -(c2d0.x - c1d0.x) * c1d1.x
        - (c2d0.y - c1d0.y) * c1d1.y
        - (c2d0.z - c1d0.z) * c1d1.z;
    let mut g2 = (c2d0.x - c1d0.x) * c2d1.x
        + (c2d0.y - c1d0.y) * c2d1.y
        + (c2d0.z - c1d0.z) * c2d1.z;
    g1 *= 2.0;
    g2 *= 2.0;
    g.set(1, g1);
    g.set(2, g2);
    true
}

/// OCCT static `_Gradient(C1, C2, X, G)` (cxx L96-119) for the 2D pair.
fn gradient_2d(c1: &dyn Adaptor2dCurve2d, c2: &dyn Adaptor2dCurve2d, x: &Vector, g: &mut Vector) -> bool {
    if x.get(1) < c1.first_parameter()
        || x.get(1) > c1.last_parameter()
        || x.get(2) < c2.first_parameter()
        || x.get(2) > c2.last_parameter()
    {
        return false;
    }

    let (c1d0, c1d1) = c1.d1(x.get(1));
    let (c2d0, c2d1) = c2.d1(x.get(2));

    let mut g1 =
        -(c2d0.x - c1d0.x) * c1d1.x - (c2d0.y - c1d0.y) * c1d1.y;
    let mut g2 = (c2d0.x - c1d0.x) * c2d1.x + (c2d0.y - c1d0.y) * c2d1.y;
    g1 *= 2.0;
    g2 *= 2.0;
    g.set(1, g1);
    g.set(2, g2);
    true
}

/// OCCT static `_Hessian(C1, C2, X, H)` (cxx L122-153) for the 3D pair.
fn hessian_3d(c1: &dyn ExtremaCurveTool, c2: &dyn ExtremaCurveTool, x: &Vector, h: &mut Matrix) -> bool {
    if x.get(1) < c1.first_parameter()
        || x.get(1) > c1.last_parameter()
        || x.get(2) < c2.first_parameter()
        || x.get(2) > c2.last_parameter()
    {
        return false;
    }

    let (c1d0, c1d1, c1d2) = c1.d2(x.get(1));
    let (c2d0, c2d1, c2d2) = c2.d2(x.get(2));

    let h11 = c1d1.x * c1d1.x
        + c1d1.y * c1d1.y
        + c1d1.z * c1d1.z
        - (c2d0.x - c1d0.x) * c1d2.x
        - (c2d0.y - c1d0.y) * c1d2.y
        - (c2d0.z - c1d0.z) * c1d2.z;
    let h12 = -c2d1.x * c1d1.x - c2d1.y * c1d1.y - c2d1.z * c1d1.z;
    let h22 = c2d1.x * c2d1.x
        + c2d1.y * c2d1.y
        + c2d1.z * c2d1.z
        + (c2d0.x - c1d0.x) * c2d2.x
        + (c2d0.y - c1d0.y) * c2d2.y
        + (c2d0.z - c1d0.z) * c2d2.z;
    h.set(1, 1, h11 * 2.0);
    h.set(1, 2, h12 * 2.0);
    h.set(2, 1, h12 * 2.0);
    h.set(2, 2, h22 * 2.0);
    true
}

/// OCCT static `_Hessian(C1, C2, X, H)` (cxx L156-185) for the 2D pair.
fn hessian_2d(c1: &dyn Adaptor2dCurve2d, c2: &dyn Adaptor2dCurve2d, x: &Vector, h: &mut Matrix) -> bool {
    if x.get(1) < c1.first_parameter()
        || x.get(1) > c1.last_parameter()
        || x.get(2) < c2.first_parameter()
        || x.get(2) > c2.last_parameter()
    {
        return false;
    }

    let (c1d0, c1d1, c1d2) = c1.d2(x.get(1));
    let (c2d0, c2d1, c2d2) = c2.d2(x.get(2));

    let h11 = c1d1.x * c1d1.x + c1d1.y * c1d1.y
        - (c2d0.x - c1d0.x) * c1d2.x
        - (c2d0.y - c1d0.y) * c1d2.y;
    let h12 = -c2d1.x * c1d1.x - c2d1.y * c1d1.y;
    let h22 = c2d1.x * c2d1.x
        + c2d1.y * c2d1.y
        + (c2d0.x - c1d0.x) * c2d2.x
        + (c2d0.y - c1d0.y) * c2d2.y;
    h.set(1, 1, h11 * 2.0);
    h.set(1, 2, h12 * 2.0);
    h.set(2, 1, h12 * 2.0);
    h.set(2, 2, h22 * 2.0);
    true
}

/// OCCT `_Value` dispatch on myType (cxx L222-232, L269-279, L337-347).
fn value_dispatch(curves: &GlobOptCurves<'_>, x: &Vector, f: &mut f64) -> bool {
    match curves {
        GlobOptCurves::Curves3d(c1, c2) => value_3d(*c1, *c2, x, f),
        GlobOptCurves::Curves2d(c1, c2) => value_2d(*c1, *c2, x, f),
    }
}

/// OCCT `_Gradient` dispatch on myType (cxx L283-293, L351-361).
fn gradient_dispatch(curves: &GlobOptCurves<'_>, x: &Vector, g: &mut Vector) -> bool {
    match curves {
        GlobOptCurves::Curves3d(c1, c2) => gradient_3d(*c1, *c2, x, g),
        GlobOptCurves::Curves2d(c1, c2) => gradient_2d(*c1, *c2, x, g),
    }
}

/// OCCT `_Hessian` dispatch on myType (cxx L377-385).
fn hessian_dispatch(curves: &GlobOptCurves<'_>, x: &Vector, h: &mut Matrix) -> bool {
    match curves {
        GlobOptCurves::Curves3d(c1, c2) => hessian_3d(*c1, *c2, x, h),
        GlobOptCurves::Curves2d(c1, c2) => hessian_2d(*c1, *c2, x, h),
    }
}

// =============================================================================
// Extrema_GlobOptFuncCCC0 (hxx L27-44, cxx L191-232)
// =============================================================================

/// OCCT Extrema_GlobOptFuncCCC0 (hxx L27-44).
pub struct GlobOptFuncCCC0<'a> {
    /// hxx L41-42: const Adaptor3d_Curve *myC1_3d, *myC2_3d; and the 2d pair.
    /// hxx L43: int myType.
    curves: GlobOptCurves<'a>,
}

impl<'a> GlobOptFuncCCC0<'a> {
    /// OCCT Extrema_GlobOptFuncCCC0(C1, C2) — the 3D form (cxx L191-199).
    pub fn new_3d(c1: &'a dyn ExtremaCurveTool, c2: &'a dyn ExtremaCurveTool) -> Self {
        GlobOptFuncCCC0 {
            curves: GlobOptCurves::Curves3d(c1, c2),
        }
    }

    /// OCCT Extrema_GlobOptFuncCCC0(C1, C2) — the 2D form (cxx L203-211).
    pub fn new_2d(c1: &'a dyn Adaptor2dCurve2d, c2: &'a dyn Adaptor2dCurve2d) -> Self {
        GlobOptFuncCCC0 {
            curves: GlobOptCurves::Curves2d(c1, c2),
        }
    }
}

impl MultipleVarFunction for GlobOptFuncCCC0<'_> {
    /// OCCT NbVariables() (cxx L215-218).
    fn nb_variables(&self) -> i32 {
        nb_variables_impl()
    }

    /// OCCT Value(X, F) (cxx L222-232).
    fn value(&mut self, x: &Vector, f: &mut f64) -> bool {
        value_dispatch(&self.curves, x, f)
    }
}

// =============================================================================
// Extrema_GlobOptFuncCCC1 (hxx L48-69, cxx L238-300)
// =============================================================================

/// OCCT Extrema_GlobOptFuncCCC1 (hxx L48-69).
pub struct GlobOptFuncCCC1<'a> {
    /// hxx L66-68: the curve pointers and myType.
    curves: GlobOptCurves<'a>,
}

impl<'a> GlobOptFuncCCC1<'a> {
    /// OCCT Extrema_GlobOptFuncCCC1(C1, C2) — the 3D form (cxx L238-246).
    pub fn new_3d(c1: &'a dyn ExtremaCurveTool, c2: &'a dyn ExtremaCurveTool) -> Self {
        GlobOptFuncCCC1 {
            curves: GlobOptCurves::Curves3d(c1, c2),
        }
    }

    /// OCCT Extrema_GlobOptFuncCCC1(C1, C2) — the 2D form (cxx L250-258).
    pub fn new_2d(c1: &'a dyn Adaptor2dCurve2d, c2: &'a dyn Adaptor2dCurve2d) -> Self {
        GlobOptFuncCCC1 {
            curves: GlobOptCurves::Curves2d(c1, c2),
        }
    }
}

impl MultipleVarFunction for GlobOptFuncCCC1<'_> {
    /// OCCT NbVariables() (cxx L262-265).
    fn nb_variables(&self) -> i32 {
        nb_variables_impl()
    }

    /// OCCT Value(X, F) (cxx L269-279).
    fn value(&mut self, x: &Vector, f: &mut f64) -> bool {
        value_dispatch(&self.curves, x, f)
    }
}

impl MultipleVarFunctionWithGradient for GlobOptFuncCCC1<'_> {
    /// OCCT Gradient(X, G) (cxx L283-293).
    fn gradient(&mut self, x: &Vector, g: &mut Vector) -> bool {
        gradient_dispatch(&self.curves, x, g)
    }

    /// OCCT Values(X, F, G) (cxx L297-300).
    fn values(&mut self, x: &Vector, f: &mut f64, g: &mut Vector) -> bool {
        self.value(x, f) && self.gradient(x, g)
    }
}

// =============================================================================
// Extrema_GlobOptFuncCCC2 (hxx L73-99, cxx L306-388)
// =============================================================================

/// OCCT Extrema_GlobOptFuncCCC2 (hxx L73-99).
pub struct GlobOptFuncCCC2<'a> {
    /// hxx L96-98: the curve pointers and myType.
    curves: GlobOptCurves<'a>,
}

impl<'a> GlobOptFuncCCC2<'a> {
    /// OCCT Extrema_GlobOptFuncCCC2(C1, C2) — the 3D form (cxx L306-314).
    pub fn new_3d(c1: &'a dyn ExtremaCurveTool, c2: &'a dyn ExtremaCurveTool) -> Self {
        GlobOptFuncCCC2 {
            curves: GlobOptCurves::Curves3d(c1, c2),
        }
    }

    /// OCCT Extrema_GlobOptFuncCCC2(C1, C2) — the 2D form (cxx L318-326).
    pub fn new_2d(c1: &'a dyn Adaptor2dCurve2d, c2: &'a dyn Adaptor2dCurve2d) -> Self {
        GlobOptFuncCCC2 {
            curves: GlobOptCurves::Curves2d(c1, c2),
        }
    }
}

impl MultipleVarFunction for GlobOptFuncCCC2<'_> {
    /// OCCT NbVariables() (cxx L330-333).
    fn nb_variables(&self) -> i32 {
        nb_variables_impl()
    }

    /// OCCT Value(X, F) (cxx L337-347).
    fn value(&mut self, x: &Vector, f: &mut f64) -> bool {
        value_dispatch(&self.curves, x, f)
    }
}

impl MultipleVarFunctionWithGradient for GlobOptFuncCCC2<'_> {
    /// OCCT Gradient(X, G) (cxx L351-361).
    fn gradient(&mut self, x: &Vector, g: &mut Vector) -> bool {
        gradient_dispatch(&self.curves, x, g)
    }

    /// OCCT Values(X, F, G) (cxx L365-368).
    fn values(&mut self, x: &Vector, f: &mut f64, g: &mut Vector) -> bool {
        self.value(x, f) && self.gradient(x, g)
    }
}

impl MultipleVarFunctionWithHessian for GlobOptFuncCCC2<'_> {
    /// OCCT Values(X, F, G, H) (cxx L372-388).
    fn values_hessian(&mut self, x: &Vector, f: &mut f64, g: &mut Vector, h: &mut Matrix) -> bool {
        let is_hessian_computed = hessian_dispatch(&self.curves, x, h);
        self.value(x, f) && self.gradient(x, g) && is_hessian_computed
    }
}

/// OCCT gp_Pnt::Distance — the 3D distance helper used by the objective.
#[allow(dead_code)]
fn pnt_distance(a: DVec3, b: DVec3) -> f64 {
    (a - b).length()
}
