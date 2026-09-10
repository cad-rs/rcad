//! OCCT IntImp package (TKGeomAlgo) — zero-of a curve-vs-surface function
//! set, solved by math_FunctionSetRoot.
//!
//! 1:1 translations:
//! - [`zer_cs_par_func::ZerCSParFunc`] — `IntImp_ZerCSParFunc.gxx`
//!   (L23-116): F(u,v,w) = S(u,v) - C(w), the function set consumed by
//!   [`int_cs::IntCS`].  The HLRBRep instantiation is
//!   HLRBRep_TheCSFunctionOfInterCSurf (ThePSurface = HLRBRep_Surface,
//!   TheCurve = HLRBRep_Curve).
//! - [`int_cs::IntCS`] — `IntImp_IntCS.gxx` (L26-198): the root solver
//!   with the 3-attempt w-restart loop.  The HLRBRep instantiation is
//!   HLRBRep_TheExactInterCSurf.

pub mod int_cs;
pub mod zer_cs_par_func;

pub use int_cs::IntCS;
pub use zer_cs_par_func::ZerCSParFunc;

use glam::DVec3;

/// OCCT `ThePSurfaceTool` template parameter of the IntImp generics —
/// static accessors over the parametrised surface.
pub trait PSurfaceTool {
    type Surface;
    /// OCCT ThePSurfaceTool::Value(S, U, V).
    fn value(s: &Self::Surface, u: f64, v: f64) -> DVec3;
    /// OCCT ThePSurfaceTool::D1(S, U, V, P, D1u, D1v).
    fn d1(s: &Self::Surface, u: f64, v: f64) -> (DVec3, DVec3, DVec3);
    /// OCCT ThePSurfaceTool::FirstUParameter(S).
    fn first_u_parameter(s: &Self::Surface) -> f64;
    /// OCCT ThePSurfaceTool::LastUParameter(S).
    fn last_u_parameter(s: &Self::Surface) -> f64;
    /// OCCT ThePSurfaceTool::FirstVParameter(S).
    fn first_v_parameter(s: &Self::Surface) -> f64;
    /// OCCT ThePSurfaceTool::LastVParameter(S).
    fn last_v_parameter(s: &Self::Surface) -> f64;
    /// OCCT ThePSurfaceTool::UResolution(S, R3d).
    fn u_resolution(s: &Self::Surface, r3d: f64) -> f64;
    /// OCCT ThePSurfaceTool::VResolution(S, R3d).
    fn v_resolution(s: &Self::Surface, r3d: f64) -> f64;
}

/// OCCT `TheCurveTool` template parameter of the IntImp generics (3D).
pub trait CurveTool3d {
    type Curve;
    /// OCCT TheCurveTool::Value(C, U).
    fn value(c: &Self::Curve, u: f64) -> DVec3;
    /// OCCT TheCurveTool::D1(C, U, P, T).
    fn d1(c: &Self::Curve, u: f64) -> (DVec3, DVec3);
    /// OCCT TheCurveTool::FirstParameter(C).
    fn first_parameter(c: &Self::Curve) -> f64;
    /// OCCT TheCurveTool::LastParameter(C).
    fn last_parameter(c: &Self::Curve) -> f64;
    /// OCCT TheCurveTool::Resolution(C, R3d).
    fn resolution(c: &Self::Curve, r3d: f64) -> f64;
}

/// OCCT Precision::Confusion() (Precision.hxx L165).
pub const PRECISION_CONFUSION: f64 = 1.0e-7;
/// OCCT Precision::SquareConfusion() (Precision.hxx L169).
pub const PRECISION_SQUARE_CONFUSION: f64 = 1.0e-14;
/// OCCT Precision::Infinite() (Precision.hxx L371).
const PRECISION_INFINITE: f64 = 2.0e100;

/// OCCT Precision::IsInfinite(R) (Precision.hxx L350-353).
fn is_infinite(r: f64) -> bool {
    r.abs() >= (0.5 * PRECISION_INFINITE)
}
