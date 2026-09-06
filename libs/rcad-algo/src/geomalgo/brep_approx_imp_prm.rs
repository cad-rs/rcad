//! OCCT BRepApprox ImpPrm chain (the implicit x parametric surfaces functor
//! of BRepApprox_Approx): 1:1 translations of
//! - [`ZerImpFunc`] — `BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox`,
//!   the `IntImp_ZerImpFunc.gxx` (L27-156) / `IntImp_ZerImpFunc.lxx`
//!   (L17-69) instantiation over the alias table below;
//! - [`ImpPrmSvSurfaces`] — `BRepApprox_TheImpPrmSvSurfacesOfApprox`, the
//!   `ApproxInt_ImpPrmSvSurfaces.gxx` instantiation (L326-1020), plus the
//!   file statics [`is_singular`] (gxx L32-57), [`singular_processing`]
//!   (gxx L80-268) and [`nonsingular_processing`] (gxx L287-322).
//!
//! Alias table of `BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox_0.cxx`
//! (L29-39) and `BRepApprox_TheImpPrmSvSurfacesOfApprox_0.cxx` (L29-44):
//! ```text
//!   ThePSurface      = BRepAdaptor_Surface    -> topalgo::brep_adaptor::surface::BRepAdaptorSurface
//!   ThePSurfaceTool  = BRepApprox_SurfaceTool -> brep_approx::surface_tool
//!   TheISurface      = IntSurf_Quadric        -> int_surf::Quadric
//!   TheISurfaceTool  = IntSurf_QuadricTool    -> int_surf::quadric::Quadric tool_* methods
//!   TheLine          = BRepApprox_ApproxLine  -> brep_approx::ApproxLine (unused by this gxx)
//! ```
//! The OCCT abstract base `ApproxInt_SvSurfaces` maps to
//! [`brep_approx::SvSurfaces`]; its base-class field `myUseSolver`
//! (ApproxInt_SvSurfaces.hxx L96) becomes the `my_use_solver` state forwarded
//! by the trait's set_use_solver / get_use_solver virtuals.

use std::cell::RefCell;

use glam::{DVec2, DVec3};
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};
use rcad_kernel::precision::{ANGULAR, APPROXIMATION, CONFUSION};

use crate::geomalgo::brep_approx::{surface_tool, SvSurfaces};
use crate::geomalgo::int_surf::{PntOn2S, Quadric, QuadricType};
use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

// OCCT IntImp_ZerImpFunc.gxx L15-17:
//   #define EpsAng 1.e-8 / EpsAng2 1.e-16 / Tolpetit 1.e-16
// (EpsAng itself is defined but never used in the gxx body — kept verbatim.)
#[allow(dead_code)]
const EPS_ANG: f64 = 1.0e-8;
const EPS_ANG2: f64 = 1.0e-16;
const TOLPETIT: f64 = 1.0e-16;

// ---------------------------------------------------------------------------
// ZerImpFunc — OCCT IntImp_ZerImpFunc (IntImp_ZerImpFunc.gxx / .lxx)
// ---------------------------------------------------------------------------

/// OCCT IntImp_ZerImpFunc (BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox)
/// — the function F(u, v) = ISurfaceTool::Value(ISurface, PSurface(u, v))
/// whose zero set is the intersection of the implicit quadric with the
/// parametric surface; derived from math_FunctionSetWithDerivatives
/// (BRepApprox_TheZerImpFuncOfTheImpPrmSvSurfacesOfApprox.hxx L38).
pub struct ZerImpFunc<'a> {
    /// OCCT void* surf (gxx L24 `#define SURF (*((ThePSurface*)(surf)))`) —
    /// the rcad typed-reference encoding of the OCCT void*.
    surf: Option<&'a BRepAdaptorSurface<'a>>,
    /// OCCT void* func (gxx L25 `#define FUNC (*((TheISurface*)(func)))`).
    func: Option<&'a Quadric>,
    /// OCCT u.
    u: f64,
    /// OCCT v.
    v: f64,
    /// OCCT tol.
    tol: f64,
    /// OCCT pntsol (gp_Pnt).
    pntsol: DVec3,
    /// OCCT valf.
    valf: f64,
    /// OCCT computed.
    computed: bool,
    /// OCCT tangent.
    tangent: bool,
    /// OCCT tgdu.
    tgdu: f64,
    /// OCCT tgdv.
    tgdv: f64,
    /// OCCT gradient (gp_Vec).
    gradient: DVec3,
    /// OCCT derived.
    derived: bool,
    /// OCCT d1u (gp_Vec).
    d1u: DVec3,
    /// OCCT d1v (gp_Vec).
    d1v: DVec3,
    /// OCCT d3d (gp_Vec).
    d3d: DVec3,
    /// OCCT d2d (gp_Dir2d — a normalized 2d direction).
    d2d: DVec2,
}

impl<'a> ZerImpFunc<'a> {
    /// OCCT IntImp_ZerImpFunc() (gxx L27-40).
    pub fn new() -> Self {
        ZerImpFunc {
            surf: None, // OCCT surf(nullptr)
            func: None, // OCCT func(nullptr)
            u: 0.0,
            v: 0.0,
            tol: 0.0,
            pntsol: DVec3::ZERO, // gp_Pnt default
            valf: 0.0,
            computed: false,
            tangent: false,
            tgdu: 0.0,
            tgdv: 0.0,
            gradient: DVec3::ZERO, // gp_Vec default
            derived: false,
            d1u: DVec3::ZERO,
            d1v: DVec3::ZERO,
            d3d: DVec3::ZERO,
            d2d: DVec2::ZERO, // gp_Dir2d default
        }
    }

    /// OCCT IntImp_ZerImpFunc(const ThePSurface& PS, const TheISurface& IS)
    /// (gxx L42-55).
    pub fn new_ps_is(ps: &'a BRepAdaptorSurface<'a>, is: &'a Quadric) -> Self {
        let mut f = ZerImpFunc::new();
        // OCCT L53-54: surf = (void*)(&PS); func = (void*)(&IS);
        f.surf = Some(ps);
        f.func = Some(is);
        f
    }

    /// OCCT IntImp_ZerImpFunc(const TheISurface& IS) (gxx L57-70).
    pub fn new_is(is: &'a Quadric) -> Self {
        let mut f = ZerImpFunc::new();
        // OCCT L69: func = (void*)(&IS);
        f.func = Some(is);
        f
    }

    /// OCCT IntImp_ZerImpFunc::Set(const ThePSurface& PS)
    /// (IntImp_ZerImpFunc.lxx L17-20).
    pub fn set_ps(&mut self, ps: &'a BRepAdaptorSurface<'a>) {
        self.surf = Some(ps);
    }

    /// OCCT IntImp_ZerImpFunc::SetImplicitSurface(const TheISurface& IS)
    /// (lxx L22-25).
    pub fn set_implicit_surface(&mut self, is: &'a Quadric) {
        self.func = Some(is);
    }

    /// OCCT IntImp_ZerImpFunc::Set(const double Tol) (lxx L27-30).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tol = tol;
    }

    /// OCCT IntImp_ZerImpFunc::Root() (lxx L32-35).
    pub fn root(&self) -> f64 {
        self.valf
    }

    /// OCCT IntImp_ZerImpFunc::Tolerance() (lxx L37-40).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// OCCT IntImp_ZerImpFunc::Point() (lxx L42-45).
    pub fn point(&self) -> DVec3 {
        self.pntsol
    }

    /// OCCT IntImp_ZerImpFunc::IsTangent() (gxx L121-150).
    pub fn is_tangent(&mut self) -> bool {
        if !self.computed {
            self.computed = true;
            if !self.derived {
                // OCCT L128: ThePSurfaceTool::D1(SURF, u, v, pntsol, d1u, d1v);
                let psurf = self.p_surface();
                let (pntsol, d1u, d1v) = surface_tool::d1(psurf, self.u, self.v);
                self.pntsol = pntsol;
                self.d1u = d1u;
                self.d1v = d1v;
                self.derived = true;
            }

            // OCCT L132-133.
            self.tgdu = self.gradient.dot(self.d1v);
            self.tgdv = -self.gradient.dot(self.d1u);
            let n2grad = self.gradient.length_squared();
            let n2grad_eps_ang2 = n2grad * EPS_ANG2;
            let n2d1u = self.d1u.length_squared();
            let n2d1v = self.d1v.length_squared();
            self.tangent = (self.tgdu * self.tgdu <= n2grad_eps_ang2 * n2d1v)
                && (self.tgdv * self.tgdv <= n2grad_eps_ang2 * n2d1u);
            if !self.tangent {
                // OCCT L141: d3d.SetLinearForm(tgdu, d1u, tgdv, d1v);
                self.d3d = self.tgdu * self.d1u + self.tgdv * self.d1v;
                // OCCT L142: d2d = gp_Dir2d(tgdu, tgdv) — the gp_Dir2d
                // constructor normalizes (and raises on the null vector);
                // rcad normalize_or_zero covers the degenerate case.
                self.d2d = DVec2::new(self.tgdu, self.tgdv).normalize_or_zero();
                if self.d3d.length() <= TOLPETIT {
                    // jag
                    self.tangent = true;
                }
            }
        }
        self.tangent
    }

    /// OCCT IntImp_ZerImpFunc::Direction3d() (lxx L47-52) — throws
    /// StdFail_UndefinedDerivative when IsTangent().
    pub fn direction_3d(&mut self) -> DVec3 {
        if self.is_tangent() {
            panic!("StdFail_UndefinedDerivative");
        }
        self.d3d
    }

    /// OCCT IntImp_ZerImpFunc::Direction2d() (lxx L54-59) — throws
    /// StdFail_UndefinedDerivative when IsTangent().
    pub fn direction_2d(&mut self) -> DVec2 {
        if self.is_tangent() {
            panic!("StdFail_UndefinedDerivative");
        }
        self.d2d
    }

    /// OCCT IntImp_ZerImpFunc::PSurface() (lxx L61-64) — the void* cast
    /// SURF; rcad panics on the null pointer OCCT would dereference blindly.
    pub fn p_surface(&self) -> &'a BRepAdaptorSurface<'a> {
        self.surf
            .expect("IntImp_ZerImpFunc: null PSurface pointer")
    }

    /// OCCT IntImp_ZerImpFunc::ISurface() (lxx L66-69) — the void* cast FUNC.
    pub fn i_surface(&self) -> &'a Quadric {
        self.func
            .expect("IntImp_ZerImpFunc: null ISurface pointer")
    }
}

impl Default for ZerImpFunc<'_> {
    fn default() -> Self {
        ZerImpFunc::new()
    }
}

impl FunctionSetWithDerivatives for ZerImpFunc<'_> {
    /// OCCT IntImp_ZerImpFunc::NbVariables() (gxx L72-75).
    fn nb_variables(&self) -> usize {
        2
    }

    /// OCCT IntImp_ZerImpFunc::NbEquations() (gxx L77-80).
    fn nb_equations(&self) -> usize {
        1
    }

    /// OCCT IntImp_ZerImpFunc::Value(X, F) (gxx L82-92).
    /// math_Vector X(1, 2) / F(1, 1): X(1)->x[0], X(2)->x[1], F(1)->f[0].
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        self.u = x[0];
        self.v = x[1];
        // OCCT L86: pntsol = ThePSurfaceTool::Value(SURF, u, v);
        let psurf = self.p_surface();
        self.pntsol = surface_tool::value(psurf, self.u, self.v);
        // OCCT L87: valf = TheISurfaceTool::Value(FUNC, pntsol.X(), .Y(), .Z());
        let func = self.i_surface();
        self.valf = func.tool_value(self.pntsol.x, self.pntsol.y, self.pntsol.z);
        f[0] = self.valf;
        self.computed = false;
        self.derived = false;
        true
    }

    /// OCCT IntImp_ZerImpFunc::Derivatives(X, D) (gxx L94-105).
    /// math_Matrix D(1, 1, 1, 2): D(1,1)->df[0][0], D(1,2)->df[0][1].
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        self.u = x[0];
        self.v = x[1];
        // OCCT L98: ThePSurfaceTool::D1(SURF, u, v, pntsol, d1u, d1v);
        let psurf = self.p_surface();
        let (pntsol, d1u, d1v) = surface_tool::d1(psurf, self.u, self.v);
        self.pntsol = pntsol;
        self.d1u = d1u;
        self.d1v = d1v;
        // OCCT L99: TheISurfaceTool::Gradient(FUNC, pntsol.X(), .Y(), .Z(), gradient);
        let func = self.i_surface();
        self.gradient = func.tool_gradient(self.pntsol.x, self.pntsol.y, self.pntsol.z);
        df[0][0] = self.d1u.dot(self.gradient);
        df[0][1] = self.d1v.dot(self.gradient);
        self.computed = false;
        self.derived = true;
        true
    }

    /// OCCT IntImp_ZerImpFunc::Values(X, F, D) (gxx L107-119).
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        self.u = x[0];
        self.v = x[1];
        // OCCT L111: ThePSurfaceTool::D1(SURF, u, v, pntsol, d1u, d1v);
        let psurf = self.p_surface();
        let (pntsol, d1u, d1v) = surface_tool::d1(psurf, self.u, self.v);
        self.pntsol = pntsol;
        self.d1u = d1u;
        self.d1v = d1v;
        // OCCT L112: TheISurfaceTool::ValueAndGradient(FUNC, ..., valf, gradient);
        let func = self.i_surface();
        let (valf, gradient) = func.tool_value_and_gradient(self.pntsol.x, self.pntsol.y, self.pntsol.z);
        self.valf = valf;
        self.gradient = gradient;
        f[0] = self.valf;
        df[0][0] = self.d1u.dot(self.gradient);
        df[0][1] = self.d1v.dot(self.gradient);
        self.computed = false;
        self.derived = true;
        true
    }
}

// ---------------------------------------------------------------------------
// ApproxInt_ImpPrmSvSurfaces.gxx file statics
// ---------------------------------------------------------------------------

/// OCCT ApproxInt_ImpPrmSvSurfaces.gxx L32-57 (static IsSingular) — returns
/// TRUE if vectors theDU || theDV or if at least one of them has
/// null-magnitude; theSqLinTol is the square of the linear tolerance.
fn is_singular(the_du: DVec3, the_dv: DVec3, the_sq_lin_tol: f64, the_ang_tol: f64) -> bool {
    let mut a_du = the_du;
    let mut a_dv = the_dv;

    let a_sq_magn_du = a_du.length_squared();
    let a_sq_magn_dv = a_dv.length_squared();

    if a_sq_magn_du < the_sq_lin_tol {
        return true;
    }

    // OCCT L44: aDU.Divide(sqrt(aSqMagnDU));
    a_du /= a_sq_magn_du.sqrt();

    if a_sq_magn_dv < the_sq_lin_tol {
        return true;
    }

    // OCCT L49: aDV.Divide(sqrt(aSqMagnDV));
    a_dv /= a_sq_magn_dv.sqrt();

    // Here aDU and aDV vectors have magnitude 1.0.

    if a_du.cross(a_dv).length_squared() < the_ang_tol * the_ang_tol {
        return true;
    }

    false
}

/// OCCT ApproxInt_ImpPrmSvSurfaces.gxx L80-268 (static SingularProcessing) —
/// computes the 2D-representation (in UV-coordinates) of theTg3D on the
/// surface in case when theDU.Crossed(theDV).Magnitude() == 0.0; theLinTol is
/// the SQUARE of the tolerance.
#[allow(clippy::too_many_arguments)]
fn singular_processing(
    the_du: DVec3,
    the_dv: DVec3,
    the_is_to_3d_tg_compute: bool,
    the_lin_tol: f64,
    the_ang_tol: f64,
    the_tg3d: &mut DVec3,
    the_tg2d: &mut DVec2,
) -> bool {
    // Attention: @ \sin theAngTol \approx theAngTol @ (for cross-product)

    // Really, vector theTg3D has to be normalized (if theIsTo3DTgCompute == FALSE).
    let a_sq_tan = the_tg3d.length_squared();

    let a_sq_magn_du = the_du.length_squared();
    let a_sq_magn_dv = the_dv.length_squared();

    // There are some reasons of singularity

    // 1.
    if (a_sq_magn_du < the_lin_tol) && (a_sq_magn_dv < the_lin_tol) {
        // For future, this case can be processed as same as in case of
        // osculating surfaces (expanding in Taylor series). Here,
        // we return only.

        return false;
    }

    // 2.
    if a_sq_magn_du < the_lin_tol {
        // In this case, theTg3D vector will be parallel with theDV.
        // Its true direction shall be precised later (the algorithm is
        // based on array of Walking-points).

        if the_is_to_3d_tg_compute {
            // theTg3D will be normalized. Its magnitude is
            let a_tg_magn = 1.0;

            // OCCT L119-121: aNorm = sqrt(aSqMagnDV); theTg3D = theDV.Divided(aNorm);
            let a_norm = a_sq_magn_dv.sqrt();
            *the_tg3d = the_dv / a_norm;
            *the_tg2d = DVec2::new(0.0, a_tg_magn / a_norm);
        } else {
            // theTg3D is already defined.
            // Here we check only, if this tangent is parallel to theDV.

            if the_dv.cross(*the_tg3d).length_squared()
                < the_ang_tol * the_ang_tol * a_sq_magn_dv * a_sq_tan
            {
                // theTg3D is parallel to theDV

                // Use sign "+" if theTg3D and theDV are codirectional
                // and sign "-" if opposite
                let a_dp = the_tg3d.dot(the_dv);
                *the_tg2d = DVec2::new(0.0, (a_sq_tan / a_sq_magn_dv).sqrt().copysign(a_dp));
            } else {
                // theTg3D is not parallel to theDV
                // It is abnormal

                return false;
            }
        }

        return true;
    }

    // 3.
    if a_sq_magn_dv < the_lin_tol {
        // In this case, theTg3D vector will be parallel with theDU.
        // Its true direction shall be precised later (the algorithm is
        // based on array of Walking-points).

        if the_is_to_3d_tg_compute {
            // theTg3D will be normalized. Its magnitude is
            let a_tg_magn = 1.0;

            // OCCT L161-163: aNorm = sqrt(aSqMagnDU); theTg3D = theDU.Divided(aNorm);
            let a_norm = a_sq_magn_du.sqrt();
            *the_tg3d = the_du / a_norm;
            *the_tg2d = DVec2::new(a_tg_magn / a_norm, 0.0);
        } else {
            // theTg3D is already defined.
            // Here we check only, if this tangent is parallel to theDU.

            if the_du.cross(*the_tg3d).length_squared()
                < the_ang_tol * the_ang_tol * a_sq_magn_du * a_sq_tan
            {
                // theTg3D is parallel to theDU

                // Use sign "+" if theTg3D and theDU are codirectional
                // and sign "-" if opposite
                let a_dp = the_tg3d.dot(the_du);
                *the_tg2d = DVec2::new((a_sq_tan / a_sq_magn_du).sqrt().copysign(a_dp), 0.0);
            } else {
                // theTg3D is not parallel to theDU
                // It is abnormal

                return false;
            }
        }

        return true;
    }

    // 4. If aSqMagnDU > 0.0 && aSqMagnDV > 0.0 but theDV || theDU.

    // OCCT L193: const double aLenU = sqrt(aSqMagnDU), aLenV = sqrt(aSqMagnDV);
    let a_len_u = a_sq_magn_du.sqrt();
    let a_len_v = a_sq_magn_dv.sqrt();

    // aLenSum > 0.0 definitely
    let a_len_sum = a_len_u + a_len_v;

    if the_dv.dot(the_du) > 0.0 {
        // Vectors theDV and theDU are codirectional.

        if the_is_to_3d_tg_compute {
            *the_tg2d = DVec2::new(1.0 / a_len_sum, 1.0 / a_len_sum);
            *the_tg3d = the_du * the_tg2d.x + the_dv * the_tg2d.y;
        } else {
            // theTg3D is already defined.
            // Here we check only, if this tangent is parallel to theDU
            //(and theDV together).

            if the_du.cross(*the_tg3d).length_squared()
                < the_ang_tol * the_ang_tol * a_sq_magn_du * a_sq_tan
            {
                // theTg3D is parallel to theDU

                let a_dp = the_tg3d.dot(the_du);
                let a_len_tg = a_sq_tan.sqrt().copysign(a_dp);
                *the_tg2d = DVec2::new(a_len_tg / a_len_sum, a_len_tg / a_len_sum);
            } else {
                // theTg3D is not parallel to theDU
                // It is abnormal

                return false;
            }
        }
    } else {
        // Vectors theDV and theDU are opposite.

        if the_is_to_3d_tg_compute {
            // Here we chose theDU as direction of theTg3D.
            // True direction shall be precised later (the algorithm is
            // based on array of Walking-points).

            *the_tg2d = DVec2::new(1.0 / a_len_sum, -1.0 / a_len_sum);
            *the_tg3d = the_du * the_tg2d.x + the_dv * the_tg2d.y;
        } else {
            // theTg3D is already defined.
            // Here we check only, if this tangent is parallel to theDU
            //(and theDV together).

            if the_du.cross(*the_tg3d).length_squared()
                < the_ang_tol * the_ang_tol * a_sq_magn_du * a_sq_tan
            {
                // theTg3D is parallel to theDU

                let a_dp = the_tg3d.dot(the_du);
                let a_len_tg = a_sq_tan.sqrt().copysign(a_dp);
                *the_tg2d = DVec2::new(a_len_tg / a_len_sum, -a_len_tg / a_len_sum);
            } else {
                // theTg3D is not parallel to theDU
                // It is abnormal

                return false;
            }
        }
    }

    true
}

/// OCCT ApproxInt_ImpPrmSvSurfaces.gxx L287-322 (static NonSingularProcessing)
/// — computes the 2D-representation (in UV-coordinates) of theTg3D on the
/// surface in case when theDU.Crossed(theDV).Magnitude() > 0.0; theLinTol is
/// the SQUARE of the tolerance.
///
/// approx_int.rs carries a reduced private copy (nonsingular_tangent_2d) for
/// the GeomInt WLineAccess semantics; that one is private to its module, so
/// this file keeps its own full-form 1:1 translation.
fn nonsingular_processing(
    the_du: DVec3,
    the_dv: DVec3,
    the_tg3d: DVec3,
    the_lin_tol: f64,
    the_ang_tol: f64,
    the_tg2d: &mut DVec2,
) -> bool {
    // OCCT L294: const gp_Vec aNormal = theDU.Crossed(theDV);
    let a_normal = the_du.cross(the_dv);
    let a_sq_magn = a_normal.length_squared();

    if is_singular(the_du, the_dv, the_lin_tol, the_ang_tol) {
        // OCCT L299-300: gp_Vec aTg(theTg3D);
        let mut a_tg = the_tg3d;
        return singular_processing(the_du, the_dv, false, the_lin_tol, the_ang_tol, &mut a_tg, the_tg2d);
    }

    // If @\vec{T}=\vec{A}*U+\vec{B}*V@ then

    //  \left\{\begin{matrix}
    //  \vec{A} \times \vec{T} = (\vec{A} \times \vec{B})*V
    //  \vec{B} \times \vec{T} = (\vec{B} \times \vec{A})*U
    //  \end{matrix}\right.

    // From here, values of U and V can be found very easily
    //(if @\left \| \vec{A} \times \vec{B} \right \| > 0.0 @,
    // else it is singular case).

    // OCCT L314-316.
    let a_tg_u = the_tg3d.cross(the_du);
    let a_tg_v = the_tg3d.cross(the_dv);
    let a_delta_u = a_tg_v.length_squared() / a_sq_magn;
    let a_delta_v = a_tg_u.length_squared() / a_sq_magn;

    *the_tg2d = DVec2::new(
        a_delta_u.sqrt().copysign(a_tg_v.dot(a_normal)),
        -a_delta_v.sqrt().copysign(a_tg_u.dot(a_normal)),
    );

    true
}

// ---------------------------------------------------------------------------
// ImpPrmSvSurfaces — OCCT ApproxInt_ImpPrmSvSurfaces
// ---------------------------------------------------------------------------

/// The OCCT member block of ApproxInt_ImpPrmSvSurfaces
/// (BRepApprox_TheImpPrmSvSurfacesOfApprox.hxx L104-122) plus the
/// ApproxInt_SvSurfaces base-class field myUseSolver (hxx L96).  The rcad
/// SvSurfaces trait methods take &self while the OCCT methods are
/// non-const, so the whole member block sits behind a RefCell in
/// [`ImpPrmSvSurfaces`].
struct ImpPrmSvSurfacesData<'a> {
    /// OCCT MyParOnS1 (gp_Pnt2d).
    my_par_on_s1: DVec2,
    /// OCCT MyParOnS2 (gp_Pnt2d).
    my_par_on_s2: DVec2,
    /// OCCT MyPnt (gp_Pnt).
    my_pnt: DVec3,
    /// OCCT MyTguv1 (gp_Vec2d).
    my_tguv1: DVec2,
    /// OCCT MyTguv2 (gp_Vec2d).
    my_tguv2: DVec2,
    /// OCCT MyTg (gp_Vec).
    my_tg: DVec3,
    /// OCCT MyIsTangent.
    my_is_tangent: bool,
    /// OCCT MyHasBeenComputed.
    my_has_been_computed: bool,
    /// OCCT MyParOnS1bis.
    my_par_on_s1_bis: DVec2,
    /// OCCT MyParOnS2bis.
    my_par_on_s2_bis: DVec2,
    /// OCCT MyPntbis.
    my_pnt_bis: DVec3,
    /// OCCT MyTguv1bis.
    my_tguv1_bis: DVec2,
    /// OCCT MyTguv2bis.
    my_tguv2_bis: DVec2,
    /// OCCT MyTgbis (gp_Vec).
    my_tg_bis: DVec3,
    /// OCCT MyIsTangentbis.
    my_is_tangent_bis: bool,
    /// OCCT MyHasBeenComputedbis.
    my_has_been_computed_bis: bool,
    /// OCCT MyImplicitFirst.
    my_implicit_first: bool,
    /// OCCT ApproxInt_SvSurfaces::myUseSolver (base-class field, hxx L96).
    my_use_solver: bool,
    /// OCCT MyZerImpFunc.
    my_zer_imp_func: ZerImpFunc<'a>,
}

impl<'a> ImpPrmSvSurfacesData<'a> {
    /// OCCT ApproxInt_ImpPrmSvSurfaces(const TheISurface& ISurf,
    /// const ThePSurface& PSurf) (gxx L326-336) — the implicit surface is
    /// the first one.
    fn new_implicit_first(i_surf: &'a Quadric, p_surf: &'a BRepAdaptorSurface<'a>) -> Self {
        ImpPrmSvSurfacesData {
            my_par_on_s1: DVec2::ZERO,  // gp_Pnt2d default
            my_par_on_s2: DVec2::ZERO,
            my_pnt: DVec3::ZERO,        // gp_Pnt default
            my_tguv1: DVec2::ZERO,      // gp_Vec2d default
            my_tguv2: DVec2::ZERO,
            my_tg: DVec3::ZERO,         // gp_Vec default
            my_is_tangent: false,
            my_has_been_computed: false,
            my_par_on_s1_bis: DVec2::ZERO,
            my_par_on_s2_bis: DVec2::ZERO,
            my_pnt_bis: DVec3::ZERO,
            my_tguv1_bis: DVec2::ZERO,
            my_tguv2_bis: DVec2::ZERO,
            my_tg_bis: DVec3::ZERO,
            my_is_tangent_bis: false,
            my_has_been_computed_bis: false,
            my_implicit_first: true,
            // OCCT MyZerImpFunc(PSurf, ISurf).
            my_zer_imp_func: ZerImpFunc::new_ps_is(p_surf, i_surf),
            my_use_solver: false, // ApproxInt_SvSurfaces() : myUseSolver(false)
        }
        // OCCT ctor body L335: SetUseSolver(true);
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces(const ThePSurface& PSurf,
    /// const TheISurface& ISurf) (gxx L340-350) — the parametric surface is
    /// the first one.
    fn new_parametric_first(p_surf: &'a BRepAdaptorSurface<'a>, i_surf: &'a Quadric) -> Self {
        ImpPrmSvSurfacesData {
            my_par_on_s1: DVec2::ZERO,
            my_par_on_s2: DVec2::ZERO,
            my_pnt: DVec3::ZERO,
            my_tguv1: DVec2::ZERO,
            my_tguv2: DVec2::ZERO,
            my_tg: DVec3::ZERO,
            my_is_tangent: false,
            my_has_been_computed: false,
            my_par_on_s1_bis: DVec2::ZERO,
            my_par_on_s2_bis: DVec2::ZERO,
            my_pnt_bis: DVec3::ZERO,
            my_tguv1_bis: DVec2::ZERO,
            my_tguv2_bis: DVec2::ZERO,
            my_tg_bis: DVec3::ZERO,
            my_is_tangent_bis: false,
            my_has_been_computed_bis: false,
            my_implicit_first: false,
            my_zer_imp_func: ZerImpFunc::new_ps_is(p_surf, i_surf),
            my_use_solver: false,
        }
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::Pnt (gxx L354-369).
    fn pnt(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, p: &mut DVec3) {
        // OCCT L360-362: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        let mut a_p = DVec3::ZERO;
        let mut a_t = DVec3::ZERO;
        let mut a_ts1 = DVec2::ZERO;
        let mut a_ts2 = DVec2::ZERO;
        // OCCT L363-366.
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        self.compute(
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        // OCCT L368: P = MyPnt;
        *p = self.my_pnt;
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::Tangency (gxx L373-389).
    fn tangency(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, t: &mut DVec3) -> bool {
        // OCCT L379-381: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        let mut a_p = DVec3::ZERO;
        let mut a_t = DVec3::ZERO;
        let mut a_ts1 = DVec2::ZERO;
        let mut a_ts2 = DVec2::ZERO;
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t_flag = self.compute(
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        // OCCT L387-388: T = MyTg; return (t);
        *t = self.my_tg;
        t_flag
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::TangencyOnSurf1 (gxx L393-409).
    fn tangency_on_surf_1(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, t: &mut DVec2) -> bool {
        // OCCT L399-401: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        let mut a_p = DVec3::ZERO;
        let mut a_t = DVec3::ZERO;
        let mut a_ts1 = DVec2::ZERO;
        let mut a_ts2 = DVec2::ZERO;
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t_flag = self.compute(
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        // OCCT L407-408: T = MyTguv1; return (t);
        *t = self.my_tguv1;
        t_flag
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::TangencyOnSurf2 (gxx L413-429).
    fn tangency_on_surf_2(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, t: &mut DVec2) -> bool {
        // OCCT L419-421: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        let mut a_p = DVec3::ZERO;
        let mut a_t = DVec3::ZERO;
        let mut a_ts1 = DVec2::ZERO;
        let mut a_ts2 = DVec2::ZERO;
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t_flag = self.compute(
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        // OCCT L427-428: T = MyTguv2; return (t);
        *t = self.my_tguv2;
        t_flag
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::Compute (gxx L437-769) — computes
    /// point on curve, 3D and 2D-tangents of a curve and parameters on the
    /// surfaces.
    #[allow(clippy::too_many_arguments)]
    fn compute(
        &mut self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        p: &mut DVec3,
        tg: &mut DVec3,
        tguv1: &mut DVec2,
        tguv2: &mut DVec2,
    ) -> bool {
        // OCCT L446-447.
        let a_q_surf: &Quadric = self.my_zer_imp_func.i_surface();
        let a_p_surf: &BRepAdaptorSurface<'a> = self.my_zer_imp_func.p_surface();
        // OCCT L448-449: gp_Vec2d& aQuadTg = MyImplicitFirst ? Tguv1 : Tguv2;
        //                    gp_Vec2d& aPrmTg  = MyImplicitFirst ? Tguv2 : Tguv1;
        // (Rust has no mutable aliases; the selection is re-made at each of
        // the four use sites below.)

        // for square
        let a_null_value = APPROXIMATION * APPROXIMATION;
        let an_ang_tol = ANGULAR;

        let tu1 = *u1;
        let tu2 = *u2;
        let tv1 = *v1;
        let tv2 = *v2;

        if self.my_has_been_computed {
            if (self.my_par_on_s1.x == *u1)
                && (self.my_par_on_s1.y == *v1)
                && (self.my_par_on_s2.x == *u2)
                && (self.my_par_on_s2.y == *v2)
            {
                return self.my_is_tangent;
            } else if self.my_has_been_computed_bis == false {
                self.my_tg_bis = self.my_tg;
                self.my_tguv1_bis = self.my_tguv1;
                self.my_tguv2_bis = self.my_tguv2;
                self.my_pnt_bis = self.my_pnt;
                self.my_par_on_s1_bis = self.my_par_on_s1;
                self.my_par_on_s2_bis = self.my_par_on_s2;
                self.my_is_tangent_bis = self.my_is_tangent;
                self.my_has_been_computed_bis = self.my_has_been_computed;
            }
        }

        if self.my_has_been_computed_bis {
            if (self.my_par_on_s1_bis.x == *u1)
                && (self.my_par_on_s1_bis.y == *v1)
                && (self.my_par_on_s2_bis.x == *u2)
                && (self.my_par_on_s2_bis.y == *v2)
            {
                // OCCT L486-492: TV/TV1/TV2/TP/TP1/TP2/TB snapshots.
                let t_v = self.my_tg;
                let t_v1 = self.my_tguv1;
                let t_v2 = self.my_tguv2;
                let t_p = self.my_pnt;
                let t_p1 = self.my_par_on_s1;
                let t_p2 = self.my_par_on_s2;
                let t_b = self.my_is_tangent;

                self.my_tg = self.my_tg_bis;
                self.my_tguv1 = self.my_tguv1_bis;
                self.my_tguv2 = self.my_tguv2_bis;
                self.my_pnt = self.my_pnt_bis;
                self.my_par_on_s1 = self.my_par_on_s1_bis;
                self.my_par_on_s2 = self.my_par_on_s2_bis;
                self.my_is_tangent = self.my_is_tangent_bis;

                self.my_tg_bis = t_v;
                self.my_tguv1_bis = t_v1;
                self.my_tguv2_bis = t_v2;
                self.my_pnt_bis = t_p;
                self.my_par_on_s1_bis = t_p1;
                self.my_par_on_s2_bis = t_p2;
                self.my_is_tangent_bis = t_b;

                return self.my_is_tangent;
            }
        }

        // OCCT L514-515: math_Vector X(1, 2), BornInf(1, 2), BornSup(1, 2),
        // Tolerance(1, 2); the rcad [f64; 2] slots map X(1)->[0], X(2)->[1].
        let mut x = [0.0f64; 2];
        let mut born_inf = [0.0f64; 2];
        let mut born_sup = [0.0f64; 2];
        let mut tolerance = [0.0f64; 2];
        //--- ThePSurfaceTool::GetResolution(aPSurf,Tolerance(1),Tolerance(2));
        tolerance[0] = 1.0e-8;
        tolerance[1] = 1.0e-8;
        let binfu = surface_tool::first_u_parameter(a_p_surf);
        let binfv = surface_tool::first_v_parameter(a_p_surf);
        let bsupu = surface_tool::last_u_parameter(a_p_surf);
        let bsupv = surface_tool::last_v_parameter(a_p_surf);
        born_inf[0] = binfu;
        born_sup[0] = bsupu;
        born_inf[1] = binfv;
        born_sup[1] = bsupv;
        let mut translation_u = 0.0;
        let mut translation_v = 0.0;

        if !self.fill_initial_vector_of_solution(
            *u1,
            *v1,
            *u2,
            *v2,
            binfu,
            bsupu,
            binfv,
            bsupv,
            &mut x,
            &mut translation_u,
            &mut translation_v,
        ) {
            // OCCT L542-543: MyIsTangent = MyIsTangentbis = false;
            //                   MyHasBeenComputed = MyHasBeenComputedbis = false;
            self.my_is_tangent = false;
            self.my_is_tangent_bis = false;
            self.my_has_been_computed = false;
            self.my_has_been_computed_bis = false;
            return false;
        }

        // OCCT L547-549.
        let mut a_rsnld_is_done = false;
        let pour_tester_u = x[0];
        let pour_tester_v = x[1];
        if self.my_use_solver {
            // OCCT L552-554: math_FunctionSetRoot Rsnld(MyZerImpFunc) — the
            // 2-arg ctor (NbIterations default 100) leaves the tolerances
            // unset; SetTolerance installs them before Perform.
            let mut rsnld = FunctionSetRoot::new(&self.my_zer_imp_func, &[0.0, 0.0], 100);
            rsnld.set_tolerance(&tolerance);
            rsnld.perform(&mut self.my_zer_imp_func, &x, &born_inf, &born_sup, false);
            a_rsnld_is_done = rsnld.is_done();
            if a_rsnld_is_done {
                // OCCT L557: Rsnld.Root(X);
                x.copy_from_slice(&rsnld.root());
            }
        }
        if a_rsnld_is_done || !self.my_use_solver {
            self.my_has_been_computed = true;

            // OCCT L563-564.
            let dist_avant_apres_u = (pour_tester_u - x[0]).abs();
            let dist_avant_apres_v = (pour_tester_v - x[1]).abs();

            // OCCT L566: MyPnt = P = ThePSurfaceTool::Value(aPSurf, X(1), X(2));
            let pnt_val = surface_tool::value(a_p_surf, x[0], x[1]);
            self.my_pnt = pnt_val;
            *p = pnt_val;

            if (dist_avant_apres_v <= 0.001) && (dist_avant_apres_u <= 0.001) {
                // OCCT L570-571: gp_Vec aD1uPrm, aD1vPrm; gp_Vec aD1uQuad, aD1vQuad;
                let a_d1u_prm: DVec3;
                let a_d1v_prm: DVec3;
                let a_d1u_quad: DVec3;
                let a_d1v_quad: DVec3;

                if self.my_implicit_first {
                    *u2 = x[0] - translation_u;
                    *v2 = x[1] - translation_v;

                    if a_q_surf.type_quadric() != QuadricType::Plane {
                        while *u1 - tu1 > std::f64::consts::PI {
                            *u1 -= std::f64::consts::PI + std::f64::consts::PI;
                        }
                        while tu1 - *u1 > std::f64::consts::PI {
                            *u1 += std::f64::consts::PI + std::f64::consts::PI;
                        }
                    }

                    self.my_par_on_s1 = DVec2::new(tu1, tv1);
                    self.my_par_on_s2 = DVec2::new(tu2, tv2);

                    // OCCT L589-591: gp_Pnt aP2;
                    //   ThePSurfaceTool::D1(aPSurf, X(1), X(2), P, aD1uPrm, aD1vPrm);
                    let (p_d1, d1u_prm, d1v_prm) = surface_tool::d1(a_p_surf, x[0], x[1]);
                    *p = p_d1;
                    a_d1u_prm = d1u_prm;
                    a_d1v_prm = d1v_prm;
                    // OCCT L592: aQSurf.D1(u1, v1, aP2, aD1uQuad, aD1vQuad);
                    let (a_p2, d1u_quad, d1v_quad) = a_q_surf.d1(*u1, *v1);
                    a_d1u_quad = d1u_quad;
                    a_d1v_quad = d1v_quad;

                    // Middle-point of P-P2 segment
                    // OCCT L595: P.BaryCenter(1.0, aP2, 1.0);
                    *p = (1.0 * *p + 1.0 * a_p2) / (1.0 + 1.0);
                } else {
                    *u1 = x[0] - translation_u;
                    *v1 = x[1] - translation_v;
                    // aQSurf.Parameters(P, u2, v2);
                    if a_q_surf.type_quadric() != QuadricType::Plane {
                        while *u2 - tu2 > std::f64::consts::PI {
                            *u2 -= std::f64::consts::PI + std::f64::consts::PI;
                        }
                        while tu2 - *u2 > std::f64::consts::PI {
                            *u2 += std::f64::consts::PI + std::f64::consts::PI;
                        }
                    }

                    self.my_par_on_s1 = DVec2::new(tu1, tv1);
                    // OCCT L611: MyParOnS2.SetCoord(tu2, tu2); — the OCCT
                    // literal (the second coordinate is tu2, not tv2).
                    self.my_par_on_s2 = DVec2::new(tu2, tu2);

                    // OCCT L613-614: gp_Pnt aP2;
                    //   ThePSurfaceTool::D1(aPSurf, X(1), X(2), P, aD1uPrm, aD1vPrm);
                    let (p_d1, d1u_prm, d1v_prm) = surface_tool::d1(a_p_surf, x[0], x[1]);
                    *p = p_d1;
                    a_d1u_prm = d1u_prm;
                    a_d1v_prm = d1v_prm;

                    // OCCT L616: aQSurf.D1(u2, v2, aP2, aD1uQuad, aD1vQuad);
                    let (a_p2, d1u_quad, d1v_quad) = a_q_surf.d1(*u2, *v2);
                    a_d1u_quad = d1u_quad;
                    a_d1v_quad = d1v_quad;

                    // Middle-point of P-P2 segment
                    // OCCT L619: P.BaryCenter(1.0, aP2, 1.0);
                    *p = (1.0 * *p + 1.0 * a_p2) / (1.0 + 1.0);
                }

                self.my_pnt = *p;

                // Normals to the surfaces
                // OCCT L625.
                let mut a_normal_prm = a_d1u_prm.cross(a_d1v_prm);
                let mut a_normal_imp = a_q_surf.normale(self.my_pnt);

                let a_sq_magn_prm = a_normal_prm.length_squared();
                let a_sq_magn_imp = a_normal_imp.length_squared();

                let mut is_prm_singular = false;
                let mut is_imp_singular = false;

                if is_singular(a_d1u_prm, a_d1v_prm, a_null_value, an_ang_tol) {
                    is_prm_singular = true;
                    // OCCT L635: SingularProcessing(aD1uPrm, aD1vPrm, true,
                    //   aNullValue, anAngTol, Tg, aPrmTg);
                    // aPrmTg = MyImplicitFirst ? Tguv2 : Tguv1
                    let a_prm_tg = if self.my_implicit_first { &mut *tguv2 } else { &mut *tguv1 };
                    if !singular_processing(
                        a_d1u_prm,
                        a_d1v_prm,
                        true,
                        a_null_value,
                        an_ang_tol,
                        tg,
                        a_prm_tg,
                    ) {
                        self.my_is_tangent = false;
                        self.my_is_tangent_bis = false;
                        self.my_has_been_computed = false;
                        self.my_has_been_computed_bis = false;
                        return false;
                    }

                    self.my_tg = *tg;
                } else {
                    // OCCT L646: aNormalPrm.Divide(sqrt(aSQMagnPrm));
                    a_normal_prm /= a_sq_magn_prm.sqrt();
                }

                // Analogically for implicit surface
                if a_sq_magn_imp < a_null_value {
                    is_imp_singular = true;

                    // OCCT L654-660: SingularProcessing(aD1uQuad, aD1vQuad,
                    //   !isPrmSingular, aNullValue, anAngTol, Tg, aQuadTg);
                    // aQuadTg = MyImplicitFirst ? Tguv1 : Tguv2
                    let a_quad_tg = if self.my_implicit_first { &mut *tguv1 } else { &mut *tguv2 };
                    if !singular_processing(
                        a_d1u_quad,
                        a_d1v_quad,
                        !is_prm_singular,
                        a_null_value,
                        an_ang_tol,
                        tg,
                        a_quad_tg,
                    ) {
                        self.my_is_tangent = false;
                        self.my_is_tangent_bis = false;
                        self.my_has_been_computed = false;
                        self.my_has_been_computed_bis = false;
                        return false;
                    }

                    self.my_tg = *tg;
                } else {
                    // OCCT L671: aNormalImp.Divide(sqrt(aSQMagnImp));
                    a_normal_imp /= a_sq_magn_imp.sqrt();
                }

                if is_imp_singular && is_prm_singular {
                    // All is OK. All abnormal cases were processed above.

                    self.my_tguv1 = *tguv1;
                    self.my_tguv2 = *tguv2;

                    self.my_is_tangent = true;
                    return self.my_is_tangent;
                } else if !(is_imp_singular || is_prm_singular) {
                    // Processing pure non-singular case
                    //(3D- and 2D-tangents are still not defined)

                    // Ask to pay attention to the fact that here
                    // aNormalImp and aNormalPrm are normalized.
                    // Therefore, @ \left \| \vec{Tg} \right \| = 0.0 @
                    // if and only if (aNormalImp || aNormalPrm).
                    *tg = a_normal_imp.cross(a_normal_prm);
                }

                let a_sq_magn_tg = tg.length_squared();

                if a_sq_magn_tg < a_null_value {
                    self.my_is_tangent = false;
                    self.my_is_tangent_bis = false;
                    self.my_has_been_computed = false;
                    self.my_has_been_computed_bis = false;
                    return false;
                }

                // Normalize Tg vector
                *tg /= a_sq_magn_tg.sqrt();
                self.my_tg = *tg;

                if !is_prm_singular {
                    // If isPrmSingular==TRUE then aPrmTg has already been computed.

                    // OCCT L713: NonSingularProcessing(aD1uPrm, aD1vPrm, Tg,
                    //   aNullValue, anAngTol, aPrmTg);
                    let a_prm_tg = if self.my_implicit_first { &mut *tguv2 } else { &mut *tguv1 };
                    if !nonsingular_processing(a_d1u_prm, a_d1v_prm, *tg, a_null_value, an_ang_tol, a_prm_tg)
                    {
                        self.my_is_tangent = false;
                        self.my_is_tangent_bis = false;
                        self.my_has_been_computed = false;
                        self.my_has_been_computed_bis = false;
                        return false;
                    }
                }

                if !is_imp_singular {
                    // If isImpSingular==TRUE then aQuadTg has already been computed.

                    // OCCT L725: NonSingularProcessing(aD1uQuad, aD1vQuad, Tg,
                    //   aNullValue, anAngTol, aQuadTg);
                    let a_quad_tg = if self.my_implicit_first { &mut *tguv1 } else { &mut *tguv2 };
                    if !nonsingular_processing(a_d1u_quad, a_d1v_quad, *tg, a_null_value, an_ang_tol, a_quad_tg)
                    {
                        self.my_is_tangent = false;
                        self.my_is_tangent_bis = false;
                        self.my_has_been_computed = false;
                        self.my_has_been_computed_bis = false;
                        return false;
                    }
                }

                self.my_tguv1 = *tguv1;
                self.my_tguv2 = *tguv2;

                self.my_is_tangent = true;

                // OCCT L738-750: the #ifdef OCCT_DEBUG dump block is not
                // translated.

                return true;
            } else {
                //-- cout<<" ApproxInt_ImpImpSvSurfaces.gxx : Distance apres recadrage Trop Grande "<<endl;

                self.my_is_tangent = false;
                self.my_is_tangent_bis = false;
                self.my_has_been_computed = false;
                self.my_has_been_computed_bis = false;
                return false;
            }
        } else {
            self.my_is_tangent = false;
            self.my_is_tangent_bis = false;
            self.my_has_been_computed = false;
            self.my_has_been_computed_bis = false;
            return false;
        }
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::SeekPoint (gxx L777-861) — computes
    /// point on curve and parameters on the surfaces.
    // The OCCT locals NewU1..NewV2 (gxx L815) are uninitialized and only
    // written inside the IsDone branch; the zero initializers are therefore
    // dead stores — kept with #[allow(unused_assignments)].
    #[allow(unused_assignments)]
    fn seek_point(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, point: &mut PntOn2S) -> bool {
        // OCCT L783-784.
        let a_q_surf: &Quadric = self.my_zer_imp_func.i_surface();
        let a_p_surf: &BRepAdaptorSurface<'a> = self.my_zer_imp_func.p_surface();

        let mut x = [0.0f64; 2];
        let mut born_inf = [0.0f64; 2];
        let mut born_sup = [0.0f64; 2];
        let mut tolerance = [0.0f64; 2];
        //--- ThePSurfaceTool::GetResolution(aPSurf,Tolerance(1),Tolerance(2));
        tolerance[0] = 1.0e-8;
        tolerance[1] = 1.0e-8;
        let binfu = surface_tool::first_u_parameter(a_p_surf);
        let binfv = surface_tool::first_v_parameter(a_p_surf);
        let bsupu = surface_tool::last_u_parameter(a_p_surf);
        let bsupv = surface_tool::last_v_parameter(a_p_surf);
        born_inf[0] = binfu;
        born_sup[0] = bsupu;
        born_inf[1] = binfv;
        born_sup[1] = bsupv;
        let mut translation_u = 0.0;
        let mut translation_v = 0.0;

        if !self.fill_initial_vector_of_solution(
            u1,
            v1,
            u2,
            v2,
            binfu,
            bsupu,
            binfv,
            bsupv,
            &mut x,
            &mut translation_u,
            &mut translation_v,
        ) {
            return false;
        }

        // OCCT L815: double NewU1, NewV1, NewU2, NewV2; — uninitialized OCCT
        // locals; neutral zero initialization.
        let mut new_u1 = 0.0;
        let mut new_v1 = 0.0;
        let mut new_u2 = 0.0;
        let mut new_v2 = 0.0;

        // OCCT L817-819: math_FunctionSetRoot Rsnld(MyZerImpFunc);
        //   Rsnld.SetTolerance(Tolerance); Rsnld.Perform(MyZerImpFunc, X, ...);
        let mut rsnld = FunctionSetRoot::new(&self.my_zer_imp_func, &[0.0, 0.0], 100);
        rsnld.set_tolerance(&tolerance);
        rsnld.perform(&mut self.my_zer_imp_func, &x, &born_inf, &born_sup, false);
        if rsnld.is_done() {
            self.my_has_been_computed = true;
            // OCCT L823: Rsnld.Root(X);
            x.copy_from_slice(&rsnld.root());

            // OCCT L825: MyPnt = ThePSurfaceTool::Value(aPSurf, X(1), X(2));
            self.my_pnt = surface_tool::value(a_p_surf, x[0], x[1]);

            if self.my_implicit_first {
                new_u2 = x[0] - translation_u;
                new_v2 = x[1] - translation_v;

                // OCCT L832: aQSurf.Parameters(MyPnt, NewU1, NewV1);
                let (q_u, q_v) = a_q_surf.parameters(self.my_pnt);
                new_u1 = q_u;
                new_v1 = q_v;
                // adjust U
                if a_q_surf.type_quadric() != QuadricType::Plane {
                    let sign: f64 = if new_u1 > u1 { -1.0 } else { 1.0 };
                    while (u1 - new_u1).abs() > std::f64::consts::PI {
                        new_u1 += sign * (std::f64::consts::PI + std::f64::consts::PI);
                    }
                }
            } else {
                new_u1 = x[0] - translation_u;
                new_v1 = x[1] - translation_v;

                // OCCT L846: aQSurf.Parameters(MyPnt, NewU2, NewV2);
                let (q_u, q_v) = a_q_surf.parameters(self.my_pnt);
                new_u2 = q_u;
                new_v2 = q_v;
                // adjust U
                if a_q_surf.type_quadric() != QuadricType::Plane {
                    let sign: f64 = if new_u2 > u2 { -1.0 } else { 1.0 };
                    while (u2 - new_u2).abs() > std::f64::consts::PI {
                        new_u2 += sign * (std::f64::consts::PI + std::f64::consts::PI);
                    }
                }
            }
        } else {
            return false;
        }

        // OCCT L859: Point.SetValue(MyPnt, NewU1, NewV1, NewU2, NewV2);
        point.set_value_all(self.my_pnt, new_u1, new_v1, new_u2, new_v2);
        true
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::FillInitialVectorOfSolution
    /// (gxx L865-1020).
    #[allow(clippy::too_many_arguments)]
    fn fill_initial_vector_of_solution(
        &mut self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        binfu: f64,
        bsupu: f64,
        binfv: f64,
        bsupv: f64,
        x: &mut [f64],
        translation_u: &mut f64,
        translation_v: &mut f64,
    ) -> bool {
        // OCCT L879: math_Vector F(1, 1); — declared in the OCCT body but
        // never used afterwards.

        // OCCT L877: const ThePSurface& aPSurf = MyZerImpFunc.PSurface();
        let a_p_surf: &BRepAdaptorSurface<'a> = self.my_zer_imp_func.p_surface();

        *translation_u = 0.0;
        *translation_v = 0.0;

        if self.my_implicit_first {
            if u2 < binfu - 0.0000000001 {
                if surface_tool::is_u_periodic(a_p_surf) {
                    let d = surface_tool::u_period(a_p_surf);
                    loop {
                        *translation_u += d;
                        if !(u2 + *translation_u < binfu) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if u2 > bsupu + 0.0000000001 {
                if surface_tool::is_u_periodic(a_p_surf) {
                    let d = surface_tool::u_period(a_p_surf);
                    loop {
                        *translation_u -= d;
                        if !(u2 + *translation_u > bsupu) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            if v2 < binfv - 0.0000000001 {
                if surface_tool::is_v_periodic(a_p_surf) {
                    let d = surface_tool::v_period(a_p_surf);
                    loop {
                        *translation_v += d;
                        if !(v2 + *translation_v < binfv) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if v2 > bsupv + 0.0000000001 {
                if surface_tool::is_v_periodic(a_p_surf) {
                    let d = surface_tool::v_period(a_p_surf);
                    loop {
                        *translation_v -= d;
                        if !(v2 + *translation_v > bsupv) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            x[0] = u2 + *translation_u;
            x[1] = v2 + *translation_v;
        } else {
            if u1 < binfu - 0.0000000001 {
                if surface_tool::is_u_periodic(a_p_surf) {
                    let d = surface_tool::u_period(a_p_surf);
                    loop {
                        *translation_u += d;
                        if !(u1 + *translation_u < binfu) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if u1 > bsupu + 0.0000000001 {
                if surface_tool::is_u_periodic(a_p_surf) {
                    let d = surface_tool::u_period(a_p_surf);
                    loop {
                        *translation_u -= d;
                        if !(u1 + *translation_u > bsupu) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            if v1 < binfv - 0.0000000001 {
                if surface_tool::is_v_periodic(a_p_surf) {
                    let d = surface_tool::v_period(a_p_surf);
                    loop {
                        *translation_v += d;
                        if !(v1 + *translation_v < binfv) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if v1 > bsupv + 0.0000000001 {
                if surface_tool::is_v_periodic(a_p_surf) {
                    let d = surface_tool::v_period(a_p_surf);
                    loop {
                        *translation_v -= d;
                        if !(v1 + *translation_v > bsupv) {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            x[0] = u1 + *translation_u;
            x[1] = v1 + *translation_v;
        }

        //=================================================================================================

        // Make a small step from boundaries in order to avoid
        // finding "outboundaried" solution (Rsnld -> NotDone).
        if self.my_use_solver {
            // OCCT L1005-1008: std::max(Precision::Confusion(), ...Resolution(aPSurf, Precision::Confusion()))
            let du = CONFUSION.max(surface_tool::u_resolution(a_p_surf, CONFUSION));
            let dv = CONFUSION.max(surface_tool::v_resolution(a_p_surf, CONFUSION));
            if x[0] - 0.0000000001 <= binfu {
                x[0] += du;
            }
            if x[0] + 0.0000000001 >= bsupu {
                x[0] -= du;
            }
            if x[1] - 0.0000000001 <= binfv {
                x[1] += dv;
            }
            if x[1] + 0.0000000001 >= bsupv {
                x[1] -= dv;
            }
        }

        true
    }
}

/// OCCT BRepApprox_TheImpPrmSvSurfacesOfApprox — the ApproxInt_ImpPrmSvSurfaces
/// instantiation over BRepAdaptor_Surface x IntSurf_Quadric.  The OCCT
/// non-const member functions mutate the member block, which sits behind a
/// RefCell so the &self [`SvSurfaces`] trait methods forward to it.
pub struct ImpPrmSvSurfaces<'a> {
    my: RefCell<ImpPrmSvSurfacesData<'a>>,
}

impl<'a> ImpPrmSvSurfaces<'a> {
    /// OCCT ApproxInt_ImpPrmSvSurfaces(const TheISurface& ISurf,
    /// const ThePSurface& PSurf) (gxx L326-336) — the implicit surface is
    /// the first one; the ctor body calls SetUseSolver(true).
    pub fn new_implicit_first(i_surf: &'a Quadric, p_surf: &'a BRepAdaptorSurface<'a>) -> Self {
        let mut data = ImpPrmSvSurfacesData::new_implicit_first(i_surf, p_surf);
        data.my_use_solver = true; // OCCT L335: SetUseSolver(true);
        ImpPrmSvSurfaces {
            my: RefCell::new(data),
        }
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces(const ThePSurface& PSurf,
    /// const TheISurface& ISurf) (gxx L340-350) — the parametric surface is
    /// the first one; the ctor body calls SetUseSolver(true).
    pub fn new_parametric_first(p_surf: &'a BRepAdaptorSurface<'a>, i_surf: &'a Quadric) -> Self {
        let mut data = ImpPrmSvSurfacesData::new_parametric_first(p_surf, i_surf);
        data.my_use_solver = true; // OCCT L349: SetUseSolver(true);
        ImpPrmSvSurfaces {
            my: RefCell::new(data),
        }
    }
}

impl SvSurfaces for ImpPrmSvSurfaces<'_> {
    /// OCCT ApproxInt_ImpPrmSvSurfaces::Compute (gxx L437-769).
    fn compute(
        &self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        pt: &mut DVec3,
        tg: &mut DVec3,
        tguv1: &mut DVec2,
        tguv2: &mut DVec2,
    ) -> bool {
        self.my.borrow_mut().compute(u1, v1, u2, v2, pt, tg, tguv1, tguv2)
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::Pnt (gxx L354-369).
    fn pnt(&self, u1: f64, v1: f64, u2: f64, v2: f64, p: &mut DVec3) {
        self.my.borrow_mut().pnt(u1, v1, u2, v2, p)
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::SeekPoint (gxx L777-861).
    fn seek_point(&self, u1: f64, v1: f64, u2: f64, v2: f64, point: &mut PntOn2S) -> bool {
        self.my.borrow_mut().seek_point(u1, v1, u2, v2, point)
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::Tangency (gxx L373-389).
    fn tangency(&self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut DVec3) -> bool {
        self.my.borrow_mut().tangency(u1, v1, u2, v2, tg)
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::TangencyOnSurf1 (gxx L393-409).
    fn tangency_on_surf_1(&self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut DVec2) -> bool {
        self.my.borrow_mut().tangency_on_surf_1(u1, v1, u2, v2, tg)
    }

    /// OCCT ApproxInt_ImpPrmSvSurfaces::TangencyOnSurf2 (gxx L413-429).
    fn tangency_on_surf_2(&self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut DVec2) -> bool {
        self.my.borrow_mut().tangency_on_surf_2(u1, v1, u2, v2, tg)
    }

    /// OCCT ApproxInt_SvSurfaces::SetUseSolver (ApproxInt_SvSurfaces.hxx L93)
    /// — writes the base-class myUseSolver field.
    fn set_use_solver(&self, the_use_sol: bool) {
        self.my.borrow_mut().my_use_solver = the_use_sol;
    }

    /// OCCT ApproxInt_SvSurfaces::GetUseSolver (hxx L95).
    fn get_use_solver(&self) -> bool {
        self.my.borrow().my_use_solver
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Line3, Plane, Surface3};
    use rcad_kernel::topo::topods::{BRepBuilder, Shape};

    /// The implicit quadric: the plane z = 0 with X/Y u/v axes.
    fn plane_quadric_z0() -> Quadric {
        Quadric::from_plane(&Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(0.0, 0.0, 1.0),
            u_dir: DVec3::new(1.0, 0.0, 0.0),
            v_dir: DVec3::new(0.0, 1.0, 0.0),
        })
    }

    /// A tilted plane face over the UV window [0, 1] x [0, 1]:
    /// P(u, v) = u * X + v * (Y + Z) / sqrt(2).  Its intersection with the
    /// plane z = 0 is the line v = 0 (3D: z = 0 along X).
    fn build_tilted_plane_brep() -> (rcad_kernel::BRep, Shape) {
        let u_dir = DVec3::new(1.0, 0.0, 0.0);
        let v_dir = DVec3::new(0.0, 1.0, 1.0).normalize();
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, rcad_kernel::geom::Point3::ZERO, 1e-7);
        let v2 = b.add_vertex(&mut brep, u_dir, 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3::new(rcad_kernel::geom::Point3::ZERO, u_dir))),
            v1,
            v2,
            [0.0, 1.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: rcad_kernel::geom::Point3::ZERO,
                normal: u_dir.cross(v_dir).normalize(),
                u_dir,
                v_dir,
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );
        (brep, face)
    }

    /// OCCT anchor: IntImp_ZerImpFunc over BRepAdaptor plane x IntSurf_Quadric
    /// plane (gxx L82-119 Value/Derivatives/Values, L121-150 IsTangent;
    /// lxx L27-69 Set/Root/Tolerance/Point/Direction3d/Direction2d).  The two
    /// planes are non-parallel: F(u, v) = v / sqrt(2), the zero set is the
    /// line v = 0 and the solver converges onto it from (0.3, 0.6).
    #[test]
    fn zer_imp_func_plane_plane_anchor() {
        let i_surf = plane_quadric_z0();
        let (brep, face) = build_tilted_plane_brep();
        let p_surf = BRepAdaptorSurface::initialize_face(&brep, &face, true);

        // OCCT ctor (const ThePSurface& PS, const TheISurface& IS) + Set(Tol).
        let mut func = ZerImpFunc::new_ps_is(&p_surf, &i_surf);
        func.set_tolerance(1.0e-12);
        assert_eq!(func.tolerance(), 1.0e-12);
        assert_eq!(func.nb_variables(), 2);
        assert_eq!(func.nb_equations(), 1);

        // Analytic Value at X = (0.3, 0.6): F = z = v / sqrt(2)
        // (gxx L82-92).
        let x = [0.3f64, 0.6];
        let mut f = [0.0f64; 1];
        assert!(func.value(&x, &mut f));
        let expected_f = 0.6 / 2.0f64.sqrt();
        assert!((f[0] - expected_f).abs() < 1e-14, "F={} vs {expected_f}", f[0]);
        assert_eq!(func.root(), f[0]); // lxx L32-35
        let expected_p = DVec3::new(0.3, expected_f, expected_f);
        assert!((func.point() - expected_p).length() < 1e-14); // lxx L42-45

        // Analytic Derivatives: D = (d1u . N, d1v . N) = (0, 1/sqrt(2))
        // (gxx L94-105).
        let mut d = [vec![0.0f64; 2]; 1];
        assert!(func.derivatives(&x, &mut d));
        assert!(d[0][0].abs() < 1e-15);
        assert!((d[0][1] - 1.0 / 2.0f64.sqrt()).abs() < 1e-14);

        // Analytic Values: F and D in one call (gxx L107-119).
        let mut f2 = [0.0f64; 1];
        let mut d2 = [vec![0.0f64; 2]; 1];
        assert!(func.values(&x, &mut f2, &mut d2));
        assert!((f2[0] - expected_f).abs() < 1e-14);
        assert!(d2[0][0].abs() < 1e-15);
        assert!((d2[0][1] - 1.0 / 2.0f64.sqrt()).abs() < 1e-14);

        // Solver convergence onto the intersection line v = 0
        // (math_FunctionSetRoot over the function set).
        let start = [0.3f64, 0.6];
        let inf = [-1.0f64, -1.0];
        let sup = [2.0f64, 2.0];
        let mut rsnld = FunctionSetRoot::new(&func, &[1.0e-12, 1.0e-12], 100);
        rsnld.perform(&mut func, &start, &inf, &sup, false);
        assert!(rsnld.is_done(), "solver not done");
        let root = rsnld.root();
        assert!((root[0] - 0.3).abs() < 1e-9, "u drifted: {}", root[0]);
        assert!(root[1].abs() < 1e-12, "v = {}", root[1]);
        assert!(func.root().abs() < 1e-12, "Root() = {}", func.root());

        // IsTangent at the root: tgdu = N . d1v = 1/sqrt(2), tgdv = -N . d1u = 0
        // -> not tangent; Direction3d = tgdu * d1u = (1/sqrt(2), 0, 0) — the
        // 3D line direction; Direction2d = the normalized (tgdu, tgdv)
        // = (1, 0) (gxx L121-150; lxx L47-59).
        assert!(!func.is_tangent());
        let d3d = func.direction_3d();
        assert!((d3d.x - 1.0 / 2.0f64.sqrt()).abs() < 1e-14);
        assert!(d3d.y.abs() < 1e-15 && d3d.z.abs() < 1e-15);
        let d2d = func.direction_2d();
        assert!((d2d.x - 1.0).abs() < 1e-15);
        assert!(d2d.y.abs() < 1e-15);
    }

    /// OCCT anchor: Compute / Pnt / Tangency / TangencyOnSurf1 /
    /// TangencyOnSurf2 / SeekPoint (gxx L354-861) over plane quadric x
    /// BRepAdaptor plane with use_solver = false.  The point (0.5, 0) lies
    /// on the intersection line of the two planes: the tangent is
    /// Tg = N_imp x N_prm (normalized) = (1, 0, 0) and both 2D tangents
    /// satisfy Tg = d1u * Tguv.x + d1v * Tguv.y.
    #[test]
    fn imp_prm_sv_surfaces_plane_plane_anchor() {
        let i_surf = plane_quadric_z0();
        let (brep, face) = build_tilted_plane_brep();
        let p_surf = BRepAdaptorSurface::initialize_face(&brep, &face, true);

        // OCCT ctor (const TheISurface&, const ThePSurface&): implicit first.
        let sv = ImpPrmSvSurfaces::new_implicit_first(&i_surf, &p_surf);
        assert!(sv.get_use_solver()); // SetUseSolver(true) in the OCCT ctor
        sv.set_use_solver(false); // the analytic path

        // (u1, v1) = (0.5, 0) on the quadric plane, (u2, v2) = (0.5, 0) on
        // the parametric plane — both map to the point (0.5, 0, 0).
        let mut u1 = 0.5f64;
        let mut v1 = 0.0f64;
        let mut u2 = 0.5f64;
        let mut v2 = 0.0f64;
        let mut p = DVec3::ZERO;
        let mut tg = DVec3::ZERO;
        let mut tguv1 = DVec2::ZERO;
        let mut tguv2 = DVec2::ZERO;
        assert!(sv.compute(&mut u1, &mut v1, &mut u2, &mut v2, &mut p, &mut tg, &mut tguv1, &mut tguv2));

        // The intersection point (MyPnt = P = midpoint of the two
        // coincident surface points, gxx L566-622).
        assert_eq!(p, DVec3::new(0.5, 0.0, 0.0));
        // Tg = N_imp x N_prm normalized = (1, 0, 0) (gxx L693-707).
        assert_eq!(tg, DVec3::new(1.0, 0.0, 0.0));
        // Out parameters: (u2, v2) = X (no translation), (u1, v1) unchanged
        // for a plane quadric (no periodic adjustment).
        assert_eq!((u1, v1, u2, v2), (0.5, 0.0, 0.5, 0.0));
        // Tguv1 = aQuadTg = (1, 0); Tguv2 = aPrmTg = (1, 0)
        // (NonSingularProcessing: Tg = d1u * 1 + d1v * 0).
        assert_eq!(tguv1, DVec2::new(1.0, 0.0));
        assert_eq!(tguv2, DVec2::new(1.0, 0.0));

        // Pnt re-enters the MyHasBeenComputed cache (same parameters) and
        // copies MyPnt (gxx L460-466, L368).
        let mut pp = DVec3::ZERO;
        sv.pnt(0.5, 0.0, 0.5, 0.0, &mut pp);
        assert_eq!(pp, DVec3::new(0.5, 0.0, 0.0));

        // Tangency / TangencyOnSurf1 / TangencyOnSurf2 read the cached
        // MyTg / MyTguv1 / MyTguv2 (gxx L387-388, L407-408, L427-428).
        let mut t = DVec3::ZERO;
        assert!(sv.tangency(0.5, 0.0, 0.5, 0.0, &mut t));
        assert_eq!(t, DVec3::new(1.0, 0.0, 0.0));
        let mut t1 = DVec2::ZERO;
        assert!(sv.tangency_on_surf_1(0.5, 0.0, 0.5, 0.0, &mut t1));
        assert_eq!(t1, DVec2::new(1.0, 0.0));
        let mut t2 = DVec2::ZERO;
        assert!(sv.tangency_on_surf_2(0.5, 0.0, 0.5, 0.0, &mut t2));
        assert_eq!(t2, DVec2::new(1.0, 0.0));

        // SeekPoint always refines with math_FunctionSetRoot (gxx L817-819);
        // from the exact initial guess it returns the same point and
        // parameters (gxx L859).
        let mut point = PntOn2S::new();
        assert!(sv.seek_point(0.5, 0.0, 0.5, 0.0, &mut point));
        assert!((point.value() - DVec3::new(0.5, 0.0, 0.0)).length() < 1e-12);
        let (su1, sv1, su2, svv2) = point.parameters();
        assert!((su1 - 0.5).abs() < 1e-12, "NewU1 = {su1}");
        assert!(sv1.abs() < 1e-12, "NewV1 = {sv1}");
        assert!((su2 - 0.5).abs() < 1e-12, "NewU2 = {su2}");
        assert!(svv2.abs() < 1e-12, "NewV2 = {svv2}");
    }

    /// OCCT anchor: the use_solver = true Compute path — the initial guess
    /// is refined by math_FunctionSetRoot (gxx L550-558) and accepted when
    /// the refinement moved both parameters by at most 0.001
    /// (gxx L568); a guess too far from the solution is rejected
    /// (gxx L754-761).
    #[test]
    fn imp_prm_sv_surfaces_solver_path_anchor() {
        let i_surf = plane_quadric_z0();
        let (brep, face) = build_tilted_plane_brep();
        let p_surf = BRepAdaptorSurface::initialize_face(&brep, &face, true);

        // The default state is use_solver = true (SetUseSolver(true)).
        let sv = ImpPrmSvSurfaces::new_implicit_first(&i_surf, &p_surf);
        assert!(sv.get_use_solver());

        // The guess v2 = 0.0005 is within the 0.001 acceptance window of the
        // root v = 0; the solver refines it onto the line.
        let mut u1 = 0.25f64;
        let mut v1 = 0.0f64;
        let mut u2 = 0.25f64;
        let mut v2 = 0.0005f64;
        let mut p = DVec3::ZERO;
        let mut tg = DVec3::ZERO;
        let mut tguv1 = DVec2::ZERO;
        let mut tguv2 = DVec2::ZERO;
        assert!(sv.compute(&mut u1, &mut v1, &mut u2, &mut v2, &mut p, &mut tg, &mut tguv1, &mut tguv2));
        assert!(v2.abs() < 1e-9, "refined v2 = {v2}");
        assert!((u2 - 0.25).abs() < 1e-9);
        // The point is the midpoint of the parametric point and the quadric
        // point (gxx L595); both lie on z = 0 up to the refinement.
        assert!(p.z.abs() < 1e-7, "P.z = {}", p.z);
        assert!((p.x - 0.25).abs() < 1e-7);
        // Tg = N_imp x N_prm = (1, 0, 0) (plane normals are constant).
        assert!((tg - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-12);
        assert!((tguv1 - DVec2::new(1.0, 0.0)).length() < 1e-9);
        assert!((tguv2 - DVec2::new(1.0, 0.0)).length() < 1e-9);

        // A guess far from the solution: the refinement moves v2 by 0.5,
        // beyond the 0.001 window -> Compute rejects (gxx L756-761).
        let sv_far = ImpPrmSvSurfaces::new_implicit_first(&i_surf, &p_surf);
        let mut u1 = 0.25f64;
        let mut v1 = 0.0f64;
        let mut u2 = 0.25f64;
        let mut v2 = 0.5f64;
        let mut p = DVec3::ZERO;
        let mut tg = DVec3::ZERO;
        let mut tguv1 = DVec2::ZERO;
        let mut tguv2 = DVec2::ZERO;
        assert!(!sv_far.compute(&mut u1, &mut v1, &mut u2, &mut v2, &mut p, &mut tg, &mut tguv1, &mut tguv2));
    }
}
