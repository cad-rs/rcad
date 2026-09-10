// OCCT GeomLProp_SLPropsBase + GeomLProp_SurfaceUtils (TKGeomBase) — the
// generic surface local-properties engine: point, first/second derivatives,
// U/V tangents, normal (CSLib), principal curvatures/directions (first and
// second fundamental forms + math_DirectPolynomialRoots), mean and Gaussian
// curvature.
//
// GeomLProp_SLProps.hxx L39-312 (the SlPropsBase template, parameterized by
// a surface type and an access policy) + GeomLProp_SurfaceUtils.hxx L34-731
// (LProp_SurfaceUtils: EvalSurfDerivatives, SetParameters, EnsureSurfDeriv,
// FindSurfTangentOrder, ComputeSurfTangent, ComputeSurfNormal (CSLib::Normal),
// ComputeSurfCurvatures, and the IsTangentU/VDefined, TangentU/V,
// IsNormalDefined, Normal, IsCurvatureDefined, RequireCurvature, IsUmbilic,
// CurvatureDirections wrappers).
//
// The OCCT template parameters (SurfaceType, Access) map to the
// [`SLPropsSurface`] accessor trait (the ToolAccess/DirectAccess policies
// collapsed; the tool's statics become the trait methods).  The
// HLRBRep_SLProps instantiation lives in rcad-algo `hlr/brep/sl_props.rs`.

use glam::DVec3;

use crate::geom::Vec3;

use crate::base::geom_lprop::LPropStatus;
use crate::core::precision::{REAL_FIRST, REAL_LAST};
use crate::math::cs_lib::normal_from_derivatives;
use crate::math::direct_polynomial_roots::DirectPolynomialRoots;

/// OCCT RealEpsilon().
const REAL_EPSILON: f64 = 2.220446049250313e-16;

/// The `LProp_SurfaceUtils::ToolAccess<Tool>` (or DirectAccess) policy for a
/// surface: the D0/D1/D2 evaluations and the parameter bounds.
pub trait SLPropsSurface {
    /// OCCT Tool::Value(S, U, V, P) / Access::D0.
    fn eval_d0(&self, u: f64, v: f64) -> DVec3;
    /// OCCT Tool::D1(S, U, V, P, D1u, D1v) / Access::D1.
    fn eval_d1(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3);
    /// OCCT Tool::D2(S, U, V, P, D1u, D1v, D2u, D2v, Duv) / Access::D2.
    fn eval_d2(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3, Vec3, Vec3, Vec3);
    /// OCCT Tool::Bounds(S, U1, V1, U2, V2) — (u1, v1, u2, v2).
    fn bounds(&self) -> (f64, f64, f64, f64);
}

/// OCCT GeomLProp_SLPropsBase.
pub struct SlPropsBase<'a, S: ?Sized> {
    /// OCCT mySurf.
    surf: Option<&'a S>,
    /// OCCT myU / myV.
    u: f64,
    v: f64,
    /// OCCT myDerOrder.
    der_order: i32,
    /// OCCT myCN.
    cn: i32,
    /// OCCT myLinTol.
    lin_tol: f64,
    /// OCCT myPnt.
    pnt: DVec3,
    /// OCCT myD1u / myD1v.
    d1u: Vec3,
    d1v: Vec3,
    /// OCCT myD2u / myD2v / myDuv.
    d2u: Vec3,
    d2v: Vec3,
    duv: Vec3,
    /// OCCT myNormal.
    normal: Vec3,
    /// OCCT myMinCurv / myMaxCurv.
    min_curv: f64,
    max_curv: f64,
    /// OCCT myDirMinCurv / myDirMaxCurv.
    dir_min_curv: Vec3,
    dir_max_curv: Vec3,
    /// OCCT myMeanCurv / myGausCurv.
    mean_curv: f64,
    gaus_curv: f64,
    /// OCCT mySignificantFirstDerivativeOrderU / V.
    significant_first_derivative_order_u: i32,
    significant_first_derivative_order_v: i32,
    /// OCCT myUTangentStatus / myVTangentStatus / myNormalStatus /
    /// myCurvatureStatus.
    u_tangent_status: LPropStatus,
    v_tangent_status: LPropStatus,
    normal_status: LPropStatus,
    curvature_status: LPropStatus,
}

impl<'a, S: SLPropsSurface + ?Sized> SlPropsBase<'a, S> {
    /// OCCT GeomLProp_SLPropsBase(N, Resolution) (hxx L94-108) — the surface
    /// is set later with SetSurface.
    pub fn new(n: i32, resolution: f64) -> Self {
        assert!((0..=2).contains(&n), "Standard_OutOfRange: SLProps(N)");
        SlPropsBase {
            surf: None,
            u: REAL_LAST, // OCCT GeomLProp_SLPropsBase.hxx: RealLast()
            v: REAL_LAST,
            der_order: n,
            cn: 0,
            lin_tol: resolution,
            pnt: DVec3::ZERO,
            d1u: Vec3::ZERO,
            d1v: Vec3::ZERO,
            d2u: Vec3::ZERO,
            d2v: Vec3::ZERO,
            duv: Vec3::ZERO,
            normal: Vec3::ZERO,
            min_curv: 0.0,
            max_curv: 0.0,
            dir_min_curv: Vec3::ZERO,
            dir_max_curv: Vec3::ZERO,
            mean_curv: 0.0,
            gaus_curv: 0.0,
            significant_first_derivative_order_u: 0,
            significant_first_derivative_order_v: 0,
            u_tangent_status: LPropStatus::Undecided,
            v_tangent_status: LPropStatus::Undecided,
            normal_status: LPropStatus::Undecided,
            curvature_status: LPropStatus::Undecided,
        }
    }

    fn s(&self) -> &'a S {
        self.surf.expect("SLProps: no surface set")
    }

    /// OCCT SetSurface(S) (hxx L112-116).
    pub fn set_surface(&mut self, s: &'a S) {
        self.surf = Some(s);
        self.cn = 4;
    }

    /// OCCT SetParameters(U, V) via LProp_SurfaceUtils::SetParameters
    /// (SurfaceUtils.hxx L485-511 + EvalSurfDerivatives L186-210).
    pub fn set_parameters(&mut self, u: f64, v: f64) {
        self.u = u;
        self.v = v;
        match self.der_order {
            0 => self.pnt = self.s().eval_d0(u, v),
            1 => {
                let (p, d1u, d1v) = self.s().eval_d1(u, v);
                self.pnt = p;
                self.d1u = d1u;
                self.d1v = d1v;
            }
            2 => {
                let (p, d1u, d1v, d2u, d2v, duv) = self.s().eval_d2(u, v);
                self.pnt = p;
                self.d1u = d1u;
                self.d1v = d1v;
                self.d2u = d2u;
                self.d2v = d2v;
                self.duv = duv;
            }
            _ => {}
        }
        self.u_tangent_status = LPropStatus::Undecided;
        self.v_tangent_status = LPropStatus::Undecided;
        self.normal_status = LPropStatus::Undecided;
        self.curvature_status = LPropStatus::Undecided;
    }

    /// OCCT Value() (hxx L141).
    pub fn value(&self) -> DVec3 {
        self.pnt
    }

    /// OCCT D1U() via EnsureSurfDeriv (SurfaceUtils.hxx L514-535).
    pub fn d1u(&mut self) -> Vec3 {
        if self.der_order < 1 {
            self.der_order = 1;
            let (p, d1u, d1v) = self.s().eval_d1(self.u, self.v);
            self.pnt = p;
            self.d1u = d1u;
            self.d1v = d1v;
        }
        self.d1u
    }

    /// OCCT D1V().
    pub fn d1v(&mut self) -> Vec3 {
        if self.der_order < 1 {
            self.der_order = 1;
            let (p, d1u, d1v) = self.s().eval_d1(self.u, self.v);
            self.pnt = p;
            self.d1u = d1u;
            self.d1v = d1v;
        }
        self.d1v
    }

    /// OCCT D2U().
    pub fn d2u(&mut self) -> Vec3 {
        if self.der_order < 2 {
            self.der_order = 2;
            let (p, d1u, d1v, d2u, d2v, duv) = self.s().eval_d2(self.u, self.v);
            self.pnt = p;
            self.d1u = d1u;
            self.d1v = d1v;
            self.d2u = d2u;
            self.d2v = d2v;
            self.duv = duv;
        }
        self.d2u
    }

    /// OCCT D2V().
    pub fn d2v(&mut self) -> Vec3 {
        if self.der_order < 2 {
            self.der_order = 2;
            let (p, d1u, d1v, d2u, d2v, duv) = self.s().eval_d2(self.u, self.v);
            self.pnt = p;
            self.d1u = d1u;
            self.d1v = d1v;
            self.d2u = d2u;
            self.d2v = d2v;
            self.duv = duv;
        }
        self.d2v
    }

    /// OCCT DUV().
    pub fn duv(&mut self) -> Vec3 {
        if self.der_order < 2 {
            self.der_order = 2;
            let (p, d1u, d1v, d2u, d2v, duv) = self.s().eval_d2(self.u, self.v);
            self.pnt = p;
            self.d1u = d1u;
            self.d1v = d1v;
            self.d2u = d2u;
            self.d2v = d2v;
            self.duv = duv;
        }
        self.duv
    }

    /// OCCT IsTangentUDefined via LProp_SurfaceUtils::IsTangentUDefined +
    /// FindSurfTangentOrder (SurfaceUtils.hxx L538-555 + L221-251).
    pub fn is_tangent_u_defined(&mut self) -> bool {
        if self.u_tangent_status == LPropStatus::Undefined {
            return false;
        }
        if self.u_tangent_status as i32 >= LPropStatus::Defined as i32 {
            return true;
        }
        let d1u = self.d1u();
        let d2u = self.d2u();
        find_surf_tangent_order(
            d1u,
            d2u,
            self.cn,
            self.lin_tol * self.lin_tol,
            &mut self.significant_first_derivative_order_u,
            &mut self.u_tangent_status,
        )
    }

    /// OCCT IsTangentVDefined (SurfaceUtils.hxx L557-575).
    pub fn is_tangent_v_defined(&mut self) -> bool {
        if self.v_tangent_status == LPropStatus::Undefined {
            return false;
        }
        if self.v_tangent_status as i32 >= LPropStatus::Defined as i32 {
            return true;
        }
        let d1v = self.d1v();
        let d2v = self.d2v();
        find_surf_tangent_order(
            d1v,
            d2v,
            self.cn,
            self.lin_tol * self.lin_tol,
            &mut self.significant_first_derivative_order_v,
            &mut self.v_tangent_status,
        )
    }

    /// OCCT TangentU(D) via TangentU + ComputeSurfTangent (SurfaceUtils.hxx
    /// L577-591 + L264-342).  Returns None for the LProp_NotDefined throw.
    pub fn tangent_u(&mut self) -> Option<Vec3> {
        if !self.is_tangent_u_defined() {
            return None; // LProp_NotDefined
        }
        if self.significant_first_derivative_order_u == 1 {
            return Some(self.d1u.normalize_or_zero()); // gp_Dir(theFirstDeriv)
        }
        let (an_uinfimum, an_vinfimum, an_usupremum, _an_vsupremum) = self.s().bounds();
        // OCCT LProp_SurfaceUtils.hxx, the U arm of ComputeSurfTangent:
        // if ((anUSupremum >= RealLast()) || (anUinfimum <= RealFirst()))
        let a_du = if an_usupremum >= REAL_LAST || an_uinfimum <= REAL_FIRST {
            0.0
        } else {
            an_usupremum - an_uinfimum
        };
        let a_delta = (a_du * 1.0e-3).max(1.0e-7);
        let mut a_v = self.d2u;
        let an_other_u = if self.u - an_uinfimum < a_delta {
            self.u + a_delta
        } else {
            self.u - a_delta
        };
        let p1 = self.s().eval_d0(self.u.min(an_other_u), self.v);
        let p2 = self.s().eval_d0(self.u.max(an_other_u), self.v);
        let a_chord = p2 - p1;
        if a_v.dot(a_chord) < 0.0 {
            a_v = -a_v;
        }
        Some(a_v.normalize_or_zero())
    }

    /// OCCT TangentV(D) (SurfaceUtils.hxx L593-607 + the V arm of
    /// ComputeSurfTangent).
    pub fn tangent_v(&mut self) -> Option<Vec3> {
        if !self.is_tangent_v_defined() {
            return None; // LProp_NotDefined
        }
        if self.significant_first_derivative_order_v == 1 {
            return Some(self.d1v.normalize_or_zero());
        }
        let (an_uinfimum, an_vinfimum, _an_usupremum, an_vsupremum) = self.s().bounds();
        // OCCT LProp_SurfaceUtils.hxx, the V arm of ComputeSurfTangent:
        // if ((anVSupremum >= RealLast()) || (anVinfimum <= RealFirst()))
        let a_dv = if an_vsupremum >= REAL_LAST || an_vinfimum <= REAL_FIRST {
            0.0
        } else {
            an_vsupremum - an_vinfimum
        };
        let a_delta = (a_dv * 1.0e-3).max(1.0e-7);
        let mut a_v = self.d2v;
        let an_other_v = if self.v - an_vinfimum < a_delta {
            self.v + a_delta
        } else {
            self.v - a_delta
        };
        let p1 = self.s().eval_d0(self.u, self.v.min(an_other_v));
        let p2 = self.s().eval_d0(self.u, self.v.max(an_other_v));
        let a_chord = p2 - p1;
        if a_v.dot(a_chord) < 0.0 {
            a_v = -a_v;
        }
        Some(a_v.normalize_or_zero())
    }

    /// OCCT IsNormalDefined via LProp_SurfaceUtils::IsNormalDefined +
    /// ComputeSurfNormal (SurfaceUtils.hxx L609-627 + L350-358; the CSLib
    /// normal with the linear tolerance passed as the sine tolerance,
    /// verbatim).
    pub fn is_normal_defined(&mut self) -> bool {
        if self.normal_status == LPropStatus::Undefined {
            return false;
        }
        if self.normal_status as i32 >= LPropStatus::Defined as i32 {
            return true;
        }
        let (normal, status) = normal_from_derivatives(self.d1u, self.d1v, self.lin_tol);
        if status == crate::math::cs_lib::DerivativeStatus::Done {
            self.normal = normal.expect("CSLib::Normal Done without a direction");
            self.normal_status = LPropStatus::Computed;
            true
        } else {
            self.normal_status = LPropStatus::Undefined;
            false
        }
    }

    /// OCCT Normal() (SurfaceUtils.hxx L629-636).  Returns None for the
    /// LProp_NotDefined throw.
    pub fn normal(&mut self) -> Option<Vec3> {
        if !self.is_normal_defined() {
            return None;
        }
        Some(self.normal)
    }

    /// OCCT IsCurvatureDefined (SurfaceUtils.hxx L638-697).
    pub fn is_curvature_defined(&mut self) -> bool {
        if self.curvature_status == LPropStatus::Undefined {
            return false;
        }
        if self.curvature_status as i32 >= LPropStatus::Defined as i32 {
            return true;
        }
        if self.cn < 2 {
            self.curvature_status = LPropStatus::Undefined;
            return false;
        }
        if !self.is_normal_defined() {
            self.curvature_status = LPropStatus::Undefined;
            return false;
        }
        if !self.is_tangent_u_defined() || !self.is_tangent_v_defined() {
            self.curvature_status = LPropStatus::Undefined;
            return false;
        }
        if self.der_order < 2 {
            self.d2u();
        }
        if compute_surf_curvatures(
            self.d1u,
            self.d1v,
            self.d2u,
            self.d2v,
            self.duv,
            self.normal,
            &mut self.min_curv,
            &mut self.max_curv,
            &mut self.dir_min_curv,
            &mut self.dir_max_curv,
            &mut self.mean_curv,
            &mut self.gaus_curv,
        ) {
            self.curvature_status = LPropStatus::Computed;
            return true;
        }
        self.curvature_status = LPropStatus::Undefined;
        false
    }

    /// OCCT IsUmbilic (SurfaceUtils.hxx L708-715).
    pub fn is_umbilic(&mut self) -> bool {
        if !self.is_curvature_defined() {
            panic!("LProp_NotDefined");
        }
        // std::abs(Epsilon(theMaxCurv)) — Epsilon(A) = A * DBL_EPSILON
        // (the gap to the next representable double at A).
        (self.max_curv - self.min_curv).abs() < (self.max_curv * f64::EPSILON).abs()
    }

    /// OCCT MaxCurvature via RequireCurvature (SurfaceUtils.hxx L699-706).
    pub fn max_curvature(&mut self) -> f64 {
        if !self.is_curvature_defined() {
            panic!("LProp_NotDefined");
        }
        self.max_curv
    }

    /// OCCT MinCurvature.
    pub fn min_curvature(&mut self) -> f64 {
        if !self.is_curvature_defined() {
            panic!("LProp_NotDefined");
        }
        self.min_curv
    }

    /// OCCT CurvatureDirections (SurfaceUtils.hxx L717-729) —
    /// (max direction, min direction).
    pub fn curvature_directions(&mut self) -> (Vec3, Vec3) {
        if !self.is_curvature_defined() {
            panic!("LProp_NotDefined");
        }
        (self.dir_max_curv, self.dir_min_curv)
    }

    /// OCCT MeanCurvature.
    pub fn mean_curvature(&mut self) -> f64 {
        if !self.is_curvature_defined() {
            panic!("LProp_NotDefined");
        }
        self.mean_curv
    }

    /// OCCT GaussianCurvature.
    pub fn gaussian_curvature(&mut self) -> f64 {
        if !self.is_curvature_defined() {
            panic!("LProp_NotDefined");
        }
        self.gaus_curv
    }
}

/// OCCT FindSurfTangentOrder (SurfaceUtils.hxx L221-251).
#[allow(clippy::needless_while_loop)]
fn find_surf_tangent_order(
    the_d1: Vec3,
    the_d2: Vec3,
    the_cn: i32,
    the_tol_sq: f64,
    the_order: &mut i32,
    the_status: &mut LPropStatus,
) -> bool {
    let a_derivs = [the_d1, the_d2];
    *the_order = 0;

    while *the_order < 3 {
        *the_order += 1;
        if the_cn >= *the_order {
            if *the_order <= 2 && a_derivs[(*the_order - 1) as usize].length_squared() > the_tol_sq
            {
                *the_status = LPropStatus::Defined;
                return true;
            }
        } else {
            *the_status = LPropStatus::Undefined;
            return false;
        }
    }

    *the_status = LPropStatus::Undefined;
    false
}

/// OCCT ComputeSurfCurvatures (SurfaceUtils.hxx L376-480) — the principal
/// curvatures/directions from the first and second fundamental forms; the
//  shape-operator root equation is solved with math_DirectPolynomialRoots.
#[allow(clippy::too_many_arguments)]
fn compute_surf_curvatures(
    the_d1u: Vec3,
    the_d1v: Vec3,
    the_d2u: Vec3,
    the_d2v: Vec3,
    the_duv: Vec3,
    the_normal: Vec3,
    the_min_curv: &mut f64,
    the_max_curv: &mut f64,
    the_dir_min: &mut Vec3,
    the_dir_max: &mut Vec3,
    the_mean_curv: &mut f64,
    the_gaus_curv: &mut f64,
) -> bool {
    let a_norm = the_normal;

    let an_e = the_d1u.length_squared();
    let an_f = the_d1u.dot(the_d1v);
    let a_g = the_d1v.length_squared();

    let a_l = a_norm.dot(the_d2u);
    let a_m = a_norm.dot(the_duv);
    let a_n = a_norm.dot(the_d2v);

    let an_a0 = an_e * a_m - an_f * a_l;
    let a_b0 = an_e * a_n - a_g * a_l;
    let a_c0 = an_f * a_n - a_g * a_m;

    let a_max_abc = an_a0.abs().max(a_b0.abs()).max(a_c0.abs());
    if a_max_abc < REAL_EPSILON {
        // Umbilic point
        if a_g < REAL_EPSILON {
            return false;
        }
        *the_min_curv = a_n / a_g;
        *the_max_curv = *the_min_curv;
        *the_dir_min = the_d1u.normalize_or_zero(); // gp_Dir(theD1u)
        *the_dir_max = the_d1u.cross(a_norm).normalize_or_zero();
        *the_mean_curv = *the_min_curv;
        *the_gaus_curv = *the_min_curv * *the_min_curv;
        return true;
    }

    let an_a = an_a0 / a_max_abc;
    let a_b = a_b0 / a_max_abc;
    let a_c = a_c0 / a_max_abc;

    let a_curv1: f64;
    let a_curv2: f64;
    let a_vect_curv1: Vec3;
    let a_vect_curv2: Vec3;

    if an_a.abs() > REAL_EPSILON {
        // math_DirectPolynomialRoots aRoot(anA, aB, aC): A x^2 + B x + C.
        let a_root = DirectPolynomialRoots::new_quadratic(an_a, a_b, a_c);
        if a_root.nb_solutions() != 2 {
            return false;
        }

        let a_root1 = a_root.value(1);
        let a_root2 = a_root.value(2);
        a_curv1 = ((a_l * a_root1 + 2.0 * a_m) * a_root1 + a_n)
            / ((an_e * a_root1 + 2.0 * an_f) * a_root1 + a_g);
        a_curv2 = ((a_l * a_root2 + 2.0 * a_m) * a_root2 + a_n)
            / ((an_e * a_root2 + 2.0 * an_f) * a_root2 + a_g);
        a_vect_curv1 = a_root1 * the_d1u + the_d1v;
        a_vect_curv2 = a_root2 * the_d1u + the_d1v;
    } else if a_c.abs() > REAL_EPSILON {
        let a_root = DirectPolynomialRoots::new_quadratic(a_c, a_b, an_a);
        if a_root.nb_solutions() != 2 {
            return false;
        }

        let a_root1 = a_root.value(1);
        let a_root2 = a_root.value(2);
        a_curv1 = ((a_n * a_root1 + 2.0 * a_m) * a_root1 + a_l)
            / ((a_g * a_root1 + 2.0 * an_f) * a_root1 + an_e);
        a_curv2 = ((a_n * a_root2 + 2.0 * a_m) * a_root2 + a_l)
            / ((a_g * a_root2 + 2.0 * an_f) * a_root2 + an_e);
        a_vect_curv1 = the_d1u + a_root1 * the_d1v;
        a_vect_curv2 = the_d1u + a_root2 * the_d1v;
    } else {
        a_curv1 = a_l / an_e;
        a_curv2 = a_n / a_g;
        a_vect_curv1 = the_d1u;
        a_vect_curv2 = the_d1v;
    }

    if a_curv1 < a_curv2 {
        *the_min_curv = a_curv1;
        *the_max_curv = a_curv2;
        *the_dir_min = a_vect_curv1.normalize_or_zero();
        *the_dir_max = a_vect_curv2.normalize_or_zero();
    } else {
        *the_min_curv = a_curv2;
        *the_max_curv = a_curv1;
        *the_dir_min = a_vect_curv2.normalize_or_zero();
        *the_dir_max = a_vect_curv1.normalize_or_zero();
    }

    let an_eg_ff = (an_e * a_g) - (an_f * an_f);
    *the_mean_curv = ((a_n * an_e) - (2.0 * a_m * an_f) + (a_l * a_g)) / (2.0 * an_eg_ff);
    *the_gaus_curv = ((a_l * a_n) - (a_m * a_m)) / an_eg_ff;
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The unit cylinder along Z: P(u,v) = (cos u, sin u, v).
    struct CylSurf;

    impl SLPropsSurface for CylSurf {
        fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
            DVec3::new(u.cos(), u.sin(), v)
        }
        fn eval_d1(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3) {
            (
                DVec3::new(u.cos(), u.sin(), v),
                Vec3::new(-u.sin(), u.cos(), 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            )
        }
        fn eval_d2(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3, Vec3, Vec3, Vec3) {
            (
                DVec3::new(u.cos(), u.sin(), v),
                Vec3::new(-u.sin(), u.cos(), 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(-u.cos(), -u.sin(), 0.0),
                Vec3::ZERO,
                Vec3::ZERO,
            )
        }
        fn bounds(&self) -> (f64, f64, f64, f64) {
            (0.0, 0.0, std::f64::consts::TAU, 2.0)
        }
    }

    /// OCCT anchor: on the unit cylinder the principal curvatures are
    /// 0 (along the axis) and 1 (around the circle), mean = 1/2, Gaussian
    /// = 0; the normal is the radial direction.
    #[test]
    fn slprops_cylinder_curvatures() {
        let mut s = SlPropsBase::<dyn SLPropsSurface>::new(2, 1e-12);
        s.set_surface(&CylSurf);
        s.set_parameters(0.3, 0.7);

        // Normal = radial.
        let n = s.normal().expect("normal defined");
        assert!((n.x - 0.3f64.cos()).abs() < 1e-12);
        assert!((n.y - 0.3f64.sin()).abs() < 1e-12);
        assert!(n.z.abs() < 1e-12);

        // Principal curvatures with the outward (radial) normal:
        // L = n . D2u = -1, N = n . D2v = 0 -> curvatures 0 (axial, max)
        // and -1 (circular, min); mean = -1/2, Gaussian = 0.
        assert!(s.max_curvature().abs() < 1e-9, "max={}", s.max_curv);
        assert!((s.min_curvature() + 1.0).abs() < 1e-9);
        assert!((s.mean_curvature() + 0.5).abs() < 1e-9);
        assert!(s.gaussian_curvature().abs() < 1e-9);

        // The max-curvature direction is the axial direction (D1v).
        let (dmax, dmin) = s.curvature_directions();
        assert!((dmax.z - 1.0).abs() < 1e-9);
        assert!((dmin.x + 0.3f64.sin()).abs() < 1e-9 && (dmin.y - 0.3f64.cos()).abs() < 1e-9);
    }

    /// OCCT anchor: the sphere is umbilic — both principal curvatures equal
    /// 1/R and IsUmbilic holds.
    #[test]
    fn slprops_sphere_umbilic() {
        /// The unit sphere: P(u,v) = (cos u sin v, sin u sin v, cos v).
        struct SphereSurf;
        impl SLPropsSurface for SphereSurf {
            fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
                DVec3::new(u.cos() * v.sin(), u.sin() * v.sin(), v.cos())
            }
            fn eval_d1(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3) {
                let sv = v.sin();
                let cv = v.cos();
                (
                    DVec3::new(u.cos() * sv, u.sin() * sv, cv),
                    Vec3::new(-u.sin() * sv, u.cos() * sv, 0.0),
                    Vec3::new(u.cos() * cv, u.sin() * cv, -sv),
                )
            }
            fn eval_d2(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3, Vec3, Vec3, Vec3) {
                let sv = v.sin();
                let cv = v.cos();
                (
                    DVec3::new(u.cos() * sv, u.sin() * sv, cv),
                    Vec3::new(-u.sin() * sv, u.cos() * sv, 0.0),
                    Vec3::new(u.cos() * cv, u.sin() * cv, -sv),
                    // D2u = (-sin u sin v, cos u sin v, 0).
                    Vec3::new(-u.cos() * sv, -u.sin() * sv, 0.0),
                    // D2v = (-sin u sin v, cos u sin v, -cos v) = -P.
                    Vec3::new(-u.cos() * sv, -u.sin() * sv, -cv),
                    // Duv = (-sin u cos v, cos u cos v, 0).
                    Vec3::new(-u.sin() * cv, u.cos() * cv, 0.0),
                )
            }
            fn bounds(&self) -> (f64, f64, f64, f64) {
                (0.0, 0.0, std::f64::consts::TAU, std::f64::consts::PI)
            }
        }

        let mut s = SlPropsBase::<dyn SLPropsSurface>::new(2, 1e-12);
        s.set_surface(&SphereSurf);
        s.set_parameters(0.9, 1.1);

        assert!(s.is_normal_defined());
        assert!(s.is_umbilic(), "the sphere point is umbilic");
        assert!((s.max_curvature() - 1.0).abs() < 1e-9, "max={}", s.max_curv);
        assert!((s.min_curvature() - 1.0).abs() < 1e-9);
        assert!((s.gaussian_curvature() - 1.0).abs() < 1e-9);
        assert!((s.mean_curvature() - 1.0).abs() < 1e-9);
    }
}
