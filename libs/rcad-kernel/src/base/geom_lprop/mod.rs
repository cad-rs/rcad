//! Curve and surface local properties (GeomLProp).
//!
//! OCCT TKGeomBase GeomLProp package: GeomLProp_CLProps (curve local properties),
//! GeomLProp_SLProps (surface local properties).

#![allow(clippy::manual_clamp)]

use glam::{DVec2, DVec3};

use crate::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};

/// Status of a local property computation.
///
/// OCCT: `LProp_Status`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LPropStatus {
    Defined,
    Undecided,
    Zero,
}

const TOL_LIN: f64 = 1e-12;

// ============================================================================
// GeomLProp_CLProps — Curve Local Properties
// ============================================================================

/// Computes local properties of a 3D curve at a parameter value:
/// point, derivatives (D1, D2, D3), tangent, curvature, normal, centre of curvature.
///
/// OCCT: `GeomLProp_CLProps`.
pub struct CLProps<'a> {
    curve: &'a Curve3,
    u: f64,
    der_order: i32,
    /// Continuity index (4 = unknown)
    cn: i32,
    lin_tol: f64,
    pnt: DVec3,
    deriv: [DVec3; 3],
    curvature: f64,
    tangent_status: LPropStatus,
    significant_first_derivative_order: i32,
}

impl<'a> CLProps<'a> {
    /// Constructor with curve, max derivative order, and linear tolerance.
    ///
    /// OCCT: `GeomLProp_CLProps(Curve, N, Resolution)`.
    pub fn new(curve: &'a Curve3, n: i32, resolution: f64) -> Self {
        assert!(n >= 0 && n <= 3, "CLProps: N must be 0, 1, 2, or 3");
        CLProps {
            curve,
            u: f64::NAN,
            der_order: n,
            cn: 4,
            lin_tol: if resolution <= 0.0 { TOL_LIN } else { resolution },
            pnt: DVec3::ZERO,
            deriv: [DVec3::ZERO, DVec3::ZERO, DVec3::ZERO],
            curvature: 0.0,
            tangent_status: LPropStatus::Undecided,
            significant_first_derivative_order: 0,
        }
    }

    /// Constructor with curve, parameter, max derivative order, and tolerance.
    ///
    /// OCCT: `GeomLProp_CLProps(Curve, U, N, Resolution)`.
    pub fn with_param(curve: &'a Curve3, u: f64, n: i32, resolution: f64) -> Self {
        let mut props = CLProps::new(curve, n, resolution);
        props.set_parameter(u);
        props
    }

    /// Set the parameter value and compute derivatives.
    ///
    /// OCCT: `SetParameter(U)`.
    pub fn set_parameter(&mut self, u: f64) {
        self.u = u;
        self.pnt = self.curve.point_at(u);
        if self.der_order >= 1 {
            self.deriv[0] = self.curve.derivative_at(u);
        }
        if self.der_order >= 2 {
            self.deriv[1] = self.curve.derivative2_at(u);
        }
        if self.der_order >= 3 {
            self.deriv[2] = self.curve.derivative3_at(u);
        }
        self.tangent_status = LPropStatus::Undecided;
    }

    /// Returns the point on the curve.
    ///
    /// OCCT: `Value()`.
    pub fn value(&self) -> DVec3 {
        self.pnt
    }

    /// Returns the first derivative (computed if not yet).
    ///
    /// OCCT: `D1()`.
    pub fn d1(&mut self) -> DVec3 {
        if self.der_order >= 1 && self.deriv[0] == DVec3::ZERO {
            self.deriv[0] = self.curve.derivative_at(self.u);
        }
        self.deriv[0]
    }

    /// Returns the second derivative.
    ///
    /// OCCT: `D2()`.
    pub fn d2(&mut self) -> DVec3 {
        if self.der_order >= 2 && self.deriv[1] == DVec3::ZERO {
            self.deriv[1] = self.curve.derivative2_at(self.u);
        }
        self.deriv[1]
    }

    /// Returns the third derivative.
    ///
    /// OCCT: `D3()`.
    pub fn d3(&mut self) -> DVec3 {
        if self.der_order >= 3 && self.deriv[2] == DVec3::ZERO {
            self.deriv[2] = self.curve.derivative3_at(self.u);
        }
        self.deriv[2]
    }

    /// Returns true if the tangent is defined.
    ///
    /// OCCT: `IsTangentDefined()`.
    pub fn is_tangent_defined(&mut self) -> bool {
        if self.tangent_status != LPropStatus::Undecided {
            return self.tangent_status == LPropStatus::Defined;
        }
        // Check if first derivative is non-zero
        let d1 = self.d1();
        if d1.length_squared() > self.lin_tol * self.lin_tol {
            self.tangent_status = LPropStatus::Defined;
            self.significant_first_derivative_order = 1;
            return true;
        }
        if self.der_order < 2 {
            self.tangent_status = LPropStatus::Zero;
            return false;
        }
        let d2 = self.d2();
        if d2.length_squared() > self.lin_tol * self.lin_tol {
            self.tangent_status = LPropStatus::Defined;
            self.significant_first_derivative_order = 2;
            return true;
        }
        if self.der_order < 3 {
            self.tangent_status = LPropStatus::Zero;
            return false;
        }
        let d3 = self.d3();
        if d3.length_squared() > self.lin_tol * self.lin_tol {
            self.tangent_status = LPropStatus::Defined;
            self.significant_first_derivative_order = 3;
            return true;
        }
        self.tangent_status = LPropStatus::Zero;
        false
    }

    /// Returns the tangent direction.
    ///
    /// OCCT: `Tangent(D)`.
    pub fn tangent(&mut self) -> Option<DVec3> {
        if !self.is_tangent_defined() {
            return None;
        }
        let idx = (self.significant_first_derivative_order - 1) as usize;
        let d = self.deriv[idx];
        Some(d.normalize_or_zero())
    }

    /// Returns the curvature.
    ///
    /// OCCT: `Curvature()`.
    pub fn curvature(&mut self) -> f64 {
        if !self.is_tangent_defined() {
            return 0.0;
        }
        let d1 = self.d1();
        let speed = d1.length();
        if speed < self.lin_tol {
            return 0.0;
        }
        let d2 = self.d2();
        self.curvature = d1.cross(d2).length() / (speed * speed * speed);
        self.curvature
    }

    /// Returns the normal direction (principal normal).
    ///
    /// OCCT: `Normal(N)`.
    pub fn normal(&mut self) -> Option<DVec3> {
        let curv = self.curvature();
        if curv < self.lin_tol {
            return None;
        }
        let d1 = self.d1();
        let d2 = self.d2();
        // N = (C' × C'') × C' / |(C' × C'') × C'|
        let cross_d1d2 = d1.cross(d2);
        let n = cross_d1d2.cross(d1);
        let len = n.length();
        if len < self.lin_tol {
            return None;
        }
        Some(n / len)
    }

    /// Returns the centre of curvature.
    ///
    /// OCCT: `CentreOfCurvature(P)`.
    pub fn centre_of_curvature(&mut self) -> Option<DVec3> {
        let curv = self.curvature();
        if curv < self.lin_tol {
            return None;
        }
        let n = self.normal()?;
        Some(self.pnt + n / curv)
    }
}

// ============================================================================
// GeomLProp_CLProps2d — Curve Local Properties (2D)
// ============================================================================

/// Computes local properties of a 2D curve at a parameter value: point,
/// first/second derivatives, tangent, curvature, normal.
///
/// OCCT: `GeomLProp_CLProps2d` (GeomLProp_CLProps.hxx — the 2D instantiation
/// of `GeomLProp_CLPropsBase`; formulas in LProp_CurveUtils.hxx):
/// - curvature = |D1 × D2| / |D1|³ (2D cross magnitude)
/// - normal = D2·|D1|² − D1·(D1·D2), normalized
pub struct ClProps2d<'a> {
    curve: &'a Curve2d,
    u: f64,
    der_order: i32,
    /// Continuity index (4 = unknown).
    cn: i32,
    lin_tol: f64,
    pnt: DVec2,
    deriv: [DVec2; 3],
    curvature: f64,
    tangent_status: LPropStatus,
    significant_first_derivative_order: i32,
}

impl<'a> ClProps2d<'a> {
    /// OCCT: `GeomLProp_CLProps2d(C, N, Resolution)`.
    pub fn new(curve: &'a Curve2d, n: i32, resolution: f64) -> Self {
        assert!(n >= 0 && n <= 3, "ClProps2d: N must be 0, 1, 2, or 3");
        ClProps2d {
            curve,
            u: f64::NAN,
            der_order: n,
            cn: 4,
            lin_tol: if resolution <= 0.0 { TOL_LIN } else { resolution },
            pnt: DVec2::ZERO,
            deriv: [DVec2::ZERO, DVec2::ZERO, DVec2::ZERO],
            curvature: 0.0,
            tangent_status: LPropStatus::Undecided,
            significant_first_derivative_order: 0,
        }
    }

    /// OCCT: `GeomLProp_CLProps2d(C, U, N, Resolution)`.
    pub fn with_param(curve: &'a Curve2d, u: f64, n: i32, resolution: f64) -> Self {
        let mut props = ClProps2d::new(curve, n, resolution);
        props.set_parameter(u);
        props
    }

    /// OCCT: `SetParameter(U)`.
    pub fn set_parameter(&mut self, u: f64) {
        self.u = u;
        self.pnt = self.curve.point_at(u);
        if self.der_order >= 1 {
            self.deriv[0] = self.curve.derivative_at(u);
        }
        if self.der_order >= 2 {
            self.deriv[1] = self.curve.derivative2_at(u);
        }
        if self.der_order >= 3 {
            // 2D third derivative not exposed by Curve2dEval; finite difference.
            let h = 1e-4;
            self.deriv[2] = (self.curve.derivative2_at(u + h) - self.curve.derivative2_at(u - h))
                / (2.0 * h);
        }
        self.tangent_status = LPropStatus::Undecided;
    }

    /// OCCT: `Value()`.
    pub fn value(&self) -> DVec2 {
        self.pnt
    }

    /// OCCT: `D1()`.
    pub fn d1(&mut self) -> DVec2 {
        if self.der_order >= 1 && self.deriv[0] == DVec2::ZERO {
            self.deriv[0] = self.curve.derivative_at(self.u);
        }
        self.deriv[0]
    }

    /// OCCT: `D2()`.
    pub fn d2(&mut self) -> DVec2 {
        if self.der_order >= 2 && self.deriv[1] == DVec2::ZERO {
            self.deriv[1] = self.curve.derivative2_at(self.u);
        }
        self.deriv[1]
    }

    /// OCCT: `IsTangentDefined()`.
    pub fn is_tangent_defined(&mut self) -> bool {
        if self.tangent_status != LPropStatus::Undecided {
            return self.tangent_status == LPropStatus::Defined;
        }
        let d1 = self.d1();
        if d1.length_squared() > self.lin_tol * self.lin_tol {
            self.tangent_status = LPropStatus::Defined;
            self.significant_first_derivative_order = 1;
            return true;
        }
        if self.der_order < 2 {
            self.tangent_status = LPropStatus::Zero;
            return false;
        }
        let d2 = self.d2();
        if d2.length_squared() > self.lin_tol * self.lin_tol {
            self.tangent_status = LPropStatus::Defined;
            self.significant_first_derivative_order = 2;
            return true;
        }
        if self.der_order < 3 {
            self.tangent_status = LPropStatus::Zero;
            return false;
        }
        let d3 = self.d3();
        if d3.length_squared() > self.lin_tol * self.lin_tol {
            self.tangent_status = LPropStatus::Defined;
            self.significant_first_derivative_order = 3;
            return true;
        }
        self.tangent_status = LPropStatus::Zero;
        false
    }

    /// OCCT: `D3()`.
    pub fn d3(&mut self) -> DVec2 {
        if self.der_order >= 3 && self.deriv[2] == DVec2::ZERO {
            let h = 1e-4;
            self.deriv[2] = (self.curve.derivative2_at(self.u + h)
                - self.curve.derivative2_at(self.u - h))
                / (2.0 * h);
        }
        self.deriv[2]
    }

    /// OCCT: `Tangent(D)`.
    pub fn tangent(&mut self) -> Option<DVec2> {
        if !self.is_tangent_defined() {
            return None;
        }
        let idx = (self.significant_first_derivative_order - 1) as usize;
        let mut d = self.deriv[idx];
        if self.significant_first_derivative_order > 1 {
            // OCCT LProp_CurveUtils::ComputeTangent — sign correction using
            // the chord near the point (the tangent of a curve whose first
            // significant derivative is of higher order).
            let dom = self.curve.default_domain();
            let (inf, sup) = (dom[0], dom[1]);
            let a_du = if sup.is_infinite() || inf.is_infinite() {
                0.0
            } else {
                sup - inf
            };
            let delta = (a_du * 1.0e-3).max(1.0e-7);
            let other_u = if self.u - inf < delta { self.u + delta } else { self.u - delta };
            let p1 = self.curve.point_at(self.u.min(other_u));
            let p2 = self.curve.point_at(self.u.max(other_u));
            let chord = p2 - p1;
            if d.dot(chord) < 0.0 {
                d = -d;
            }
        }
        Some(d.normalize_or_zero())
    }

    /// OCCT: `Curvature()`. Returns |D1 × D2| / |D1|³.
    pub fn curvature(&mut self) -> f64 {
        if !self.is_tangent_defined() {
            return 0.0;
        }
        let d1 = self.d1();
        let d2 = self.d2();
        let a_dd1 = d1.length_squared();
        let a_dd2 = d2.length_squared();
        if a_dd2 <= self.lin_tol * self.lin_tol {
            return 0.0;
        }
        let a_n = (d1.x * d2.y - d1.y * d2.x).powi(2); // CrossSquareMagnitude
        let a_t = a_n / a_dd1 / a_dd2;
        if a_t <= self.lin_tol * self.lin_tol {
            return 0.0;
        }
        self.curvature = a_n.sqrt() / a_dd1 / a_dd1.sqrt();
        self.curvature
    }

    /// OCCT: `Normal(N)` — D2·|D1|² − D1·(D1·D2), normalized.
    pub fn normal(&mut self) -> Option<DVec2> {
        let curv = self.curvature();
        if curv.abs() <= self.lin_tol {
            return None;
        }
        let d1 = self.d1();
        let d2 = self.d2();
        let n = d2 * d1.length_squared() - d1 * d1.dot(d2);
        let len = n.length();
        if len < self.lin_tol {
            return None;
        }
        Some(n / len)
    }
}

// ============================================================================
// GeomLProp_SLProps — Surface Local Properties
// ============================================================================

/// Computes local properties of a 3D surface at (u, v):
/// point, first/second derivatives, tangents, normal, curvature analysis.
///
/// OCCT: `GeomLProp_SLProps`.
pub struct SLProps<'a> {
    surface: &'a Surface3,
    u: f64,
    v: f64,
    der_order: i32,
    cn: i32,
    lin_tol: f64,
    pnt: DVec3,
    d1u: DVec3,
    d1v: DVec3,
    d2u: DVec3,
    d2v: DVec3,
    duv: DVec3,
    normal: DVec3,
    min_curv: f64,
    max_curv: f64,
    dir_min_curv: DVec3,
    dir_max_curv: DVec3,
    mean_curv: f64,
    gaus_curv: f64,
    u_tangent_status: LPropStatus,
    v_tangent_status: LPropStatus,
    normal_status: LPropStatus,
    curvature_status: LPropStatus,
    significant_first_derivative_order_u: i32,
    significant_first_derivative_order_v: i32,
}

impl<'a> SLProps<'a> {
    /// Constructor with surface, parameter (u, v), max derivative order, and tolerance.
    ///
    /// OCCT: `GeomLProp_SLProps(Surface, U, V, N, Resolution)`.
    pub fn new(surface: &'a Surface3, u: f64, v: f64, n: i32, resolution: f64) -> Self {
        assert!(n >= 0 && n <= 2, "SLProps: N must be 0, 1, or 2");
        let mut props = SLProps {
            surface,
            u: f64::NAN,
            v: f64::NAN,
            der_order: n,
            cn: 4,
            lin_tol: if resolution <= 0.0 { TOL_LIN } else { resolution },
            pnt: DVec3::ZERO,
            d1u: DVec3::ZERO,
            d1v: DVec3::ZERO,
            d2u: DVec3::ZERO,
            d2v: DVec3::ZERO,
            duv: DVec3::ZERO,
            normal: DVec3::Z,
            min_curv: 0.0,
            max_curv: 0.0,
            dir_min_curv: DVec3::ZERO,
            dir_max_curv: DVec3::ZERO,
            mean_curv: 0.0,
            gaus_curv: 0.0,
            u_tangent_status: LPropStatus::Undecided,
            v_tangent_status: LPropStatus::Undecided,
            normal_status: LPropStatus::Undecided,
            curvature_status: LPropStatus::Undecided,
            significant_first_derivative_order_u: 0,
            significant_first_derivative_order_v: 0,
        };
        props.set_parameters(u, v);
        props
    }

    /// Constructor with surface only (call set_parameters later).
    ///
    /// OCCT: `GeomLProp_SLProps(Surface, N, Resolution)`.
    pub fn from_surface(surface: &'a Surface3, n: i32, resolution: f64) -> Self {
        assert!(n >= 0 && n <= 2, "SLProps: N must be 0, 1, or 2");
        SLProps {
            surface,
            u: f64::NAN,
            v: f64::NAN,
            der_order: n,
            cn: 4,
            lin_tol: if resolution <= 0.0 { TOL_LIN } else { resolution },
            pnt: DVec3::ZERO,
            d1u: DVec3::ZERO,
            d1v: DVec3::ZERO,
            d2u: DVec3::ZERO,
            d2v: DVec3::ZERO,
            duv: DVec3::ZERO,
            normal: DVec3::Z,
            min_curv: 0.0,
            max_curv: 0.0,
            dir_min_curv: DVec3::ZERO,
            dir_max_curv: DVec3::ZERO,
            mean_curv: 0.0,
            gaus_curv: 0.0,
            u_tangent_status: LPropStatus::Undecided,
            v_tangent_status: LPropStatus::Undecided,
            normal_status: LPropStatus::Undecided,
            curvature_status: LPropStatus::Undecided,
            significant_first_derivative_order_u: 0,
            significant_first_derivative_order_v: 0,
        }
    }

    /// Set the surface.
    ///
    /// OCCT: `SetSurface(S)`.
    pub fn set_surface(&mut self, surface: &'a Surface3) {
        self.surface = surface;
        self.cn = 4;
    }

    /// Set the parameter values and compute derivatives.
    ///
    /// OCCT: `SetParameters(U, V)`.
    pub fn set_parameters(&mut self, u: f64, v: f64) {
        self.u = u;
        self.v = v;
        self.pnt = self.surface.point_at(u, v);
        if self.der_order >= 1 {
            let (_, pu, pv) = self.surface.derivatives(u, v);
            self.d1u = pu;
            self.d1v = pv;
        }
        if self.der_order >= 2 {
            let (_, pu, pv, puu, puv, pvv) = self.surface.derivatives2(u, v);
            self.d1u = pu;
            self.d1v = pv;
            self.d2u = puu;
            self.d2v = pvv;
            self.duv = puv;
        }
        self.u_tangent_status = LPropStatus::Undecided;
        self.v_tangent_status = LPropStatus::Undecided;
        self.normal_status = LPropStatus::Undecided;
        self.curvature_status = LPropStatus::Undecided;
    }

    /// Returns the point.
    ///
    /// OCCT: `Value()`.
    pub fn value(&self) -> DVec3 {
        self.pnt
    }

    /// Returns the first U derivative.
    ///
    /// OCCT: `D1U()`.
    pub fn d1u(&self) -> DVec3 {
        self.d1u
    }

    /// Returns the first V derivative.
    ///
    /// OCCT: `D1V()`.
    pub fn d1v(&self) -> DVec3 {
        self.d1v
    }

    /// Returns the second U derivative.
    ///
    /// OCCT: `D2U()`.
    pub fn d2u(&self) -> DVec3 {
        self.d2u
    }

    /// Returns the second V derivative.
    ///
    /// OCCT: `D2V()`.
    pub fn d2v(&self) -> DVec3 {
        self.d2v
    }

    /// Returns the mixed UV derivative.
    ///
    /// OCCT: `DUV()`.
    pub fn duv(&self) -> DVec3 {
        self.duv
    }

    /// Returns true if the U tangent is defined.
    ///
    /// OCCT: `IsTangentUDefined()`.
    pub fn is_tangent_u_defined(&self) -> bool {
        self.d1u.length_squared() > self.lin_tol * self.lin_tol
    }

    /// Returns true if the V tangent is defined.
    ///
    /// OCCT: `IsTangentVDefined()`.
    pub fn is_tangent_v_defined(&self) -> bool {
        self.d1v.length_squared() > self.lin_tol * self.lin_tol
    }

    /// Returns the U tangent direction.
    ///
    /// OCCT: `TangentU(D)`.
    pub fn tangent_u(&self) -> Option<DVec3> {
        if !self.is_tangent_u_defined() {
            return None;
        }
        Some(self.d1u.normalize_or_zero())
    }

    /// Returns the V tangent direction.
    ///
    /// OCCT: `TangentV(D)`.
    pub fn tangent_v(&self) -> Option<DVec3> {
        if !self.is_tangent_v_defined() {
            return None;
        }
        Some(self.d1v.normalize_or_zero())
    }

    /// Returns true if the normal is defined.
    ///
    /// OCCT: `IsNormalDefined()`.
    pub fn is_normal_defined(&self) -> bool {
        self.d1u.cross(self.d1v).length_squared() > self.lin_tol * self.lin_tol
    }

    /// Returns the unit normal.
    ///
    /// OCCT: `Normal()`.
    pub fn normal(&mut self) -> Option<DVec3> {
        if self.normal_status != LPropStatus::Undecided && self.normal_status == LPropStatus::Defined {
            return Some(self.normal);
        }
        let n = self.d1u.cross(self.d1v);
        let len = n.length();
        if len < self.lin_tol {
            self.normal_status = LPropStatus::Zero;
            return None;
        }
        self.normal = n / len;
        self.normal_status = LPropStatus::Defined;
        Some(self.normal)
    }

    /// Returns true if curvature is defined.
    ///
    /// OCCT: `IsCurvatureDefined()`.
    pub fn is_curvature_defined(&mut self) -> bool {
        if self.curvature_status != LPropStatus::Undecided {
            return self.curvature_status == LPropStatus::Defined;
        }
        if !self.is_normal_defined() || self.der_order < 2 {
            self.curvature_status = LPropStatus::Zero;
            return false;
        }
        self.compute_curvature();
        self.curvature_status == LPropStatus::Defined
    }

    /// Returns the maximum curvature.
    ///
    /// OCCT: `MaxCurvature()`.
    pub fn max_curvature(&mut self) -> f64 {
        self.compute_curvature();
        self.max_curv
    }

    /// Returns the minimum curvature.
    ///
    /// OCCT: `MinCurvature()`.
    pub fn min_curvature(&mut self) -> f64 {
        self.compute_curvature();
        self.min_curv
    }

    /// Returns the curvature directions (max and min).
    ///
    /// OCCT: `CurvatureDirections(MaxD, MinD)`.
    pub fn curvature_directions(&mut self) -> Option<(DVec3, DVec3)> {
        if !self.is_curvature_defined() {
            return None;
        }
        Some((self.dir_max_curv, self.dir_min_curv))
    }

    /// Returns the mean curvature.
    ///
    /// OCCT: `MeanCurvature()`.
    pub fn mean_curvature(&mut self) -> f64 {
        self.compute_curvature();
        self.mean_curv
    }

    /// Returns the Gaussian curvature.
    ///
    /// OCCT: `GaussianCurvature()`.
    pub fn gaussian_curvature(&mut self) -> f64 {
        self.compute_curvature();
        self.gaus_curv
    }

    fn compute_curvature(&mut self) {
        if self.curvature_status != LPropStatus::Undecided {
            return;
        }
        if self.der_order < 2 || !self.is_normal_defined() {
            self.curvature_status = LPropStatus::Zero;
            return;
        }

        let n = self.d1u.cross(self.d1v).normalize_or_zero();
        let nu = n.length_squared();

        // First fundamental form coefficients
        let e = self.d1u.dot(self.d1u);
        let f = self.d1u.dot(self.d1v);
        let g = self.d1v.dot(self.d1v);

        let denom = e * g - f * f;
        if denom.abs() < self.lin_tol * self.lin_tol {
            self.curvature_status = LPropStatus::Zero;
            return;
        }

        // Second fundamental form coefficients
        let l = self.d2u.dot(n);
        let m = self.duv.dot(n);
        let n_val = self.d2v.dot(n);

        // Shape operator (Weingarten map) eigenvalues
        // K = (L*G - M*F) / (E*G - F^2), M' = (M*E - L*F) / (E*G - F^2)
        let k1 = (l * g - m * f) / denom;
        let k2 = (m * e - l * f) / denom;
        let k3 = (m * g - n_val * f) / denom;
        let k4 = (n_val * e - m * f) / denom;

        // Mean curvature H = (E*N + G*L - 2*F*M) / (2*(E*G - F^2))
        self.mean_curv = (e * n_val + g * l - 2.0 * f * m) / (2.0 * denom);

        // Gaussian curvature K = (L*N - M^2) / (E*G - F^2)
        self.gaus_curv = (l * n_val - m * m) / denom;

        // Principal curvatures: k = H ± sqrt(H^2 - K)
        let disc = self.mean_curv * self.mean_curv - self.gaus_curv;
        if disc > 0.0 {
            let sqrt_disc = disc.sqrt();
            self.max_curv = self.mean_curv + sqrt_disc;
            self.min_curv = self.mean_curv - sqrt_disc;

            // Principal direction for max curvature
            // (k1 - K1)*du + k2*dv = 0 → direction in UV space
            if k2.abs() > self.lin_tol {
                let du = 1.0;
                let dv = (self.max_curv - k1) / k2;
                let dir = self.d1u * du + self.d1v * dv;
                if dir.length_squared() > self.lin_tol * self.lin_tol {
                    self.dir_max_curv = dir.normalize();
                } else {
                    self.dir_max_curv = self.d1u.normalize_or_zero();
                }
            } else if k1.abs() > self.lin_tol {
                let dv = 1.0;
                let du = (self.max_curv - k4) / k1;
                let dir = self.d1u * du + self.d1v * dv;
                if dir.length_squared() > self.lin_tol * self.lin_tol {
                    self.dir_max_curv = dir.normalize();
                } else {
                    self.dir_max_curv = self.d1u.normalize_or_zero();
                }
            } else {
                self.dir_max_curv = self.d1u.normalize_or_zero();
            }

            // Min direction is orthogonal to max in the tangent plane
            let n_norm = n.normalize_or_zero();
            self.dir_min_curv = n_norm.cross(self.dir_max_curv).normalize_or_zero();

            self.curvature_status = LPropStatus::Defined;
        } else if disc.abs() < self.lin_tol {
            // Umbilic point
            self.max_curv = self.mean_curv;
            self.min_curv = self.mean_curv;
            self.dir_max_curv = self.d1u.normalize_or_zero();
            self.dir_min_curv = self.d1v.normalize_or_zero();
            self.curvature_status = LPropStatus::Defined;
        } else {
            self.curvature_status = LPropStatus::Zero;
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn test_clprops_line() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let mut props = CLProps::new(&line, 2, 1e-12);
        props.set_parameter(5.0);
        assert!((props.value() - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-12);
        let d1 = props.d1();
        assert!((d1 - DVec3::X).length() < 1e-12);
        assert!(props.is_tangent_defined());
        assert!((props.curvature() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_clprops_circle() {
        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let mut props = CLProps::new(&circle, 2, 1e-12);
        props.set_parameter(0.0);
        let p = props.value();
        assert!((p.length() - 5.0).abs() < 1e-10);
        let curv = props.curvature();
        assert!((curv - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_slprops_plane() {
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let mut props = SLProps::new(&plane, 0.0, 0.0, 2, 1e-12);
        assert!((props.value() - DVec3::ZERO).length() < 1e-12);
        assert!((props.d1u() - DVec3::X).length() < 1e-10);
        let n = props.normal();
        assert!(n.is_some());
        assert!((n.unwrap() - DVec3::Z).length() < 1e-10);
        let gc = props.gaussian_curvature();
        assert!(gc.abs() < 1e-12);
    }

    #[test]
    fn test_slprops_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: DVec3::X,
        });
        // OCCT ElSLib::SphereValue: V = latitude [-pi/2, pi/2], V=0 is the
        // equator. Evaluate there: the parametrization degenerates at the
        // poles (Pu = 0) where curvature is undefined (matching OCCT).
        let mut props = SLProps::new(&sphere, 0.0, 0.0, 2, 1e-12);
        let gc = props.gaussian_curvature();
        assert!((gc - 1.0).abs() < 1e-7); // Unit sphere: K = 1
    }
}
