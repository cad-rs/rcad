//! IntPatch_TheSurfFunction — F(u,v) = Q(P(u,v))
//!
//! OCCT IntPatch_TheSurfFunction.hxx / .cxx
//!
//! math_FunctionSetWithDerivatives: Nv=2 (u,v), Ne=1 (F=0).
//! F(u,v) = algebraic distance from P(u,v) on parametric surface to quadric Q.
//! Used by IntStart_SearchInside (interior point search) and
//! IntStart_SearchOnBoundaries (tangency computation).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};
use crate::tolerance::TOLERANCE_LEN_SQ_DIV_SAFE;
use super::super::super::inttools::int_surf_quadric::Quadric;

/// IntPatch_TheSurfFunction
///
/// Fields (hxx:78-95):
///   surf:     void* — parametric surface (Handle(Adaptor3d_Surface))
///   func:     void* — quadric function evaluator
///   u, v:     double — last evaluated UV
///   tol:      double — tolerance
///   pntsol:   gp_Pnt — last computed 3D point
///   valf:     double — last F value
///   computed: bool — whether Value has been called
///   tangent:  bool — whether gradient is near zero
///   tgdu,tgdv: double — tangent direction in UV
///   gradient: gp_Vec — gradient of F (3D)
///   derived:  bool — whether Derivatives has been called
///   d1u,d1v:  gp_Vec — surface derivatives
///   d3d:      gp_Vec — direction 3d
///   d2d:      gp_Dir2d — direction 2d
pub struct SurfFunction {
    // Parametric surface
    surf: Surface3,
    // Quadric (implicit surface)
    quad: Quadric,
    // Tolerance
    tol: f64,
    // Last evaluated UV
    u: f64,
    v: f64,
    // Last F value
    valf: f64,
    // Last computed 3D point
    pntsol: DVec3,
    // State flags
    computed: bool,
    tangent: bool,
    derived: bool,
    // Tangent direction in UV
    tgdu: f64,
    tgdv: f64,
    // Gradient of F (3D) = grad(Q) at pntsol
    gradient: DVec3,
    // Surface derivatives at last evaluated point
    d1u: DVec3,
    d1v: DVec3,
    // Direction 3D (tangent to intersection curve)
    d3d: DVec3,
    // Direction 2D (tangent in UV space)
    d2d: DVec2,
}

impl SurfFunction {
    /// OCCT L37: default constructor
    pub fn new() -> Self {
        Self {
            surf: Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
            quad: Quadric::new(),
            tol: 1e-7,
            u: 0.0, v: 0.0,
            valf: 0.0,
            pntsol: DVec3::ZERO,
            computed: false,
            tangent: false,
            derived: false,
            tgdu: 0.0, tgdv: 0.0,
            gradient: DVec3::ZERO,
            d1u: DVec3::ZERO, d1v: DVec3::ZERO,
            d3d: DVec3::ZERO,
            d2d: DVec2::ZERO,
        }
    }

    /// OCCT L40: constructor with parametric surface + quadric
    pub fn with_surface(surf: Surface3, quad: Quadric) -> Self {
        Self {
            surf,
            quad,
            ..Self::new()
        }
    }

    /// OCCT L42: constructor with quadric only
    pub fn with_quadric(quad: Quadric) -> Self {
        let mut s = Self::new();
        s.quad = quad;
        s
    }

    /// OCCT L44: Set(Handle(Adaptor3d_Surface))
    pub fn set_surface(&mut self, surf: Surface3) {
        self.surf = surf;
    }

    /// OCCT L46: SetImplicitSurface(IntSurf_Quadric)
    pub fn set_implicit_surface(&mut self, quad: Quadric) {
        self.quad = quad;
    }

    /// OCCT L48: Set(Tolerance)
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tol = tol;
    }

    /// OCCT L50: NbVariables() → 2 (u, v)
    pub fn nb_variables(&self) -> i32 { 2 }

    /// OCCT L52: NbEquations() → 1 (F=0)
    pub fn nb_equations(&self) -> i32 { 1 }

    /// OCCT L54: Value(X, F)
    /// X = [u, v], F[0] = Q(P(u,v))
    pub fn value(&mut self, x: &[f64; 2]) -> Option<f64> {
        use rcad_kernel::geom::SurfaceEval;
        self.u = x[0];
        self.v = x[1];
        self.pntsol = self.surf.point_at(self.u, self.v);
        if !self.pntsol.is_finite() {
            return None;
        }
        self.valf = self.quad.distance(self.pntsol);
        self.computed = true;
        self.derived = false;
        Some(self.valf)
    }

    /// OCCT L56: Derivatives(X, D)
    /// D[0][0] = dF/du, D[0][1] = dF/dv
    /// dF/du = grad(Q) · dP/du, dF/dv = grad(Q) · dP/dv
    pub fn derivatives(&mut self, x: &[f64; 2]) -> Option<[f64; 2]> {
        use rcad_kernel::geom::SurfaceEval;
        // Ensure Value has been called
        if !self.computed || (self.u != x[0] || self.v != x[1]) {
            self.value(x)?;
        }
        let (_, d1u, d1v) = self.surf.derivatives(self.u, self.v);
        self.d1u = d1u;
        self.d1v = d1v;
        self.gradient = self.quad.gradient(self.pntsol);
        let df_du = self.gradient.dot(self.d1u);
        let df_dv = self.gradient.dot(self.d1v);
        self.derived = true;
        Some([df_du, df_dv])
    }

    /// OCCT L58: Values(X, F, D)
    pub fn values(&mut self, x: &[f64; 2]) -> Option<(f64, [f64; 2])> {
        let f = self.value(x)?;
        let d = self.derivatives(x)?;
        Some((f, d))
    }

    /// OCCT L60: Root() → F value
    pub fn root(&self) -> f64 { self.valf }

    /// OCCT L64: Tolerance()
    pub fn tolerance(&self) -> f64 { self.tol }

    /// OCCT L66: Point() → last 3D point
    pub fn point(&self) -> &DVec3 { &self.pntsol }

    /// OCCT L68: IsTangent() — true if |grad(F)| is near zero
    pub fn is_tangent(&mut self) -> bool {
        if !self.derived {
            // Compute gradient
            let x = [self.u, self.v];
            let _ = self.derivatives(&x);
        }
        let g = self.gradient.length();
        self.tangent = g < 1e-10;
        self.tangent
    }

    /// OCCT L70: Direction3d() → tangent direction in 3D
    /// The 3D direction of the intersection curve tangent.
    pub fn direction_3d(&mut self) -> DVec3 {
        if !self.derived {
            let x = [self.u, self.v];
            let _ = self.derivatives(&x);
        }
        // Tangent = dP/du * dF/dv - dP/dv * dF/du
        // i.e., cross product of surface normal and gradient projected to 3D
        let df_du = self.gradient.dot(self.d1u);
        let df_dv = self.gradient.dot(self.d1v);
        let tg = self.d1u * df_dv - self.d1v * df_du;
        let len = tg.length();
        if len > TOLERANCE_LEN_SQ_DIV_SAFE {
            self.d3d = tg / len;
        }
        self.d3d
    }

    /// OCCT L72: Direction2d() → UV-space direction of intersection
    pub fn direction_2d(&mut self) -> DVec2 {
        if !self.derived {
            let x = [self.u, self.v];
            let _ = self.derivatives(&x);
        }
        let df_du = self.gradient.dot(self.d1u);
        let df_dv = self.gradient.dot(self.d1v);
        // 2D tangent = (dF/dv, -dF/du) — perpendicular to gradient in UV
        let len = (df_du * df_du + df_dv * df_dv).sqrt();
        if len > TOLERANCE_LEN_SQ_DIV_SAFE {
            self.d2d = DVec2::new(df_dv / len, -df_du / len);
        }
        self.d2d
    }

    /// OCCT L74: PSurface() → parametric surface
    pub fn p_surface(&self) -> &Surface3 { &self.surf }

    /// OCCT L76: ISurface() → quadric
    pub fn i_surface(&self) -> &Quadric { &self.quad }
}
