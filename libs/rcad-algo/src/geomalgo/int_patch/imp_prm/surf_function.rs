// OCCT IntImp_ZerImpFunc.gxx / .lxx (IntPatch_TheSurfFunction) 1:1 Rust
// translation — the algebraic distance function F(u,v) = Q(P(u,v)) used by
// IntStart_SearchInside and IntWalk_IWalking.
//
// math_FunctionSetWithDerivatives: Nv = 2 (u,v), Ne = 1 (F = 0).
//
// OCCT IntImp_ZerImpFunc.gxx L1-131, IntImp_ZerImpFunc.lxx L20-69.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};

use crate::geomalgo::int_surf::quadric::Quadric;

/// OCCT IntPatch_TheSurfFunction (IntImp_ZerImpFunc) — algebraic distance
/// from the parametric-surface point P(u,v) to the implicit quadric Q.
#[derive(Clone)]
pub struct SurfFunction {
    // OCCT: surf (void* — ThePSurface), func (void* — TheISurface).
    surf: Surface3,
    quad: Quadric,
    // OCCT: u, v — last evaluated parameters.
    u: f64,
    v: f64,
    // OCCT: tol — the tolerance below which F is considered null.
    tol: f64,
    // OCCT: pntsol — last computed surface point.
    pntsol: DVec3,
    // OCCT: valf — last F value.
    valf: f64,
    // OCCT: computed — one-shot guard for IsTangent.
    computed: bool,
    // OCCT: tangent — cached IsTangent result.
    tangent: bool,
    // OCCT: tgdu, tgdv — tangent direction in the parametric space.
    tgdu: f64,
    tgdv: f64,
    // OCCT: gradient — gradient of F in 3D.
    gradient: DVec3,
    // OCCT: derived — whether d1u/d1v/gradient are fresh.
    derived: bool,
    // OCCT: d1u, d1v — surface first derivatives.
    d1u: DVec3,
    d1v: DVec3,
    // OCCT: d3d — intersection tangent in 3D.
    d3d: DVec3,
    // OCCT: d2d — intersection tangent in 2D.
    d2d: DVec2,
}

impl SurfFunction {
    /// OCCT IntImp_ZerImpFunc() (L29-40).
    pub fn new() -> Self {
        SurfFunction {
            surf: Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
            quad: Quadric::new(),
            u: 0.0,
            v: 0.0,
            tol: 0.0,
            pntsol: DVec3::ZERO,
            valf: 0.0,
            computed: false,
            tangent: false,
            tgdu: 0.0,
            tgdv: 0.0,
            gradient: DVec3::ZERO,
            derived: false,
            d1u: DVec3::ZERO,
            d1v: DVec3::ZERO,
            d3d: DVec3::ZERO,
            d2d: DVec2::ZERO,
        }
    }

    /// OCCT IntImp_ZerImpFunc(PS, IS) (L42-53).
    #[allow(dead_code)]
    pub fn with_surface(surf: Surface3, quad: Quadric) -> Self {
        let mut s = SurfFunction::new();
        s.surf = surf;
        s.quad = quad;
        s
    }

    /// OCCT IntImp_ZerImpFunc(IS) (L55-63).
    pub fn with_quadric(quad: Quadric) -> Self {
        let mut s = SurfFunction::new();
        s.quad = quad;
        s
    }

    /// OCCT Set(PS) (lxx L26-30).
    pub fn set_surface(&mut self, surf: Surface3) {
        self.surf = surf;
    }

    /// OCCT SetImplicitSurface(IS) (lxx L32-36).
    pub fn set_implicit_surface(&mut self, quad: Quadric) {
        self.quad = quad;
    }

    /// OCCT Set(Tol) (lxx L38-41).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tol = tol;
    }

    /// OCCT NbVariables() — always 2.
    pub fn nb_variables(&self) -> i32 {
        2
    }

    /// OCCT NbEquations() — always 1.
    pub fn nb_equations(&self) -> i32 {
        1
    }

    /// OCCT Value(X, F) (gxx L77-89): F(1) = Q(P(u,v)).
    pub fn value(&mut self, x: &[f64; 2]) -> Option<f64> {
        self.u = x[0];
        self.v = x[1];
        self.pntsol = self.surf.point_at(self.u, self.v);
        self.valf = self.quad.distance(self.pntsol);
        self.computed = false;
        self.derived = false;
        Some(self.valf)
    }

    /// OCCT Derivatives(X, D) (gxx L91-103): D(1,1)=d1u·grad, D(1,2)=d1v·grad.
    pub fn derivatives(&mut self, x: &[f64; 2]) -> Option<[f64; 2]> {
        self.u = x[0];
        self.v = x[1];
        let (p, d1u, d1v) = self.surf.derivatives(self.u, self.v);
        self.pntsol = p;
        self.d1u = d1u;
        self.d1v = d1v;
        self.gradient = self.quad.gradient(self.pntsol);
        self.computed = false;
        self.derived = true;
        Some([self.d1u.dot(self.gradient), self.d1v.dot(self.gradient)])
    }

    /// OCCT Values(X, F, D) (gxx L105-120): D1 + ValueAndGradient.
    pub fn values(&mut self, x: &[f64; 2]) -> Option<(f64, [f64; 2])> {
        self.u = x[0];
        self.v = x[1];
        let (p, d1u, d1v) = self.surf.derivatives(self.u, self.v);
        self.pntsol = p;
        self.d1u = d1u;
        self.d1v = d1v;
        let (valf, gradient) = self.quad.val_and_grad(self.pntsol);
        self.valf = valf;
        self.gradient = gradient;
        self.computed = false;
        self.derived = true;
        Some((self.valf, [self.d1u.dot(self.gradient), self.d1v.dot(self.gradient)]))
    }

    /// OCCT IsTangent() (gxx L122-141).
    pub fn is_tangent(&mut self) -> bool {
        const EPS_ANG2: f64 = 1e-16;
        const TOLPETIT: f64 = 1e-16;
        if !self.computed {
            self.computed = true;
            if !self.derived {
                let (p, d1u, d1v) = self.surf.derivatives(self.u, self.v);
                self.pntsol = p;
                self.d1u = d1u;
                self.d1v = d1v;
                self.derived = true;
            }
            self.tgdu = self.gradient.dot(self.d1v);
            self.tgdv = -self.gradient.dot(self.d1u);
            let n2grad = self.gradient.length_squared();
            let n2grad_eps_ang2 = n2grad * EPS_ANG2;
            let n2d1u = self.d1u.length_squared();
            let n2d1v = self.d1v.length_squared();
            self.tangent = (self.tgdu * self.tgdu <= n2grad_eps_ang2 * n2d1v)
                && (self.tgdv * self.tgdv <= n2grad_eps_ang2 * n2d1u);
            if !self.tangent {
                // d3d.SetLinearForm(tgdu, d1u, tgdv, d1v) = tgdu*d1u + tgdv*d1v.
                self.d3d = self.d1u * self.tgdu + self.d1v * self.tgdv;
                // d2d = gp_Dir2d(tgdu, tgdv) — normalized.
                self.d2d = DVec2::new(self.tgdu, self.tgdv).normalize_or_zero();
                if self.d3d.length() <= TOLPETIT {
                    self.tangent = true;
                }
            }
        }
        self.tangent
    }

    /// OCCT Root() (lxx L43-45) — last F value.
    pub fn root(&self) -> f64 {
        self.valf
    }

    /// OCCT Tolerance() (lxx L47-49).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// OCCT Point() (lxx L51-53).
    pub fn point(&self) -> DVec3 {
        self.pntsol
    }

    /// OCCT Direction3d() (lxx L55-59) — throws when tangent; returns d3d.
    pub fn direction_3d(&mut self) -> DVec3 {
        if self.is_tangent() {
            // OCCT throws StdFail_UndefinedDerivative.  rcad: return zero.
            return DVec3::ZERO;
        }
        self.d3d
    }

    /// OCCT Direction2d() (lxx L61-65).
    pub fn direction_2d(&mut self) -> DVec2 {
        if self.is_tangent() {
            return DVec2::ZERO;
        }
        self.d2d
    }

    /// OCCT GetStateNumber() (math_FunctionSetWithDerivatives).
    pub fn get_state_number(&mut self) -> i32 {
        0
    }

    /// OCCT PSurface() (lxx L67-70).
    pub fn p_surface(&self) -> &Surface3 {
        &self.surf
    }

    /// OCCT ISurface() (lxx L72-75).
    pub fn i_surface(&self) -> &Quadric {
        &self.quad
    }
}
