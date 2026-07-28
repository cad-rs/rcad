//! IntPatch_ArcFunction — F(t) = Q(C(t)).
//!
//! OCCT IntPatch_ArcFunction.hxx / .cxx
//!
//! Given a quadric surface Q and a 2D curve C(t) on a parametric surface S,
//! evaluates F(t) = distance from Q to S(C(t)) — i.e., the algebraic distance
//! from the 3D point on the parametric surface to the quadric.
//! Used by IntStart_SearchOnBoundaries to find zero crossings on boundary curves.

use super::super::super::int_tools::int_surf_quadric::Quadric;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, SurfaceEval};

/// IntPatch_ArcFunction = math_FunctionWithDerivative
///
/// Fields (hxx:62-67):
///   myArc:  Handle(Adaptor2d_Curve2d) — boundary curve
///   mySurf: Handle(Adaptor3d_Surface) — parametric surface
///   myQuad: IntSurf_Quadric          — quadric (analytic surface)
///   ptsol:  gp_Pnt                   — last computed point
///   seqpt:  Sequence of gp_Pnt       — sample points
pub struct ArcFunction {
    // rcad: curve on parametric surface, as a 2D curve in the surface's parameter space
    arc: Curve2d,
    // rcad: the parametric surface (Surface3 stores both geometry and domain)
    surf: Surface3,
    // rcad: the quadric (analytic surface)
    quad: Quadric,
    // Last computed 3D point
    ptsol: DVec3,
    // Sample points from evaluation
    seqpt: Vec<DVec3>,
}

impl ArcFunction {
    /// OCCT L32: default constructor
    pub fn new() -> Self {
        Self {
            arc: Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::X,
            }),
            surf: Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
            quad: Quadric::new(),
            ptsol: DVec3::ZERO,
            seqpt: Vec::new(),
        }
    }

    /// OCCT L34: SetQuadric
    pub fn set_quadric(&mut self, q: Quadric) {
        self.quad = q;
    }

    /// OCCT L38: Set(Handle(Adaptor3d_Surface))
    pub fn set_surface(&mut self, surf: Surface3) {
        self.surf = surf;
    }

    /// OCCT L36: Set(Handle(Adaptor2d_Curve2d))
    pub fn set_arc(&mut self, arc: Curve2d) {
        self.arc = arc;
    }

    /// OCCT L40: Value(X, F) → bool
    /// F(t) = Q(P(C(t))) where P maps 2D→3D on the parametric surface
    pub fn value(&mut self, x: f64) -> Option<f64> {
        let p2d = self.arc.point_at(x);
        let p3d = self.surf.point_at(p2d.x, p2d.y);
        if !p3d.is_finite() {
            return None;
        }
        self.ptsol = p3d;
        let f = self.quad.distance(p3d);
        self.seqpt.push(p3d);
        Some(f)
    }

    /// OCCT L42: Derivative(X, D) → bool
    /// D = dF/dt = grad(Q) · dP/dt  where dP/dt = dP/du * du/dt + dP/dv * dv/dt
    pub fn derivative(&self, x: f64) -> Option<f64> {
        let p2d = self.arc.point_at(x);
        let (_, d1u, d1v) = self.surf.derivatives(p2d.x, p2d.y);
        // 2D tangent via finite difference since Curve2d has no tangent_at
        let eps = 1e-7;
        let p2d_eps = self.arc.point_at(x + eps);
        let d2d = (p2d_eps - p2d) / eps;
        // dP/dt = dP/du * du/dt + dP/dv * dv/dt
        let dp_dt = d1u * d2d.x + d1v * d2d.y;
        let grad_q = self.quad.gradient(self.ptsol);
        Some(grad_q.dot(dp_dt))
    }

    /// OCCT L44: Values(X, F, D) → bool
    pub fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let f = self.value(x)?;
        let d = self.derivative(x)?;
        Some((f, d))
    }

    /// OCCT L46: NbSamples() → int
    pub fn nb_samples(&self) -> i32 {
        10 // OCTC default
    }

    /// OCCT L48: GetStateNumber
    pub fn get_state_number(&mut self) -> i32 {
        0
    }

    /// OCCT L50: Valpoint(Index) → const gp_Pnt&
    pub fn valpoint(&self, index: usize) -> &DVec3 {
        &self.seqpt[index]
    }

    /// OCCT L52: Quadric()
    pub fn quadric(&self) -> &Quadric {
        &self.quad
    }

    /// OCCT L54: Arc()
    pub fn arc(&self) -> &Curve2d {
        &self.arc
    }

    /// OCCT L56: Surface()
    pub fn surface(&self) -> &Surface3 {
        &self.surf
    }

    /// OCCT L60: LastComputedPoint()
    pub fn last_computed_point(&self) -> &DVec3 {
        &self.ptsol
    }
}
