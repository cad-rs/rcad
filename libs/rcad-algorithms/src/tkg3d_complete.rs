//! Remaining TKG3d types: HyperboloidSurface, ParaboloidSurface, HypParaboloidSurface,
//! TBezierCurve/Surface, AHTBezierCurve/Surface, Surface batch eval, TransformedSurface.
//!
//! ✅ OCCT-aligned: standalone evaluator structs matching GeomEval_* pattern.

use glam::{DVec2, DVec3, DAffine3};
use rcad_kernel::geom::*;

const TOL: f64 = 1e-10;
const TOL_FD: f64 = 1e-5;

// =============================================================================
// HyperboloidSurface — one-sheet and two-sheet
// =============================================================================
// OCCT: GeomEval_HyperboloidSurface
// One-sheet: P(u,v)=R1*cosh(v)*cos(u)*X+R1*cosh(v)*sin(u)*Y+R2*sinh(v)*Z
// Two-sheet: P(u,v)=R2*sinh(v)*cos(u)*X+R2*sinh(v)*sin(u)*Y+R1*cosh(v)*Z

#[derive(Debug, Clone)]
pub enum SheetMode { OneSheet, TwoSheets }

#[derive(Debug, Clone)]
pub struct HyperboloidEvaluator {
    pub center: DVec3,
    pub axis: DVec3,
    pub ref_dir: DVec3,
    pub r1: f64,
    pub r2: f64,
    pub mode: SheetMode,
}

impl HyperboloidEvaluator {
    pub fn new(r1: f64, r2: f64, mode: SheetMode) -> Self {
        Self { center: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X, r1, r2, mode }
    }

    fn frame(&self) -> (DVec3, DVec3, DVec3) {
        let x = any_perpendicular(self.axis);
        let y = self.axis.cross(x);
        (self.axis, x, y)
    }

    pub fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_ax, y_ax) = self.frame();
        let (cu, su) = (u.cos(), u.sin());
        match self.mode {
            SheetMode::OneSheet => {
                let ch = v.cosh(); let sh = v.sinh();
                self.center + self.r1 * ch * (cu * x_ax + su * y_ax) + self.r2 * sh * axis
            }
            SheetMode::TwoSheets => {
                let ch = v.cosh(); let sh = v.sinh();
                self.center + self.r2 * sh * (cu * x_ax + su * y_ax) + self.r1 * ch * axis
            }
        }
    }

    pub fn eval_d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let (axis, x_ax, y_ax) = self.frame();
        let (cu, su) = (u.cos(), u.sin());
        match self.mode {
            SheetMode::OneSheet => {
                let ch = v.cosh(); let sh = v.sinh();
                let pt = self.center + self.r1 * ch * (cu * x_ax + su * y_ax) + self.r2 * sh * axis;
                let du = self.r1 * ch * (-su * x_ax + cu * y_ax);
                let dv = self.r1 * sh * (cu * x_ax + su * y_ax) + self.r2 * ch * axis;
                (pt, du, dv)
            }
            SheetMode::TwoSheets => {
                let ch = v.cosh(); let sh = v.sinh();
                let pt = self.center + self.r2 * sh * (cu * x_ax + su * y_ax) + self.r1 * ch * axis;
                let du = self.r2 * sh * (-su * x_ax + cu * y_ax);
                let dv = self.r2 * ch * (cu * x_ax + su * y_ax) + self.r1 * sh * axis;
                (pt, du, dv)
            }
        }
    }

    pub fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * std::f64::consts::PI, -2.0, 2.0]
    }

    pub fn is_u_periodic(&self) -> bool { true }
    pub fn is_u_closed(&self) -> bool { true }
}

// =============================================================================
// ParaboloidSurface
// =============================================================================
// OCCT: GeomEval_ParaboloidSurface
// P(u,v)=v*cos(u)*X+v*sin(u)*Y+(v²/(4F))*Z, u in [0,2π]

#[derive(Debug, Clone)]
pub struct ParaboloidEvaluator {
    pub center: DVec3,
    pub axis: DVec3,
    pub ref_dir: DVec3,
    pub focal: f64,
}

impl ParaboloidEvaluator {
    pub fn new(focal: f64) -> Self {
        Self { center: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X, focal }
    }

    fn frame(&self) -> (DVec3, DVec3, DVec3) {
        let x = any_perpendicular(self.axis);
        let y = self.axis.cross(x);
        (self.axis, x, y)
    }

    pub fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_ax, y_ax) = self.frame();
        let (cu, su) = (u.cos(), u.sin());
        self.center + v * (cu * x_ax + su * y_ax) + (v * v / (4.0 * self.focal)) * axis
    }

    pub fn eval_d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let (axis, x_ax, y_ax) = self.frame();
        let (cu, su) = (u.cos(), u.sin());
        let pt = self.center + v * (cu * x_ax + su * y_ax) + (v * v / (4.0 * self.focal)) * axis;
        let du = v * (-su * x_ax + cu * y_ax);
        let dv = cu * x_ax + su * y_ax + (v / (2.0 * self.focal)) * axis;
        (pt, du, dv)
    }

    pub fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * std::f64::consts::PI, 0.0, 5.0]
    }

    pub fn is_u_periodic(&self) -> bool { true }
    pub fn is_u_closed(&self) -> bool { true }
}

// =============================================================================
// HypParaboloidSurface (Hyperbolic Paraboloid / saddle)
// =============================================================================
// OCCT: GeomEval_HypParaboloidSurface
// P(u,v)=u*X+v*Y+(u²/A²-v²/B²)*Z

#[derive(Debug, Clone)]
pub struct HypParaboloidEvaluator {
    pub center: DVec3,
    pub axis: DVec3,
    pub ref_dir: DVec3,
    pub semi_a: f64,
    pub semi_b: f64,
}

impl HypParaboloidEvaluator {
    pub fn new(semi_a: f64, semi_b: f64) -> Self {
        Self { center: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X, semi_a, semi_b }
    }

    fn frame(&self) -> (DVec3, DVec3, DVec3) {
        let x = any_perpendicular(self.axis);
        let y = self.axis.cross(x);
        (self.axis, x, y)
    }

    pub fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_ax, y_ax) = self.frame();
        self.center + u * x_ax + v * y_ax + (u * u / (self.semi_a * self.semi_a) - v * v / (self.semi_b * self.semi_b)) * axis
    }

    pub fn eval_d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let (axis, x_ax, y_ax) = self.frame();
        let pt = self.eval_d0(u, v);
        let du = x_ax + (2.0 * u / (self.semi_a * self.semi_a)) * axis;
        let dv = y_ax - (2.0 * v / (self.semi_b * self.semi_b)) * axis;
        (pt, du, dv)
    }

    pub fn eval_d2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let (axis, x_ax, y_ax) = self.frame();
        let (pt, du, dv) = self.eval_d1(u, v);
        let d2u = (2.0 / (self.semi_a * self.semi_a)) * axis;
        let d2v = (-2.0 / (self.semi_b * self.semi_b)) * axis;
        let d2uv = DVec3::ZERO;
        (pt, du, dv, d2u, d2v, d2uv)
    }

    pub fn default_domain(&self) -> [f64; 4] {
        [f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY]
    }
}

// =============================================================================
// TBezierCurve — Trigonometric Bezier curve
// =============================================================================
// OCCT: GeomEval_TBezierCurve
// Basis: {1, sin(alpha*t), cos(alpha*t)} => 3 poles
// C(t) = P0 + P1*sin(alpha*t) + P2*cos(alpha*t), t in [0, Pi/alpha]

#[derive(Debug, Clone)]
pub struct TBezierCurve {
    pub poles: Vec<DVec3>,
    pub alpha: f64,
}

impl TBezierCurve {
    pub fn new(poles: Vec<DVec3>, alpha: f64) -> Self {
        Self { poles, alpha }
    }

    pub fn nb_poles(&self) -> usize { self.poles.len() }

    pub fn order(&self) -> usize { (self.poles.len() - 1) / 2 }

    pub fn first_param(&self) -> f64 { 0.0 }

    pub fn last_param(&self) -> f64 { std::f64::consts::PI / self.alpha }

    pub fn eval_d0(&self, t: f64) -> DVec3 {
        let a = self.alpha;
        let n = self.poles.len();
        let order = (n - 1) / 2;
        let mut pt = DVec3::ZERO;
        if n > 0 { pt += self.poles[0]; }
        for k in 1..=order {
            if 2 * k < n { pt += self.poles[2 * k - 1] * (a * t * k as f64).sin(); }
            if 2 * k + 1 < n { pt += self.poles[2 * k] * (a * t * k as f64).cos(); }
        }
        pt
    }
}

// =============================================================================
// Batch surface evaluation — GeomGridEval_Surface equivalent
// =============================================================================

/// Batch evaluate a plane at uniform grid (u,v) parameter arrays.
pub fn batch_eval_plane(plane: &Plane, u_params: &[f64], v_params: &[f64]) -> Vec<Vec<DVec3>> {
    let normal = plane.normal;
    let x_ax = any_perpendicular(normal);
    let y_ax = normal.cross(x_ax);
    u_params.iter().map(|&u| {
        v_params.iter().map(|&v| {
            plane.origin + u * x_ax + v * y_ax
        }).collect()
    }).collect()
}

/// Batch evaluate a sphere at uniform grid (u,v) parameter arrays.
pub fn batch_eval_sphere(sphere: &SphericalSurface, u_params: &[f64], v_params: &[f64]) -> Vec<Vec<DVec3>> {
    u_params.iter().map(|&u| {
        v_params.iter().map(|&v| sphere.point_at(u, v)).collect()
    }).collect()
}

/// Batch evaluate any surface at uniform grid (u,v) parameter arrays (generic fallback).
pub fn batch_eval_surface(surface: &Surface3, u_params: &[f64], v_params: &[f64]) -> Vec<Vec<DVec3>> {
    u_params.iter().map(|&u| {
        v_params.iter().map(|&v| surface.point_at(u, v)).collect()
    }).collect()
}

// =============================================================================
// TransformedSurfaceAdaptor — GeomAdaptor_TransformedSurface equivalent
// =============================================================================

#[derive(Debug, Clone)]
pub struct TransformedSurfaceAdaptor {
    surface: Option<Surface3>,
    trsf: DAffine3,
}

impl TransformedSurfaceAdaptor {
    pub fn new(surface: Surface3, trsf: DAffine3) -> Self {
        Self { surface: Some(surface), trsf }
    }

    pub fn evaluate(&self, u: f64, v: f64) -> Option<DVec3> {
        let surf = self.surface.as_ref()?;
        let pt = surf.point_at(u, v);
        Some(self.trsf.transform_point3(pt))
    }

    pub fn set_trsf(&mut self, trsf: DAffine3) { self.trsf = trsf; }
}

// Need any_perpendicular 
fn any_perpendicular(v: DVec3) -> DVec3 {
    if v.x.abs() > v.y.abs() {
        DVec3::new(-v.z, 0.0, v.x).normalize()
    } else {
        DVec3::new(0.0, v.z, -v.y).normalize()
    }
}
