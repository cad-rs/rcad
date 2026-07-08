//! Remaining TKG3d types: AHTBezierCurve/Surface, TBezierSurface,
//! TransformedSurfaceAdaptor, GridEval surfaces.
//!
//! ✅ OCCT-aligned: standalone evaluator structs matching GeomEval_* pattern.

use glam::{DVec2, DVec3, DVec4, DAffine3};
use rcad_kernel::geom::*;

const TOL: f64 = 1e-10;
const TOL_FD: f64 = 1e-5;

// =============================================================================
// HyperboloidSurface
// =============================================================================

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
        let pt = self.eval_d0(u, v);
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
// HypParaboloidSurface (saddle)
// =============================================================================

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
        let z = u*u/(self.semi_a*self.semi_a) - v*v/(self.semi_b*self.semi_b);
        self.center + u * x_ax + v * y_ax + z * axis
    }

    pub fn eval_d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let (axis, x_ax, y_ax) = self.frame();
        let pt = self.eval_d0(u, v);
        let du = x_ax + (2.0*u/(self.semi_a*self.semi_a)) * axis;
        let dv = y_ax - (2.0*v/(self.semi_b*self.semi_b)) * axis;
        (pt, du, dv)
    }

    pub fn eval_d2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let (_pt, du, dv) = self.eval_d1(u, v);
        (self.eval_d0(u, v), du, dv,
         (2.0/(self.semi_a*self.semi_a))*self.axis,
         (-2.0/(self.semi_b*self.semi_b))*self.axis,
         DVec3::ZERO)
    }

    pub fn default_domain(&self) -> [f64; 4] {
        [f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY]
    }
}

// =============================================================================
// TBezierCurve — Trigonometric Bezier
// =============================================================================
// Basis: {1, sin(αt), cos(αt), sin(2αt), cos(2αt), ...}
// Domain: [0, π/α] (default α=1 → [0, π])

#[derive(Debug, Clone)]
pub struct TBezierCurve {
    pub poles: Vec<DVec3>,
    pub weights: Vec<f64>,
    pub alpha: f64,
}

impl TBezierCurve {
    pub fn new(poles: Vec<DVec3>, alpha: f64) -> Self {
        let w = vec![1.0; poles.len()];
        Self { poles, weights: w, alpha }
    }

    pub fn new_rational(poles: Vec<DVec3>, weights: Vec<f64>, alpha: f64) -> Self {
        Self { poles, weights, alpha }
    }

    pub fn nb_poles(&self) -> usize { self.poles.len() }
    pub fn order(&self) -> usize { (self.poles.len() - 1) / 2 }
    pub fn is_rational(&self) -> bool { self.weights.iter().any(|&w| (w - 1.0).abs() > TOL) }
    pub fn first_param(&self) -> f64 { 0.0 }
    pub fn last_param(&self) -> f64 { std::f64::consts::PI / self.alpha }

    /// EvalD0: C(t) = sum_k P_k * φ_k(t)
    /// φ_0 = 1, φ_{2k-1} = sin(kαt), φ_{2k} = cos(kαt)
    pub fn eval_d0(&self, t: f64) -> DVec3 {
        let a = self.alpha;
        let n = self.poles.len();
        let order = (n - 1) / 2;
        let mut num = DVec3::ZERO;
        let mut den = 0.0;
        if n > 0 { num += self.poles[0] * self.weights[0]; den += self.weights[0]; }
        for k in 1..=order {
            let sk = (a * t * k as f64).sin();
            let ck = (a * t * k as f64).cos();
            if 2*k-1 < n { num += self.poles[2*k-1] * self.weights[2*k-1] * sk; den += self.weights[2*k-1] * sk; }
            if 2*k < n   { num += self.poles[2*k]   * self.weights[2*k]   * ck; den += self.weights[2*k]   * ck; }
        }
        if den > TOL { num / den } else { num }
    }

    /// EvalD1: C'(t) = sum_k P_k * φ'_k(t)
    /// φ'_0 = 0, φ'_{2k-1} = kα*cos(kαt), φ'_{2k} = -kα*sin(kαt)
    pub fn eval_d1(&self, t: f64) -> (DVec3, DVec3) {
        let a = self.alpha;
        let n = self.poles.len();
        let order = (n - 1) / 2;
        let pt = self.eval_d0(t);
        let mut d1 = DVec3::ZERO;
        for k in 1..=order {
            let kf = k as f64;
            let ck = (a * t * kf).cos();
            let sk = (a * t * kf).sin();
            if 2*k-1 < n { d1 += self.poles[2*k-1] * (kf * a * ck); }
            if 2*k < n   { d1 += self.poles[2*k]   * (-kf * a * sk); }
        }
        (pt, d1)
    }
}

// =============================================================================
// TBezierSurface — Tensor-product trigonometric Bezier
// =============================================================================
// S(u,v) = sum_i sum_j P_ij * φ_i(u) * φ_j(v)
// Domain: [0, π/αu] × [0, π/αv]

#[derive(Debug, Clone)]
pub struct TBezierSurface {
    pub poles: Vec<Vec<DVec3>>,
    pub weights: Vec<Vec<f64>>,
    pub alpha_u: f64,
    pub alpha_v: f64,
}

impl TBezierSurface {
    pub fn new(poles: Vec<Vec<DVec3>>, alpha_u: f64, alpha_v: f64) -> Self {
        let w = poles.iter().map(|r| vec![1.0; r.len()]).collect();
        Self { poles, weights: w, alpha_u, alpha_v }
    }

    pub fn nb_u_poles(&self) -> usize { self.poles.len() }
    pub fn nb_v_poles(&self) -> usize { self.poles.first().map(|r| r.len()).unwrap_or(0) }
    pub fn order_u(&self) -> usize { (self.nb_u_poles() - 1) / 2 }
    pub fn order_v(&self) -> usize { (self.nb_v_poles() - 1) / 2 }

    pub fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        // Evaluate TBezier basis in U, then in V (tensor product)
        let a_u = self.alpha_u; let a_v = self.alpha_v;
        let nu = self.nb_u_poles(); let nv = self.nb_v_poles();
        let ou = (nu - 1) / 2; let ov = (nv - 1) / 2;

        // Compute basis values: bu[i] = φ_i(u), bv[j] = φ_j(v)
        let mut bu = vec![0.0; nu];
        if nu > 0 { bu[0] = 1.0; }
        for k in 1..=ou {
            if 2*k-1 < nu { bu[2*k-1] = (a_u * u * k as f64).sin(); }
            if 2*k   < nu { bu[2*k]   = (a_u * u * k as f64).cos(); }
        }
        let mut bv = vec![0.0; nv];
        if nv > 0 { bv[0] = 1.0; }
        for k in 1..=ov {
            if 2*k-1 < nv { bv[2*k-1] = (a_v * v * k as f64).sin(); }
            if 2*k   < nv { bv[2*k]   = (a_v * v * k as f64).cos(); }
        }

        let mut num = DVec3::ZERO;
        let mut den = 0.0;
        for i in 0..nu {
            for j in 0..nv {
                let b = bu[i] * bv[j];
                num += self.poles[i][j] * self.weights[i][j] * b;
                den += self.weights[i][j] * b;
            }
        }
        if den > TOL { num / den } else { num }
    }

    pub fn bounds(&self) -> [f64; 4] {
        [0.0, std::f64::consts::PI / self.alpha_u,
         0.0, std::f64::consts::PI / self.alpha_v]
    }
}

// =============================================================================
// AHTBezierCurve — Arbitrary Homogeneous Technique Bezier
// =============================================================================
// Basis: {1, t, t², ..., t^D, sinh(αt), cosh(αt), sin(βt), cos(βt)}
// D = algDegree, domain [0, 1]
// NbPoles = (D+1) + 2 + 2 = D + 5 when α>0 and β>0
// For polynomial only (α=β=0): NbPoles = D+1

#[derive(Debug, Clone)]
pub struct AHTBezierCurve {
    pub poles: Vec<DVec3>,
    pub weights: Vec<f64>,
    pub alg_degree: i32,     // polynomial degree D
    pub alpha: f64,
    pub beta: f64,
}

impl AHTBezierCurve {
    pub fn new(poles: Vec<DVec3>, alg_degree: i32, alpha: f64, beta: f64) -> Self {
        let w = vec![1.0; poles.len()];
        Self { poles, weights: w, alg_degree, alpha, beta }
    }

    pub fn new_rational(poles: Vec<DVec3>, weights: Vec<f64>, alg_degree: i32, alpha: f64, beta: f64) -> Self {
        Self { poles, weights, alg_degree, alpha, beta }
    }

    pub fn nb_poles(&self) -> usize { self.poles.len() }
    pub fn is_rational(&self) -> bool { self.weights.iter().any(|&w| (w - 1.0).abs() > TOL) }

    fn eval_basis(&self, t: f64) -> Vec<f64> {
        let mut b = Vec::new();
        // Polynomial terms: 1, t, t², ..., t^D
        let d = self.alg_degree as usize;
        for k in 0..=d { b.push(t.powi(k as i32)); }

        // Hyperbolic terms: sinh(αt), cosh(αt)
        if self.alpha > TOL {
            b.push((self.alpha * t).sinh());
            b.push((self.alpha * t).cosh());
        }
        // Trig terms: sin(βt), cos(βt)
        if self.beta > TOL {
            b.push((self.beta * t).sin());
            b.push((self.beta * t).cos());
        }
        b
    }

    pub fn eval_d0(&self, t: f64) -> DVec3 {
        let basis = self.eval_basis(t);
        let mut num = DVec3::ZERO;
        let mut den = 0.0;
        for i in 0..basis.len().min(self.poles.len()) {
            num += self.poles[i] * self.weights[i] * basis[i];
            den += self.weights[i] * basis[i];
        }
        if den > TOL { num / den } else { num }
    }
}

// =============================================================================
// AHTBezierSurface — Tensor-product AHT Bezier
// =============================================================================

#[derive(Debug, Clone)]
pub struct AHTBezierSurface {
    pub poles: Vec<Vec<DVec3>>,
    pub weights: Vec<Vec<f64>>,
    pub alg_degree_u: i32,
    pub alg_degree_v: i32,
    pub alpha_u: f64, pub beta_u: f64,
    pub alpha_v: f64, pub beta_v: f64,
}

impl AHTBezierSurface {
    pub fn new(poles: Vec<Vec<DVec3>>, du: i32, dv: i32, au: f64, bu: f64, av: f64, bv: f64) -> Self {
        let w = poles.iter().map(|r| vec![1.0; r.len()]).collect();
        Self { poles, weights: w, alg_degree_u: du, alg_degree_v: dv, alpha_u: au, beta_u: bu, alpha_v: av, beta_v: bv }
    }

    pub fn nb_poles_u(&self) -> usize { self.poles.len() }
    pub fn nb_poles_v(&self) -> usize { self.poles.first().map(|r| r.len()).unwrap_or(0) }

    pub fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        let nu = self.nb_poles_u(); let nv = self.nb_poles_v();

        // Build basis arrays
        let mut bu = Vec::new();
        for k in 0..=self.alg_degree_u as usize { bu.push(u.powi(k as i32)); }
        if self.alpha_u > TOL { bu.push((self.alpha_u * u).sinh()); bu.push((self.alpha_u * u).cosh()); }
        if self.beta_u > TOL { bu.push((self.beta_u * u).sin()); bu.push((self.beta_u * u).cos()); }

        let mut bv = Vec::new();
        for k in 0..=self.alg_degree_v as usize { bv.push(v.powi(k as i32)); }
        if self.alpha_v > TOL { bv.push((self.alpha_v * v).sinh()); bv.push((self.alpha_v * v).cosh()); }
        if self.beta_v > TOL { bv.push((self.beta_v * v).sin()); bv.push((self.beta_v * v).cos()); }

        let mut num = DVec3::ZERO;
        let mut den = 0.0;
        for i in 0..nu.min(bu.len()) {
            for j in 0..nv.min(bv.len()) {
                let b = bu[i] * bv[j];
                num += self.poles[i][j] * self.weights[i][j] * b;
                den += self.weights[i][j] * b;
            }
        }
        if den > TOL { num / den } else { num }
    }

    pub fn bounds(&self) -> [f64; 4] { [0.0, 1.0, 0.0, 1.0] }
}

// =============================================================================
// TransformedSurfaceAdaptor — GeomAdaptor_TransformedSurface
// =============================================================================

#[derive(Debug, Clone)]
pub struct TransformedSurfaceAdaptor {
    surface: Surface3,
    trsf: DAffine3,
    transformed_surface: Option<Surface3>,
}

impl TransformedSurfaceAdaptor {
    pub fn new(surface: Surface3, trsf: DAffine3) -> Self {
        let transformed = if trsf == DAffine3::IDENTITY { Some(surface.clone()) } else { None };
        Self { surface, trsf, transformed_surface: transformed }
    }

    pub fn geom_surface_original(&self) -> &Surface3 { &self.surface }
    pub fn trsf(&self) -> DAffine3 { self.trsf }

    pub fn set_trsf(&mut self, trsf: DAffine3) {
        self.trsf = trsf;
        self.transformed_surface = None; // invalidate cache
    }

    pub fn evaluate(&self, u: f64, v: f64) -> DVec3 {
        let pt = self.surface.point_at(u, v);
        self.trsf.transform_point3(pt)
    }
}

// =============================================================================
// GridEval for surfaces
// =============================================================================

/// Batch evaluate a plane at uniform grid.
pub fn batch_eval_plane_grid(plane: &Plane, u_params: &[f64], v_params: &[f64]) -> Vec<Vec<DVec3>> {
    let normal = plane.normal;
    let x_ax = any_perpendicular(normal);
    let y_ax = normal.cross(x_ax);
    u_params.iter().map(|&u| {
        v_params.iter().map(|&v| plane.origin + u * x_ax + v * y_ax).collect()
    }).collect()
}

/// Batch evaluate a sphere at uniform grid.
pub fn batch_eval_sphere_grid(sphere: &SphericalSurface, u_params: &[f64], v_params: &[f64]) -> Vec<Vec<DVec3>> {
    let c = sphere.center;
    let a = sphere.axis;
    let x = sphere.ref_dir;
    let y = a.cross(x);
    let r = sphere.radius;
    u_params.iter().map(|&u| {
        let (cu, su) = (u.cos(), u.sin());
        v_params.iter().map(|&v| {
            let (cv, sv) = (v.cos(), v.sin());
            c + r * (cv * (cu * x + su * y) + sv * a)
        }).collect()
    }).collect()
}

/// Batch evaluate any surface at uniform grid (generic fallback).
pub fn batch_eval_surface_grid(surface: &Surface3, u_params: &[f64], v_params: &[f64]) -> Vec<Vec<DVec3>> {
    u_params.iter().map(|&u| {
        v_params.iter().map(|&v| surface.point_at(u, v)).collect()
    }).collect()
}

fn any_perpendicular(v: DVec3) -> DVec3 {
    if v.x.abs() > v.y.abs() {
        DVec3::new(-v.z, 0.0, v.x).normalize()
    } else {
        DVec3::new(0.0, v.z, -v.y).normalize()
    }
}
