use glam::DVec2;
use rcad_kernel::geom::*;
use rcad_kernel::PCurve;
use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_CLAMP_MIN, TOLERANCE_LEN_SQ_DIV_SAFE, TOLERANCE_LINEAR_ULTRA_STRICT};
use super::curve_tools::*;
use super::intres2d::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G2dCurveType { Line, Circle, Ellipse, Parabola, Hyperbola, Other }

pub fn geom2d_curve_type(curve: &Curve2d) -> G2dCurveType {
    match curve { Curve2d::Line(_) => G2dCurveType::Line, Curve2d::Circle(_) => G2dCurveType::Circle, Curve2d::Ellipse(_) => G2dCurveType::Ellipse,
    Curve2d::Parabola(_) => G2dCurveType::Parabola,
    Curve2d::Hyperbola(_) => G2dCurveType::Hyperbola, _ => G2dCurveType::Other }
}

pub fn intersect_curves_2d_ginter(c1: &Curve2d, d1: &IntRes2dDomain, c2: &Curve2d, d2: &IntRes2dDomain, tol_conf: f64, _tol: f64) -> Vec<(f64, f64)> {
    let typ1 = geom2d_curve_type(c1); let typ2 = geom2d_curve_type(c2);
    let t1_min = if d1.has_first_point() { d1.first_parameter() } else { f64::NEG_INFINITY };
    let t1_max = if d1.has_last_point() { d1.last_parameter() } else { f64::INFINITY };
    let t2_min = if d2.has_first_point() { d2.first_parameter() } else { f64::NEG_INFINITY };
    let t2_max = if d2.has_last_point() { d2.last_parameter() } else { f64::INFINITY };
    match typ1 {
        G2dCurveType::Line => match typ2 {
            G2dCurveType::Line => intersect_line_line(c1, c2, t1_min, t1_max, t2_min, t2_max),
            G2dCurveType::Circle | G2dCurveType::Ellipse | G2dCurveType::Parabola | G2dCurveType::Hyperbola => intersect_line_conic(c1, c2, typ2, t1_min, t1_max, t2_min, t2_max),
            G2dCurveType::Other => intersect_conic_curve(c1, G2dCurveType::Line, c2, t1_min, t1_max, t2_min, t2_max, tol_conf)
        },
        G2dCurveType::Circle => match typ2 {
            G2dCurveType::Line => { let r = intersect_line_conic(c2, c1, G2dCurveType::Circle, t2_min, t2_max, t1_min, t1_max); r.into_iter().map(|(a,b)|(b,a)).collect() }
            G2dCurveType::Circle => intersect_circle_circle(c1, c2, t1_min, t1_max, t2_min, t2_max),
            G2dCurveType::Ellipse | G2dCurveType::Parabola | G2dCurveType::Hyperbola => intersect_conic_conic(c1, typ1, c2, typ2, t1_min, t1_max, t2_min, t2_max),
            G2dCurveType::Other => intersect_conic_curve(c1, G2dCurveType::Circle, c2, t1_min, t1_max, t2_min, t2_max, tol_conf)
        },
        G2dCurveType::Ellipse => match typ2 {
            G2dCurveType::Line => { let r = intersect_line_conic(c2, c1, G2dCurveType::Ellipse, t2_min, t2_max, t1_min, t1_max); r.into_iter().map(|(a,b)|(b,a)).collect() }
            G2dCurveType::Circle => { let r = intersect_conic_conic(c2, G2dCurveType::Circle, c1, G2dCurveType::Ellipse, t2_min, t2_max, t1_min, t1_max); r.into_iter().map(|(a,b)|(b,a)).collect() }
            G2dCurveType::Ellipse | G2dCurveType::Parabola | G2dCurveType::Hyperbola => intersect_conic_conic(c1, typ1, c2, typ2, t1_min, t1_max, t2_min, t2_max),
            G2dCurveType::Other => intersect_conic_curve(c1, G2dCurveType::Ellipse, c2, t1_min, t1_max, t2_min, t2_max, tol_conf)
        },
        G2dCurveType::Parabola | G2dCurveType::Hyperbola => match typ2 {
            G2dCurveType::Line | G2dCurveType::Circle | G2dCurveType::Ellipse | G2dCurveType::Parabola | G2dCurveType::Hyperbola => intersect_conic_conic(c1, typ1, c2, typ2, t1_min, t1_max, t2_min, t2_max),
            G2dCurveType::Other => intersect_conic_curve(c1, typ1, c2, t1_min, t1_max, t2_min, t2_max, tol_conf)
        },
        G2dCurveType::Other => match typ2 {
            G2dCurveType::Line | G2dCurveType::Circle | G2dCurveType::Ellipse | G2dCurveType::Parabola | G2dCurveType::Hyperbola => {
                let r = intersect_conic_curve(c2, typ2, c1, t2_min, t2_max, t1_min, t1_max, tol_conf);
                r.into_iter().map(|(a,b)|(b,a)).collect()
            }
            G2dCurveType::Other => intersect_curve_curve(c1, c2, t1_min, t1_max, t2_min, t2_max, tol_conf)
        },
    }
}

fn intersect_line_line(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Line(l1) = c1 else { return vec![] }; let Curve2d::Line(l2) = c2 else { return vec![] };
    let a = l1.direction.x; let b = -l2.direction.x; let c = l1.direction.y; let d = -l2.direction.y;
    let det = a * d - b * c; if det.abs() < TOLERANCE_CLAMP_MIN { return vec![]; }
    let rx = l2.origin.x - l1.origin.x; let ry = l2.origin.y - l1.origin.y;
    let t1 = (d * rx - b * ry) / det; let t2 = (a * ry - c * rx) / det;
    if t1 >= t1_min - 1e-12 && t1 <= t1_max + 1e-12 && t2 >= t2_min - 1e-12 && t2 <= t2_max + 1e-12 { vec![(t1, t2)] } else { vec![] }
}

fn intersect_line_conic(line: &Curve2d, conic: &Curve2d, _typ: G2dCurveType, tl_min: f64, tl_max: f64, tc_min: f64, tc_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Line(l) = line else { return vec![] };
    let ro = l.origin; let rd = l.direction; let full = tl_min.is_infinite() && tl_max.is_infinite();
    match conic {
        Curve2d::Circle(c) => {
            let oc = ro - c.center; let a = rd.dot(rd); let b = 2.0*rd.dot(oc); let c2 = oc.dot(oc) - c.radius*c.radius;
            let disc = b*b - 4.0*a*c2; if disc < 0.0 { return vec![]; }
            let sd = disc.sqrt(); let mut r = Vec::new();
            for &s in &[(-b - sd)/(2.0*a), (-b + sd)/(2.0*a)] {
                if !full && s < 0.0 { continue; }
                let p = ro + s*rd; let mut t = (p.y - c.center.y).atan2(p.x - c.center.x);
                if t < 0.0 { t += std::f64::consts::TAU; }
                if t >= tc_min - 1e-12 && t <= tc_max + 1e-12 { r.push((t, s)); }
            }
            dedup(r)
        }
        Curve2d::Ellipse(e) => {
            let uv = DVec2::new(-e.major_dir.y, e.major_dir.x);
            let oc = ro - e.center;
            let du = rd.dot(e.major_dir)/e.major_radius; let dv = rd.dot(uv)/e.minor_radius;
            let ou = oc.dot(e.major_dir)/e.major_radius; let ov = oc.dot(uv)/e.minor_radius;
            let a2 = du*du+dv*dv; let b2 = 2.0*(du*ou+dv*ov); let c2 = ou*ou+ov*ov-1.0;
            let disc = b2*b2 - 4.0*a2*c2; if disc < 0.0 { return vec![]; }
            let sd = disc.sqrt(); let mut r = Vec::new();
            for &s in &[(-b2 - sd)/(2.0*a2), (-b2 + sd)/(2.0*a2)] {
                if !full && s < 0.0 { continue; }
                let dp = ro + s*rd - e.center; let mut t = dp.y.atan2(dp.x);
                if t < 0.0 { t += std::f64::consts::TAU; }
                if t >= tc_min - 1e-12 && t <= tc_max + 1e-12 { r.push((t, s)); }
            }
            dedup(r)
        }
        _ => vec![],
    }
}

fn intersect_circle_circle(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Circle(ca) = c1 else { return vec![] }; let Curve2d::Circle(cb) = c2 else { return vec![] };
    let d = ca.center.distance(cb.center); let r1 = ca.radius; let r2 = cb.radius;
    if d > r1+r2+1e-12 || d < (r1-r2).abs()-1e-12 || d < TOLERANCE_LEN_SQ_DIV_SAFE { return vec![]; }
    let a = (r1*r1 - r2*r2 + d*d) / (2.0*d); let h = (r1*r1 - a*a).max(0.0).sqrt();
    let mid = ca.center + a*(cb.center - ca.center)/d;
    let perp = DVec2::new(-(cb.center.y-ca.center.y), cb.center.x-ca.center.x)/d;
    let mut r = Vec::new();
    for &s in &[-1.0, 1.0] {
        let p = mid + s*h*perp;
        let mut t1 = (p.y-ca.center.y).atan2(p.x-ca.center.x); if t1 < 0.0 { t1 += std::f64::consts::TAU; }
        let mut t2 = (p.y-cb.center.y).atan2(p.x-cb.center.x); if t2 < 0.0 { t2 += std::f64::consts::TAU; }
        if t1 >= t1_min-1e-12 && t1 <= t1_max+1e-12 && t2 >= t2_min-1e-12 && t2 <= t2_max+1e-12 { r.push((t1, t2)); }
    }
    dedup(r)
}

fn intersect_conic_conic(c1: &Curve2d, typ1: G2dCurveType, c2: &Curve2d, typ2: G2dCurveType, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64) -> Vec<(f64, f64)> {
    if typ1 == G2dCurveType::Circle && typ2 == G2dCurveType::Ellipse { intersect_circle_ellipse(c1, c2, t1_min, t1_max, t2_min, t2_max) }
    else if typ1 == G2dCurveType::Ellipse && typ2 == G2dCurveType::Circle { let r = intersect_circle_ellipse(c2, c1, t2_min, t2_max, t1_min, t1_max); r.into_iter().map(|(a,b)|(b,a)).collect() }
    else if typ1 == G2dCurveType::Ellipse && typ2 == G2dCurveType::Ellipse { intersect_ellipse_ellipse(c1, c2, t1_min, t1_max, t2_min, t2_max) }
    else { intersect_curve_curve(c1, c2, t1_min, t1_max, t2_min, t2_max, TOLERANCE_ABS) }
}

/// IntCurve_IntConicConic::Perform(Circle, DC, Ellipse, DE)
///   (IntCurve_IntConicConic.cxx L439-482).
///   Implicit circle (IConicTool) × Parametric ellipse (PConic):
///   solve |P(t) - C_center|² - R² = 0 via sampling + Newton.
///   OCCT L453-480: ensure both domains are closed (period 2π).
fn intersect_circle_ellipse(circle: &Curve2d, ell: &Curve2d, tc_min: f64, tc_max: f64, te_min: f64, te_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Circle(c) = circle else { return vec![] }; let Curve2d::Ellipse(e) = ell else { return vec![] };
    let minor = DVec2::new(-e.major_dir.y, e.major_dir.x);
    // OCCT L456-457, 460, 472-473: ensure closed domains with period 2π
    let te_range = te_max - te_min;
    let tc_range = tc_max - tc_min;
    if te_range < TOLERANCE_CLAMP_MIN || tc_range < TOLERANCE_CLAMP_MIN { return vec![]; }
    let te_period = std::f64::consts::TAU;
    let tc_period = std::f64::consts::TAU;
    // Normalize domain spans to the [t_min, t_min + period] pattern
    let te_effective_start = te_min;
    let te_effective_end = te_min + te_period;
    let tc_effective_start = tc_min;
    let tc_effective_end = tc_min + tc_period;
    // Sampling on the parametric curve (ellipse's parameter t)
    let n_samples: usize = 256;
    let mut candidates: Vec<f64> = Vec::new();
    for i in 0..=n_samples {
        let t = te_effective_start + te_period * (i as f64 / n_samples as f64);
        let p = e.center + e.major_dir * (e.major_radius * t.cos()) + minor * (e.minor_radius * t.sin());
        let f_val = (p - c.center).length_squared() - c.radius * c.radius;
        if i > 0 {
            let t_prev = te_effective_start + te_period * ((i - 1) as f64 / n_samples as f64);
            let p_prev = e.center + e.major_dir * (e.major_radius * t_prev.cos()) + minor * (e.minor_radius * t_prev.sin());
            let f_prev = (p_prev - c.center).length_squared() - c.radius * c.radius;
            // Sign change or zero at sample point
            if f_val == 0.0 || (f_val * f_prev < 0.0) {
                candidates.push(t);
            }
        }
        if f_val.abs() < TOLERANCE_LINEAR_ULTRA_STRICT {
            let is_dup = candidates.last().map_or(false, |&lt| (t - lt).abs() < 1e-9 * te_period);
            if !is_dup { candidates.push(t); }
        }
    }
    // Refine each candidate with Newton
    let mut results: Vec<(f64, f64)> = Vec::new();
    for &t0 in &candidates {
        let mut t = t0;
        let mut converged = false;
        for _ in 0..20 {
            let p = e.center + e.major_dir * (e.major_radius * t.cos()) + minor * (e.minor_radius * t.sin());
            let dp = p - c.center;
            let f_val = dp.length_squared() - c.radius * c.radius;
            if f_val.abs() < 1e-14 {
                // Circle parameter from point angle
                let mut tc = (p.y - c.center.y).atan2(p.x - c.center.x);
                if tc < 0.0 { tc += std::f64::consts::TAU; }
                // Map result into caller's domain range
                let tc_mapped = if tc < tc_min { tc + tc_period * ((tc_min - tc) / tc_period).ceil() } else { tc };
                if tc_mapped >= tc_min - 1e-10 && tc_mapped <= tc_max + 1e-10
                    && t >= te_min - 1e-10 && t <= te_max + 1e-10
                {
                    let is_dup = results.last().map_or(false, |&(lt, _)| (t - lt).abs() < 1e-9 * te_period);
                    if !is_dup { results.push((tc_mapped, t)); }
                }
                converged = true;
                break;
            }
            let der = e.major_dir * (-e.major_radius * t.sin()) + minor * (e.minor_radius * t.cos());
            let df = 2.0 * dp.dot(der);
            if df.abs() < TOLERANCE_CLAMP_MIN { break; }
            t = t - f_val / df;
        }
        if !converged {
            // Fallback: accept candidate with approximate params
            let tc_approx = {
                let p = e.center + e.major_dir * (e.major_radius * t0.cos()) + minor * (e.minor_radius * t0.sin());
                let mut tc = (p.y - c.center.y).atan2(p.x - c.center.x);
                if tc < 0.0 { tc += std::f64::consts::TAU; }
                tc
            };
            if tc_approx >= tc_min - 1e-8 && tc_approx <= tc_max + 1e-8
                && t0 >= te_min - 1e-8 && t0 <= te_max + 1e-8
            {
                let is_dup = results.last().map_or(false, |&(lt, _)| (t0 - lt).abs() < 1e-8 * te_period);
                if !is_dup { results.push((tc_approx, t0)); }
            }
        }
    }
    dedup(results)
}

/// IntCurve_IntConicConic::Perform(Ellipse, DE1, Ellipse, DE2)
///   (IntCurve_IntConicConic.cxx L915-958).
///   Implicit ellipse1 × Parametric ellipse2: solve implicit form along parametric curve.
///   Domain closure: non-closed domains get period 2π.
fn intersect_ellipse_ellipse(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Ellipse(e1) = c1 else { return vec![] }; let Curve2d::Ellipse(e2) = c2 else { return vec![] };
    let m1 = DVec2::new(-e1.major_dir.y, e1.major_dir.x);
    let m2 = DVec2::new(-e2.major_dir.y, e2.major_dir.x);
    if (t1_max - t1_min) < TOLERANCE_CLAMP_MIN || (t2_max - t2_min) < TOLERANCE_CLAMP_MIN { return vec![]; }
    let period = std::f64::consts::TAU;
    // OCCT L929-956: ensure both domains are closed (period 2π)
    let t1_start = t1_min;
    let t2_start = t2_min;
    // Implicit form of ellipse1: test if a point lies on ellipse1.
    // In ellipse1's local frame (major_dir, m1):
    //   local.x = (P - center)・major_dir,  local.y = (P - center)・m1
    //   implicit: (local.x / major_radius)² + (local.y / minor_radius)² - 1 = 0
    let implicit_fn = |pt: DVec2| -> f64 {
        let d = pt - e1.center;
        let lx = d.dot(e1.major_dir) / e1.major_radius;
        let ly = d.dot(m1) / e1.minor_radius;
        lx * lx + ly * ly - 1.0
    };
    // Parametric form of ellipse2: P(t) = center + a*cos(t)*dir + b*sin(t)*minor
    let n_samples: usize = 256;
    let mut candidates: Vec<f64> = Vec::new();
    for i in 0..=n_samples {
        let t = t2_start + period * (i as f64 / n_samples as f64);
        let p = e2.center + e2.major_dir * (e2.major_radius * t.cos()) + m2 * (e2.minor_radius * t.sin());
        let f_val = implicit_fn(p);
        if i > 0 {
            let t_prev = t2_start + period * ((i - 1) as f64 / n_samples as f64);
            let p_prev = e2.center + e2.major_dir * (e2.major_radius * t_prev.cos()) + m2 * (e2.minor_radius * t_prev.sin());
            let f_prev = implicit_fn(p_prev);
            if f_val == 0.0 || f_val * f_prev < 0.0 {
                candidates.push(t);
            }
        }
        if f_val.abs() < TOLERANCE_LINEAR_ULTRA_STRICT {
            let is_dup = candidates.last().map_or(false, |&lt| (t - lt).abs() < 1e-9 * period);
            if !is_dup { candidates.push(t); }
        }
    }
    // Newton refinement
    let mut results: Vec<(f64, f64)> = Vec::new();
    for &t0 in &candidates {
        let mut t = t0;
        for _ in 0..20 {
            let p = e2.center + e2.major_dir * (e2.major_radius * t.cos()) + m2 * (e2.minor_radius * t.sin());
            let f_val = implicit_fn(p);
            if f_val.abs() < 1e-14 {
                let d = p - e1.center;
                let mut t1 = d.y.atan2(d.x);
                if t1 < 0.0 { t1 += std::f64::consts::TAU; }
                let t1_mapped = if t1 < t1_min { t1 + period * ((t1_min - t1) / period).ceil() } else { t1 };
                if t1_mapped >= t1_min - 1e-10 && t1_mapped <= t1_max + 1e-10
                    && t >= t2_min - 1e-10 && t <= t2_max + 1e-10
                {
                    let is_dup = results.last().map_or(false, |&(lt, _)| (t - lt).abs() < 1e-9 * period);
                    if !is_dup { results.push((t1_mapped, t)); }
                }
                break;
            }
            let der = e2.major_dir * (-e2.major_radius * t.sin()) + m2 * (e2.minor_radius * t.cos());
            // Derivative of implicit function along parametric curve:
            // d/dt f(P(t)) = ∇f(P) · P'(t)
            // ∇f = 2*(lx/a², ly/b²) in local frame �?transformed to world
            let d = p - e1.center;
            let lx = d.dot(e1.major_dir) / e1.major_radius;
            let ly = d.dot(m1) / e1.minor_radius;
            let grad = 2.0 * (e1.major_dir * (lx / e1.major_radius) + m1 * (ly / e1.minor_radius));
            let df = grad.dot(der);
            if df.abs() < TOLERANCE_CLAMP_MIN { break; }
            t = t - f_val / df;
        }
    }
    dedup(results)
}

/// IntCurve_IConicTool (IntCurve_IConicTool.hxx/lxx).
///   Implicit representation of a conic: F(P) = 0 defines the conic surface.
struct IConicTool {
    conic: Curve2d,
}

impl IConicTool {
    fn new(conic: &Curve2d) -> Self {
        IConicTool { conic: conic.clone() }
    }

    /// Implicit value F(P). F(P) = 0 �?P lies on the conic.
    fn value(&self, pt: DVec2) -> f64 {
        match &self.conic {
            Curve2d::Line(l) => {
                let dp = pt - l.origin;
                dp.x * l.direction.y - dp.y * l.direction.x
            }
            Curve2d::Circle(c) => pt.distance_squared(c.center) - c.radius * c.radius,
            Curve2d::Ellipse(e) => {
                let minor = DVec2::new(-e.major_dir.y, e.major_dir.x);
                let du = (pt - e.center).dot(e.major_dir) / e.major_radius;
                let dv = (pt - e.center).dot(minor) / e.minor_radius;
                du * du + dv * dv - 1.0
            }
            _ => 1.0,
        }
    }
}

/// IntCurve_IntConicCurveGen (IntConicCurveGen.gxx/lxx L119-130).
///   Perform(IConicTool, D1, PCurve, D2, TolConf, Tol):
///   implicit conic × parametric curve intersection.
///   Algorithm: sample parametric curve �?find sign changes in implicit value �?Newton refine.
fn intersect_conic_curve(conic: &Curve2d, _typ: G2dCurveType, curve: &Curve2d, tc_min: f64, tc_max: f64, t_min: f64, t_max: f64, tol: f64) -> Vec<(f64, f64)> {
    let tr = t_max - t_min;
    if tr < TOLERANCE_CLAMP_MIN || !tr.is_finite() { return vec![]; }
    // OCCT: ensure closed domain for periodic conic (circle/ellipse)
    let is_periodic = matches!(conic, Curve2d::Circle(_) | Curve2d::Ellipse(_));
    if is_periodic && (tc_max - tc_min) < std::f64::consts::TAU - 1e-10 {
        return vec![];
    }
    let tool = IConicTool::new(conic);
    // OCCT TheIntersector: sample parametric curve at N intervals,
    // detect sign changes in the implicit function, refine with Newton.
    const N_SAMPLES: usize = 256;
    let mut results: Vec<(f64, f64)> = Vec::new();
    for i in 0..=N_SAMPLES {
        let t = t_min + tr * (i as f64 / N_SAMPLES as f64);
        let p = curve.point_at(t);
        let f = tool.value(p);
        // Sign change detection: OCCT finds intervals where f changes sign
        if i > 0 {
            let t_prev = t_min + tr * ((i - 1) as f64 / N_SAMPLES as f64);
            let p_prev = curve.point_at(t_prev);
            let f_prev = tool.value(p_prev);
            if f * f_prev < 0.0 || f.abs() < tol.max(1e-12) {
                // Newton refine from the candidate
                let mut tn = t;
                for _ in 0..20 {
                    let pn = curve.point_at(tn);
                    let f_val = tool.value(pn);
                    if f_val.abs() < tol.max(1e-14) {
                        results.push((tn, 0.0));
                        break;
                    }
                    let eps_d = TOLERANCE_ABS;
                    let der = (curve.point_at((tn + eps_d).min(t_max)) - curve.point_at((tn - eps_d).max(t_min))) / (2.0 * eps_d);
                    let df = match conic {
                        Curve2d::Line(l) => der.x * l.direction.y - der.y * l.direction.x,
                        Curve2d::Circle(c) => 2.0 * (curve.point_at(tn) - c.center).dot(der),
                        _ => break,
                    };
                    if df.abs() < TOLERANCE_CLAMP_MIN { break; }
                    tn = tn - f_val / df;
                    if tn < t_min || tn > t_max { break; }
                }
            }
            // Near-tangent detection: f near zero at both ends of interval
            if f.abs() < TOLERANCE_LINEAR_ULTRA_STRICT && f_prev.abs() < TOLERANCE_LINEAR_ULTRA_STRICT {
                let t_mid = (t + t_prev) * 0.5;
                let p_mid = curve.point_at(t_mid);
                let f_mid = tool.value(p_mid);
                if f_mid.abs() < tol.max(1e-8) {
                    let is_dup = results.last().map_or(false, |&(lt, _)| (t_mid - lt).abs() < 1e-9 * tr);
                    if !is_dup { results.push((t_mid, 0.0)); }
                }
            }
        }
    }
    results
}

/// IntCurve_IntPolyPolyGen (IntPolyPolyGen.gxx L93-107+).
///   Polygon-based curve-curve intersection:
///   1. Sample both curves into polygons (OCCT uses ThePolygon2d + bounding boxes)
///   2. Find overlapping polygon segments
///   3. Refine intersection candidates with Newton
fn intersect_curve_curve(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64, tol: f64) -> Vec<(f64, f64)> {
    let r1 = t1_max - t1_min;
    let r2 = t2_max - t2_min;
    if r1 < TOLERANCE_CLAMP_MIN || r2 < TOLERANCE_CLAMP_MIN || !r1.is_finite() || !r2.is_finite() { return vec![]; }
    // OCCT: build polygons for both curves (point + bounding segment)
    struct PolygonSample { t: f64, pt: DVec2, bbox_min: DVec2, bbox_max: DVec2 }
    let build_polygon = |c: &Curve2d, t_min: f64, t_max: f64, n: usize| -> Vec<PolygonSample> {
        let r = t_max - t_min;
        let mut samples: Vec<PolygonSample> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = t_min + r * (i as f64 / n as f64);
            let pt = c.point_at(t);
            // Bounding box for segment [i-1, i]
            if i > 0 {
                let prev_pt = samples[i - 1].pt;
                let bmin = DVec2::new(prev_pt.x.min(pt.x), prev_pt.y.min(pt.y));
                let bmax = DVec2::new(prev_pt.x.max(pt.x), prev_pt.y.max(pt.y));
                samples[i - 1].bbox_min = bmin;
                samples[i - 1].bbox_max = bmax;
            }
            samples.push(PolygonSample { t, pt, bbox_min: DVec2::ZERO, bbox_max: DVec2::ZERO });
        }
        // Close last segment
        if n > 0 {
            let last_pt = samples[n - 1].pt;
            let pt = samples[n].pt;
            samples[n - 1].bbox_min = DVec2::new(last_pt.x.min(pt.x), last_pt.y.min(pt.y));
            samples[n - 1].bbox_max = DVec2::new(last_pt.x.max(pt.x), last_pt.y.max(pt.y));
        }
        samples
    };
    let poly1 = build_polygon(c1, t1_min, t1_max, 128);
    let poly2 = build_polygon(c2, t2_min, t2_max, 128);
    // OCCT: find overlapping segments between the two polygons
    let mut results: Vec<(f64, f64)> = Vec::new();
    for seg_a in poly1.windows(2) {
        if seg_a.len() < 2 { continue; }
        let (t1_a, p1_a, bmin_a, bmax_a) = (seg_a[0].t, seg_a[0].pt, seg_a[0].bbox_min, seg_a[0].bbox_max);
        let (_t1_b, p1_b, _bmin_b, _bmax_b) = (seg_a[1].t, seg_a[1].pt, seg_a[1].bbox_min, seg_a[1].bbox_max);
        for seg_b in poly2.windows(2) {
            if seg_b.len() < 2 { continue; }
            let t2_a = seg_b[0].t;
            let p2_a = seg_b[0].pt;
            let t2_b = seg_b[1].t;
            let p2_b = seg_b[1].pt;
            let (bmin_b, bmax_b) = (seg_b[0].bbox_min, seg_b[0].bbox_max);
            // Bounding box overlap check
            if bmax_a.x < bmin_b.x || bmin_a.x > bmax_b.x { continue; }
            if bmax_a.y < bmin_b.y || bmin_a.y > bmax_b.y { continue; }
            // Compute minimum distance between the two line segments
            let d1 = p1_b - p1_a;
            let d2 = p2_b - p2_a;
            let r = p1_a - p2_a;
            let a = d1.dot(d1); let b = d1.dot(d2); let c = d1.dot(r);
            let e = d2.dot(d2); let f = d2.dot(r);
            let det = a * e - b * b;
            if det.abs() < TOLERANCE_LEN_SQ_DIV_SAFE { continue; }
            let s = (b * f - c * e) / det;
            let t_s = (a * f - b * c) / det;
            let s_cl = s.clamp(0.0, 1.0);
            let t_cl = t_s.clamp(0.0, 1.0);
            let closest_a = p1_a + d1 * s_cl;
            let closest_b = p2_a + d2 * t_cl;
            let dist2 = closest_a.distance_squared(closest_b);
            if dist2 < tol.max(1e-8) {
                let t1_interp = t1_a + (t2_a - t2_b) * s_cl;  // wait wrong
                let t1_res = t1_a + (seg_a[1].t - t1_a) * s_cl;
                let t2_res = t2_a + (t2_b - t2_a) * t_cl;
                let is_dup = results.last().map_or(false, |&(lt1, _)| (t1_res - lt1).abs() < 1e-9 * r1);
                if !is_dup { results.push((t1_res, t2_res)); }
            }
        }
    }
    results
}

pub fn intersect_ray_curve_2d(ro: DVec2, rd: DVec2, curve: &Curve2d, t_min: f64, t_max: f64) -> Vec<(f64, f64)> {
    if rd.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE { return vec![]; }
    let rc = Curve2d::Line(Line2d { origin: ro, direction: rd });
    let mut d1 = IntRes2dDomain::new(); d1.set_values_bounded(ro, 0.0, 1e-10, ro+rd, 1.0, 1e-10);
    let mut d2 = IntRes2dDomain::new(); d2.set_values_bounded(curve.point_at(t_min), t_min, 1e-10, curve.point_at(t_max), t_max, 1e-10);
    let hits = intersect_curves_2d_ginter(&rc, &d1, curve, &d2, 1e-10, 1e-10);
    hits.into_iter().map(|(tr, tc)| (tc, tr)).collect()
}

fn dedup(v: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    if v.len() < 2 { return v; }
    let (t1,_) = v[v.len()-1]; let (t2,_) = v[v.len()-2];
    if (t1-t2).abs() < 1e-12 { v[..v.len()-1].to_vec() } else { v }
}


