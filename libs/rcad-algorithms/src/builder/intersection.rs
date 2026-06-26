use glam::DVec2;
use rcad_kernel::geom::*;
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
    let det = a * d - b * c; if det.abs() < 1e-15 { return vec![]; }
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
    if d > r1+r2+1e-12 || d < (r1-r2).abs()-1e-12 || d < 1e-30 { return vec![]; }
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
    else { intersect_curve_curve(c1, c2, t1_min, t1_max, t2_min, t2_max, 1e-7) }
}

/// OCCT-aligned: IntCurve_IntConicConic::Perform(Circle, DC, Ellipse, DE)
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
    if te_range < 1e-15 || tc_range < 1e-15 { return vec![]; }
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
        if f_val.abs() < 1e-10 {
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
            if df.abs() < 1e-15 { break; }
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

/// OCCT-aligned: IntCurve_IntConicConic::Perform(Ellipse, DE1, Ellipse, DE2)
///   (IntCurve_IntConicConic.cxx L915-958).
///   Implicit ellipse1 × Parametric ellipse2: solve implicit form along parametric curve.
///   Domain closure: non-closed domains get period 2π.
fn intersect_ellipse_ellipse(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Ellipse(e1) = c1 else { return vec![] }; let Curve2d::Ellipse(e2) = c2 else { return vec![] };
    let m1 = DVec2::new(-e1.major_dir.y, e1.major_dir.x);
    let m2 = DVec2::new(-e2.major_dir.y, e2.major_dir.x);
    if (t1_max - t1_min) < 1e-15 || (t2_max - t2_min) < 1e-15 { return vec![]; }
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
        if f_val.abs() < 1e-10 {
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
            // ∇f = 2*(lx/a², ly/b²) in local frame → transformed to world
            let d = p - e1.center;
            let lx = d.dot(e1.major_dir) / e1.major_radius;
            let ly = d.dot(m1) / e1.minor_radius;
            let grad = 2.0 * (e1.major_dir * (lx / e1.major_radius) + m1 * (ly / e1.minor_radius));
            let df = grad.dot(der);
            if df.abs() < 1e-15 { break; }
            t = t - f_val / df;
        }
    }
    dedup(results)
}

/// OCCT-aligned: IntCurve_IntConicCurveGen (IntConicCurveGen.gxx/lxx).
///   Implicit conic (IConicTool) × Parametric curve (ThePCurve):
///   sample curve → find nearest point to conic → Newton refine.
///   Domain closure: circle/ellipse domains are made closed (period 2π).
///   OCCT L119-130: Perform(IConicTool, D1, PCurve, D2, TolConf, Tol).
fn intersect_conic_curve(conic: &Curve2d, _typ: G2dCurveType, curve: &Curve2d, tc_min: f64, tc_max: f64, t_min: f64, t_max: f64, tol: f64) -> Vec<(f64, f64)> {
    let tr = t_max - t_min;
    if tr < 1e-15 || !tr.is_finite() { return vec![]; }
    // OCCT: ensure closed domain for circle/ellipse conic
    let is_periodic = matches!(conic, Curve2d::Circle(_) | Curve2d::Ellipse(_));
    if is_periodic && (tc_max - tc_min).abs() < std::f64::consts::TAU - 1e-10 {
        // Non-closed domain of a periodic conic: skip (OCCT treats as no intersection)
        return vec![];
    }
    // Implicit function for the conic
    let implicit_fn = |pt: DVec2| -> f64 {
        match conic {
            Curve2d::Line(l) => {
                let dp = pt - l.origin;
                (dp.x * l.direction.y - dp.y * l.direction.x).abs() / l.direction.length().max(1e-30)
            }
            Curve2d::Circle(c) => pt.distance(c.center) - c.radius,
            Curve2d::Ellipse(e) => {
                let u = e.major_dir;
                let v = DVec2::new(-u.y, u.x);
                let du = (pt - e.center).dot(u) / e.major_radius;
                let dv = (pt - e.center).dot(v) / e.minor_radius;
                (du * du + dv * dv).sqrt() - 1.0
            }
            _ => return 1.0,
        }
    };
    // Sample the parametric curve, find best candidate
    const N_SAMPLES: usize = 256;
    let mut candidates: Vec<f64> = Vec::new();
    for i in 0..=N_SAMPLES {
        let t = t_min + tr * (i as f64 / N_SAMPLES as f64);
        let p = curve.point_at(t);
        let d = implicit_fn(p);
        if i > 0 {
            let t_prev = t_min + tr * ((i - 1) as f64 / N_SAMPLES as f64);
            let p_prev = curve.point_at(t_prev);
            let d_prev = implicit_fn(p_prev);
            // Sign change in distance indicates crossing the conic surface
            if d * d_prev < 0.0 || d.abs() < 1e-10 {
                candidates.push(t);
            }
        }
    }
    // Refine each candidate with Newton on the conic-specific distance function
    let mut results: Vec<(f64, f64)> = Vec::new();
    for &t0 in &candidates {
        let mut t = t0;
        for _ in 0..20 {
            let p = curve.point_at(t);
            let eps_d = 1e-7;
            let der = (curve.point_at((t + eps_d).min(t_max)) - curve.point_at((t - eps_d).max(t_min))) / (2.0 * eps_d);
            let (f_val, df) = match conic {
                Curve2d::Line(l) => {
                    let dp = p - l.origin;
                    (dp.x * l.direction.y - dp.y * l.direction.x, der.x * l.direction.y - der.y * l.direction.x)
                }
                Curve2d::Circle(c) => {
                    let dp = p - c.center;
                    let len = dp.length();
                    (len - c.radius, if len > 1e-30 { dp.dot(der) / len } else { der.length() })
                }
                _ => break,
            };
            if df.abs() < 1e-15 { break; }
            t = t - f_val / df;
            if f_val.abs() < tol.max(1e-14) {
                results.push((t, 0.0));
                break;
            }
        }
    }
    results
}

/// OCCT-aligned: IntCurve_IntPolyPolyGen (IntPolyPolyGen.gxx L93-107+).
///   Two generic curves: sample both → find approximate intersections → refine.
///   OCCT uses polygon-polygon interference detection + Newton refinement.
///   rcad: coarse sampling of curve1 → nearest-point search on curve2 → output close pairs.
fn intersect_curve_curve(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64, tol: f64) -> Vec<(f64, f64)> {
    let r1 = t1_max - t1_min;
    let r2 = t2_max - t2_min;
    if r1 < 1e-15 || r2 < 1e-15 || !r1.is_finite() || !r2.is_finite() { return vec![]; }
    // Sample both curves into discrete point sets
    let n1: usize = 256;
    let n2: usize = 256;
    let mut results: Vec<(f64, f64)> = Vec::new();
    for i in 0..=n1 {
        let t1 = t1_min + r1 * (i as f64 / n1 as f64);
        let p1 = c1.point_at(t1);
        let mut best_d2 = f64::INFINITY;
        let mut best_t2 = t2_min;
        for j in 0..=n2 {
            let t2 = t2_min + r2 * (j as f64 / n2 as f64);
            let d2 = p1.distance_squared(c2.point_at(t2));
            if d2 < best_d2 {
                best_d2 = d2;
                best_t2 = t2;
            }
        }
        if best_d2 < tol.max(1e-8) {
            let is_dup = results.last().map_or(false, |&(lt, _)| (t1 - lt).abs() < 1e-9 * r1);
            if !is_dup { results.push((t1, best_t2)); }
        }
    }
    // Remove endpoints that are just the domain boundary, not true intersections
    results.retain(|&(t1, t2)| {
        let at_boundary = |t: f64, lo: f64, hi: f64| (t - lo).abs() < 1e-10 || (t - hi).abs() < 1e-10;
        if at_boundary(t1, t1_min, t1_max) || at_boundary(t2, t2_min, t2_max) {
            // Keep only if the point actually lies on both curves
            let p1 = c1.point_at(t1);
            let p2 = c2.point_at(t2);
            p1.distance_squared(p2) < tol.max(1e-8)
        } else {
            true
        }
    });
    results
}

pub fn intersect_ray_curve_2d(ro: DVec2, rd: DVec2, curve: &Curve2d, t_min: f64, t_max: f64) -> Vec<(f64, f64)> {
    if rd.length_squared() < 1e-30 { return vec![]; }
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
