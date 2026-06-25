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

fn intersect_circle_ellipse(circle: &Curve2d, ell: &Curve2d, tc_min: f64, tc_max: f64, te_min: f64, te_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Circle(c) = circle else { return vec![] }; let Curve2d::Ellipse(e) = ell else { return vec![] };
    let minor = DVec2::new(-e.major_dir.y, e.major_dir.x);
    let tr = te_max - te_min; if tr < 1e-15 { return vec![]; }
    let mut cand = Vec::new();
    for i in 0..=128 {
        let t = te_min + tr*(i as f64)/128.0;
        let p = e.center + e.major_dir*(e.major_radius*t.cos()) + minor*(e.minor_radius*t.sin());
        let v = p - c.center; let f = v.length_squared() - c.radius*c.radius;
        if i > 0 { let t2 = te_min + tr*((i-1)as f64)/128.0;
            let p2 = e.center + e.major_dir*(e.major_radius*t2.cos()) + minor*(e.minor_radius*t2.sin());
            let f2 = (p2 - c.center).length_squared() - c.radius*c.radius;
            if f == 0.0 || f*f2 < 0.0 { cand.push(t); }
        }
        if f.abs() < 1e-8 { let dup = cand.last().map_or(false, |&lt| (t-lt).abs() < 1e-9*tr); if !dup { cand.push(t); } }
    }
    let mut res = Vec::new();
    for &t0 in &cand {
        let mut t = t0.clamp(te_min, te_max);
        for _ in 0..20 {
            let p = e.center + e.major_dir*(e.major_radius*t.cos()) + minor*(e.minor_radius*t.sin());
            let dp = p - c.center; let f = dp.length_squared() - c.radius*c.radius;
            if f.abs() < 1e-14 {
                let mut tc = (p.y - c.center.y).atan2(p.x - c.center.x); if tc < 0.0 { tc += std::f64::consts::TAU; }
                if tc >= tc_min-1e-10 && tc <= tc_max+1e-10 && t >= te_min-1e-10 && t <= te_max+1e-10 {
                    let dup = res.last().map_or(false, |&(lt,_): &(f64,f64)| (t-lt).abs() < 1e-9*tr); if !dup { res.push((tc, t)); }
                }
                break;
            }
            let der = e.major_dir*(-e.major_radius*t.sin()) + minor*(e.minor_radius*t.cos());
            let df = 2.0*dp.dot(der); if df.abs() < 1e-15 { break; }
            t = (t - f/df).clamp(te_min, te_max); if f.abs() < 1e-14 { break; }
        }
    }
    dedup(res)
}

fn intersect_ellipse_ellipse(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64) -> Vec<(f64, f64)> {
    let Curve2d::Ellipse(e1) = c1 else { return vec![] }; let Curve2d::Ellipse(e2) = c2 else { return vec![] };
    let m1 = DVec2::new(-e1.major_dir.y, e1.major_dir.x); let m2 = DVec2::new(-e2.major_dir.y, e2.major_dir.x);
    let r1 = t1_max - t1_min; if r1 < 1e-15 { return vec![]; }
    let mut res = Vec::new();
    for i in 0..=128 {
        let t1 = t1_min + r1*(i as f64)/128.0;
        let p1 = e1.center + e1.major_dir*(e1.major_radius*t1.cos()) + m1*(e1.minor_radius*t1.sin());
        let r2 = t2_max - t2_min; if r2 < 1e-15 { continue; }
        let mut bd2 = f64::INFINITY; let mut bt2 = t2_min;
        for j in 0..=256 { let t2 = t2_min + r2*(j as f64)/256.0;
            let p2 = e2.center + e2.major_dir*(e2.major_radius*t2.cos()) + m2*(e2.minor_radius*t2.sin());
            let d2 = p1.distance_squared(p2); if d2 < bd2 { bd2 = d2; bt2 = t2; }
        }
        if bd2 < 1e-8 && res.last().map_or(true, |&(lt,_):&(f64,f64)| (t1-lt).abs() > 1e-9*r1) { res.push((t1, bt2)); }
    }
    res
}

fn intersect_conic_curve(conic: &Curve2d, _typ: G2dCurveType, curve: &Curve2d, _tc_min: f64, _tc_max: f64, t_min: f64, t_max: f64, _tol: f64) -> Vec<(f64, f64)> {
    let tr = t_max - t_min; if tr < 1e-15 || !tr.is_finite() { return vec![]; }
    let mut best: Option<(f64,f64,f64)> = None;
    for i in 0..=256 {
        let t = t_min + tr*(i as f64)/256.0; let p = curve.point_at(t);
        let d = if let Curve2d::Line(l) = conic { let dp = p - l.origin; (dp.x*l.direction.y - dp.y*l.direction.x).abs()/l.direction.length().max(1e-30) }
        else if let Curve2d::Circle(c) = conic { (p.distance(c.center) - c.radius).abs() }
        else if let Curve2d::Ellipse(e) = conic { let u = e.major_dir; let v = DVec2::new(-u.y,u.x); let du = (p-e.center).dot(u)/e.major_radius; let dv = (p-e.center).dot(v)/e.minor_radius; (du*du + dv*dv).sqrt() - 1.0 }
        else { continue; };
        if best.map_or(true, |(bd,_,_)| d < bd) { best = Some((d, t, 0.0)); }
    }
    if let Some((dist, t0, _)) = best {
        if dist > 1e-4 { return vec![]; }
        let mut t = t0;
        for _ in 0..20 {
            let p = curve.point_at(t);
            let der = (curve.point_at((t+1e-7).clamp(t_min,t_max)) - curve.point_at((t-1e-7).clamp(t_min,t_max)))/(2.0*1e-7);
            let (f, df) = if let Curve2d::Line(l) = conic { let dp = p - l.origin; (dp.x*l.direction.y - dp.y*l.direction.x, der.x*l.direction.y - der.y*l.direction.x) }
            else if let Curve2d::Circle(c) = conic { let dp = p - c.center; let len = dp.length(); (len - c.radius, if len > 1e-30 { dp.dot(der)/len } else { der.length() }) }
            else { break; };
            if df.abs() < 1e-15 { break; }
            t = (t - f/df).clamp(t_min, t_max); if f.abs() < 1e-14 { return vec![(t, 0.0)]; }
        }
    }
    vec![]
}

fn intersect_curve_curve(c1: &Curve2d, c2: &Curve2d, t1_min: f64, t1_max: f64, t2_min: f64, t2_max: f64, _tol: f64) -> Vec<(f64, f64)> {
    let r1 = t1_max - t1_min; let r2 = t2_max - t2_min;
    if r1 < 1e-15 || r2 < 1e-15 || !r1.is_finite() || !r2.is_finite() { return vec![]; }
    let mut res = Vec::new();
    for i in 0..=64 {
        let t1 = t1_min + r1*(i as f64)/64.0; let p1 = c1.point_at(t1);
        let mut bd2 = f64::INFINITY; let mut bt2 = t2_min;
        for j in 0..=128 { let t2 = t2_min + r2*(j as f64)/128.0; let d2 = p1.distance_squared(c2.point_at(t2)); if d2 < bd2 { bd2 = d2; bt2 = t2; } }
        if bd2 < 1e-8 && res.last().map_or(true, |&(lt,_):&(f64,f64)| (t1-lt).abs() > 1e-9*r1) { res.push((t1, bt2)); }
    }
    res
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
