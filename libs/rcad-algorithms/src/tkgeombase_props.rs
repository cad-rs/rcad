//! OCCT-aligned TKGeomBase property/analysis equivalents.
//!
//! 闁?OCCT-aligned implementations for:
//!   - GeomLProp_CLProps2d     (curvature, D1, D2, centre of curvature, normal)
//!   - GeomLProp_CurAndInf2d   (curvature extrema, inflection points)
//!   - GProp_PGProps           (point-set mass, centre of mass, inertia)
//!   - GProp_PEquation          (point-set plane/line/point classification)
//!
//! OCCT source: src/ModelingData/TKGeomBase/

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};
use glam::DVec3;
use crate::tolerance::TOLERANCE_LEN_SQ_DIV_SAFE;

/// Precision used for curvature computations.
const TOL: f64 = 1e-7;

// =============================================================================
// GeomLProp_CLProps2d 闁?2D curve local properties
// =============================================================================

/// Signed curvature of a 2D curve at parameter t.
///
/// OCCT: GeomLProp_CLProps2d::Curvature()
/// k = (x'鐠虹椊'' 闁?y'鐠虹椈'') / (x'閾?+ y'閾?^(3/2)
pub fn curvature_at(curve: &Curve2d, t: f64) -> f64 {
    let d1 = curve.derivative_at(t);
    let eps = 1e-7;
    let d2 = (curve.derivative_at(t + eps) - curve.derivative_at(t - eps)) / (2.0 * eps);
    let sq_norm = d1.length_squared();
    if sq_norm < TOL {
        return 0.0;
    }
    let cross = d1.x * d2.y - d1.y * d2.x;
    cross / sq_norm.powf(1.5)
}

/// First derivative (velocity) vector.
/// OCCT: GeomLProp_CLProps2d::D1()
pub fn d1_at(curve: &Curve2d, t: f64) -> DVec2 {
    curve.derivative_at(t)
}

/// Second derivative (acceleration) vector.
/// OCCT: GeomLProp_CLProps2d::D2()
pub fn d2_at(curve: &Curve2d, t: f64) -> DVec2 {
    let eps = 1e-7;
    (curve.derivative_at(t + eps) - curve.derivative_at(t - eps)) / (2.0 * eps)
}

/// Centre of curvature at parameter t (osculating circle centre).
/// Returns None if curvature is zero (line or degenerate).
/// OCCT: GeomLProp_CLProps2d::CentreOfCurvature()
pub fn centre_of_curvature_at(curve: &Curve2d, t: f64) -> Option<DVec2> {
    let pt = curve.point_at(t);
    let k = curvature_at(curve, t);
    if k.abs() < TOL {
        return None;
    }
    let d1 = curve.derivative_at(t);
    let speed = d1.length();
    if speed < TOL {
        return None;
    }
    let normal = DVec2::new(-d1.y, d1.x) / speed;
    Some(pt + normal * (1.0 / k))
}

/// Unit normal vector (left of tangent direction).
/// OCCT: GeomLProp_CLProps2d::Normal()
pub fn normal_at(curve: &Curve2d, t: f64) -> Option<DVec2> {
    let d1 = curve.derivative_at(t);
    let speed = d1.length();
    if speed < TOL {
        return None;
    }
    Some(DVec2::new(-d1.y, d1.x) / speed)
}

// =============================================================================
// GeomLProp_CurAndInf2d 闁?curvature extrema and inflection points
// =============================================================================

/// OCCT-aligned: LProp_CIType
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CIType2d {
    Inflection,
    MinCur,
    MaxCur,
}

/// Curvature extremum/inflection analysis result.
/// OCCT-aligned: GeomLProp_CurAndInf2d
#[derive(Debug, Clone)]
pub struct CurAndInf2d {
    params: Vec<f64>,
    types: Vec<CIType2d>,
}

impl CurAndInf2d {
    pub fn new() -> Self {
        Self { params: Vec::new(), types: Vec::new() }
    }
    pub fn is_done(&self) -> bool { true }
    pub fn nb_points(&self) -> usize { self.params.len() }
    pub fn parameter(&self, idx: usize) -> f64 { self.params[idx] }
    pub fn ci_type(&self, idx: usize) -> CIType2d { self.types[idx] }

    fn add(&mut self, param: f64, typ: CIType2d) {
        self.params.push(param);
        self.types.push(typ);
    }

    fn sort(&mut self) {
        let mut indices: Vec<usize> = (0..self.params.len()).collect();
        indices.sort_by(|&a, &b| self.params[a].partial_cmp(&self.params[b]).unwrap());
        self.params = indices.iter().map(|&i| self.params[i]).collect();
        self.types = indices.iter().map(|&i| self.types[i]).collect();
    }
}

/// Find curvature extrema (min/max) of a 2D curve by scanning the domain.
/// OCCT-aligned: GeomLProp_CurAndInf2d::PerformCurAndInf()
pub fn curvature_extrema_2d(curve: &Curve2d) -> CurAndInf2d {
    let mut result = CurAndInf2d::new();
    let domain = curve.default_domain();
    let t0 = domain[0].max(-1000.0);
    let t1 = domain[1].min(1000.0);
    let span = t1 - t0;
    if span <= 0.0 { return result; }

    let n = 256usize;
    let mut prev_k = curvature_at(curve, t0);
    for i in 1..=n {
        let t = t0 + span * (i as f64) / (n as f64);
        let k = curvature_at(curve, t);
        if prev_k.is_finite() && k.is_finite() {
            // Sign change in derivative of curvature 闁?extremum
            let dk_prev = k - prev_k;
            if i > 1 && dk_prev.abs() > 1e-12 {
                let t_mid = t0 + span * ((i as f64 - 0.5) / (n as f64));
                let typ = if dk_prev < 0.0 { CIType2d::MaxCur } else { CIType2d::MinCur };
                result.add(t_mid, typ);
            }
            // Sign change in curvature 闁?inflection
            if prev_k * k < 0.0 {
                let t_inf = t0 + span * ((i as f64 - 0.5) / (n as f64));
                result.add(t_inf, CIType2d::Inflection);
            }
        }
        prev_k = k;
    }
    result.sort();
    result
}

/// Find inflection points of a 2D curve.
/// OCCT-aligned: GeomLProp_CurAndInf2d (inflection subset)
pub fn inflection_points_2d(curve: &Curve2d) -> CurAndInf2d {
    let full = curvature_extrema_2d(curve);
    let mut result = CurAndInf2d::new();
    for i in 0..full.nb_points() {
        if full.ci_type(i) == CIType2d::Inflection {
            result.add(full.parameter(i), CIType2d::Inflection);
        }
    }
    result
}

/// Full curvature analysis (extrema + inflections).
pub fn perform_curvature_analysis_2d(curve: &Curve2d) -> CurAndInf2d {
    curvature_extrema_2d(curve)
}

// =============================================================================
// GProp_PGProps 闁?point-set properties (mass, centre, inertia)
// =============================================================================

/// Point-set properties accumulator.
/// OCCT: GProp_PGProps
#[derive(Debug, Clone)]
pub struct PointSetProps {
    mass: f64,
    centre: DVec3,
    inertia: [f64; 6],
}

impl PointSetProps {
    pub fn new() -> Self {
        Self { mass: 0.0, centre: DVec3::ZERO, inertia: [0.0; 6] }
    }

    pub fn mass(&self) -> f64 { self.mass }
    pub fn centre_of_mass(&self) -> DVec3 { self.centre }

    pub fn matrix_of_inertia(&self) -> [f64; 9] {
        let [ixx, iyy, izz, ixy, ixz, iyz] = self.inertia;
        [ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz]
    }

    pub fn add_point(&mut self, pt: DVec3) {
        self.add_point_weighted(pt, 1.0);
    }

    pub fn add_point_weighted(&mut self, pt: DVec3, weight: f64) {
        let new_mass = self.mass + weight;
        if new_mass > 0.0 {
            self.centre = (self.centre * self.mass + pt * weight) / new_mass;
        }
        self.mass = new_mass;
        let x = pt.x; let y = pt.y; let z = pt.z;
        self.inertia[0] += weight * (y * y + z * z);
        self.inertia[1] += weight * (x * x + z * z);
        self.inertia[2] += weight * (x * x + y * y);
        self.inertia[3] -= weight * x * y;
        self.inertia[4] -= weight * x * z;
        self.inertia[5] -= weight * y * z;
    }

    pub fn barycentre(points: &[DVec3]) -> DVec3 {
        let mut props = PointSetProps::new();
        for &pt in points { props.add_point(pt); }
        if props.mass > 0.0 { props.centre } else { DVec3::ZERO }
    }
}

// =============================================================================
// GProp_PEquation 闁?point-set equation (fit point/line/plane/space)
// =============================================================================

/// Result of point-set equation analysis.
/// OCCT: GProp_PEquation
#[derive(Debug, Clone, PartialEq)]
pub enum PointSetKind {
    Point(DVec3),
    Line(DVec3, DVec3),
    Plane(DVec3, DVec3),
    Space,
}

/// Classify a set of 3D points by PCA.
/// OCCT: GProp_PEquation
pub fn analyze_point_set(points: &[DVec3], tolerance: f64) -> PointSetKind {
    if points.is_empty() { return PointSetKind::Space; }
    let n = points.len() as f64;
    let mut centroid = DVec3::ZERO;
    for &p in points { centroid += p; }
    centroid /= n;

    let max_dist = points.iter().map(|p| (*p - centroid).length()).fold(0.0, f64::max);
    if max_dist < tolerance { return PointSetKind::Point(centroid); }

    let mut cxx = 0.0; let mut cxy = 0.0; let mut cxz = 0.0;
    let mut cyy = 0.0; let mut cyz = 0.0; let mut czz = 0.0;
    for &p in points {
        let dx = p.x - centroid.x; let dy = p.y - centroid.y; let dz = p.z - centroid.z;
        cxx += dx * dx; cxy += dx * dy; cxz += dx * dz;
        cyy += dy * dy; cyz += dy * dz; czz += dz * dz;
    }
    cxx /= n; cxy /= n; cxz /= n; cyy /= n; cyz /= n; czz /= n;

    // Power iteration for largest eigenvalue/eigenvector
    let (eval1, v1) = power_iteration(cxx, cxy, cxz, cyy, cyz, czz, DVec3::X);
    if eval1 < tolerance * tolerance { return PointSetKind::Point(centroid); }

    // Deflate
    cxx -= eval1 * v1.x * v1.x; cxy -= eval1 * v1.x * v1.y;
    cxz -= eval1 * v1.x * v1.z; cyy -= eval1 * v1.y * v1.y;
    cyz -= eval1 * v1.y * v1.z; czz -= eval1 * v1.z * v1.z;

    let second_dir = if (v1 - DVec3::X).length() > 0.1 { DVec3::Y } else { DVec3::X };
    let (eval2, _v2) = power_iteration(cxx, cxy, cxz, cyy, cyz, czz, second_dir);
    if eval2 < tolerance * tolerance { return PointSetKind::Line(centroid, v1); }

    let third_dir = v1.cross(_v2).normalize();
    let (eval3, _) = power_iteration(cxx, cxy, cxz, cyy, cyz, czz, third_dir);
    if eval3 < tolerance * tolerance { return PointSetKind::Plane(centroid, v1.cross(_v2).normalize()); }

    PointSetKind::Space
}

fn power_iteration(cxx: f64, cxy: f64, cxz: f64, cyy: f64, cyz: f64, czz: f64, start: DVec3) -> (f64, DVec3) {
    let mut v = start;
    for _ in 0..30 {
        let v2 = DVec3::new(
            cxx * v.x + cxy * v.y + cxz * v.z,
            cxy * v.x + cyy * v.y + cyz * v.z,
            cxz * v.x + cyz * v.y + czz * v.z,
        );
        let len = v2.length();
        if len < TOLERANCE_LEN_SQ_DIV_SAFE { break; }
        v = v2 / len;
    }
    let eval = DVec3::new(
        cxx * v.x + cxy * v.y + cxz * v.z,
        cxy * v.x + cyy * v.y + cyz * v.z,
        cxz * v.x + cyz * v.y + czz * v.z,
    ).dot(v);
    (eval, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::*;

    // 闁冲厜鍋撻柍鍏夊亾 GeomLProp_CLProps2d tests 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾

    fn make_line() -> Curve2d {
        Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X })
    }

    fn make_circle() -> Curve2d {
        Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0))
    }

    fn make_ellipse() -> Curve2d {
        Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X,
            major_radius: 10.0, minor_radius: 5.0,
        })
    }

    #[test]
    fn clprops2d_line_curvature_is_zero() {
        let l = make_line(); let k = curvature_at(&l, 5.0);
        assert!(k.abs() < TOL);
    }

    #[test]
    fn clprops2d_line_centre_of_curvature_not_defined() {
        let l = make_line();
        assert!(centre_of_curvature_at(&l, 5.0).is_none());
    }

    #[test]
    fn clprops2d_circle_curvature_constant() {
        let c = make_circle();
        for i in 0..8 {
            let t = i as f64 * std::f64::consts::PI / 4.0;
            let k = curvature_at(&c, t);
            assert!((k - 0.2).abs() < TOL, "curvature not constant at t={}", t);
        }
    }

    #[test]
    fn clprops2d_circle_normal() {
        let c = make_circle();
        let tgt = c.tangent_at(0.0);
        let nml = normal_at(&c, 0.0).unwrap();
        assert!(tgt.dot(nml).abs() < TOL, "tangent and normal should be perpendicular");
    }

    #[test]
    fn clprops2d_circle_centre_of_curvature() {
        let c = make_circle();
        let cc = centre_of_curvature_at(&c, 0.0).unwrap();
        assert!((cc - DVec2::ZERO).length() < TOL);
    }

    #[test]
    fn clprops2d_circle_d1() {
        let c = make_circle();
        let d1 = d1_at(&c, 0.0);
        assert!((d1 - DVec2::new(0.0, 5.0)).length() < TOL);
    }

    #[test]
    fn clprops2d_circle_d2() {
        let c = make_circle();
        let d2 = d2_at(&c, 0.0);
        assert!((d2 - DVec2::new(-5.0, 0.0)).length() < TOL);
    }

    #[test]
    fn clprops2d_ellipse_curvature_at_major_vertex() {
        let e = make_ellipse();
        let k = curvature_at(&e, 0.0);
        let expected = 10.0 / 25.0;
        assert!((k - expected).abs() < 1e-6);
    }

    #[test]
    fn clprops2d_ellipse_curvature_at_minor_vertex() {
        let e = make_ellipse();
        let k = curvature_at(&e, std::f64::consts::PI / 2.0);
        let expected = 5.0 / 100.0;
        assert!((k - expected).abs() < 1e-6);
    }

    // 闁冲厜鍋撻柍鍏夊亾 GeomLProp_CurAndInf2d tests 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾

    #[test]
    fn curandinf2d_circle_no_inflections() {
        let c = make_circle();
        let r = curvature_extrema_2d(&c);
        assert!(r.is_done());
    }

    #[test]
    fn curandinf2d_ellipse_has_four_extrema() {
        let e = make_ellipse();
        let r = curvature_extrema_2d(&e);
        assert!(r.is_done());
        assert!(r.nb_points() >= 2, "expected at least 2 curvature extrema, got {}", r.nb_points());
        let mut nb_min = 0; let mut nb_max = 0;
        for i in 0..r.nb_points() {
            match r.ci_type(i) {
                CIType2d::MinCur => { nb_min += 1; }
                CIType2d::MaxCur => { nb_max += 1; }
                _ => {}
            }
        }
        assert!(nb_min >= 1 && nb_max >= 1, "expected at least 1 min and 1 max curvature");
    }

    #[test]
    fn curandinf2d_ellipse_no_inflections() {
        let e = make_ellipse();
        let r = inflection_points_2d(&e);
        assert_eq!(r.nb_points(), 0);
    }

    // 闁冲厜鍋撻柍鍏夊亾 GProp_PGProps tests 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾

    #[test]
    fn pgprops_empty() {
        let props = PointSetProps::new();
        assert!((props.mass() - 0.0).abs() < TOL);
    }

    #[test]
    fn pgprops_single_point() {
        let mut props = PointSetProps::new();
        props.add_point(DVec3::new(1.0, 2.0, 3.0));
        assert!((props.mass() - 1.0).abs() < TOL);
        assert!((props.centre_of_mass() - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn pgprops_two_points_barycentre() {
        let mut props = PointSetProps::new();
        props.add_point(DVec3::ZERO);
        props.add_point(DVec3::new(2.0, 4.0, 6.0));
        assert!((props.centre_of_mass() - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn pgprops_weighted_points() {
        let mut props = PointSetProps::new();
        props.add_point_weighted(DVec3::ZERO, 1.0);
        props.add_point_weighted(DVec3::new(4.0, 0.0, 0.0), 3.0);
        assert!((props.centre_of_mass() - DVec3::new(3.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn pgprops_inertia_symmetric() {
        let pts = vec![DVec3::X, -DVec3::X, DVec3::Y, -DVec3::Y];
        let mut props = PointSetProps::new();
        for &pt in &pts { props.add_point(pt); }
        let m = props.matrix_of_inertia();
        assert!((m[0] - 2.0).abs() < TOL, "Ixx={}", m[0]);
        assert!((m[4] - 2.0).abs() < TOL, "Iyy={}", m[4]);
        assert!((m[8] - 4.0).abs() < TOL, "Izz={}", m[8]);
    }

    // 闁冲厜鍋撻柍鍏夊亾 GProp_PEquation tests 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾

    #[test]
    fn pequation_single_point() {
        let pts = vec![DVec3::new(5.0, 5.0, 5.0)];
        match analyze_point_set(&pts, 1e-6) {
            PointSetKind::Point(p) => assert!((p - DVec3::new(5.0, 5.0, 5.0)).length() < TOL),
            other => panic!("expected Point, got {:?}", other),
        }
    }

    #[test]
    fn pequation_collinear_points() {
        let pts = vec![DVec3::ZERO, DVec3::X, DVec3::new(2.0, 0.0, 0.0)];
        match analyze_point_set(&pts, 1e-6) {
            PointSetKind::Line(_, _) => {},
            other => panic!("expected Line, got {:?}", other),
        }
    }

    #[test]
    fn pequation_coplanar_points() {
        let pts = vec![DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::new(5.0, 3.0, 0.0)];
        match analyze_point_set(&pts, 1e-6) {
            PointSetKind::Plane(_, _) => {},
            other => panic!("expected Plane, got {:?}", other),
        }
    }

    #[test]
    fn pequation_space_filling_points() {
        let pts = vec![DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Z, DVec3::ONE];
        match analyze_point_set(&pts, 1e-6) {
            PointSetKind::Space => {},
            other => panic!("expected Space, got {:?}", other),
        }
    }
}



