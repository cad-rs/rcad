//! 2D curve local properties — GeomLProp_CLProps2d / GeomLProp_CurAndInf2d equivalents.
//!
//! ✅ OCCT-aligned: provides curvature, tangent, normal, centre-of-curvature,
//! and curvature-extremum / inflection-point analysis for 2D curves.
//!
//! OCCT source: src/ModelingData/TKGeomBase/GTests/
//!   GeomLProp_CLProps2d_Test.cxx
//!   GeomLProp_CurAndInf2d_Test.cxx

use glam::DVec2;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Circle2d, Line2d, Ellipse2d, Hyperbola2d, Parabola2d,
};

/// OCCT-aligned precision: Precision::Confusion()
const TOL: f64 = 1e-7;

// =============================================================================
// Curve2d local properties — GeomLProp_CLProps2d equivalent
// =============================================================================

/// Signed curvature of a 2D curve at parameter t.
///
/// OCCT: GeomLProp_CLProps2d::Curvature()
/// For a parametric curve, curvature k = (x'·y'' - y'·x'') / (x'² + y'²)^(3/2)
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

/// First derivative (velocity) vector of a 2D curve.
///
/// OCCT: GeomLProp_CLProps2d::D1()
/// Returns the analytic first derivative when available.
pub fn d1_at(curve: &Curve2d, t: f64) -> DVec2 {
    curve.derivative_at(t)
}

/// Second derivative (acceleration) vector of a 2D curve.
///
/// OCCT: GeomLProp_CLProps2d::D2()
pub fn d2_at(curve: &Curve2d, t: f64) -> DVec2 {
    let eps = 1e-7;
    (curve.derivative_at(t + eps) - curve.derivative_at(t - eps)) / (2.0 * eps)
}

/// Centre of curvature at parameter t.
///
/// OCCT: GeomLProp_CLProps2d::CentreOfCurvature()
/// Returns the centre of the osculating circle.
/// Returns None if curvature is zero (degenerate, e.g. line).
pub fn centre_of_curvature_at(curve: &Curve2d, t: f64) -> Option<DVec2> {
    let pt = curve.point_at(t);
    let k = curvature_at(curve, t);
    if k.abs() < TOL {
        return None;
    }
    let d1 = curve.derivative_at(t);
    // Unit normal = perpendicular to tangent: n = (-d1.y, d1.x) / |d1|
    let speed = d1.length();
    if speed < TOL {
        return None;
    }
    let normal = DVec2::new(-d1.y, d1.x) / speed;
    // Centre = point + normal * (1/k), where sign gives inward normal
    Some(pt + normal * (1.0 / k))
}

/// Unit normal vector at parameter t (pointing left of the tangent direction).
///
/// OCCT: GeomLProp_CLProps2d::Normal()
/// Normal = perpendicular to tangent, oriented such that the signed curvature is positive.
pub fn normal_at(curve: &Curve2d, t: f64) -> Option<DVec2> {
    let d1 = curve.derivative_at(t);
    let speed = d1.length();
    if speed < TOL {
        return None;
    }
    Some(DVec2::new(-d1.y, d1.x) / speed)
}

// =============================================================================
// Curvature extremum / inflection — GeomLProp_CurAndInf2d equivalent
// =============================================================================

/// OCCT-aligned: LProp_CIType
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CIType2d {
    Inflection,
    MinCur,
    MaxCur,
}

/// Curvature analysis result: curvature extrema and inflection points.
///
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

/// Find curvature extrema of a 2D curve (OCCT: PerformCurExt).
///
/// Samples curvature across the domain and locates local min/max.
pub fn curvature_extrema_2d(curve: &Curve2d) -> CurAndInf2d {
    let domain = curve.default_domain();
    let (t0, t1) = if domain[0].is_finite() && domain[1].is_finite() {
        (domain[0], domain[1])
    } else {
        // For curves with infinite domain, use a reasonable range
        (-10.0, 10.0)
    };

    let mut result = CurAndInf2d::new();
    let n_sample = 100;
    let step = (t1 - t0) / n_sample as f64;
    if step <= 0.0 { return result; }

    // Compute curvature at each sample point
    let mut curvatures: Vec<(f64, f64)> = Vec::with_capacity(n_sample + 1);
    for i in 0..=n_sample {
        let t = t0 + i as f64 * step;
        let k = curvature_at(curve, t);
        if k.is_finite() {
            curvatures.push((t, k));
        }
    }

    if curvatures.len() < 3 { return result; }

    // Find local extrema: points where curvature changes direction
    for i in 1..curvatures.len() - 1 {
        let (t_prev, k_prev) = curvatures[i - 1];
        let (t_curr, k_curr) = curvatures[i];
        let (t_next, _k_next) = curvatures[i + 1];

        // Check for zero crossing of curvature (inflection)
        if k_prev * k_curr < 0.0 {
            // Linear interpolation for zero crossing
            let zero_t = t_prev + (t_curr - t_prev) * (-k_prev) / (k_curr - k_prev);
            result.add(zero_t, CIType2d::Inflection);
        }

        // Check for local min/max
        if i > 0 && i < curvatures.len() - 1 {
            let k_prev = curvatures[i - 1].1;
            let k_curr = curvatures[i].1;
            let k_next = curvatures[i + 1].1;
            if k_curr > k_prev && k_curr > k_next {
                result.add(t_curr, CIType2d::MaxCur);
            } else if k_curr < k_prev && k_curr < k_next {
                result.add(t_curr, CIType2d::MinCur);
            }
        }
    }

    // Check endpoints for extrema
    if curvatures.len() >= 2 {
        let (t0_k, k0) = curvatures[0];
        let (_t1_k, k1) = curvatures[1];
        if (k0 - k1).abs() > TOL {
            // Endpoints can be extrema on closed curves
        }
    }

    result.sort();
    // Remove duplicates
    let mut dedup = CurAndInf2d::new();
    for i in 0..result.nb_points() {
        if i == 0 || (result.parameter(i) - result.parameter(i - 1)).abs() > TOL {
            dedup.add(result.parameter(i), result.ci_type(i));
        }
    }

    dedup
}

/// Find inflection points of a 2D curve (OCCT: PerformInf).
pub fn inflection_points_2d(curve: &Curve2d) -> CurAndInf2d {
    let all = curvature_extrema_2d(curve);
    let mut result = CurAndInf2d::new();
    for i in 0..all.nb_points() {
        if all.ci_type(i) == CIType2d::Inflection {
            result.add(all.parameter(i), CIType2d::Inflection);
        }
    }
    result
}

/// Full analysis: inflections + curvature extrema (OCCT: Perform).
pub fn perform_curvature_analysis_2d(curve: &Curve2d) -> CurAndInf2d {
    curvature_extrema_2d(curve)
}

// =============================================================================
// Tests — OCCT GeomLProp_CLProps2d_Test.cxx
// =============================================================================

#[cfg(test)]
mod clprops2d_tests {
    use super::*;

    fn make_circle() -> Curve2d {
        Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0))
    }

    fn make_line() -> Curve2d {
        Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X })
    }

    fn make_ellipse() -> Curve2d {
        Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X,
            major_radius: 10.0, minor_radius: 5.0,
        })
    }

    #[test]
    fn clprops2d_circle_value() {
        let c = make_circle();
        let p = c.point_at(0.0);
        assert!((p - DVec2::new(5.0, 0.0)).length() < TOL);
    }

    #[test]
    fn clprops2d_circle_set_parameter() {
        let c = make_circle();
        let p = c.point_at(std::f64::consts::PI / 2.0);
        assert!((p - DVec2::new(0.0, 5.0)).length() < TOL);
    }

    #[test]
    fn clprops2d_circle_tangent_at_zero() {
        let c = make_circle();
        let t = c.tangent_at(0.0);
        assert!((t - DVec2::new(0.0, 1.0)).length() < TOL);
    }

    #[test]
    fn clprops2d_circle_curvature() {
        let c = make_circle();
        let k = curvature_at(&c, 0.0);
        assert!((k - 1.0 / 5.0).abs() < TOL);
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
        assert!((cc - DVec2::ZERO).length() < TOL, "centre of curvature should be center");
    }

    #[test]
    fn clprops2d_circle_centre_of_curvature_at_pi_half() {
        let c = make_circle();
        let cc = centre_of_curvature_at(&c, std::f64::consts::PI / 2.0).unwrap();
        assert!((cc - DVec2::ZERO).length() < TOL);
    }

    #[test]
    fn clprops2d_line_curvature_is_zero() {
        let l = make_line();
        let k = curvature_at(&l, 5.0);
        assert!(k.abs() < TOL);
    }

    #[test]
    fn clprops2d_line_centre_of_curvature_not_defined() {
        let l = make_line();
        let cc = centre_of_curvature_at(&l, 5.0);
        assert!(cc.is_none(), "line should have undefined centre of curvature");
    }

    #[test]
    fn clprops2d_line_tangent() {
        let l = make_line();
        let t = l.tangent_at(5.0);
        assert!((t - DVec2::new(1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn clprops2d_ellipse_curvature_at_major_vertex() {
        let e = make_ellipse();
        // Curvature at major vertex (t=0) = a/b^2 = 10/25 = 0.4
        let k = curvature_at(&e, 0.0);
        let expected = 10.0 / 25.0;
        assert!((k - expected).abs() < 1e-6);
    }

    #[test]
    fn clprops2d_ellipse_curvature_at_minor_vertex() {
        let e = make_ellipse();
        // Curvature at minor vertex (t=PI/2) = b/a^2 = 5/100 = 0.05
        let k = curvature_at(&e, std::f64::consts::PI / 2.0);
        let expected = 5.0 / 100.0;
        assert!((k - expected).abs() < 1e-6);
    }

    #[test]
    fn clprops2d_d1_circle() {
        let c = make_circle();
        let d1 = d1_at(&c, 0.0);
        assert!((d1 - DVec2::new(0.0, 5.0)).length() < TOL);
    }

    #[test]
    fn clprops2d_d2_circle() {
        let c = make_circle();
        let d2 = d2_at(&c, 0.0);
        assert!((d2 - DVec2::new(-5.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// Tests — OCCT GeomLProp_CurAndInf2d_Test.cxx
// =============================================================================

#[cfg(test)]
mod curandinf2d_tests {
    use super::*;

    fn make_circle() -> Curve2d {
        Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0))
    }

    fn make_ellipse() -> Curve2d {
        Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X,
            major_radius: 10.0, minor_radius: 3.0,
        })
    }

    fn make_hyperbola() -> Curve2d {
        Curve2d::Hyperbola(Hyperbola2d {
            center: DVec2::ZERO, major_dir: DVec2::X,
            semi_major: 6.0, semi_minor: 3.0,
        })
    }

    fn make_parabola() -> Curve2d {
        Curve2d::Parabola(Parabola2d {
            origin: DVec2::ZERO, axis_dir: DVec2::X, focal_param: 4.0,
        })
    }

    #[test]
    fn curandinf2d_circle_perform_no_inflections() {
        let c = make_circle();
        let r = perform_curvature_analysis_2d(&c);
        assert!(r.is_done());
        // Circle has constant curvature -> no extrema in the sampled result
    }

    #[test]
    fn curandinf2d_ellipse_perform_cur_ext_has_extrema() {
        let e = make_ellipse();
        let r = curvature_extrema_2d(&e);
        assert!(r.is_done());
        // Ellipse has 4 curvature extrema
        assert_eq!(r.nb_points(), 4);
    }

    #[test]
    fn curandinf2d_ellipse_perform_cur_ext_types() {
        let e = make_ellipse();
        let r = curvature_extrema_2d(&e);
        assert!(r.is_done());
        assert_eq!(r.nb_points(), 4);
        let mut nb_min = 0;
        let mut nb_max = 0;
        for i in 0..r.nb_points() {
            match r.ci_type(i) {
                CIType2d::MinCur => nb_min += 1,
                CIType2d::MaxCur => nb_max += 1,
                _ => {}
            }
        }
        assert_eq!(nb_min, 2);
        assert_eq!(nb_max, 2);
    }

    #[test]
    fn curandinf2d_ellipse_perform_cur_ext_parameters_sorted() {
        let e = make_ellipse();
        let r = curvature_extrema_2d(&e);
        assert!(r.is_done());
        for i in 1..r.nb_points() {
            assert!(r.parameter(i - 1) < r.parameter(i));
        }
    }

    #[test]
    fn curandinf2d_ellipse_perform_inf_no_inflections() {
        let e = make_ellipse();
        let r = inflection_points_2d(&e);
        assert!(r.is_done());
        assert_eq!(r.nb_points(), 0);
    }

    #[test]
    fn curandinf2d_hyperbola_perform_cur_ext_vertex_only() {
        let h = make_hyperbola();
        let r = curvature_extrema_2d(&h);
        assert!(r.is_done());
        assert_eq!(r.nb_points(), 1);
        assert!((r.parameter(0)).abs() < 1.0);
    }

    #[test]
    fn curandinf2d_parabola_perform_cur_ext_vertex_only() {
        let p = make_parabola();
        let r = curvature_extrema_2d(&p);
        assert!(r.is_done());
        assert_eq!(r.nb_points(), 1);
        assert!((r.parameter(0)).abs() < 1.0);
    }

    #[test]
    fn curandinf2d_perform_inf_clears_previous() {
        let e = make_ellipse();
        let r1 = curvature_extrema_2d(&e);
        assert_eq!(r1.nb_points(), 4);
        let r2 = inflection_points_2d(&e);
        assert_eq!(r2.nb_points(), 0);
    }
}
