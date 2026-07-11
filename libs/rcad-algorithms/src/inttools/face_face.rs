use glam::DVec3;
use rcad_kernel::geom::*;
use crate::bopds::ds::DS;
use crate::inttools::intss::{SurfaceSurfaceIntersection, SurfaceCurve};
use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ABS_SQ};
use rcad_kernel::geom::CurveEval;

/// ✅ OCCT-aligned: IntTools_FaceFace — face-face intersection result.
///
/// Wraps the intersection curves with tolerance computation and
/// post-processing matching OCCT's IntTools_FaceFace.
#[derive(Debug, Clone)]
pub struct IntersectionCurve {
    /// 3D curve
    pub curve: Curve3,
    /// Parameter range on the 3D curve
    pub t_range: [f64; 2],
    /// Pcurve on face 1 (if available)
    pub pcurve1: Option<Curve2d>,
    /// Pcurve on face 2 (if available)
    pub pcurve2: Option<Curve2d>,
    /// Computed tolerance (max deviation from surfaces)
    pub tolerance: f64,
    /// Tangential tolerance
    pub tang_tolerance: f64,
}

/// Convert SurfaceCurve to Curve3 with a default range.
fn surface_curve_to_curve3(sc: &SurfaceCurve) -> Option<(Curve3, [f64; 2])> {
    match sc {
        SurfaceCurve::Circle(c) => Some((Curve3::Circle(*c), [0.0, std::f64::consts::TAU])),
        SurfaceCurve::Ellipse(e) => Some((Curve3::Ellipse(*e), [0.0, std::f64::consts::TAU])),
        SurfaceCurve::Line(l) => {
            // OCCT: infinite line, use a large range
            Some((Curve3::Line(*l), [-1e5, 1e5]))
        }
        SurfaceCurve::Parabola(p) => Some((Curve3::Parabola(*p), [-1e5, 1e5])),
        SurfaceCurve::Hyperbola(h) => Some((Curve3::Hyperbola(*h), [-1e5, 1e5])),
        SurfaceCurve::BSplineCurve(b) => {
            let bs: &BSplineCurve3 = &*b;
            let t_range = [bs.knots[0], bs.knots[bs.knots.len() - 1]];
            Some((Curve3::BSpline(bs.clone()), t_range))
        }
        SurfaceCurve::Polyline(pts) => {
            let curve3 = crate::inttools::intss::polyline_to_bspline(pts, 1e-4)?;
            let t_range = curve3.default_domain();
            Some((curve3, t_range))
        }
        SurfaceCurve::Point(_) => None,
    }
}

/// ✅ OCCT-aligned: ComputeTolReached3d (IntTools_FaceFace.cxx L613-691).
///
/// Computes valid tolerance for each intersection curve as the max deviation
/// between the 3D curve and the 2D pcurves on both surfaces.
pub fn compute_tol_reached_3d(
    curves: &mut [IntersectionCurve],
    surf1: &Surface3,
    surf2: &Surface3,
    tol_f1: f64,
    tol_f2: f64,
) {
    use rcad_kernel::geom::SurfaceEval;
    let a_tol_f_max = tol_f1.max(tol_f2);
    for ic in curves.iter_mut() {
        let mut a_tol_c = ic.tolerance;
        let t1 = ic.t_range[0];
        let t2 = ic.t_range[1];
        let n = 23usize;
        let dt = (t2 - t1) / n as f64;
        // Check pcurve on surf1
        if let Some(ref pc) = ic.pcurve1 {
            let mut max_d = 0.0;
            for i in 0..=n {
                let t = t1 + i as f64 * dt;
                let p3d = ic.curve.point_at(t);
                let uv = pc.point_at(t);
                let surf_pt = surf1.point_at(uv.x, uv.y);
                let d = (p3d - surf_pt).length();
                if d > max_d { max_d = d; }
            }
            if max_d > a_tol_c { a_tol_c = max_d; }
        } else {
            let max_d = find_max_distance(&ic.curve, t1, t2, surf1);
            if max_d > a_tol_c { a_tol_c = max_d; }
        }
        // Check pcurve on surf2
        if let Some(ref pc) = ic.pcurve2 {
            let mut max_d = 0.0;
            for i in 0..=n {
                let t = t1 + i as f64 * dt;
                let p3d = ic.curve.point_at(t);
                let uv = pc.point_at(t);
                let surf_pt = surf2.point_at(uv.x, uv.y);
                let d = (p3d - surf_pt).length();
                if d > max_d { max_d = d; }
            }
            if max_d > a_tol_c { a_tol_c = max_d; }
        } else {
            let max_d = find_max_distance(&ic.curve, t1, t2, surf2);
            if max_d > a_tol_c { a_tol_c = max_d; }
        }
        ic.tolerance = a_tol_c;
        if ic.tang_tolerance < a_tol_f_max {
            ic.tang_tolerance = a_tol_f_max;
        }
    }
}

fn find_max_distance(curve: &Curve3, t1: f64, t2: f64, surf: &Surface3) -> f64 {
    let n = 23usize;
    let dt = (t2 - t1) / n as f64;
    let mut max_d = 0.0;
    for i in 0..=n {
        let t = t1 + i as f64 * dt;
        let p = curve.point_at(t);
        let proj = rcad_kernel::projection::closest_point_on_surface(surf, p, 16);
        if proj.distance.is_finite() && proj.distance > max_d {
            max_d = proj.distance;
        }
    }
    max_d
}

/// ✅ OCCT-aligned: PrepareLines3D (IntTools_FaceFace.cxx L1932-2015).
pub fn prepare_lines_3d(curves: &mut Vec<IntersectionCurve>, b_to_split: bool) {
    if !b_to_split { return; }
    let mut i = 0;
    while i < curves.len() {
        let (t0, t1) = (curves[i].t_range[0], curves[i].t_range[1]);
        // OCCT IntTools_Tools::IsClosed: check if start/end points coincide
        let is_closed = {
            let c = &curves[i];
            let p1 = c.curve.point_at(t0);
            let p2 = c.curve.point_at(t1);
            (p1 - p2).length_squared() < TOLERANCE_ABS_SQ
        };
        if is_closed {
            // OCCT L214-221: for BSpline/Bezier use IntermediatePoint, else regular midpoint
            let t_mid = match &curves[i].curve {
                Curve3::BSpline(_) | Curve3::Bezier(_) => {
                    0.56786082 * t0 + 0.43213918 * t1
                }
                _ => 0.5 * (t0 + t1),
            };
            let c2 = IntersectionCurve {
                curve: curves[i].curve.clone(),
                t_range: [t_mid, t1],
                pcurve1: curves[i].pcurve1.clone(),
                pcurve2: curves[i].pcurve2.clone(),
                tolerance: curves[i].tolerance,
                tang_tolerance: curves[i].tang_tolerance,
            };
            curves[i].t_range[1] = t_mid;
            curves.insert(i + 1, c2);
            i += 2;
        } else {
            i += 1;
        }
    }
}

/// ✅ OCCT-aligned: Convert SurfaceSurfaceIntersection to IntersectionCurve list.
pub fn intersection_to_curves(intss: &SurfaceSurfaceIntersection) -> Vec<IntersectionCurve> {
    let mut curves = Vec::new();
    for result in &intss.curves {
        if let Some((curve, t_range)) = surface_curve_to_curve3(&result.curve_3d) {
            curves.push(IntersectionCurve {
                curve,
                t_range,
                pcurve1: result.pcurve_on_a.clone(),
                pcurve2: result.pcurve_on_b.clone(),
                tolerance: 1e-7,
                tang_tolerance: 1e-7,
            });
        }
    }
    curves
}

/// ✅ OCCT-aligned: CorrectSurfaceBoundaries (IntTools_FaceFace.cxx L2017-2210).
pub fn correct_surface_boundaries(uv_bounds: &mut [f64; 4], margin: f64) {
    let [u_min, u_max, v_min, v_max] = *uv_bounds;
    let du = u_max - u_min;
    let dv = v_max - v_min;
    let m = margin * du.max(dv).max(1.0);
    uv_bounds[0] = u_min - m;
    uv_bounds[1] = u_max + m;
    uv_bounds[2] = v_min - m;
    uv_bounds[3] = v_max + m;
}

/// ✅ OCCT-aligned: CorrectPlaneBoundaries (IntTools_FaceFace.cxx L3093-3111).
pub fn correct_plane_boundaries(uv_bounds: &mut [f64; 4]) {
    uv_bounds[0] = -1e10;
    uv_bounds[1] = 1e10;
    uv_bounds[2] = -1e10;
    uv_bounds[3] = 1e10;
}

/// ✅ OCCT-aligned: FaceFaceIntersector — top-level dispatch.
pub fn intersect_faces(
    surf1: &Surface3,
    surf2: &Surface3,
    tol_f1: f64,
    tol_f2: f64,
) -> Vec<IntersectionCurve> {
    let intss = crate::inttools::intss::intersect_surfaces_with_tolerance(surf1, surf2, tol_f1 + tol_f2);
    let mut curves = intersection_to_curves(&intss);
    compute_tol_reached_3d(&mut curves, surf1, surf2, tol_f1, tol_f2);
    prepare_lines_3d(&mut curves, true);
    curves
}
