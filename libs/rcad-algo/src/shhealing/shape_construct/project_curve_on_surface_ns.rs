//! OCCT ShapeConstruct_ProjectCurveOnSurface.cxx anonymous namespace
//! (L68-637) plus the file-local helpers, translated 1:1.
//!
//! Contents:
//! - `THE_NCONTROL` (L71)
//! - `SurfaceProjectorWithCache` (L78-224)
//! - `extractBSplineCurve` (L232-243)
//! - `adjustSecondToFirstPoint` (L253-273)
//! - `fixPeriodicityTroubles` (L286-415)
//! - `isBSplineCurveInvalid` (L428-517)
//! - `generateCurvePoints` (L531-586)
//! - `projectDegeneratedPoints` (L600-635)

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{BSplineCurve3, BSplineSurface, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::el::in_period;
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};

use super::array1::Array1;
use super::gap_deps::{
    shape_analysis_adjust_to_period, surface_u_period, surface_v_period, ShapeAnalysisSurface,
};

/// OCCT L71: default number of control points for discretization.
pub const THE_NCONTROL: i32 = 23;

/// OCCT L60: `CachePoint = std::pair<gp_Pnt, gp_Pnt2d>`.
#[derive(Debug, Clone, Copy)]
pub struct CachePoint {
    pub first: DVec3,
    pub second: DVec2,
}

// OCCT's Geom_BSpline(Curve/Surface) carry (knots, multiplicities); the rcad
// BSpline encoding stores the flat expanded knot vector (arch. diff.).  The
// two helpers below re-derive the OCCT-form quantities from the flat vector.

/// Distinct knot values of a flat knot vector (`Geom_BSplineCurve::Knot(i)`).
pub fn bspline_distinct_knots(flat_knots: &[f64]) -> Vec<f64> {
    let mut knots: Vec<f64> = Vec::new();
    for &k in flat_knots {
        if knots.last().map(|last: &f64| *last != k).unwrap_or(true) {
            knots.push(k);
        }
    }
    knots
}

/// First-knot multiplicity of a flat knot vector.
fn flat_first_mult(flat_knots: &[f64]) -> usize {
    let mut m = 1usize;
    while m < flat_knots.len() && flat_knots[m] == flat_knots[0] {
        m += 1;
    }
    m
}

/// Last-knot multiplicity of a flat knot vector.
fn flat_last_mult(flat_knots: &[f64]) -> usize {
    let mut m = 1usize;
    while m < flat_knots.len() && flat_knots[flat_knots.len() - m] == flat_knots[flat_knots.len() - 1]
    {
        m += 1;
    }
    m
}

/// OCCT L78-224: utility class for projecting points onto a surface with a
/// B-spline corner cache optimization.  For clamped B-spline surfaces, caches
/// corner pole positions and their exact UV parameters to avoid expensive
/// ValueOfUV calls when projected points coincide with surface corners.
pub struct SurfaceProjectorWithCache<'a> {
    /// OCCT mySurf — surface to project on.
    my_surf: Option<&'a ShapeAnalysisSurface>,
    /// OCCT myCorners3d — 3D positions of B-spline surface corners.
    my_corners3d: Vec<DVec3>,
    /// OCCT myCorners2d — UV parameters of B-spline surface corners.
    my_corners2d: Vec<DVec2>,
}

impl<'a> SurfaceProjectorWithCache<'a> {
    /// OCCT L84-98: constructor — initializes the projector with a surface.
    /// For B-spline surfaces, builds the corner cache automatically.
    pub fn new(the_surf: Option<&'a ShapeAnalysisSurface>) -> Self {
        let mut prj = SurfaceProjectorWithCache {
            my_surf: the_surf,
            my_corners3d: Vec::new(),
            my_corners2d: Vec::new(),
        };
        if let Some(a_surf) = the_surf {
            if let Surface3::BSpline(a_bspline_surf) = surf_surface(a_surf) {
                prj.build_corner_cache(a_bspline_surf);
            }
        }
        prj
    }

    /// OCCT L106-115: projects a 3D point onto the surface.  First checks the
    /// B-spline corner cache, then falls back to ValueOfUV.
    pub fn value_of_uv(&self, the_point: DVec3, the_tol: f64, the_tol_sq: f64) -> DVec2 {
        let surf = match self.my_surf {
            Some(s) => s,
            None => return DVec2::ZERO,
        };
        let mut a_result = DVec2::ZERO;
        if self.find_in_corner_cache(the_point, the_tol_sq, &mut a_result) {
            // Corner UV is exact, but refine for numerical stability
            return surf.next_value_of_uv(a_result, the_point, the_tol, the_tol);
        }
        surf.value_of_uv(the_point, the_tol)
    }

    /// OCCT L125-138: projects a 3D point onto the surface using a hint from a
    /// previous projection.  First checks the B-spline corner cache, then
    /// falls back to NextValueOfUV.
    pub fn next_value_of_uv(
        &self,
        the_hint: DVec2,
        the_point: DVec3,
        the_tol: f64,
        the_tol_sq: f64,
        the_step: f64,
    ) -> DVec2 {
        let surf = match self.my_surf {
            Some(s) => s,
            None => return DVec2::ZERO,
        };
        let mut a_result = DVec2::ZERO;
        if self.find_in_corner_cache(the_point, the_tol_sq, &mut a_result) {
            // Corner UV is exact, but refine for numerical stability
            return surf.next_value_of_uv(a_result, the_point, the_tol, the_tol);
        }
        surf.next_value_of_uv(the_hint, the_point, the_tol, the_step)
    }

    /// OCCT L141: returns the gap from the last projection.
    pub fn gap(&self) -> f64 {
        match self.my_surf {
            Some(s) => s.gap(),
            None => 0.0,
        }
    }

    /// OCCT L148-200: builds the cache from B-spline surface corner poles.
    /// For clamped B-splines (multiplicity = degree + 1 at ends), the surface
    /// passes exactly through corner poles, so their UV parameters can be
    /// used directly.
    fn build_corner_cache(&mut self, the_surface: &BSplineSurface) {
        let a_nb_poles_u = the_surface.control_points.len();
        let a_nb_poles_v = the_surface.control_points.first().map(|r| r.len()).unwrap_or(0);
        let a_degree_u = the_surface.degree_u;
        let a_degree_v = the_surface.degree_v;

        // Check if surface is clamped (end multiplicities = degree + 1)
        let a_first_u_mult = flat_first_mult(&the_surface.knots_u);
        let a_last_u_mult = flat_last_mult(&the_surface.knots_u);
        let a_first_v_mult = flat_first_mult(&the_surface.knots_v);
        let a_last_v_mult = flat_last_mult(&the_surface.knots_v);

        let is_u_first_clamped = a_first_u_mult >= a_degree_u + 1;
        let is_u_last_clamped = a_last_u_mult >= a_degree_u + 1;
        let is_v_first_clamped = a_first_v_mult >= a_degree_v + 1;
        let is_v_last_clamped = a_last_v_mult >= a_degree_v + 1;

        // Get parameter bounds
        let a_u_first = the_surface.knots_u[0];
        let a_u_last = the_surface.knots_u[the_surface.knots_u.len() - 1];
        let a_v_first = the_surface.knots_v[0];
        let a_v_last = the_surface.knots_v[the_surface.knots_v.len() - 1];

        // Cache corner poles where surface passes through them
        // Corner (1, 1) - UFirst, VFirst
        if is_u_first_clamped && is_v_first_clamped {
            self.my_corners3d.push(the_surface.control_points[0][0]);
            self.my_corners2d.push(DVec2::new(a_u_first, a_v_first));
        }

        // Corner (NbPolesU, 1) - ULast, VFirst
        if is_u_last_clamped && is_v_first_clamped {
            self.my_corners3d
                .push(the_surface.control_points[a_nb_poles_u - 1][0]);
            self.my_corners2d.push(DVec2::new(a_u_last, a_v_first));
        }

        // Corner (1, NbPolesV) - UFirst, VLast
        if is_u_first_clamped && is_v_last_clamped {
            self.my_corners3d
                .push(the_surface.control_points[0][a_nb_poles_v - 1]);
            self.my_corners2d.push(DVec2::new(a_u_first, a_v_last));
        }

        // Corner (NbPolesU, NbPolesV) - ULast, VLast
        if is_u_last_clamped && is_v_last_clamped {
            self.my_corners3d
                .push(the_surface.control_points[a_nb_poles_u - 1][a_nb_poles_v - 1]);
            self.my_corners2d.push(DVec2::new(a_u_last, a_v_last));
        }
    }

    /// OCCT L207-218: finds the closest corner to a 3D point within tolerance.
    fn find_in_corner_cache(&self, the_point: DVec3, the_tol_sq: f64, the_uv: &mut DVec2) -> bool {
        for i in 0..self.my_corners3d.len() {
            if self.my_corners3d[i].distance_squared(the_point) < the_tol_sq {
                *the_uv = self.my_corners2d[i];
                return true;
            }
        }
        false
    }
}

/// OCCT `theSurf->Surface()` member read through the carrier.
fn surf_surface(the_surf: &ShapeAnalysisSurface) -> &Surface3 {
    the_surf.surface()
}

/// OCCT L232-243: extracts a B-spline curve from a possibly nested trimmed
/// curve.  Recursively unwraps trimmed curves to find the B-spline basis.
pub fn extract_bspline_curve(the_curve: &Curve3) -> Option<BSplineCurve3> {
    let mut a_curve = the_curve.clone();

    // Recursively unwrap trimmed curves
    while let Curve3::Trimmed(tc) = &a_curve {
        a_curve = tc.curve.as_ref().clone();
    }

    match &a_curve {
        Curve3::BSpline(bs) => Some(bs.clone()),
        _ => None,
    }
}

/// OCCT L253-273: adjusts the second point to the first point considering
/// surface periodicity.  For periodic surfaces, adjusts the second point
/// coordinates to be within half a period of the first point.
pub fn adjust_second_to_first_point(
    the_first_point: DVec2,
    the_second_point: &mut DVec2,
    the_surf: &Surface3,
) {
    if the_surf.is_u_periodic() {
        let an_u_period = surface_u_period(the_surf);
        let a_new_u = in_period(
            the_second_point.x,
            the_first_point.x - an_u_period / 2.0,
            the_first_point.x + an_u_period / 2.0,
        );
        the_second_point.x = a_new_u;
    }
    if the_surf.is_v_periodic() {
        let a_v_period = surface_v_period(the_surf);
        let a_new_v = in_period(
            the_second_point.y,
            the_first_point.y - a_v_period / 2.0,
            the_first_point.y + a_v_period / 2.0,
        );
        the_second_point.y = a_new_v;
    }
}

/// Coordinate access by OCCT's 1-based `gp_Pnt2d::Coord(theIdx)`
/// (theIdx = 1 for X, 2 for Y).
pub(crate) fn p2d_coord(the_p: DVec2, the_idx: usize) -> f64 {
    if the_idx == 1 {
        the_p.x
    } else {
        the_p.y
    }
}

/// Coordinate write by OCCT's 1-based `gp_Pnt2d::SetCoord(theIdx, val)`.
pub(crate) fn p2d_set_coord(the_p: &mut DVec2, the_idx: usize, the_val: f64) {
    if the_idx == 1 {
        the_p.x = the_val;
    } else {
        the_p.y = the_val;
    }
}

/// OCCT L286-415: fixes possible period jumps in an array of 4 points.
/// Handles the walking period parameter to ensure smooth transition across
/// periodic boundaries.  Returns true if a period jump was detected and fixed.
pub fn fix_periodicity_troubles(
    the_pnt: &mut [DVec2; 4],
    the_idx: usize,
    the_period: f64,
    the_saved_point: i32,
    the_saved_param: f64,
) -> bool {
    let a_saved_param: f64;
    let a_saved_point: i32;
    let mut a_min_param = 0.0f64;
    let mut a_max_param = the_period;

    if the_saved_point < 0 {
        a_saved_param = 0.5 * the_period;
        a_saved_point = 0;
    } else {
        a_saved_param = the_saved_param;
        a_saved_point = the_saved_point;
        while a_min_param > a_saved_param {
            a_min_param -= the_period;
            a_max_param -= the_period;
        }
        while a_max_param < a_saved_param {
            a_min_param += the_period;
            a_max_param += the_period;
        }
    }

    let mut a_fix_iso_param = a_min_param;
    let mut is_iso_line = false;
    if a_max_param - a_saved_param < PCONFUSION || a_saved_param - a_min_param < PCONFUSION {
        a_fix_iso_param = a_saved_param;
        is_iso_line = true;
    }

    // Normalize all coordinates to [aMinParam, aMaxParam)
    for i in 0..4 {
        let mut a_param = p2d_coord(the_pnt[i], the_idx);
        let a_shift = shape_analysis_adjust_to_period(a_param, a_min_param, a_max_param);
        a_param += a_shift;

        if is_iso_line {
            if a_max_param - a_param < PCONFUSION || a_param - a_min_param < PCONFUSION {
                a_param = a_fix_iso_param;
            }
        } else {
            if a_max_param - a_param < PCONFUSION {
                a_param = a_max_param;
            }
            if a_param - a_min_param < PCONFUSION {
                a_param = a_min_param;
            }
        }

        p2d_set_coord(&mut the_pnt[i], the_idx, a_param);
    }

    // Find possible period jump and increasing or decreasing coordinates vector
    let mut is_jump = false;
    let mut a_prev_diff = 0.0f64;
    let mut a_sum_diff = 1.0f64;
    for i in 0..3 {
        let a_diff = p2d_coord(the_pnt[i + 1], the_idx) - p2d_coord(the_pnt[i], the_idx);
        if a_diff < -PCONFUSION {
            a_sum_diff *= -1.0;
        }
        if a_diff * a_prev_diff < -PCONFUSION {
            is_jump = true;
        }
        a_prev_diff = a_diff;
    }

    if !is_jump {
        return false;
    }

    if a_sum_diff > 0.0 {
        let mut i = a_saved_point;
        while i > 0 {
            if p2d_coord(the_pnt[i as usize], the_idx)
                > p2d_coord(the_pnt[(i - 1) as usize], the_idx)
            {
                let v = p2d_coord(the_pnt[(i - 1) as usize], the_idx) + the_period;
                p2d_set_coord(&mut the_pnt[(i - 1) as usize], the_idx, v);
            }
            i -= 1;
        }
        let mut i = a_saved_point;
        while i < 3 {
            if p2d_coord(the_pnt[i as usize], the_idx) < p2d_coord(the_pnt[(i + 1) as usize], the_idx)
            {
                let v = p2d_coord(the_pnt[(i + 1) as usize], the_idx) - the_period;
                p2d_set_coord(&mut the_pnt[(i + 1) as usize], the_idx, v);
            }
            i += 1;
        }
    } else {
        let mut i = a_saved_point;
        while i > 0 {
            if p2d_coord(the_pnt[i as usize], the_idx)
                < p2d_coord(the_pnt[(i - 1) as usize], the_idx)
            {
                let v = p2d_coord(the_pnt[(i - 1) as usize], the_idx) - the_period;
                p2d_set_coord(&mut the_pnt[(i - 1) as usize], the_idx, v);
            }
            i -= 1;
        }
        let mut i = a_saved_point;
        while i < 3 {
            if p2d_coord(the_pnt[i as usize], the_idx) > p2d_coord(the_pnt[(i + 1) as usize], the_idx)
            {
                let v = p2d_coord(the_pnt[(i + 1) as usize], the_idx) + the_period;
                p2d_set_coord(&mut the_pnt[(i + 1) as usize], the_idx, v);
            }
            i += 1;
        }
    }

    true
}

/// OCCT L428-517: checks if a B-spline curve has uneven parameterization
/// requiring special handling.  Computes the ratio of maximum to minimum
/// parameterization speed across knot intervals.  If this ratio exceeds a
/// threshold, the curve should be projected using ProjLib instead of the
/// standard approximation approach.
pub fn is_bspline_curve_invalid(
    the_curve: &Curve3,
    the_first: f64,
    the_last: f64,
    the_bspline: &mut Option<BSplineCurve3>,
) -> bool {
    *the_bspline = extract_bspline_curve(the_curve);
    let the_bs = match the_bspline {
        Some(bs) => bs,
        None => return false,
    };

    // Compute parametrization speed on each knot interval inside [theFirst,
    // theLast].  If quotient = (MaxSpeed / MinSpeed) >= aMaxQuotientCoeff then
    // use PerformByProjLib.
    let mut a_first_param = the_first;
    let mut a_last_param = the_last;

    // OCCT iterates over the distinct knot values (Geom_BSplineCurve::Knot(i));
    // rcad stores flat knots (arch. diff.), so the distinct knots are derived.
    let a_knots = bspline_distinct_knots(&the_bs.knots);
    let a_nb_knots = a_knots.len() as i32;

    // Find first knot index
    let mut an_idx = 1i32;
    while an_idx <= a_nb_knots && a_first_param < the_last {
        if a_knots[(an_idx - 1) as usize] > the_first {
            break;
        }
        an_idx += 1;
    }

    let mut a_min_par_speed = f64::INFINITY; // OCCT Precision::Infinite()
    let mut a_knot_coeffs: Vec<f64> = Vec::new();

    while an_idx <= a_nb_knots && a_first_param < the_last {
        // Fill current knot interval
        a_last_param = the_last.min(a_knots[(an_idx - 1) as usize]);
        let mut a_nb_int_pnts = THE_NCONTROL;

        // Adapt number of inner points according to the length of the interval
        // to avoid a lot of calculations on small range of parameters.
        if an_idx > 1 {
            let a_len_thres = 1.0e-2;
            let a_len_ratio = (a_last_param - a_first_param)
                / (a_knots[(an_idx - 1) as usize] - a_knots[(an_idx - 2) as usize]);
            if a_len_ratio < a_len_thres {
                a_nb_int_pnts = (a_len_ratio / a_len_thres * a_nb_int_pnts as f64) as i32;
                if a_nb_int_pnts < 2 {
                    a_nb_int_pnts = 2;
                }
            }
        }

        let a_step = (a_last_param - a_first_param) / (a_nb_int_pnts - 1) as f64;
        let mut p3d1 = DVec3::ZERO;
        let mut p3d2 = DVec3::ZERO;

        // Start filling from first point
        p3d1 = the_curve.point_at(a_first_param);

        let mut a_length3d = 0.0f64;
        for an_int_idx in 1..a_nb_int_pnts {
            let a_param = a_first_param + a_step * an_int_idx as f64;
            p3d2 = the_curve.point_at(a_param);
            let a_dist = p3d2.distance(p3d1);

            a_length3d += a_dist;
            p3d1 = p3d2;

            a_min_par_speed = a_min_par_speed.min(a_dist / a_step);
        }

        let a_coeff = a_length3d / (a_last_param - a_first_param);
        if a_coeff.abs() > GP_RESOLUTION {
            a_knot_coeffs.push(a_coeff);
        }
        a_first_param = a_last_param;
        an_idx += 1;
    }

    let mut an_evenly_coeff = 0.0f64;
    if !a_knot_coeffs.is_empty() {
        let max_c = a_knot_coeffs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_c = a_knot_coeffs.iter().cloned().fold(f64::INFINITY, f64::min);
        an_evenly_coeff = max_c / min_c;
    }

    let a_max_quotient_coeff = 1500.0;
    an_evenly_coeff > a_max_quotient_coeff && a_min_par_speed > CONFUSION
}

/// OCCT L531-586: generates discretization points for a curve.  Uses uniform
/// distribution for general curves, and adjusts the number of points for
/// B-splines based on their knot structure to ensure adequate sampling.
/// Returns the actual number of generated points.
pub fn generate_curve_points(
    the_curve: &Curve3,
    the_first: f64,
    the_last: f64,
    the_nb_control_points: i32,
    the_points: &mut Array1<DVec3>,
    the_params: &mut Array1<f64>,
) -> i32 {
    let a_bspline = extract_bspline_curve(the_curve);
    let mut a_nb_pini = the_nb_control_points;

    if let Some(bs) = &a_bspline {
        let a_knots = bspline_distinct_knots(&bs.knots);
        let a_nb_knots = a_knots.len() as i32;
        let mut a_used_knots = 0i32;
        for i in 1..a_nb_knots {
            if a_knots[i as usize] > the_first && a_knots[(i - 1) as usize] < the_last {
                a_used_knots += 1;
            }
        }
        let a_min_pnt = a_used_knots * (bs.degree as i32 + 1);
        while a_nb_pini < a_min_pnt {
            a_nb_pini += THE_NCONTROL - 1;
        }
    }

    the_points.resize(1, a_nb_pini as usize, DVec3::ZERO);
    the_params.resize(1, a_nb_pini as usize, 0.0);

    let a_delta_param = (the_last - the_first) / (a_nb_pini - 1) as f64;

    for i in 1..=a_nb_pini {
        let a_param = if i == 1 {
            the_first
        } else if i == a_nb_pini {
            the_last
        } else {
            the_first + (i - 1) as f64 * a_delta_param
        };

        let a_point = the_curve.point_at(a_param);
        the_points.set_value(i as usize, a_point);
        the_params.set_value(i as usize, a_param);
    }

    a_nb_pini
}

/// OCCT L600-635: wrapper for ShapeAnalysis_Surface::ProjectDegenerated.
/// Converts NCollection_Array1 containers to sequences, performs projection,
/// then copies results back.  Required because ShapeAnalysis_Surface uses
/// sequence containers.
pub fn project_degenerated_points(
    the_surf: &ShapeAnalysisSurface,
    the_nb_pnt: i32,
    the_points: &Array1<DVec3>,
    the_points2d: &mut Array1<DVec2>,
    the_preci: f64,
    the_direct: bool,
) {
    // Convert arrays to sequences for ShapeAnalysis_Surface::ProjectDegenerated
    let mut a_points3d: Vec<DVec3> = Vec::new();
    let mut a_points2d: Vec<DVec2> = Vec::new();

    for i in 1..=the_nb_pnt {
        a_points3d.push(*the_points.value(i as usize));
        a_points2d.push(*the_points2d.value(i as usize));
    }

    // Call the method that expects sequences
    the_surf.project_degenerated(the_nb_pnt, &a_points3d, &mut a_points2d, the_preci, the_direct);

    // Copy results back to array
    for i in 1..=the_nb_pnt {
        the_points2d.set_value(i as usize, a_points2d[(i - 1) as usize]);
    }
}
