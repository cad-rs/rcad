//! OCCT HelixGeom_Tools (TKHelix HelixGeom package).
//!
//! 1:1 translation of `HelixGeom_Tools.cxx` — the `HelixGeom_Tools_Eval`
//! evaluator class (L29-91) and the `ApprCurve3D` / `ApprHelix`
//! approximations (L95-190).

use super::helix_curve::HelixCurve;
use rcad_kernel::geom::BSplineCurve3;
use rcad_kernel::math::adv_approx::{ApproxAFunction, DichoCutting, EvaluatorFunction};
use rcad_kernel::math::GeomAbsShape;

/// OCCT class HelixGeom_Tools_Eval — evaluator for the approximation.
struct ToolsEval<'a> {
    fonct: &'a HelixCurve,
}

impl<'a> EvaluatorFunction for ToolsEval<'a> {
    /// OCCT HelixGeom_Tools_Eval::Evaluate (L48-91).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        param: f64,
        order: i32,
        result: &mut [f64],
    ) -> i32 {
        let _ = start_end;
        let par = param;
        let _par = par;

        match order {
            0 => {
                let pnt = self.fonct.eval_d0(par);
                result[0] = pnt.x;
                result[1] = pnt.y;
                result[2] = pnt.z;
                0
            }
            1 => {
                let (_p, v1) = self.fonct.eval_d1(par);
                result[0] = v1.x;
                result[1] = v1.y;
                result[2] = v1.z;
                0
            }
            2 => {
                let (_p, _v1, v2) = self.fonct.eval_d2(par);
                result[0] = v2.x;
                result[1] = v2.y;
                result[2] = v2.z;
                0
            }
            _ => {
                result[0] = 0.0;
                result[1] = 0.0;
                result[2] = 0.0;
                3
            }
        }
    }
}

/// OCCT HelixGeom_Tools::ApprCurve3D (L95-156) — approximates a helix
/// adaptor curve by a BSpline.  Returns `(error_code, bspline, max_error)`;
/// error code 0 on success.
pub fn appr_curve3d(
    the_hc: &HelixCurve,
    the_tol: f64,
    the_cont: GeomAbsShape,
    the_max_seg: i32,
    the_max_deg: i32,
) -> (i32, Option<BSplineCurve3>, f64) {
    let first = the_hc.first_parameter();
    let last = the_hc.last_parameter();
    // Setup approximation dimensions and tolerances: Num3DSS = 1.
    let three_d_tol = [the_tol];

    // Setup approximation function and perform approximation.
    let mut ev = ToolsEval { fonct: the_hc };
    let a_approx = ApproxAFunction::new(
        0,
        0,
        1,
        None,
        None,
        Some(&three_d_tol),
        first,
        last,
        the_cont,
        the_max_deg,
        the_max_seg,
        &mut ev,
    );

    // Check if approximation was successful.
    if !a_approx.is_done() {
        return (1, None, 0.0);
    }
    // Initialize error and check for results.
    if !a_approx.has_result() {
        return (2, None, 0.0);
    }
    // Extract B-spline curve data from approximation.
    let degree = a_approx.degree() as usize;
    let knots = a_approx.knots_vec().to_vec();
    let mults = a_approx.multiplicities_vec().to_vec();
    let poles = a_approx.poles_flat(1);
    let max_error = a_approx.max_error_at(3, 1);
    let control_points: Vec<glam::DVec3> = (0..poles.len() / 3)
        .map(|i| glam::DVec3::new(poles[i * 3], poles[i * 3 + 1], poles[i * 3 + 2]))
        .collect();
    let bspl = BSplineCurve3::from_knots_mults(degree, knots, mults, control_points);
    (0, Some(bspl), max_error)
}

/// OCCT HelixGeom_Tools::ApprHelix (L160-190).
#[allow(clippy::too_many_arguments)]
pub fn appr_helix(
    a_t1: f64,
    a_t2: f64,
    a_pitch: f64,
    a_r_start: f64,
    a_taper_angle: f64,
    a_is_cw: bool,
    the_tol: f64,
) -> (i32, Option<BSplineCurve3>, f64) {
    // Load helix parameters and create adaptor.
    let mut a_adaptor = HelixCurve::new();
    a_adaptor.load(a_t1, a_t2, a_pitch, a_r_start, a_taper_angle, a_is_cw);
    // Set default approximation parameters.
    let a_cont = GeomAbsShape::C2;
    let a_max_degree = 8;
    let a_max_seg = 150;
    // Perform curve approximation.
    appr_curve3d(
        &a_adaptor,
        the_tol,
        a_cont,
        a_max_seg,
        a_max_degree,
    )
}

/// OCCT AdvApprox_DichoCutting default construction site kept local for the
/// builder (the engine ctor takes the cut tool by reference in OCCT).
#[allow(dead_code)]
fn _cut_tool_witness() -> DichoCutting {
    DichoCutting
}
