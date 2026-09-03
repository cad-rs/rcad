//! OCCT GeomFill_CorrectedFrenet (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_CorrectedFrenet.hxx (L31-118) + GeomFill_CorrectedFrenet.cxx
//! (whole file L26-1064), together with the GeomFill_Trihedron enum
//! (GeomFill_Trihedron.hxx L18-30) and the static ComputeTorsion /
//! smoothlaw / FindPlane helpers.
//!
//! Architecture mappings: `Adaptor3d_Curve` -> rcad `Curve3`; the
//! `myTrimmed` view is `Curve3::Trimmed` evaluated at unchanged parameters;
//! `BndLib_Add3dCurve::Add(SnglrFunc-as-curve, ...)` ->
//! `bnd_lib::curve_box_range_fn` with the SnglrFunc evaluation closure
//! (LengthMin = GetGap() * 1e-4, the gap being the added tolerance).
//! `TLaw = EvolAroundT` copies the composite in Rust (in OCCT the two
//! handles alias one object; the alias is not mutated after Init except
//! through the TLaw path itself, and `Prepare` re-locates its cached
//! function when out of range, so the copy is behaviorally neutral).

use std::cell::RefCell;
use std::f64::consts::PI;
use std::rc::Rc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::bnd_lib::curve_box_range_fn;
use rcad_kernel::base::geom_lib::{axe_of_inertia, fuse_intervals};
use rcad_kernel::geom::{Curve3, CurveEval, Plane};
use rcad_kernel::math::gp::Ax2;
use rcad_kernel::math::GeomAbsShape;

use super::frenet::{gp_vec_angle, law_d3, Frenet};
use crate::geomalgo::law::{
    LawBSpFunc, LawBSpline, LawComposite, LawConstant, LawFunction, LawFunctionHandle,
    LawInterpolate,
};
use super::sngrl_func::SnglrFunc;
use super::trihedron_law::{
    curve_first_parameter, curve_last_parameter, TrihedronLaw, TrihedronLawBase,
};

// OCCT Precision values used below.
const P_CONFUSION: f64 = 1e-12;
const SQUARE_CONFUSION: f64 = 1e-14;
const ANGULAR: f64 = 1e-12;
const CONFUSION: f64 = 1e-7;
// OCCT gp::Resolution().
const GP_RESOLUTION: f64 = 2.2250738585072014e-308;

/// OCCT GeomFill_Trihedron (GeomFill_Trihedron.hxx L18-30).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trihedron {
    IsCorrectedFrenet,
    IsFixed,
    IsFrenet,
    IsConstantNormal,
    IsDarboux,
    IsGuideAC,
    IsGuidePlan,
    IsGuideACWithContact,
    IsGuidePlanWithContact,
    IsDiscreteTrihedron,
}

/// OCCT static ComputeTorsion (L52-72).
fn compute_torsion(param: f64, curve: &Curve3) -> f64 {
    let (_, dc1, dc2, dc3) = law_d3(curve, param);
    let dc1crossdc2 = dc1.cross(dc2);
    let norm_dc1crossdc2 = dc1crossdc2.length();
    let dc1dc2dc3 = dc1crossdc2.dot(dc3); // mixed product
    let tol = GP_RESOLUTION;
    let square_norm_dc1crossdc2 = norm_dc1crossdc2 * norm_dc1crossdc2;
    if square_norm_dc1crossdc2 <= tol {
        0.0
    } else {
        dc1dc2dc3 / square_norm_dc1crossdc2
    }
}

/// OCCT static smoothlaw (L83-190) — smooth a law: reduce the number of
/// knots (two RemoveKnot passes with an interpolation control between them).
fn smoothlaw(law: &mut LawBSpline, points: &[f64], param: &[f64], tol: f64) {
    let tol_orig = tol;
    let mut bs = law.clone();
    let nbk = bs.nb_knots();
    let mut tol = tol / 10.0;
    let mut ok = false;
    // Une premiere passe tolerance serres
    let mut ii = nbk as i32 - 1;
    while ii > 1 {
        let b = bs.remove_knot(ii as usize, 0, tol);
        if b {
            ok = true;
        }
        ii -= 1;
    }
    if ok {
        // controle
        let mut dev = 0.0f64;
        for ii in 0..param.len() {
            let d = (bs.value(param[ii]) - points[ii]).abs();
            if d > dev {
                dev = d;
            }
            ok = dev <= tol_orig;
            if !ok {
                break;
            }
        }
        if ok {
            tol = (tol_orig - dev) / 2.0;
        } else {
            // Echec
            return;
        }
    } else {
        tol = tol_orig / 2.0;
    }
    if ok {
        *law = bs.clone();
    }
    // Une deuxieme passe tolerance desserre
    ok = false;
    let nbk = bs.nb_knots();
    let mut ii = nbk as i32 - 1;
    while ii > 1 {
        let b = bs.remove_knot(ii as usize, 0, tol);
        if b {
            ok = true;
        }
        ii -= 1;
    }
    if ok {
        // controle
        let mut dev = 0.0f64;
        for ii in 0..param.len() {
            let d = (bs.value(param[ii]) - points[ii]).abs();
            if d > dev {
                dev = d;
            }
            ok = dev <= tol_orig;
            if !ok {
                break;
            }
        }
    }
    if ok {
        *law = bs;
    }
}

/// OCCT static FindPlane (L189-290) — find an average plane for the curve.
fn find_plane(curve: &Curve3) -> Option<Plane> {
    let mut found = true;
    let tab_p: Option<Vec<DVec3>>;
    match curve {
        Curve3::Line(_) => {
            found = false;
            tab_p = None;
        }
        Curve3::Circle(c) => {
            return Some(Plane::new(c.center, c.normal));
        }
        Curve3::Ellipse(e) => {
            return Some(Plane::new(e.center, e.normal));
        }
        Curve3::Hyperbola(h) => {
            return Some(Plane::new(h.center, h.normal));
        }
        Curve3::Parabola(p) => {
            return Some(Plane::new(p.vertex, p.normal));
        }
        Curve3::Bezier(gc) => {
            let nbp = gc.control_points.len();
            if nbp < 2 || nbp == 2 {
                found = false;
                tab_p = None;
            } else {
                tab_p = Some(gc.control_points.clone());
            }
        }
        Curve3::BSpline(gc) => {
            let nbp = gc.control_points.len();
            if nbp < 2 || nbp == 2 {
                found = false;
                tab_p = None;
            } else {
                tab_p = Some(gc.control_points.clone());
            }
        }
        _ => {
            // On utilise un echantillonnage
            let nbp = 15 + super::frenet::curve_nb_intervals(curve, GeomAbsShape::C3);
            let f = curve_first_parameter(curve);
            let l = curve_last_parameter(curve);
            let inv = 1.0 / (nbp - 1) as f64;
            let mut samples = Vec::with_capacity(nbp);
            for ii in 1..=nbp {
                let mut t = f * ((nbp - ii) as f64) + l * ((ii - 1) as f64);
                t *= inv;
                samples.push(curve.point_at(t));
            }
            tab_p = Some(samples);
        }
    }
    if let Some(points) = tab_p {
        // Recherche d'un plan moyen et controle
        let mut inertia = Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X);
        let mut issingular = false;
        // OCCT: GeomLib::AxeOfInertia(TabP->Array1(), inertia, issingular)
        // with the default Tol = 1.0e-7.
        axe_of_inertia(&points, &mut inertia, &mut issingular, 1.0e-7);
        if issingular {
            found = false;
        }
        if found {
            let the_p = Plane::new(inertia.location, inertia.direction);
            // theP->Coefficients(a, b, c, d): N = (a, b, c),
            // d = -N . Location.
            let (a, b, c) = (the_p.normal.x, the_p.normal.y, the_p.normal.z);
            let d = -the_p.normal.dot(the_p.origin);
            for xyz in &points {
                let dist = a * xyz.x + b * xyz.y + c * xyz.z + d;
                found = dist.abs() <= CONFUSION;
                if !found {
                    break;
                }
            }
            if found {
                return Some(the_p);
            }
            return None;
        }
    }
    if found {
        None
    } else {
        None
    }
}

/// OCCT static corr2PI_PI (L678-680).
fn corr_2pi_pi(ang: f64) -> f64 {
    if ang < PI {
        ang
    } else {
        ang - 2.0 * PI
    }
}

/// OCCT static diffAng (L682-685).
fn diff_ang(a: f64, ao: f64) -> f64 {
    let mut da = (a - ao) - ((a - ao) / 2.0 / PI).floor() * 2.0 * PI;
    da = if da >= 0.0 {
        corr_2pi_pi(da)
    } else {
        -corr_2pi_pi(-da)
    };
    da
}

/// OCCT gp_Vec2d::Angle — the angle between two 2D vectors in (-PI, PI].
fn gp_vec2d_angle(a: DVec2, b: DVec2) -> f64 {
    let mut angle = f64::atan2(a.x * b.y - a.y * b.x, a.dot(b));
    if angle <= -PI {
        angle += 2.0 * PI;
    }
    angle
}

/// OCCT gp_Vec::IsOpposite (gp_Vec.hxx L130-134).
fn vec_is_opposite(a: DVec3, b: DVec3, angular_tolerance: f64) -> bool {
    PI - gp_vec_angle(a, b) <= angular_tolerance
}

/// OCCT CalcAngleAT (L702-732) — OCC78: the angle of rotation of the
/// trihedron normal between two positions.
fn calc_angle_at(tangent: DVec3, normal: DVec3, prev_tangent: DVec3, prev_normal: DVec3) -> f64 {
    let angle = gp_vec_angle(tangent, prev_tangent);
    let normal_rot = if angle.abs() > ANGULAR && angle.abs() < PI - ANGULAR {
        let cross = tangent.cross(prev_tangent).normalize();
        normal
            + angle.sin() * cross.cross(normal)
            + (1.0 - angle.cos()) * cross.cross(cross.cross(normal))
    } else {
        normal
    };
    let mut angle_at = gp_vec_angle(normal_rot, prev_normal);
    if angle_at > ANGULAR
        && PI - angle_at > ANGULAR
        && vec_is_opposite(normal_rot.cross(prev_normal), prev_tangent, ANGULAR)
    {
        angle_at = -angle_at;
    }
    angle_at
}

/// OCCT GeomFill_CorrectedFrenet.
#[derive(Clone)]
pub struct CorrectedFrenet {
    pub(crate) base: TrihedronLawBase,
    frenet: Frenet,
    evol_around_t: Option<LawComposite>,
    tlaw: Option<LawFunctionHandle>,
    at: DVec3,
    an: DVec3,
    is_frenet: bool,
    my_for_evaluation: bool,
    harr_poles: Option<Vec<f64>>,
    harr_angle: Option<Vec<f64>>,
    harr_tangent: Option<Vec<DVec3>>,
    harr_normal: Option<Vec<DVec3>>,
}

impl CorrectedFrenet {
    /// OCCT GeomFill_CorrectedFrenet() (L292-297).
    pub fn new() -> Self {
        CorrectedFrenet {
            base: TrihedronLawBase::default(),
            frenet: Frenet::new(),
            evol_around_t: None,
            tlaw: None,
            at: DVec3::ZERO,
            an: DVec3::ZERO,
            is_frenet: false,
            my_for_evaluation: false,
            harr_poles: None,
            harr_angle: None,
            harr_tangent: None,
            harr_normal: None,
        }
    }

    /// OCCT GeomFill_CorrectedFrenet(const bool ForEvaluation) (L299-305).
    pub fn new_for_evaluation(for_evaluation: bool) -> Self {
        let mut law = Self::new();
        law.my_for_evaluation = for_evaluation;
        law
    }

    /// OCCT Init (L354-425) — builds the Normal angle evolution law.
    fn init(&mut self) {
        self.evol_around_t = Some(LawComposite::new());
        let nb_i = TrihedronLaw::nb_intervals(&self.frenet, GeomAbsShape::C0);
        let mut t = vec![0.0f64; nb_i + 1];
        TrihedronLaw::intervals(&self.frenet, &mut t, GeomAbsShape::C0);
        let mut func: Option<LawFunctionHandle> = None;
        // OCC78
        let mut seq_poles: Vec<f64> = Vec::new();
        let mut seq_angle: Vec<f64> = Vec::new();
        let mut seq_tangent: Vec<DVec3> = Vec::new();
        let mut seq_normal: Vec<DVec3> = Vec::new();
        let mut tangent = DVec3::ZERO;
        let mut normal = DVec3::ZERO;
        let mut bn = DVec3::ZERO;
        let trimmed_first;
        let trimmed_last;
        let trimmed_is_periodic;
        {
            let trimmed = self.base.my_trimmed.as_ref().unwrap();
            trimmed_first = curve_first_parameter(trimmed);
            trimmed_last = curve_last_parameter(trimmed);
            trimmed_is_periodic = trimmed.is_periodic();
        }
        TrihedronLaw::d0(
            &self.frenet,
            trimmed_first,
            &mut tangent,
            &mut normal,
            &mut bn,
        );
        let mut start_ang = 0.0f64;
        let av_step = (trimmed_last - trimmed_first) / 10.0;
        // AT/AN are the OCCT members, passed by reference into InitInterval.
        let mut at = DVec3::ZERO;
        let mut an = DVec3::ZERO;
        for i in 1..=nb_i {
            let nb_step = (((t[i] - t[i - 1]) / av_step) as i32).max(3);
            let step = (t[i] - t[i - 1]) / nb_step as f64;
            let ok = self.init_interval(
                t[i - 1],
                t[i],
                step,
                &mut start_ang,
                &mut tangent,
                &mut normal,
                &mut at,
                &mut an,
                &mut func,
                &mut seq_poles,
                &mut seq_angle,
                &mut seq_tangent,
                &mut seq_normal,
            );
            if !ok && self.is_frenet {
                self.is_frenet = false;
            }
            self.evol_around_t
                .as_mut()
                .unwrap()
                .change_laws()
                .push(func.clone().unwrap());
        }
        if trimmed_is_periodic {
            self.evol_around_t.as_mut().unwrap().set_periodic();
        }
        // TLaw = EvolAroundT (handle alias in OCCT; a value copy here).
        self.tlaw = Some(Rc::new(RefCell::new(self.evol_around_t.clone().unwrap())));
        // OCC78
        let i_end = seq_poles.len();
        if i_end != 0 {
            self.harr_poles = Some(seq_poles);
            self.harr_angle = Some(seq_angle);
            self.harr_tangent = Some(seq_tangent);
            self.harr_normal = Some(seq_normal);
        }
        self.at = at;
        self.an = an;
    }

    /// OCCT InitInterval (L432-662) — computes the angle law on one span.
    /// (OCCT declares it const and mutates the shared frenet through the
    /// handle; the rcad port takes &mut self.)
    #[allow(clippy::too_many_arguments)]
    fn init_interval(
        &mut self,
        first: f64,
        last: f64,
        step: f64,
        start_ang: &mut f64,
        prev_tangent: &mut DVec3,
        prev_normal: &mut DVec3,
        a_t: &mut DVec3,
        a_n: &mut DVec3,
        func_int: &mut Option<LawFunctionHandle>,
        seq_poles: &mut Vec<f64>,
        seq_angle: &mut Vec<f64>,
        seq_tangent: &mut Vec<DVec3>,
        seq_normal: &mut Vec<DVec3>,
    ) -> bool {
        let mut tangent = DVec3::ZERO;
        let mut normal = DVec3::ZERO;
        let mut bn = DVec3::ZERO;
        let mut parameters: Vec<f64> = Vec::new();
        let _ = (&mut tangent, &mut normal, &mut bn);
        let mut evol_at: Vec<f64> = Vec::new();
        let mut param = first;
        let mut is_zero = true;
        let mut is_const = true;
        // frenet->SetInterval(First, Last) — to have right evaluation at
        // bounds.
        TrihedronLaw::set_interval(&mut self.frenet, first, last);
        let my_curve = self.base.my_curve.as_ref().unwrap().clone();
        let cs = SnglrFunc::new(my_curve.clone());
        let eval = |u: f64| cs.eval_d0(u);
        let (_box, gap) = curve_box_range_fn(&eval, first, last, 1.0e-2);
        let length_min = gap * 1.0e-4;
        *a_t = DVec3::ZERO;
        *a_n = DVec3::ZERO;
        let mut angle_at = 0.0f64;
        let mut curr_step = step;
        let mut is_planar = false;
        if !self.my_for_evaluation {
            is_planar = find_plane(&my_curve).is_some();
        }
        let mut i = 1usize;
        #[allow(unused_mut)]
        let mut curr_param = param;
        let d_last = last - P_CONFUSION;
        while param < last {
            if curr_param > d_last {
                if (d_last - param).abs() < SQUARE_CONFUSION {
                    param = curr_param;
                }
                curr_step = d_last - param;
                curr_param = last;
            }
            if is_planar {
                curr_param = last;
            }
            TrihedronLaw::d0(&self.frenet, curr_param, &mut tangent, &mut normal, &mut bn);
            if gp_vec_angle(*prev_tangent, tangent) < PI / 3.0 || i == 1 {
                parameters.push(curr_param);
                // OCC78
                seq_poles.push(param);
                seq_angle.push(if i > 1 { evol_at[i - 2] } else { *start_ang });
                seq_tangent.push(*prev_tangent);
                seq_normal.push(*prev_normal);
                angle_at = calc_angle_at(tangent, normal, *prev_tangent, *prev_normal);
                if is_const && i > 1 && angle_at.abs() > P_CONFUSION {
                    is_const = false;
                }
                angle_at += if i > 1 { evol_at[i - 2] } else { *start_ang };
                evol_at.push(angle_at);
                *prev_normal = normal;
                if is_zero && angle_at.abs() > P_CONFUSION {
                    is_zero = false;
                }
                *a_t += tangent;
                let cross = tangent.cross(normal);
                *a_n = angle_at.sin() * cross
                    + (1.0 - angle_at.cos()) * tangent.cross(cross)
                    + (normal + *a_n);
                *prev_tangent = tangent;
                param = curr_param;
                i += 1;
                // Evaluate the Next step
                let (pon_c, d1) = cs.eval_d1(param);
                let l = (pon_c.length() / 2.0).max(length_min);
                let mut norm = d1.length();
                if norm < CONFUSION {
                    norm = CONFUSION;
                }
                curr_step = l / norm;
                if curr_step > step {
                    curr_step = step; // default value
                }
            } else {
                curr_step /= 2.0; // Step too long!
            }
            curr_param = param + curr_step;
        }
        if !is_planar {
            *a_t /= (parameters.len() - 1) as f64;
            *a_n /= (parameters.len() - 1) as f64;
        }
        *start_ang = angle_at;
        // Interpolation
        if is_const || is_planar {
            let mut constant = LawConstant::new();
            constant.set(angle_at, first, last);
            *func_int = Some(Rc::new(RefCell::new(constant)));
        } else {
            let pararr = parameters;
            let angle_at_arr = evol_at;
            let mut law_at = LawInterpolate::with_parameters(
                angle_at_arr.clone(),
                pararr.clone(),
                false,
                P_CONFUSION,
            );
            law_at.perform();
            let bs = law_at.curve();
            smoothlaw(&mut bs.borrow_mut(), &angle_at_arr, &pararr, 0.1);
            *func_int = Some(Rc::new(RefCell::new(LawBSpFunc::with_curve(
                bs, first, last,
            ))));
        }
        is_zero
    }

    /// OCCT GetAngleAT (L736-794) — OCC78.
    fn get_angle_at(&self, param: f64) -> f64 {
        // Search index of low margin from poles of TLaw by bisection method
        let harr_poles = self.harr_poles.as_ref().unwrap();
        let harr_angle = self.harr_angle.as_ref().unwrap();
        let harr_tangent = self.harr_tangent.as_ref().unwrap();
        let harr_normal = self.harr_normal.as_ref().unwrap();
        let mut i_b = 1usize;
        let mut i_e = harr_poles.len();
        let mut i_c = (i_e + i_b) / 2;
        if param == harr_poles[i_b - 1] {
            return self.tlaw.as_ref().unwrap().borrow_mut().value(param);
        }
        if param > harr_poles[i_e - 1] {
            i_c = i_e;
        }
        if i_c < i_e {
            while harr_poles[i_c - 1] > param || param > harr_poles[i_c] {
                if harr_poles[i_c - 1] < param {
                    i_b = i_c;
                } else {
                    i_e = i_c;
                }
                i_c = (i_e + i_b) / 2;
            }
            if harr_poles[i_c - 1] == param || param == harr_poles[i_c] {
                return self.tlaw.as_ref().unwrap().borrow_mut().value(param);
            }
        }
        // Calculate differentiation between approximated and local values
        // of AngleAT
        let ang_p = self.tlaw.as_ref().unwrap().borrow_mut().value(param);
        let ang_po = harr_angle[i_c - 1];
        let d_ang = ang_p - ang_po;
        let mut tangent = DVec3::ZERO;
        let mut normal = DVec3::ZERO;
        let mut bn = DVec3::ZERO;
        TrihedronLaw::d0(&self.frenet, param, &mut tangent, &mut normal, &mut bn);
        let d_ang_local =
            calc_angle_at(tangent, normal, harr_tangent[i_c - 1], harr_normal[i_c - 1]);
        let da = diff_ang(d_ang_local, d_ang);
        // The correction (there is core of OCC78 bug)
        if da.abs() > PI / 2.0 {
            return ang_po + d_ang_local;
        }
        ang_p
    }
}

impl Default for CorrectedFrenet {
    fn default() -> Self {
        Self::new()
    }
}

impl TrihedronLaw for CorrectedFrenet {
    fn my_curve(&self) -> &Option<Curve3> {
        &self.base.my_curve
    }

    fn my_trimmed(&self) -> &Option<Curve3> {
        &self.base.my_trimmed
    }

    fn set_my_curve(&mut self, c: Curve3) {
        self.base.my_curve = Some(c);
    }

    fn set_my_trimmed(&mut self, c: Option<Curve3>) {
        self.base.my_trimmed = c;
    }

    fn copy_law(&self) -> Box<dyn TrihedronLaw> {
        let mut copy = CorrectedFrenet::new();
        if let Some(curve) = &self.base.my_curve {
            TrihedronLaw::set_curve(&mut copy, curve.clone());
        }
        Box::new(copy)
    }

    /// OCCT SetCurve (L332-356).
    #[allow(unconditional_recursion)]
    fn set_curve(&mut self, c: Curve3) -> bool {
        // Base-class initialization (the lint sees the same-name method).
        #[allow(unconditional_recursion)]
        fn base_set_curve(this: &mut CorrectedFrenet, c: Curve3) -> bool {
            TrihedronLaw::set_curve(this, c)
        }
        base_set_curve(self, c.clone());
        // frenet->SetCurve(C)
        TrihedronLaw::set_curve(&mut self.frenet, c.clone());
        let analytic = matches!(
            c,
            Curve3::Line(_)
                | Curve3::Circle(_)
                | Curve3::Ellipse(_)
                | Curve3::Hyperbola(_)
                | Curve3::Parabola(_)
        );
        if analytic {
            // No probleme isFrenet
            self.is_frenet = true;
        } else {
            // We have to search singularities
            self.is_frenet = true;
            self.init();
        }
        self.is_frenet
    }

    /// OCCT D0 (L746-768).
    fn d0(
        &self,
        param: f64,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
    ) -> bool {
        TrihedronLaw::d0(&self.frenet, param, tangent, normal, binormal);
        if self.is_frenet {
            return true;
        }
        // angleAT = TLaw->Value(Param);
        let angle_at = self.get_angle_at(param); // OCC78
        // rotation around Tangent
        let cross = tangent.cross(*normal);
        *normal = angle_at.sin() * cross + (1.0 - angle_at.cos()) * tangent.cross(cross) + *normal;
        *binormal = tangent.cross(*normal);
        true
    }

    /// OCCT D1 (L770-845).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
    ) -> bool {
        TrihedronLaw::d1(
            &self.frenet,
            param,
            tangent,
            dtangent,
            normal,
            dnormal,
            binormal,
            dbinormal,
        );
        if self.is_frenet {
            return true;
        }
        let mut tl_angle = 0.0f64;
        let mut d_angle_at = 0.0f64;
        self.tlaw
            .as_ref()
            .unwrap()
            .borrow_mut()
            .d1(param, &mut tl_angle, &mut d_angle_at);
        let angle_at = self.get_angle_at(param); // OCC78
        let sina = angle_at.sin();
        let cosa = angle_at.cos();
        let cross = tangent.cross(*normal);
        let dcross = dtangent.cross(*normal) + tangent.cross(*dnormal);
        let tcross = tangent.cross(cross);
        let dtcross = dtangent.cross(cross) + tangent.cross(dcross);
        // aux = sina*dcross + cosa*d_angleAT*cross;
        let aux = sina * dcross + cosa * d_angle_at * cross;
        // aux = (1-cosa)*dtcross + sina*d_angleAT*tcross + aux;
        let aux = (1.0 - cosa) * dtcross + sina * d_angle_at * tcross + aux;
        *dnormal += aux;
        *normal = sina * cross + (1.0 - cosa) * tcross + *normal;
        *binormal = tangent.cross(*normal);
        *dbinormal = dtangent.cross(*normal) + tangent.cross(*dnormal);
        true
    }

    /// OCCT D2 (L847-995).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        d2tangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        d2normal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
        d2binormal: &mut DVec3,
    ) -> bool {
        TrihedronLaw::d2(
            &self.frenet,
            param,
            tangent,
            dtangent,
            d2tangent,
            normal,
            dnormal,
            d2normal,
            binormal,
            dbinormal,
            d2binormal,
        );
        if self.is_frenet {
            return true;
        }
        let mut tl_angle = 0.0f64;
        let mut d_angle_at = 0.0f64;
        let mut d2_angle_at = 0.0f64;
        self.tlaw
            .as_ref()
            .unwrap()
            .borrow_mut()
            .d2(param, &mut tl_angle, &mut d_angle_at, &mut d2_angle_at);
        let angle_at = self.get_angle_at(param); // OCC78
        let sina = angle_at.sin();
        let cosa = angle_at.cos();
        let cross = tangent.cross(*normal);
        let dcross = dtangent.cross(*normal) + tangent.cross(*dnormal);
        let d2cross =
            d2tangent.cross(*normal) + 2.0 * dtangent.cross(*dnormal) + tangent.cross(*d2normal);
        let tcross = tangent.cross(cross);
        let dtcross = dtangent.cross(cross) + tangent.cross(dcross);
        let d2tcross = d2tangent.cross(cross) + 2.0 * dtangent.cross(dcross) + tangent.cross(d2cross);
        // aux = sina*d2cross + 2*cosa*d_angleAT*dcross
        //     + (cosa*d2_angleAT - sina*d_angleAT^2)*cross;
        let aux1 = sina * d2cross
            + 2.0 * cosa * d_angle_at * dcross
            + (cosa * d2_angle_at - sina * d_angle_at * d_angle_at) * cross;
        // aux = (1-cosa)*d2tcross + 2*sina*d_angleAT*dtcross
        //     + (cosa*d_angleAT^2 + sina*d2_angleAT)*tcross + aux;
        let aux = (1.0 - cosa) * d2tcross
            + 2.0 * sina * d_angle_at * dtcross
            + (cosa * d_angle_at * d_angle_at + sina * d2_angle_at) * tcross
            + aux1;
        *d2normal += aux;
        // aux = sina*dcross + cosa*d_angleAT*cross;
        // aux = (1-cosa)*dtcross + sina*d_angleAT*tcross + aux;
        let aux = sina * dcross + cosa * d_angle_at * cross;
        let aux = (1.0 - cosa) * dtcross + sina * d_angle_at * tcross + aux;
        *dnormal += aux;
        *normal = sina * cross + (1.0 - cosa) * tcross + *normal;
        *binormal = tangent.cross(*normal);
        *dbinormal = dtangent.cross(*normal) + tangent.cross(*dnormal);
        *d2binormal =
            d2tangent.cross(*normal) + 2.0 * dtangent.cross(*dnormal) + tangent.cross(*d2normal);
        true
    }

    /// OCCT NbIntervals (L997-1028).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let nb_frenet = TrihedronLaw::nb_intervals(&self.frenet, s);
        if self.is_frenet {
            return nb_frenet;
        }
        let nb_law = self
            .evol_around_t
            .as_ref()
            .map(|c| LawFunction::nb_intervals(c, s))
            .unwrap_or(0);
        if nb_frenet == 1 {
            return nb_law;
        }
        let frenet_int = {
            let mut t = vec![0.0f64; nb_frenet + 1];
            TrihedronLaw::intervals(&self.frenet, &mut t, s);
            t
        };
        let law_int = {
            let mut t = vec![0.0f64; nb_law + 1];
            if let Some(evol) = &self.evol_around_t {
                LawFunction::intervals(evol, &mut t, s);
            }
            t
        };
        let mut fusion: Vec<f64> = Vec::new();
        fuse_intervals(&frenet_int, &law_int, &mut fusion, P_CONFUSION, true);
        fusion.len() - 1
    }

    /// OCCT Intervals (L934-965) — note the OCCT fall-through when
    /// NbFrenet == 1 (the law intervals are overwritten by the fusion).
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        if self.is_frenet {
            TrihedronLaw::intervals(&self.frenet, t, s);
            return;
        }
        let nb_frenet = TrihedronLaw::nb_intervals(&self.frenet, s);
        if nb_frenet == 1 {
            if let Some(evol) = &self.evol_around_t {
                LawFunction::intervals(evol, t, s);
            }
        }
        let nb_law = self
            .evol_around_t
            .as_ref()
            .map(|c| LawFunction::nb_intervals(c, s))
            .unwrap_or(0);
        let frenet_int = {
            let mut t = vec![0.0f64; nb_frenet + 1];
            TrihedronLaw::intervals(&self.frenet, &mut t, s);
            t
        };
        let law_int = {
            let mut t = vec![0.0f64; nb_law + 1];
            if let Some(evol) = &self.evol_around_t {
                LawFunction::intervals(evol, &mut t, s);
            }
            t
        };
        let mut fusion: Vec<f64> = Vec::new();
        fuse_intervals(&frenet_int, &law_int, &mut fusion, P_CONFUSION, true);
        for i in 1..=fusion.len() {
            t[i - 1] = fusion[i - 1];
        }
    }

    /// OCCT SetInterval (L967-976).  Calls the base implementation
    /// explicitly (the lint sees the same-name method as recursion).
    #[allow(unconditional_recursion)]
    fn set_interval(&mut self, first: f64, last: f64) {
        TrihedronLaw::set_interval(self, first, last);
        TrihedronLaw::set_interval(&mut self.frenet, first, last);
        if !self.is_frenet {
            if let Some(evol) = &mut self.evol_around_t {
                let trimmed = LawFunction::trim(evol, first, last, P_CONFUSION / 2.0);
                self.tlaw = Some(trimmed);
            }
        }
    }

    /// OCCT GetAverageLaw (L1045-1057).
    fn get_average_law(&self, atangent: &mut DVec3, anormal: &mut DVec3, abinormal: &mut DVec3) {
        if self.is_frenet {
            TrihedronLaw::get_average_law(&self.frenet, atangent, anormal, abinormal);
        } else {
            *atangent = self.at;
            *anormal = self.an;
            *abinormal = *atangent;
            *abinormal = abinormal.cross(*anormal);
        }
    }

    /// OCCT IsConstant (L1059-1062).
    fn is_constant(&self) -> bool {
        matches!(self.base.my_curve, Some(Curve3::Line(_)))
    }

    /// OCCT IsOnlyBy3dCurve (L1064-1067).
    fn is_only_by3d_curve(&self) -> bool {
        true
    }
}

impl CorrectedFrenet {
    /// OCCT EvaluateBestMode (L978-1043).
    pub fn evaluate_best_mode(&self) -> Trihedron {
        let Some(evol) = &self.evol_around_t else {
            return Trihedron::IsFrenet; // Frenet
        };
        let max_angle = 3.0 * PI / 4.0;
        let max_torsion = 100.0;
        let nb_int = LawFunction::nb_intervals(evol, GeomAbsShape::CN);
        let mut int = vec![0.0f64; nb_int + 1];
        LawFunction::intervals(evol, &mut int, GeomAbsShape::CN);
        let mut old = DVec2::ZERO;
        let mut prev_vec = DVec2::ZERO;
        let nb_samples = 10usize;
        let mut k = 1usize;
        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        for i in 1..=nb_int {
            let tmin = int[i - 1];
            let tmax = int[i];
            let torsion = compute_torsion(tmin, trimmed);
            if torsion.abs() > max_torsion {
                return Trihedron::IsDiscreteTrihedron; // DiscreteTrihedron
            }
            let trimmedlaw = LawFunction::trim(evol, tmin, tmax, P_CONFUSION / 2.0);
            let step = (int[i] - int[i - 1]) / nb_samples as f64;
            for j in 0..=nb_samples {
                let u = tmin + j as f64 * step;
                let v = trimmedlaw.borrow_mut().value(u);
                let point2d = DVec2::new(u, v);
                if j != 0 {
                    let a_vec = point2d - old;
                    if k > 2 {
                        let the_angle = gp_vec2d_angle(prev_vec, a_vec);
                        if the_angle.abs() > max_angle {
                            return Trihedron::IsDiscreteTrihedron; // DiscreteTrihedron
                        }
                    }
                    prev_vec = a_vec;
                }
                old = point2d;
                k += 1;
            }
        }
        Trihedron::IsCorrectedFrenet // CorrectedFrenet
    }

}
