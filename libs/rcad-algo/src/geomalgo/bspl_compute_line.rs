// OCCT Approx_BSplComputeLine (TKGeomBase Approx) — 1:1 Rust translation of
// `Approx_BSplComputeLine.gxx` (L139-1458) as instantiated by
// `AppDef_BSplineCompute_0.cxx`:
//   MultiLine            = AppDef_MultiLine   (rcad app_def::MultiLine)
//   LineTool             = AppDef_MyLineTool  (rcad app_def::my_line_tool)
//   Approx_MyBSplGradient = AppDef_MyBSplGradientOfBSplineCompute
//     (= AppParCurves_BSpGradient.gxx, rcad app_par_curves_bsp::BSpGradient)
//   Approx_MyGradientbis  = AppDef_MyGradientbisOfBSplineCompute
//     (= AppParCurves_Gradient.gxx, rcad app_par_curves::Gradient)
//   Approx_BSpParLeastSquareOfMyBSplGradient
//     (= AppParCurves_BSpParLeastSquare.gxx, the BSP constructor family of
//      rcad app_par_curves::LeastSquare)
//
// The struct keeps the AppDef instantiation name; the gxx is shared with the
// Approx_BSplComputeLine and BRepApprox_TheComputeLineOfApprox bindings,
// which will extract the generic engine when their consumers land (the
// int_curve_curve_gen precedent: concrete first, genericize at the second
// instantiation).
//
// Documented deferral: `ComputeCurve` is declared in AppDef_BSplineCompute.hxx
// but has no definition anywhere in OCCT (a stale declaration copied from
// Approx_ComputeLine) — nothing calls it; rcad omits it.

use rcad_kernel::math::math_matrix::Vector as RVector;
use rcad_kernel::math::VecD;

use super::app_def::{my_line_tool, MultiLine};
use super::app_par_curves::{Gradient, LeastSquare};
use super::app_par_curves_bsp::BSpGradient;
use super::approx_int::{
    AppParConstraint, ApproxParamType, ConstraintCouple, MultiBSpCurve, MultiCurve,
};

/// OCCT Standard_Real RealLast() (Standard_Real.hxx).
const REAL_LAST: f64 = f64::MAX;

/// OCCT AppDef_BSplineCompute.
pub struct BSplineCompute {
    /// OCCT TheMultiBSpCurve.
    the_multi_bsp_curve: MultiBSpCurve,
    alldone: bool,
    tolreached: bool,
    /// OCCT Par.
    par: ApproxParamType,
    /// OCCT myParameters (HArray1(firstpt, lastpt)); None until Perform.
    my_parameters: Option<RVector>,
    /// OCCT myfirstParam.
    myfirst_param: Option<RVector>,
    /// OCCT myknots.
    myknots: Option<Vec<f64>>,
    /// OCCT mymults.
    mymults: Option<Vec<i32>>,
    myhasknots: bool,
    myhasmults: bool,
    /// OCCT myConstraints (HArray1(1, 2)).
    my_constraints: [ConstraintCouple; 2],
    mydegremin: i32,
    mydegremax: i32,
    mytol3d: f64,
    mytol2d: f64,
    currenttol3d: f64,
    currenttol2d: f64,
    mycut: bool,
    mysquares: bool,
    myitermax: i32,
    myfirstc: AppParConstraint,
    mylastc: AppParConstraint,
    realfirstc: AppParConstraint,
    reallastc: AppParConstraint,
    mycont: i32,
    mylambda1: f64,
    mylambda2: f64,
    my_periodic: bool,
}

impl BSplineCompute {
    /// OCCT AppDef_BSplineCompute(degmin..Squares) (gxx L573-597) — the
    /// empty constructor with the parametrization.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
        squares: bool,
    ) -> Self {
        BSplineCompute {
            the_multi_bsp_curve: MultiBSpCurve::new_nbpol(1),
            alldone: false,
            tolreached: false,
            par: parametrization,
            my_parameters: None,
            myfirst_param: None,
            myknots: None,
            mymults: None,
            myhasknots: false,
            myhasmults: false,
            my_constraints: [ConstraintCouple { index: 0, constraint: AppParConstraint::NoConstraint }, ConstraintCouple { index: 0, constraint: AppParConstraint::NoConstraint }],
            mydegremin: degreemin,
            mydegremax: degreemax,
            mytol3d: tolerance3d,
            mytol2d: tolerance2d,
            currenttol3d: REAL_LAST,
            currenttol2d: REAL_LAST,
            mycut: cutting,
            mysquares: squares,
            myitermax: nb_iterations,
            myfirstc: AppParConstraint::TangencyPoint,
            mylastc: AppParConstraint::TangencyPoint,
            realfirstc: AppParConstraint::NoConstraint,
            reallastc: AppParConstraint::NoConstraint,
            mycont: -1,
            mylambda1: 0.0,
            mylambda2: 0.0,
            my_periodic: false,
        }
    }

    /// OCCT AppDef_BSplineCompute(Line, degmin..Squares) (gxx L599-623).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_line(
        line: &MultiLine,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
        squares: bool,
    ) -> Self {
        let mut r = BSplineCompute::new(
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            parametrization,
            squares,
        );
        r.perform(line);
        r
    }

    /// OCCT AppDef_BSplineCompute(Line, Parameters, degmin..Squares)
    /// (gxx L529-552) — Par is forced to Approx_IsoParametric.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_line_and_parameters(
        line: &MultiLine,
        parameters: &RVector,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        squares: bool,
    ) -> Self {
        let mut r = BSplineCompute {
            the_multi_bsp_curve: MultiBSpCurve::new_nbpol(1),
            alldone: false,
            tolreached: false,
            par: ApproxParamType::IsoParametric,
            my_parameters: None,
            myfirst_param: Some(parameters.clone()),
            myknots: None,
            mymults: None,
            myhasknots: false,
            myhasmults: false,
            my_constraints: [ConstraintCouple { index: 0, constraint: AppParConstraint::NoConstraint }, ConstraintCouple { index: 0, constraint: AppParConstraint::NoConstraint }],
            mydegremin: degreemin,
            mydegremax: degreemax,
            mytol3d: tolerance3d,
            mytol2d: tolerance2d,
            currenttol3d: REAL_LAST,
            currenttol2d: REAL_LAST,
            mycut: cutting,
            mysquares: squares,
            myitermax: nb_iterations,
            myfirstc: AppParConstraint::TangencyPoint,
            mylastc: AppParConstraint::TangencyPoint,
            realfirstc: AppParConstraint::NoConstraint,
            reallastc: AppParConstraint::NoConstraint,
            mycont: -1,
            mylambda1: 0.0,
            mylambda2: 0.0,
            my_periodic: false,
        };
        // myfirstParam = new HArray1(Parameters.Lower(), Upper()); fill.
        let mut fp = RVector::new(parameters.lower(), parameters.upper());
        for i in parameters.lower()..=parameters.upper() {
            fp.set(i, parameters.get(i));
        }
        r.myfirst_param = Some(fp);
        r.perform(line);
        r
    }

    /// OCCT AppDef_BSplineCompute(Parameters, degmin..Squares) (gxx
    /// L554-571).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_parameters(
        parameters: &RVector,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        squares: bool,
    ) -> Self {
        let mut r = BSplineCompute {
            the_multi_bsp_curve: MultiBSpCurve::new_nbpol(1),
            alldone: false,
            tolreached: false,
            par: ApproxParamType::IsoParametric,
            my_parameters: None,
            myfirst_param: Some(parameters.clone()),
            myknots: None,
            mymults: None,
            myhasknots: false,
            myhasmults: false,
            my_constraints: [ConstraintCouple { index: 0, constraint: AppParConstraint::NoConstraint }, ConstraintCouple { index: 0, constraint: AppParConstraint::NoConstraint }],
            mydegremin: degreemin,
            mydegremax: degreemax,
            mytol3d: tolerance3d,
            mytol2d: tolerance2d,
            currenttol3d: REAL_LAST,
            currenttol2d: REAL_LAST,
            mycut: cutting,
            mysquares: squares,
            myitermax: nb_iterations,
            myfirstc: AppParConstraint::TangencyPoint,
            mylastc: AppParConstraint::TangencyPoint,
            realfirstc: AppParConstraint::NoConstraint,
            reallastc: AppParConstraint::NoConstraint,
            mycont: -1,
            mylambda1: 0.0,
            mylambda2: 0.0,
            my_periodic: false,
        };
        r
    }

    /// OCCT FirstTangencyVector(Line, index, V) (gxx L139-283).
    fn first_tangency_vector(&self, line: &MultiLine, index: i32, v: &mut RVector) {
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_v = vec![glam::DVec3::ZERO; mynb_p3d as usize];
        let mut tab_v2d = vec![glam::DVec2::ZERO; mynb_p2d as usize];

        let ok = if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::tangency_3d_2d(line, index as usize, &mut tab_v, &mut tab_v2d)
        } else if nb_p2d != 0 {
            my_line_tool::tangency_2d(line, index as usize, &mut tab_v2d)
        } else {
            my_line_tool::tangency_3d(line, index as usize, &mut tab_v)
        };

        if ok {
            if nb_p3d != 0 {
                let mut j = 1;
                for t in tab_v.iter().take(nb_p3d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    v.set(j + 2, t.z);
                    j += 3;
                }
            }
            if nb_p2d != 0 {
                let mut j = nb_p3d * 3 + 1;
                for t in tab_v2d.iter().take(nb_p2d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    j += 2;
                }
            }
        } else {
            // Search for a tangent vector by construction of a parabola:
            let first_c = AppParConstraint::PassPoint;
            let last_c = AppParConstraint::PassPoint;
            let nbpoles = 3;
            let mut mypar = RVector::new(index, index + 2);
            self.parameters_compute(line, index, index + 2, &mut mypar);
            let mut lsq = LeastSquare::new(
                line,
                index,
                index + 2,
                first_c,
                last_c,
                &rvector_to_vecd(&mypar),
                nbpoles,
            );
            let c = lsq.bezier_value();

            let mut j = 1;
            for i in 1..=nb_p3d {
                let (_p, v1) = c.d1(i as usize, 0.0);
                v.set(j, v1.x);
                v.set(j + 1, v1.y);
                v.set(j + 2, v1.z);
                j += 3;
            }
            let mut j = nb_p3d * 3 + 1;
            for i in (nb_p3d + 1)..=(nb_p3d + nb_p2d) {
                let (_p2d, v2d) = c.d1(i as usize, 0.0);
                v.set(j, v2d.x);
                v.set(j + 1, v2d.y);
                j += 2;
            }
        }
    }

    /// OCCT LastTangencyVector(Line, index, V) (gxx L286-390) — same shape
    /// with the D1 evaluated at parameter 1.0 on the last parabola segment.
    fn last_tangency_vector(&self, line: &MultiLine, index: i32, v: &mut RVector) {
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_v = vec![glam::DVec3::ZERO; mynb_p3d as usize];
        let mut tab_v2d = vec![glam::DVec2::ZERO; mynb_p2d as usize];

        let ok = if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::tangency_3d_2d(line, index as usize, &mut tab_v, &mut tab_v2d)
        } else if nb_p2d != 0 {
            my_line_tool::tangency_2d(line, index as usize, &mut tab_v2d)
        } else {
            my_line_tool::tangency_3d(line, index as usize, &mut tab_v)
        };

        if ok {
            if nb_p3d != 0 {
                let mut j = 1;
                for t in tab_v.iter().take(nb_p3d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    v.set(j + 2, t.z);
                    j += 3;
                }
            }
            if nb_p2d != 0 {
                let mut j = nb_p3d * 3 + 1;
                for t in tab_v2d.iter().take(nb_p2d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    j += 2;
                }
            }
        } else {
            // Search for a tangent vector by construction of a parabola:
            let first_c = AppParConstraint::PassPoint;
            let last_c = AppParConstraint::PassPoint;
            let nbpoles = 3;
            let mut mypar = RVector::new(index - 2, index);
            self.parameters_compute(line, index - 2, index, &mut mypar);
            let mut lsq = LeastSquare::new(
                line,
                index - 2,
                index,
                first_c,
                last_c,
                &rvector_to_vecd(&mypar),
                nbpoles,
            );
            let c = lsq.bezier_value();

            let mut j = 1;
            for i in 1..=nb_p3d {
                let (p, v1) = c.d1(i as usize, 1.0);
                let _ = p;
                v.set(j, v1.x);
                v.set(j + 1, v1.y);
                v.set(j + 2, v1.z);
                j += 3;
            }
            let mut j = nb_p3d * 3 + 1;
            for i in (nb_p3d + 1)..=(nb_p3d + nb_p2d) {
                let (p2d, v2d) = c.d1(i as usize, 1.0);
                v.set(j, v2d.x);
                v.set(j + 1, v2d.y);
                j += 2;
            }
        }
    }

    /// OCCT SearchFirstLambda(Line, aPar, Theknots, V, index) (gxx
    /// L392-461).
    fn search_first_lambda(
        &self,
        line: &MultiLine,
        a_par: &RVector,
        theknots: &[f64],
        v: &RVector,
        index: i32,
    ) -> f64 {
        // dq/dw = lambda* V = (p2-p1)/(u2-u1)
        let nb_p2d = my_line_tool::nb_p2d(line);
        let nb_p3d = my_line_tool::nb_p3d(line);
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_p1 = vec![glam::DVec3::ZERO; mynb_p3d];
        let mut tab_p2 = vec![glam::DVec3::ZERO; mynb_p3d];
        let mut tab_p12d = vec![glam::DVec2::ZERO; mynb_p2d];
        let mut tab_p22d = vec![glam::DVec2::ZERO; mynb_p2d];

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, index as usize, &mut tab_p1, &mut tab_p12d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, index as usize, &mut tab_p12d);
        } else {
            my_line_tool::value_3d(line, index as usize, &mut tab_p1);
        }

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, (index + 1) as usize, &mut tab_p2, &mut tab_p22d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, (index + 1) as usize, &mut tab_p22d);
        } else {
            my_line_tool::value_3d(line, (index + 1) as usize, &mut tab_p2);
        }

        let u1 = a_par.get(index);
        let u2 = a_par.get(index + 1);
        let low = 1; // OCCT V.Lower().
        let nbknots = theknots.len() as i32;

        let (lambda, s);
        if nb_p3d != 0 {
            let p1 = tab_p1[0];
            let p2 = tab_p2[0];
            let p1p2 = p2 - p1;
            let my_v = glam::DVec3::new(v.get(low), v.get(low + 1), v.get(low + 2));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        } else {
            let p12d = tab_p12d[0];
            let p22d = tab_p22d[0];
            let p1p2 = p22d - p12d;
            let my_v = glam::DVec2::new(v.get(low), v.get(low + 1));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        }
        (s * lambda) * (theknots[1] - theknots[0]) / (theknots[(nbknots - 1) as usize] - theknots[0])
    }

    /// OCCT SearchLastLambda(Line, aPar, Theknots, V, index) (gxx
    /// L464-527).
    fn search_last_lambda(
        &self,
        line: &MultiLine,
        a_par: &RVector,
        theknots: &[f64],
        v: &RVector,
        index: i32,
    ) -> f64 {
        // dq/dw = lambda* V = (p2-p1)/(u2-u1)
        let nb_p2d = my_line_tool::nb_p2d(line);
        let nb_p3d = my_line_tool::nb_p3d(line);
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_p = vec![glam::DVec3::ZERO; mynb_p3d];
        let mut tab_p2 = vec![glam::DVec3::ZERO; mynb_p3d];
        let mut tab_p2d = vec![glam::DVec2::ZERO; mynb_p2d];
        let mut tab_p22d = vec![glam::DVec2::ZERO; mynb_p2d];

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, (index - 1) as usize, &mut tab_p, &mut tab_p2d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, (index - 1) as usize, &mut tab_p2d);
        } else {
            my_line_tool::value_3d(line, (index - 1) as usize, &mut tab_p);
        }

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, index as usize, &mut tab_p2, &mut tab_p22d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, index as usize, &mut tab_p22d);
        } else {
            my_line_tool::value_3d(line, index as usize, &mut tab_p2);
        }

        let u1 = a_par.get(index - 1);
        let u2 = a_par.get(index);
        let low = 1;
        let nbknots = theknots.len() as i32;

        let (lambda, s);
        if nb_p3d != 0 {
            let p1 = tab_p[0];
            let p2 = tab_p2[0];
            let p1p2 = p2 - p1;
            let my_v = glam::DVec3::new(v.get(low), v.get(low + 1), v.get(low + 2));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        } else {
            let p12d = tab_p2d[0];
            let p22d = tab_p22d[0];
            let p1p2 = p22d - p12d;
            let my_v = glam::DVec2::new(v.get(low), v.get(low + 1));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        }

        (s * lambda) * (theknots[(nbknots - 1) as usize] - theknots[(nbknots - 2) as usize])
            / (theknots[(nbknots - 1) as usize] - theknots[0])
    }

    /// OCCT Perform(Line) (gxx L625-780).
    pub fn perform(&mut self, line: &MultiLine) {
        let mut finish = false;
        let mut begin = true;

        // Search for the actual constraints given by the Line:
        self.find_real_constraints(line);

        let thefirstpt = my_line_tool::first_point(line) as i32;
        let thelastpt = my_line_tool::last_point(line) as i32;
        let myfirstpt = thefirstpt;
        let mylastpt = thelastpt;

        self.my_constraints[0] = ConstraintCouple { index: myfirstpt, constraint: self.realfirstc };
        self.my_constraints[1] = ConstraintCouple { index: mylastpt, constraint: self.reallastc };

        let mut the_param = RVector::new_init(thefirstpt, thelastpt, 0.0);
        match &self.myfirst_param {
            None => {
                self.parameters_compute(line, thefirstpt, thelastpt, &mut the_param);
            }
            Some(fp) => {
                for i in fp.lower()..=fp.upper() {
                    the_param.set(i + thefirstpt - 1, fp.get(i));
                }
            }
        }

        let mut my_parameters = RVector::new(the_param.lower(), the_param.upper());
        for i in the_param.lower()..=the_param.upper() {
            my_parameters.set(i, the_param.get(i));
        }
        self.my_parameters = Some(my_parameters);
        let mut nbknots = 2i32;
        self.alldone = false;

        if !self.mycut {
            // Case where no additional knots are desired.
            // ==================================================
            if !self.myhasknots {
                let theknots = [0.0f64, 1.0];
                let mut themults = [0i32; 2];
                let a_params = rvector_to_vecd(&the_param);
                self.alldone = self.compute(
                    line,
                    myfirstpt,
                    mylastpt,
                    &mut the_param,
                    &theknots,
                    &mut themults,
                );
            } else if !self.myhasmults {
                let knots = self.myknots.clone().unwrap();
                let mut themults = vec![0i32; knots.len()];
                let a_params = rvector_to_vecd(&the_param);
                self.alldone = self.compute(
                    line,
                    myfirstpt,
                    mylastpt,
                    &mut the_param,
                    &knots,
                    &mut themults,
                );
            } else {
                let knots = self.myknots.clone().unwrap();
                let mut mults = self.mymults.clone().unwrap();
                let a_params = rvector_to_vecd(&the_param);
                self.alldone =
                    self.compute(line, myfirstpt, mylastpt, &mut the_param, &knots, &mut mults);
            }
        } else {
            // cas ou on va iterer a partir de noeuds donnes par
            // l''utilisateur ou a partir d''une bezier.
            // ================================================================
            while !finish {
                self.currenttol3d = REAL_LAST;
                self.currenttol2d = REAL_LAST;

                if self.myhasknots && begin {
                    if !self.myhasmults {
                        // 1st case: the user provides starting knots but
                        // it is up to us to set the multiplicities according to the
                        // desired continuity.
                        // ========================================================
                        let knots = self.myknots.clone().unwrap();
                        let mut themults = vec![0i32; knots.len()];
                        let a_params = rvector_to_vecd(&the_param);
                        self.alldone = self.compute(
                            line,
                            myfirstpt,
                            mylastpt,
                            &mut the_param,
                            &knots,
                            &mut themults,
                        );
                    } else {
                        // 2nd case: the user provides starting knots
                        // with their multiplicities.
                        // ===================================================
                        let knots = self.myknots.clone().unwrap();
                        let mut mults = self.mymults.clone().unwrap();
                        let a_params = rvector_to_vecd(&the_param);
                        self.alldone = self.compute(
                            line,
                            myfirstpt,
                            mylastpt,
                            &mut the_param,
                            &knots,
                            &mut mults,
                        );
                    }
                    begin = false;
                } else {
                    // 3eme cas: l''utilisateur ne donne aucun noeuds de depart
                    // ========================================================
                    let mut theknots = vec![0.0f64; nbknots as usize];
                    let mut themults = vec![0i32; nbknots as usize];
                    theknots[0] = 0.0;
                    theknots[(nbknots - 1) as usize] = 1.0;
                    for i in 2..=(nbknots - 1) {
                        let l = (mylastpt - myfirstpt) as f64 * (i - 1) as f64
                            / (nbknots - 1) as f64;
                        let ll = l as i32;
                        let a = l - ll as f64;
                        let p1 = the_param.get(ll + myfirstpt);
                        let p2 = the_param.get(ll + 1 + myfirstpt);
                        theknots[(i - 1) as usize] = (1.0 - a) * p1 + a * p2;
                    }

                    let a_params = rvector_to_vecd(&the_param);
                    self.alldone = self.compute(
                        line,
                        myfirstpt,
                        mylastpt,
                        &mut the_param,
                        &theknots,
                        &mut themults,
                    );
                }

                if !self.alldone {
                    nbknots += 1;
                } else {
                    finish = true;
                }
            }
        }
    }

    /// OCCT Parameters() (gxx L785-791) — the parameters accessor.
    pub fn parameters(&self) -> Option<&RVector> {
        self.my_parameters.as_ref()
    }

    /// OCCT Value() (gxx L797-800).
    pub fn value(&self) -> &MultiBSpCurve {
        &self.the_multi_bsp_curve
    }

    /// OCCT ChangeValue() (gxx L805-808).
    pub fn change_value(&mut self) -> &mut MultiBSpCurve {
        &mut self.the_multi_bsp_curve
    }

    /// OCCT Parameters(Line, firstP, LastP, TheParameters) (gxx L811-865).
    fn parameters_compute(
        &self,
        line: &MultiLine,
        first_p: i32,
        last_p: i32,
        the_parameters: &mut RVector,
    ) {
        let a_nbp = last_p - first_p + 1;

        // The first parameter should always be zero according to all the logic
        // below, so division by any value will give zero anyway, so it should
        // never be scaled to avoid case when there is only one parameter in
        // the array thus division by zero happens.
        the_parameters.set(first_p, 0.0);
        if a_nbp == 2 {
            the_parameters.set(last_p, 1.0);
        } else if self.par == ApproxParamType::ChordLength
            || self.par == ApproxParamType::Centripetal
        {
            let nb_p3d = my_line_tool::nb_p3d(line) as i32;
            let nb_p2d = my_line_tool::nb_p2d(line) as i32;
            let mut mynb_p3d = nb_p3d;
            let mut mynb_p2d = nb_p2d;
            if nb_p3d == 0 {
                mynb_p3d = 1;
            }
            if nb_p2d == 0 {
                mynb_p2d = 1;
            }

            let mut tab_p = vec![glam::DVec3::ZERO; mynb_p3d as usize];
            let mut tab_pp = vec![glam::DVec3::ZERO; mynb_p3d as usize];
            let mut tab_p2d = vec![glam::DVec2::ZERO; mynb_p2d as usize];
            let mut tab_pp2d = vec![glam::DVec2::ZERO; mynb_p2d as usize];

            for i in (first_p + 1)..=last_p {
                if nb_p3d != 0 && nb_p2d != 0 {
                    my_line_tool::value_3d_2d(line, (i - 1) as usize, &mut tab_p, &mut tab_p2d);
                } else if nb_p2d != 0 {
                    my_line_tool::value_2d(line, (i - 1) as usize, &mut tab_p2d);
                } else {
                    my_line_tool::value_3d(line, (i - 1) as usize, &mut tab_p);
                }

                if nb_p3d != 0 && nb_p2d != 0 {
                    my_line_tool::value_3d_2d(line, i as usize, &mut tab_pp, &mut tab_pp2d);
                } else if nb_p2d != 0 {
                    my_line_tool::value_2d(line, i as usize, &mut tab_pp2d);
                } else {
                    my_line_tool::value_3d(line, i as usize, &mut tab_pp);
                }
                let mut dist = 0.0;
                for j in 1..=nb_p3d {
                    let a_p1 = tab_p[(j - 1) as usize];
                    let a_p2 = tab_pp[(j - 1) as usize];
                    dist += a_p2.distance_squared(a_p1);
                }
                for j in 1..=nb_p2d {
                    let a_p12d = tab_p2d[(j - 1) as usize];
                    let a_p22d = tab_pp2d[(j - 1) as usize];
                    dist += a_p22d.distance_squared(a_p12d);
                }

                dist = dist.sqrt();
                if self.par == ApproxParamType::ChordLength {
                    the_parameters.set(i, the_parameters.get(i - 1) + dist);
                } else {
                    // Par == Approx_Centripetal
                    the_parameters.set(i, the_parameters.get(i - 1) + dist.sqrt());
                }
            }
            for i in (first_p + 1)..=last_p {
                the_parameters.set(i, the_parameters.get(i) / the_parameters.get(last_p));
            }
        } else {
            for i in (first_p + 1)..=last_p {
                the_parameters
                    .set(i, (i as f64 - first_p as f64) / (last_p as f64 - first_p as f64));
            }
        }
    }

    /// OCCT Compute(Line, fpt, lpt, Para, Knots, Mults) (gxx L868-1085).
    #[allow(clippy::too_many_arguments)]
    fn compute(
        &mut self,
        line: &MultiLine,
        fpt: i32,
        lpt: i32,
        para: &mut RVector,
        knots: &[f64],
        mults: &mut [i32],
    ) -> bool {
        let mut mydone;
        let mut the_tol3d = 0.0;
        let mut the_tol2d = 0.0;
        let nbp = lpt - fpt + 1;
        self.mylambda1 = 0.0;
        self.mylambda2 = 0.0;

        let mut a_params = rvector_to_vecd(para);

        for deg in self.mydegremin..=self.mydegremax {
            a_params = rvector_to_vecd(para);

            let nbpoles;
            if !self.myhasmults {
                // OCCT Mults.Lower() = 0 (the rcad slices are 0-based).
                let low = 0usize;
                let up = mults.len() - 1;
                mults[low] = deg + 1;
                mults[up] = deg + 1;
                let mut n = deg + 1;
                let multinter = if self.mycont == -1 {
                    1
                } else {
                    (deg - self.mycont).max(1)
                };
                for m in mults.iter_mut().take(up).skip(low + 1) {
                    *m = multinter;
                    n += multinter;
                }
                nbpoles = n;
            } else {
                let mut n = -deg - 1;
                for m in mults.iter() {
                    n += *m;
                }
                nbpoles = n;
            }

            let mut nbpolestocompare = nbpoles;
            if self.realfirstc == AppParConstraint::TangencyPoint {
                nbpolestocompare += 1;
            }
            if self.reallastc == AppParConstraint::TangencyPoint {
                nbpolestocompare += 1;
            }
            if self.realfirstc == AppParConstraint::CurvaturePoint {
                nbpolestocompare += 1;
            }
            if self.reallastc == AppParConstraint::CurvaturePoint {
                nbpolestocompare += 1;
            }
            if nbpolestocompare > nbp {
                self.interpol(line);
                self.tolreached = true;
                return true;
            }

            let mut my_scu = MultiBSpCurve::new_nbpol(nbpoles as usize);

            if self.mysquares {
                let mut sq = LeastSquare::new_bsp(
                    line,
                    knots,
                    mults,
                    fpt,
                    lpt,
                    self.realfirstc,
                    self.reallastc,
                    &a_params,
                    nbpoles,
                );
                mydone = sq.is_done();
                if mydone {
                    my_scu = sq.bspline_value().clone();
                    let mut fv = 0.0;
                    sq.error(&mut fv, &mut the_tol3d, &mut the_tol2d);
                } else {
                    continue;
                }
            } else if nbpoles != deg + 1 {
                if deg == self.mydegremin
                    && (self.realfirstc >= AppParConstraint::TangencyPoint
                        || self.reallastc >= AppParConstraint::TangencyPoint)
                {
                    let thefitt = LeastSquare::new_bsp(
                        line,
                        knots,
                        mults,
                        fpt,
                        lpt,
                        self.realfirstc,
                        self.reallastc,
                        &a_params,
                        nbpoles,
                    );
                    self.mylambda1 = thefitt.first_lambda() * deg as f64;
                    self.mylambda2 = thefitt.last_lambda() * deg as f64;
                }
                let l1 = self.mylambda1 / deg as f64;
                let l2 = self.mylambda2 / deg as f64;

                let mut grad = BSpGradient::new_with_lambdas(
                    line,
                    fpt,
                    lpt,
                    &self.my_constraints,
                    &mut a_params,
                    knots,
                    mults,
                    deg,
                    self.mytol3d,
                    self.mytol2d,
                    self.myitermax,
                    l1,
                    l2,
                );

                mydone = grad.is_done();
                if mydone {
                    my_scu = grad.value();
                    the_tol3d = grad.max_error_3d();
                    the_tol2d = grad.max_error_2d();
                } else {
                    continue;
                }
            } else {
                let mut grad2 = Gradient::new(
                    line,
                    fpt,
                    lpt,
                    &self.my_constraints,
                    &mut a_params,
                    deg,
                    self.mytol3d,
                    self.mytol2d,
                    self.myitermax,
                );
                mydone = grad2.is_done();
                if mydone {
                    if grad2.value().nb_curves() == 0 {
                        continue;
                    }
                    my_scu = MultiBSpCurve::from_bezier(
                        &grad2.value(),
                        knots.to_vec(),
                        mults.iter().map(|m| *m as usize).collect(),
                    );
                    the_tol3d = grad2.max_error_3d();
                    the_tol2d = grad2.max_error_2d();
                } else {
                    continue;
                }
            }
            let mut save = true;

            for i in 1..=(a_params.len() as i32) {
                let v = a_params.get(i as usize);
                if v <= -0.000001 || v >= 1.000001 {
                    save = false;
                    break;
                }
            }

            if mydone && the_tol3d <= self.mytol3d && the_tol2d <= self.mytol2d {
                // Stockage de la multicurve approximee.
                self.tolreached = true;
                self.the_multi_bsp_curve = my_scu;
                self.currenttol3d = the_tol3d;
                self.currenttol2d = the_tol2d;
                if save {
                    if let Some(mp) = self.my_parameters.as_mut() {
                        for i in mp.lower()..=mp.upper() {
                            mp.set(i, a_params.get((i - mp.lower() + 1) as usize));
                        }
                    }
                }
                return true;
            }

            if the_tol3d <= self.currenttol3d && the_tol2d <= self.currenttol2d {
                self.the_multi_bsp_curve = my_scu;
                self.currenttol3d = the_tol3d;
                self.currenttol2d = the_tol2d;
                if save {
                    if let Some(mp) = self.my_parameters.as_mut() {
                        for i in mp.lower()..=mp.upper() {
                            mp.set(i, a_params.get((i - mp.lower() + 1) as usize));
                        }
                    }
                }
            }
        }

        false
    }

    /// OCCT SetParameters(ThePar) (gxx L1088-1097).
    pub fn set_parameters(&mut self, the_par: &RVector) {
        let mut fp = RVector::new(the_par.lower(), the_par.upper());
        for i in the_par.lower()..=the_par.upper() {
            fp.set(i, the_par.get(i));
        }
        self.myfirst_param = Some(fp);
    }

    /// OCCT SetKnots(Knots) (gxx L1099-1108).
    pub fn set_knots(&mut self, knots: &[f64]) {
        self.myhasknots = true;
        self.myknots = Some(knots.to_vec());
    }

    /// OCCT SetKnotsAndMultiplicities(Knots, Mults) (gxx L1110-1130).
    pub fn set_knots_and_multiplicities(&mut self, knots: &[f64], mults: &[i32]) {
        self.myhasknots = true;
        self.myhasmults = true;
        self.myknots = Some(knots.to_vec());
        self.mymults = Some(mults.to_vec());
    }

    /// OCCT Init(degreemin..Squares) (gxx L1132-1145).
    #[allow(clippy::too_many_arguments)]
    pub fn init(
        &mut self,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
        squares: bool,
    ) {
        self.mydegremin = degreemin;
        self.mydegremax = degreemax;
        self.mytol3d = tolerance3d;
        self.mytol2d = tolerance2d;
        self.par = parametrization;
        self.mysquares = squares;
        self.mycut = cutting;
        self.myitermax = nb_iterations;
    }

    /// OCCT SetDegrees(degreemin, degreemax) (gxx L1147-1152).
    pub fn set_degrees(&mut self, degreemin: i32, degreemax: i32) {
        self.mydegremin = degreemin;
        self.mydegremax = degreemax;
    }

    /// OCCT SetTolerances(Tolerance3d, Tolerance2d) (gxx L1154-1159).
    pub fn set_tolerances(&mut self, tolerance3d: f64, tolerance2d: f64) {
        self.mytol3d = tolerance3d;
        self.mytol2d = tolerance2d;
    }

    /// OCCT SetConstraints(FirstC, LastC) (gxx L1161-1166).
    pub fn set_constraints(&mut self, first_c: AppParConstraint, last_c: AppParConstraint) {
        self.myfirstc = first_c;
        self.mylastc = last_c;
    }

    /// OCCT SetPeriodic(thePeriodic) (gxx L1168-1171).
    pub fn set_periodic(&mut self, the_periodic: bool) {
        self.my_periodic = the_periodic;
    }

    /// OCCT IsAllApproximated() (gxx L1173-1176).
    pub fn is_all_approximated(&self) -> bool {
        self.alldone
    }

    /// OCCT IsToleranceReached() (gxx L1178-1181).
    pub fn is_tolerance_reached(&self) -> bool {
        self.tolreached
    }

    /// OCCT Error(tol3d, tol2d) (gxx L1183-1188).
    pub fn error(&self) -> (f64, f64) {
        (self.currenttol3d, self.currenttol2d)
    }

    /// OCCT SetContinuity(C) (gxx L1190-1193).
    pub fn set_continuity(&mut self, c: i32) {
        self.mycont = c;
    }

    /// OCCT FindRealConstraints(Line) (gxx L1226-1317).
    fn find_real_constraints(&mut self, line: &MultiLine) {
        self.realfirstc = self.myfirstc;
        self.reallastc = self.mylastc;
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let mut tab_v = vec![glam::DVec3::ZERO; nb_p3d.max(1) as usize];
        let mut tab_v2d = vec![glam::DVec2::ZERO; nb_p2d.max(1) as usize];
        let thefirstpt = my_line_tool::first_point(line) as i32;
        let thelastpt = my_line_tool::last_point(line) as i32;

        if self.myfirstc >= AppParConstraint::TangencyPoint {
            let ok = if nb_p3d != 0 && nb_p2d != 0 {
                my_line_tool::tangency_3d_2d(line, thefirstpt as usize, &mut tab_v, &mut tab_v2d)
            } else if nb_p2d != 0 {
                my_line_tool::tangency_2d(line, thefirstpt as usize, &mut tab_v2d)
            } else {
                my_line_tool::tangency_3d(line, thefirstpt as usize, &mut tab_v)
            };

            self.realfirstc = AppParConstraint::PassPoint;
            if ok {
                self.realfirstc = AppParConstraint::TangencyPoint;
                if self.myfirstc == AppParConstraint::CurvaturePoint {
                    let ok2 = if nb_p3d != 0 && nb_p2d != 0 {
                        my_line_tool::tangency_3d_2d(
                            line,
                            thefirstpt as usize,
                            &mut tab_v,
                            &mut tab_v2d,
                        )
                    } else if nb_p2d != 0 {
                        my_line_tool::tangency_2d(line, thefirstpt as usize, &mut tab_v2d)
                    } else {
                        my_line_tool::tangency_3d(line, thefirstpt as usize, &mut tab_v)
                    };
                    if ok2 {
                        self.realfirstc = AppParConstraint::CurvaturePoint;
                    }
                }
            }
        }

        if self.mylastc >= AppParConstraint::TangencyPoint {
            let ok = if nb_p3d != 0 && nb_p2d != 0 {
                my_line_tool::tangency_3d_2d(line, thelastpt as usize, &mut tab_v, &mut tab_v2d)
            } else if nb_p2d != 0 {
                my_line_tool::tangency_2d(line, thelastpt as usize, &mut tab_v2d)
            } else {
                my_line_tool::tangency_3d(line, thelastpt as usize, &mut tab_v)
            };

            self.reallastc = AppParConstraint::PassPoint;
            if ok {
                self.reallastc = AppParConstraint::TangencyPoint;
                if self.mylastc == AppParConstraint::CurvaturePoint {
                    let ok2 = if nb_p3d != 0 && nb_p2d != 0 {
                        my_line_tool::tangency_3d_2d(
                            line,
                            thelastpt as usize,
                            &mut tab_v,
                            &mut tab_v2d,
                        )
                    } else if nb_p2d != 0 {
                        my_line_tool::tangency_2d(line, thelastpt as usize, &mut tab_v2d)
                    } else {
                        my_line_tool::tangency_3d(line, thelastpt as usize, &mut tab_v)
                    };
                    if ok2 {
                        self.reallastc = AppParConstraint::CurvaturePoint;
                    }
                }
            }
        }
    }

    /// OCCT Interpol(Line) (gxx L1319-1429) — the C2 degree-3 interpolation.
    pub fn interpol(&mut self, line: &MultiLine) {
        let deg = 3i32;
        self.mycont = 2;
        let thefirstpt = my_line_tool::first_point(line) as i32;
        let thelastpt = my_line_tool::last_point(line) as i32;
        let mut the_param = RVector::new_init(thefirstpt, thelastpt, 0.0);
        // Par = Approx_ChordLength; (OCCT comment)
        match &self.myfirst_param {
            None => {
                self.parameters_compute(line, thefirstpt, thelastpt, &mut the_param);
            }
            Some(fp) => {
                for i in fp.lower()..=fp.upper() {
                    the_param.set(i + thefirstpt - 1, fp.get(i));
                }
            }
        }
        let mut cons = AppParConstraint::TangencyPoint;

        // Recherche du nombre de noeuds.
        let nbpoints = thelastpt - thefirstpt + 1;

        if nbpoints == 2 {
            cons = AppParConstraint::NoConstraint;
            let mydeg = 1;
            let mut lsq = LeastSquare::new(
                line,
                thefirstpt,
                thelastpt,
                cons,
                cons,
                &rvector_to_vecd(&the_param),
                mydeg + 1,
            );
            self.alldone = lsq.is_done();
            let the_knots = [the_param.get(thefirstpt), the_param.get(thelastpt)];
            let the_mults = [mydeg + 1, mydeg + 1];
            self.the_multi_bsp_curve = MultiBSpCurve::from_bezier(
                &lsq.bezier_value(),
                the_knots.to_vec(),
                the_mults.iter().map(|m| *m as usize).collect(),
            );
            let mut fv = 0.0;
            let mut e3 = 0.0;
            let mut e2 = 0.0;
            let mut lsq = lsq;
            lsq.error(&mut fv, &mut e3, &mut e2);
            self.currenttol3d = e3;
            self.currenttol2d = e2;
        } else {
            let nbpoles = nbpoints + 2;
            let nbknots = nbpoints;

            // Resolution:
            let mut theknots = vec![0.0f64; nbknots as usize];
            theknots[0] = the_param.get(thefirstpt);
            theknots[(nbknots - 1) as usize] = the_param.get(thelastpt);
            let mut themults = vec![0i32; nbknots as usize];
            themults[0] = deg + 1;
            themults[(nbknots - 1) as usize] = deg + 1;

            for i in 2..=(nbknots - 1) {
                theknots[(i - 1) as usize] = the_param.get(i + thefirstpt - 1);
                themults[(i - 1) as usize] = 1;
            }

            let nb_p = (3 * my_line_tool::nb_p3d(line) + 2 * my_line_tool::nb_p2d(line)) as i32;
            let mut v1 = RVector::new_init(1, nb_p, 0.0);
            let mut v2 = RVector::new_init(1, nb_p, 0.0);

            let mut lambda1;
            let mut lambda2;
            if nbpoints == 3 || nbpoints == 4 {
                self.first_tangency_vector(line, thefirstpt, &mut v1);
                lambda1 =
                    self.search_first_lambda(line, &the_param, &theknots, &v1, thefirstpt);

                self.last_tangency_vector(line, thelastpt, &mut v2);
                lambda2 = self.search_last_lambda(line, &the_param, &theknots, &v2, thelastpt);

                lambda1 /= deg as f64;
                lambda2 /= deg as f64;
            } else {
                let nnp = nbpoints.min(9);
                let nnpol = nnp;
                let lastp = thelastpt.min(thefirstpt + nnp - 1);
                let mut sq1 = LeastSquare::new_no_params(
                    line,
                    thefirstpt,
                    lastp,
                    cons,
                    cons,
                    nnpol,
                );

                let mut p1 = VecD::new((lastp - thefirstpt + 1) as usize);
                for i in thefirstpt..=lastp {
                    p1.set((i - thefirstpt + 1) as usize, the_param.get(i));
                }
                sq1.perform(&p1);
                let c1 = sq1.bezier_value();
                let u = 0.0;
                self.tangency_vector(line, &c1, u, &mut v1);

                let firstp = std::cmp::max(thefirstpt, thelastpt - nnp + 1);

                if firstp == thefirstpt && lastp == thelastpt {
                    let u = 1.0;
                    self.tangency_vector(line, &c1, u, &mut v2);
                } else {
                    let mut sq2 = LeastSquare::new_no_params(
                        line,
                        firstp,
                        thelastpt,
                        cons,
                        cons,
                        nnpol,
                    );

                    let mut p2 = VecD::new((thelastpt - firstp + 1) as usize);
                    for i in firstp..=thelastpt {
                        p2.set((i - firstp + 1) as usize, the_param.get(i));
                    }
                    sq2.perform(&p2);
                    let c2 = sq2.bezier_value();
                    let u = 1.0;
                    self.tangency_vector(line, &c2, u, &mut v2);
                }

                lambda1 = 1.0 / deg as f64;
                lambda1 *= (theknots[1] - theknots[0])
                    / (theknots[(nbknots - 1) as usize] - theknots[0]);
                lambda2 = 1.0 / deg as f64;
                lambda2 *= (theknots[(nbknots - 1) as usize] - theknots[(nbknots - 2) as usize])
                    / (theknots[(nbknots - 1) as usize] - theknots[0]);
            }

            if self.my_periodic {
                let avg_v1: Vec<f64> = (1..=nb_p)
                    .map(|i| 0.5 * (v1.get(i as i32) + v2.get(i as i32)))
                    .collect();
                for (i, val) in avg_v1.iter().enumerate() {
                    v1.set(i as i32 + 1, *val);
                    v2.set(i as i32 + 1, *val);
                }
            }

            let mut sq = LeastSquare::new_bsp_no_params(
                line,
                &theknots,
                &themults,
                thefirstpt,
                thelastpt,
                cons,
                cons,
                nbpoles,
            );

            sq.perform_v1tv2t(
                &rvector_to_vecd(&the_param),
                &rvector_to_vecd(&v1),
                &rvector_to_vecd(&v2),
                lambda1,
                lambda2,
            );
            self.alldone = sq.is_done();
            self.the_multi_bsp_curve = sq.bspline_value().clone();
            let mut fv = 0.0;
            let mut e3 = 0.0;
            let mut e2 = 0.0;
            sq.error(&mut fv, &mut e3, &mut e2);
            self.currenttol3d = e3;
            self.currenttol2d = e2;
            self.tolreached = true;
        }
        let mut my_parameters = RVector::new(the_param.lower(), the_param.upper());
        for i in the_param.lower()..=the_param.upper() {
            my_parameters.set(i, the_param.get(i));
        }
        self.my_parameters = Some(my_parameters);
    }

    /// OCCT TangencyVector(Line, C, U, V) (gxx L1431-1458).
    fn tangency_vector(&self, line: &MultiLine, c: &MultiCurve, u: f64, v: &mut RVector) {
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;

        let mut j = 1;
        for i in 1..=nb_p3d {
            let (p, v1) = c.d1(i as usize, u);
            let _ = p;
            v.set(j, v1.x);
            v.set(j + 1, v1.y);
            v.set(j + 2, v1.z);
            j += 3;
        }
        let mut j = nb_p3d * 3 + 1;
        for i in (nb_p3d + 1)..=(nb_p3d + nb_p2d) {
            let (p2d, v2d) = c.d1(i as usize, u);
            v.set(j, v2d.x);
            v.set(j + 1, v2d.y);
            j += 2;
        }
    }
}

/// OCCT math_Vector -> VecD conversion for the LeastSquare calls (the rcad
/// LeastSquare takes the flat 1-based VecD).
fn rvector_to_vecd(v: &RVector) -> VecD {
    let mut r = VecD::new(v.length() as usize);
    for i in 1..=v.length() {
        r.set(i as usize, v.get(i));
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    /// OCCT anchor: the quarter unit circle sampled at 17 points (3D-only
    /// MultiLine), approximated with cutting from degree 3 to 8 — Perform
    /// reaches the tolerances, the resulting MultiBSpCurve carries the
    /// poles through the end points, and its degree stays inside the
    /// requested window (Perform/Compute, gxx L625-1085).
    #[test]
    fn bspline_compute_quarter_circle() {
        let n = 17;
        let pts: Vec<DVec3> = (0..n)
            .map(|i| {
                let t = 0.5 * std::f64::consts::FRAC_PI_2 * (i as f64) / ((n - 1) as f64);
                DVec3::new(t.cos(), t.sin(), 0.0)
            })
            .collect();
        let line = MultiLine::new_tab_p3d(&pts);

        let mut c = BSplineCompute::new(
            3,
            8,
            1.0e-6,
            1.0e-6,
            20,
            true,
            ApproxParamType::ChordLength,
            false,
        );
        c.perform(&line);

        assert!(c.is_all_approximated(), "alldone");
        assert!(c.is_tolerance_reached(), "tolreached");
        let (tol3d, tol2d) = c.error();
        assert!(tol3d <= 1.0e-6 + 1.0e-12, "tol3d={}", tol3d);
        assert!(tol2d <= 1.0e-6 + 1.0e-12, "tol2d={}", tol2d);

        let cu = c.value();
        // NbCurves() = the per-pole multipoint coordinate count (1 = 3D-only).
        assert_eq!(cu.poles[0].nb_points(), 1);
        assert!((3..=8).contains(&cu.degree), "degree={}", cu.degree);
        // The BSpline carries knots with degree+1 multiplicities at both
        // ends (the Compute mults assignment, gxx L892-893).
        assert_eq!(cu.knots.len(), cu.mults.len());
        // The end poles pass through the sampled end points (the PassPoint
        // constraint at both extremities).
        let first = cu.poles[0].point(1);
        let last = cu.poles[cu.poles.len() - 1].point(1);
        assert!(first.distance(pts[0]) < 1.0e-9, "first={:?}", first);
        assert!(
            last.distance(pts[n - 1]) < 1.0e-9,
            "last={:?}",
            last
        );
        // The stored parameters are the post-approximation values in (0,1)
        // (myParameters updated inside Compute with the save guard).
        let pars = c.parameters().expect("myParameters");
        assert!(pars.get(1) >= 0.0 && pars.get(pars.length()) <= 1.0 + 1.0e-9);
    }

    /// OCCT anchor: a two-point MultiLine takes the Interpol degenerate
    /// branch (gxx L1345-1368) — a degree-1 two-pole result whose knots are
    /// the parameter bounds.
    #[test]
    fn bspline_compute_two_points_interpol() {
        let pts = vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 2.0, 3.0)];
        let line = MultiLine::new_tab_p3d(&pts);

        let mut c = BSplineCompute::new(
            3,
            8,
            1.0e-6,
            1.0e-6,
            20,
            true,
            ApproxParamType::ChordLength,
            false,
        );
        c.perform(&line);

        assert!(c.is_all_approximated(), "alldone");
        assert!(c.is_tolerance_reached(), "tolreached");
        let cu = c.value();
        assert_eq!(cu.degree, 1, "degree={}", cu.degree);
        assert_eq!(cu.nb_poles(), 2);
        let p0 = cu.poles[0].point(1);
        let p1 = cu.poles[1].point(1);
        assert_eq!(p0, pts[0]);
        assert_eq!(p1, pts[1]);
    }
}
