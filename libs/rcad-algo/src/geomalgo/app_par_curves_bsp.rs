// OCCT AppParCurves (TKGeomBase) — the BSP variants of the gradient
// machinery, 1:1 Rust translation:
//   - AppParCurves_BSpFunction.gxx (L21-337)  -> [`BSpParFunction`]
//   - AppParCurves_BSpGradient.gxx (L22-412)  -> [`BSpGradient`]
//
// The non-BSP twins (AppParCurves_Function.gxx / AppParCurves_Gradient.gxx)
// are the landed `ParFunction` / `Gradient` of `app_par_curves.rs`; the
// AppDef instantiations AppDef_BSpParFunctionOfMyBSplGradientOfBSplineCompute
// / AppDef_MyBSplGradientOfBSplineCompute are these same gxx compiled over
// MultiLine = AppDef_MultiLine, ToolLine = AppDef_MyLineTool (the _0.cxx
// alias tables), which is the binding used here.
//
// The BSpGradient_BFGS member (AppDef_BSpGradient_BFGSOfMyBSplGradientOfBSplineCompute)
// is the landed generic [`GradientBfgs`] shell over math_BFGS —
// BSpParFunction implements [`GradientFunction`], so the same shell serves.

use glam::{DVec2, DVec3};

use rcad_kernel::math::math_matrix::{IntegerVector as IVector, Matrix, Vector as RVector};
use rcad_kernel::math::VecD;

use super::app_def::{my_line_tool, MultiLine};
use super::app_par_curves::{
    GradientBfgs, GradientFunction, LeastSquare,
};
use super::approx_int::{AppParConstraint, ConstraintCouple, MultiBSpCurve, MultiPoint};

/// OCCT AppParCurves_BSpFunction — the F = sum(||C(ui) - Ptli||) function
/// over the BSpline least square, with its gradient.
pub struct BSpParFunction {
    /// OCCT MyMultiLine.
    my_multi_line: MultiLine,
    /// OCCT MyMultiBSpCurve.
    my_multi_bsp_curve: MultiBSpCurve,
    /// OCCT myParameters.
    my_parameters: RVector,
    /// OCCT ValGrad_F.
    val_grad_f: RVector,
    /// OCCT MyF.
    my_f: Matrix,
    /// OCCT PTLX / PTLY / PTLZ (filled only when intermediate points are
    /// constrained).
    ptlx: Matrix,
    ptly: Matrix,
    ptlz: Matrix,
    /// OCCT A / DA.
    a: Matrix,
    da: Matrix,
    /// OCCT MyLeastSquare.
    my_least_square: LeastSquare,
    first_p: i32,
    last_p: i32,
    nb_p: i32,
    a_deb: i32,
    a_fin: i32,
    nbpoles: i32,
    nbcu: i32,
    contraintes: bool,
    /// OCCT myConstraints.
    my_constraints: Vec<ConstraintCouple>,
    /// OCCT tabdim (HArray1(0, NbCu-1)).
    tabdim: Vec<i32>,
    /// OCCT FVal / ERR3d / ERR2d / Done.
    f_val: f64,
    err3d: f64,
    err2d: f64,
    done: bool,
    /// OCCT mylambda1 / mylambda2.
    mylambda1: f64,
    mylambda2: f64,
}

impl BSpParFunction {
    /// OCCT AppParCurves_BSpFunction(SSP, FirstPoint, LastPoint,
    /// TheConstraints, Parameters, Knots, Mults, NbPol) (gxx ctor L21-118).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &VecD,
        knots: &[f64],
        mults: &[i32],
        nb_pol: i32,
    ) -> Self {
        let nb3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb2d = my_line_tool::nb_p2d(ssp) as i32;
        let nb_pts = nb3d + nb2d;
        let mut f = BSpParFunction {
            my_multi_line: ssp.clone(),
            my_multi_bsp_curve: MultiBSpCurve::new_nbpol(nb_pol as usize),
            my_parameters: RVector::new(1, parameters.len() as i32),
            val_grad_f: RVector::new(first_point, last_point),
            my_f: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            ptlx: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            ptly: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            ptlz: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            a: Matrix::new(first_point, last_point, 1, nb_pol),
            da: Matrix::new(first_point, last_point, 1, nb_pol),
            my_least_square: LeastSquare::new_bsp_no_params(
                ssp,
                knots,
                mults,
                first_point,
                last_point,
                Self::first_constraint(the_constraints, first_point),
                Self::last_constraint(the_constraints, last_point),
                nb_pol,
            ),
            first_p: 0,
            last_p: 0,
            nb_p: 0,
            a_deb: 0,
            a_fin: 0,
            nbpoles: 0,
            nbcu: 0,
            contraintes: false,
            my_constraints: the_constraints.to_vec(),
            tabdim: Vec::new(),
            f_val: 0.0,
            err3d: 0.0,
            err2d: 0.0,
            done: false,
            mylambda1: 0.0,
            mylambda2: 0.0,
        };
        for i in 1..=(parameters.len() as i32) {
            let v = parameters.get(i as usize);
            f.my_parameters.set(i, v);
        }
        // MyMultiBSpCurve.SetKnots(Knots); SetMultiplicities(Mults).
        f.my_multi_bsp_curve.set_knots(knots);
        f.my_multi_bsp_curve.set_multiplicities_i32(mults);
        f.first_p = first_point;
        f.last_p = last_point;
        f.nb_p = f.last_p - f.first_p + 1;
        f.a_deb = f.first_p;
        f.a_fin = f.last_p;
        f.nbpoles = nb_pol;
        f.contraintes = false;
        for couple in the_constraints {
            let cons = couple.constraint;
            let myindex = couple.index;
            if myindex == f.first_p {
                if cons as i32 >= 1 {
                    f.a_deb += 1;
                }
            } else if myindex == f.last_p {
                if cons as i32 >= 1 {
                    f.a_fin -= 1;
                }
            } else if cons as i32 >= 1 {
                f.contraintes = true;
            }
        }
        let mut mynb3d = nb3d;
        let mut mynb2d = nb2d;
        if nb3d == 0 {
            mynb3d = 1;
        }
        if nb2d == 0 {
            mynb2d = 1;
        }
        f.nbcu = nb3d + nb2d;
        f.tabdim = vec![0i32; f.nbcu as usize];
        if f.contraintes {
            for i in 1..=f.nbcu {
                if i <= nb3d {
                    f.tabdim[(i - 1) as usize] = 3;
                } else {
                    f.tabdim[(i - 1) as usize] = 2;
                }
            }
            let mut tab_p = vec![DVec3::ZERO; mynb3d as usize];
            let mut tab_p2d = vec![DVec2::ZERO; mynb2d as usize];
            for i in f.first_p..=f.last_p {
                if nb3d != 0 && nb2d != 0 {
                    my_line_tool::value_3d_2d(ssp, i as usize, &mut tab_p, &mut tab_p2d);
                } else if nb3d != 0 {
                    my_line_tool::value_3d(ssp, i as usize, &mut tab_p);
                } else {
                    my_line_tool::value_2d(ssp, i as usize, &mut tab_p2d);
                }
                for j in 1..=f.nbcu {
                    if f.tabdim[(j - 1) as usize] == 3 {
                        let p = tab_p[(j - 1) as usize];
                        f.ptlx.set(i, j, p.x);
                        f.ptly.set(i, j, p.y);
                        f.ptlz.set(i, j, p.z);
                    } else {
                        let p = tab_p2d[(j - 1) as usize];
                        f.ptlx.set(i, j, p.x);
                        f.ptly.set(i, j, p.y);
                    }
                }
            }
        }
        f
    }

    /// OCCT FirstConstraint(TheConstraints, FirstPoint) (gxx L120-142).
    pub fn first_constraint(
        the_constraints: &[ConstraintCouple],
        first_point: i32,
    ) -> AppParConstraint {
        let mut cons = AppParConstraint::NoConstraint;
        for couple in the_constraints {
            cons = couple.constraint;
            let myindex = couple.index;
            if myindex == first_point {
                break;
            }
        }
        cons
    }

    /// OCCT LastConstraint(TheConstraints, LastPoint) (gxx L144-166).
    pub fn last_constraint(
        the_constraints: &[ConstraintCouple],
        last_point: i32,
    ) -> AppParConstraint {
        let mut cons = AppParConstraint::NoConstraint;
        for couple in the_constraints {
            cons = couple.constraint;
            let myindex = couple.index;
            if myindex == last_point {
                break;
            }
        }
        cons
    }

    /// OCCT Value(X, F) (gxx L168-198).
    pub fn value(&mut self, x: &VecD, f: &mut f64) -> bool {
        // myParameters = X.
        for i in 1..=self.my_parameters.length() {
            let v = x.get(i as usize);
            self.my_parameters.set(i, v);
        }
        // Resolution moindres carres:
        // ===========================
        self.my_least_square.perform_l1l2(
            &{
                let mut v = VecD::new(self.my_parameters.length() as usize);
                for i in 1..=self.my_parameters.length() {
                    v.set(i as usize, self.my_parameters.get(i));
                }
                v
            },
            self.mylambda1,
            self.mylambda2,
        );
        if !self.my_least_square.is_done() {
            self.done = false;
            return false;
        }
        if !self.contraintes {
            let mut fval = 0.0;
            let mut e3 = 0.0;
            let mut e2 = 0.0;
            self.my_least_square.error(&mut fval, &mut e3, &mut e2);
            self.f_val = fval;
            self.err3d = e3;
            self.err2d = e2;
            *f = self.f_val;
        }
        // Resolution avec contraintes:
        // ============================
        // (the OCCT constrained branch body is empty for the BSP function)
        true
    }

    /// OCCT Perform(X) (gxx L200-228).
    pub fn perform(&mut self, x: &VecD) {
        // myParameters = X.
        for i in 1..=self.my_parameters.length() {
            let v = x.get(i as usize);
            self.my_parameters.set(i, v);
        }
        // Resolution moindres carres:
        // ===========================
        self.my_least_square.perform_l1l2(
            &{
                let mut v = VecD::new(self.my_parameters.length() as usize);
                for i in 1..=self.my_parameters.length() {
                    v.set(i as usize, self.my_parameters.get(i));
                }
                v
            },
            self.mylambda1,
            self.mylambda2,
        );
        if !self.my_least_square.is_done() {
            self.done = false;
            return;
        }
        for j in 1..=self.val_grad_f.length() {
            self.val_grad_f.set(j, 0.0);
        }
        if !self.contraintes {
            let mut fval = 0.0;
            let mut e3 = 0.0;
            let mut e2 = 0.0;
            let mut grad = VecD::new(self.val_grad_f.length() as usize);
            self.my_least_square.error_gradient(&mut grad, &mut fval, &mut e3, &mut e2);
            self.f_val = fval;
            self.err3d = e3;
            self.err2d = e2;
            for i in 1..=self.val_grad_f.length() {
                self.val_grad_f.set(i, grad.get(i as usize));
            }
        }
        // Resolution avec contraintes:
        // ============================
        // (the OCCT constrained branch body is empty for the BSP function)
    }

    /// OCCT SetFirstLambda(l1) (gxx L230-233).
    pub fn set_first_lambda(&mut self, l1: f64) {
        self.mylambda1 = l1;
    }

    /// OCCT SetLastLambda(l2) (gxx L235-238).
    pub fn set_last_lambda(&mut self, l2: f64) {
        self.mylambda2 = l2;
    }

    /// OCCT NbVariables() (gxx L240-243).
    pub fn nb_variables(&self) -> i32 {
        self.nb_p
    }

    /// OCCT Gradient(X, G) (gxx L245-250).
    pub fn gradient(&mut self, x: &VecD, g: &mut RVector) -> bool {
        self.perform(x);
        for i in 1..=self.val_grad_f.length() {
            g.set(i, self.val_grad_f.get(i));
        }
        true
    }

    /// OCCT Values(X, F, G) (gxx L252-258).
    pub fn values(&mut self, x: &VecD, f: &mut f64, g: &mut RVector) -> bool {
        self.perform(x);
        *f = self.f_val;
        for i in 1..=self.val_grad_f.length() {
            g.set(i, self.val_grad_f.get(i));
        }
        true
    }

    /// OCCT CurveValue() (gxx L278-283).
    pub fn curve_value(&mut self) -> MultiBSpCurve {
        if !self.contraintes {
            self.my_multi_bsp_curve = self.my_least_square.bspline_value().clone();
        }
        self.my_multi_bsp_curve.clone()
    }

    /// OCCT Error(IPoint, CurveIndex) (gxx L285-293).
    pub fn error(&mut self, ipoint: i32, curve_index: i32) -> f64 {
        let d = self.my_least_square.distance().get(ipoint as usize, curve_index as usize);
        if !self.contraintes {
            d
        } else {
            self.my_f.get(ipoint, curve_index).sqrt()
        }
    }

    /// OCCT MaxError3d() (gxx L295-298).
    pub fn max_error_3d(&self) -> f64 {
        self.err3d
    }

    /// OCCT MaxError2d() (gxx L300-303).
    pub fn max_error_2d(&self) -> f64 {
        self.err2d
    }

    /// OCCT NewParameters() (gxx L305-308).
    pub fn new_parameters(&self) -> &RVector {
        &self.my_parameters
    }

    /// OCCT FunctionMatrix() (gxx L310-313).
    pub fn function_matrix(&self) -> &Matrix {
        self.my_least_square.function_matrix()
    }

    /// OCCT DerivativeFunctionMatrix() (gxx L315-318).
    pub fn derivative_function_matrix(&self) -> &Matrix {
        self.my_least_square.derivative_function_matrix()
    }

    /// OCCT Index() (gxx L320-323).
    pub fn index(&self) -> &IVector {
        self.my_least_square.k_index()
    }

    /// OCCT MyConstraints accessor for the ResolConstraint consumers (the
    /// myConstraints member).
    pub(crate) fn constraints(&self) -> &[ConstraintCouple] {
        &self.my_constraints
    }
}

/// OCCT AppParCurves_BSpGradient — the BSpline parameter-minimization
/// engine (the MyBSplGradient template body).
pub struct BSpGradient {
    /// OCCT SCU.
    scu: MultiBSpCurve,
    /// OCCT ParError.
    par_error: RVector,
    /// OCCT AvError / MError3d / MError2d / Done.
    av_error: f64,
    m_error3d: f64,
    m_error2d: f64,
    done: bool,
}

impl BSpGradient {
    /// OCCT AppParCurves_BSpGradient(SSP, FirstPoint, LastPoint,
    /// TheConstraints, Parameters, Knots, Mults, Deg, Tol3d, Tol2d,
    /// NbIterations) (gxx ctor L62-82) — the lambda-undefined form.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut VecD,
        knots: &[f64],
        mults: &[i32],
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut g = BSpGradient {
            scu: MultiBSpCurve::new_nbpol(1),
            par_error: RVector::new_init(first_point, last_point, 0.0),
            av_error: 0.0,
            m_error3d: 0.0,
            m_error2d: 0.0,
            done: false,
        };
        g.perform(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            knots,
            mults,
            deg,
            tol3d,
            tol2d,
            nb_iterations,
            None,
        );
        g
    }

    /// OCCT AppParCurves_BSpGradient(..., lambda1, lambda2) (gxx ctor
    /// L84-104) — the lambda-defined form.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_lambdas(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut VecD,
        knots: &[f64],
        mults: &[i32],
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
        lambda1: f64,
        lambda2: f64,
    ) -> Self {
        let mut g = BSpGradient {
            scu: MultiBSpCurve::new_nbpol(1),
            par_error: RVector::new_init(first_point, last_point, 0.0),
            av_error: 0.0,
            m_error3d: 0.0,
            m_error2d: 0.0,
            done: false,
        };
        g.perform(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            knots,
            mults,
            deg,
            tol3d,
            tol2d,
            nb_iterations,
            Some((lambda1, lambda2)),
        );
        g
    }

    /// OCCT Perform(...) (gxx L106-380) — the lambdas parameter carries the
    /// two ctor variants (None = myIsLambdaDefined false).
    #[allow(clippy::too_many_arguments)]
    fn perform(
        &mut self,
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut VecD,
        knots: &[f64],
        mults: &[i32],
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
        lambdas: Option<(f64, f64)>,
    ) {
        let (mut mylambda1, mut mylambda2, my_is_lambda_defined) = match lambdas {
            Some((l1, l2)) => (l1, l2, true),
            None => (0.0, 0.0, false),
        };
        let nb_p3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb_p2d = my_line_tool::nb_p2d(ssp) as i32;
        let nb_p = nb_p3d + nb_p2d;
        self.done = false;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_p = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p2d = vec![DVec2::ZERO; mynb_p2d as usize];
        let _ = (&mut tab_p, &mut tab_p2d);

        // Calculation of the function F= sum(||C(ui)-Ptli||2):
        // ================================================================
        let mut nbpoles = -deg - 1;
        for m in mults {
            nbpoles += *m;
        }
        let mut the_poles = vec![DVec3::ZERO; (nbpoles as usize) * (mynb_p3d as usize)];
        let mut the_poles2d = vec![DVec2::ZERO; (nbpoles as usize) * (mynb_p2d as usize)];
        let first_cons = BSpParFunction::first_constraint(the_constraints, first_point);
        let last_cons = BSpParFunction::last_constraint(the_constraints, last_point);
        let mut my_f = BSpParFunction::new(
            ssp,
            first_point,
            last_point,
            the_constraints,
            parameters,
            knots,
            mults,
            nbpoles,
        );
        if first_cons >= AppParConstraint::TangencyPoint
            || last_cons >= AppParConstraint::TangencyPoint
        {
            if !my_is_lambda_defined {
                let thefitt = LeastSquare::new_bsp(
                    ssp,
                    knots,
                    mults,
                    first_point,
                    last_point,
                    first_cons,
                    last_cons,
                    parameters,
                    nbpoles,
                );
                if first_cons >= AppParConstraint::TangencyPoint {
                    mylambda1 = thefitt.first_lambda();
                    my_f.set_first_lambda(mylambda1);
                }
                if last_cons >= AppParConstraint::TangencyPoint {
                    mylambda2 = thefitt.last_lambda();
                    my_f.set_last_lambda(mylambda2);
                }
            } else {
                my_f.set_first_lambda(mylambda1);
                my_f.set_last_lambda(mylambda2);
            }
        }
        let mut fval = 0.0;
        let _ = my_f.value(parameters, &mut fval);
        self.m_error3d = my_f.max_error_3d();
        self.m_error2d = my_f.max_error_2d();
        self.scu = my_f.curve_value();
        if self.m_error3d > tol3d || self.m_error2d > tol2d {
            // Storage of curve poles for projection:
            // ============================================
            let mut i2 = 0i32;
            for k in 1..=nb_p3d {
                let mut tab_pole: Vec<DVec3> = Vec::new();
                self.scu.curve(k as usize, &mut tab_pole);
                for j in 1..=nbpoles {
                    the_poles[(j - 1 + i2) as usize] = tab_pole[(j - 1) as usize];
                }
                i2 += nbpoles;
            }
            let mut i2 = 0i32;
            for k in 1..=nb_p2d {
                let mut tab_pole2d: Vec<DVec2> = Vec::new();
                self.scu.curve2d(k as usize, &mut tab_pole2d);
                for j in 1..=nbpoles {
                    the_poles2d[(j - 1 + i2) as usize] = tab_pole2d[(j - 1) as usize];
                }
                i2 += nbpoles;
            }
            //  Une iteration rapide de projection est faite par la methode de
            //  Rogers & Fog 89, methode equivalente a Hoschek 88 qui ne
            //  necessite pas le calcul de D2.
            let a = my_f.function_matrix().clone();
            let da = my_f.derivative_function_matrix().clone();
            let index = my_f.index().clone();
            for j in (first_point + 1)..=(last_point - 1) {
                let uf = parameters.get(j as usize);
                if nb_p3d != 0 && nb_p2d != 0 {
                    my_line_tool::value_3d_2d(ssp, j as usize, &mut tab_p, &mut tab_p2d);
                } else if nb_p2d != 0 {
                    my_line_tool::value_2d(ssp, j as usize, &mut tab_p2d);
                } else {
                    my_line_tool::value_3d(ssp, j as usize, &mut tab_p);
                }
                let mut fu = 0.0;
                let mut dfu = 0.0;
                let mut i2 = 0i32;
                let indexdeb = index.get(j) + 1;
                let indexfin = indexdeb + deg;
                for k in 1..=nb_p3d {
                    let mut a_acc = 0.0;
                    let mut b_acc = 0.0;
                    let mut c_acc = 0.0;
                    let mut d_acc = 0.0;
                    let mut e_acc = 0.0;
                    let mut f_acc = 0.0;
                    for l in indexdeb..=indexfin {
                        let pt = the_poles[(l - 1 + i2) as usize];
                        let px = pt.x;
                        let py = pt.y;
                        let pz = pt.z;
                        let aa = a.get(j, l);
                        let daa = da.get(j, l);
                        a_acc += aa * px;
                        d_acc += daa * px;
                        b_acc += aa * py;
                        e_acc += daa * py;
                        c_acc += aa * pz;
                        f_acc += daa * pz;
                    }
                    let pt = DVec3::new(a_acc, b_acc, c_acc);
                    let v1 = DVec3::new(d_acc, e_acc, f_acc);
                    i2 += nbpoles;
                    let my_v = tab_p[(k - 1) as usize] - pt;
                    fu += my_v.dot(v1);
                    dfu += v1.length_squared();
                }
                let mut i2 = 0i32;
                for k in 1..=nb_p2d {
                    let mut a_acc = 0.0;
                    let mut b_acc = 0.0;
                    let mut d_acc = 0.0;
                    let mut e_acc = 0.0;
                    for l in indexdeb..=indexfin {
                        let pt2d = the_poles2d[(l - 1 + i2) as usize];
                        let px = pt2d.x;
                        let py = pt2d.y;
                        let aa = a.get(j, l);
                        let daa = da.get(j, l);
                        a_acc += aa * px;
                        d_acc += daa * px;
                        b_acc += aa * py;
                        e_acc += daa * py;
                    }
                    let pt2d = DVec2::new(a_acc, b_acc);
                    let v12d = DVec2::new(d_acc, e_acc);
                    i2 += nbpoles;
                    let my_v2d = tab_p2d[(k - 1) as usize] - pt2d;
                    fu += my_v2d.dot(v12d);
                    dfu += v12d.length_squared();
                }
                // OCCT RealEpsilon().
                if dfu >= f64::EPSILON {
                    let mut du = fu / dfu;
                    du = du.abs().min(5.0e-02) * du.signum(); // copysign(min(5e-2, |DU|), DU)
                    let uf = uf + du;
                    parameters.set(j as usize, uf);
                }
            }
            let mut fval = 0.0;
            let _ = my_f.value(parameters, &mut fval);
            self.m_error3d = my_f.max_error_3d();
            self.m_error2d = my_f.max_error_2d();
        }
        if self.m_error3d <= tol3d && self.m_error2d <= tol2d {
            self.done = true;
        } else if nb_iterations != 0 {
            // NbIterations de gradient conjugue:
            // =================================
            let eps = 1.0e-07;
            // AppParCurves_BSpGradient_BFGS FResol(MyF, Parameters, Tol3d,
            // Tol2d, Eps, NbIterations).
            let mut parameters_copy = parameters.clone();
            let _f_resol = GradientBfgs::new(
                &mut my_f,
                &parameters_copy,
                tol3d,
                tol2d,
                eps,
                nb_iterations,
            );
            let _ = &mut parameters_copy;
        }
        self.scu = my_f.curve_value();
        self.av_error = 0.0;
        for j in first_point..=last_point {
            // Recherche des erreurs maxi et moyenne a un index donne:
            for k in 1..=nb_p {
                let e = my_f.error(j, k);
                let v = self.par_error.get(j).max(e);
                self.par_error.set(j, v);
            }
            let v = self.av_error + self.par_error.get(j);
            self.av_error = v;
        }
        self.av_error = self.av_error / (last_point - first_point + 1) as f64;
        self.m_error3d = my_f.max_error_3d();
        self.m_error2d = my_f.max_error_2d();
        if self.m_error3d <= tol3d && self.m_error2d <= tol2d {
            self.done = true;
        }
    }

    /// OCCT Value() (gxx L382-385).
    pub fn value(&self) -> MultiBSpCurve {
        self.scu.clone()
    }

    /// OCCT IsDone() (gxx L387-390).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Error(Index) (gxx L392-395).
    pub fn error(&self, index: i32) -> f64 {
        self.par_error.get(index)
    }

    /// OCCT AverageError() (gxx L397-400).
    pub fn average_error(&self) -> f64 {
        self.av_error
    }

    /// OCCT MaxError3d() (gxx L402-405).
    pub fn max_error_3d(&self) -> f64 {
        self.m_error3d
    }

    /// OCCT MaxError2d() (gxx L407-410).
    pub fn max_error_2d(&self) -> f64 {
        self.m_error2d
    }
}

// The MultiPoint import keeps the parity with the MultiBSpCurve pole
// storage (AppParCurves_MultiPoint).
#[allow(unused)]
fn _multipoint_parity(_mp: &MultiPoint) {}

// ---------------------------------------------------------------------------
// The math_MultipleVarFunctionWithGradient interface of BSpParFunction (the
// math_BFGS consumption) — the same delegation shape as the landed
// ParFunction.
// ---------------------------------------------------------------------------

impl rcad_kernel::math::math_bfgs::MultipleVarFunction for BSpParFunction {
    fn nb_variables(&self) -> i32 {
        BSpParFunction::nb_variables(self)
    }
    fn value(&mut self, x: &rcad_kernel::math::math_matrix::Vector, f: &mut f64) -> bool {
        let mut xv = VecD::new(x.length() as usize);
        for i in 1..=x.length() {
            xv.set(i as usize, x.get(i));
        }
        BSpParFunction::value(self, &xv, f)
    }
}
impl rcad_kernel::math::math_bfgs::MultipleVarFunctionWithGradient for BSpParFunction {
    fn gradient(&mut self, x: &rcad_kernel::math::math_matrix::Vector, g: &mut rcad_kernel::math::math_matrix::Vector) -> bool {
        let mut xv = VecD::new(x.length() as usize);
        for i in 1..=x.length() {
            xv.set(i as usize, x.get(i));
        }
        let mut gv = RVector::new(g.lower, g.upper());
        let ok = BSpParFunction::gradient(self, &xv, &mut gv);
        for i in g.lower..=g.upper() {
            let v = gv.get(i);
            g.set(i, v);
        }
        ok
    }
    fn values(&mut self, x: &rcad_kernel::math::math_matrix::Vector, f: &mut f64, g: &mut rcad_kernel::math::math_matrix::Vector) -> bool {
        let mut xv = VecD::new(x.length() as usize);
        for i in 1..=x.length() {
            xv.set(i as usize, x.get(i));
        }
        let mut gv = RVector::new(g.lower, g.upper());
        let ok = BSpParFunction::values(self, &xv, f, &mut gv);
        for i in g.lower..=g.upper() {
            let v = gv.get(i);
            g.set(i, v);
        }
        ok
    }
}
impl GradientFunction for BSpParFunction {
    fn max_error_3d(&self) -> f64 {
        BSpParFunction::max_error_3d(self)
    }
    fn max_error_2d(&self) -> f64 {
        BSpParFunction::max_error_2d(self)
    }
}
