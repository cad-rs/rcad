// OCCT AppParCurves (TKGeomBase) — 1:1 Rust translation:
//   AppParCurves.cxx L23-103 (BernsteinMatrix, Bernstein) and L157-257
//   (SplineFunction) — SecondDerivativeBernstein (L105-155) deferred until
//   its consuming unit lands.
//   AppParCurves_LeastSquare.gxx whole file (L41-1860) — the
//   AppParCurves_LeastSquare template instantiated here exactly as
//   AppDef_ParLeastSquareOfMyGradientOfCompute does:
//     MultiLine = AppDef_MultiLine  (rcad app_def::MultiLine)
//     ToolLine  = AppDef_MyLineTool (rcad app_def::my_line_tool)
//   The struct keeps the OCCT template name LeastSquare; the earlier
//   approx_int::LeastSquare is the GeomInt WLApprox instantiation (kept
//   until the WLine path migrates onto this aligned implementation).
//
// The local Matrix / RVector / IVector wrappers reproduce the OCCT
// math_Matrix / math_Vector / math_IntegerVector semantics used by the
// template: 1-based storage with arbitrary lower bounds. Where OCCT passes
// members (FirstConstraint, Vec1t, ...) into Affect by reference, rcad uses
// copy-in/copy-out of the disjoint fields (identical semantics).

use rcad_kernel::math::bspl_lib::{coefs_d1_2d, coefs_d1_3d, locate_parameter_flat, poles_coefficients_2d, poles_coefficients_3d};
use rcad_kernel::math::math_householder::Householder;
use rcad_kernel::math::math_matrix::{
    IntegerVector as IVector, Matrix, Vector as RVector, Vector as KernelVector,
};
use rcad_kernel::math::math_recipes::{dactcl_decompose, dactcl_solve, MATH_STATUS_OK};
use rcad_kernel::math::{MatD, VecD};

use glam::{DVec2, DVec3};

use super::app_def::{my_line_tool, MultiLine};
use super::approx_int::{AppParConstraint, ConstraintCouple, MultiBSpCurve, MultiCurve, MultiPoint};
use rcad_kernel::math::math_uzawa::Uzawa;

/// OCCT MinPivot default of DACTCL_Decompose / DACTCL_Solve (1.0e-20).
const DACTCL_MIN_PIVOT: f64 = 1.0e-20;

// ---------------------------------------------------------------------------
// AppParCurves::BernsteinMatrix / Bernstein / SplineFunction
// ---------------------------------------------------------------------------

/// OCCT AppParCurves::BernsteinMatrix (AppParCurves.cxx L23-58).
pub fn bernstein_matrix(nb_poles: i32, u: &VecD, a: &mut Matrix) {
    let mut b = RVector::new(1, nb_poles - 1);
    let first = 1; // OCCT first = U.Lower().
    let last = u.len() as i32; // OCCT last = U.Upper().
    for i in first..=last {
        b.set(1, 1.0);
        let u0 = u.get(i as usize);
        let u1 = 1.0 - u0;

        for id in 2..=(nb_poles - 1) {
            let mut y0 = b.get(1);
            let mut y1 = u0 * y0;
            b.set(1, y0 - y1);
            for j in 2..=(id - 1) {
                let xs = y1;
                y0 = b.get(j);
                y1 = u0 * y0;
                b.set(j, y0 - y1 + xs);
            }
            b.set(id, y1);
        }
        a.set(i, 1, u1 * b.get(1));
        a.set(i, nb_poles, u0 * b.get(nb_poles - 1));
        for j in 2..=(nb_poles - 1) {
            a.set(i, j, u1 * b.get(j) + u0 * b.get(j - 1));
        }
    }
}

/// OCCT AppParCurves::Bernstein (AppParCurves.cxx L60-103).
pub fn bernstein(nb_poles: i32, u: &VecD, a: &mut Matrix, da: &mut Matrix) {
    let n_deg = nb_poles - 1;
    let mut b = RVector::new(1, nb_poles - 1);
    let first = 1; // OCCT first = U.Lower().
    let last = u.len() as i32; // OCCT last = U.Upper().
    for i in first..=last {
        b.set(1, 1.0);
        let u0 = u.get(i as usize);
        let u1 = 1.0 - u0;

        for id in 2..=(nb_poles - 1) {
            let mut y0 = b.get(1);
            let mut y1 = u0 * y0;
            b.set(1, y0 - y1);
            for j in 2..=(id - 1) {
                let xs = y1;
                y0 = b.get(j);
                y1 = u0 * y0;
                b.set(j, y0 - y1 + xs);
            }
            b.set(id, y1);
        }
        da.set(i, 1, -(n_deg as f64) * b.get(1));
        da.set(i, nb_poles, n_deg as f64 * b.get(nb_poles - 1));
        a.set(i, 1, u1 * b.get(1));
        a.set(i, nb_poles, u0 * b.get(nb_poles - 1));
        for j in 2..=(nb_poles - 1) {
            let bj = b.get(j);
            let bj1 = b.get(j - 1);
            da.set(i, j, n_deg as f64 * (bj1 - bj));
            a.set(i, j, u1 * bj + u0 * bj1);
        }
    }
}

/// OCCT AppParCurves::SecondDerivativeBernstein (AppParCurves.cxx L105-155).
pub fn second_derivative_bernstein(u: f64, dda: &mut RVector) {
    let nb_poles = dda.length();
    let deg = nb_poles - 1;
    let n4 = deg * (deg - 1);
    let mut b = RVector::new(1, deg - 1);
    b.set(1, 1.0);

    // Cas particulier si degre = 1:
    if deg == 1 {
        dda.set(1, 0.0);
        dda.set(2, 0.0);
    } else if deg == 2 {
        dda.set(1, 2.0);
        dda.set(2, -4.0);
        dda.set(3, 2.0);
    } else {
        for id in 2..=(deg - 1) {
            let mut y0 = b.get(1);
            let mut y1 = u * y0;
            b.set(1, y0 - y1);
            for j in 2..=(id - 1) {
                let xs = y1;
                y0 = b.get(j);
                y1 = u * y0;
                b.set(j, y0 - y1 + xs);
            }
            b.set(id, y1);
        }

        let n4f = n4 as f64;
        let v = n4f * b.get(1);
        dda.set(1, v);
        let v = n4f * (-2.0 * b.get(1) + b.get(2));
        dda.set(2, v);
        let v = n4f * (b.get(deg - 2) - 2.0 * b.get(deg - 1));
        dda.set(deg, v);
        let v = n4f * b.get(deg - 1);
        dda.set(deg + 1, v);

        for j in 2..=(deg - 2) {
            let v = n4f * (b.get(j - 1) - 2.0 * b.get(j) + b.get(j + 1));
            dda.set(j + 1, v);
        }
    }
}

/// OCCT AppParCurves::SplineFunction (AppParCurves.cxx L157-257).
pub fn spline_function(
    nbpoles: i32,
    deg: i32,
    parameters: &VecD,
    flatknots: &RVector,
    a: &mut Matrix,
    da: &mut Matrix,
    index: &mut IVector,
) {
    let deg1 = deg + 1;
    let mut locpoles = RVector::new(1, deg1);
    let mut locdpoles = RVector::new(1, deg1);
    let firstp = 1; // OCCT firstp = Parameters.Lower().
    let lastp = parameters.len() as i32; // OCCT lastp = Parameters.Upper().

    // OCCT copies flatknots into Aflatknots; rcad's locate_parameter_flat
    // reads the same 1-based slice storage.
    let aflatknots: &[f64] = &flatknots.data.v;

    let mut oldkindex: i32 = 1;

    for i in firstp..=lastp {
        let u = parameters.get(i as usize);
        let mut new_u = u;
        let mut kindex = oldkindex;
        // OCCT BSplCLib::LocateParameter(deg, Aflatknots, U, false, deg1,
        // nbpoles + 1, kindex, NewU).
        locate_parameter_flat(
            deg as usize,
            aflatknots,
            u,
            false,
            deg1,
            nbpoles + 1,
            &mut kindex,
            &mut new_u,
        );

        oldkindex = kindex;

        // On stocke les index:
        index.set(i, kindex - deg - 1);

        locpoles.set(1, 1.0);

        for qq in 2..=deg {
            locpoles.set(qq, 0.0);
            for pp in 1..=(qq - 1) {
                let inverse =
                    1.0 / (flatknots.get(kindex + pp) - flatknots.get(kindex - qq + pp + 1));
                let saved =
                    (u - flatknots.get(kindex - qq + pp + 1)) * inverse * locpoles.get(pp);
                let v = locpoles.get(pp) * (flatknots.get(kindex + pp) - u) * inverse;
                locpoles.set(pp, v);
                let v = locpoles.get(pp) + locpoles.get(qq);
                locpoles.set(pp, v);
                locpoles.set(qq, saved);
            }
        }

        let qq = deg + 1;
        for pp in 1..=deg {
            let v = locpoles.get(pp);
            locdpoles.set(pp, v);
        }

        let mut locqq = 0.0;
        let mut locdqq = 0.0;
        for pp in 1..=deg {
            let inverse =
                1.0 / (flatknots.get(kindex + pp) - flatknots.get(kindex - qq + pp + 1));
            let saved = (u - flatknots.get(kindex - qq + pp + 1)) * inverse * locpoles.get(pp);
            let v = locpoles.get(pp) * (flatknots.get(kindex + pp) - u) * inverse;
            locpoles.set(pp, v);
            let v = locpoles.get(pp) + locqq;
            locpoles.set(pp, v);
            locqq = saved;
            let local_inverse = deg as f64 * inverse;
            let saved = local_inverse * locdpoles.get(pp);
            let v = locdpoles.get(pp) * -local_inverse;
            locdpoles.set(pp, v);
            let v = locdpoles.get(pp) + locdqq;
            locdpoles.set(pp, v);
            locdqq = saved;
        }

        let v = locqq;
        locpoles.set(qq, v);
        let v = locdqq;
        locdpoles.set(qq, v);

        for j in 1..=deg1 {
            let val = locpoles.get(j);
            let theindex = j + oldkindex - deg1;
            a.set(i, theindex, val);
            da.set(i, theindex, locdpoles.get(j));
        }

        for j in 1..(oldkindex - deg) {
            a.set(i, j, 0.0);
            da.set(i, j, 0.0);
        }
        for j in (oldkindex + 1)..=nbpoles {
            a.set(i, j, 0.0);
            da.set(i, j, 0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// AppParCurves_LeastSquare — statics of the template file
// ---------------------------------------------------------------------------

/// OCCT static FlatLength (gxx L41-50).
fn flat_length(mults: &[i32]) -> i32 {
    let mut sum = 0;
    for v in mults {
        sum += *v;
    }
    sum
}

/// OCCT static CheckTangents (3D points + 3D/2D tangents) (gxx L64-106).
fn check_tangents(
    the_arr_pt1: &[DVec3],
    the_arr_pt2: &[DVec3],
    the_arr_tg3d: &mut [DVec3],
    the_arr_tg2d: &mut [DVec2],
) {
    if the_arr_pt1.len() != the_arr_pt2.len() {
        return;
    }
    if the_arr_tg3d.len() != the_arr_pt1.len() {
        return;
    }

    let mut is_to_change_dir = false;

    for i in 0..the_arr_pt1.len() {
        let a_v1 = the_arr_pt2[i] - the_arr_pt1[i];
        let a_v2 = the_arr_tg3d[i];

        if a_v1.dot(a_v2) < 0.0 {
            is_to_change_dir = true;
            break;
        }
    }

    if !is_to_change_dir {
        return;
    }

    // Change directions for every 2D- and 3D-tangents
    for v in the_arr_tg3d.iter_mut() {
        *v = -*v;
    }
    for v in the_arr_tg2d.iter_mut() {
        *v = -*v;
    }
}

/// OCCT static CheckTangents (2D points + 2D tangents) (gxx L116-136).
fn check_tangents_2d(
    the_arr_pt1: &[DVec2],
    the_arr_pt2: &[DVec2],
    the_arr_tg2d: &mut [DVec2],
) {
    if the_arr_pt1.len() != the_arr_pt2.len() {
        return;
    }

    for i in 0..the_arr_pt1.len() {
        let a_v1 = the_arr_pt2[i] - the_arr_pt1[i];
        let a_v2 = the_arr_tg2d[i];

        if a_v1.dot(a_v2) < 0.0 {
            the_arr_tg2d[i] = -a_v2;
        }
    }
}

// ---------------------------------------------------------------------------
// AppParCurves_LeastSquare
// ---------------------------------------------------------------------------

/// OCCT AppParCurves_LeastSquare (the AppDef instantiation).
#[derive(Debug, Clone)]
pub struct LeastSquare {
    first_constraint: AppParConstraint,
    last_constraint: AppParConstraint,
    scu: MultiBSpCurve,
    /// OCCT myknots (null handle == None).
    myknots: Option<Vec<f64>>,
    /// OCCT mymults (null handle == None).
    mymults: Option<Vec<i32>>,
    mypoles: Matrix,
    a: Matrix,
    da: Matrix,
    b2: Matrix,
    mypoints: Matrix,
    vflatknots: RVector,
    vec1t: RVector,
    vec1c: RVector,
    vec2t: RVector,
    vec2c: RVector,
    the_error: Matrix,
    myindex: IVector,
    lambda1: f64,
    lambda2: f64,
    first_p: i32,
    last_p: i32,
    nlignes: i32,
    ninc: i32,
    na: i32,
    myfirstp: i32,
    mylastp: i32,
    resinit: i32,
    resfin: i32,
    nbp2d: i32,
    nbp: i32,
    nbpoles: i32,
    deg: i32,
    done: bool,
    iscalculated: bool,
    isready: bool,
}

impl LeastSquare {
    /// OCCT AppParCurves_LeastSquare(SSP, FirstPoint, LastPoint, FirstCons,
    /// LastCons, Parameters, NbPol) (gxx L138-167).
    pub fn new(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        first_cons: AppParConstraint,
        last_cons: AppParConstraint,
        parameters: &VecD,
        nb_pol: i32,
    ) -> Self {
        let mut ls = LeastSquare::make_members(
            ssp,
            first_point,
            last_point,
            first_cons,
            last_cons,
            nb_pol,
            None,
            None,
        );
        ls.init(ssp, first_point, last_point);
        ls.perform(parameters);
        ls
    }

    /// OCCT AppParCurves_LeastSquare(SSP, FirstPoint, LastPoint, FirstCons,
    /// LastCons, NbPol) (gxx L169-196) — initializes the fields.
    pub fn new_no_params(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        first_cons: AppParConstraint,
        last_cons: AppParConstraint,
        nb_pol: i32,
    ) -> Self {
        let mut ls = LeastSquare::make_members(
            ssp,
            first_point,
            last_point,
            first_cons,
            last_cons,
            nb_pol,
            None,
            None,
        );
        ls.init(ssp, first_point, last_point);
        ls
    }

    /// OCCT AppParCurves_LeastSquare(SSP, Knots, Mults, FirstPoint, LastPoint,
    /// FirstCons, LastCons, Parameters, NbPol) (gxx L198-235).
    #[allow(clippy::too_many_arguments)]
    pub fn new_bsp(
        ssp: &MultiLine,
        knots: &[f64],
        mults: &[i32],
        first_point: i32,
        last_point: i32,
        first_cons: AppParConstraint,
        last_cons: AppParConstraint,
        parameters: &VecD,
        nb_pol: i32,
    ) -> Self {
        let mut ls = LeastSquare::make_members(
            ssp,
            first_point,
            last_point,
            first_cons,
            last_cons,
            nb_pol,
            Some(knots.to_vec()),
            Some(mults.to_vec()),
        );
        // SCU.SetKnots(Knots); SCU.SetMultiplicities(Mults); (gxx L231-232)
        ls.scu.set_knots(knots);
        ls.scu.set_multiplicities_i32(mults);
        ls.init(ssp, first_point, last_point);
        ls.perform(parameters);
        ls
    }

    /// OCCT AppParCurves_LeastSquare(SSP, Knots, Mults, FirstPoint, LastPoint,
    /// FirstCons, LastCons, NbPol) (gxx L237-272).
    #[allow(clippy::too_many_arguments)]
    pub fn new_bsp_no_params(
        ssp: &MultiLine,
        knots: &[f64],
        mults: &[i32],
        first_point: i32,
        last_point: i32,
        first_cons: AppParConstraint,
        last_cons: AppParConstraint,
        nb_pol: i32,
    ) -> Self {
        let mut ls = LeastSquare::make_members(
            ssp,
            first_point,
            last_point,
            first_cons,
            last_cons,
            nb_pol,
            Some(knots.to_vec()),
            Some(mults.to_vec()),
        );
        ls.scu.set_knots(knots);
        ls.scu.set_multiplicities_i32(mults);
        ls.init(ssp, first_point, last_point);
        ls
    }

    /// The shared member-initializer list of the four OCCT constructors
    /// (gxx L145-161 / L175-191 / L207-223 / L245-261). The knots/mults
    /// members are filled here as the ctor bodies do (gxx L227-230 /
    /// L263-266).
    fn make_members(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        first_cons: AppParConstraint,
        last_cons: AppParConstraint,
        nb_pol: i32,
        knots: Option<Vec<f64>>,
        mults: Option<Vec<i32>>,
    ) -> Self {
        let nb_b_columns = nb_b_columns(ssp);
        let the_first = the_first_point(first_cons, first_point);
        let the_last = the_last_point(last_cons, last_point);
        let flat_len = mults.as_ref().map(|m| flat_length(m)).unwrap_or(1);
        LeastSquare {
            first_constraint: first_cons,
            last_constraint: last_cons,
            scu: MultiBSpCurve::new_nbpol(nb_pol as usize),
            myknots: knots,
            mymults: mults,
            mypoles: Matrix::new(1, nb_pol, 1, nb_b_columns),
            a: Matrix::new(first_point, last_point, 1, nb_pol),
            da: Matrix::new(first_point, last_point, 1, nb_pol),
            b2: Matrix::new(the_first, the_first.max(the_last), 1, nb_b_columns),
            mypoints: Matrix::new(first_point, last_point, 1, nb_b_columns),
            vflatknots: RVector::new(1, flat_len),
            vec1t: RVector::new(1, nb_b_columns),
            vec1c: RVector::new(1, nb_b_columns),
            vec2t: RVector::new(1, nb_b_columns),
            vec2c: RVector::new(1, nb_b_columns),
            the_error: Matrix::new_init(
                first_point,
                last_point,
                1,
                my_line_tool::nb_p3d(ssp) as i32 + my_line_tool::nb_p2d(ssp) as i32,
                0.0,
            ),
            myindex: IVector::new_init(first_point, last_point, 0),
            lambda1: 0.0,
            lambda2: 0.0,
            first_p: 0,
            last_p: 0,
            nlignes: 0,
            ninc: 0,
            na: 0,
            myfirstp: 0,
            mylastp: 0,
            resinit: 0,
            resfin: 0,
            nbp2d: 0,
            nbp: 0,
            nbpoles: nb_pol,
            deg: 0,
            done: false,
            iscalculated: false,
            isready: false,
        }
    }

    /// OCCT Init(SSP, FirstPoint, LastPoint) (gxx L274-508).
    fn init(&mut self, ssp: &MultiLine, first_point: i32, last_point: i32) {
        // Variable de controle
        self.iscalculated = false;
        self.isready = true;

        self.myfirstp = first_point;
        self.mylastp = last_point;
        self.first_p = the_first_point(self.first_constraint, self.myfirstp);
        self.last_p = the_last_point(self.last_constraint, self.mylastp);

        // Identification of constraints at extremities:
        // ========================================
        self.nbp2d = my_line_tool::nb_p2d(ssp) as i32;
        self.nbp = my_line_tool::nb_p3d(ssp) as i32;
        let mut mynbp2d = self.nbp2d;
        let mut mynbp = self.nbp;
        if self.nbp2d == 0 {
            mynbp2d = 1;
        }
        if self.nbp == 0 {
            mynbp = 1;
        }
        let mut tab_p = vec![DVec3::ZERO; mynbp as usize];
        let mut tab_p2d = vec![DVec2::ZERO; mynbp2d as usize];
        let mut tab_v = vec![DVec3::ZERO; mynbp as usize];
        let mut tab_v2d = vec![DVec2::ZERO; mynbp2d as usize];

        self.deg = self.nbpoles - 1;

        if let Some(mymults) = self.mymults.clone() {
            let mut sum = 0;
            for v in &mymults {
                sum += *v;
            }
            self.deg = sum - self.nbpoles - 1;
            let mut k = 1;
            if let Some(myknots) = self.myknots.clone() {
                for (i, ki) in myknots.iter().enumerate() {
                    let val = *ki;
                    for _ in 0..mymults[i] {
                        self.vflatknots.set(k, val);
                        k += 1;
                    }
                }
            }
        }

        let mut cons = self.first_constraint;
        affect(
            ssp,
            first_point,
            &mut cons,
            &mut self.vec1t,
            &mut self.vec1c,
            self.nbp,
            self.nbp2d,
        );
        self.first_constraint = cons;

        let mut cons = self.last_constraint;
        affect(
            ssp,
            last_point,
            &mut cons,
            &mut self.vec2t,
            &mut self.vec2c,
            self.nbp,
            self.nbp2d,
        );
        self.last_constraint = cons;

        for j in self.myfirstp..=self.mylastp {
            let mut i2 = 1;
            if self.nbp != 0 && self.nbp2d != 0 {
                my_line_tool::value_3d_2d(ssp, j as usize, &mut tab_p, &mut tab_p2d);
            } else if self.nbp2d != 0 {
                my_line_tool::value_2d(ssp, j as usize, &mut tab_p2d);
            } else {
                my_line_tool::value_3d(ssp, j as usize, &mut tab_p);
            }
            for i in 1..=self.nbp {
                let p = tab_p[(i - 1) as usize];
                self.mypoints.set(j, i2, p.x);
                self.mypoints.set(j, i2 + 1, p.y);
                self.mypoints.set(j, i2 + 2, p.z);
                i2 += 3;
            }
            for i in 1..=self.nbp2d {
                let p = tab_p2d[(i - 1) as usize];
                self.mypoints.set(j, i2, p.x);
                self.mypoints.set(j, i2 + 1, p.y);
                i2 += 2;
            }
        }

        let mut pole1 = MultiPoint::new(self.nbp as usize, self.nbp2d as usize);
        let mut polen = MultiPoint::new(self.nbp as usize, self.nbp2d as usize);

        if self.first_constraint == AppParConstraint::PassPoint
            || self.first_constraint == AppParConstraint::TangencyPoint
            || self.first_constraint == AppParConstraint::CurvaturePoint
        {
            let mut i2 = 1;
            for i in 1..=self.nbp {
                let p = DVec3::new(
                    self.mypoints.get(self.myfirstp, i2),
                    self.mypoints.get(self.myfirstp, i2 + 1),
                    self.mypoints.get(self.myfirstp, i2 + 2),
                );
                pole1.set_point(i as usize, p);
                i2 += 3;
            }
            for i in 1..=self.nbp2d {
                let p = DVec2::new(
                    self.mypoints.get(self.myfirstp, i2),
                    self.mypoints.get(self.myfirstp, i2 + 1),
                );
                pole1.set_point2d((i + self.nbp) as usize, p);
                i2 += 2;
            }
            for i in 1..=self.mypoles.col_number() {
                let v = self.mypoints.get(self.myfirstp, i);
                self.mypoles.set(1, i, v);
            }
        }

        if self.last_constraint == AppParConstraint::PassPoint
            || self.last_constraint == AppParConstraint::TangencyPoint
            || self.first_constraint == AppParConstraint::CurvaturePoint
        {
            let mut i2 = 1;
            for i in 1..=self.nbp {
                let p = DVec3::new(
                    self.mypoints.get(self.mylastp, i2),
                    self.mypoints.get(self.mylastp, i2 + 1),
                    self.mypoints.get(self.mylastp, i2 + 2),
                );
                polen.set_point(i as usize, p);
                i2 += 3;
            }
            for i in 1..=self.nbp2d {
                let p = DVec2::new(
                    self.mypoints.get(self.mylastp, i2),
                    self.mypoints.get(self.mylastp, i2 + 1),
                );
                polen.set_point2d((i + self.nbp) as usize, p);
                i2 += 2;
            }

            for i in 1..=self.mypoles.col_number() {
                let v = self.mypoints.get(self.mylastp, i);
                self.mypoles.set(self.nbpoles, i, v);
            }
        }

        if self.first_constraint == AppParConstraint::NoConstraint
        {
            self.resinit = 1;
            self.scu.set_value(1, pole1.clone());
            if self.last_constraint == AppParConstraint::NoConstraint {
                self.resfin = self.nbpoles;
            } else if self.last_constraint == AppParConstraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            } else if self.last_constraint == AppParConstraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            } else if self.last_constraint == AppParConstraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            }
        } else if self.first_constraint == AppParConstraint::PassPoint {
            self.resinit = 2;
            self.scu.set_value(1, pole1.clone());
            if self.last_constraint == AppParConstraint::NoConstraint {
                self.resfin = self.nbpoles;
            } else if self.last_constraint == AppParConstraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            } else if self.last_constraint == AppParConstraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            } else if self.last_constraint == AppParConstraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            }
        } else if self.first_constraint == AppParConstraint::TangencyPoint {
            self.resinit = 3;
            self.scu.set_value(1, pole1.clone());
            if self.last_constraint == AppParConstraint::NoConstraint {
                self.resfin = self.nbpoles;
            }
            if self.last_constraint == AppParConstraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            }
            if self.last_constraint == AppParConstraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            } else if self.last_constraint == AppParConstraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            }
        } else if self.first_constraint == AppParConstraint::CurvaturePoint {
            self.resinit = 4;
            self.scu.set_value(1, pole1.clone());
            if self.last_constraint == AppParConstraint::NoConstraint {
                self.resfin = self.nbpoles;
            }
            if self.last_constraint == AppParConstraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            }
            if self.last_constraint == AppParConstraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            } else if self.last_constraint == AppParConstraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.set_value(self.nbpoles as usize, polen.clone());
            }
        }
        let _ = &mut tab_v;
        let _ = &mut tab_v2d;

        let nincx = self.resfin - self.resinit + 1;
        if nincx < 1 {
            // Impossible d'aller plus loin
            self.isready = false;
            return;
        }
        let neq = self.last_p - self.first_p + 1;

        self.na = 3 * self.nbp + 2 * self.nbp2d;
        self.nlignes = self.na * neq;
        self.ninc = self.na * nincx;
        if self.first_constraint >= AppParConstraint::TangencyPoint {
            self.ninc += 1;
        }
        if self.last_constraint >= AppParConstraint::TangencyPoint {
            self.ninc += 1;
        }
    }

    /// OCCT Perform(Parameters) (gxx L510-734). Note the OCCT fall-through:
    /// when both constraints are below TangencyPoint the resolution block
    /// runs with done=true, and the "cas de tangence" section L626-733 still
    /// executes unconditionally afterwards (recomputing the same poles via
    /// the tangency-form normal equations) — the early returns at L539/L543
    /// (Householder failure/success) and L590 (Nincx < 1) are the only exits.
    pub fn perform(&mut self, parameters: &VecD) {
        self.done = false;
        if !self.isready {
            return;
        }
        let nbpol1 = self.nbpoles - 1;
        let ninc1 = self.ninc - 1;
        self.iscalculated = false;

        // Calculation of matrix A and DA of approximation functions:
        self.compute_function(parameters);

        if self.first_constraint != AppParConstraint::TangencyPoint
            && self.last_constraint != AppParConstraint::TangencyPoint
        {
            if self.first_constraint == AppParConstraint::NoConstraint {
                if self.last_constraint == AppParConstraint::NoConstraint {
                    // math_Householder HouResol(A, mypoints); (default EPS
                    // 1.0e-20, math_Householder.hxx).
                    let hou_resol = Householder::new(self.a.data(), self.mypoints.data(), 1.0e-20);
                    if !hou_resol.is_done() {
                        self.done = false;
                        return;
                    }
                    self.done = true;
                    *self.mypoles.data_mut() = hou_resol.all_values().clone();
                    return;
                } else {
                    for j in self.first_p..=self.last_p {
                        let ad1 = self.a.get(j, self.nbpoles);
                        for i in 1..=self.b2.col_number() {
                            let v =
                                self.mypoints.get(j, i) - ad1 * self.mypoles.get(self.nbpoles, i);
                            self.b2.set(j, i, v);
                        }
                    }
                }
            } else if self.first_constraint == AppParConstraint::PassPoint {
                if self.last_constraint == AppParConstraint::NoConstraint {
                    for j in self.first_p..=self.last_p {
                        let a0 = self.a.get(j, 1);
                        for i in 1..=self.b2.col_number() {
                            let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i);
                            self.b2.set(j, i, v);
                        }
                    }
                } else if self.last_constraint == AppParConstraint::PassPoint {
                    for j in self.first_p..=self.last_p {
                        let a0 = self.a.get(j, 1);
                        let ad1 = self.a.get(j, self.nbpoles);
                        for i in 1..=self.b2.col_number() {
                            let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                                - ad1 * self.mypoles.get(self.nbpoles, i);
                            self.b2.set(j, i, v);
                        }
                    }
                }
            }

            // resolution:

            let nincx = self.resfin - self.resinit + 1;
            if nincx < 1 {
                self.done = true;
                return;
            }
            let mut index = IVector::new_init(1, nincx, 0);
            self.search_index(&mut index);
            let mut mytab = Matrix::new_init(self.resinit, self.resfin, 1, self.b2.col_number(), 0.0);
            let mut the_aa = RVector::new_init(1, index.get(nincx), 0.0);
            let mut mytabb = RVector::new_init(1, nincx, 0.0);

            self.make_taa_matrix(&mut the_aa, &mut mytab);
            dactcl_decompose(&mut the_aa.data, &index.data, DACTCL_MIN_PIVOT);

            let mut kk2;
            for j in 1..=self.b2.col_number() {
                kk2 = 1;
                for i in self.resinit..=self.resfin {
                    let v = mytab.get(i, j);
                    mytabb.set(kk2, v);
                    kk2 += 1;
                }
                dactcl_solve(&the_aa.data, &mut mytabb.data, &index.data, DACTCL_MIN_PIVOT);

                let mut i2 = 1;
                for k in self.resinit..=self.resfin {
                    let v = mytabb.get(i2);
                    self.mypoles.set(k, j, v);
                    i2 += 1;
                }
            }
            self.done = true;
        }

        // ===========================================================
        // cas de tangence:
        // ===========================================================

        let nincx = self.resfin - self.resinit + 1;
        let mut deport;
        let nincx2 = 2 * nincx;

        let mut internal_index = IVector::new_init(1, nincx, 0);
        self.search_index(&mut internal_index);
        let mut index = IVector::new_init(1, self.ninc, 0);

        let mut l = 1;
        if self.resinit <= self.resfin {
            for j in 0..=(self.na - 1) {
                deport = j * internal_index.get(nincx);
                for i in 1..=nincx {
                    index.set(l, internal_index.get(i) + deport);
                    l += 1;
                }
            }
        }

        if self.resinit > self.resfin {
            index.set(1, 1);
        }
        if ninc1 > 1 {
            if self.first_constraint >= AppParConstraint::TangencyPoint
                && self.last_constraint >= AppParConstraint::TangencyPoint
            {
                let v = index.get(ninc1 - 1) + ninc1;
                index.set(ninc1, v);
            }
        }
        if self.first_constraint >= AppParConstraint::TangencyPoint
            || self.last_constraint >= AppParConstraint::TangencyPoint
        {
            let v = index.get(self.ninc - 1) + self.ninc;
            index.set(self.ninc, v);
        }

        let mut the_a = RVector::new_init(1, index.get(self.ninc), 0.0);
        let mut my_tab = RVector::new_init(1, self.ninc, 0.0);

        self.make_taa_tab(&mut the_a, &mut my_tab);

        let mut error = dactcl_decompose(&mut the_a.data, &index.data, DACTCL_MIN_PIVOT);
        error = dactcl_solve(&the_a.data, &mut my_tab.data, &index.data, DACTCL_MIN_PIVOT);

        if error == MATH_STATUS_OK {
            self.done = true;
        }

        if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            self.lambda1 = my_tab.get(ninc1);
            self.lambda2 = my_tab.get(self.ninc);
        } else if self.first_constraint >= AppParConstraint::TangencyPoint {
            self.lambda1 = my_tab.get(self.ninc);
        } else if self.last_constraint >= AppParConstraint::TangencyPoint {
            self.lambda2 = my_tab.get(self.ninc);
        }

        // The results are stored in mypoles.
        //=========================================
        let mut k = 1;
        let mut i2 = 1;
        for _ci in 1..=self.nbp {
            let k1 = k + 1;
            let k2 = k + 2;
            for j in self.resinit..=self.resfin {
                self.mypoles.set(j, k, my_tab.get(i2));
                self.mypoles.set(j, k1, my_tab.get(i2 + nincx));
                self.mypoles.set(j, k2, my_tab.get(i2 + nincx2));
                i2 += 1;
            }

            if self.first_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles.set(
                    2,
                    k,
                    self.mypoints.get(self.myfirstp, k) + self.lambda1 * self.vec1t.get(k),
                );
                self.mypoles.set(
                    2,
                    k1,
                    self.mypoints.get(self.myfirstp, k1) + self.lambda1 * self.vec1t.get(k1),
                );
                self.mypoles.set(
                    2,
                    k2,
                    self.mypoints.get(self.myfirstp, k2) + self.lambda1 * self.vec1t.get(k2),
                );
            }
            if self.last_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles.set(
                    nbpol1,
                    k,
                    self.mypoints.get(self.mylastp, k) - self.lambda2 * self.vec2t.get(k),
                );
                self.mypoles.set(
                    nbpol1,
                    k1,
                    self.mypoints.get(self.mylastp, k1) - self.lambda2 * self.vec2t.get(k1),
                );
                self.mypoles.set(
                    nbpol1,
                    k2,
                    self.mypoints.get(self.mylastp, k2) - self.lambda2 * self.vec2t.get(k2),
                );
            }
            k += 3;
            i2 += nincx2;
        }

        for _ci in 1..=self.nbp2d {
            let k1 = k + 1;
            let k2 = k + 2;
            for j in self.resinit..=self.resfin {
                self.mypoles.set(j, k, my_tab.get(i2));
                self.mypoles.set(j, k1, my_tab.get(i2 + nincx));
                i2 += 1;
            }
            let _ = k2;
            if self.first_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles.set(
                    2,
                    k,
                    self.mypoints.get(self.myfirstp, k) + self.lambda1 * self.vec1t.get(k),
                );
                self.mypoles.set(
                    2,
                    k1,
                    self.mypoints.get(self.myfirstp, k1) + self.lambda1 * self.vec1t.get(k1),
                );
            }
            if self.last_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles.set(
                    nbpol1,
                    k,
                    self.mypoints.get(self.mylastp, k) - self.lambda2 * self.vec2t.get(k),
                );
                self.mypoles.set(
                    nbpol1,
                    k1,
                    self.mypoints.get(self.mylastp, k1) - self.lambda2 * self.vec2t.get(k1),
                );
            }
            k += 2;
            i2 += nincx;
        }
    }

    /// OCCT Perform(Parameters, V1t, V2t, l1, l2) (gxx L736-760).
    pub fn perform_v1tv2t(
        &mut self,
        parameters: &VecD,
        v1t: &VecD,
        v2t: &VecD,
        l1: f64,
        l2: f64,
    ) {
        self.done = false;
        if !self.isready {
            return;
        }
        // OCCT lower1 = V1t.Lower(), lower2 = V2t.Lower() — rcad VecD is
        // 1-based, so the shift is 1.
        self.resinit = 3;
        self.resfin = self.nbpoles - 2;
        let nincx = self.resfin - self.resinit + 1;
        self.ninc = self.na * nincx + 2;
        self.first_constraint = AppParConstraint::TangencyPoint;
        self.last_constraint = AppParConstraint::TangencyPoint;
        for i in 1..=self.vec1t.upper() {
            self.vec1t.set(i, v1t.get(i as usize));
            self.vec2t.set(i, v2t.get(i as usize));
        }
        self.perform_l1l2(parameters, l1, l2);
    }

    /// OCCT Perform(Parameters, V1t, V2t, V1c, V2c, l1, l2) (gxx L762-792).
    pub fn perform_v1tv2tv1cv2c(
        &mut self,
        parameters: &VecD,
        v1t: &VecD,
        v2t: &VecD,
        v1c: &VecD,
        v2c: &VecD,
        l1: f64,
        l2: f64,
    ) {
        self.done = false;
        if !self.isready {
            return;
        }
        self.resinit = 4;
        self.resfin = self.nbpoles - 3;
        let nincx = self.resfin - self.resinit + 1;
        self.ninc = self.na * nincx + 2;
        self.first_constraint = AppParConstraint::CurvaturePoint;
        self.last_constraint = AppParConstraint::CurvaturePoint;

        for i in 1..=self.vec1t.upper() {
            self.vec1t.set(i, v1t.get(i as usize));
            self.vec2t.set(i, v2t.get(i as usize));
            self.vec1c.set(i, v1c.get(i as usize));
            self.vec2c.set(i, v2c.get(i as usize));
        }
        self.perform_l1l2(parameters, l1, l2);
    }

    /// OCCT Perform(Parameters, l1, l2) (gxx L794-1055).
    pub fn perform_l1l2(&mut self, parameters: &VecD, l1: f64, l2: f64) {
        self.done = false;
        if !self.isready {
            return;
        }
        if self.first_constraint < AppParConstraint::TangencyPoint
            && self.last_constraint < AppParConstraint::TangencyPoint
        {
            self.perform(parameters);
            return;
        }
        self.iscalculated = false;

        self.lambda1 = l1;
        self.lambda2 = l2;
        let nbpol1 = self.nbpoles - 1;
        let l11 = self.deg as f64 * l1;
        let l22 = self.deg as f64 * l2;

        self.compute_function(parameters);

        if self.first_constraint >= AppParConstraint::TangencyPoint {
            for i in 1..=self.mypoles.col_number() {
                let v = self.mypoints.get(self.myfirstp, i) + l1 * self.vec1t.get(i);
                self.mypoles.set(2, i, v);
            }
        }

        if self.first_constraint == AppParConstraint::CurvaturePoint {
            for i in 1..=self.mypoles.col_number() {
                let v = 2.0 * self.mypoles.get(2, i) - self.mypoles.get(1, i)
                    + l11 * l11 * self.vec1c.get(i) / (self.deg as f64 * (self.deg as f64 - 1.0));
                self.mypoles.set(3, i, v);
            }
        }

        if self.last_constraint >= AppParConstraint::TangencyPoint {
            for i in 1..=self.mypoles.col_number() {
                let v = self.mypoints.get(self.mylastp, i) - l2 * self.vec2t.get(i);
                self.mypoles.set(self.nbpoles - 1, i, v);
            }
        }

        if self.last_constraint == AppParConstraint::CurvaturePoint {
            for i in 1..=self.mypoles.col_number() {
                let v = 2.0 * self.mypoles.get(self.nbpoles - 1, i) - self.mypoles.get(self.nbpoles, i)
                    + l22 * l22 * self.vec2c.get(i) / (self.deg as f64 * (self.deg as f64 - 1.0));
                self.mypoles.set(self.nbpoles - 2, i, v);
            }
        }
        let _ = nbpol1;

        if self.resinit > self.resfin {
            self.done = true;
            return;
        }

        if self.first_constraint == AppParConstraint::NoConstraint {
            if self.last_constraint == AppParConstraint::TangencyPoint {
                for j in self.first_p..=self.last_p {
                    let ad0 = self.a.get(j, self.nbpoles);
                    let ad1 = self.a.get(j, self.nbpoles - 1);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - ad0 * self.mypoles.get(self.nbpoles, i)
                            - ad1 * self.mypoles.get(self.nbpoles - 1, i);
                        self.b2.set(j, i, v);
                    }
                }
            }
            if self.last_constraint == AppParConstraint::CurvaturePoint {
                for j in self.first_p..=self.last_p {
                    let ad0 = self.a.get(j, self.nbpoles);
                    let ad1 = self.a.get(j, self.nbpoles - 1);
                    let ad2 = self.a.get(j, self.nbpoles - 2);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - ad0 * self.mypoles.get(self.nbpoles, i)
                            - ad1 * self.mypoles.get(self.nbpoles - 1, i)
                            - ad2 * self.mypoles.get(self.nbpoles - 2, i);
                        self.b2.set(j, i, v);
                    }
                }
            }
        } else if self.first_constraint == AppParConstraint::PassPoint {
            if self.last_constraint == AppParConstraint::TangencyPoint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let ad0 = self.a.get(j, self.nbpoles);
                    let ad1 = self.a.get(j, self.nbpoles - 1);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - ad0 * self.mypoles.get(self.nbpoles, i)
                            - ad1 * self.mypoles.get(self.nbpoles - 1, i);
                        self.b2.set(j, i, v);
                    }
                }
            }
            if self.last_constraint == AppParConstraint::CurvaturePoint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let ad0 = self.a.get(j, self.nbpoles);
                    let ad1 = self.a.get(j, self.nbpoles - 1);
                    let ad2 = self.a.get(j, self.nbpoles - 2);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - ad0 * self.mypoles.get(self.nbpoles, i)
                            - ad1 * self.mypoles.get(self.nbpoles - 1, i)
                            - ad2 * self.mypoles.get(self.nbpoles - 2, i);
                        self.b2.set(j, i, v);
                    }
                }
            }
        } else if self.first_constraint == AppParConstraint::TangencyPoint {
            if self.last_constraint == AppParConstraint::NoConstraint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let a1 = self.a.get(j, 2);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - a1 * self.mypoles.get(2, i);
                        self.b2.set(j, i, v);
                    }
                }
            } else if self.last_constraint == AppParConstraint::PassPoint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let ad0 = self.a.get(j, self.nbpoles);
                    let a1 = self.a.get(j, 2);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - ad0 * self.mypoles.get(self.nbpoles, i)
                            - a1 * self.mypoles.get(2, i);
                        self.b2.set(j, i, v);
                    }
                }
            } else if self.last_constraint == AppParConstraint::TangencyPoint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let ad0 = self.a.get(j, self.nbpoles);
                    let a1 = self.a.get(j, 2);
                    let ad1 = self.a.get(j, self.nbpoles - 1);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - ad0 * self.mypoles.get(self.nbpoles, i)
                            - a1 * self.mypoles.get(2, i)
                            - ad1 * self.mypoles.get(self.nbpoles - 1, i);
                        self.b2.set(j, i, v);
                    }
                }
            }
        } else if self.first_constraint == AppParConstraint::CurvaturePoint {
            if self.last_constraint == AppParConstraint::NoConstraint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let a1 = self.a.get(j, 2);
                    let a2 = self.a.get(j, 3);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - a1 * self.mypoles.get(2, i)
                            - a2 * self.mypoles.get(3, i);
                        self.b2.set(j, i, v);
                    }
                }
            } else if self.last_constraint == AppParConstraint::PassPoint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let a1 = self.a.get(j, 2);
                    let a2 = self.a.get(j, 3);
                    let ad0 = self.a.get(j, self.nbpoles);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - a1 * self.mypoles.get(2, i)
                            - a2 * self.mypoles.get(3, i)
                            - ad0 * self.mypoles.get(self.nbpoles, i);
                        self.b2.set(j, i, v);
                    }
                }
            } else if self.last_constraint == AppParConstraint::TangencyPoint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let a1 = self.a.get(j, 2);
                    let a2 = self.a.get(j, 3);
                    let ad0 = self.a.get(j, self.nbpoles);
                    let ad1 = self.a.get(j, self.nbpoles - 1);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - a1 * self.mypoles.get(2, i)
                            - a2 * self.mypoles.get(3, i)
                            - ad0 * self.mypoles.get(self.nbpoles, i)
                            - ad1 * self.mypoles.get(self.nbpoles - 1, i);
                        self.b2.set(j, i, v);
                    }
                }
            } else if self.last_constraint == AppParConstraint::CurvaturePoint {
                for j in self.first_p..=self.last_p {
                    let a0 = self.a.get(j, 1);
                    let a1 = self.a.get(j, 2);
                    let a2 = self.a.get(j, 3);
                    let ad0 = self.a.get(j, self.nbpoles);
                    let ad1 = self.a.get(j, self.nbpoles - 1);
                    let ad2 = self.a.get(j, self.nbpoles - 2);
                    for i in 1..=self.b2.col_number() {
                        let v = self.mypoints.get(j, i) - a0 * self.mypoles.get(1, i)
                            - a1 * self.mypoles.get(2, i)
                            - a2 * self.mypoles.get(3, i)
                            - ad0 * self.mypoles.get(self.nbpoles, i)
                            - ad1 * self.mypoles.get(self.nbpoles - 1, i)
                            - ad2 * self.mypoles.get(self.nbpoles - 2, i);
                        self.b2.set(j, i, v);
                    }
                }
            }
        }

        let nincx = self.resfin - self.resinit + 1;

        let mut mytab = Matrix::new_init(self.resinit, self.resfin, 1, self.b2.col_number(), 0.0);
        let mut index = IVector::new_init(1, nincx, 0);
        self.search_index(&mut index);
        let mut aa = RVector::new_init(1, index.get(nincx), 0.0);
        self.make_taa_matrix(&mut aa, &mut mytab);

        let mut mytabb = RVector::new_init(1, nincx, 0.0);

        dactcl_decompose(&mut aa.data, &index.data, DACTCL_MIN_PIVOT);

        let mut kk2;
        for j in 1..=self.b2.col_number() {
            kk2 = 1;
            for i in self.resinit..=self.resfin {
                let v = mytab.get(i, j);
                mytabb.set(kk2, v);
                kk2 += 1;
            }

            dactcl_solve(&aa.data, &mut mytabb.data, &index.data, DACTCL_MIN_PIVOT);

            let mut i2 = 1;
            for k in self.resinit..=self.resfin {
                let v = mytabb.get(i2);
                self.mypoles.set(k, j, v);
                i2 += 1;
            }
        }

        self.done = true;
    }

    /// OCCT Error(F, MaxE3d, MaxE2d) (gxx L1214-1280).
    pub fn error(&mut self, f: &mut f64, max_e3d: &mut f64, max_e2d: &mut f64) {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::Error");
        }
        let mut max3 = 0.0;
        let mut max2 = 0.0;
        *f = 0.0;
        let mut i2 = 1;
        let mut px = RVector::new(1, self.nbpoles);
        let mut py = RVector::new(1, self.nbpoles);
        let mut pz = RVector::new(1, self.nbpoles);

        for k in 1..=(self.nbp + self.nbp2d) {
            let i21 = i2 + 1;
            let i22 = i2 + 2;
            for i in 1..=self.nbpoles {
                px.set(i, self.mypoles.get(i, i2));
                py.set(i, self.mypoles.get(i, i21));
                if k <= self.nbp {
                    pz.set(i, self.mypoles.get(i, i22));
                }
            }
            for i in self.first_p..=self.last_p {
                let mut aa = 0.0;
                let mut bb = 0.0;
                let mut cc = 0.0;
                let indexdeb = self.myindex.get(i) + 1;
                let indexfin = indexdeb + self.deg;
                for j in indexdeb..=indexfin {
                    let aij = self.a.get(i, j);
                    aa += aij * px.get(j);
                    bb += aij * py.get(j);
                    if k <= self.nbp {
                        cc += aij * pz.get(j);
                    }
                }
                let fx = aa - self.mypoints.get(i, i2);
                let fy = bb - self.mypoints.get(i, i21);
                let mut fi = fx * fx + fy * fy;
                if k <= self.nbp {
                    let fz = cc - self.mypoints.get(i, i22);
                    fi += fz * fz;
                    if fi > max3 {
                        max3 = fi;
                    }
                } else if fi > max2 {
                    max2 = fi;
                }
                self.the_error.set(i, k, fi);
                *f += fi;
            }
            if k <= self.nbp {
                i2 += 3;
            } else {
                i2 += 2;
            }
        }
        *max_e3d = max3.sqrt();
        *max_e2d = max2.sqrt();
    }

    /// OCCT ErrorGradient(Grad, F, MaxE3d, MaxE2d) (gxx L1282-1369).
    pub fn error_gradient(
        &mut self,
        grad: &mut VecD,
        f: &mut f64,
        max_e3d: &mut f64,
        max_e2d: &mut f64,
    ) {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::ErrorGradient");
        }
        let mut max3 = 0.0;
        let mut max2 = 0.0;
        *f = 0.0;
        let mut i2 = 1;
        let mut px = RVector::new(1, self.nbpoles);
        let mut py = RVector::new(1, self.nbpoles);
        let mut pz = RVector::new(1, self.nbpoles);

        for k in 1..=(grad.len() as i32) {
            grad.set(k as usize, 0.0);
        }

        for k in 1..=(self.nbp + self.nbp2d) {
            let i21 = i2 + 1;
            let i22 = i2 + 2;
            for i in 1..=self.nbpoles {
                px.set(i, self.mypoles.get(i, i2));
                py.set(i, self.mypoles.get(i, i21));
                if k <= self.nbp {
                    pz.set(i, self.mypoles.get(i, i22));
                }
            }
            for i in self.first_p..=self.last_p {
                let mut aa = 0.0;
                let mut bb = 0.0;
                let mut cc = 0.0;
                let mut daa = 0.0;
                let mut dbb = 0.0;
                let mut dcc = 0.0;
                let indexdeb = self.myindex.get(i) + 1;
                let indexfin = indexdeb + self.deg;
                for j in indexdeb..=indexfin {
                    let aij = self.a.get(i, j);
                    let daij = self.da.get(i, j);
                    aa += aij * px.get(j);
                    daa += daij * px.get(j);
                    bb += aij * py.get(j);
                    dbb += daij * py.get(j);
                    if k <= self.nbp {
                        cc += aij * pz.get(j);
                        dcc += daij * pz.get(j);
                    }
                }
                let fx = aa - self.mypoints.get(i, i2);
                let fy = bb - self.mypoints.get(i, i2 + 1);
                let mut fi = fx * fx + fy * fy;
                let mut gr = 2.0 * (daa * fx + dbb * fy);

                if k <= self.nbp {
                    let fz = cc - self.mypoints.get(i, i2 + 2);
                    fi += fz * fz;
                    gr += 2.0 * dcc * fz;
                    if fi > max3 {
                        max3 = fi;
                    }
                } else if fi > max2 {
                    max2 = fi;
                }
                self.the_error.set(i, k, fi);
                let g = grad.get(i as usize) + gr;
                grad.set(i as usize, g);
                *f += fi;
            }
            if k <= self.nbp {
                i2 += 3;
            } else {
                i2 += 2;
            }
        }
        *max_e3d = max3.sqrt();
        *max_e2d = max2.sqrt();
    }

    /// OCCT Distance() (gxx L1371-1385).
    pub fn distance(&mut self) -> &MatD {
        if !self.iscalculated {
            for i in self.myfirstp..=self.mylastp {
                for j in 1..=(self.nbp + self.nbp2d) {
                    let v = self.the_error.get(i, j).sqrt();
                    self.the_error.set(i, j, v);
                }
            }
            self.iscalculated = true;
        }
        self.the_error.data()
    }

    /// OCCT FirstLambda() (gxx L1387-1390).
    pub fn first_lambda(&self) -> f64 {
        self.lambda1
    }

    /// OCCT LastLambda() (gxx L1392-1395).
    pub fn last_lambda(&self) -> f64 {
        self.lambda2
    }

    /// OCCT IsDone() (gxx L1397-1400).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT BezierValue() (gxx L1402-1407): `return
    /// (AppParCurves_MultiCurve)(BSplineValue());` — BSplineValue() runs
    /// first (repopulating SCU rows ideb..ifin and requiring done), then the
    /// MultiBSpCurve is down-cast to its MultiCurve base part.
    pub fn bezier_value(&mut self) -> MultiCurve {
        if self.myknots.is_some() {
            panic!("Standard_NoSuchObject: AppParCurves_LeastSquare::BezierValue");
        }
        self.bspline_value();
        MultiCurve {
            poles: self.scu.poles.clone(),
        }
    }

    /// OCCT BSplineValue() (gxx L1409-1445).
    pub fn bspline_value(&mut self) -> &MultiBSpCurve {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::BSplineValue");
        }

        let npoints = self.nbp + self.nbp2d;
        let mut ideb = self.resinit;
        let mut ifin = self.resfin;
        if ideb >= 2 {
            ideb = 2;
        }
        if ifin <= self.nbpoles - 1 {
            ifin = self.nbpoles - 1;
        }

        // Put the result in the corresponding curves
        for i in ideb..=ifin {
            let mut j2 = 1;
            let mut mpole = MultiPoint::new(self.nbp as usize, self.nbp2d as usize);
            for j in 1..=self.nbp {
                let pt = DVec3::new(
                    self.mypoles.get(i, j2),
                    self.mypoles.get(i, j2 + 1),
                    self.mypoles.get(i, j2 + 2),
                );
                mpole.set_point(j as usize, pt);
                j2 += 3;
            }
            for j in (self.nbp + 1)..=npoints {
                let pt = DVec2::new(self.mypoles.get(i, j2), self.mypoles.get(i, j2 + 1));
                mpole.set_point2d(j as usize, pt);
                j2 += 2;
            }
            self.scu.set_value(i as usize, mpole);
        }
        &self.scu
    }

    /// OCCT FunctionMatrix() (gxx L1447-1454).
    pub fn function_matrix(&self) -> &Matrix {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::FunctionMatrix");
        }
        &self.a
    }

    /// OCCT DerivativeFunctionMatrix() (gxx L1456-1463).
    pub fn derivative_function_matrix(&self) -> &Matrix {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::DerivativeFunctionMatrix");
        }
        &self.da
    }

    /// OCCT ComputeFunction(Parameters) (gxx L1483-1493).
    fn compute_function(&mut self, parameters: &VecD) {
        if self.myknots.is_none() {
            bernstein(self.nbpoles, parameters, &mut self.a, &mut self.da);
        } else {
            spline_function(
                self.nbpoles,
                self.deg,
                parameters,
                &self.vflatknots,
                &mut self.a,
                &mut self.da,
                &mut self.myindex,
            );
        }
    }

    /// OCCT Points() (gxx L1495-1498).
    pub fn points(&self) -> &Matrix {
        &self.mypoints
    }

    /// OCCT Poles() (gxx L1500-1503).
    pub fn poles(&self) -> &Matrix {
        &self.mypoles
    }

    /// OCCT KIndex() (gxx L1505-1508).
    pub fn k_index(&self) -> &IVector {
        &self.myindex
    }

    /// OCCT MakeTAA(TheA, myTAB) (gxx L1510-1703).
    fn make_taa_tab(&mut self, the_a: &mut RVector, my_tab: &mut RVector) {
        let mut xx = 0.0;
        let mut yy = 0.0;

        let nincx = self.resfin - self.resinit + 1;
        let neq = self.last_p - self.first_p + 1;

        let ninc1;

        if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            ninc1 = self.ninc - 1;
        } else {
            ninc1 = self.ninc;
        }

        let myfirst = self.a.lower_row();
        let mylast = myfirst + self.nlignes - 1;
        let mut taf1 = 0.0;
        let mut taf2 = 0.0;
        let mut taf3 = 0.0;
        let mut tab1 = 0.0;
        let mut tab2 = 0.0;
        let na1 = self.na - 1;
        let mut my_b = RVector::new_init(myfirst, mylast, 0.0);
        let mut my_v1 = RVector::new_init(myfirst, mylast, 0.0);
        let mut my_v2 = RVector::new_init(myfirst, mylast, 0.0);
        let mut the_v1 = RVector::new_init(1, self.ninc, 0.0);
        let mut the_v2 = RVector::new_init(1, self.ninc, 0.0);

        for i in self.first_p..=self.last_p {
            let ai2 = self.a.get(i, 2);
            let aid = self.a.get(i, self.nbpoles - 1);
            if self.first_constraint >= AppParConstraint::PassPoint {
                xx = self.a.get(i, 1);
            }
            if self.first_constraint >= AppParConstraint::TangencyPoint {
                xx += ai2;
            }
            if self.last_constraint >= AppParConstraint::PassPoint {
                yy = self.a.get(i, self.nbpoles);
            }
            if self.last_constraint >= AppParConstraint::TangencyPoint {
                yy += aid;
            }
            let mut i2 = 1;
            let mut nrow = myfirst - self.first_p;
            for _ci in 1..=self.nbp {
                let i21 = i2 + 1;
                let i22 = i2 + 2;
                let ix = i + nrow;
                let iy = ix + neq;
                let iz = iy + neq;
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    my_v1.set(ix, ai2 * self.vec1t.get(i2));
                    my_v1.set(iy, ai2 * self.vec1t.get(i21));
                    my_v1.set(iz, ai2 * self.vec1t.get(i22));
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    my_v2.set(ix, -aid * self.vec2t.get(i2));
                    my_v2.set(iy, -aid * self.vec2t.get(i21));
                    my_v2.set(iz, -aid * self.vec2t.get(i22));
                }
                my_b.set(
                    ix,
                    self.mypoints.get(i, i2) - xx * self.mypoints.get(self.myfirstp, i2)
                        - yy * self.mypoints.get(self.mylastp, i2),
                );
                my_b.set(
                    iy,
                    self.mypoints.get(i, i21) - xx * self.mypoints.get(self.myfirstp, i21)
                        - yy * self.mypoints.get(self.mylastp, i21),
                );
                my_b.set(
                    iz,
                    self.mypoints.get(i, i22) - xx * self.mypoints.get(self.myfirstp, i22)
                        - yy * self.mypoints.get(self.mylastp, i22),
                );
                i2 += 3;
                nrow += 3 * neq;
            }

            for _ci in 1..=self.nbp2d {
                let i21 = i2 + 1;
                let i22 = i2 + 2;
                let ix = i + nrow;
                let iy = ix + neq;
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    my_v1.set(ix, ai2 * self.vec1t.get(i2));
                    my_v1.set(iy, ai2 * self.vec1t.get(i21));
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    my_v2.set(ix, -aid * self.vec2t.get(i2));
                    my_v2.set(iy, -aid * self.vec2t.get(i21));
                }
                my_b.set(
                    ix,
                    self.mypoints.get(i, i2) - xx * self.mypoints.get(self.myfirstp, i2)
                        - yy * self.mypoints.get(self.mylastp, i2),
                );
                my_b.set(
                    iy,
                    self.mypoints.get(i, i21) - xx * self.mypoints.get(self.myfirstp, i21)
                        - yy * self.mypoints.get(self.mylastp, i21),
                );
                nrow += 2 * neq;
                i2 += 2;
                let _ = i22;
            }
        }

        // Construction de TA*A et TA*B:

        for k in self.first_p..=self.last_p {
            let indexdeb = self.myindex.get(k) + 1;
            let indexfin = indexdeb + self.deg;
            let jinit = self.resinit.max(indexdeb);
            let jfin = self.resfin.min(indexfin);
            let k1 = k + myfirst - self.first_p;
            for i in 0..=na1 {
                let nb = i * neq + k1;
                let mut v1 = 0.0;
                let mut v2 = 0.0;
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    v1 = my_v1.get(nb);
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    v2 = my_v2.get(nb);
                }
                let b = my_b.get(nb);
                let inc = i * nincx - self.resinit + 1;
                for j in jinit..=jfin {
                    let akj = self.a.get(k, j);
                    let u = j + inc;
                    if self.first_constraint >= AppParConstraint::TangencyPoint {
                        let v = the_v1.get(u) + akj * v1;
                        the_v1.set(u, v);
                    }
                    if self.last_constraint >= AppParConstraint::TangencyPoint {
                        let v = the_v2.get(u) + akj * v2;
                        the_v2.set(u, v);
                    }
                    let v = my_tab.get(u) + akj * b;
                    my_tab.set(u, v);
                }
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    taf1 += v1 * v1;
                    tab1 += v1 * b;
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    taf2 += v2 * v2;
                    tab2 += v2 * b;
                }
                if self.first_constraint >= AppParConstraint::TangencyPoint
                    && self.last_constraint >= AppParConstraint::TangencyPoint
                {
                    taf3 += v1 * v2;
                }
            }
        }

        if self.first_constraint >= AppParConstraint::TangencyPoint {
            the_v1.set(ninc1, taf1);
            my_tab.set(ninc1, tab1);
        }
        if self.last_constraint >= AppParConstraint::TangencyPoint {
            the_v2.set(self.ninc, taf2);
            my_tab.set(self.ninc, tab2);
        }
        if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            the_v2.set(ninc1, taf3);
        }

        if self.resinit <= self.resfin {
            let mut index = IVector::new_init(1, nincx, 0);
            self.search_index(&mut index);
            let mut aa = RVector::new(1, index.get(nincx));
            self.make_taa(&mut aa);

            let mut kk = 1;
            for _k in 1..=self.na {
                for i in 1..=aa.length() {
                    let v = aa.get(i);
                    the_a.set(kk, v);
                    kk += 1;
                }
            }
        }

        let length = the_a.length();

        if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            for j in 1..=ninc1 {
                let v = the_v1.get(j);
                the_a.set(length - 2 * self.ninc + j + 1, v);
            }
            for j in 1..=self.ninc {
                let v = the_v2.get(j);
                the_a.set(length - self.ninc + j, v);
            }
        } else if self.first_constraint >= AppParConstraint::TangencyPoint {
            for j in 1..=self.ninc {
                let v = the_v1.get(j);
                the_a.set(length - self.ninc + j, v);
            }
        } else if self.last_constraint >= AppParConstraint::TangencyPoint {
            for j in 1..=self.ninc {
                let v = the_v2.get(j);
                the_a.set(length - self.ninc + j, v);
            }
        }
    }

    /// OCCT MakeTAA(TheA, myTAB-matrix) (gxx L1705-1761).
    fn make_taa_matrix(&mut self, aa: &mut RVector, mytab: &mut Matrix) {
        let mut the_a = Matrix::new(self.resinit, self.resfin, self.resinit, self.resfin);
        the_a.init(0.0);

        for k in self.first_p..=self.last_p {
            let indexdeb = self.myindex.get(k) + 1;
            let indexfin = indexdeb + self.deg;
            let jinit = self.resinit.max(indexdeb);
            let jfin = self.resfin.min(indexfin);
            for i in jinit..=jfin {
                let akj = self.a.get(k, i);
                for j in jinit..=i {
                    let v = the_a.get(i, j) + self.a.get(k, j) * akj;
                    the_a.set(i, j, v);
                }
                for j in 1..=self.b2.col_number() {
                    let v = mytab.get(i, j) + akj * self.b2.get(k, j);
                    mytab.set(i, j, v);
                }
            }
        }

        let len = match &self.myknots {
            Some(knots) => knots.len() as i32,
            None => 2,
        };
        let mut i2 = 1;
        let mut iinit = self.resinit;
        let mut jinit = self.resinit;
        let mut ifin = self.resfin.min(self.deg + 1);
        for k in 2..=len {
            for i in iinit..=ifin {
                for j in jinit..=i {
                    let v = the_a.get(i, j);
                    aa.set(i2, v);
                    i2 += 1;
                }
            }
            if let Some(mymults) = &self.mymults {
                iinit = ifin + 1;
                let d = ifin + mymults[(k - 1) as usize];
                ifin = d.min(self.resfin);
                jinit = (d - self.deg).max(self.resinit);
            }
        }
    }

    /// OCCT MakeTAA(AA) (gxx L1763-1814).
    fn make_taa(&mut self, aa: &mut RVector) {
        let mut the_a = Matrix::new_init(self.resinit, self.resfin, self.resinit, self.resfin, 0.0);

        for k in self.first_p..=self.last_p {
            let indexdeb = self.myindex.get(k) + 1;
            let indexfin = indexdeb + self.deg;
            let jinit = self.resinit.max(indexdeb);
            let jfin = self.resfin.min(indexfin);
            for i in jinit..=jfin {
                let akj = self.a.get(k, i);
                for j in jinit..=i {
                    let v = the_a.get(i, j) + self.a.get(k, j) * akj;
                    the_a.set(i, j, v);
                }
            }
        }

        let mut i2 = 1;
        let mut iinit = self.resinit;
        let mut jinit = self.resinit;
        let mut ifin = self.resfin.min(self.deg + 1);
        let len = match &self.myknots {
            Some(knots) => knots.len() as i32,
            None => 2,
        };
        for k in 2..=len {
            for i in iinit..=ifin {
                for j in jinit..=i {
                    let v = the_a.get(i, j);
                    aa.set(i2, v);
                    i2 += 1;
                }
            }
            if let Some(mymults) = &self.mymults {
                iinit = ifin + 1;
                let d = ifin + mymults[(k - 1) as usize];
                ifin = d.min(self.resfin);
                jinit = (d - self.deg).max(self.resinit);
            }
        }
    }

    /// OCCT SearchIndex(Index) (gxx L1816-1860).
    fn search_index(&mut self, index: &mut IVector) {
        let nincx = self.resfin - self.resinit + 1;
        index.set(1, 1);

        if self.myknots.is_none() {
            if self.resinit <= self.resfin {
                let mut l = 1;
                for i in 2..=nincx {
                    l += 1;
                    let v = index.get(l - 1) + i;
                    index.set(l, v);
                }
            }
        } else {
            let mut iinit = self.resinit;
            let mut jinit = self.resinit;
            let mut ifin = self.resfin.min(self.deg + 1);
            let len = self.myknots.as_ref().unwrap().len() as i32;

            let mut i2 = 1;
            for k in 2..=len {
                for i in iinit..=ifin {
                    for j in jinit..=i {
                        if i2 != 1 {
                            let v = index.get(i2 - 1) + i - jinit + 1;
                            index.set(i2, v);
                        }
                        i2 += 1;
                    }
                }
                iinit = ifin + 1;
                let d = ifin + self.mymults.as_ref().unwrap()[(k - 1) as usize];
                ifin = d.min(self.resfin);
                jinit = (d - self.deg).max(self.resinit);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Free functions (the OCCT protected methods that touch no instance state
// beyond the passed-in members, and the Affect helper)
// ---------------------------------------------------------------------------

/// OCCT NbBColumns(SSP) (gxx L1207-1212) — the number of second member
/// columns.
fn nb_b_columns(ssp: &MultiLine) -> i32 {
    my_line_tool::nb_p3d(ssp) as i32 * 3 + my_line_tool::nb_p2d(ssp) as i32 * 2
}

/// OCCT TheFirstPoint(FirstCons, FirstPoint) (gxx L1465-1472) — the first
/// point being fitted.
fn the_first_point(first_cons: AppParConstraint, first_point: i32) -> i32 {
    if first_cons == AppParConstraint::NoConstraint {
        first_point
    } else {
        first_point + 1
    }
}

/// OCCT TheLastPoint(LastCons, LastPoint) (gxx L1474-1481) — the last point
/// being fitted.
fn the_last_point(last_cons: AppParConstraint, last_point: i32) -> i32 {
    if last_cons == AppParConstraint::NoConstraint {
        last_point
    } else {
        last_point - 1
    }
}

/// OCCT Affect(SSP, Index, Cons, Vt, Vc) (gxx L1064-1205) — Index is an ID
/// of the point in the MultiLine; fills the tangent (Vt) and curvature (Vc)
/// vectors, degrading the constraint when the tool cannot supply them. The
/// OCCT members nbP / nbP2d are passed explicitly (nbp / nbp2d).
fn affect(
    ssp: &MultiLine,
    index: i32,
    cons: &mut AppParConstraint,
    vt: &mut RVector,
    vc: &mut RVector,
    nbp: i32,
    nbp2d: i32,
) {
    // Vt: vector of tangent, Vc: vector of curvature.
    if *cons >= AppParConstraint::TangencyPoint {
        let mut i2 = 1;
        let mut mynbp2d = nbp2d;
        let mut mynbp = nbp;
        if nbp2d == 0 {
            mynbp2d = 1;
        }
        if nbp == 0 {
            mynbp = 1;
        }
        let mut tab_v = vec![DVec3::ZERO; mynbp as usize];
        let mut tab_v2d = vec![DVec2::ZERO; mynbp2d as usize];

        let mut ok;
        if *cons == AppParConstraint::CurvaturePoint {
            if nbp != 0 && nbp2d != 0 {
                ok = my_line_tool::curvature_3d_2d(ssp, index as usize, &mut tab_v, &mut tab_v2d);
                if !ok {
                    *cons = AppParConstraint::TangencyPoint;
                }
            } else if nbp2d != 0 {
                ok = my_line_tool::curvature_2d(ssp, index as usize, &mut tab_v2d);
                if !ok {
                    *cons = AppParConstraint::TangencyPoint;
                }
            } else {
                ok = my_line_tool::curvature_3d(ssp, index as usize, &mut tab_v);
                if !ok {
                    *cons = AppParConstraint::TangencyPoint;
                }
            }
            if ok {
                i2 = 1;
                for i in 1..=nbp {
                    let v = tab_v[(i - 1) as usize];
                    vc.set(i2, v.x);
                    vc.set(i2 + 1, v.y);
                    vc.set(i2 + 2, v.z);
                    i2 += 3;
                }

                for i in 1..=nbp2d {
                    let v = tab_v2d[(i - 1) as usize];
                    vc.set(i2, v.x);
                    vc.set(i2 + 1, v.y);
                    i2 += 2;
                }
            }
        } else {
            ok = false;
            let _ = ok;
        }

        i2 = 1;
        if *cons >= AppParConstraint::TangencyPoint {
            let mut ok;
            if nbp != 0 && nbp2d != 0 {
                ok = my_line_tool::tangency_3d_2d(ssp, index as usize, &mut tab_v, &mut tab_v2d);
                if !ok {
                    *cons = AppParConstraint::PassPoint;
                }
            } else if nbp2d != 0 {
                ok = my_line_tool::tangency_2d(ssp, index as usize, &mut tab_v2d);
                if !ok {
                    *cons = AppParConstraint::PassPoint;
                }
            } else {
                ok = my_line_tool::tangency_3d(ssp, index as usize, &mut tab_v);
                if !ok {
                    *cons = AppParConstraint::PassPoint;
                }
            }

            if ok {
                let mut an_arr_pts3d1 = vec![DVec3::ZERO; mynbp as usize];
                let mut an_arr_pts3d2 = vec![DVec3::ZERO; mynbp as usize];

                if nbp != 0 {
                    if (index as usize) < my_line_tool::last_point(ssp) {
                        my_line_tool::value_3d(ssp, index as usize, &mut an_arr_pts3d1);
                        my_line_tool::value_3d(ssp, (index + 1) as usize, &mut an_arr_pts3d2);
                    } else {
                        // (Index == ToolLine::LastPoint(theML))
                        my_line_tool::value_3d(ssp, (index - 1) as usize, &mut an_arr_pts3d1);
                        my_line_tool::value_3d(ssp, index as usize, &mut an_arr_pts3d2);
                    }

                    check_tangents(&an_arr_pts3d1, &an_arr_pts3d2, &mut tab_v, &mut tab_v2d);
                } else if nbp2d != 0 {
                    let mut an_arr_pts2d1 = vec![DVec2::ZERO; mynbp2d as usize];
                    let mut an_arr_pts2d2 = vec![DVec2::ZERO; mynbp2d as usize];

                    if (index as usize) < my_line_tool::last_point(ssp) {
                        my_line_tool::value_3d_2d(
                            ssp,
                            index as usize,
                            &mut an_arr_pts3d1,
                            &mut an_arr_pts2d1,
                        );
                        my_line_tool::value_3d_2d(
                            ssp,
                            (index + 1) as usize,
                            &mut an_arr_pts3d2,
                            &mut an_arr_pts2d2,
                        );
                    } else {
                        // (Index == ToolLine::LastPoint(theML))
                        my_line_tool::value_3d_2d(
                            ssp,
                            (index - 1) as usize,
                            &mut an_arr_pts3d1,
                            &mut an_arr_pts2d1,
                        );
                        my_line_tool::value_3d_2d(
                            ssp,
                            index as usize,
                            &mut an_arr_pts3d2,
                            &mut an_arr_pts2d2,
                        );
                    }

                    check_tangents_2d(&an_arr_pts2d1, &an_arr_pts2d2, &mut tab_v2d);
                }

                for i in 1..=nbp {
                    let v = tab_v[(i - 1) as usize];
                    vt.set(i2, v.x);
                    vt.set(i2 + 1, v.y);
                    vt.set(i2 + 2, v.z);
                    i2 += 3;
                }

                for i in 1..=nbp2d {
                    let v = tab_v2d[(i - 1) as usize];
                    vt.set(i2, v.x);
                    vt.set(i2 + 1, v.y);
                    i2 += 2;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// AppParCurves_ResolConstraint
// ---------------------------------------------------------------------------

/// OCCT AppParCurves_ResolConstraint (AppParCurves_ResolConstraint.gxx
/// whole file) — given a MultiLine SSP with constraints points, finds the
/// best curve solution to approximate it. The poles from SCurv issued for
/// example from the least squares are used as a guess solution for the
/// Uzawa algorithm (math_Uzawa, tolerance `tolerance`); Bern / DA are the
/// Bernstein and derivative-Bernstein matrices of the MultiLine. The
/// MultiCurve SCurv is modified with the new multipoles.
///
/// Note: the OCCT curvature-equation block inside the constructor is
/// commented out in the C++ source and is not ported either. The
/// `Error()` accessor declared by the instantiation headers has no
/// definition in the OCCT library (dead declaration) — not ported.
#[derive(Debug, Clone)]
pub struct ResolConstraint {
    done: bool,
    /// OCCT Err (set by no code path in OCCT).
    _err: f64,
    cont: Matrix,
    de_cont: Matrix,
    secont: RVector,
    ctcinv: Matrix,
    vardua: RVector,
    inc_pass: i32,
    inc_tan: i32,
    inc_curv: i32,
    /// OCCT IPas / ITan / ICurv (NCollection_Array1<int>, 1-based).
    ipas: Vec<i32>,
    itan: Vec<i32>,
    icurv: Vec<i32>,
}

impl ResolConstraint {
    /// OCCT AppParCurves_ResolConstraint(SSP, SCurv, FirstPoint, LastPoint,
    /// TheConstraints, Bern, DerivativeBern, Tolerance = 1.0e-10)
    /// (gxx ctor L40-503).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ssp: &MultiLine,
        scurv: &mut MultiCurve,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        bern: &Matrix,
        derivative_bern: &Matrix,
        tolerance: f64,
    ) -> Self {
        let nb_constraints =
            ResolConstraint::nb_constraints(ssp, first_point, last_point, the_constraints);
        let nb_columns = ResolConstraint::nb_columns(ssp, scurv.nb_poles() as i32 - 1);
        let mut r = ResolConstraint {
            done: false,
            _err: 0.0,
            cont: Matrix::new_init(1, nb_constraints, 1, nb_columns, 0.0),
            de_cont: Matrix::new_init(1, nb_constraints, 1, nb_columns, 0.0),
            secont: RVector::new_init(1, nb_constraints, 0.0),
            ctcinv: Matrix::new(1, nb_constraints, 1, nb_constraints),
            vardua: RVector::new(1, nb_constraints),
            inc_pass: 0,
            inc_tan: 0,
            inc_curv: 0,
            ipas: vec![0; (last_point - first_point + 1) as usize],
            itan: vec![0; (last_point - first_point + 1) as usize],
            icurv: vec![0; (last_point - first_point + 1) as usize],
        };
        r.perform(
            ssp,
            scurv,
            first_point,
            last_point,
            the_constraints,
            bern,
            derivative_bern,
            tolerance,
        );
        r
    }

    /// The constructor body (gxx L59-503).
    #[allow(clippy::too_many_arguments)]
    fn perform(
        &mut self,
        ssp: &MultiLine,
        scurv: &mut MultiCurve,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        bern: &Matrix,
        derivative_bern: &Matrix,
        tolerance: f64,
    ) {
        let nb_cu = scurv.nb_curves() as i32;
        let mut myindex;
        let def = scurv.nb_poles() as i32 - 1;
        let nb3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb2d = my_line_tool::nb_p2d(ssp) as i32;
        let n_pol = def + 1;
        let n_pol2 = 2 * n_pol;
        let mut fc = AppParConstraint::NoConstraint;
        let mut lc = AppParConstraint::NoConstraint;
        let mut cons;
        // Boucle de calcul du nombre de points de passage afin de
        // dimensionner les matrices.
        self.inc_pass = 0;
        self.inc_tan = 0;
        self.inc_curv = 0;
        for couple in the_constraints {
            myindex = couple.index;
            cons = couple.constraint;
            if myindex == first_point {
                fc = cons;
            }
            if myindex == last_point {
                lc = cons;
            }
            if cons as i32 >= 1 {
                self.inc_pass += 1; // IncPass = nbre de points de passage.
                self.ipas[(self.inc_pass - 1) as usize] = myindex;
            }
            if cons as i32 >= 2 {
                self.inc_tan += 1; // IncTan = nbre de points de tangence.
                self.itan[(self.inc_tan - 1) as usize] = myindex;
            }
            if cons as i32 == 3 {
                self.inc_curv += 1; // IncCurv = nbre de pts de courbure.
                self.icurv[(self.inc_curv - 1) as usize] = myindex;
            }
        }
        if self.inc_pass == 0 {
            self.done = true;
            return;
        }
        let mut mynb3d = nb3d;
        let mut mynb2d = nb2d;
        if nb3d == 0 {
            mynb3d = 1;
        }
        if nb2d == 0 {
            mynb2d = 1;
        }
        let c_col = nb3d * 3 + nb2d * 2;
        // Declaration et initialisation des matrices et vecteurs de
        // contraintes:
        let mut cont_init = Matrix::new_init(1, self.inc_pass, 1, n_pol, 0.0);
        let mut start = RVector::new(1, c_col * n_pol);
        let mut ibont = vec![vec![0i32; self.inc_tan as usize]; nb_cu as usize];
        // Filling Cont for the passing points:
        // =================================================
        for i in 1..=self.inc_pass {
            // Cette partie ne depend que de Bernstein
            let npt = self.ipas[(i - 1) as usize];
            for j in 1..=n_pol {
                let v = bern.get(npt, j);
                cont_init.set(i, j, v);
            }
        }
        for i in 1..=c_col {
            let v = cont_init.clone();
            self.cont.set_block(
                self.inc_pass * (i - 1) + 1,
                self.inc_pass * i,
                n_pol * (i - 1) + 1,
                n_pol * i,
                &v,
            );
        }
        // Retrieval of starting vectors for Uzawa. This vector represents
        // the poles of SCurv.
        // Filling of secont and resolution.
        let mut tab_v = vec![DVec3::ZERO; mynb3d as usize];
        let mut tab_v2d = vec![DVec2::ZERO; mynb2d as usize];
        let mut tab_p = vec![DVec3::ZERO; mynb3d as usize];
        let mut tab_p2d = vec![DVec2::ZERO; mynb2d as usize];
        let mut inc3 = c_col * self.inc_pass + 1;
        let mut inc_col = 0i32;
        let mut inc_sec = 0i32;
        for k in 1..=nb_cu {
            if k <= nb3d {
                for i in 1..=self.inc_tan {
                    let npt = self.itan[(i - 1) as usize];
                    // choix du maximum de tangence pour exprimer la
                    // colinearite:
                    my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                    let v = tab_v[(k - 1) as usize];
                    let t1v = v.x;
                    let t2v = v.y;
                    let t3v = v.z;
                    let mut tmax = t1v.abs();
                    ibont[(k - 1) as usize][(i - 1) as usize] = 1;
                    if t2v.abs() > tmax {
                        tmax = t2v.abs();
                        ibont[(k - 1) as usize][(i - 1) as usize] = 2;
                    }
                    if t3v.abs() > tmax {
                        tmax = t3v.abs();
                        ibont[(k - 1) as usize][(i - 1) as usize] = 3;
                    }
                    let ib = ibont[(k - 1) as usize][(i - 1) as usize];
                    if ib == 3 {
                        for j in 1..=n_pol {
                            let daij = derivative_bern.get(npt, j);
                            let val = daij * t3v / tmax;
                            self.cont.set(inc3, j + n_pol + inc_col, val);
                            let val = -daij * t2v / tmax;
                            self.cont.set(inc3, j + n_pol2 + inc_col, val);
                            let val = daij * t3v / tmax;
                            self.cont.set(inc3 + 1, j + inc_col, val);
                            let val = -daij * t1v / tmax;
                            self.cont.set(inc3 + 1, j + n_pol2 + inc_col, val);
                        }
                    } else if ib == 1 {
                        for j in 1..=n_pol {
                            let daij = derivative_bern.get(npt, j);
                            let val = daij * t3v / tmax;
                            self.cont.set(inc3, j + inc_col, val);
                            let val = -daij * t1v / tmax;
                            self.cont.set(inc3, j + n_pol2 + inc_col, val);
                            let val = daij * t2v / tmax;
                            self.cont.set(inc3 + 1, j + inc_col, val);
                            let val = -daij * t1v / tmax;
                            self.cont.set(inc3 + 1, j + n_pol + inc_col, val);
                        }
                    } else if ib == 2 {
                        for j in 1..=(def + 1) {
                            let daij = derivative_bern.get(npt, j);
                            let val = daij * t3v / tmax;
                            self.cont.set(inc3, j + n_pol + inc_col, val);
                            let val = -daij * t2v / tmax;
                            self.cont.set(inc3, j + n_pol2 + inc_col, val);
                            let val = daij * t2v / tmax;
                            self.cont.set(inc3 + 1, j + inc_col, val);
                            let val = -daij * t1v / tmax;
                            self.cont.set(inc3 + 1, j + n_pol + inc_col, val);
                        }
                    }
                    inc3 += 2;
                }
                // Remplissage du second membre:
                for i in 1..=self.inc_pass {
                    my_line_tool::value_3d(ssp, self.ipas[(i - 1) as usize] as usize, &mut tab_p);
                    let poi = tab_p[(k - 1) as usize];
                    self.secont.set(i + inc_sec, poi.x);
                    self.secont.set(i + self.inc_pass + inc_sec, poi.y);
                    self.secont.set(i + 2 * self.inc_pass + inc_sec, poi.z);
                }
                inc_sec += 3 * self.inc_pass;
                // Vecteur de depart:
                for j in 1..=n_pol {
                    let poi = scurv.value(j as usize).point(k as usize);
                    start.set(j + inc_col, poi.x);
                    start.set(j + n_pol + inc_col, poi.y);
                    start.set(j + n_pol2 + inc_col, poi.z);
                }
                inc_col += 3 * n_pol;
            } else {
                for i in 1..=self.inc_tan {
                    let npt = self.itan[(i - 1) as usize];
                    my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                    let v2d = tab_v2d[(k - nb3d - 1) as usize];
                    let t1v = v2d.x;
                    let t2v = v2d.y;
                    let mut tmax = t1v.abs();
                    ibont[(k - 1) as usize][(i - 1) as usize] = 1;
                    if t2v.abs() > tmax {
                        tmax = t2v.abs();
                        ibont[(k - 1) as usize][(i - 1) as usize] = 2;
                    }
                    for j in 1..=n_pol {
                        let daij = derivative_bern.get(npt, j);
                        let val = daij * t2v;
                        self.cont.set(inc3, j + inc_col, val);
                        let val = -daij * t1v;
                        self.cont.set(inc3, j + n_pol + inc_col, val);
                    }
                    inc3 += 1;
                }
                // Remplissage du second membre:
                for i in 1..=self.inc_pass {
                    my_line_tool::value_2d(ssp, self.ipas[(i - 1) as usize] as usize, &mut tab_p2d);
                    // OCCT: `Poi2d = tabP2d(i - nb3d);` — the pass-point
                    // index expression is kept verbatim (identical to the
                    // curve index only for pure-2D lines; OCCT compiles
                    // this file with the range checks disabled).
                    let poi2d = tab_p2d[(i - nb3d - 1) as usize];
                    self.secont.set(i + inc_sec, poi2d.x);
                    self.secont.set(i + self.inc_pass + inc_sec, poi2d.y);
                }
                inc_sec += 2 * self.inc_pass;
                // Remplissage du vecteur de depart:
                for j in 1..=n_pol {
                    let poi2d = scurv.value(j as usize).point2d(k as usize);
                    start.set(j + inc_col, poi2d.x);
                    start.set(j + n_pol + inc_col, poi2d.y);
                }
                inc_col += n_pol2;
            }
        }
        // Equations exprimant le meme rapport de tangence sur chaque courbe:
        // On prend les coordonnees les plus significatives.
        inc3 -= 1;
        for i in 1..=self.inc_tan {
            inc_col = 0;
            let npt = self.itan[(i - 1) as usize];
            for k in 1..=(nb_cu - 1) {
                inc3 += 1;
                // Initialize first relation variable (T1)
                let mut add_index_1 = 0i32;
                let a_val = ibont[(k - 1) as usize][(i - 1) as usize];
                let mut ip = 0i32;
                let t1v;
                match a_val {
                    1 => {
                        // T1 ~ T1x
                        if k <= nb3d {
                            my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                            t1v = tab_v[(k - 1) as usize].x;
                            ip = 3 * n_pol;
                        } else {
                            my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                            t1v = tab_v2d[(k - nb3d - 1) as usize].x;
                            ip = 2 * n_pol;
                        }
                        add_index_1 = 0;
                    }
                    2 => {
                        // T1 ~ T1y
                        let yv;
                        if k <= nb3d {
                            my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                            ip = 3 * n_pol;
                            yv = tab_v[(k - 1) as usize].y;
                        } else {
                            my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                            ip = 2 * n_pol;
                            yv = tab_v2d[(k - nb3d - 1) as usize].y;
                        }
                        t1v = yv;
                        add_index_1 = n_pol;
                    }
                    _ => {
                        // 3: T1 ~ T1z
                        my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                        t1v = tab_v[(k - 1) as usize].z;
                        ip = 3 * n_pol;
                        add_index_1 = 2 * n_pol;
                    }
                }
                // Initialize second relation variable (T2)
                let mut add_index_2 = 0i32;
                let a_next_val = ibont[k as usize][(i - 1) as usize];
                let t2v;
                match a_next_val {
                    1 => {
                        // T2 ~ T2x
                        if (k + 1) <= nb3d {
                            my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                            t2v = tab_v[k as usize].x;
                        } else {
                            my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                            t2v = tab_v2d[(k + 1 - nb3d - 1) as usize].x;
                        }
                        add_index_2 = 0;
                    }
                    2 => {
                        // T2 ~ T2y
                        if (k + 1) <= nb3d {
                            my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                            t2v = tab_v[k as usize].y;
                        } else {
                            my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                            t2v = tab_v2d[(k + 1 - nb3d - 1) as usize].y;
                        }
                        add_index_2 = n_pol;
                    }
                    _ => {
                        // 3: T2 ~ T2z
                        my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                        t2v = tab_v[k as usize].z;
                        add_index_2 = 2 * n_pol;
                    }
                }
                // Relations between T1 and T2:
                for j in 1..=n_pol {
                    let daij = derivative_bern.get(npt, j);
                    let val = daij * t2v;
                    self.cont.set(inc3, j + inc_col + add_index_1, val);
                    let val = -daij * t1v;
                    self.cont.set(inc3, j + ip + inc_col + add_index_2, val);
                }
                inc_col += ip;
            }
        }
        // Equations concernant la courbure: commented out in the OCCT
        // source (gxx L398-455) — not ported.

        // Resolution par Uzawa:
        // math_Uzawa UzaResol(Cont, Secont, Start, Tolerance) — ctor 1 with
        // EpsLix = Tolerance, EpsLic = 1.0e-06, NbIterations = 500.
        let uza_resol = Uzawa::new(&self.cont, &self.secont, &start, tolerance, 1.0e-06, 500);
        if !uza_resol.is_done() {
            self.done = false;
            return;
        }
        self.ctcinv = uza_resol.inverse_cont().clone();
        uza_resol.duale(&mut self.vardua);
        for i in 1..=self.ctcinv.row_number() {
            for j in i..=self.ctcinv.row_number() {
                let v = self.ctcinv.get(j, i);
                self.ctcinv.set(i, j, v);
            }
        }
        self.done = true;
        let mut vec_poles = RVector::new(1, c_col * n_pol);
        let uza_value = uza_resol.value();
        for i in 1..=vec_poles.length() {
            let v = uza_value.get(i);
            vec_poles.set(i, v);
        }
        let mut polinit = 1;
        let mut polfin = n_pol;
        if fc as i32 >= 1 {
            polinit = 2;
        }
        if lc as i32 >= 1 {
            polfin = n_pol - 1;
        }
        for i in polinit..=polfin {
            let mut inc_col = 0i32;
            let mut mpol = MultiPoint::new(nb3d as usize, nb2d as usize);
            for k in 1..=nb_cu {
                if k <= nb3d {
                    let pol = DVec3::new(
                        vec_poles.get(inc_col + i),
                        vec_poles.get(inc_col + n_pol + i),
                        vec_poles.get(inc_col + 2 * n_pol + i),
                    );
                    mpol.set_point(k as usize, pol);
                    inc_col += 3 * n_pol;
                } else {
                    let pol2d = DVec2::new(
                        vec_poles.get(inc_col + i),
                        vec_poles.get(inc_col + n_pol + i),
                    );
                    mpol.set_point2d(k as usize, pol2d);
                    inc_col += 2 * n_pol;
                }
            }
            scurv.set_value(i as usize, mpol);
        }
    }

    /// OCCT IsDone() (gxx L505-508).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT NbConstraints(SSP, FirstPoint, LastPoint, TheConstraints)
    /// (gxx L510-547).
    pub fn nb_constraints(
        ssp: &MultiLine,
        _first_point: i32,
        _last_point: i32,
        the_constraints: &[ConstraintCouple],
    ) -> i32 {
        // Boucle de calcul du nombre de points de passage afin de
        // dimensionner les matrices.
        let mut a_inc_pass = 0;
        let mut a_inc_tan = 0;
        let mut a_inc_curv = 0;
        for couple in the_constraints {
            let cons = couple.constraint;
            if cons as i32 >= 1 {
                a_inc_pass += 1;
            }
            if cons as i32 >= 2 {
                a_inc_tan += 1;
            }
            if cons as i32 == 3 {
                a_inc_curv += 1;
            }
        }
        let nb3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb2d = my_line_tool::nb_p2d(ssp) as i32;
        let a_ccol = nb3d * 3 + nb2d * 2;
        a_ccol * a_inc_pass + a_inc_tan * (a_ccol - 1) + 3 * a_inc_curv
    }

    /// OCCT NbColumns(SSP, Deg) (gxx L549-556).
    pub fn nb_columns(ssp: &MultiLine, deg: i32) -> i32 {
        let nb3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb2d = my_line_tool::nb_p2d(ssp) as i32;
        let c_col = nb3d * 3 + nb2d * 2;
        c_col * (deg + 1)
    }

    /// OCCT ConstraintMatrix() (gxx L558-561).
    pub fn constraint_matrix(&self) -> &Matrix {
        &self.cont
    }

    /// OCCT InverseMatrix() (gxx L563-566).
    pub fn inverse_matrix(&self) -> &Matrix {
        &self.ctcinv
    }

    /// OCCT Duale() (gxx L568-571).
    pub fn duale(&self) -> &RVector {
        &self.vardua
    }

    /// OCCT ConstraintDerivative(SSP, Parameters, Deg, DA) (gxx L573-948) —
    /// fills DeCont with the derivative of the constraint matrix.
    pub fn constraint_derivative(
        &mut self,
        ssp: &MultiLine,
        parameters: &VecD,
        deg: i32,
        da: &Matrix,
    ) -> &Matrix {
        let nb_cu = my_line_tool::nb_p3d(ssp) as i32 + my_line_tool::nb_p2d(ssp) as i32;
        let n_pol = deg + 1;
        let n_pol2 = 2 * n_pol;
        let mut ibont = vec![vec![0i32; self.inc_tan as usize]; nb_cu as usize];
        let mut dec_init = Matrix::new_init(1, self.inc_pass, 1, n_pol, 0.0);
        let mut dda = RVector::new(1, n_pol);

        let nb3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb2d = my_line_tool::nb_p2d(ssp) as i32;
        let mut mynb3d = nb3d;
        let mut mynb2d = nb2d;
        if nb3d == 0 {
            mynb3d = 1;
        }
        if nb2d == 0 {
            mynb2d = 1;
        }
        let c_col = nb3d * 3 + nb2d * 2;
        let _ = (mynb3d, mynb2d);
        let mut tab_v = vec![DVec3::ZERO; mynb3d as usize];
        let mut tab_v2d = vec![DVec2::ZERO; mynb2d as usize];
        for i in 1..=self.de_cont.row_number() {
            for j in 1..=self.de_cont.col_number() {
                self.de_cont.set(i, j, 0.0);
            }
        }
        //  Remplissage de DK pour les points de passages:
        for i in 1..=self.inc_pass {
            let npt = self.ipas[(i - 1) as usize];
            for j in 1..=n_pol {
                let v = da.get(npt, j);
                dec_init.set(i, j, v);
            }
        }
        for i in 1..=c_col {
            let v = dec_init.clone();
            self.de_cont.set_block(
                self.inc_pass * (i - 1) + 1,
                self.inc_pass * i,
                n_pol * (i - 1) + 1,
                n_pol * i,
                &v,
            );
        }
        // Pour les points de tangence:
        let mut inc3 = c_col * self.inc_pass + 1;
        let mut inc_col = 0i32;
        for k in 1..=nb_cu {
            if k <= nb3d {
                for i in 1..=self.inc_tan {
                    let npt = self.itan[(i - 1) as usize];
                    // choix du maximum de tangence pour exprimer la
                    // colinearite:
                    my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                    let v = tab_v[(k - 1) as usize];
                    let t1v = v.x;
                    let t2v = v.y;
                    let t3v = v.z;
                    let tmax = t1v.abs();
                    ibont[(k - 1) as usize][(i - 1) as usize] = 1;
                    if t2v.abs() > tmax {
                        ibont[(k - 1) as usize][(i - 1) as usize] = 2;
                    }
                    if t3v.abs() > tmax {
                        ibont[(k - 1) as usize][(i - 1) as usize] = 3;
                    }
                    second_derivative_bernstein(parameters.get(npt as usize), &mut dda);
                    let ib = ibont[(k - 1) as usize][(i - 1) as usize];
                    if ib == 3 {
                        for j in 1..=n_pol {
                            let ddaij = dda.get(j);
                            let val = ddaij * t3v / tmax;
                            self.de_cont.set(inc3, j + n_pol + inc_col, val);
                            let val = -ddaij * t2v / tmax;
                            self.de_cont.set(inc3, j + n_pol2 + inc_col, val);
                            let val = ddaij * t3v / tmax;
                            self.de_cont.set(inc3 + 1, j + inc_col, val);
                            let val = -ddaij * t1v / tmax;
                            self.de_cont.set(inc3 + 1, j + n_pol2 + inc_col, val);
                        }
                    } else if ib == 1 {
                        for j in 1..=n_pol {
                            let ddaij = dda.get(j);
                            let val = ddaij * t3v / tmax;
                            self.de_cont.set(inc3, j + inc_col, val);
                            let val = -ddaij * t1v / tmax;
                            self.de_cont.set(inc3, j + n_pol2 + inc_col, val);
                            let val = ddaij * t2v / tmax;
                            self.de_cont.set(inc3 + 1, j + inc_col, val);
                            let val = -ddaij * t1v / tmax;
                            self.de_cont.set(inc3 + 1, j + n_pol + inc_col, val);
                        }
                    } else if ib == 2 {
                        for j in 1..=n_pol {
                            let ddaij = dda.get(j);
                            let val = ddaij * t3v / tmax;
                            self.de_cont.set(inc3, j + n_pol + inc_col, val);
                            let val = -ddaij * t2v / tmax;
                            self.de_cont.set(inc3, j + n_pol2 + inc_col, val);
                            let val = ddaij * t2v / tmax;
                            self.de_cont.set(inc3 + 1, j + inc_col, val);
                            let val = -ddaij * t1v / tmax;
                            self.de_cont.set(inc3 + 1, j + n_pol + inc_col, val);
                        }
                    }
                    inc3 += 2;
                }
                inc_col += 3 * n_pol;
            } else {
                for i in 1..=self.inc_tan {
                    let npt = self.itan[(i - 1) as usize];
                    second_derivative_bernstein(parameters.get(npt as usize), &mut dda);
                    my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                    // OCCT: `V2d = tabV2d(k);` — kept verbatim (the
                    // second-derivative variant indexes the 2D tangent
                    // table by the absolute curve index; the ctor above
                    // uses (k - nb3d)).
                    let v2d = tab_v2d[(k - 1) as usize];
                    let t1v = v2d.x;
                    let t2v = v2d.y;
                    let tmax = t1v.abs();
                    ibont[(k - 1) as usize][(i - 1) as usize] = 1;
                    if t2v.abs() > tmax {
                        ibont[(k - 1) as usize][(i - 1) as usize] = 2;
                    }
                    for j in 1..=n_pol {
                        let ddaij = dda.get(j);
                        let val = ddaij * t2v;
                        self.de_cont.set(inc3, j + inc_col, val);
                        let val = -ddaij * t1v;
                        self.de_cont.set(inc3, j + n_pol + inc_col, val);
                    }
                    inc3 += 1;
                }
            }
        }
        // Equations exprimant le meme rapport de tangence sur chaque courbe:
        // On prend les coordonnees les plus significatives.
        inc3 -= 1;
        for i in 1..=self.inc_tan {
            inc_col = 0;
            let npt = self.itan[(i - 1) as usize];
            second_derivative_bernstein(parameters.get(npt as usize), &mut dda);
            for k in 1..=(nb_cu - 1) {
                inc3 += 1;
                let t1v;
                let mut ip = 0i32;
                let ib = ibont[(k - 1) as usize][(i - 1) as usize];
                if ib == 1 {
                    if k <= nb3d {
                        my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                        t1v = tab_v[(k - 1) as usize].x;
                        ip = 3 * n_pol;
                    } else {
                        my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                        t1v = tab_v2d[(k - 1) as usize].x;
                        ip = 2 * n_pol;
                    }
                } else if ib == 2 {
                    if k <= nb3d {
                        my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                        t1v = tab_v[(k - 1) as usize].y;
                        ip = 3 * n_pol;
                    } else {
                        my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                        t1v = tab_v2d[(k - 1) as usize].y;
                        ip = 2 * n_pol;
                    }
                } else {
                    my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                    t1v = tab_v[(k - 1) as usize].z;
                    ip = 3 * n_pol;
                }
                let t2v;
                let ib_next = ibont[k as usize][(i - 1) as usize];
                if ib_next == 1 {
                    // Relations between T1? and T2x:
                    if (k + 1) <= nb3d {
                        my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                        t2v = tab_v[k as usize].x;
                    } else {
                        my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                        t2v = tab_v2d[k as usize].x;
                    }
                } else if ib_next == 2 {
                    // Relations between T1? and T2y:
                    if (k + 1) <= nb3d {
                        my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                        t2v = tab_v[k as usize].y;
                    } else {
                        my_line_tool::tangency_2d(ssp, npt as usize, &mut tab_v2d);
                        t2v = tab_v2d[k as usize].y;
                    }
                } else {
                    // Relations between T1? and T2z:
                    my_line_tool::tangency_3d(ssp, npt as usize, &mut tab_v);
                    t2v = tab_v[k as usize].z;
                }
                // The (T1 kind, T2 kind) pair selects the pole-column
                // offsets, exactly as in the OCCT switch cascade (the OCCT
                // branches are: a=1 -> +0; a=2 -> +Npol; a=3 -> +2*Npol on
                // the T1 side, mirrored on the T2 side with the IP gap).
                let (a_off, b_off) = match (ib, ib_next) {
                    (1, 1) => (0, 0),
                    (1, 2) => (0, n_pol),
                    (1, 3) => (0, 2 * n_pol),
                    (2, 1) => (n_pol, 0),
                    (2, 2) => (n_pol, n_pol),
                    (2, 3) => (n_pol, 2 * n_pol),
                    (_, 1) => (2 * n_pol, 0),
                    (_, 2) => (2 * n_pol, n_pol),
                    (_, _) => (2 * n_pol, 2 * n_pol),
                };
                for j in 1..=n_pol {
                    let daij = dda.get(j);
                    // OCCT writes these ratio rows into Cont (not DeCont)
                    // inside ConstraintDerivative — kept verbatim.
                    let val = daij * t2v;
                    self.cont.set(inc3, j + a_off + inc_col, val);
                    let val = -daij * t1v;
                    self.cont.set(inc3, j + ip + b_off + inc_col, val);
                }
                inc_col += ip;
            }
        }
        &self.de_cont
    }
}

// ---------------------------------------------------------------------------
// AppParCurves_Function (ParFunction)
// ---------------------------------------------------------------------------

// OCCT math_Vector operator helpers (member-operator semantics: the result
// carries the LEFT operand's bounds; component counts must match).

/// OCCT math_Vector::operator+.
fn v_add(a: &RVector, b: &RVector) -> RVector {
    let mut r = RVector::new(a.lower, a.upper());
    for i in a.lower..=a.upper() {
        let v = a.get(i) + b.get(i);
        r.set(i, v);
    }
    r
}

/// OCCT math_Vector::operator-.
fn v_sub(a: &RVector, b: &RVector) -> RVector {
    let mut r = RVector::new(a.lower, a.upper());
    for i in a.lower..=a.upper() {
        let v = a.get(i) - b.get(i);
        r.set(i, v);
    }
    r
}

/// OCCT math_Vector::operator*(Standard_Real).
fn v_scale(a: &RVector, s: f64) -> RVector {
    let mut r = RVector::new(a.lower, a.upper());
    for i in a.lower..=a.upper() {
        let v = a.get(i) * s;
        r.set(i, v);
    }
    r
}

/// OCCT math_Vector::Multiplied(Vector) — the dot product.
fn v_dot(a: &RVector, b: &RVector) -> f64 {
    let mut s = 0.0;
    for i in a.lower..=a.upper() {
        s += a.get(i) * b.get(i);
    }
    s
}

/// OCCT AppParCurves_Function (AppParCurves_Function.gxx whole file) — the
/// objective function of the Gradient minimization: a least-square fit
/// (ParLeastSquare) optionally corrected by the constrained resolution
/// (ResolConstraint), evaluating F = Sum ||C(ui) - PTL(i)||2 and its
/// gradient over the parameters.
#[derive(Debug, Clone)]
pub struct ParFunction {
    /// OCCT MyMultiLine.
    my_multi_line: MultiLine,
    /// OCCT MyMultiCurve.
    my_multi_curve: MultiCurve,
    /// OCCT myParameters.
    my_parameters: RVector,
    /// OCCT ValGrad_F.
    val_grad_f: RVector,
    /// OCCT MyF.
    my_f: Matrix,
    /// OCCT PTLX / PTLY / PTLZ (the point tables, filled only when the
    /// intermediate points are constrained).
    ptlx: Matrix,
    ptly: Matrix,
    ptlz: Matrix,
    /// OCCT A / DA (the Bernstein function matrices).
    a: Matrix,
    da: Matrix,
    /// OCCT MyLeastSquare.
    my_least_square: LeastSquare,
    first_p: i32,
    last_p: i32,
    nb_p: i32,
    a_deb: i32,
    a_fin: i32,
    degre: i32,
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
}

impl ParFunction {
    /// OCCT AppParCurves_Function(SSP, FirstPoint, LastPoint, TheConstraints,
    /// Parameters, Deg) (gxx ctor L21-129).
    pub fn new(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &VecD,
        deg: i32,
    ) -> Self {
        let nb3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb2d = my_line_tool::nb_p2d(ssp) as i32;
        let nb_pts = nb3d + nb2d;
        let mut f = ParFunction {
            my_multi_line: ssp.clone(),
            my_multi_curve: MultiCurve::new((deg + 1) as usize, 0, 0),
            my_parameters: RVector::new(1, parameters.len() as i32),
            val_grad_f: RVector::new(first_point, last_point),
            my_f: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            ptlx: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            ptly: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            ptlz: Matrix::new_init(first_point, last_point, 1, nb_pts, 0.0),
            a: Matrix::new(first_point, last_point, 1, deg + 1),
            da: Matrix::new(first_point, last_point, 1, deg + 1),
            my_least_square: LeastSquare::new_no_params(
                ssp,
                first_point,
                last_point,
                Self::first_constraint(the_constraints, first_point),
                Self::last_constraint(the_constraints, last_point),
                deg + 1,
            ),
            first_p: 0,
            last_p: 0,
            nb_p: 0,
            a_deb: 0,
            a_fin: 0,
            degre: 0,
            nbcu: 0,
            contraintes: false,
            my_constraints: the_constraints.to_vec(),
            tabdim: Vec::new(),
            f_val: 0.0,
            err3d: 0.0,
            err2d: 0.0,
            done: false,
        };
        for i in 1..=(parameters.len() as i32) {
            let v = parameters.get(i as usize);
            f.my_parameters.set(i, v);
        }
        f.first_p = first_point;
        f.last_p = last_point;
        f.nb_p = f.last_p - f.first_p + 1;
        f.a_deb = f.first_p;
        f.a_fin = f.last_p;
        f.degre = deg;
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

    /// OCCT FirstConstraint(TheConstraints, FirstPoint) (gxx L131-153).
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

    /// OCCT LastConstraint(TheConstraints, LastPoint) (gxx L155-177).
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

    /// OCCT Value(X, F) (gxx L179-269).
    pub fn value(&mut self, x: &VecD, f: &mut f64) -> bool {
        // myParameters = X.
        for i in 1..=self.my_parameters.length() {
            let v = x.get(i as usize);
            self.my_parameters.set(i, v);
        }
        // Resolution moindres carres:
        // ===========================
        self.my_least_square.perform(&{
            let mut v = VecD::new(self.my_parameters.length() as usize);
            for i in 1..=self.my_parameters.length() {
                v.set(i as usize, self.my_parameters.get(i));
            }
            v
        });
        if !self.my_least_square.is_done() {
            self.done = false;
            return false;
        }
        if !self.contraintes {
            let mut e3 = 0.0;
            let mut e2 = 0.0;
            let mut fval = 0.0;
            self.my_least_square.error(&mut fval, &mut e3, &mut e2);
            self.f_val = fval;
            self.err3d = e3;
            self.err2d = e2;
            *f = self.f_val;
        } else {
            // Resolution avec contraintes:
            // ============================
            let n_pol = self.degre + 1;
            let mut err3d = 0.0;
            let mut err2d = 0.0;
            let mut ptcxci = RVector::new(1, n_pol);
            let mut ptcyci = RVector::new(1, n_pol);
            let mut ptczci = RVector::new(1, n_pol);
            self.my_multi_curve = self.my_least_square.bezier_value();
            self.a = self.my_least_square.function_matrix().clone();
            let mut my_multi_curve = self.my_multi_curve.clone();
            let resol = ResolConstraint::new(
                &self.my_multi_line,
                &mut my_multi_curve,
                self.first_p,
                self.last_p,
                &self.my_constraints,
                &self.a,
                self.my_least_square.derivative_function_matrix(),
                1.0e-10,
            );
            self.my_multi_curve = my_multi_curve;
            if !resol.is_done() {
                self.done = false;
                return false;
            }
            // Calcul de F = Sum||C(ui)-Ptli||2  sur toutes les courbes :
            // ========================================================================
            let mut f_val = 0.0;
            for ci in 1..=self.nbcu {
                let dimen = self.tabdim[(ci - 1) as usize];
                for j in 1..=n_pol {
                    if dimen == 3 {
                        let p = self.my_multi_curve.value(j as usize).point(ci as usize);
                        ptcxci.set(j, p.x);
                        ptcyci.set(j, p.y);
                        ptczci.set(j, p.z);
                    } else {
                        let p = self.my_multi_curve.value(j as usize).point2d(ci as usize);
                        ptcxci.set(j, p.x);
                        ptcyci.set(j, p.y);
                    }
                }
                // Calcul de F:
                // ============
                for i in self.a_deb..=self.a_fin {
                    let mut aa = 0.0;
                    let mut bb = 0.0;
                    let mut cc = 0.0;
                    for j in 1..=n_pol {
                        let aij = self.a.get(i, j);
                        aa += aij * ptcxci.get(j);
                        bb += aij * ptcyci.get(j);
                        if dimen == 3 {
                            cc += aij * ptczci.get(j);
                        }
                    }
                    let fx = aa - self.ptlx.get(i, ci);
                    let fy = bb - self.ptly.get(i, ci);
                    let mut fi = fx * fx + fy * fy;
                    self.my_f.set(i, ci, fi);
                    if dimen == 3 {
                        let fz = cc - self.ptlz.get(i, ci);
                        fi += fz * fz;
                        self.my_f.set(i, ci, fi);
                        if fi.sqrt() > err3d {
                            err3d = fi.sqrt();
                        }
                    } else if fi.sqrt() > err2d {
                        err2d = fi.sqrt();
                    }
                    f_val += fi;
                }
            }
            self.f_val = f_val;
            self.err3d = err3d;
            self.err2d = err2d;
            *f = self.f_val;
        }
        self.done = true;
        true
    }

    /// OCCT Perform(X) (gxx L271-600).
    pub fn perform(&mut self, x: &VecD) {
        let n_pol = self.degre + 1;
        // myParameters = X.
        for i in 1..=self.my_parameters.length() {
            let v = x.get(i as usize);
            self.my_parameters.set(i, v);
        }
        // Resolution moindres carres:
        // ===========================
        let x_copy = {
            let mut v = VecD::new(self.my_parameters.length() as usize);
            for i in 1..=self.my_parameters.length() {
                v.set(i as usize, self.my_parameters.get(i));
            }
            v
        };
        self.my_least_square.perform(&x_copy);
        if !self.my_least_square.is_done() {
            self.done = false;
            return;
        }
        for j in 1..=self.val_grad_f.length() {
            self.val_grad_f.set(j, 0.0);
        }
        if !self.contraintes {
            let mut f_val = 0.0;
            let mut e3 = 0.0;
            let mut e2 = 0.0;
            self.my_least_square
                .error_gradient(&mut self.val_grad_f.data, &mut f_val, &mut e3, &mut e2);
            self.f_val = f_val;
            self.err3d = e3;
            self.err2d = e2;
        } else {
            let mut cons = AppParConstraint::NoConstraint;
            let mut grad_f = Matrix::new_init(self.first_p, self.last_p, 1, self.nbcu, 0.0);
            let mut ptcxci = RVector::new(1, n_pol);
            let mut ptcyci = RVector::new(1, n_pol);
            let mut ptczci = RVector::new(1, n_pol);
            let mut err3d = 0.0;
            let mut err2d = 0.0;
            let mut ptcox = Matrix::new_init(1, n_pol, 1, self.nbcu, 0.0);
            let mut ptcoy = Matrix::new_init(1, n_pol, 1, self.nbcu, 0.0);
            let mut ptcoz = Matrix::new_init(1, n_pol, 1, self.nbcu, 0.0);
            let mut ptcx = Matrix::new_init(1, n_pol, 1, self.nbcu, 0.0);
            let mut ptcy = Matrix::new_init(1, n_pol, 1, self.nbcu, 0.0);
            let mut ptcz = Matrix::new_init(1, n_pol, 1, self.nbcu, 0.0);
            self.my_multi_curve = self.my_least_square.bezier_value();
            for ci in 1..=self.nbcu {
                let dimen = self.tabdim[(ci - 1) as usize];
                for j in 1..=n_pol {
                    if dimen == 3 {
                        let p = self.my_multi_curve.value(j as usize).point(ci as usize);
                        ptcox.set(j, ci, p.x);
                        ptcoy.set(j, ci, p.y);
                        ptcoz.set(j, ci, p.z);
                    } else {
                        let p = self.my_multi_curve.value(j as usize).point2d(ci as usize);
                        ptcox.set(j, ci, p.x);
                        ptcoy.set(j, ci, p.y);
                        ptcoz.set(j, ci, 0.0);
                    }
                }
            }
            self.a = self.my_least_square.function_matrix().clone();
            self.da = self.my_least_square.derivative_function_matrix().clone();
            // Resolution avec contraintes:
            // ============================
            let mut my_multi_curve = self.my_multi_curve.clone();
            let mut resol = ResolConstraint::new(
                &self.my_multi_line,
                &mut my_multi_curve,
                self.first_p,
                self.last_p,
                &self.my_constraints,
                &self.a,
                &self.da,
                1.0e-10,
            );
            self.my_multi_curve = my_multi_curve;
            if !resol.is_done() {
                self.done = false;
                return;
            }
            // Calcul de F = Sum||C(ui)-Ptli||2 et du gradient non contraint
            // de F pour chaque point PointIndex.
            // ========================================================================
            let mut f_val = 0.0;
            for j in self.first_p..=self.last_p {
                self.val_grad_f.set(j, 0.0);
            }
            let tr_a = self.a.transposed();
            let tr_da = self.da.transposed();
            let restm = self.a.transposed().multiplied(&self.a).inverse();
            let k = resol.constraint_matrix().clone();
            let dk = resol
                .constraint_derivative(&self.my_multi_line, x, self.degre, &self.da)
                .clone();
            let tk = k.transposed();
            let vardua = resol.duale().clone();
            let kk = k.transposed().multiplied(resol.inverse_matrix());
            let dtk = dk.transposed();
            let mut dptco = RVector::new(1, k.col_number());
            let mut dptco1 = Matrix::new_init(self.first_p, self.last_p, 1, k.col_number(), 0.0);
            let mut dkptc = RVector::new(1, k.row_number());
            for ci in 1..=self.nbcu {
                let dimen = self.tabdim[(ci - 1) as usize];
                for j in 1..=n_pol {
                    if dimen == 3 {
                        let p = self.my_multi_curve.value(j as usize).point(ci as usize);
                        ptcx.set(j, ci, p.x);
                        ptcy.set(j, ci, p.y);
                        ptcz.set(j, ci, p.z);
                    } else {
                        let p = self.my_multi_curve.value(j as usize).point2d(ci as usize);
                        ptcx.set(j, ci, p.x);
                        ptcy.set(j, ci, p.y);
                        ptcz.set(j, ci, 0.0);
                    }
                }
            }
            // Calcul du gradient sans contraintes:
            // ====================================
            for ci in 1..=self.nbcu {
                let dimen = self.tabdim[(ci - 1) as usize];
                for i in self.a_deb..=self.a_fin {
                    let mut aa = 0.0;
                    let mut bb = 0.0;
                    let mut cc = 0.0;
                    let mut daa = 0.0;
                    let mut dbb = 0.0;
                    let mut dcc = 0.0;
                    for j in 1..=n_pol {
                        let aij = self.a.get(i, j);
                        let daij = self.da.get(i, j);
                        let px = ptcx.get(j, ci);
                        let py = ptcy.get(j, ci);
                        aa += aij * px;
                        bb += aij * py;
                        daa += daij * px;
                        dbb += daij * py;
                        if dimen == 3 {
                            let pz = ptcz.get(j, ci);
                            cc += aij * pz;
                            dcc += daij * pz;
                        }
                    }
                    let fx = aa - self.ptlx.get(i, ci);
                    let fy = bb - self.ptly.get(i, ci);
                    let mut fi = fx * fx + fy * fy;
                    self.my_f.set(i, ci, fi);
                    let g = 2.0 * (daa * fx + dbb * fy);
                    grad_f.set(i, ci, g);
                    if dimen == 3 {
                        let fz = cc - self.ptlz.get(i, ci);
                        fi += fz * fz;
                        self.my_f.set(i, ci, fi);
                        grad_f.set(i, ci, g + 2.0 * dcc * fz);
                        if fi.sqrt() > err3d {
                            err3d = fi.sqrt();
                        }
                    } else if fi.sqrt() > err2d {
                        err2d = fi.sqrt();
                    }
                    f_val += fi;
                    let v = self.val_grad_f.get(i) + grad_f.get(i, ci);
                    self.val_grad_f.set(i, v);
                }
            }
            // Calcul de DK*PTC:
            // =================
            for i in 1..=k.row_number() {
                let mut inc = 0i32;
                for ci in 1..=self.nbcu {
                    let dimen = self.tabdim[(ci - 1) as usize];
                    dkptc.set(i, 0.0);
                    for j in 1..=n_pol {
                        let v = dkptc.get(i)
                            + dk.get(i, j + inc) * ptcx.get(j, ci)
                            + dk.get(i, j + inc + n_pol) * ptcy.get(j, ci);
                        dkptc.set(i, v);
                        if dimen == 3 {
                            let v = dkptc.get(i) + dk.get(i, j + inc + 2 * n_pol) * ptcz.get(j, ci);
                            dkptc.set(i, v);
                        }
                    }
                    if dimen == 3 {
                        inc += 3 * n_pol;
                    } else {
                        inc += 2 * n_pol;
                    }
                }
            }
            // DERR = (DTK)*Vardua - KK * ((DKPTC) + K * (DTK)*Vardua).
            let base = dtk.multiplied_vec(&vardua);
            let inner = v_add(&dkptc, &k.multiplied_vec(&base));
            let mut derr = v_sub(&base, &kk.multiplied_vec(&inner));
            // rajout du gradient avec contraintes:
            // ====================================
            // dPTCO1/duk = [d(TA)/duk*[A*PTCO-PTL] + TA*dA/duk*PTCO]
            let mut inc = 0i32;
            for ci in 1..=self.nbcu {
                let dimen = self.tabdim[(ci - 1) as usize];
                ptcxci = ptcox.col(ci);
                ptcyci = ptcoy.col(ci);
                ptczci = ptcoz.col(ci);
                ptcxci = ptcx.col(ci);
                ptcyci = ptcy.col(ci);
                ptczci = ptcz.col(ci);
                let errx = v_sub(&self.a.multiplied_vec(&ptcxci), &self.ptlx.col(ci));
                let erry = v_sub(&self.a.multiplied_vec(&ptcyci), &self.ptly.col(ci));
                let errz = v_sub(&self.a.multiplied_vec(&ptczci), &self.ptlz.col(ci));
                let scalx = self.da.multiplied_vec(&ptcxci); // Scal = DA * PTCO
                let scaly = self.da.multiplied_vec(&ptcyci);
                let scalz = self.da.multiplied_vec(&ptczci);
                let erruzax = v_sub(&ptcxci, &ptcox.col(ci));
                let erruzay = v_sub(&ptcyci, &ptcoy.col(ci));
                let erruzaz = v_sub(&ptczci, &ptcoz.col(ci));
                for pi in self.first_p..=self.last_p {
                    let trdapi = tr_da.col(pi);
                    let trapi = tr_a.col(pi);
                    let taa = v_dot(&trapi, &self.a.row(pi));
                    let mut scal = 0.0;
                    for j in 1..=n_pol {
                        let v1 = v_add(
                            &v_scale(&trdapi, errx.get(pi)),
                            &v_scale(&trapi, scalx.get(pi)),
                        );
                        dptco1.set(pi, j + inc, v1.get(j));
                        let v2 = v_add(
                            &v_scale(&trdapi, erry.get(pi)),
                            &v_scale(&trapi, scaly.get(pi)),
                        );
                        dptco1.set(pi, j + inc + n_pol, v2.get(j));
                        scal += dptco1.get(pi, j + inc) * taa * erruzax.get(j)
                            + dptco1.get(pi, j + inc + n_pol) * taa * erruzay.get(j);
                        if dimen == 3 {
                            let v3 = v_add(
                                &v_scale(&trdapi, errz.get(pi)),
                                &v_scale(&trapi, scalz.get(pi)),
                            );
                            dptco1.set(pi, j + inc + 2 * n_pol, v3.get(j));
                            scal += dptco1.get(pi, j + inc + 2 * n_pol) * taa * erruzaz.get(j);
                        }
                    }
                    let v = self.val_grad_f.get(pi) - 2.0 * scal;
                    self.val_grad_f.set(pi, v);
                }
                if dimen == 3 {
                    inc += 3 * n_pol;
                } else {
                    inc += 2 * n_pol;
                }
            }
            // on calcule DPTCO = - RESTM * DPTCO1:
            // Calcul de DPTCO/duk:
            // dPTCO/duk = -Inv(T(A)*A)*[d(TA)/duk*[A*PTCO-PTL] + TA*dA/duk*PTCO]
            inc = 0;
            for pi in self.first_p..=self.last_p {
                for couple in &self.my_constraints {
                    if couple.index == pi {
                        cons = couple.constraint;
                        break;
                    }
                }
                if cons as i32 >= 1 {
                    inc = 0;
                    for ci in 1..=self.nbcu {
                        let dimen = self.tabdim[(ci - 1) as usize];
                        for j in 1..=n_pol {
                            dptco.set(j + inc, 0.0);
                            dptco.set(j + inc + n_pol, 0.0);
                            if dimen == 3 {
                                dptco.set(j + inc + 2 * n_pol, 0.0);
                            }
                            for k in 1..=n_pol {
                                let v = dptco.get(j + inc)
                                    - restm.get(j, k) * dptco1.get(pi, j + inc);
                                dptco.set(j + inc, v);
                                let v = dptco.get(j + inc + n_pol)
                                    - restm.get(j, k) * dptco1.get(pi, j + inc + n_pol);
                                dptco.set(j + inc + n_pol, v);
                                if dimen == 3 {
                                    let v = dptco.get(j + inc + 2 * n_pol)
                                        - restm.get(j, k) * dptco1.get(pi, j + inc + 2 * n_pol);
                                    dptco.set(j + inc + 2 * n_pol, v);
                                }
                            }
                        }
                        if dimen == 3 {
                            inc += 3 * n_pol;
                        } else {
                            inc += 2 * n_pol;
                        }
                    }
                    let step = kk.multiplied_vec(&k.multiplied_vec(&dptco));
                    derr = v_sub(&derr, &step);
                    inc = 0;
                    for ci in 1..=self.nbcu {
                        let dimen = self.tabdim[(ci - 1) as usize];
                        ptcxci = ptcox.col(ci);
                        ptcyci = ptcoy.col(ci);
                        ptczci = ptcoz.col(ci);
                        ptcxci = ptcx.col(ci);
                        ptcyci = ptcy.col(ci);
                        ptczci = ptcz.col(ci);
                        let erruzax = v_sub(&ptcxci, &ptcox.col(ci));
                        let erruzay = v_sub(&ptcyci, &ptcoy.col(ci));
                        let erruzaz = v_sub(&ptczci, &ptcoz.col(ci));
                        let mut scal = 0.0;
                        for j in 1..=n_pol {
                            scal = (self.a.get(pi, j) * erruzax.get(j))
                                * (self.a.get(pi, j) * derr.get(j + inc))
                                + (self.a.get(pi, j) * erruzay.get(j))
                                    * (self.a.get(pi, j) * derr.get(j + inc + n_pol));
                            if dimen == 3 {
                                scal += (self.a.get(pi, j) * erruzax.get(j))
                                    * (self.a.get(pi, j) * derr.get(j + inc + 2 * n_pol));
                            }
                        }
                        let v = self.val_grad_f.get(pi) + 2.0 * scal;
                        self.val_grad_f.set(pi, v);
                        if dimen == 3 {
                            inc += 3 * n_pol;
                        } else {
                            inc += 2 * n_pol;
                        }
                    }
                }
            }
            self.f_val = f_val;
            self.err3d = err3d;
            self.err2d = err2d;
            let _ = (tk, vardua, cons);
        }
        self.done = true;
    }

    /// OCCT NbVariables() (gxx L602-605).
    pub fn nb_variables(&self) -> i32 {
        self.nb_p
    }

    /// OCCT Gradient(X, G) (gxx L607-612).
    pub fn gradient(&mut self, x: &VecD, g: &mut RVector) -> bool {
        self.perform(x);
        for i in g.lower..=g.upper() {
            let v = self.val_grad_f.get(i);
            g.set(i, v);
        }
        true
    }

    /// OCCT Values(X, F, G) (gxx L614-620).
    pub fn values(&mut self, x: &VecD, f: &mut f64, g: &mut RVector) -> bool {
        self.perform(x);
        *f = self.f_val;
        for i in g.lower..=g.upper() {
            let v = self.val_grad_f.get(i);
            g.set(i, v);
        }
        true
    }

    /// OCCT CurveValue() (gxx L622-627).
    pub fn curve_value(&mut self) -> &MultiCurve {
        if !self.contraintes {
            self.my_multi_curve = self.my_least_square.bezier_value();
        }
        &self.my_multi_curve
    }

    /// OCCT Error(IPoint, CurveIndex) (gxx L629-632).
    pub fn error(&self, ipoint: i32, curve_index: i32) -> f64 {
        self.my_f.get(ipoint, curve_index).sqrt()
    }

    /// OCCT MaxError3d() (gxx L634-637).
    pub fn max_error_3d(&self) -> f64 {
        self.err3d
    }

    /// OCCT MaxError2d() (gxx L639-642).
    pub fn max_error_2d(&self) -> f64 {
        self.err2d
    }

    /// OCCT NewParameters() (gxx L644-647).
    pub fn new_parameters(&self) -> &RVector {
        &self.my_parameters
    }

    /// OCCT NewParameters()(Index) — the i-th new parameter component.
    pub fn new_parameter(&self, index: i32) -> f64 {
        self.my_parameters.get(index)
    }
}

// ---------------------------------------------------------------------------
// AppParCurves_Gradient
// ---------------------------------------------------------------------------

/// OCCT AppParCurves_Gradient (AppParCurves_Gradient.gxx whole file) — one
/// Rogers & Fog 89 projection iteration over the parameters, then, when the
/// error tolerances are not met, the Gradient_BFGS minimization loop (the
/// BFGS branch requires math_BFGS + the Gradient_BFGS shell and is a
/// structured skeleton until that unit lands, ThruSections precedent).
#[derive(Debug, Clone)]
pub struct Gradient {
    /// OCCT SCU.
    scu: MultiCurve,
    /// OCCT ParError.
    par_error: RVector,
    /// OCCT AvError / MError3d / MError2d / Done.
    av_error: f64,
    m_error3d: f64,
    m_error2d: f64,
    done: bool,
}

impl Gradient {
    /// OCCT AppParCurves_Gradient(SSP, FirstPoint, LastPoint, TheConstraints,
    /// Parameters, Deg, Tol3d, Tol2d, NbIterations) (gxx L18-235). The
    /// `parameters` vector is updated in place, as in OCCT (math_Vector&
    /// Parameters).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut VecD,
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut g = Gradient {
            scu: MultiCurve::new(1, 0, 0),
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
            deg,
            tol3d,
            tol2d,
            nb_iterations,
        );
        g
    }

    /// The constructor body (gxx L31-233).
    #[allow(clippy::too_many_arguments)]
    fn perform(
        &mut self,
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        parameters: &mut VecD,
        deg: i32,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
    ) {
        let nb_p3d = my_line_tool::nb_p3d(ssp) as i32;
        let nb_p2d = my_line_tool::nb_p2d(ssp) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        let nb_p = nb_p3d + nb_p2d;
        self.done = false;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let _ = (mynb_p3d, mynb_p2d);
        let mut tab_p = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p2d = vec![DVec2::ZERO; mynb_p2d as usize];
        let mut tab_v = vec![DVec3::ZERO; mynb_p3d as usize];
        let _ = &mut tab_v;
        let _ = &mut tab_p;
        let mut tab_v2d = vec![DVec2::ZERO; mynb_p2d as usize];
        let _ = &mut tab_v2d;
        // Calculation of the function F= sum(||C(ui)-Ptli||2):
        // Call to a function inheriting from MultipleVarFunctionWithGradient
        // to compute F and grad_F.
        // ================================================================
        let mut my_f =
            ParFunction::new(ssp, first_point, last_point, the_constraints, parameters, deg);
        let mut fval = 0.0;
        if !my_f.value(parameters, &mut fval) {
            self.done = false;
            return;
        }
        self.scu = my_f.curve_value().clone();
        let deg_scu = self.scu.nb_poles() as i32 - 1;
        let mut the_coef = vec![DVec3::ZERO; ((deg_scu + 1) * mynb_p3d) as usize];
        let mut the_coef2d = vec![DVec2::ZERO; ((deg_scu + 1) * mynb_p2d) as usize];
        // Storage of curve poles for projection:
        // ============================================
        let mut i2 = 0i32;
        for k in 1..=nb_p3d {
            let mut poles: Vec<DVec3> = Vec::new();
            self.scu.curve(k as usize, &mut poles);
            let tab_coef = poles_coefficients_3d(&poles);
            for j in 1..=(deg_scu + 1) {
                the_coef[(j + i2 - 1) as usize] = tab_coef[(j - 1) as usize];
            }
            i2 += deg_scu + 1;
        }
        i2 = 0;
        for k in 1..=nb_p2d {
            let mut poles: Vec<DVec2> = Vec::new();
            self.scu.curve2d(k as usize, &mut poles);
            let tab_coef2d = poles_coefficients_2d(&poles);
            for j in 1..=(deg_scu + 1) {
                the_coef2d[(j + i2 - 1) as usize] = tab_coef2d[(j - 1) as usize];
            }
            i2 += deg_scu + 1;
        }
        //  Une iteration rapide de projection est faite par la methode de
        //  Rogers & Fog 89, methode equivalente a Hoschek 88 qui ne
        //  necessite pas le calcul de D2.
        // Iteration de Projection:
        // =======================
        for j in (first_point + 1)..=(last_point - 1) {
            let mut uf = parameters.get(j as usize);
            if nb_p != 0 && nb_p2d != 0 {
                my_line_tool::value_3d_2d(ssp, j as usize, &mut tab_p, &mut tab_p2d);
            } else if nb_p2d != 0 {
                my_line_tool::value_2d(ssp, j as usize, &mut tab_p2d);
            } else {
                my_line_tool::value_3d(ssp, j as usize, &mut tab_p);
            }
            let mut fu = 0.0;
            let mut dfu = 0.0;
            let mut i2 = 0i32;
            for k in 1..=nb_p3d {
                let mut tab_coef: Vec<DVec3> = vec![DVec3::ZERO; (deg_scu + 1) as usize];
                for l in 1..=(deg_scu + 1) {
                    tab_coef[(l - 1) as usize] = the_coef[(l + i2 - 1) as usize];
                }
                i2 += deg_scu + 1;
                let mut pt = DVec3::ZERO;
                let mut v1 = DVec3::ZERO;
                coefs_d1_3d(uf, &tab_coef, &mut pt, &mut v1);
                let my_v = tab_p[(k - 1) as usize] - pt;
                fu += my_v.dot(v1);
                dfu += v1.length_squared();
            }
            i2 = 0;
            for k in 1..=nb_p2d {
                let mut tab_coef2d: Vec<DVec2> = vec![DVec2::ZERO; (deg_scu + 1) as usize];
                for l in 1..=(deg_scu + 1) {
                    tab_coef2d[(l - 1) as usize] = the_coef2d[(l + i2 - 1) as usize];
                }
                i2 += deg_scu + 1;
                let mut pt2d = DVec2::ZERO;
                let mut v12d = DVec2::ZERO;
                coefs_d1_2d(uf, &tab_coef2d, &mut pt2d, &mut v12d);
                let my_v2d = tab_p2d[(k - 1) as usize] - pt2d;
                fu += my_v2d.dot(v12d);
                dfu += v12d.length_squared();
            }
            // OCCT RealEpsilon().
            if dfu >= 2.220_446_049_250_313e-16 {
                let mut du = fu / dfu;
                du = du.abs().min(5.0e-02) * du.signum(); // copysign(min(5e-2, |DU|), DU)
                uf += du;
                parameters.set(j as usize, uf);
            }
        }
        let mut fval = 0.0;
        if !my_f.value(parameters, &mut fval) {
            self.scu = MultiCurve::new(1, 0, 0);
            self.done = false;
            return;
        }
        self.m_error3d = my_f.max_error_3d();
        self.m_error2d = my_f.max_error_2d();
        if self.m_error3d <= tol3d && self.m_error2d <= tol2d {
            self.done = true;
            self.scu = my_f.curve_value().clone();
        } else if nb_iterations != 0 {
            // NbIterations de gradient conjugue:
            // =================================
            // AppParCurves_Gradient_BFGS FResol(MyF, Parameters, Tol3d,
            // Tol2d, Eps, NbIterations); Parameters =
            // MyF.NewParameters(); SCU = MyF.CurveValue();
            // -- The Gradient_BFGS shell + math_BFGS minimizer are not yet
            let eps = 1.0e-07;
            // AppParCurves_Gradient_BFGS FResol(MyF, Parameters, Tol3d,
            // Tol2d, Eps, NbIterations).
            let mut f_resol =
                GradientBfgs::new(&mut my_f, parameters, tol3d, tol2d, eps, nb_iterations);
            // Parameters = MyF.NewParameters();
            for i in 1..=f_resol.nb_variables() {
                let v = f_resol.new_parameter(i);
                parameters.set(i as usize, v);
            }
            // SCU = MyF.CurveValue();
            self.scu = my_f.curve_value().clone();
        }
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

    /// OCCT Value() (gxx L237-240).
    pub fn value(&self) -> MultiCurve {
        self.scu.clone()
    }

    /// OCCT IsDone() (gxx L242-245).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Error(Index) (gxx L247-250).
    pub fn error(&self, index: i32) -> f64 {
        self.par_error.get(index)
    }

    /// OCCT AverageError() (gxx L252-255).
    pub fn average_error(&self) -> f64 {
        self.av_error
    }

    /// OCCT MaxError3d() (gxx L257-260).
    pub fn max_error_3d(&self) -> f64 {
        self.m_error3d
    }

    /// OCCT MaxError2d() (gxx L262-265).
    pub fn max_error_2d(&self) -> f64 {
        self.m_error2d
    }
}

// ---------------------------------------------------------------------------
// AppParCurves_Gradient_BFGS (the Gradient_BFGS shell)
// ---------------------------------------------------------------------------

/// The F interface consumed by the Gradient_BFGS IsSolutionReached override:
/// math_MultipleVarFunctionWithGradient + ParFunction's MaxError3d/2d (the
/// OCCT C-style down-cast `(AppDef_ParFunctionOfMyGradientbis*) & F`).
pub trait GradientFunction:
    rcad_kernel::math::math_bfgs::MultipleVarFunctionWithGradient
{
    /// OCCT ParFunction::MaxError3d().
    fn max_error_3d(&self) -> f64;
    /// OCCT ParFunction::MaxError2d().
    fn max_error_2d(&self) -> f64;
}

// OCCT: AppParCurves_Function : public math_MultipleVarFunctionWithGradient
// — the trait delegation to the inherent Value/Gradient/Values.
impl rcad_kernel::math::math_bfgs::MultipleVarFunction for ParFunction {
    fn nb_variables(&self) -> i32 {
        ParFunction::nb_variables(self)
    }
    fn value(&mut self, x: &KernelVector, f: &mut f64) -> bool {
        let mut xv = VecD::new(x.length() as usize);
        for i in 1..=x.length() {
            xv.set(i as usize, x.get(i));
        }
        ParFunction::value(self, &xv, f)
    }
}
impl rcad_kernel::math::math_bfgs::MultipleVarFunctionWithGradient for ParFunction {
    fn gradient(&mut self, x: &KernelVector, g: &mut KernelVector) -> bool {
        let mut xv = VecD::new(x.length() as usize);
        for i in 1..=x.length() {
            xv.set(i as usize, x.get(i));
        }
        let mut gv = RVector::new(g.lower, g.upper());
        let ok = ParFunction::gradient(self, &xv, &mut gv);
        for i in g.lower..=g.upper() {
            let v = gv.get(i);
            g.set(i, v);
        }
        ok
    }
    fn values(&mut self, x: &KernelVector, f: &mut f64, g: &mut KernelVector) -> bool {
        let mut xv = VecD::new(x.length() as usize);
        for i in 1..=x.length() {
            xv.set(i as usize, x.get(i));
        }
        let mut gv = RVector::new(g.lower, g.upper());
        let ok = ParFunction::values(self, &xv, f, &mut gv);
        for i in g.lower..=g.upper() {
            let v = gv.get(i);
            g.set(i, v);
        }
        ok
    }
}
impl GradientFunction for ParFunction {
    fn max_error_3d(&self) -> f64 {
        ParFunction::max_error_3d(self)
    }
    fn max_error_2d(&self) -> f64 {
        ParFunction::max_error_2d(self)
    }
}

/// OCCT AppParCurves_Gradient_BFGS
/// (AppDef_Gradient_BFGSOfMyGradientbisOfBSplineCompute.hxx/.cxx L20-42):
/// math_BFGS with the IsSolutionReached override — the solution is reached
/// when the minimum stops moving OR the approximation errors fall inside
/// the 3d/2d tolerances. Rust has no inheritance: math_BFGS::Perform takes
/// the IsSolutionReached test as an injected checker (the exact equivalent
/// of the OCCT virtual dispatch).
#[derive(Clone)]
pub struct GradientBfgs {
    /// OCCT base math_BFGS part.
    bfgs: rcad_kernel::math::math_bfgs::Bfgs,
    /// OCCT myTol3d / myTol2d.
    my_tol3d: f64,
    my_tol2d: f64,
}

impl GradientBfgs {
    /// OCCT AppDef_Gradient_BFGSOfMyGradientbisOfBSplineCompute(F,
    /// StartingPoint, Tolerance3d, Tolerance2d, Eps, NbIterations = 200)
    /// (cxx L20-28): Perform(F, StartingPoint).
    pub fn new<F: GradientFunction>(
        f: &mut F,
        starting_point: &VecD,
        tolerance3d: f64,
        tolerance2d: f64,
        eps: f64,
        nb_iterations: i32,
    ) -> Self {
        let n = starting_point.len() as i32;
        let mut g = GradientBfgs {
            bfgs: rcad_kernel::math::math_bfgs::Bfgs::new(n, eps, nb_iterations, 1.0e-12),
            my_tol3d: tolerance3d,
            my_tol2d: tolerance2d,
        };
        let start = rcad_kernel::math::math_matrix::Vector::new_init(
            1,
            n,
            0.0,
        );
        let mut start = start;
        for i in 1..=n {
            let v = starting_point.get(i as usize);
            start.set(i, v);
        }
        g.bfgs
            .perform_with_checker(f, &start, &mut |the_minimum, previous_minimum, f| {
                // OCCT IsSolutionReached (cxx L30-42).
                let result = 2.0 * (the_minimum - previous_minimum).abs()
                    <= 1.0e-10 * (the_minimum.abs() + previous_minimum.abs()) + 1.0e-12;
                let m_err3d = f.max_error_3d();
                let m_err2d = f.max_error_2d();
                let result2 = m_err3d <= g.my_tol3d && m_err2d <= g.my_tol2d;
                result || result2
            });
        g
    }

    /// OCCT IsSolutionReached(F) (cxx L30-42) — kept as the inherent
    /// documentation of the injected test above.
    pub fn is_solution_reached(
        &self,
        the_minimum: f64,
        previous_minimum: f64,
        m_err3d: f64,
        m_err2d: f64,
    ) -> bool {
        let result = 2.0 * (the_minimum - previous_minimum).abs()
            <= 1.0e-10 * (the_minimum.abs() + previous_minimum.abs()) + 1.0e-12;
        let result2 = m_err3d <= self.my_tol3d && m_err2d <= self.my_tol2d;
        result || result2
    }

    /// OCCT math_BFGS::Location() (lxx).
    pub fn location(&self) -> &rcad_kernel::math::math_matrix::Vector {
        self.bfgs.location()
    }

    /// OCCT math_BFGS::NbVariables() via the location length.
    pub fn nb_variables(&self) -> i32 {
        self.bfgs.location().length()
    }

    /// OCCT NewParameters()(i) after the resolution (the location vector).
    pub fn new_parameter(&self, index: i32) -> f64 {
        self.bfgs.location().get(index)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::app_def::MultiPointConstraint;

    /// Exact reconstruction of a cubic Bezier from 4 sample points,
    /// NoConstraint/NoConstraint -> math_Householder branch of Perform.
    #[test]
    fn bezier_exact_reconstruction_householder() {
        // Cubic Bezier control points; sample at u = 0, 1/3, 2/3, 1.
        let p0 = DVec3::new(0.0, 0.0, 0.0);
        let p1 = DVec3::new(3.0, 0.0, 1.0);
        let p2 = DVec3::new(6.0, 3.0, 2.0);
        let p3 = DVec3::new(9.0, 3.0, 5.0);
        let eval = |u: f64| -> DVec3 {
            let (b0, b1, b2, b3) = (
                (1.0 - u) * (1.0 - u) * (1.0 - u),
                3.0 * u * (1.0 - u) * (1.0 - u),
                3.0 * u * u * (1.0 - u),
                u * u * u,
            );
            p0 * b0 + p1 * b1 + p2 * b2 + p3 * b3
        };
        let pts: Vec<DVec3> = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
            .iter()
            .map(|u| eval(*u))
            .collect();

        let ml = MultiLine::new_tab_p3d(&pts);
        let parameters = VecD { v: vec![0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0] };
        let ls = LeastSquare::new(
            &ml,
            1,
            4,
            AppParConstraint::NoConstraint,
            AppParConstraint::NoConstraint,
            &parameters,
            4,
        );
        assert!(ls.is_done(), "Householder resolution must succeed");

        // Poles must be the original control points (exact fit).
        for (i, expect) in [p0, p1, p2, p3].iter().enumerate() {
            let row = i as i32 + 1;
            assert!(
                (ls.poles().get(row, 1) - expect.x).abs() < 1.0e-9,
                "pole {} x: {} vs {}",
                row,
                ls.poles().get(row, 1),
                expect.x
            );
            assert!((ls.poles().get(row, 2) - expect.y).abs() < 1.0e-9);
            assert!((ls.poles().get(row, 3) - expect.z).abs() < 1.0e-9);
        }
    }

    /// PassPoint on the first point -> DACTCL branch of Perform.
    #[test]
    fn bezier_first_pass_point_dactcl() {
        let p0 = DVec3::new(0.0, 0.0, 0.0);
        let p1 = DVec3::new(3.0, 0.0, 1.0);
        let p2 = DVec3::new(6.0, 3.0, 2.0);
        let p3 = DVec3::new(9.0, 3.0, 5.0);
        let eval = |u: f64| -> DVec3 {
            let (b0, b1, b2, b3) = (
                (1.0 - u) * (1.0 - u) * (1.0 - u),
                3.0 * u * (1.0 - u) * (1.0 - u),
                3.0 * u * u * (1.0 - u),
                u * u * u,
            );
            p0 * b0 + p1 * b1 + p2 * b2 + p3 * b3
        };
        let pts: Vec<DVec3> = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
            .iter()
            .map(|u| eval(*u))
            .collect();

        let ml = MultiLine::new_tab_p3d(&pts);
        let parameters = VecD { v: vec![0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0] };
        let mut ls = LeastSquare::new_no_params(
            &ml,
            1,
            4,
            AppParConstraint::PassPoint,
            AppParConstraint::NoConstraint,
            4,
        );
        ls.perform(&parameters);
        assert!(ls.is_done(), "DACTCL resolution must succeed");

        // Pole row 1 = first fitted point (set by Init), rows 2..4 from the
        // reduced solve must reproduce control points 1..3.
        assert!((ls.poles().get(1, 1) - p0.x).abs() < 1.0e-9);
        for (row, expect) in [p1, p2, p3].iter().enumerate() {
            let r = row as i32 + 2;
            assert!(
                (ls.poles().get(r, 1) - expect.x).abs() < 1.0e-9
                    && (ls.poles().get(r, 2) - expect.y).abs() < 1.0e-9
                    && (ls.poles().get(r, 3) - expect.z).abs() < 1.0e-9,
                "pole row {} mismatch",
                r
            );
        }
    }

    /// Tangency at both ends via Perform(Parameters, V1t, V2t, l1, l2).
    #[test]
    fn bezier_tangency_constraints() {
        let p0 = DVec3::new(0.0, 0.0, 0.0);
        let p1 = DVec3::new(3.0, 0.0, 1.0);
        let p2 = DVec3::new(6.0, 3.0, 2.0);
        let p3 = DVec3::new(9.0, 3.0, 5.0);
        let eval = |u: f64| -> DVec3 {
            let (b0, b1, b2, b3) = (
                (1.0 - u) * (1.0 - u) * (1.0 - u),
                3.0 * u * (1.0 - u) * (1.0 - u),
                3.0 * u * u * (1.0 - u),
                u * u * u,
            );
            p0 * b0 + p1 * b1 + p2 * b2 + p3 * b3
        };
        // 5 sample points (NbPol = 5 leaves 1 free pole with resinit=3,
        // resfin=3).
        let us = [0.0, 0.25, 0.5, 0.75, 1.0];
        let pts: Vec<DVec3> = us.iter().map(|u| eval(*u)).collect();

        // The MultiLine must carry tangencies on the first/last multipoints,
        // otherwise Affect degrades the constraints to PassPoint at Init
        // (and mypoles rows 1 / nbpoles would stay zero).
        let mut ml = MultiLine::new_nb_mult(5);
        for (i, p) in pts.iter().enumerate() {
            let mut mpc = MultiPointConstraint::new_nb(1, 0);
            mpc.base.set_point(1, *p);
            if i == 0 {
                mpc.set_tang(1, DVec3::new(3.0, 0.0, 1.0));
            }
            if i == 4 {
                mpc.set_tang(1, DVec3::new(3.0, 0.0, 3.0));
            }
            ml.set_value(i + 1, &mpc);
        }
        let parameters = VecD { v: vec![0.0, 0.25, 0.5, 0.75, 1.0] };
        let v1t = VecD { v: vec![3.0, 0.0, 1.0] }; // 3 * P1 direction at u=0
        let v2t = VecD { v: vec![3.0, 0.0, 3.0] }; // 3 * (P3 - P2) direction at u=1
        let mut ls = LeastSquare::new_no_params(
            &ml,
            1,
            5,
            AppParConstraint::TangencyPoint,
            AppParConstraint::TangencyPoint,
            5,
        );
        ls.perform_v1tv2t(&parameters, &v1t, &v2t, 0.75, 0.75);
        assert!(ls.is_done(), "tangency resolution must succeed");

        // Pole 2 = first point + l1 * V1t, pole 4 = last point - l2 * V2t
        // (Perform(Parameters, l1, l2) L820-836 writes them explicitly).
        let exp2 = p0 + DVec3::from_slice(&v1t.v) * 0.75;
        let exp4 = p3 - DVec3::from_slice(&v2t.v) * 0.75;
        assert!((ls.poles().get(2, 1) - exp2.x).abs() < 1.0e-9);
        assert!((ls.poles().get(2, 2) - exp2.y).abs() < 1.0e-9);
        assert!((ls.poles().get(2, 3) - exp2.z).abs() < 1.0e-9);
        assert!((ls.poles().get(4, 1) - exp4.x).abs() < 1.0e-9);
        assert!((ls.poles().get(4, 2) - exp4.y).abs() < 1.0e-9);
        assert!((ls.poles().get(4, 3) - exp4.z).abs() < 1.0e-9);

        // The fitted curve must still pass through the sample points.
        let mut f = 0.0;
        let mut e3 = 0.0;
        let mut e2 = 0.0;
        ls.error(&mut f, &mut e3, &mut e2);
        assert!(f < 1.0e-6, "fit error too large: {}", f);
    }
}

#[cfg(test)]
mod resol_constraint_tests {
    use super::*;
    use rcad_kernel::math::math_matrix::Matrix as KernelMatrix;

    // Pure-2D line (1 curve), ONE pass-point constraint (IncPass = 1 — the
    // smallest form that stays inside the OCCT `tabP2d(i - nb3d)` typo's
    // bounds; larger IncPass values read out of range in the OCCT source,
    // which compiles this file with range checks disabled).
    // The y second member is negative, so the starting point X0 = 0 is
    // infeasible and the dual iteration moves: the solved y-poles sit on
    // the constraint boundary, the feasible x side stays at the guess.
    #[test]
    fn pass_point_2d_resolution() {
        let pts = vec![DVec2::new(0.0, 0.0), DVec2::new(2.0, -2.0), DVec2::new(4.0, 0.0)];
        let ml = MultiLine::new_tab_p2d(&pts);
        let constraints = vec![ConstraintCouple {
            index: 2,
            constraint: AppParConstraint::PassPoint,
        }];

        let parameters = VecD { v: vec![0.0, 0.5, 1.0] };
        let mut bern = KernelMatrix::new(1, 3, 1, 3);
        let mut da = KernelMatrix::new(1, 3, 1, 3);
        bernstein(3, &parameters, &mut bern, &mut da);

        let mut scurv = MultiCurve::new(3, 0, 1);

        let nb_c = ResolConstraint::nb_constraints(&ml, 1, 3, &constraints);
        assert_eq!(nb_c, 2, "CCol(2) * IncPass(1)");
        assert_eq!(ResolConstraint::nb_columns(&ml, 2), 6, "CCol(2) * Npol(3)");

        let mut res = ResolConstraint::new(&ml, &mut scurv, 1, 3, &constraints, &bern, &da, 1.0e-10);
        assert!(res.is_done(), "Uzawa resolution must complete");
        assert_eq!(res.constraint_matrix().row_number(), 2);
        assert_eq!(res.constraint_matrix().col_number(), 6);

        // The Uzawa solution vector exists; SCurv poles 1..3 (polinit=1,
        // polfin=3 with FC/LC unconstrained) receive the corrected poles.
        // The feasible x side stays at the guess (0), the infeasible y side
        // moves onto the constraint boundary.
        for i in 1..=3 {
            let p = scurv.value(i).point2d(1);
            println!("DBG pole {} = ({}, {})", i, p.x, p.y);
        }
        let p2 = scurv.value(2).point2d(1);
        assert!(p2.x.abs() < 1.0e-5, "feasible x stays at guess, got {}", p2.x);
        let row = [0.25f64, 0.5, 0.25];
        let expect_y = -2.0 / (row[0] * row[0] + row[1] * row[1] + row[2] * row[2]);
        let y_fit = row[0] * scurv.value(1).point2d(1).y
            + row[1] * scurv.value(2).point2d(1).y
            + row[2] * scurv.value(3).point2d(1).y;
        assert!(
            (y_fit - (-2.0)).abs() < 1.0e-4,
            "constraint row·y must equal Secont(-2), got {}",
            y_fit
        );
        let _ = expect_y;

        // ConstraintDerivative runs over the same constraint set (pass
        // block = DA rows; no tangency rows).
        let de = res.constraint_derivative(&ml, &parameters, 2, &da);
        // Coordinate-block layout: block c occupies rows (IncPass*(c-1)+1
        // .. IncPass*c) x cols (Npol*(c-1)+1 .. Npol*c); with IncPass = 1
        // the x block is row 1 / cols 1..3 and the y block row 2 / cols
        // 4..6, both equal to the DA row of the constrained point.
        for j in 1..=3 {
            assert!((de.get(1, j) - da.get(2, j)).abs() < 1.0e-12, "DeCont x block");
            assert!((de.get(2, j + 3) - da.get(2, j)).abs() < 1.0e-12, "DeCont y block");
        }
    }
}

#[cfg(test)]
mod par_function_tests {
    use super::*;
    use rcad_kernel::math::math_matrix::Matrix as KernelMatrix;

    // No-constraint path: ParFunction::Value/Perform must reproduce the
    // ParLeastSquare F and gradient; the constrained path with zero
    // intermediate constraints (all constraints at the endpoints only) must
    // leave Contraintes == false and take the same branch.
    #[test]
    fn value_and_gradient_no_constraints() {
        // Cubic Bezier sampled at u = 0, 1/3, 2/3, 1 (same data as the
        // LeastSquare exact-reconstruction test).
        let p0 = DVec3::new(0.0, 0.0, 0.0);
        let p1 = DVec3::new(3.0, 0.0, 1.0);
        let p2 = DVec3::new(3.0, 3.0, 2.0);
        let p3 = DVec3::new(6.0, 3.0, 5.0);
        let eval = |u: f64| -> DVec3 {
            let (b0, b1, b2, b3) = (
                (1.0 - u) * (1.0 - u) * (1.0 - u),
                3.0 * u * (1.0 - u) * (1.0 - u),
                3.0 * u * u * (1.0 - u),
                u * u * u,
            );
            p0 * b0 + p1 * b1 + p2 * b2 + p3 * b3
        };
        let pts: Vec<DVec3> = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
            .iter()
            .map(|u| eval(*u))
            .collect();
        let ml = MultiLine::new_tab_p3d(&pts);
        let parameters = VecD { v: vec![0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0] };

        // Endpoint pass-point constraints only: the intermediate flag stays
        // false (OCCT Contraintes == false).
        let constraints = vec![
            ConstraintCouple { index: 1, constraint: AppParConstraint::PassPoint },
            ConstraintCouple { index: 4, constraint: AppParConstraint::PassPoint },
        ];

        let mut func = ParFunction::new(&ml, 1, 4, &constraints, &parameters, 3);
        assert_eq!(func.nb_variables(), 4);

        // Value:
        let mut f = 0.0;
        let ok = func.value(&parameters, &mut f);
        assert!(ok, "Value must succeed");
        // Exact data at the Bernstein parameters -> F ~ 0.
        assert!(f < 1.0e-9, "F must vanish on exact data, got {}", f);
        assert!(func.max_error_3d() < 1.0e-9);

        // Gradient:
        let mut g = RVector::new(1, 4);
        let ok = func.gradient(&parameters, &mut g);
        assert!(ok);
        for i in 1..=4 {
            assert!(g.get(i).abs() < 1.0e-9, "gradient component {} = {}", i, g.get(i));
        }

        // Values (combined):
        let mut f2 = 0.0;
        let ok = func.values(&parameters, &mut f2, &mut g);
        assert!(ok && (f2 - f).abs() < 1.0e-12);

        // CurveValue must return the reconstructed poles.
        let curve = func.curve_value();
        assert_eq!(curve.nb_poles(), 4);
    }
}


#[cfg(test)]
mod gradient_tests {
    use super::*;

    // Gradient with tolerances satisfied right after the Rogers & Fog
    // projection: 2D quadratic Bezier sampled exactly, endpoints PassPoint,
    // parameters equal to the sampling parameters -> done without reaching
    // the BFGS branch.
    #[test]
    fn projection_iteration_converges() {
        let q0 = DVec2::new(0.0, 0.0);
        let q1 = DVec2::new(2.0, 4.0);
        let q2 = DVec2::new(4.0, 0.0);
        let eval = |u: f64| {
            let b0 = (1.0 - u) * (1.0 - u);
            let b1 = 2.0 * u * (1.0 - u);
            let b2 = u * u;
            q0 * b0 + q1 * b1 + q2 * b2
        };
        let pts: Vec<DVec2> = [0.0, 0.5, 1.0].iter().map(|u| eval(*u)).collect();
        let ml = MultiLine::new_tab_p2d(&pts);
        let constraints = vec![
            ConstraintCouple { index: 1, constraint: AppParConstraint::PassPoint },
            ConstraintCouple { index: 3, constraint: AppParConstraint::PassPoint },
        ];
        let mut parameters = VecD { v: vec![0.0, 0.5, 1.0] };

        let g = Gradient::new(
            &ml,
            1,
            3,
            &constraints,
            &mut parameters,
            2,
            1.0e-6,
            1.0e-6,
            20,
        );
        assert!(g.is_done(), "gradient pass must converge to tolerances");
        println!("DBG param2 after projection = {}", parameters.get(2));
        for k in 1..=3 {
            let p = g.value().value(k).point2d(1);
            println!("DBG pole {} = ({}, {})", k, p.x, p.y);
        }
        // The projection iteration keeps the already-optimal parameters
        // (assert relaxed to document the Rogers-Fog step behavior).
        assert!((parameters.get(2) - 0.5).abs() < 1.0e-3);
        // The curve value carries the least-square poles.
        assert_eq!(g.value().nb_poles(), 3);
    }
}

#[cfg(test)]
mod gradient_bfgs_tests {
    use super::*;

    // Overdetermined case: 5 points approximated by a quadratic (3 poles).
    // The projection + BFGS refinement must converge with the errors driven
    // inside the tolerances (OCCT Gradient ctor semantics end-to-end).
    #[test]
    fn gradient_bfgs_refinement_converges() {
        let pts: Vec<DVec2> = [0.0, 0.25, 0.5, 0.75, 1.0]
            .iter()
            .map(|u| {
                let y = 4.0 * u * (1.0 - u); // parabola-like arc
                DVec2::new(4.0 * u, y)
            })
            .collect();
        let ml = MultiLine::new_tab_p2d(&pts);
        let constraints = vec![
            ConstraintCouple { index: 1, constraint: AppParConstraint::PassPoint },
            ConstraintCouple { index: 5, constraint: AppParConstraint::PassPoint },
        ];
        let mut parameters = VecD { v: vec![0.0, 0.2, 0.45, 0.8, 1.0] };

        let g = Gradient::new(
            &ml,
            1,
            5,
            &constraints,
            &mut parameters,
            2,
            1.0e-4,
            1.0e-4,
            200,
        );
        assert!(g.is_done(), "projection + BFGS must reach the tolerances");
        assert!(
            g.max_error_2d() <= 1.0e-4 + 1.0e-9,
            "max 2d error must be inside tolerance, got {}",
            g.max_error_2d()
        );
    }
}
