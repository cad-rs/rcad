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

use rcad_kernel::math::bspl_lib::locate_parameter_flat;
use rcad_kernel::math::math_householder::Householder;
use rcad_kernel::math::math_matrix::{IntegerVector as IVector, Matrix, Vector as RVector};
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

    /// OCCT BezierValue() (gxx L1402-1407).
    pub fn bezier_value(&self) -> MultiCurve {
        if self.myknots.is_some() {
            panic!("Standard_NoSuchObject: AppParCurves_LeastSquare::BezierValue");
        }
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
