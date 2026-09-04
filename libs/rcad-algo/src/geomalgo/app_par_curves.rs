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
use rcad_kernel::math::math_recipes::{dactcl_decompose, dactcl_solve, MATH_STATUS_OK};
use rcad_kernel::math::{IntVec, MatD, VecD};

use glam::{DVec2, DVec3};

use super::app_def::{my_line_tool, MultiLine};
use super::approx_int::{AppParConstraint, MultiBSpCurve, MultiCurve, MultiPoint};

/// OCCT MinPivot default of DACTCL_Decompose / DACTCL_Solve (1.0e-20).
const DACTCL_MIN_PIVOT: f64 = 1.0e-20;

// ---------------------------------------------------------------------------
// math_Matrix / math_Vector / math_IntegerVector wrappers
// ---------------------------------------------------------------------------

/// OCCT math_Matrix — 1-based storage with arbitrary LowerRow/LowerCol.
#[derive(Debug, Clone)]
struct Matrix {
    data: MatD,
    lower_row: i32,
    lower_col: i32,
}

impl Matrix {
    /// OCCT math_Matrix(I1, I2, J1, J2) (zero-initialized storage).
    fn new(i1: i32, i2: i32, j1: i32, j2: i32) -> Self {
        assert!(i2 >= i1 && j2 >= j1, "math_Matrix: bad range");
        Matrix {
            data: MatD::new((i2 - i1 + 1) as usize, (j2 - j1 + 1) as usize),
            lower_row: i1,
            lower_col: j1,
        }
    }

    /// OCCT math_Matrix(I1, I2, J1, J2, InitValue).
    fn new_init(i1: i32, i2: i32, j1: i32, j2: i32, init: f64) -> Self {
        let mut m = Matrix::new(i1, i2, j1, j2);
        m.init(init);
        m
    }

    /// OCCT Init(Value).
    fn init(&mut self, init: f64) {
        for r in 1..=self.row_number() {
            for c in 1..=self.col_number() {
                self.data.m[(r - 1) as usize][(c - 1) as usize] = init;
            }
        }
    }

    #[inline]
    fn get(&self, i: i32, j: i32) -> f64 {
        self.data.get(
            (i - self.lower_row + 1).max(1) as usize,
            (j - self.lower_col + 1).max(1) as usize,
        )
    }

    #[inline]
    fn set(&mut self, i: i32, j: i32, v: f64) {
        self.data.set(
            (i - self.lower_row + 1).max(1) as usize,
            (j - self.lower_col + 1).max(1) as usize,
            v,
        );
    }

    /// OCCT RowNumber().
    fn row_number(&self) -> i32 {
        self.data.n_rows() as i32
    }

    /// OCCT ColNumber().
    fn col_number(&self) -> i32 {
        self.data.n_cols() as i32
    }

    /// OCCT LowerRow().
    fn lower_row(&self) -> i32 {
        self.lower_row
    }

    /// Normalized (1-based) storage view, for the rcad Householder whose
    /// entry points take MatD.
    fn data(&self) -> &MatD {
        &self.data
    }

    fn data_mut(&mut self) -> &mut MatD {
        &mut self.data
    }
}

/// OCCT math_Vector — 1-based storage with arbitrary lower bound.
#[derive(Debug, Clone)]
struct RVector {
    data: VecD,
    lower: i32,
}

impl RVector {
    /// OCCT math_Vector(I1, I2).
    fn new(i1: i32, i2: i32) -> Self {
        assert!(i2 >= i1, "math_Vector: bad range");
        RVector {
            data: VecD::new((i2 - i1 + 1) as usize),
            lower: i1,
        }
    }

    /// OCCT math_Vector(I1, I2, InitValue).
    fn new_init(i1: i32, i2: i32, init: f64) -> Self {
        let mut v = RVector::new(i1, i2);
        for r in 1..=v.data.len() {
            v.data.set(r, init);
        }
        v
    }

    #[inline]
    fn get(&self, i: i32) -> f64 {
        self.data.get((i - self.lower + 1).max(1) as usize)
    }

    #[inline]
    fn set(&mut self, i: i32, v: f64) {
        self.data.set((i - self.lower + 1).max(1) as usize, v);
    }

    /// OCCT Length().
    fn length(&self) -> i32 {
        self.data.len() as i32
    }

    /// OCCT Upper().
    fn upper(&self) -> i32 {
        self.lower + self.length() - 1
    }
}

/// OCCT math_IntegerVector — 1-based storage with arbitrary lower bound.
#[derive(Debug, Clone)]
struct IVector {
    data: IntVec,
    lower: i32,
}

impl IVector {
    /// OCCT math_IntegerVector(I1, I2, InitValue).
    fn new_init(i1: i32, i2: i32, init: i32) -> Self {
        assert!(i2 >= i1, "math_IntegerVector: bad range");
        let mut v = IVector {
            data: IntVec::new((i2 - i1 + 1) as usize),
            lower: i1,
        };
        for r in 1..=v.data.len() {
            v.data.set(r, init);
        }
        v
    }

    #[inline]
    fn get(&self, i: i32) -> i32 {
        self.data.get((i - self.lower + 1).max(1) as usize)
    }

    #[inline]
    fn set(&mut self, i: i32, v: i32) {
        self.data.set((i - self.lower + 1).max(1) as usize, v);
    }
}

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
